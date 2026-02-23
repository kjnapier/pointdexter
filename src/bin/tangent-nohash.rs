

use pointdexter::cli::{Config, Cli};
use pointdexter::Exposure;
use pointdexter::load_detections;
use pointdexter::xyz_to_proj_matrix;

use spacerocks::{SpiceKernel};
use pointdexter::Detection;
use pointdexter::TangentPlaneExposure;
use nalgebra::Vector3;
use clap::Parser;

use spacerocks::coordinates::Origin;
use spacerocks::time::Time;
use spacerocks::constants::SPEED_OF_LIGHT;

// use std::collections::HashMap;
use hashbrown::HashMap;
use ahash::RandomState;



use std::ops::Range;

#[inline]
fn pack(i: u32, j: u32) -> u64 {
    ((i as u64) << 32) | (j as u64)
}

#[inline]
fn unpack(key: u64) -> (u32, u32) {
    ((key >> 32) as u32, (key & 0xFFFF_FFFF) as u32)
}

#[derive(Clone, Copy, Debug)]
pub struct Cell {
    pub key: u64,   // packed (i,j)
    pub count: u32, // occupancy
}

/// Build a hashless sparse histogram:
/// - returns sorted unique occupied cells with counts
/// - plus WCS params + bin dims
pub fn build_grid_sparse_sorted(
    alphas: &[f64],
    betas: &[f64],
    pixel_scale: f64,
) -> (Vec<Cell>, f64, f64, f64, u32, u32) {
    debug_assert_eq!(alphas.len(), betas.len());
    let n = alphas.len();
    if n == 0 {
        return (Vec::new(), 0.0, 0.0, pixel_scale, 0, 0);
    }

    let mut alpha_min = alphas[0];
    let mut alpha_max = alphas[0];
    let mut beta_min  = betas[0];
    let mut beta_max  = betas[0];

    for i in 1..n {
        let a = alphas[i];
        let b = betas[i];
        if a < alpha_min { alpha_min = a; }
        if a > alpha_max { alpha_max = a; }
        if b < beta_min  { beta_min  = b; }
        if b > beta_max  { beta_max  = b; }
    }

    let inv = 1.0 / pixel_scale;

    // Important: if you want the max value to land in-bounds when it hits the edge,
    // consider +1 here. Your current code "ceil()" often works, but can edge-case.
    let alpha_bins = ((alpha_max - alpha_min) * inv).ceil().max(1.0) as u32;
    let beta_bins  = ((beta_max  - beta_min)  * inv).ceil().max(1.0) as u32;

    let mut keys: Vec<u64> = Vec::with_capacity(n);

    for (&a, &b) in alphas.iter().zip(betas.iter()) {
        // Using floor (via as u32 after clamp) is usually what you want for bins.
        let mut ai = ((a - alpha_min) * inv).floor() as i64;
        let mut bi = ((b - beta_min)  * inv).floor() as i64;

        // clamp
        if ai < 0 { ai = 0; }
        if bi < 0 { bi = 0; }
        if ai >= alpha_bins as i64 { ai = alpha_bins.saturating_sub(1) as i64; }
        if bi >= beta_bins  as i64 { bi = beta_bins.saturating_sub(1)  as i64; }

        keys.push(pack(ai as u32, bi as u32));
    }

    keys.sort_unstable();

    // run-length encode
    let mut cells: Vec<Cell> = Vec::with_capacity(keys.len().min(1024));
    let mut cur = keys[0];
    let mut cnt: u32 = 1;

    for &k in keys.iter().skip(1) {
        if k == cur {
            cnt += 1;
        } else {
            cells.push(Cell { key: cur, count: cnt });
            cur = k;
            cnt = 1;
        }
    }
    cells.push(Cell { key: cur, count: cnt });

    (cells, alpha_min, beta_min, pixel_scale, alpha_bins, beta_bins)
}

/// Row index for fast "look in row i for column j" without hashing.
/// Stores only rows that exist: (row_i, range into `cells`)
#[derive(Debug)]
pub struct RowIndex {
    pub rows: Vec<(u32, Range<usize>)>,
}

pub fn build_row_index(cells: &[Cell]) -> RowIndex {
    let mut rows: Vec<(u32, Range<usize>)> = Vec::new();
    if cells.is_empty() {
        return RowIndex { rows };
    }

    let (mut cur_i, _) = unpack(cells[0].key);
    let mut start = 0usize;

    for (idx, c) in cells.iter().enumerate().skip(1) {
        let (i, _) = unpack(c.key);
        if i != cur_i {
            rows.push((cur_i, start..idx));
            cur_i = i;
            start = idx;
        }
    }
    rows.push((cur_i, start..cells.len()));
    RowIndex { rows }
}

#[inline]
fn row_range<'a>(ri: &'a RowIndex, i: u32) -> Option<Range<usize>> {
    // rows is sorted by i
    match ri.rows.binary_search_by_key(&i, |(row_i, _)| *row_i) {
        Ok(pos) => Some(ri.rows[pos].1.clone()),
        Err(_) => None,
    }
}

#[inline]
fn get_count_in_row(cells: &[Cell], range: Range<usize>, j: u32) -> u32 {
    // In this row slice, keys are sorted by (i,j), so by j within fixed i.
    // We'll binary_search by packed key with same i (i is constant inside row).
    // Instead of repacking, just search by unpacking j with comparator.
    let slice = &cells[range];
    match slice.binary_search_by(|c| {
        let (_, cj) = unpack(c.key);
        cj.cmp(&j)
    }) {
        Ok(pos) => slice[pos].count,
        Err(_) => 0,
    }
}

/// Hashless local maxima:
/// A cell is a peak if no neighbor has strictly larger count.
/// Ties allowed.
pub fn find_local_maxima_sparse_sorted(
    cells: &[Cell],
    row_index: &RowIndex,
    alpha_bins: u32,
    beta_bins: u32,
    min_count: u32,
) -> Vec<(u32, u32)> {
    let mut peaks = Vec::new();

    for c in cells {
        if c.count < min_count {
            continue;
        }
        let (i, j) = unpack(c.key);

        let i0 = i.saturating_sub(1);
        let j0 = j.saturating_sub(1);
        let i1 = (i + 1).min(alpha_bins.saturating_sub(1));
        let j1 = (j + 1).min(beta_bins.saturating_sub(1));

        let mut is_peak = true;

        'nbrs: for ni in i0..=i1 {
            if let Some(rr) = row_range(row_index, ni) {
                for nj in j0..=j1 {
                    if ni == i && nj == j {
                        continue;
                    }
                    let v = get_count_in_row(cells, rr.clone(), nj);
                    if v > c.count {
                        is_peak = false;
                        break 'nbrs;
                    }
                }
            }
        }

        if is_peak {
            peaks.push((i, j));
        }
    }

    peaks
}

/// Hashless 3x3 sum around each peak
pub fn convolve_peaks_sparse_sorted(
    cells: &[Cell],
    row_index: &RowIndex,
    peaks: &[(u32, u32)],
    alpha_bins: u32,
    beta_bins: u32,
) -> Vec<u32> {
    let mut out = Vec::with_capacity(peaks.len());

    for &(i, j) in peaks {
        let i0 = i.saturating_sub(1);
        let j0 = j.saturating_sub(1);
        let i1 = (i + 1).min(alpha_bins.saturating_sub(1));
        let j1 = (j + 1).min(beta_bins.saturating_sub(1));

        let mut sum = 0u32;

        for ni in i0..=i1 {
            if let Some(rr) = row_range(row_index, ni) {
                for nj in j0..=j1 {
                    sum += get_count_in_row(cells, rr.clone(), nj);
                }
            }
        }

        out.push(sum);
    }

    out
}


pub fn sort_detections_into_exposures(detections: &Vec<Detection>) -> Vec<Exposure> {
    // sort the detections by epoch
    let mut detections = detections.clone();
    detections.sort_by(|a, b| a.epoch.partial_cmp(&b.epoch).unwrap());
    let mut exposure_groups: Vec<Vec<Detection>> = Vec::new();
    let mut current_group: Vec<Detection> = Vec::new();
    let mut current_epoch = detections[0].epoch;
    for detection in detections.iter() {
        if (detection.epoch - current_epoch).abs() < 1e-6 {
            current_group.push(detection.clone());
        } else {
            exposure_groups.push(current_group);
            current_group = vec![detection.clone()];
            current_epoch = detection.epoch;
        }
    }
    if !current_group.is_empty() {
        exposure_groups.push(current_group);
    }

    let mut exposures: Vec<Exposure> = Vec::with_capacity(exposure_groups.len());
    for (idx, group) in exposure_groups.iter().enumerate() {
        let exposure = Exposure {
            id: format!("exposure_{}", idx),
            epoch: group[0].epoch,
            filter: None,
            detections: group.iter().map(|d| d.rho_hat).collect(),
            observer_position: group[0].observer_position,
            observer_velocity: group[0].observer_velocity,
            reference_plane: group[0].reference_plane.clone(),
        };
        exposures.push(exposure);
    }
    exposures
}

pub fn transform_exposures_to_tangent_plane(exposures: &Vec<Exposure>, center: Vector3<f64>) -> Vec<TangentPlaneExposure> {
    // do stuff
    exposures.iter().map(|e| e.transform_to_tangent_plane(center)).collect()
}

pub fn sync_detections(exposures: &Vec<TangentPlaneExposure>, mu: f64, gamma: f64, gdot: f64, adot: f64, bdot: f64, ref_epoch: f64, ndets: usize) -> (Vec<f64>, Vec<f64>) {

    let mut alphas = Vec::with_capacity(ndets);
    let mut betas = Vec::with_capacity(ndets); 

    let z0 = 1.0/gamma;
    
    let mm_sqr = (mu * gamma * gamma * gamma);
    for exposure in exposures.iter() {

        // Calculate the mean theta_x and theta_y across all detections in this exposure.
        let mut mean_theta_x = 0.0;
        let mut mean_theta_y = 0.0;
        let mut count = 0;
        for tx in exposure.theta_x.iter() {
            mean_theta_x += tx;
            count += 1;
        }
        for ty in exposure.theta_y.iter() {
            mean_theta_y += ty;
        }
        mean_theta_x /= count as f64;
        mean_theta_y /= count as f64;

        let mut synced_exposure: Vec<(f64, f64)> = Vec::new();
        let t = exposure.epoch - ref_epoch;

        let rho2 = (1.0 + mean_theta_x * mean_theta_x + mean_theta_y * mean_theta_y)*(z0 - exposure.xyz_e.z).powi(2);
        let rho = rho2.sqrt();

        // Improve light time correction at some point.
        let dt = rho/SPEED_OF_LIGHT;

        let t = exposure.epoch - ref_epoch;
        let tp = t - dt;
        let f = (1.0 - 0.5 * mm_sqr * tp * tp);
        let g = tp;
        let fac = (1.0 + g/f * gdot - gamma/f * exposure.xyz_e.z);

        for (theta_x, theta_y) in exposure.theta_x.iter().zip(exposure.theta_y.iter()) {
            let phi_x = theta_x * fac + gamma/f * exposure.xyz_e.x;
            let phi_y = theta_y * fac + gamma/f * exposure.xyz_e.y;
            let alpha = phi_x - g/f * adot;
            let beta = phi_y - g/f * bdot;
            //println!("Exposure {}: alpha = {}, beta = {}", exposure.id, alpha, beta);
            alphas.push(alpha);
            betas.push(beta);
        }
    }
    (alphas, betas)
}




pub fn main() ->  Result<(), Box<dyn std::error::Error>> {
    // We need to make this more flexible.
    let args = Cli::try_parse()?;
    let f = std::fs::File::open(args.config)?;
    let config: Config = serde_yaml::from_reader(f)?;

    let mut kernel = SpiceKernel::new();
    kernel.load_spk(format!("{}/sb441-n16.bsp", config.spice_path).as_str())?;
    kernel.load_spk(format!("{}/de440s.bsp", config.spice_path).as_str())?;
    kernel.load_bpc(format!("{}/earth_1962_240827_2124_combined.bpc", config.spice_path).as_str())?;

    let spacerock_origin = Origin::from_str("ssb")?;
    let mu = spacerock_origin.mu();
    //println!("Loaded SPICE kernels. Solar system barycenter mu: {}", mu);

    let mut detections = load_detections(&config.detection_catalog, &config.orbit_reference_plane, &kernel)?;
    let exposures = sort_detections_into_exposures(&detections);
    //println!("Loaded {} exposures", exposures.len());


    //let ref_vec = exposures[0].detections[0];
    // Calculate the reference vector as the mean of the detections in all exposures.
    let mut ref_vec = Vector3::zeros();
    let mut count = 0;
    for exposure in exposures.iter() {
        for det in exposure.detections.iter() {
            ref_vec += det;
            count += 1;
        }
    }
    ref_vec /= count as f64;
    ref_vec = ref_vec.normalize();

    let mut tangent_exposures = Vec::new();

    for exposure in exposures.iter() {
        let tangent_exposure = exposure.transform_to_tangent_plane(ref_vec);
        tangent_exposures.push(tangent_exposure);
    }

    let epsilon = config.epsilon_arcsec;
    let pixel_scale = epsilon / 3.0 * std::f64::consts::PI / (180.0 * 3600.0); // convert arcsec to radians
    
    let ref_epoch = config.reference_epoch;
    let gamma = 1.0/40.0; 
    let gdot = 0.0; 
    let adot = 24.0 * (0.5/3600.0) * std::f64::consts::PI/180.0; // 0.5 arcsec/day in radians/day 
    let bdot = 24.0 * (0.1/3600.0) * std::f64::consts::PI/180.0; // 0.1 arcsec/day in radians/day
    
    println!("Calculating alpha and beta for each detection...");
    
    let start_time = std::time::Instant::now();
    let (alphas, betas) = sync_detections(&tangent_exposures, mu, gamma, gdot, adot, bdot, ref_epoch, count);

    // Now sync the detections to the reference epoch and calculate alpha and beta for each detection.
    let (cells, alpha_min, beta_min, pixel_scale, alpha_bins, beta_bins) =
        build_grid_sparse_sorted(&alphas, &betas, pixel_scale);

    let row_index = build_row_index(&cells);

    let min_count = 3;
    let peaks = find_local_maxima_sparse_sorted(
        &cells,
        &row_index,
        alpha_bins,
        beta_bins,
        min_count,
    );

    let convolved_counts = convolve_peaks_sparse_sorted(
        &cells,
        &row_index,
        &peaks,
        alpha_bins,
        beta_bins,
    );

    let threshold = 12;
    let significant_peaks: Vec<((u32, u32), u32)> = peaks
        .iter()
        .zip(convolved_counts.iter())
        .filter(|&(_, &c)| c >= threshold)
        .map(|(&p, &c)| (p, c))
        .collect();


    let duration = start_time.elapsed();    
    println!("completed in {:?}", duration);
    


    
    
   
    Ok(())

    

}