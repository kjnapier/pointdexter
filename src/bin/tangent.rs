use pointdexter::cli::{Config, Cli};
use pointdexter::Exposure;
use pointdexter::load_detections;
use pointdexter::xyz_to_proj_matrix;
use pointdexter::load_initial_conditions;

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




pub fn build_grid(alphas: &Vec<f64>, betas: &Vec<f64>, pixel_scale: f64) -> Vec<Vec<usize>> {
    let mut alpha_min = std::f64::INFINITY;
    let mut alpha_max = std::f64::NEG_INFINITY;
    let mut beta_min = std::f64::INFINITY;;
    let mut beta_max = std::f64::NEG_INFINITY;
    for (alpha, beta) in alphas.iter().zip(betas.iter()) {
        if *alpha < alpha_min {
            alpha_min = *alpha;
        }
        if *alpha > alpha_max {
            alpha_max = *alpha;
        }
        if *beta < beta_min {
            beta_min = *beta;
        }
        if *beta > beta_max {
            beta_max = *beta;
        }
    }
    let alpha_bins = ((alpha_max - alpha_min)/pixel_scale).ceil() as usize;
    let beta_bins = ((beta_max - beta_min)/pixel_scale).ceil() as usize;
    let mut grid = vec![vec![0; beta_bins]; alpha_bins];
    for (alpha, beta) in alphas.iter().zip(betas.iter()) {
        let alpha_idx = ((*alpha - alpha_min)/pixel_scale).floor() as usize;
        let beta_idx = ((*beta - beta_min)/pixel_scale).floor() as usize;
        //println!("alpha: {}, beta: {}, alpha_idx: {}, beta_idx: {}", alpha, beta, alpha_idx, beta_idx);
        grid[alpha_idx][beta_idx] += 1;
    }
    //println!("Grid size: {} x {}", alpha_bins, beta_bins);
    grid
}

// Return grid and the parameters needed to convert back to alpha and beta values.
pub fn build_grid_and_wcs(alphas: &Vec<f64>, betas: &Vec<f64>, pixel_scale: f64) -> (Vec<Vec<usize>>, f64, f64, f64) {
    let mut alpha_min = std::f64::INFINITY;
    let mut alpha_max = std::f64::NEG_INFINITY;
    let mut beta_min = std::f64::INFINITY;;
    let mut beta_max = std::f64::NEG_INFINITY;
    for (alpha, beta) in alphas.iter().zip(betas.iter()) {
        if *alpha < alpha_min {
            alpha_min = *alpha;
        }
        if *alpha > alpha_max {
            alpha_max = *alpha;
        }
        if *beta < beta_min {
            beta_min = *beta;
        }
        if *beta > beta_max {
            beta_max = *beta;
        }
    }
    let alpha_bins = ((alpha_max - alpha_min)/pixel_scale).ceil() as usize;
    let beta_bins = ((beta_max - beta_min)/pixel_scale).ceil() as usize;
    let mut grid = vec![vec![0; beta_bins]; alpha_bins];
    for (alpha, beta) in alphas.iter().zip(betas.iter()) {
        let alpha_idx = ((*alpha - alpha_min)/pixel_scale).floor() as usize;
        let beta_idx = ((*beta - beta_min)/pixel_scale).floor() as usize;
        //println!("alpha: {}, beta: {}, alpha_idx: {}, beta_idx: {}", alpha, beta, alpha_idx, beta_idx);
        grid[alpha_idx][beta_idx] += 1;
    }
    //println!("Grid size: {} x {}", alpha_bins, beta_bins);
    (grid, alpha_min, beta_min, pixel_scale)
}

pub type CountsMap = HashMap<u64, u32, RandomState>;
// Return sparse counts + WCS params + bin counts
pub fn build_grid_sparse(alphas: &[f64], betas: &[f64], pixel_scale: f64) -> (CountsMap, f64, f64, f64, usize, usize) {
    debug_assert_eq!(alphas.len(), betas.len());
    let n = alphas.len();
    if n == 0 {
        return (CountsMap::with_hasher(RandomState::new()), 0.0, 0.0, pixel_scale, 0, 0);
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
    let alpha_bins = ((alpha_max - alpha_min) * inv).ceil() as usize;
    let beta_bins  = ((beta_max  - beta_min)  * inv).ceil() as usize;

    
    // let mut counts: HashMap<u64, u32> = HashMap::with_capacity(n * 2);
    let mut counts: CountsMap = CountsMap::with_capacity_and_hasher(n * 2, RandomState::new());


    for (&a, &b) in alphas.iter().zip(betas.iter()) {
        let mut ai = ((a - alpha_min) * inv) as isize; // trunc toward 0; OK since >=0 typically
        let mut bi = ((b - beta_min)  * inv) as isize;

        // Clamp just in case of edge/rounding issues
        if ai < 0 { ai = 0; }
        if bi < 0 { bi = 0; }
        if ai as usize >= alpha_bins { ai = alpha_bins.saturating_sub(1) as isize; }
        if bi as usize >= beta_bins  { bi = beta_bins.saturating_sub(1)  as isize; }

        let key = pack(ai as u32, bi as u32);
        *counts.entry(key).or_insert(0) += 1;
    }

    (counts, alpha_min, beta_min, pixel_scale, alpha_bins, beta_bins)
}


#[inline]
fn pack(i: u32, j: u32) -> u64 {
    ((i as u64) << 32) | (j as u64)
}

#[inline]
fn unpack(key: u64) -> (u32, u32) {
    ((key >> 32) as u32, (key & 0xFFFF_FFFF) as u32)
}

/// Sum counts in the 3x3 neighborhood around each peak (including the center).
pub fn convolve_peaks_sparse(
    counts: &CountsMap,
    peaks: &[(u32, u32)],
    alpha_bins: u32,
    beta_bins: u32,
) -> Vec<u32> {
    let mut out = Vec::with_capacity(peaks.len());

    for &(i, j) in peaks {
        let mut sum: u32 = 0;

        // neighborhood bounds (avoid per-neighbor checks where possible)
        let i0 = i.saturating_sub(1);
        let j0 = j.saturating_sub(1);
        let i1 = (i + 1).min(alpha_bins.saturating_sub(1));
        let j1 = (j + 1).min(beta_bins.saturating_sub(1));

        for ni in i0..=i1 {
            for nj in j0..=j1 {
                if let Some(&v) = counts.get(&pack(ni, nj)) {
                    sum += v;
                }
            }
        }

        out.push(sum);
    }

    out
}

/// Find local maxima among occupied cells, skipping any cell with count <= min_count.
/// A cell is a peak if no neighbor in its 8-neighborhood has a strictly larger count.
/// (Ties are allowed.)
pub fn find_local_maxima_sparse(
    counts: &CountsMap,
    alpha_bins: u32,
    beta_bins: u32,
    min_count: u32, // pass 3 to "skip <= 2"
) -> Vec<(u32, u32)> {
    let mut peaks: Vec<(u32, u32)> = Vec::new();

    for (&key, &count) in counts.iter() {
        if count < min_count {
            continue;
        }

        let (i, j) = unpack(key);

        let i0 = i.saturating_sub(1);
        let j0 = j.saturating_sub(1);
        let i1 = (i + 1).min(alpha_bins.saturating_sub(1));
        let j1 = (j + 1).min(beta_bins.saturating_sub(1));

        let mut is_peak = true;

        'nbrs: for ni in i0..=i1 {
            for nj in j0..=j1 {
                if ni == i && nj == j {
                    continue;
                }
                if let Some(&v) = counts.get(&pack(ni, nj)) {
                    if v > count {
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
pub fn save_grid_png(grid: &Vec<Vec<usize>>, path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use image::{GrayImage, Luma};
    let alpha_bins = grid.len();
    let beta_bins = grid[0].len();
    let mut img = GrayImage::new(beta_bins as u32, alpha_bins as u32);
    for (i, row) in grid.iter().enumerate() {
        for (j, count) in row.iter().enumerate() {
            let intensity = (*count as f64 / 10.0).min(1.0) * 255.0; // scale counts to [0, 255], adjust as needed
            img.put_pixel(j as u32, i as u32, Luma([intensity as u8]));
        }
    }
    img.save(path)?;
    Ok(())
}

pub fn convolve_peaks(grid: &Vec<Vec<usize>>, peaks: &Vec<(usize, usize)>) -> Vec<usize> {
    let mut counts: Vec<usize> = Vec::with_capacity(peaks.len());
    for peak in peaks.iter() {
        let i = peak.0;
        let j = peak.1;
        let mut count_sum = 0;
        for di in -1..=1 {
            for dj in -1..=1 {
                let ni = i as isize + di;
                let nj = j as isize + dj;
                if ni < 0 || ni >= grid.len() as isize || nj < 0 || nj >= grid[i].len() as isize {
                    continue;
                }
                count_sum += grid[ni as usize][nj as usize];
            }
        }
        counts.push(count_sum);
    } 
    counts
}


pub fn find_local_maxima(grid: &Vec<Vec<usize>>) -> Vec<(usize, usize)> {
    let mut peaks: Vec<(usize, usize)> = Vec::new();
    for i in 0..grid.len() {
        for j in 0..grid[i].len() {
            let count = grid[i][j];
            if count == 0 {
                continue;
            }
            let mut is_peak = true;
            for di in -1..=1 {
                for dj in -1..=1 {
                    if di == 0 && dj == 0 {
                        continue;
                    }   
                    let ni = i as isize + di;
                    let nj = j as isize + dj;
                    if ni < 0 || ni >= grid.len() as isize || nj < 0 || nj >= grid[i].len() as isize {
                        continue;
                    }
                    if grid[ni as usize][nj as usize] > count {
                        is_peak = false;
                        break;
                    }
                }
                if !is_peak {
                    break;
                }
            }
            if is_peak {
                peaks.push((i, j));
            }
        }
    }
    peaks
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
        // This can be calculated up front. 
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

    // Let's read a csv file.
    let mut rdr = csv::Reader::from_path(&config.initial_conditions_file)?;

    // Iterate over each record and save the values to variables. The csv file should have columns: r, vr, vo, psi, id
    let mut initial_conditions: Vec<(f64, f64, f64, f64, String)> = Vec::new();
    for result in rdr.records() {
        // The iterator yields Result<StringRecord, Error>, so we check the error
        let record = result?;
        let r = record[0].parse::<f64>()?;
        let vr = record[1].parse::<f64>()?;
        let adot = record[2].parse::<f64>()?;
        let bdot = record[3].parse::<f64>()?;
        let id = record[4].to_string();        
        initial_conditions.push((r, vr, adot, bdot, id));
    }

    // // Print the initial conditions to check that we read them correctly.
    // println!("Initial conditions:");
    // for (r, vr, vo, psi, id) in initial_conditions.iter() {
    //     println!("id: {}, r: {}, vr: {}, vo: {}, psi: {}", id, r, vr, vo, psi);
    // }
    
    let epsilon = config.epsilon_arcsec;
    let pixel_scale = epsilon / 3.0 * std::f64::consts::PI / (180.0 * 3600.0); // convert arcsec to radians
    
    let ref_epoch = config.reference_epoch;
    let time_span = 60.0; // days; TODO: make this configurable    
    let gamma = 1.0/30.0; 
    let gdot = 0.0; 


    let start_time = std::time::Instant::now();

    //let adot = 24.0 * (0.5/3600.0) * std::f64::consts::PI/180.0; // 0.5 arcsec/day in radians/day 
    //let bdot = 24.0 * (0.1/3600.0) * std::f64::consts::PI/180.0; // 0.1 arcsec/day in radians/day

    // Now we can loop over our grid of adot and bdot values, 
    // sync the detections to the reference epoch for each pair, 
    // build a grid in alpha-beta space, and count the number of detections that fall into each cell.
    
    let min_count = 3;

    for (r, vr, adot, bdot, id) in initial_conditions.iter() {
        let gamma = 1.0/r;
        let gdot = vr/r;
            
        let (alphas, betas) = sync_detections(&tangent_exposures, mu, gamma, gdot, *adot, *bdot, ref_epoch, count);
        let (counts, alpha_min, beta_min, pixel_scale, alpha_bins, beta_bins) = build_grid_sparse(&alphas, &betas, pixel_scale);
        let peaks = find_local_maxima_sparse(&counts, alpha_bins as u32, beta_bins as u32, min_count);
        let convolved_counts = convolve_peaks_sparse(&counts, &peaks, alpha_bins as u32, beta_bins as u32);
        let mut peak_counts: Vec<((u32, u32), u32)> = peaks.into_iter().zip(convolved_counts.into_iter()).collect();
        peak_counts.sort_by(|a, b| b.1.cmp(&a.1)); // sort in descending order of convolved count
        let threshold = 12;
        let significant_peaks = peak_counts.iter().filter(|(_, count)| *count >= threshold).count();
        if significant_peaks > 0 {
            // Print the adot, bdot, and number of significant peaks for this pair.
            // Units should be in arcsec/day for readability.
            let adot_arcsec_per_day = adot * 3600.0 * 180.0 / std::f64::consts::PI;
            let adot_arcsec_per_hour = adot_arcsec_per_day / 24.0;
            let bdot_arcsec_per_day = bdot * 3600.0 * 180.0 / std::f64::consts::PI;
            let bdot_arcsec_per_hour = bdot_arcsec_per_day / 24.0;
            // Find the top peak for this pair and print its alpha, beta, and convolved count.
            let top_peak = peak_counts.iter().max_by_key(|(_, count)| *count).unwrap();
            let alpha = (top_peak.0 .0 as f64 + 0.5) * pixel_scale + alpha_min; // add 0.5 to get the center of the pixel
            let beta = (top_peak.0 .1 as f64 + 0.5) * pixel_scale + beta_min; // add 0.5 to get the center of the pixel
            println!("adot: {:.3} arcsec/hour, bdot: {:.3} arcsec/hour, top peak: (alpha: {:.6}, beta: {:.6}, convolved count: {}), significant peaks: {}", adot_arcsec_per_hour, bdot_arcsec_per_hour, alpha, beta, top_peak.1, significant_peaks);
            //println!("adot: {:.3} arcsec/hour, bdot: {:.3} arcsec/hour, significant peaks: {}", adot_arcsec_per_hour, bdot_arcsec_per_hour, significant_peaks);

        }
    } 

    let duration = start_time.elapsed();    
    println!("completed in {:?}", duration);

    // // Now sync the detections to the reference epoch and calculate alpha and beta for each detection.
    // let (alphas, betas) = sync_detections(&tangent_exposures, mu, gamma, gdot, adot, bdot, ref_epoch, count);

    // let (counts, alpha_min, beta_min, pixel_scale, alpha_bins, beta_bins) = build_grid_sparse(&alphas, &betas, pixel_scale);
    
    // let peaks = find_local_maxima_sparse(&counts, alpha_bins as u32, beta_bins as u32, min_count);
    // let convolved_counts = convolve_peaks_sparse(&counts, &peaks, alpha_bins as u32, beta_bins as u32);
    
    // // save_grid_png(&grid, "tangent_plane.png")?;
    

    // // Sort peaks by convolved counts and print the top 10
    // let mut peak_counts: Vec<((u32, u32), u32)> = peaks.into_iter().zip(convolved_counts.into_iter()).collect();
    // peak_counts.sort_by(|a, b| b.1.cmp(&a.1)); // sort in descending order of convolved count
    
    //  println!("Top 10 peaks (alpha, beta, convolved count):");
    // for (peak, count) in peak_counts.iter().take(10) {
    //     let alpha = (peak.0 as f64 + 0.5) * pixel_scale + alpha_min; // add 0.5 to get the center of the pixel
    //     let beta = (peak.1 as f64 + 0.5) * pixel_scale + beta_min; // add 0.5 to get the center of the pixel
    //     println!("({}, {}, {})", alpha, beta, count);
    // }

    // // Count the number of peaks with convolved count above a threshold (e.g. 5)
    // let threshold = 12;
    // let significant_peaks = peak_counts.iter().filter(|(_, count)| *count >= threshold).count();
    // println!("Number of peaks with convolved count >= {}: {}", threshold, significant_peaks);

    // let duration = start_time.elapsed();    
    // println!("completed in {:?}", duration);

    // Now convert the significant peaks back to alpha and beta values and print them
    // println!("Significant peaks (alpha, beta, convolved count):");
    // for (peak, count) in peak_counts.iter().filter(|(_, count)| *count >= threshold) {
    //     let alpha = (peak.0 as f64 + 0.5) * pixel_scale + alpha_min; // add 0.5 to get the center of the pixel
    //     let beta = (peak.1 as f64 + 0.5) * pixel_scale + beta_min; // add 0.5 to get the center of the pixel
    //     println!("({}, {}, {})", alpha, beta, count);
    // }

    
    
   
    Ok(())

    

}