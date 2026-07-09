use pointdexter::cli::{Config, Cli};
use pointdexter::Exposure;
use pointdexter::load_detections;

use spacerocks::{SpiceKernel};
use pointdexter::Detection;
use pointdexter::TangentPlaneExposure;
use nalgebra::Vector3;
use clap::Parser;
use rayon::prelude::*;

use spacerocks::coordinates::Origin;
use spacerocks::constants::SPEED_OF_LIGHT;

// Shared back-end reused from matt-0 (all re-exported by pointdexter::utils):
// the coupled αβγ tangent-plane fit + the ecliptic→equatorial output helpers.
use pointdexter::{iterative_reject, ecliptic_to_equatorial, ecliptic_to_equatorial_vector};
use std::io::Write;
use std::collections::HashSet;
use std::f64::consts::PI;

// use std::collections::HashMap;
use hashbrown::HashMap;
use ahash::RandomState;




pub fn build_grid(alphas: &Vec<f64>, betas: &Vec<f64>, pixel_scale: f64) -> Vec<Vec<usize>> {
    let mut alpha_min = std::f64::INFINITY;
    let mut alpha_max = std::f64::NEG_INFINITY;
    let mut beta_min = std::f64::INFINITY;
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
    let mut beta_min = std::f64::INFINITY;
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

pub type MembersMap = HashMap<u64, Vec<u32>, RandomState>;

// Same binning as build_grid_sparse, but also records which entry indices
// (positions in the alphas/betas arrays) fall in each cell. This is the
// peak->detection membership the back-end needs to fit/emit a candidate.
pub fn build_grid_sparse_members(alphas: &[f64], betas: &[f64], pixel_scale: f64)
    -> (CountsMap, MembersMap, f64, f64, f64, usize, usize) {
    debug_assert_eq!(alphas.len(), betas.len());
    let n = alphas.len();
    if n == 0 {
        return (CountsMap::with_hasher(RandomState::new()),
                MembersMap::with_hasher(RandomState::new()),
                0.0, 0.0, pixel_scale, 0, 0);
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

    let mut counts: CountsMap = CountsMap::with_capacity_and_hasher(n * 2, RandomState::new());
    let mut members: MembersMap = MembersMap::with_capacity_and_hasher(n * 2, RandomState::new());

    for (idx, (&a, &b)) in alphas.iter().zip(betas.iter()).enumerate() {
        let mut ai = ((a - alpha_min) * inv) as isize;
        let mut bi = ((b - beta_min)  * inv) as isize;
        if ai < 0 { ai = 0; }
        if bi < 0 { bi = 0; }
        if ai as usize >= alpha_bins { ai = alpha_bins.saturating_sub(1) as isize; }
        if bi as usize >= beta_bins  { bi = beta_bins.saturating_sub(1)  as isize; }
        let key = pack(ai as u32, bi as u32);
        *counts.entry(key).or_insert(0) += 1;
        members.entry(key).or_insert_with(Vec::new).push(idx as u32);
    }

    (counts, members, alpha_min, beta_min, pixel_scale, alpha_bins, beta_bins)
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

// Return the exposures, the mean time of the detections, and the epoch-sorted
// detections. The exposures are contiguous slices of the sorted vec in the same
// order sync_phi emits, so the flat entry index (position in phi/alphas/betas)
// equals the index into this returned vec -- that is how the back-end recovers
// each detection's identity (detid/ast_ucty/mag/obs) for the fit and emit.
pub fn sort_detections_into_exposures(detections: &Vec<Detection>) -> (Vec<Exposure>, f64, Vec<Detection>) {
    // sort the detections by epoch
    let mut detections = detections.clone();
    detections.sort_by(|a, b| a.epoch.partial_cmp(&b.epoch).unwrap());

    // Calculate the mean epoch of all the detections to use as a reference epoch for the exposures.
    let mut mean_epoch = 0.0;
    for det in detections.iter() {
        mean_epoch += det.epoch;
    }
    mean_epoch /= detections.len() as f64;

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
        // Calculate the central rho vector for this exposure as the mean of the detections in the group.
        let mut central_rho = Vector3::zeros();
        for det in group.iter() {
            central_rho += det.rho_hat;
        }
        central_rho /= group.len() as f64;
        // Normalize the central rho vector so that it is a unit vector.
        central_rho = central_rho.normalize();
        let exposure = Exposure {
            id: format!("exposure_{}", idx),
            epoch: group[0].epoch,
            filter: None,
            central_rho,
            detections: group.iter().map(|d| d.rho_hat).collect(),
            observer_position: group[0].observer_position,
            observer_velocity: group[0].observer_velocity,
            reference_plane: group[0].reference_plane.clone(),
        };
        exposures.push(exposure);
    }
    (exposures, mean_epoch, detections)
}

pub fn transform_exposures_to_tangent_plane(exposures: &Vec<Exposure>, center: Vector3<f64>) -> Vec<TangentPlaneExposure> {
    // do stuff
    exposures.iter().map(|e| e.transform_to_tangent_plane(center)).collect()
}

// Compute the (gamma, gdot)-only part of the sync: per-detection phi_x, phi_y and
// the per-exposure coefficient c = g/f. These do NOT depend on (adot, bdot), so
// this is done ONCE per outer iteration; the inner (adot, bdot) loop then only
// shifts: alpha = phi_x - c*adot, beta = phi_y - c*bdot. Outputs are in the same
// entry order as the rest of the pipeline (exposures, then within-exposure).
pub fn sync_phi(exposures: &Vec<TangentPlaneExposure>, mu: f64, gamma: f64, gdot: f64, ref_epoch: f64, ndets: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {

    let mut phi_xs = Vec::with_capacity(ndets);
    let mut phi_ys = Vec::with_capacity(ndets);
    let mut cs     = Vec::with_capacity(ndets);   // c = g/f, per entry (per exposure)

    let z0 = 1.0/gamma;
    let z0dot = z0 * gdot;

    // Mean motion squared, used in the light time correction.
    let mm2 = mu * gamma * gamma * gamma;
    let mut dt = 0.0;
    for exposure in exposures.iter() {

        // Get the central theta_x0 and theta_y0 for this exposure.
        let theta_x0 = exposure.theta_x0;
        let theta_y0 = exposure.theta_y0;

        let t = exposure.epoch - ref_epoch;

        // Iterate to solve for the light time correction (only depends on gamma,gdot).
        // Re-use the last dt as the initial guess to speed up convergence.
        for _ in 0..2 {
            let tp = t - dt;
            let z = z0 + z0dot * tp - 0.5 * mm2 * tp * tp;
            let rho2 = (1.0 + theta_x0 * theta_x0 + theta_y0 * theta_y0)*(z - exposure.ze).powi(2);
            let rho = rho2.sqrt();
            dt = rho/SPEED_OF_LIGHT;
        }

        let tp = t - dt;
        let f = 1.0 - 0.5 * mm2 * tp * tp;
        let g = tp;
        let fac = 1.0 + g/f * gdot - gamma/f * exposure.ze;
        let c = g/f;

        for (theta_x, theta_y) in exposure.theta_x.iter().zip(exposure.theta_y.iter()) {
            phi_xs.push(theta_x * fac + gamma/f * exposure.xe);
            phi_ys.push(theta_y * fac + gamma/f * exposure.ye);
            cs.push(c);
        }
    }
    (phi_xs, phi_ys, cs)
}

fn linear_fit_with_residuals_and_max_residual(x: &[f64], y: &[f64]) -> Option<(f64, f64, Vec<f64>, f64)> {
    if x.len() != y.len() || x.len() < 2 {
        return None;
    }

    let n = x.len() as f64;
    let sum_x = x.iter().sum::<f64>();
    let sum_y = y.iter().sum::<f64>();
    let sum_xy = x.iter().zip(y).map(|(xi, yi)| xi * yi).sum::<f64>();
    let sum_x2 = x.iter().map(|xi| xi * xi).sum::<f64>();

    let denom = n * sum_x2 - sum_x * sum_x;
    if denom.abs() < 1e-12 {
        return None; // Prevent division by zero
    }

    let a = (n * sum_xy - sum_x * sum_y) / denom;
    let b = (sum_y * sum_x2 - sum_x * sum_xy) / denom;

    let residuals: Vec<f64> = x.iter().zip(y).map(|(xi, yi)| yi - (a * xi + b)).collect();
    let max_residual = residuals.iter().map(|r| r.abs()).fold(0./0., f64::max);

    Some((a, b, residuals, max_residual)) // slope, intercept, residuals, max residual
} 

pub fn analyze_tangent_exposures(tangent_exposures: &Vec<TangentPlaneExposure>, mean_epoch: f64) -> (f64, f64, f64, f64, f64, f64, f64) {

    // Collect the t, xe, ye, t^2, t*xe, t*ye values from the 
    // tangent exposures, to check for nonlinearity in the calculated
    // observed tangent plane locations.
    let mut xe_values = Vec::new();
    let mut ye_values = Vec::new();
    let mut t_values = Vec::new();
    let mut t2_values = Vec::new();
    let mut t_xe_values = Vec::new();
    let mut t_ye_values = Vec::new();
    for exposure in tangent_exposures.iter() {
        xe_values.push(exposure.xe);
        ye_values.push(exposure.ye);
        let t = exposure.epoch - mean_epoch;
        t_values.push(t);
        t2_values.push(t * t);
        t_xe_values.push(t * exposure.xe);
        t_ye_values.push(t * exposure.ye);
    }

    // Fit a line vs time to the xe and ye values. This is just for diagnostics and to check for nonlinearity.
    let (_xe_slope, _xe_intercept, _xe_residuals, max_xe_residual) = linear_fit_with_residuals_and_max_residual(&t_values, &xe_values).unwrap_or((0.0, 0.0, Vec::new(), 0.0));
    let (_ye_slope, _ye_intercept, _ye_residuals, max_ye_residual) = linear_fit_with_residuals_and_max_residual(&t_values, &ye_values).unwrap_or((0.0, 0.0, Vec::new(), 0.0));

    // Fit a line vs time to the t*xe and t*ye values. This is just for diagnostics and to check for nonlinearity.
    let (_t_xe_slope, _t_xe_intercept, _t_xe_residuals, max_t_xe_residual) = linear_fit_with_residuals_and_max_residual(&t_values, &t_xe_values).unwrap_or((0.0, 0.0, Vec::new(), 0.0));
    let (_t_ye_slope, _t_ye_intercept, _t_ye_residuals, max_t_ye_residual) = linear_fit_with_residuals_and_max_residual(&t_values, &t_ye_values).unwrap_or((0.0, 0.0, Vec::new(), 0.0));

    // Fit a line vs time to the t^2 values. This is just for diagnostics and to check for nonlinearity.
    let (_t2_slope, _t2_intercept, _t2_residuals, max_t2_residual) = linear_fit_with_residuals_and_max_residual(&t_values, &t2_values).unwrap_or((0.0, 0.0, Vec::new(), 0.0));

    // Calculate the min and max time values to check the time span of the exposures.
    let t_min = t_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let t_max = t_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    
    (t_min, t_max, max_xe_residual, max_t_xe_residual, max_ye_residual, max_t_ye_residual, max_t2_residual)
}

// TODO:
// 0. Add a calculation of the central time for all the exposures.
// 1. Exposures should have a central RA and Dec.
// 2. The initial conditions should be organized to be more efficient to loop over.
// 3. Work on generating the grid.
// 4. Think about the 1-d approach.
// 

// Fit one peak's member detections with the shared coupled αβγ tangent-plane fit
// (`iterative_reject`) and append matt-0-format rows to `out`. `entry_ids` are
// positions in the alphas/betas arrays (== indices into `dets`, the sorted
// detections); the e_* arrays hold the IC-independent per-entry fit inputs. The
// tangent-plane coords are used as-is (no re-projection). Returns true if a fit
// was emitted.
// --- SYNCED-coordinate fit (prototype, env TANGENT_SYNCED_FIT=1) --------------------------
// Gather sync identity: phi = alpha0 + c*rate per axis (c = g/f). Fitting the SYNCED coords --
// which already carry the 2nd-order two-body f,g + exact parallax -- matches the gather's order,
// so a true clump leaves only noise and the arc-end points are NOT clipped (the raw-theta
// 1st-order fit clips them). Two weighted linear regressions phi_x~c, phi_y~c + worst-|pull| reject.
fn wls_phi_c(keep: &[usize], c: &[f64], phi: &[f64], sig: &[f64]) -> Option<(f64, f64)> {
    let (mut sw, mut swc, mut swcc, mut swp, mut swcp) = (0.0f64, 0.0, 0.0, 0.0, 0.0);
    for &i in keep {
        let w = 1.0 / (sig[i] * sig[i]);
        sw += w; swc += w * c[i]; swcc += w * c[i] * c[i];
        swp += w * phi[i]; swcp += w * c[i] * phi[i];
    }
    let d = sw * swcc - swc * swc;
    if d.abs() < 1e-300 { return None; }
    Some(((swcc * swp - swc * swcp) / d, (sw * swcp - swc * swp) / d))
}

fn synced_reject(
    t: &[f64], c: &[f64], phi_x: &[f64], phi_y: &[f64], sig: &[f64],
    sigma_threshold: f64, min_unique_times: usize, min_nights: usize,
) -> Option<([f64; 4], f64, f64, Vec<usize>, Vec<usize>)> {
    let mut keep: Vec<usize> = (0..t.len()).collect();
    let mut rejected: Vec<usize> = Vec::new();
    loop {
        let mut tn: Vec<i64> = keep.iter().map(|&i| t[i].round() as i64).collect();
        tn.sort_unstable(); tn.dedup();
        if tn.len() < min_nights { return None; }
        let (a, adot) = wls_phi_c(&keep, c, phi_x, sig)?;
        let (b0, bdot) = wls_phi_c(&keep, c, phi_y, sig)?;
        let mut worst = 0.0f64;
        let mut worst_i: Option<usize> = None;
        let mut chi2 = 0.0f64;
        for &i in &keep {
            let px = (phi_x[i] - (a + adot * c[i])) / sig[i];
            let py = (phi_y[i] - (b0 + bdot * c[i])) / sig[i];
            chi2 += px * px + py * py;
            let pm = px.abs().max(py.abs());
            if pm > worst { worst = pm; worst_i = Some(i); }
        }
        if worst < sigma_threshold {
            let dof = (2 * keep.len()).saturating_sub(4) as f64;
            return Some(([a, adot, b0, bdot], chi2, dof, keep.clone(), rejected.clone()));
        }
        match worst_i {
            Some(wi) => { rejected.push(wi); keep.retain(|&x| x != wi); }
            None => return None,
        }
        if keep.len() < min_unique_times { return None; }
    }
}

fn fit_and_emit_peak(
    out: &mut String,
    entry_ids: &[u32],
    e_t: &[f64], e_tx: &[f64], e_ty: &[f64],
    e_mxe: &[f64], e_mye: &[f64], e_ast: &[f64],
    phi_x_g: &[f64], phi_y_g: &[f64], cc_g: &[f64],
    dets: &Vec<Detection>,
    config: &Config,
    obliq_rad: f64,
    ic_id: &str,
) -> bool {
    let m = entry_ids.len();
    if m < config.min_points_per_cell {
        return false;
    }
    let t:   Vec<f64> = entry_ids.iter().map(|&i| e_t[i as usize]).collect();
    let tx:  Vec<f64> = entry_ids.iter().map(|&i| e_tx[i as usize]).collect();
    let ty:  Vec<f64> = entry_ids.iter().map(|&i| e_ty[i as usize]).collect();
    let mxe: Vec<f64> = entry_ids.iter().map(|&i| e_mxe[i as usize]).collect();
    let mye: Vec<f64> = entry_ids.iter().map(|&i| e_mye[i as usize]).collect();
    let ast: Vec<f64> = entry_ids.iter().map(|&i| e_ast[i as usize]).collect();

    // -- SYNCED-fit path (prototype, env TANGENT_SYNCED_FIT=1): fit the gather's synced phi
    //    (alpha0 + c*rate) instead of raw theta, matching the gather order so arc-end points
    //    are not clipped. Same TANGENT_REJECT_DUMP diagnostic as the raw path. --
    if std::env::var("TANGENT_SYNCED_FIT").is_ok() {
        let c_g:  Vec<f64> = entry_ids.iter().map(|&i| cc_g[i as usize]).collect();
        let px_g: Vec<f64> = entry_ids.iter().map(|&i| phi_x_g[i as usize]).collect();
        let py_g: Vec<f64> = entry_ids.iter().map(|&i| phi_y_g[i as usize]).collect();
        let sr = synced_reject(&t, &c_g, &px_g, &py_g, &ast,
            config.sigma_threshold, config.min_unique_times, config.min_nights);
        if let Ok(dump_path) = std::env::var("TANGENT_REJECT_DUMP") {
            let rej_min: usize = std::env::var("TANGENT_REJECT_MIN").ok()
                .and_then(|s| s.parse().ok()).unwrap_or(6);
            if m >= rej_min {
                use std::io::Write;
                let gathered: Vec<String> = (0..m)
                    .map(|i| dets[entry_ids[i] as usize].detid.clone().unwrap_or_default())
                    .collect();
                let (n_rej, rej_ids): (String, String) = match &sr {
                    Some((_, _, _, _, rejected)) => {
                        let r: Vec<String> = rejected.iter()
                            .map(|&i| dets[entry_ids[i] as usize].detid.clone().unwrap_or_default())
                            .collect();
                        (r.len().to_string(), r.join(","))
                    }
                    None => ("ERR".to_string(), String::new()),
                };
                if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&dump_path) {
                    let _ = writeln!(f, "{}\tn_gathered={}\tn_rejected={}\tgathered={}\trejected={}",
                        ic_id, m, n_rej, gathered.join(","), rej_ids);
                }
            }
        }
        let (params, chi2, dof, kept, _rejected) = match sr {
            Some(v) => v,
            None => return false,
        };
        out.push_str(&format!("# Candidate (IC {}): params {:?} chi2 {:.3} dof {:.0}\n",
            ic_id, params, chi2, dof));
        out.push_str("# detid             epoch          t(days)      c          pull_x    pull_y   (synced)\n");
        for &i in &kept {
            let det = &dets[entry_ids[i] as usize];
            let detid = det.detid.clone().unwrap_or_default();
            let px = (px_g[i] - (params[0] + params[1] * c_g[i])) / ast[i];
            let py = (py_g[i] - (params[2] + params[3] * c_g[i])) / ast[i];
            out.push_str(&format!("{:14} {:13.7} {:12.6} {:11.4} {:8.3} {:8.3}\n",
                detid, det.epoch, t[i], c_g[i], px, py));
        }
        out.push_str(&format!("# Kept {} detections in synced fit (IC {})\n", kept.len(), ic_id));
        return true;
    }

    // x = a + b*t + c*mxe ; y = d + e*t + c*mye  (mxe=-xe, mye=-ye), sigma-clipped.
    let reject_result = iterative_reject(
        &t, &mxe, &mye, &tx, &ty, &ast, &ast,
        config.sigma_threshold,
        config.min_points_per_cell,
        config.min_unique_times,
        config.min_nights,
    );
    // -- diagnostic (env TANGENT_REJECT_DUMP=path, gate TANGENT_REJECT_MIN default 6):
    //    per emitted candidate, the PRE-reject gathered detids + which iterative_reject
    //    rejected. Separates "never gathered" from "gathered-then-clipped" for the 8/10
    //    recovery cap. Inert unless TANGENT_REJECT_DUMP is set. --
    if let Ok(dump_path) = std::env::var("TANGENT_REJECT_DUMP") {
        let rej_min: usize = std::env::var("TANGENT_REJECT_MIN").ok()
            .and_then(|s| s.parse().ok()).unwrap_or(6);
        if m >= rej_min {
            use std::io::Write;
            let gathered: Vec<String> = (0..m)
                .map(|i| dets[entry_ids[i] as usize].detid.clone().unwrap_or_default())
                .collect();
            let (n_rej, rej_ids): (String, String) = match &reject_result {
                Ok((_, _, rejected)) => {
                    let r: Vec<String> = (0..m).filter(|&i| rejected.contains(&i))
                        .map(|i| dets[entry_ids[i] as usize].detid.clone().unwrap_or_default())
                        .collect();
                    (r.len().to_string(), r.join(","))
                }
                Err(_) => ("ERR".to_string(), String::new()),
            };
            if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&dump_path) {
                let _ = writeln!(f, "{}\tn_gathered={}\tn_rejected={}\tgathered={}\trejected={}",
                    ic_id, m, n_rej, gathered.join(","), rej_ids);
            }
        }
    }
    let (fit, _kept, rejected) = match reject_result {
        Ok(v) => v,
        Err(_) => return false,
    };
    let p = fit.params;
    out.push_str(&format!("# Candidate (IC {}): params {:?} chi2 {:.3} dof {:.0}\n",
        ic_id, p, fit.chi2, fit.dof));
    out.push_str("# detid             epoch          RA(deg)      sig(\")   Dec(deg)     obs_x          obs_y          obs_z      obscode  mag mag_sig  filt   intid     t(days)   theta_x(\")  theta_y(\")  x_model(\")  y_model(\")  pull_x    pull_y\n");
    let mut kept_count = 0;
    for i in 0..m {
        if rejected.contains(&i) {
            continue;
        }
        kept_count += 1;
        let x_model = p[0] + p[1] * t[i] + p[2] * mxe[i];
        let y_model = p[3] + p[4] * t[i] + p[2] * mye[i];
        let pull_x = (tx[i] - x_model) / ast[i];
        let pull_y = (ty[i] - y_model) / ast[i];

        let det = &dets[entry_ids[i] as usize];
        let rho = det.rho_hat;                       // ecliptic unit vector
        let lam = rho[1].atan2(rho[0]);
        let bet = rho[2].asin();
        let (ra, dec) = ecliptic_to_equatorial(lam, bet, obliq_rad);
        let obs2 = ecliptic_to_equatorial_vector(det.observer_position, obliq_rad);
        let detid = det.detid.clone().unwrap_or_default();
        let filt = det.filter.clone().unwrap_or_default();
        let obscode = det.obscode.clone().unwrap_or_default();
        let mag = det.mag.unwrap_or(0.0);
        let mag_sig = det.mag_ucty.unwrap_or(0.0);

        out.push_str(&format!(
            "{:14} {:13.7} {:11.7}   {:6.3}   {:11.7}  {:13.10}  {:13.10}  {:13.10}   {}   {:.3} {:.3} {:6}  {:8} {:12.6}  {:9.4}  {:9.4}  {:9.4}  {:9.4}  {:8.3}  {:8.3}\n",
            detid, det.epoch, ra * 180.0 / PI, ast[i] * 206265.0, dec * 180.0 / PI,
            obs2[0], obs2[1], obs2[2], obscode, mag, mag_sig, filt,
            entry_ids[i], t[i], tx[i] * 206265.0, ty[i] * 206265.0,
            x_model * 206265.0, y_model * 206265.0, pull_x, pull_y));
    }
    out.push_str(&format!("# Kept {} detections in fit (IC {})\n", kept_count, ic_id));
    true
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

    let detections = load_detections(&config.detection_catalog, &config.orbit_reference_plane, &kernel)?;
    let (exposures, mean_epoch, detections_sorted) = sort_detections_into_exposures(&detections);

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

    let (t_min, t_max,
        _max_xe_residual, max_t_xe_residual,
        _max_ye_residual, max_t_ye_residual,
        max_t2_residual) = analyze_tangent_exposures(&tangent_exposures, mean_epoch);

    let epsilon = config.epsilon_arcsec;
    let pixel_scale = epsilon / 3.0 * std::f64::consts::PI / (180.0 * 3600.0); // convert arcsec to radians

    let epsilon_radians = epsilon * std::f64::consts::PI / (180.0 * 3600.0); // convert arcsec to radians
    let max_t_xeye_residual = max_t_xe_residual.max(max_t_ye_residual);

    let dgamma = epsilon_radians / max_t_xeye_residual;
    let gamma = 1.0/35.0;
    let mm2_max = 2.0 * mu * gamma * gamma * gamma;
    let angular_residual = max_t2_residual * mm2_max;
    let dgamma_v2 = epsilon_radians / angular_residual;


    //println!("Max residuals in observer position: xe: {:.3} au, ye: {:.3} au", max_xe_residual, max_ye_residual);
    println!("Gamma step size corresponding to max t*xe/ye residual: {:.3e} au^-1", dgamma);
    println!("Gamma step size corresponding to max t^2 residuals: {:.3e} au^-1", dgamma_v2);
    println!("Time span of exposures: t_min: {:.3} days, t_max: {:.3} days", t_min, t_max);

    
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
    
    
    let ref_epoch = mean_epoch;
    let time_span = t_max - t_min;
    println!("Reference epoch: {:.3}, time span: {:.3} days", ref_epoch, time_span);

    // Precompute the IC-independent per-entry fit inputs, in the SAME order
    // sync_phi emits phi (exposures, then within-exposure), so the entry index
    // aligns with phi/alphas/betas and detections_sorted. t/observer
    // terms are per-exposure; theta_x/theta_y are the observed tangent coords used
    // as-is; ast_ucty is per detection (radians, from load_detections).
    let default_ast = 0.2 * (1.0 / 3600.0) * PI / 180.0;
    let mut e_t:   Vec<f64> = Vec::with_capacity(count);
    let mut e_tx:  Vec<f64> = Vec::with_capacity(count);
    let mut e_ty:  Vec<f64> = Vec::with_capacity(count);
    let mut e_mxe: Vec<f64> = Vec::with_capacity(count);
    let mut e_mye: Vec<f64> = Vec::with_capacity(count);
    let mut e_ast: Vec<f64> = Vec::with_capacity(count);
    let mut gidx = 0usize;
    for te in tangent_exposures.iter() {
        let t = te.epoch - ref_epoch;
        let mxe = -te.xe;
        let mye = -te.ye;
        for k in 0..te.theta_x.len() {
            e_t.push(t);
            e_tx.push(te.theta_x[k]);
            e_ty.push(te.theta_y[k]);
            e_mxe.push(mxe);
            e_mye.push(mye);
            e_ast.push(detections_sorted[gidx].ast_ucty.unwrap_or(default_ast));
            gidx += 1;
        }
    }

    let obliq_rad = (84381.448 / 3600.0) * PI / 180.0;
    let mut outfile = std::fs::File::create(&config.output_path)?;
    writeln!(outfile, "# tangent_v2 search output")?;
    writeln!(outfile, "# detection_catalog: {}", config.detection_catalog)?;
    writeln!(outfile, "# reference_epoch: {:.6}  time_span_days: {:.3}", ref_epoch, time_span)?;
    writeln!(outfile, "# epsilon_arcsec: {}  sigma_threshold: {}  min_points_per_cell: {}  min_unique_times: {}  min_nights: {}",
        config.epsilon_arcsec, config.sigma_threshold, config.min_points_per_cell, config.min_unique_times, config.min_nights)?;

    // --- DIAGNOSTIC (env-gated): dump per-detection synced phi at each (gamma,gdot).
    // Set TANGENT_PHI_DUMP=/path to enable. With a 1-IC grid at the true (gamma,gdot)
    // this exposes the irreducible synced residual at truth: fit (adot,bdot) in post
    // (alpha = phi_x - c*adot is linear in adot) and plot residual vs time. Sub-arcsec
    // residual at truth => model is fine, 7/10 cap is grid under-sampling; |t|-growing
    // residual at truth => a model/convention bug (light-time/aberration/observer).
    let phi_dump = std::env::var("TANGENT_PHI_DUMP").ok();
    if let Some(ref p) = phi_dump {
        let mut df = std::fs::File::create(p)?;
        writeln!(df, "detid,t,epoch,phi_x,phi_y,c,gamma,gdot")?;
    }

    let start_time = std::time::Instant::now();

    // Now we can loop over our grid of adot and bdot values,
    // sync the detections to the reference epoch for each pair, 
    // build a grid in alpha-beta space, and count the number of detections that fall into each cell.
    
    // Group ICs by (gamma,gdot) == (1/r, vr/r), i.e. by (r,vr), so sync_phi runs
    // ONCE per group and its phi/c are reused across the group's (adot,bdot).
    initial_conditions.sort_by(|a, b|
        a.0.partial_cmp(&b.0).unwrap().then(a.1.partial_cmp(&b.1).unwrap()));
    let mut groups: Vec<(f64, f64, Vec<(f64, f64, String)>)> = Vec::new();
    for (r, vr, adot, bdot, id) in initial_conditions.into_iter() {
        match groups.last_mut() {
            Some(g) if g.0 == r && g.1 == vr => g.2.push((adot, bdot, id)),
            _ => groups.push((r, vr, vec![(adot, bdot, id)])),
        }
    }
    let n_ics: usize = groups.iter().map(|g| g.2.len()).sum();
    println!("Grouped {} ICs into {} (gamma,gdot) groups", n_ics, groups.len());

    let min_count = 3;
    // Per-(gamma,gdot) group work: independent -> parallelizable. Each group runs sync_phi
    // once, sweeps its (adot,bdot) motions, and returns its output text + candidate count.
    // Per-motion (adot,bdot) work: independent within a group (reads the group's synced
    // phi/cc read-only, emits its own output). Parallelize the ~1370 motions/group across all
    // cores -> full utilization, no per-group load-imbalance tail (the outer-group version's
    // problem). Returns (output text, candidate count) for one motion.
    let process_motion = |adot: f64, bdot: f64, id: &String,
                          phi_x: &[f64], phi_y: &[f64], cc: &[f64],
                          gamma: f64, gdot: f64| -> (String, usize) {
        let mut alphas = Vec::with_capacity(count);
        let mut betas  = Vec::with_capacity(count);
        for i in 0..phi_x.len() {
            alphas.push(phi_x[i] - cc[i] * adot);
            betas.push(phi_y[i] - cc[i] * bdot);
        }
        let (counts, cell_members, _ax, _bx, _ps, alpha_bins, beta_bins) =
            build_grid_sparse_members(&alphas, &betas, pixel_scale);
        let peaks = find_local_maxima_sparse(&counts, alpha_bins as u32, beta_bins as u32, min_count);
        let convolved_counts = convolve_peaks_sparse(&counts, &peaks, alpha_bins as u32, beta_bins as u32);

        let mut out = String::new();
        let mut mc = 0usize;
        let mut seen: HashSet<Vec<u32>> = HashSet::new();
        for (peak, conv) in peaks.iter().zip(convolved_counts.iter()) {
            if (*conv as usize) < config.min_points_per_cell {
                continue;
            }
            let (pi, pj) = (peak.0, peak.1);
            let mut ids: Vec<u32> = Vec::new();
            for di in -1i64..=1 {
                for dj in -1i64..=1 {
                    let ni = pi as i64 + di;
                    let nj = pj as i64 + dj;
                    if ni < 0 || nj < 0 { continue; }
                    if let Some(v) = cell_members.get(&pack(ni as u32, nj as u32)) {
                        ids.extend_from_slice(v);
                    }
                }
            }
            ids.sort_unstable();
            if !seen.insert(ids.clone()) {
                continue;
            }
            let ic_id = format!("{} (g={:.4},gdot={:.2e},adot={:.3e},bdot={:.3e})",
                                id, gamma, gdot, adot, bdot);
            if fit_and_emit_peak(&mut out, &ids, &e_t, &e_tx, &e_ty, &e_mxe, &e_mye,
                                 &e_ast, phi_x, phi_y, cc, &detections_sorted, &config, obliq_rad, &ic_id) {
                mc += 1;
            }
        }
        (out, mc)
    };

    // Outer (gamma,gdot) loop stays SERIAL (one sync_phi each + ordered file writes); inner
    // motions run in parallel (or serial when a diagnostic dump appends to shared files).
    let dumps_on = phi_dump.is_some() || std::env::var("TANGENT_REJECT_DUMP").is_ok();
    let mut n_candidates = 0usize;
    for (r, vr, motions) in groups.iter() {
        let gamma = 1.0 / r;
        let gdot = vr / r;
        let (phi_x, phi_y, cc) = sync_phi(&tangent_exposures, mu, gamma, gdot, ref_epoch, count);

        if let Some(ref p) = phi_dump {
            if let Ok(mut df) = std::fs::OpenOptions::new().append(true).open(p) {
                for i in 0..phi_x.len() {
                    let detid = detections_sorted[i].detid.clone().unwrap_or_default();
                    let _ = writeln!(df, "{},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.9e},{:.9e}",
                        detid, e_t[i], detections_sorted[i].epoch,
                        phi_x[i], phi_y[i], cc[i], gamma, gdot);
                }
            }
        }

        let results: Vec<(String, usize)> = if dumps_on {
            motions.iter()
                .map(|(adot, bdot, id)| process_motion(*adot, *bdot, id, &phi_x, &phi_y, &cc, gamma, gdot))
                .collect()
        } else {
            motions.par_iter()
                .map(|(adot, bdot, id)| process_motion(*adot, *bdot, id, &phi_x, &phi_y, &cc, gamma, gdot))
                .collect()
        };
        for (s, c) in &results {
            n_candidates += *c;
            if !s.is_empty() {
                write!(outfile, "{}", s)?;
            }
        }
    }
    outfile.flush()?;

    let duration = start_time.elapsed();
    println!("completed in {:?}: {} candidates -> {}", duration, n_candidates, config.output_path);

    
    
   
    Ok(())

    

}