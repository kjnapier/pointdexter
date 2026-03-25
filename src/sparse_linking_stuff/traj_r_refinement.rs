use rayon::prelude::*;
use std::fmt;
use crate::sparse_linking_stuff::orbit::{Orbit, find_max_separation_cheat};
use crate::sparse_linking_stuff::InitialCondition;
use crate::sparse_linking_stuff::InputFormat;

use ordered_float::OrderedFloat;
use std::collections::HashSet;
use spacerocks::Time;
use std::sync::Arc;
use std::cmp::Ordering;
use std::sync::Mutex;
use std::collections::HashMap;
use crate::sparse_linking_stuff::Detection;
use crate::sparse_linking_stuff::trajectory::{Trajectory, angular_separation};
use spacerocks::SpaceRock;
use spacerocks::constants::{SPEED_OF_LIGHT};
use spacerocks::transforms::{solve_for_universal_anomaly, stumpff_c, stumpff_s};
use nalgebra::Vector3;
use crate::sparse_linking_stuff::grid::*;
use std::f64::consts::PI;

const ARCSEC_PER_RAD: f64 = 3600.0 * 180.0 / std::f64::consts::PI;
const RADS_TO_ARCSEC: f64 = 206265.0;

#[derive(Debug, Clone)]
pub struct MetricEntry {
    pub phi_vals: Vec<f64>,
    pub theta_vals: Vec<f64>,
    pub epochs: Vec<f64>,
}

fn cost_function(det: &Detection, r: f64, rho: f64) -> f64 {                                                                                  
    r.powi(2)
        - rho.powi(2)
        - 2.0 * rho * det.pointing.dot(&det.observer.position)
        - det.observer.position.dot(&det.observer.position)
}

fn gradient(det: &Detection, rho: f64) -> f64 {
    -2.0 * (rho + det.pointing.dot(&det.observer.position))
}

fn optimize_rho(det: &Detection, r: f64, mut rho: f64) -> f64 {
    const TOL: f64 = 1.0e-10;
    const MAX_ITER: usize = 50;

    let mut f = cost_function(det, r, rho);
    for _ in 0..MAX_ITER {
        if f.abs() < TOL {
            break;
        }
        let g = gradient(det, rho);
        if g.abs() < 1e-16 {
            break; // avoid div-by-zero
        }
        rho -= f / g;
        f = cost_function(det, r, rho);
    }
    rho
}

pub fn construct_orbit_from_r(detection: &Detection, ic: &InitialCondition, r: f64) -> Option<([f64; 3])> {

    // Solve for the distance from the observer to the object numerically. More accurate.
    let rho = optimize_rho(detection, r, ic.r0 - 1.0);

    // Get the actual epoch of the detection after light travel time correction.
    // This is the epoch at which the object was actually at the position observed by the observer
    let light_corrected_epoch = detection.epoch.clone() - rho / SPEED_OF_LIGHT;
    let dt = ic.epoch.epoch - light_corrected_epoch.utc().jd();

    // Calculate the position of the object at the time of detection in barycentric coordinates and its distance 
    let r_vec = detection.observer.position + rho * detection.pointing();
    let r = r_vec.norm();

    // Get the angular coordinates of the object in spherical coordinates
    let phi = r_vec.y.atan2(r_vec.x);
    let theta = (r_vec.z / r).asin();

    // Use angular coordinates to set up Kev's basis reference frame 
    let ahat = Vector3::new(-phi.sin(), phi.cos(), 0.0);
    let dhat = Vector3::new(-theta.sin() * phi.cos(), -theta.sin() * phi.sin(), theta.cos());

    // Assert that psi = inclination. This is an assertion that will get worse as the object gets to larger latitudes.
    let mut inclination = ic.psi.abs();
    let kappa = ic.psi.signum();

    let theta_abs = theta.abs();

    if inclination < theta_abs {
        inclination = match (ic.p4_min.map(|v| v.abs()), ic.p4_max.map(|v| v.abs())) {
            (Some(min), Some(max)) => {
                if min >= theta_abs && max >= theta_abs {
                    min.min(max)
                } else if min >= theta_abs {
                    min
                } else if max >= theta_abs {
                    max
                } else {
                    // println!("Early exit: no valid p4_min/max ≥ |theta| = {}, min = {}, max = {}: Det orbit id = {}, ic_id = {}", theta_abs, min, max, detection.orbit_id.unwrap(), ic.id);
                    return None;
                }
            }
            (Some(min), None) if min >= theta_abs => min,
            (None, Some(max)) if max >= theta_abs => max,
            _ => {
                println!("Early exit: no valid p4 bounds available");
                return None;
            }
        };
    }


    // Get the velocity components of the object
    let vo = ic.h / r;
    let vsq = 2.0 * (ic.energy + ic.mu / r);
    let vrsq = vsq - vo.powi(2);
    let mut vr = if vrsq < 0.0 { 0.0 } else { vrsq.sqrt() };
    let mut mean_anomaly = (ic.mean_anomaly + ic.n() * dt) % (2.0 * PI);
    if mean_anomaly > PI {
        vr = -vr;
    }



    // Get the  full velocity vector
    let cos_psi = inclination.cos() / theta.cos();
    let mut sin_psi = kappa * (1.0 - cos_psi.powi(2)).sqrt();
    let mut v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);

    
    // Using r and v, propagate the orbit to the time of the initial condition w/ the Lagrange coefficients.
    let alpha = -2.0 * ic.energy / ic.mu;
    
    let chi = match solve_for_universal_anomaly(r, vr, alpha, ic.mu, dt, 1e-10, 100) {
        Ok(val) => val,
        Err(e) => {
            println!("Error solving for universal anomaly: {}", e);
            return None
        }
    };


    let stumpff_s_val = stumpff_s(alpha * chi.powi(2));
    let stumpff_c_val = stumpff_c(alpha * chi.powi(2));

    let f = 1.0 - chi.powi(2) * stumpff_c_val / r;
    let g = dt - chi.powi(3) * stumpff_s_val / ic.mu.sqrt();
    let mut new_position = f * r_vec + g * v_vec;

    
    let new_pointing = new_position.normalize();

    Some([new_pointing.x, new_pointing.y, new_pointing.z])

}



fn linspace(start: &f64, end: &f64, num: usize) -> Vec<f64> {
    if num == 0 {
        return Vec::new();
    }
    if num == 1 {
        return vec![*start];
    }
    let step = if num > 1 {
        (end - start) / ((num - 1) as f64)
    } else {
        0.0
    };
    (0..num).map(|i| start + (i as f64) * step).collect()
}

fn pointing_to_phi_theta(x: f64, y: f64, z: f64) -> Option<(f64, f64)> {
    let phi = y.atan2(x);
    let theta = z.asin();
    Some((phi, theta))
}

pub fn mean(xs: &[f64]) -> f64 {
    if xs.is_empty() { return f64::NAN; }
    xs.iter().sum::<f64>() / xs.len() as f64
}
pub fn variance(xs: &[f64]) -> f64 {
    if xs.len() < 2 { return 0.0; }
    let m = mean(xs);
    xs.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (xs.len() as f64)
}

/// R² of a simple linear regression y ~ a*t + b
pub fn r2_linear(y: &[f64], t: &[f64]) -> f64 {
    let n = y.len();
    if n < 2 || t.len() != n {
        return 0.0;
    }
    // standardize t (helps numerics with JD-scale times)
    let tm = mean(t);
    let tsd = (t.iter().map(|v| (v - tm).powi(2)).sum::<f64>() / n as f64).sqrt();
    let t_scaled: Vec<f64> = if tsd == 0.0 { t.iter().map(|_| 0.0).collect() }
                             else { t.iter().map(|v| (v - tm) / tsd).collect() };

    let y_mean = mean(y);
    let sst = y.iter().map(|v| (v - y_mean).powi(2)).sum::<f64>();
    if sst == 0.0 { return 0.0; }

    // closed-form slope/intercept on standardized t
    let sx = t_scaled.iter().sum::<f64>();
    let sy = y.iter().sum::<f64>();
    let sxx = t_scaled.iter().map(|x| x * x).sum::<f64>();
    let sxy = t_scaled.iter().zip(y.iter()).map(|(x, yy)| x * yy).sum::<f64>();
    let n_f = n as f64;
    let denom = n_f * sxx - sx * sx + 1e-12;
    let a = (n_f * sxy - sx * sy) / denom;
    let b = (sy - a * sx) / n_f;

    let ssr = t_scaled.iter()
        .zip(y.iter())
        .map(|(x, yy)| {
            let yhat = a * *x + b;
            (yy - yhat).powi(2)
        })
        .sum::<f64>();

    (1.0 - ssr / sst).clamp(0.0, 1.0)
}

fn norm01_mut(v: &mut [f64]) {
    if v.is_empty() { return; }
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for &x in v.iter() {
        if x.is_finite() {
            if x < lo { lo = x; }
            if x > hi { hi = x; }
        }
    }
    if !lo.is_finite() || !hi.is_finite() || hi == lo {
        // collapse to zeros (no influence)
        for x in v.iter_mut() { *x = 0.0; }
        return;
    }
    for x in v.iter_mut() {
        *x = (*x - lo) / (hi - lo);
    }
}



/// Recenters the r-bounds around `best_r` and tightens them until the
/// separation along ONLY the r (p1) axis is < `epsilon`.
///
/// Returns `Some(new_ic)` when successful, or `None` if we couldn't get
/// below `epsilon` before hitting `min_half_width` or `max_iter`.
fn refactor_home_cell(
    ic: &InitialCondition,
    best_r: f64,
    epsilon: f64,
    t_bounds: (f64, f64),
) -> Option<InitialCondition> {
    // NOTE: axis 0 == p1; for KEV this is r.
    // If using QEF/KEP, this still shrinks p1 (not "r").
    let t = ic.epoch.epoch;
    let mu = ic.mu;
    let format = ic.format.clone();

    // Start with the current half-width in r
    let rmin = ic.p1_min?;
    let rmax = ic.p1_max?;
    let mut half_width = 0.5 * (rmax - rmin);

    // Shrink controls
    let shrink_factor = 0.9;     // tighten by 10% each try (
    let min_half_width = 1.0e-6;  // AU (or whatever units p1 uses); guard against over-shrinking
    let max_iter = 64;

    for _ in 0..max_iter {
        let new_rmin = best_r - half_width;
        let new_rmax = best_r + half_width;

        let mut new_ic = ic.clone();
        new_ic.p1_min = Some(new_rmin);
        new_ic.p1_max = Some(new_rmax);
        new_ic.p1 = best_r;

        // Build a cell and measure separation along ONLY the r-axis
        let cell = match Cell::from_initial_condition(&new_ic) {
            Some(c) => c,
            None => return None,
        };

        let sep_r_arcsec = calculate_axis_separation(
            &cell,
            0,          // axis 0 == p1 (r for KEV)
            t,
            t_bounds,
            mu,
            format.clone(),
        );

        if sep_r_arcsec < epsilon {
            return Some(new_ic);
        }

        // Tighten and try again
        half_width *= shrink_factor;
        if half_width < min_half_width {
            break;
        }
    }

    None
}




pub fn optimize_r_and_refine_ic(
    ic: &InitialCondition,
    trajectory: &Trajectory,
    det_map: &HashMap<u64, Detection>,
    num_r: usize,
    w_spread: f64,
    w_linearity: f64,
    epsilon: f64,
    t_bounds: (f64, f64),
) -> Option<(Trajectory, InitialCondition)> {
    // ----- 1) Gather and sort detections -----
    let mut dets: Vec<&Detection> = trajectory
        .detection_ids
        .iter()
        .filter_map(|id| det_map.get(&(*id as u64)))
        .collect();

    dets.sort_by(|a, b| {
        a.epoch.epoch
            .partial_cmp(&b.epoch.epoch)
            .unwrap_or(Ordering::Equal)
    });

    // Precompute epochs to align with phi/theta lengths
    let epochs: Vec<f64> = dets.iter().map(|d| d.epoch.epoch).collect();

    // ----- 2) Grid search over r to pick best_r by your metric -----
    let rmin = ic.p1_min?;
    let rmax = ic.p1_max?;
    let r_vals = linspace(&rmin, &rmax, num_r);

    let mut spreads: Vec<f64> = Vec::with_capacity(r_vals.len());
    let mut linearities: Vec<f64> = Vec::with_capacity(r_vals.len());

    for &trial_r in &r_vals {
        let mut phi: Vec<f64> = Vec::with_capacity(dets.len());
        let mut theta: Vec<f64> = Vec::with_capacity(dets.len());

        for det in dets.iter() {
            if let Some([x, y, z]) = construct_orbit_from_r(det, ic, trial_r) {
                if let Some((ph, th)) = pointing_to_phi_theta(x, y, z) {
                    phi.push(ph);
                    theta.push(th);
                }
            }
        }

        if phi.len() < 3 {
            spreads.push(f64::INFINITY);
            linearities.push(0.0);
            continue;
        }

        let sp = variance(&phi) + variance(&theta);
        let r2_phi = r2_linear(&phi, &epochs[..phi.len()]);
        let r2_theta = r2_linear(&theta, &epochs[..theta.len()]);
        let lin = 0.5 * (r2_phi + r2_theta);

        spreads.push(sp);
        linearities.push(lin);
    }

    // Normalize & combine
    let mut spreads_n = spreads.clone();
    norm01_mut(&mut spreads_n);

    let mut lin_n = linearities.clone();
    norm01_mut(&mut lin_n);
    let inv_lin: Vec<f64> = lin_n.iter().map(|v| 1.0 - *v).collect();

    let (mut best_idx, mut best_score) = (0usize, f64::INFINITY);
    for i in 0..r_vals.len() {
        let s = w_spread * spreads_n[i] + w_linearity * inv_lin[i];
        if s < best_score {
            best_score = s;
            best_idx = i;
        }
    }

    let best_r = r_vals[best_idx];

    // ----- 3) Recenter & tighten r bounds until axis-0 separation < epsilon -----
    // Uses the helper we discussed earlier.
    let new_ic = refactor_home_cell(ic, best_r, epsilon, t_bounds)?;

    Some((trajectory.clone(), new_ic))
}