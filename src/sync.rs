use nalgebra::Vector3;
use crate::detection::Detection;
use crate::initial_condition::InitialCondition;

use spacerocks::data::{SPEED_OF_LIGHT};

// --- Orbit Sync ---
#[inline(never)]
fn cost_and_gradient(det: &Detection, ic: &InitialCondition, rho: f64) -> (f64, f64, f64) {
    let ltt = rho / SPEED_OF_LIGHT;
    let ep = det.epoch - ltt;
    let r = ic.r_at_epoch(ep);

    let cost = r * r
        - rho * rho
        - 2.0 * rho * det.rho_hat_dot_observer_position
        - det.observer_distance_squared;
    let grad = -2.0 * (rho + det.rho_hat_dot_observer_position);

    (cost, grad, r)
}

#[inline(never)]
fn optimize_rho(det: &Detection, ic: &InitialCondition, rho0: f64) -> (f64, f64) {
    // tunable parameter
    let tol = 1e-8;
    let mut rho = rho0;
    let (mut cost, mut grad, mut r) = cost_and_gradient(det, ic, rho);
    while cost.abs() > tol {
        if grad == 0.0 {
            break;
        }
        rho -= cost / grad;
        let (c, g, new_r) = cost_and_gradient(det, ic, rho);
        cost = c;
        grad = g;
        r = new_r;
    }
    (rho, r)
}

#[inline(never)]
pub fn sync_detection_to_orbit(det: &Detection, ic: &InitialCondition) -> Option<[f64; 3]> {
    let (rho, r) = optimize_rho(det, ic, ic.r - 1.0);

    let r_vec = det.observer_position + rho * det.rho_hat;

    let light_corrected_epoch = det.epoch - rho / SPEED_OF_LIGHT;
    if ((r_vec[2].abs() / r).asin()).abs() > ic.inc {
        return None;
    }

    let (x, y, z) = (r_vec[0], r_vec[1], r_vec[2]);
    let xy_norm = (x * x + y * y).sqrt();

    let ahat = Vector3::new(-y / xy_norm, x / xy_norm, 0.0);
    let dhat = Vector3::new(-z * x / (r * xy_norm), -z * y / (r * xy_norm), xy_norm / r);

    let cos_theta = xy_norm / r;

    let vo = ic.h / r;
    let vr = ic.vr_at_epoch(light_corrected_epoch);
    let cos_psi = ic.inc.cos() / cos_theta;
    let sin_psi = ic.kappa as f64 * (1.0 - cos_psi * cos_psi).sqrt();

    let v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);

    let f = ic.f_at_epoch(ic.epoch + (ic.epoch - light_corrected_epoch));
    let g = ic.g_at_epoch(ic.epoch + (ic.epoch - light_corrected_epoch));
    let new_position = r_vec * f + v_vec * g;

    let new_pointing = new_position / ic.r;
    Some([new_pointing[0], new_pointing[1], new_pointing[2]])
}


pub fn sync_detection_to_orbit_with_time(det: &Detection, ic: &InitialCondition) -> Option<([f64; 3], f64)> {
    let (rho, r) = optimize_rho(det, ic, ic.r - 1.0); // Keep in mind that there can be multiple solutions, but not for TNOs.

    let r_vec = det.observer_position + rho * det.rho_hat;

    let light_corrected_epoch = det.epoch - rho / SPEED_OF_LIGHT;
        
    if (r_vec[2] / r).abs() > ic.inc.sin().abs() {
        return None;
    }

    let (x, y, z) = (r_vec[0], r_vec[1], r_vec[2]);
    let xy_norm = (x * x + y * y).sqrt();

    let ahat = Vector3::new(-y / xy_norm, x / xy_norm, 0.0);
    let dhat = Vector3::new(-z * x / (r * xy_norm), -z * y / (r * xy_norm), xy_norm / r);

    let cos_theta = xy_norm / r;
    
    let vo = ic.h / r;
    let vr = ic.vr_at_epoch(light_corrected_epoch);
    let cos_psi = ic.inc.cos() / cos_theta;
    let sin_psi = ic.kappa as f64 * (1.0 - cos_psi * cos_psi).sqrt();
    
    let v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);
    
    let f = ic.f_at_epoch(ic.epoch + (ic.epoch - light_corrected_epoch));
    let g = ic.g_at_epoch(ic.epoch + (ic.epoch - light_corrected_epoch));
    let new_position = r_vec * f + v_vec * g;

    let new_pointing = new_position / ic.r;
    Some(([new_pointing[0], new_pointing[1], new_pointing[2]], light_corrected_epoch))
}
