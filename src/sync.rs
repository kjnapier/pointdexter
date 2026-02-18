use nalgebra::Vector3;
use crate::detection::Detection;
use crate::initial_condition::InitialCondition;


use spacerocks::data::{SPEED_OF_LIGHT};


pub fn optimize_rho_v2(det: &Detection, ic: &InitialCondition, mut rho0: f64, sgn: f64) -> Option<(f64, f64)> {
    let tol = 1e-6;         
    let max_iter = 100;
    let inv_c = 1.0 / SPEED_OF_LIGHT;

    let dot = det.rho_hat_dot_observer_position;
    let dot2 = dot * dot;
    let obs2 = det.observer_distance_squared;

    for _ in 0..max_iter {
        let dt = rho0 * inv_c;
        let r = ic.r_at_epoch(det.epoch - dt);

        // disc = r*r + dot2 - obs2
        // mul_add can be a wash, but often maps nicely to FMA where available
        let disc = r.mul_add(r, dot2 - obs2);
        if disc < 0.0 { return None; }

        let rho = -dot + sgn * disc.sqrt();

        let d = rho - rho0;
        if d.abs() < tol {
            return Some((rho, r));
        }
        rho0 = rho;
    }
    None
}




#[inline(never)]
pub fn sync_detection_to_orbit(det: &Detection, ic: &InitialCondition) -> Option<[f64; 3]> {
    // let (rho, r) = optimize_rho(det, ic, ic.r - 1.0);
    let rho_opt = optimize_rho_v2(det, ic, ic.r - 1.0, 1.0); // Keep in mind that there can be multiple solutions, but not for TNOs.
    let (rho, r) = match rho_opt {
        Some((rho_val, r_val)) => (rho_val, r_val),
        None => return None,
    };
    
    let r_vec = det.observer_position + rho * det.rho_hat;
    

    let light_corrected_epoch = det.epoch - rho / SPEED_OF_LIGHT;
    let (f, g) = ic.fg_at_epoch(ic.epoch + (ic.epoch - light_corrected_epoch));

    let (x, y, z) = (r_vec[0], r_vec[1], r_vec[2]);
    let xy_norm = (x * x + y * y).sqrt();

    let ahat = Vector3::new(-y / xy_norm, x / xy_norm, 0.0);
    let dhat = Vector3::new(-z * x / (r * xy_norm), -z * y / (r * xy_norm), xy_norm / r);

    let cos_theta = xy_norm / r;

    let vo = ic.h / r;
    let vsq = (ic.energy + ic.mu / r).max(0.0) * 2.0; // ensure non-negative vsq
    let vr = (vsq - vo * vo).max(0.0).sqrt(); // ensure non-negative argument for sqrt
    
    if (r_vec[2] / r).abs() > ic.sin_latitude_threshold {
        return None;
    }

    let cos_psi = ic.cos_inc / cos_theta;
    let sin_psi = ic.kappa as f64 * (1.0 - cos_psi * cos_psi).sqrt();

    let v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);

    
    let new_position = r_vec * f + v_vec * g;

    let new_pointing = new_position / ic.r;
    Some([new_pointing[0], new_pointing[1], new_pointing[2]])
}