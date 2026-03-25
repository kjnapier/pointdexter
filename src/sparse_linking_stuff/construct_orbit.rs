use crate::sparse_linking_stuff::initial_condition::{InitialCondition, nice_acos};
use crate::sparse_linking_stuff::Detection;

use spacerocks::constants::{SPEED_OF_LIGHT};
use spacerocks::transforms::{solve_for_universal_anomaly, stumpff_c, stumpff_s};
use spacerocks::SpaceRock;
use spacerocks::time::Time;

use nalgebra::Vector3;
use std::f64::consts::PI;

pub fn cost_function(det: &Detection, ic: &InitialCondition, rho: f64) -> Option<f64> {
    let ltt = rho / SPEED_OF_LIGHT;                      
    let ep = det.epoch.epoch - ltt;                            
    let r = ic.calculate_r(Time::new(ep, "utc", "jd").expect("Failed to create Time"))?;
    Some(r.powi(2)
        - rho.powi(2)
        - 2.0 * rho * det.pointing.dot(&det.observer.position)
        - det.observer.position.dot(&det.observer.position))
}

pub fn gradient(det: &Detection, _ic: &InitialCondition, rho: f64) -> f64 {
    -2.0 * (rho + det.pointing.dot(&det.observer.position))
}

pub fn optimize_rho(det: &Detection, ic: &InitialCondition, mut rho: f64) -> Option<f64> {
    const TOL: f64 = 1.0e-10;
    const MAX_ITER: usize = 50;

    let mut f = cost_function(det, ic, rho)?;
    for _ in 0..MAX_ITER {
        if f.abs() < TOL {
            break;
        }
        let g = gradient(det, ic, rho);
        if g.abs() < 1e-16 {
            break;
        }
        rho -= f / g;
        f = cost_function(det, ic, rho)?;
    }
    Some(rho)
}

pub fn construct_orbit(detection: &Detection, ic: &InitialCondition, tracklet: Option<bool>) -> Option<([f64; 3], f64)> {
    // Approximate analytic solution for the distance rho from the observer to the object at the time of detection.

    // let xi = detection.pointing().dot(&detection.observer.position);
    // let term_a = ic.r0 * ic.vr / SPEED_OF_LIGHT - xi;
    // let term_b = (term_a.powi(2) - ((ic.vr.powi(2) / SPEED_OF_LIGHT.powi(2) - 1.0) * (ic.r0.powi(2) - detection.observer.position.norm().powi(2)))).sqrt();
    // let term_c = ic.vr.powi(2) / SPEED_OF_LIGHT.powi(2) - 1.0;

    // let mut rho = term_a - term_b / term_c;
    // if rho < 0.0 {
    //     rho = term_a + term_b / term_c
    // }

    // Solve for the distance from the observer to the object numerically. More accurate.
    let rho = optimize_rho(detection, ic, ic.r0 - 1.0)?;

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

    // Only perform rate validation if tracklet is Some(true). I'm not using this for now but I'm gonna leave it here for future use.
    // if tracklet.unwrap_or(false) {
    //     let observed_ra_rate = match detection.ra_rate {
    //         Some(rate) => rate,
    //         None => {
    //             eprintln!("Warning: tracklet=true but ra_rate is None for detection {}", detection.objid);
    //             return None;
    //         }
    //     };
        
    //     let observed_dec_rate = match detection.dec_rate {
    //         Some(rate) => rate,
    //         None => {
    //             eprintln!("Warning: tracklet=true but dec_rate is None for detection {}", detection.objid);
    //             return None;
    //         }
    //     };
        
    //     let mut rock = SpaceRock::from_spherical(
    //         &detection.objid.clone(), 
    //         phi,
    //         theta,
    //         ic.r0,
    //         ic.vr, 
    //         ic.vo,
    //         ic.psi,
    //         detection.epoch.clone(),
    //         "J2000",
    //         "SSB"
    //     ).expect("Failed to create SpaceRock from spherical coordinates");

    //     rock.analytic_propagate(&detection.observer.epoch); 

    //     let obs = match rock.observe(&detection.observer) {
    //         Ok(observation) => observation,
    //         Err(e) => {println!("Error observing rock: {}", e);
    //             return None;
    //         }
    //     };
        
    //     let theoretical_ra_rate = obs.ra_rate()?;
    //     let theoretical_dec_rate = obs.dec_rate()?;

    //     // println!("Theoretical RA rate: {}, Observed RA rate: {}", theoretical_ra_rate, observed_ra_rate);
    //     // println!("Theoretical Dec rate: {}, Observed Dec rate: {}", theoretical_dec_rate, observed_dec_rate);

    //     // Validate against detection's observed rates (within 10%)
    //     let ra_rate_match = rates_within_tolerance(
    //         observed_ra_rate, 
    //         theoretical_ra_rate, 
    //         0.25 // 25% tolerance
    //     );
        
    //     let dec_rate_match = rates_within_tolerance(
    //         observed_dec_rate, 
    //         theoretical_dec_rate, 
    //         0.25 // 25% tolerance
    //     );
        
    //     // Both rates must match within tolerance
    //     if !ra_rate_match || !dec_rate_match {
    //         return None;
    //     }
    // }


    // Assert that psi = inclination. This is an assertion that will get worse as the object gets to larger latitudes.
    let mut inclination = ic.psi.abs();
    let kappa = ic.psi.signum();


    // // Eary exit if we are at an unphysical theta given the inclination.
    // if inclination < theta.abs() {
    //     // println!("Early exit: inclination {} is less than theta {}", inclination, theta);
    //     inclination = ic.p4_min?;
    //     // println!("Early exit: inclination {} is less than theta {}", inclination, theta);
    //     // return None;
    // }

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

    // println!("Mean anomaly: {}, vr: {}, ic: {:?}", mean_anomaly, vr, ic.id);



    // Get the  full velocity vector
    let cos_psi = inclination.cos() / theta.cos();
    let mut sin_psi = kappa * (1.0 - cos_psi.powi(2)).sqrt();
    let mut v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);

    // println!("cos_psi: {}, sin_psi: {}, v_vec: {:?}, ic: {:?}", cos_psi, sin_psi, v_vec, ic.id);
    // println!("v_vec: {:?}, ic: {:?}", v_vec, ic.id);

    // Using r and v, propagate the orbit to the time of the initial condition w/ the Lagrange coefficients.
    let alpha = -2.0 * ic.energy / ic.mu;
    // let chi = solve_for_universal_anomaly(r, vr, alpha, ic.mu, dt, 1e-10, 100).expect("Failed to solve for universal anomaly");

    let chi = match solve_for_universal_anomaly(r, vr, alpha, ic.mu, dt, 1e-6, 100) {
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
    // let fdot = chi * ic.mu.sqrt() / (r * ic.r0) * (alpha * chi.powi(2) * stumpff_s_val - 1.0);
    // let gdot = 1.0 - chi.powi(2) * stumpff_c_val / r;
    let mut new_position = f * r_vec + g * v_vec;

    // let new_velocity = fdot * r_vec + gdot * v_vec;

    // if new_velocity[2].signum() != kappa {
    //     sin_psi = -sin_psi;
    //     v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);
    //     new_position = f * r_vec + g * v_vec;
    // }
    // let v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);
    // let new_position = f * r_vec + g * v_vec;

    let new_pointing = new_position.normalize();
    // println!("New pointing: {:?} for ic {:?}", new_pointing, ic.id);

    Some(([new_pointing.x, new_pointing.y, new_pointing.z], theta))

}

fn rates_within_tolerance(observed: f64, theoretical: f64, rel_tolerance: f64) -> bool {
    // For non-zero rates, check relative difference
    let relative_diff = (observed - theoretical).abs() / theoretical.abs().max(observed.abs());
    relative_diff <= rel_tolerance
}

pub fn detection_to_spacerock<'a>(
    detection: &Detection,
    fakes: &'a [SpaceRock],
) -> Option<SpaceRock> {
    let fakeid = detection.fakeid.as_ref()?;
    fakes.iter().find(|rock| rock.name == *fakeid).cloned()
}




pub fn construct_orbit_from_fake(detection: &Detection, ic: &InitialCondition, fakes: &[SpaceRock]) -> Option<([f64; 3], f64)> {
    // Approximate analytic solution for the distance rho from the observer to the object at the time of detection.

    // let xi = detection.pointing().dot(&detection.observer.position);
    // let term_a = ic.r0 * ic.vr / SPEED_OF_LIGHT - xi;
    // let term_b = (term_a.powi(2) - ((ic.vr.powi(2) / SPEED_OF_LIGHT.powi(2) - 1.0) * (ic.r0.powi(2) - detection.observer.position.norm().powi(2)))).sqrt();
    // let term_c = ic.vr.powi(2) / SPEED_OF_LIGHT.powi(2) - 1.0;

    // let mut rho = term_a - term_b / term_c;
    // if rho < 0.0 {
    //     rho = term_a + term_b / term_c
    // }

    // Solve for the distance from the observer to the object numerically. More accurate.
    let rho = optimize_rho(detection, ic, ic.r0 - 1.0)?;

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

    
    // Get the true inclination if the detection is associated with a fake by calling detection_to_spacerock
    let mut fake_rock = detection_to_spacerock(detection, fakes)?;
    fake_rock.change_reference_plane("J2000");

    // Assert that psi = inclination. This is an assertion that will get worse as the object gets to larger latitudes.
    let mut inclination = fake_rock.inc().abs();
    let kappa = ic.psi.signum();


    // // Eary exit if we are at an unphysical theta given the inclination.
    // if inclination < theta.abs() {
    //     // println!("Early exit: inclination {} is less than theta {}", inclination, theta);
    //     inclination = ic.p4_min?;
    //     // println!("Early exit: inclination {} is less than theta {}", inclination, theta);
    //     // return None;
    // }

    let theta_abs = theta.abs();

    if inclination < theta_abs {
        inclination = match (ic.p4_min, ic.p4_max) {
            (Some(min), Some(max)) => {
                if min >= theta_abs && max >= theta_abs {
                    min.min(max)
                } else if min >= theta_abs {
                    min
                } else if max >= theta_abs {
                    max
                } else {
                    // println!("Early exit: no valid p4_min/max ≥ |theta| = {}, min = {}, max = {}: Det orbit id = {}, ic_id = {}, inc = {}", theta_abs, min, max, detection.orbit_id.unwrap(), ic.id, inclination);
                    inclination += 0.001; // Adjust inclination slightly to avoid early exit
                    if inclination < theta_abs {
                        // println!("Adjusted inclination still less than theta, returning None");
                        return None;
                    }
                    // println!("Adjusted inclination to {}", inclination);
                    // If we reach here, we can return the adjusted inclination
                    // but we should still return None to avoid further processing
                    // This is a workaround to avoid early exit in the current logic
                    // but it may need to be revisited in the future.
                    inclination
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

    // println!("Mean anomaly: {}, vr: {}, ic: {:?}", mean_anomaly, vr, ic.id);



    // Get the  full velocity vector
    let cos_psi = inclination.cos() / theta.cos();
    let mut sin_psi = kappa * (1.0 - cos_psi.powi(2)).sqrt();
    let mut v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);

    // println!("cos_psi: {}, sin_psi: {}, v_vec: {:?}, ic: {:?}", cos_psi, sin_psi, v_vec, ic.id);
    // println!("v_vec: {:?}, ic: {:?}", v_vec, ic.id);

    // Using r and v, propagate the orbit to the time of the initial condition w/ the Lagrange coefficients.
    let alpha = -2.0 * ic.energy / ic.mu;
    // let chi = solve_for_universal_anomaly(r, vr, alpha, ic.mu, dt, 1e-10, 100).expect("Failed to solve for universal anomaly");

    let chi = match solve_for_universal_anomaly(r, vr, alpha, ic.mu, dt, 1e-6, 100) {
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
    // let fdot = chi * ic.mu.sqrt() / (r * ic.r0) * (alpha * chi.powi(2) * stumpff_s_val - 1.0);
    // let gdot = 1.0 - chi.powi(2) * stumpff_c_val / r;
    let mut new_position = f * r_vec + g * v_vec;

    // let new_velocity = fdot * r_vec + gdot * v_vec;

    // if new_velocity[2].signum() != kappa {
    //     sin_psi = -sin_psi;
    //     v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);
    //     new_position = f * r_vec + g * v_vec;
    // }
    // let v_vec = vr * r_vec / r + vo * (cos_psi * ahat + sin_psi * dhat);
    // let new_position = f * r_vec + g * v_vec;

    let new_pointing = new_position.normalize();
    // println!("New pointing: {:?} for ic {:?}", new_pointing, ic.id);

    Some(([new_pointing.x, new_pointing.y, new_pointing.z], theta))

}