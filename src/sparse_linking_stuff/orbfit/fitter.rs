use crate::detection::Detection;

use spacerocks::Time;
use spacerocks::SpaceRock;
use spacerocks::coordinates::{Origin, ReferencePlane};

use nalgebra::{DMatrix, DVector};
use serde::{Serialize, Deserialize};


pub fn residuals(detections: &Vec<&Detection>, theta: &Vec<f64>, epoch: &Time) -> Option<DVector<f64>> {

    // let state = StateVector::new(theta[0], theta[1], theta[2], theta[3], theta[4], theta[5]);
    // let trial = SpaceRock::from_state("rock", state, epoch.clone(), "J2000", "SSB");
    // let trial = SpaceRock::from_state("rock", state, epoch.clone(), &ReferencePlane::J2000, &Origin::SSB);
    let trial = SpaceRock::from_xyz("rock", theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], epoch.clone(), "J2000", "SSB").expect("Error in constructing SpaceRock");

    let mut residuals: DVector<f64> = DVector::zeros(detections.len());
    let mut rock = trial.clone();
    
    for (idx, detection) in detections.iter().enumerate() {

        if let Err(_) = rock.analytic_propagate(&Time::new(detection.epoch, "utc", "jd").expect("Failed to create Time for propagation")) {
            return None; // Return None if propagation fails
        }

        // rock.analytic_propagate(&detection.epoch);
        let astro = rock.observe(detection.observer.as_ref()?).unwrap();

        let ra_residual = (astro.ra() - detection.ra?) / detection.ra_ucty?;
        let dec_residual = (astro.dec() - detection.dec?) / detection.dec_ucty?;

        residuals[idx] = (ra_residual.powi(2) + dec_residual.powi(2)).sqrt();
    }
    Some(residuals)
}

pub fn radec_residuals(detections: &Vec<&Detection>, theta: &Vec<f64>, epoch: &Time) -> Option<Vec<(f64, f64)>> {

    // let state = StateVector::new(theta[0], theta[1], theta[2], theta[3], theta[4], theta[5]);
    // let trial = SpaceRock::from_state("rock", state, epoch.clone(), "J2000", "SSB");
    // let trial = SpaceRock::from_state("rock", state, epoch.clone(), &ReferencePlane::J2000, &Origin::SSB);
    let trial = SpaceRock::from_xyz("rock", theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], epoch.clone(), "J2000", "SSB").expect("Error in constructing SpaceRock");

    let mut residuals: Vec<(f64, f64)> = Vec::new();
    let mut rock = trial.clone();
    
    for detection in detections.iter() {

        // rock.analytic_propagate(&detection.epoch);
        if let Err(_) = rock.analytic_propagate(&Time::new(detection.epoch, "utc", "jd").expect("Failed to create Time for propagation")) {
            return None; // Return None if propagation fails
        }
        let astro = rock.observe(detection.observer.as_ref()?).unwrap();


        let ra_residual = (astro.ra() - detection.ra?) / detection.ra_ucty?;
        let dec_residual = (astro.dec() - detection.dec?) / detection.dec_ucty?;

        residuals.push((ra_residual, dec_residual));
    }
   Some(residuals)
}

// compute the jacobian of the residuals with respect to the state vector
pub fn jacobian(detections: &Vec<&Detection>, theta: &Vec<f64>, epoch: &Time) -> Option<DMatrix<f64>> {

    let mut jac = DMatrix::zeros(detections.len(), 6);
    let mut theta_plus = theta.clone();
    let mut theta_minus = theta.clone();
    let eps = (f64::EPSILON).sqrt();

    for i in 0..6 {
        theta_plus[i] += eps;
        theta_minus[i] -= eps;
        let res_plus = residuals(detections, &theta_plus, epoch)?;
        let res_minus = residuals(detections, &theta_minus, epoch)?;
        let deriv = (res_plus - res_minus) / (2.0 * eps);
        for j in 0..detections.len() {
            jac[(j, i)] = deriv[j]
        }
        theta_plus[i] -= eps;
        theta_minus[i] += eps;
    }

    Some(jac)
}

// Right now, I am just using the sum of the squares of the raw residuals as the cost function to get chi-squared
// In spacerocks, I use the Mahalanobis distance... This should be fine for now
pub fn cost(detections: &Vec<&Detection>, theta: &Vec<f64>, epoch: &Time) -> Option<f64> {

    // let state = StateVector::new(theta[0], theta[1], theta[2], theta[3], theta[4], theta[5]);
    // let trial = SpaceRock::from_state("rock", state, epoch.clone(), &ReferencePlane::J2000, &Origin::SSB);
    let trial = SpaceRock::from_xyz("rock", theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], epoch.clone(), "J2000", "SSB").expect("Error in constructing SpaceRock");

    let res = residuals(detections, theta, epoch)?;
    let csq = res.dot(&res);
    Some(csq)
}



#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FitResult {
    pub chisq: f64,
    pub rock: SpaceRock,
    pub niter: usize,
    pub dof: f64,
    pub residuals: Vec<f64>,
    pub ra_residuals: Vec<f64>,
    pub dec_residuals: Vec<f64>,
    pub covariance: DMatrix<f64>,
}

pub fn fit_orbit_lm(detections: &Vec<&Detection>, initial_guess: Vec<f64>, epoch: &Time) -> Option<FitResult> {    
    
    let csq_tol = 1e-1;
    let grad_tol = 1e-3;
    let theta_tol = 1e-4;
    let rho_accept = 0.0;

    // let csq_tol = 1e-3;
    // let grad_tol = 1e-15;
    // let theta_tol = 1e-15;
    // let rho_accept = 0.0;

    let maxiter = 100;
    let mut niter = 0;

    let dof = detections.len() as f64 * 2.0 - 6.0;
    let mut theta = initial_guess.clone();
    
    let mut csq = cost(detections, &theta, epoch)? / dof;
    // println!("Initial chisq: {}", csq);
    let mut new_csq: f64 = csq;
    
    let mut csq_change: f64 = 1000.0;
    
    let mut lambda: f64 = 0.0001;

    while csq_change.abs() > csq_tol {
    // while true {

        let res = residuals(detections, &theta, epoch)?;
        let mut j = jacobian(detections, &theta, epoch)?;

        let mut a = j.clone().transpose() * j.clone();
        for idx in 0..6 {
            a[(idx, idx)] += lambda;
        }
        
        let a_inv = match a.try_inverse() {
            Some(a_inv) => a_inv,
            None => return None,
        };

        //check if the gradient has converged
        let grad = j.clone().transpose() * res.clone();
        let mut max_grad = 0.0;
        for idx in 0..6 {
            if grad[idx].abs() > max_grad {
                max_grad = grad[idx].abs();
            }
        }
        if max_grad < grad_tol {
            // println!("Gradient converged");
            break;
        }


        let h = a_inv * grad.clone();

        // if no parameter is being changed by more than theta_tol, break
        let mut max_h = 0.0;
        for idx in 0..6 {
            let perturbation = (h[idx] / theta[idx]).abs();
            if perturbation > max_h {
                max_h = perturbation;
            }
        }
        if max_h < theta_tol {
            // println!("Parameters converged");
            break;
        }

        for idx in 0..6 {
            theta[idx] -= h[idx];
        }

        new_csq = cost(detections, &theta, epoch)? / dof;
        let csq_change = new_csq - csq;

        let rho = (csq - new_csq) / (h.clone().transpose() * (lambda * h.clone() + grad.clone())).norm();
        if rho > rho_accept {
            // accept the step and shrink lambda
            lambda *= 0.5;
        } else {
            // reject the step, increase lambda, and reset the parameters
            lambda *= 2.0;
            for idx in 0..6 {
                theta[idx] += h[idx];
            }
        }

        csq = new_csq;

        niter += 1;
        if niter >= maxiter {
            break;
        }
        
    }

    // let rock = SpaceRock::from_state("rock", StateVector::new(theta[0], theta[1], theta[2], theta[3], theta[4], theta[5]), epoch.clone(), "J2000", "SSB");
    // let rock = SpaceRock::from_state("rock", StateVector::new(theta[0], theta[1], theta[2], theta[3], theta[4], theta[5]), epoch.clone(), &ReferencePlane::J2000, &Origin::SSB);
    let rock = SpaceRock::from_xyz("rock", theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], epoch.clone(), "J2000", "SSB").expect("Error in constructing SpaceRock");

    let j = jacobian(detections, &theta, epoch)?;
    let a = j.clone().transpose() * j.clone();
    let cov = a.pseudo_inverse(1e-10).unwrap();

    let rd_resid = radec_residuals(detections, &theta, epoch)?;

    return Some(FitResult {
        chisq: new_csq * dof,
        dof: dof,
        residuals: residuals(detections, &theta, epoch)?.iter().map(|r| *r).collect::<Vec<_>>(),
        ra_residuals: rd_resid.iter().map(|(ra, _)| *ra).collect::<Vec<_>>(),
        dec_residuals: rd_resid.iter().map(|(_, dec)| *dec).collect::<Vec<_>>(),
        rock: rock,
        niter: niter,
        covariance: cov
    });


}