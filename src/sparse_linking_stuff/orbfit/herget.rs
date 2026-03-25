/*

Kevin Napier, 1 October 2023

Implementation of the Herget (1965) algorithm for orbit determination.

1. Guess the topocentric distance at two observations. 
2. Calculate the dt between the two observations.
3. Use Lambert's method to calculate the orbit passing through 
   both positions, given dt.
4. Refine the orbit by perturbing the range guesses 
   to minimize the residuals of the other observations. 

The advantage of this method is that the optimization is done in
two dimensions, rather than six.

*/

use crate::detection::Detection;
use nalgebra::Vector3;

pub fn ahat(phi: f64, theta: f64) -> Vector3<f64> {
    return Vector3::new(-phi.sin(), phi.cos(), 0.0);
}

pub fn dhat(phi: f64, theta: f64) -> Vector3<f64> {
    return Vector3::new(-theta.sin() * phi.cos(), -theta.sin() * phi.sin(), theta.cos());
}

pub fn lambert(d1: Vector3<f64>, d2: Vector3<f64>, dt: f64) -> Vector3<f64> {
    return (d1 + d2) / (d1.norm() + d2.norm()) * dt;
}

pub fn herget_residuals() {
    None
}

pub fn herget(detections: &Vec<Detection>, rho_1_guess: f64, rho_n_guess: f64) -> f64 {
    let d_1 = detections[0];
    let d_n = detections[detections.len() - 1];
    let dt = d_n.epoch.epoch - d_1.epoch.epoch;

    let position_1 = d_1.observer.position + rho_1_guess * d_1.pointing();
    let position_n = d_n.observer.position + rho_n_guess * d_n.pointing();

    let orbit = lambert(position_1, position_n, dt);
}