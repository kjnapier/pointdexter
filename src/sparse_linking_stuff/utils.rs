use rand::seq::IndexedRandom;
use rand::rng;
use std::collections::HashMap;
use nalgebra::{DMatrix, DVector};

use crate::initial_condition::InitialCondition;
use crate::detection::Detection;
use crate::sparse_linking_stuff::orbfit::fitter::FitResult as HsFit;

use spacerocks::orbfit::fitter::FitResult as SrFit;
use spacerocks::nbody::Simulation;
use spacerocks::SpaceRock;
use spacerocks::SpiceKernel;
use spacerocks::Time;

pub const ARCSEC_PER_RAD: f64 = 3600.0 * 180.0 / std::f64::consts::PI;

pub fn get_n_random_orbits<'a>(ics: &'a [InitialCondition], n: usize) -> Vec<&'a InitialCondition> {
    // Check if we have enough points to sample
    if ics.len() < n {
        return ics.iter().collect();
    }

    let mut rng = rng();
    // Using choose_multiple with copied references
    ics.choose_multiple(&mut rng, n).collect()
}


pub fn to_sparse_linking_stuff_fit(sr: SrFit) -> HsFit {
    HsFit {
        chisq:         sr.chisq,
        rock:          sr.rock,
        niter:         sr.niter,
        dof:           sr.dof,
        residuals:     sr.residuals,
        ra_residuals:  sr.ra_residuals,
        dec_residuals: sr.dec_residuals,
        covariance:    sr.covariance,
    }
}

pub fn index_detections_by_expnum<'a>(
    dets_all: &'a [Detection]
) -> HashMap<i64, Vec<&'a Detection>> {
    let mut map: HashMap<i64, Vec<&Detection>> = HashMap::new();
    for d in dets_all {
        if let Some(expnum) = d.expnum {
            map.entry(expnum).or_default().push(d);
        }
    }
    map
}

/// Same here, I need to move this somewhere else. Just a slightly modified version of the SR residuals function, this takes a spacerock
/// while the SR function takes a 7 element state vector.
pub fn residuals(detections: &Vec<&Detection>, rock: &mut SpaceRock, kernel: &SpiceKernel) 
    -> Result<(DVector<f64>, Vec<f64>, Vec<f64>), Box<dyn std::error::Error>> {

    let mut sim = Simulation::giants(&rock.epoch, "J2000", "SSB", kernel)?;

    let trial: SpaceRock = rock.clone();
    sim.integrate(&trial.epoch);
    sim.add(trial)?;

    let mut residuals: DVector<f64> = DVector::zeros(detections.len());

    // Create vectors to store the raw RA and Dec residuals
    let mut ra_residuals: Vec<f64> = Vec::new();
    let mut dec_residuals: Vec<f64> = Vec::new();
    
    for (idx, detection) in detections.iter().enumerate() {

        let observation = detection.to_observation().unwrap();
        let observed_parameters = DVector::from_vec(vec![observation.ra(), observation.dec()]);
          

        // Get the rock to the epoch of the detection
        sim.integrate(&Time::new(detection.epoch, "utc", "jd")?);
        let mut rock = sim.get_particle("rock")?.clone();

        // calculate the model observations
        let astro = rock.observe(detection.observer.as_ref().ok_or("Detection missing observer")?)?;        
        let model_parameters = DVector::from_vec(vec![astro.ra(), astro.dec()]);
           

        let d = observed_parameters - model_parameters;
        // Store the raw differences 
        ra_residuals.push(d[0]);   // RA difference in radians
        dec_residuals.push(d[1]);  // Dec difference in radians

        // Convert 0.15 arcseconds to radians -- this is just a standard uncertainty for now
        let smear = 0.15 / 3600.0 * (std::f64::consts::PI / 180.0);

        // Create covariance matrix as a 2x2 diagonal matrix
        let cov = DMatrix::from_diagonal(&DVector::from_vec(vec![smear.powi(2), smear.powi(2)]));

        // Calculate inverse of covariance matrix
        let inverse_covariance = cov.try_inverse().expect("Matrix should be invertible");

        // This is the mahalanobis distance calculation
        let m = &d.transpose() * &inverse_covariance * &d;
        residuals[idx] = m[0].sqrt();
    }
    // Return the residuals and the raw RA and Dec residuals
    Ok((residuals, ra_residuals, dec_residuals))
}