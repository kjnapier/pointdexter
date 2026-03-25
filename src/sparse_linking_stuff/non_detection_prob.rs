use crate::sparse_linking_stuff::ccd::compute_chip;
use crate::detection::Detection;
use spacerocks::{SpaceRock, Simulation, observing::Observatory};
use spacerocks::time::Time;
use spacerocks::SpiceKernel;
use statrs::distribution::{Binomial, DiscreteCDF};
use std::collections::{HashMap, HashSet};
use ordered_float::OrderedFloat;
use crate::sparse_linking_stuff::config::cli::Config;
use crate::sparse_linking_stuff::ExposureRow;

// Detection efficiency model given some fit parameters using the fakes
pub fn detection_efficiency(mag: f64, eta0: f64, m50: f64, sigma: f64) -> f64 {
    eta0 / (1.0 + 10f64.powf((mag - m50) / sigma))
}


// Estimate the angular distance between two points on the sky in degrees
fn angular_distance_deg(ra1: f64, dec1: f64, ra2: f64, dec2: f64) -> f64 {
    let d_ra = ra1 - ra2;
    let d_dec = dec1 - dec2;
    ((d_dec.powi(2) + (d_ra * dec2.cos()).powi(2)).sqrt()).to_degrees()
}


// Returns the cdf of the binomial distribution for the number of detections k
// given n trials and detection probability p. This is the probability of
// detecting k or fewer times in n trials.
pub fn non_detection_prob_pass(
    rock: &SpaceRock,
    dets_by_expnum: &HashMap<i64, Vec<&Detection>>, 
    exposures: &[ExposureRow],
    kernel: &SpiceKernel,
    config: &Config,
) -> f64 {

    // Set up a simulation and add the best fit rock
    let name = rock.name.clone();
    let mut sim = Simulation::giants(&rock.epoch, "J2000", "SSB", kernel).unwrap();
    sim.add(rock.clone()).unwrap();

    // Set up the observatory
    let observatory = Observatory::from_obscode("W84").unwrap();

    let mut n_trials = 0usize;
    let mut k_hits   = 0usize;

    // collect mags from matched detections 
    let mut matched_mags: Vec<f64> = Vec::new();

    // Run through exposures in in order
    let mut exps_sorted = exposures.to_vec();
    exps_sorted.sort_by(|a, b| a.epoch.total_cmp(&b.epoch));
    for exp in &exps_sorted {

        // 1) propagate to the exposure epoch
        let t = Time::new(exp.epoch, "utc", "jd").unwrap();
        let observer = observatory.at(&t, "J2000", "SSB", kernel).unwrap();
        sim.integrate(&t);
        let mut body = sim.get_particle(&name).unwrap().clone();
        let obs = body.observe(&observer).unwrap();

        // predicted RA/Dec at this epoch
        let ra_obs = obs.ra();
        let dec_obs = obs.dec();


        // 2) does it land on a valid DECam CCD for this exposure? (Trial)
        let (_chip, ccd_id) = compute_chip(ra_obs, dec_obs, exp.radeg.to_radians(), exp.decdeg.to_radians());

        // Check if this exposure AND this specific CCD are valid 
        let exposure_valid = !config.invalid_expnums.contains(&exp.expnum);
        let ccd_valid = !config.invalid_ccds.contains(&ccd_id); 

        if !exposure_valid || !ccd_valid {
            continue;
        }
        // If the expsosure didn't get filtered out and the CCD is working, count this as a trial
        n_trials += 1;


        // 3) any detection in this expnum within 5" ? 
        if let Some(cands) = dets_by_expnum.get(&exp.expnum) {
            let mut hit = false;
            for d in cands {
                if let (Some(ra), Some(dec)) = (d.ra, d.dec) {
                    let dist = angular_distance_deg(ra, dec, ra_obs, dec_obs);
                    if dist < (5.0 / 3600.0) {
                        hit = true;
                        if let Some(m) = d.mag {
                            matched_mags.push(m);
                        }
                    }
                }
            }
            if hit { k_hits += 1; }
        }
    }

    // If k_hits is >= 30, return a cdf of 1
    if k_hits >= 20 {
        println!("Non-detection test: k={}, n={} is very high, returning p=1.0", k_hits, n_trials);
        return 1.0;
    }

    if n_trials == 0 {
        // No exposures where the object would have been seen on a valid chip.
        println!("Non-detection test: no valid exposures (n=0).");
        return 0.0;
    }

    // 4) Estimate p from matched mags by taking the mean mag from the detections, as estimated by our source extraction code... This needs to be
    // modified to account for the fact that the magnitudes used for the detection efficiency fits are the injected mags on Ed's side,
    // not the recovered mags after SE. For now, we will just use the mean of the matched mags as a proxy.
    let mean_mag = if matched_mags.is_empty() {
        println!("Non-detection test: no matches found (k=0).");
        return 0.0;
    } else {
        matched_mags.iter().sum::<f64>() / (matched_mags.len() as f64)
    };

    if mean_mag > 28.0 {
        println!("Non-detection test: Too faint to be real (mean_mag > 28).");
        return 0.0;
    }

    // Add in 0.2 mag to be conservative
    let eta0 = config.eta0;
    let m50  = config.m50;
    let sigma = config.sigma;
    let prob_threshold = config.prob_threshold;

    let test_mag = matched_mags.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    //let p = detection_efficiency(mean_mag + 0.3, eta0, m50, sigma);
    let p = detection_efficiency(test_mag+0.2, eta0, m50, sigma);

    // 5) Binomial CDF P[K <= k | n, p]
    let binom = match Binomial::new(p, n_trials as u64) {
        Ok(b) => b,
        Err(e) => {
            println!(
                "Non-detection test: invalid binomial params (p={}, n={}): {:?}",
                p, n_trials, e
            );
            return 0.0;
        }
    };
    
    let cdf = binom.cdf(k_hits as u64);

    println!(
        "Non-detection test: k={}, n={}, mean_mag={:.2}, p={:.3}, CDF={:.6}, threshold={}",
        k_hits, n_trials, mean_mag, p, cdf, prob_threshold
    );

    cdf
}
