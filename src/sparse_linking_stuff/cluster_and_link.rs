use std::time::Instant;
use std::collections::{HashMap, HashSet};
use std::sync::{Mutex, Arc};
use std::sync::atomic::{AtomicUsize, Ordering, AtomicU64};
use std::io::BufWriter;
// use itertools::Itertools;
use std::mem; 
use ordered_float::OrderedFloat;

use clap::Parser;
use indicatif::ProgressIterator;
use rayon::prelude::*;
use rayon::iter::IntoParallelRefIterator;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use nalgebra::{DVector, DMatrix};
// use lazy_static::lazy_static;
use std::collections::BTreeSet;


use flate2::write::GzEncoder;
use flate2::Compression;
use spice;
use kiddo::KdTree;
use kiddo::immutable::float::kdtree::ImmutableKdTree;
use kiddo::SquaredEuclidean;

use serde_yaml;
use serde::{Serialize, Deserialize};
use crate::initial_condition::InitialCondition;
use crate::detection::Detection;
// use Clustering-and-Linking::detection::Detection;
use crate::sync::sync_detection_to_orbit as construct_orbit;
use crate::sparse_linking_stuff::config::{Config, Cli};
use crate::sparse_linking_stuff::ccd::*;
use crate::sparse_linking_stuff::non_detection_prob::*;
use crate::sparse_linking_stuff::io::ExposureRow;
use crate::sparse_linking_stuff::grid::GridCacheKey;
use crate::sparse_linking_stuff::utils::*;

use spacerocks::StateVector;
// use spacerocks::orbfit::fitter::residuals as residuals_sr;
use bloomfilter::Bloom;
use std::fs::File;
use std::path::PathBuf;
use crate::sparse_linking_stuff::Link;
use spacerocks::nbody::Simulation;
use spacerocks::orbfit::gauss::gauss;
use spacerocks::SpiceKernel;

use spacerocks::orbfit::fitter::{FitResult as FitResultSr};
use crate::sparse_linking_stuff::orbfit::fitter::{FitResult};

use spacerocks::orbfit::fitter::fit_orbit_lm;
use spacerocks::orbfit::fitter::FitResult as SrFit;
use crate::sparse_linking_stuff::orbfit::fitter::FitResult as HsFit;

use spacerocks::observing::observation::ObservationType;
use spacerocks::observing::{Observer, Observation};
use spacerocks::SpaceRock;
use itertools::Itertools;
use spacerocks::constants::MU_BARY;
use spacerocks::Time;

use std::hash::{Hash, Hasher};
use crate::sparse_linking_stuff::grid::{Cell, IcBoundsKey, refine_with_grid_cache_sparse};
use crate::sparse_linking_stuff::trajectory::*;
use crate::sparse_linking_stuff::utils::*;
use crate::sparse_linking_stuff::pipeline_metrics::PipelineMetrics;



pub fn run_ic_cluster_and_trajectory_search(
    ic: &InitialCondition,
    detections: &[Detection],
    det_map: &HashMap<u64, Detection>,
    epsilon_arcsec: f64,
    config: &Config,
    pipeline_metrics: &Arc<PipelineMetrics>,
) -> Vec<Trajectory> {

    let mut valid_trajs = Vec::new();
    let eps = epsilon_arcsec / ARCSEC_PER_RAD;
    let cluster_radius = 2.0 * (1.0 - eps.cos());

    let mut point_map: HashMap<u64, [f64; 3]> = HashMap::new();
    let mut center_list = Vec::new();

    // --- 1: collect all valid points first ---
    let mut points: Vec<[f64; 3]> = Vec::new();
    let mut intids: Vec<u64> = Vec::new();

    for det in detections {
        if let Some(pt) = construct_orbit(det, ic) {
            if pt.iter().any(|v| !v.is_finite()) {
                println!("Skipping invalid KD point: {:?} from det {}", pt, det.intid);
                continue;
            }
            points.push(pt);
            intids.push(det.intid as u64);
            point_map.insert(det.intid as u64, pt);
            center_list.push((det.intid as u64, pt));
        }
    }

    if points.is_empty() {
        return valid_trajs;
    }

    // --- 2: build the immutable tree from the collected points ---
    type Tree = ImmutableKdTree<f64, u64, 3, 32>;
    let tree: Tree = ImmutableKdTree::new_from_slice(&points);

    // --- 3: cluster search ---
    for &(center_id, center_pt) in &center_list {
        let center_epoch = det_map.get(&center_id).unwrap().epoch;
        let all_neighbors = tree.within_unsorted::<SquaredEuclidean>(&center_pt, cluster_radius);

        // Find all the neihbors that satisfy our combined spatial and temporal criteria. 
        // The spatial criteria is a time-dependent cone, where the radius of the cone grows linearly with time from the center detection.
        // The temporal criteria is just that the detections must be within t_max/2 of the center detection.
        let neighbors: Vec<_> = all_neighbors
            .into_iter()
            .filter(|nb| {
                let intid = intids[nb.item as usize];
                let det = det_map.get(&intid).unwrap();
                let dt = (det.epoch - center_epoch).abs();

                if dt > config.t_max / 2.0 + 1.0 {
                    return false;
                }

                let pt2 = *point_map.get(&intid).unwrap();
                let sep_arcsec = angular_separation(center_pt, pt2)
                    .iter()
                    .map(|x| x.powi(2))
                    .sum::<f64>()
                    .sqrt();

                if sep_arcsec < 0.1 * config.epsilon_arcsec {
                    return true;
                }

                sep_arcsec < config.epsilon_arcsec * dt / (config.t_max / 2.0)
            })
            .collect();

        if neighbors.len() < config.min_detections || neighbors.len() > config.max_detections {
            continue;
        }

        let mut nights = HashSet::new();
        for nb in &neighbors {
            let intid = intids[nb.item as usize];
            let d = det_map.get(&intid).unwrap();
            nights.insert(d.nite);
        }

        if nights.len() < config.min_nites {
            continue;
        }

        let mut expnums = HashSet::new();
        for nb in &neighbors {
            let intid = intids[nb.item as usize];
            let d = det_map.get(&intid).unwrap();
            if let Some(expnum) = d.expnum {
                expnums.insert(expnum);
            }
        }

        if expnums.len() < config.min_detections {
            continue;
        }

        let times: Vec<f64> = neighbors.iter()
            .map(|nb| det_map.get(&intids[nb.item as usize]).unwrap().epoch)
            .collect();
        let min_t = *times.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        let max_t = *times.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        if (max_t - min_t) < config.min_duration {
            continue;
        }

        let cluster_ids: Vec<usize> = neighbors.iter()
            .map(|nb| intids[nb.item as usize] as usize)
            .collect();
        pipeline_metrics.track_cluster(&cluster_ids);

        let raw_vec: Vec<RawDetection> = neighbors
            .iter()
            .map(|nb| {
                let intid = intids[nb.item as usize];
                let det = det_map.get(&intid).unwrap();
                let pt2 = *point_map.get(&intid).unwrap();
                let [dphi, dtheta] = angular_separation(center_pt, pt2);
                let mag = det.mag.unwrap_or(0.0);
                RawDetection {
                    id: intid as usize,
                    delta_phi: dphi,
                    delta_theta: dtheta,
                    epoch: det.epoch,
                    magnitude: mag,
                    night: det.nite.unwrap_or(0),
                    expnum: det.expnum.unwrap_or(0),
                }
            })
            .collect();

        let meta_pts = create_meta_points(&raw_vec, 1.0, 0.2);
        let trajs = graph_based_trajectory_search(&meta_pts, 0.6, 45., &config);

        for traj in trajs {
            let ids: Vec<usize> = traj.iter()
                .flat_map(|mp| mp.original_detections.iter().map(|rd| rd.id))
                .collect();
            valid_trajs.push(Trajectory {
                detection_ids: ids,
                center_id,
            });
        }
    }

    valid_trajs
}




/// Function that takes a traj and tries to orbit fit. Also included the non-detection-prob test. If successful, returns a Link object.
pub fn try_link(
    traj: &Trajectory,
    ic: &InitialCondition,
    trajectory_filter: &Arc<Mutex<TrajectoryFilter>>,
    det_map: &HashMap<u64, Detection>,
    exposures: &[ExposureRow],
    dets_by_expnum: &HashMap<i64, Vec<&Detection>>,
    config: &Config,
    kernel: &SpiceKernel,
    pipeline_metrics: &Arc<PipelineMetrics>,
) -> Option<Link> {

    let ids: HashSet<usize> = traj.detection_ids.iter().copied().collect();
    let mut tried_detection_sets: HashSet<BTreeSet<usize>> = HashSet::new();


    /// Check if this trajectory either:
    /// a. has already been accepted as a link 
    /// b. has failed (this should only be if it fails for the same IC I think)
    if trajectory_filter.lock().unwrap().should_skip(&ids) {
        return None;
    }


    /// Get the detections from in the trajecory, and get the set of all possible triplets
    let cluster_detections: Vec<&Detection> = ids.iter().map(|&id| det_map.get(&(id as u64)).unwrap()).collect();
    let triplets = cluster_detections.clone().into_iter().combinations(3);

    for mut triplet in triplets {
        triplet.sort_by(|a, b| a.epoch.partial_cmp(&b.epoch).unwrap());

        /// Ensure the triplet is from at least 3 different nights
        let nights: HashSet<_> = triplet.iter().map(|d| d.nite).collect();
        if nights.len() != 3 {
            continue;
        }

        /// Ensure the triplet duration is long enough
        let span = triplet[2].epoch - triplet[0].epoch;
        if span < config.min_duration {
            continue;
        }

        /// Convert the detections to SpaceRocks Observation objects for the gauss fit
        /// NOTE: I can make a Hashmap of intids and observation objects so that I don't have to keep recomputing these
        let (o1, o2, o3) = match (
            triplet[0].to_observation(),
            triplet[1].to_observation(),
            triplet[2].to_observation(),
        ) {
            (Some(o1), Some(o2), Some(o3)) => (o1, o2, o3),
            _ => continue,
        };

        let orbits = match gauss(&o1, &o2, &o3, 10.0) {
            Some(v) if v.len() == 1 => v, /// Do I need this condition? Am I safe to reject multiple gauss solutions?
            _ => continue,
        };

        /// Skip very hyperbolic orbits
        let mut orbit = orbits.into_iter().next().unwrap();
        if orbit.e() > 3.0 {
            continue;
        }


        /// Find the mahalanobis residuals for this gauss fit
        let (resids, _, _) = match residuals(&cluster_detections, &mut orbit, kernel) {
            Ok(r) => r,
            Err(_) => continue,
        };

        /// ONLY KEEP detections less than 7 sigma from the gauss fit (I was previously using 5? I think I changed it to 7 for some fake.
        /// Need to nail this down, but 7 works for now, probably just permits a few too many FPs). If there aren't enough detections, skip.
        /// Also, if the detections don't match the IDs in the trajectory, skip (Shouldn't ever happen)
        let dets: Vec<_> = cluster_detections
            .iter()
            .zip(resids.iter())
            .filter(|(_, r)| **r < 7.0)
            .map(|(&d, _)| d.clone())
            .collect();

        let surviving_ids: HashSet<usize> = dets.iter().map(|d| d.intid as usize).collect();
        if dets.len() < config.min_detections || !ids.iter().all(|id| surviving_ids.contains(id)) {
            continue;
        }

        /// Ensure the detections span a long enough time
        let t0 = dets.iter().map(|d| d.epoch).fold(f64::INFINITY, f64::min);
        let t1 = dets.iter().map(|d| d.epoch).fold(f64::NEG_INFINITY, f64::max);
        if (t1 - t0) < config.min_duration {
            continue;
        }

        /// Get the gauss fit parameters into the correct format for the SpaceRocks Levenberg-Marquardt fitting
        let theta0 = {
            let tdb_epoch = orbit.epoch.to_tdb().epoch; /// The LM fitter expects the epoch in TDB
            [
                orbit.position.x,
                orbit.position.y,
                orbit.position.z,
                orbit.velocity.x,
                orbit.velocity.y,
                orbit.velocity.z,
                tdb_epoch,
            ]
        };

        /// Make observations from the remaining detections
        let observations: Vec<Observation> = dets.iter().filter_map(|d| d.to_observation()).collect();
        let obs_refs: Vec<&Observation> = observations.iter().collect();


        let surviving_id_set: BTreeSet<usize> = surviving_ids.iter().copied().collect();
        if !tried_detection_sets.insert(surviving_id_set.clone()) {
            continue;  // skip duplicate set if I arrived there from a different triplet
        }

        /// Setup a simulation object
        let sim0 = match Simulation::giants(&orbit.epoch, "J2000", "SSB", kernel) {
            Ok(s) => s,
            Err(_) => continue,
        };

        /// Fit the orbit with SpaceRocks Levenberg-Marquardt fitter, using the gauss fit as a guess
        let sr_fit = match fit_orbit_lm(&obs_refs, &theta0, sim0) {
            Ok(Some(f)) => f,
            _ => continue,
        };

        /// If ANY of the residuals^2 are > 25 (SR residuals function returns the square of the mahalnobis distance),
        /// the fit fails.
        if sr_fit.residuals.iter().any(|&r| r >= 25.0) {
            continue;
        }

        /// Convert the SpaceRocks fit to a sparse_linking_stuff fit -- Just haven't made the effort to make the SR fit serializable
        /// Check if the chisq/dof is too high or if semi-major axis is non-physical
        let fit = to_sparse_linking_stuff_fit(sr_fit);

        /// NOTE: The semi-major axis cut is disabled for now, because it was a couple of marginal real things. The cost is an increased number of FPs. 
        /// For now, just ensuring that our chisq/dof is not too high. NEED TO CHECK chi2 in SR, im not posiitive it's being calculated correctly.
        // if fit.chisq / fit.dof > 5.0 || fit.rock.a() < 0.0 { 
        if fit.chisq / fit.dof > 5.0 {
            continue;
        }

        // If we make it here, we have a successful orbit fit. Track this in the pipeline metrics
        pipeline_metrics.track_post_orbitfit_trajectory(&surviving_ids);


        /// Ok, we have some set of detections for which an orbit has successfully been fit. Now we do our virtual non-detection stacking a la Pedro.
        let should_accept = if config.non_detection {
            let pass = non_detection_prob_pass(
                &fit.rock,
                dets_by_expnum,
                exposures,
                kernel,
                config,
            );
            pass >= config.prob_threshold
        } else {
            true
        };

        if should_accept {
            trajectory_filter.lock().unwrap().add_accepted(&ids);
            return Some(Link::new(
                0,
                dets.iter().filter_map(|d| d.objid.clone()).collect(),
                fit,
                ic.q,
                ic.e,
                ic.inc,
                ic.true_anomaly,
            ));
        } else {
            trajectory_filter.lock().unwrap().add_failed_minimal(&ids);
            continue;
        }
    }

    None
}
