// src/sparse_linking.rs

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
use kiddo::SquaredEuclidean;

use serde_yaml;
use serde::{Serialize, Deserialize};
use pointdexter::sparse_linking_stuff::initial_condition::InitialCondition;
use pointdexter::sparse_linking_stuff::detection::Detection;
// use Clustering-and-Linking::detection::Detection;
use pointdexter::sparse_linking_stuff::construct_orbit::{construct_orbit, construct_orbit_from_fake};
use pointdexter::sparse_linking_stuff::config::{Config, Cli};
use pointdexter::sparse_linking_stuff::ccd::*;
use pointdexter::sparse_linking_stuff::non_detection_prob::*;
use pointdexter::sparse_linking_stuff::io::ExposureRow;
use pointdexter::sparse_linking_stuff::grid::GridCacheKey;

use spacerocks::StateVector;
// use spacerocks::orbfit::fitter::residuals as residuals_sr;
use bloomfilter::Bloom;
use std::fs::File;
use std::path::PathBuf;
use pointdexter::sparse_linking_stuff::Link;
use spacerocks::nbody::Simulation;
use spacerocks::orbfit::gauss::gauss;
use spacerocks::SpiceKernel;

use spacerocks::orbfit::fitter::{FitResult as FitResultSr};
use pointdexter::sparse_linking_stuff::orbfit::fitter::{FitResult};

use spacerocks::orbfit::fitter::fit_orbit_lm;
use spacerocks::orbfit::fitter::FitResult as SrFit;
use pointdexter::sparse_linking_stuff::orbfit::fitter::FitResult as HsFit;

use spacerocks::observing::observation::ObservationType;
use spacerocks::observing::{Observer, Observation};
use spacerocks::SpaceRock;
use itertools::Itertools;
use spacerocks::constants::MU_BARY;

use std::hash::{Hash, Hasher};
use pointdexter::sparse_linking_stuff::grid::{Cell, IcBoundsKey, refine_with_grid_cache_sparse};
use pointdexter::sparse_linking_stuff::trajectory::*;



/// Adding in a struct to see how the data is refined through the pipeline
#[derive(Debug)]
pub struct PipelineMetrics {
    cluster_count: AtomicUsize,
    initial_trajectory_count: AtomicUsize,
    post_epsilon_trajectory_count: AtomicUsize,
    post_orbitfit_trajectory_count: AtomicUsize,
    final_trajectory_count: AtomicUsize,
}

impl PipelineMetrics {
    pub fn new() -> Self {
        Self {
            cluster_count: AtomicUsize::new(0),
            initial_trajectory_count: AtomicUsize::new(0),
            post_epsilon_trajectory_count: AtomicUsize::new(0),
            post_orbitfit_trajectory_count: AtomicUsize::new(0),
            final_trajectory_count: AtomicUsize::new(0),
        }
    }
    
    // Just increment counters - no deduplication
    pub fn track_cluster(&self, _detection_ids: &[usize]) {
        self.cluster_count.fetch_add(1, Ordering::Relaxed);
    }
    
    pub fn track_initial_trajectory(&self, _detection_ids: &[usize]) {
        self.initial_trajectory_count.fetch_add(1, Ordering::Relaxed);
    }
    
    pub fn track_post_epsilon_trajectory(&self, _detection_ids: &[usize]) {
        self.post_epsilon_trajectory_count.fetch_add(1, Ordering::Relaxed);
    }

    pub fn track_post_orbitfit_trajectory(&self, _detection_ids: &HashSet<usize>) {
        self.post_orbitfit_trajectory_count.fetch_add(1, Ordering::Relaxed);
    }
    
    pub fn track_final_trajectory(&self, _detection_ids: &HashSet<usize>) {
        self.final_trajectory_count.fetch_add(1, Ordering::Relaxed);
    }
    
    pub fn print_stage_summary(&self, stage_name: &str) {
        println!("\n{}", stage_name);
        println!("1. Total clusters: {}", self.cluster_count.load(Ordering::Relaxed));
        println!("2. Total initial trajectories: {}", self.initial_trajectory_count.load(Ordering::Relaxed));
        println!("3. Total post-epsilon cascade trajectories (1\"): {}", self.post_epsilon_trajectory_count.load(Ordering::Relaxed));
        println!("4. Total post-orbit-fit trajectories: {}", self.post_orbitfit_trajectory_count.load(Ordering::Relaxed));
        println!("5. Total final trajectories (post orbit-fit & non-det prob): {}\n", self.final_trajectory_count.load(Ordering::Relaxed));
    }
}




/// I need to move this somewhere else. The reason for this is because FitResult from SR isn't serializable. Look into modifying this in
/// spacerocks, or just move this somewhere else.
fn to_sparse_linking_stuff_fit(sr: SrFit) -> HsFit {
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

fn index_detections_by_expnum<'a>(
    dets_all: &'a [Detection]
) -> HashMap<i64, Vec<&'a Detection>> {
    let mut map: HashMap<i64, Vec<&Detection>> = HashMap::new();
    for d in dets_all {
        map.entry(d.expnum).or_default().push(d);
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

        let observation = detection.to_observation();
        let observed_parameters = DVector::from_vec(vec![observation.ra(), observation.dec()]);
          

        // Get the rock to the epoch of the detection
        sim.integrate(&detection.epoch);
        let mut rock = sim.get_particle("rock")?.clone();

        // calculate the model observations
        let astro = rock.observe(&detection.observer)?;
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


/// Function that takes a traj and tries to orbit fit. Also included the non-detection-prob test. If successful, returns a Link object.
fn try_link(
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
        triplet.sort_by(|a, b| a.epoch.epoch.partial_cmp(&b.epoch.epoch).unwrap());

        /// Ensure the triplet is from at least 3 different nights
        let nights: HashSet<_> = triplet.iter().map(|d| d.nite).collect();
        if nights.len() != 3 {
            continue;
        }

        /// Ensure the triplet duration is long enough
        let span = triplet[2].epoch.epoch - triplet[0].epoch.epoch;
        if span < config.min_duration {
            continue;
        }

        /// Convert the detections to SpaceRocks Observation objects for the gauss fit
        /// NOTE: I can make a Hashmap of intids and observation objects so that I don't have to keep recomputing these
        let (o1, o2, o3) = (
            &triplet[0].to_observation(),
            &triplet[1].to_observation(),
            &triplet[2].to_observation(),
        );

        let orbits = match gauss(o1, o2, o3, 10.0) {
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
        let t0 = dets.iter().map(|d| d.epoch.epoch).fold(f64::INFINITY, f64::min);
        let t1 = dets.iter().map(|d| d.epoch.epoch).fold(f64::NEG_INFINITY, f64::max);
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
        let observations: Vec<Observation> = dets.iter().map(|d| d.to_observation()).collect();
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
        let pass = non_detection_prob_pass(
            &fit.rock,
            // &dets,
            // dets_all,
            // epoch_expnum,
            // expnum_pos,
            dets_by_expnum,
            exposures, 
            // // The following 3 params comes from my detection probability model
            // 0.808, 
            // 24.435,
            // 0.609,
            // 0.1, // This is the minimal cdf probability for a link to be accepted
            kernel,
            config,
        );


        if pass >= config.prob_threshold {
            trajectory_filter.lock().unwrap().add_accepted(&ids);
            return Some(Link::new(
                0,
                dets.iter().map(|d| d.objid.clone()).collect(),
                fit,
                ic.q,
                ic.e,
                ic.psi,
                ic.true_anomaly,
            ));

        } else {
            trajectory_filter.lock().unwrap().add_failed_minimal(&ids);
            continue;
        }
    }

    None
}

fn run_ic_cluster_and_trajectory_search(
    ic: &InitialCondition,
    detections: &[Detection],
    det_map: &HashMap<u64, Detection>,
    epsilon_arcsec: f64,
    config: &Config,
    arcsec_per_rad: f64,
    pipeline_metrics: &Arc<PipelineMetrics>,
    // fakes: &[SpaceRock]
) -> Vec<Trajectory> {

    let mut valid_trajs = Vec::new();
    let eps = epsilon_arcsec / arcsec_per_rad;
    let cluster_radius = 2.0 * (1.0 - eps.cos());

    let mut tree = KdTree::new();
    let mut point_map = HashMap::new();
    let mut center_list = Vec::new();

    for det in detections {
        // If the det.orbit_id != ic.id, skip it. Note: THIS IS JUST FOR MAKING PURE CLUSTERS/TRAJS WITH FAKE CATALOGS
        // if let Some(orbit_id) = det.orbit_id {
        //     if orbit_id as usize != ic.id {
        //         continue;
        //     }
        // }


        /// Iterate over detections, attempting to seed a cluster with each detection
        if let Some((pt, theta)) = construct_orbit(det, ic, Some(false)) {
            if pt.iter().any(|v| !v.is_finite()) {
                println!("Skipping invalid KD point: {:?} from det {}", pt, det.intid);
                continue;
            }

            /// Add all the pointing vectors to the KdTree
            tree.add(&pt, det.intid as u64);
            point_map.insert(det.intid as u64, pt);


            // // Only use detections with fakeids as cluster centers. NOTE: ONLY UNCOMMENT FOR TESTING WITH FAKE CATALOGS
            // if let Some(fakeid) = &det.fakeid {
            //     center_list.push((det.intid as u64, pt, theta));
            // }
            center_list.push((det.intid as u64, pt, theta));
        }
    }


    // // Ensure that we are not building near-duplicate clusters from the same IC --> Deprecated for now
    // let used_centers = Mutex::new(KdTree::new());
    // let one_arcsec = 1.0 / arcsec_per_rad;
    // let min_center_dist = 2.0 * (1.0 - one_arcsec.cos());

    for &(center_id, center_pt, theta) in &center_list {

        // {
        //     let uc = used_centers.lock().unwrap();
        //     if uc.size() > 0 {
        //         if let Some(n) = uc.nearest_n::<SquaredEuclidean>(&center_pt, 1).get(0) {
        //             if n.distance < min_center_dist {
        //                 continue;
        //             }
        //         }
        //     }
        // }


        /// Get all cluster neighbors that are witin +- config.t_max/2.0 of the center detection's epoch, AND that are within config.epsilon_arcsec * dt / config.t_max

        // For testing, make sure clusters are pure

        // let center_fakeid = det_map.get(&center_id).unwrap().fakeid.clone().unwrap_or_default();
        // let all_neighbors = tree.within_unsorted::<SquaredEuclidean>(&center_pt, cluster_radius);
        // let center_epoch = det_map.get(&center_id).unwrap().epoch.epoch;
        // let cluster = tree.within_unsorted::<SquaredEuclidean>(&center_pt, cluster_radius);

        let center_epoch = det_map.get(&center_id).unwrap().epoch.epoch;
        let all_neighbors = tree.within_unsorted::<SquaredEuclidean>(&center_pt, cluster_radius);

        let neighbors: Vec<_> = all_neighbors
            .into_iter()
            .filter(|nb| {
                let det = det_map.get(&nb.item).unwrap();
                let dt = (det.epoch.epoch - center_epoch).abs();

                // Reject if too far in time
                if dt > config.t_max / 2.0 + 1.0 {
                    return false;
                }

                let pt2 = *point_map.get(&nb.item).unwrap();
                let sep_arcsec = angular_separation(center_pt, pt2)
                    .iter()
                    .map(|x| x.powi(2))
                    .sum::<f64>()
                    .sqrt();

                // Case 1: very close neighbor — always accept
                if sep_arcsec < 0.1 * config.epsilon_arcsec {
                    return true;
                }

                // Case 2: consistent motion — velocity-based acceptance
                sep_arcsec <  config.epsilon_arcsec * dt / (config.t_max / 2.0)
            })
            .collect();


        /// Skip if we don't have enough detections in the cluster
        if neighbors.len() < config.min_detections || neighbors.len() > config.max_detections {
            continue;
        }

        /// Ensure we have enough nights, and that the duration is long enough
        let mut nights = HashSet::new();
        for nb in &neighbors {
            let d = det_map.get(&nb.item).unwrap();
            
            nights.insert(d.nite);
        }

        if nights.len() < config.min_nites {
            continue;
        }

        /// count the number of unique expnums in det.expnum and ensure we have at least 6
        let mut expnums = HashSet::new();
        for nb in &neighbors {
            let d = det_map.get(&nb.item).unwrap();
            expnums.insert(d.expnum);
        }

        if expnums.len() < config.min_detections {
            continue;
        }   

        let times: Vec<f64> = neighbors.iter().map(|nb| det_map.get(&nb.item).unwrap().epoch.epoch).collect();
        let min_t = *times.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        let max_t = *times.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        if (max_t - min_t) < config.min_duration {
            continue;
        }

        /// Track this unique cluster in the pipeline metrics
        let cluster_ids: Vec<usize> = neighbors.iter().map(|nb| nb.item as usize).collect();
        pipeline_metrics.track_cluster(&cluster_ids);

        /// Make the raw detection vector for the trajectory search, which has the separation in d_phi, d_theta space
        let raw_vec: Vec<RawDetection> = neighbors
            .iter()
            .map(|nb| {
                let det = det_map.get(&nb.item).unwrap();
                let pt2 = *point_map.get(&nb.item).unwrap();
                let [dphi, dtheta] = angular_separation(center_pt, pt2);
                let mag = det.mag.clone().and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
                RawDetection {
                    id: nb.item as usize,
                    delta_phi: dphi,
                    delta_theta: dtheta,
                    epoch: det.epoch.epoch,
                    magnitude: mag,
                    night: det.nite,
                    expnum: det.expnum,
                }
            })
            .collect();

        /// Make my meta points and run the trajectory search. Skip if I don't find any valid trajectories
        /// Values are currently hardcoded here, need to make configurable later. 
        // 1.0 and 0.2 are the spatial and temporal separations for MPs
        // 0.6 and 45.0 are the max absolute velocity difference and max angle for the trajectory search

        let meta_pts = create_meta_points(&raw_vec, 1.0, 0.2);
        let trajs = graph_based_trajectory_search(&meta_pts, 0.6, 45., &config);

        // Pull out the individual detections from the MPs and save these as Trajectory objects
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

    /// Return all valid trajectories found in this IC
    valid_trajs
}



fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start = Instant::now();
    let args = Cli::parse();
    let f = std::fs::File::open(&args.config)?;
    let config: Config = serde_yaml::from_reader(f)?;
    let mu = MU_BARY;


    // Shared pipeline metrics across threads
    let pipeline_metrics = Arc::new(PipelineMetrics::new());


    // Set up constants and epsilon ladder
    const ARCSEC_PER_RAD: f64 = 3600.0 * 180.0 / std::f64::consts::PI;
    let epsilon_values = vec![15.0, 5.0, 3.0, 1.0]; // CURRENT USE
    let t_bounds = (0.0, config.t_max);

    let mut kernel = SpiceKernel::new();
    kernel.load_spk(&format!("{}/sb441-n16.bsp", config.spice_path))?;
    kernel.load_spk(&format!("{}/de440s.bsp", config.spice_path))?;
    kernel.load_bpc(&format!("{}/earth_1962_240827_2124_combined.bpc", config.spice_path))?;


    // Read in my catalogs
    let ics: Vec<InitialCondition> = pointdexter::sparse_linking_stuff::io::read_initial_conditions_bounded(&config.initial_conditions, MU_BARY)?;
    let detections: Vec<Detection> = pointdexter::sparse_linking_stuff::io::catalog_to_detections(&config.detection_catalog, &kernel)?;
    let all_x3_exps: Vec<ExposureRow> = pointdexter::sparse_linking_stuff::io::read_exposure_metadata(&config.all_x3_exps_file)?;

    let dets_by_expnum = index_detections_by_expnum(&detections);

    // Make a HashMap of detections by intid
    let mut det_map: HashMap<u64, Detection> = HashMap::new();
    for d in &detections {
        det_map.insert(d.intid as u64, d.clone());
    }

    let trajectory_filter = Arc::new(Mutex::new(TrajectoryFilter::new()));

    // Comment out when not orbit-fitting
    let mut final_links = HashSet::new();
    let mut all_saved_trajectories = Vec::new(); // Store all saved trajectories

    let batch_size = 5_000; // Size of this can be adjusted based on memory and performance, just don't want to have too many trajectories in memory at once
    let num_batches = (ics.len() + batch_size - 1) / batch_size;
    println!("Processing {} batches of size {}", num_batches, batch_size);

    for (batch_index, ic_batch) in ics.chunks(batch_size).enumerate() {
        println!("Processing batch {} of {}", batch_index + 1, num_batches);
    
        // Shared bucket for trajectories found in this batch
        let shared_traj_bucket: Arc<Mutex<Vec<(Trajectory, InitialCondition)>>> = Arc::new(Mutex::new(Vec::new()));

        // Iterate over initial conditions in parallel, collecting trajectories
        ic_batch.par_iter().progress_count(ic_batch.len() as u64).for_each(|ic| {
            let mut valid_trajs = run_ic_cluster_and_trajectory_search(
                ic, &detections, &det_map, config.epsilon_arcsec, &config, ARCSEC_PER_RAD, &pipeline_metrics,
            );

            for traj in valid_trajs.drain(..) {
                shared_traj_bucket.lock().unwrap().push((traj, ic.clone()));
            }
        });

        let mut current_trajs: Vec<(Trajectory, InitialCondition)> = shared_traj_bucket.lock().unwrap().drain(..).collect();
        
        /// Comment out for analyzing initial trajectories
        let grid_cache = Arc::new(Mutex::new(HashMap::new()));


        // Loop through current_trajs to see if any trajs are in trajectory_filter
        current_trajs.retain(|(traj, _)| {
            let ids: HashSet<usize> = traj.detection_ids.iter().copied().collect();
            !trajectory_filter.lock().unwrap().should_skip(&ids)
        });


        for (traj, _ic) in &current_trajs {
            pipeline_metrics.track_initial_trajectory(&traj.detection_ids);
        }

        println!("Initial trajectory count: {}", current_trajs.len());

        /// Remove for analyzing initial trajectories
        for &epsilon in &epsilon_values {
            let eps_rad = epsilon / ARCSEC_PER_RAD;
            let cluster_radius = 2.0 * (1.0 - eps_rad.cos());

            current_trajs = refine_with_grid_cache_sparse(
                &current_trajs,
                &grid_cache,
                &det_map,
                cluster_radius,
                epsilon,
                t_bounds,
                2, // The number of surviving sub-cells to keep. I.e, keep the 2 'tightest' clusters of each trajectory
                // &fakes,
            );

            

            println!("Remaining trajectories after epsilon = {} arcsec: {}", epsilon, current_trajs.len());
        }

        println!("Number of current trajectories: {}", current_trajs.len());

        let mut seen: HashSet<BTreeSet<usize>> = HashSet::new();
        let mut deduped_trajs = Vec::new();

        // Deduplicate trajectories based on their detection IDs before orbit fitting
        for (traj, ic) in current_trajs.into_iter() {
            // Collect detection IDs into a sorted set
            let id_set: BTreeSet<usize> = traj.detection_ids.iter().copied().collect();

            // Insert into HashSet; if it's new, keep the trajectory
            if seen.insert(id_set.clone()) {
                deduped_trajs.push((traj, ic));
            }
        }
        

        println!("Number of deduplicated trajectories: {}", deduped_trajs.len());

        for (traj, _ic) in &deduped_trajs {
            pipeline_metrics.track_post_epsilon_trajectory(&traj.detection_ids);
        }
      

        deduped_trajs.sort_by_key(|(traj, _)| std::cmp::Reverse(traj.detection_ids.len()));

        // Convert deduplicated trajectories to saved format
        for (traj, ic) in &deduped_trajs {
            let dets: Vec<Detection> = traj.detection_ids.iter()
                .filter_map(|&id| det_map.get(&(id as u64)).cloned())
                .collect();
            
            let detection_objids: Vec<String> = dets.iter()
                .map(|d| d.objid.clone())
                .collect();

            let saved_traj = SavedTrajectory {
                detection_objids,
                ic_params: [ic.r0, ic.vr, ic.vo, ic.psi],
            };

            all_saved_trajectories.push(saved_traj);
        }

        /// Comment out when not orbit-fitting
        let collected_links: HashSet<Link> = deduped_trajs
            .par_iter()
            .progress_count(deduped_trajs.len() as u64)
            .filter_map(|(traj, ic)| {
                // Create new SpiceKernel per trajectory so that I can use in parallel
                let mut local_kernel = SpiceKernel::new();
                if local_kernel.load_spk(&format!("{}/sb441-n16.bsp", config.spice_path)).is_err() { return None; }
                if local_kernel.load_spk(&format!("{}/de440s.bsp", config.spice_path)).is_err() { return None; }
                if local_kernel.load_bpc(&format!("{}/earth_1962_240827_2124_combined.bpc", config.spice_path)).is_err() { return None; }

                try_link(
                    traj,
                    ic,
                    &trajectory_filter, 
                    &det_map,
                    &all_x3_exps,
                    &dets_by_expnum,
                    &config,
                    &local_kernel,
                    &pipeline_metrics,
                )
            })
            .collect();
        
        final_links.extend(collected_links);
    }


    for link in &final_links {
        let ids: HashSet<usize> = link.objids.iter()
            .filter_map(|objid| {
                detections.iter()
                    .find(|d| &d.objid == objid)
                    .map(|d| d.intid as usize)
            })
            .collect();
        pipeline_metrics.track_final_trajectory(&ids);
    }
    
    pipeline_metrics.print_stage_summary("Final Results (Post Orbit-Fit & Non-Det Prob)");

    // deduplicate all saved trajectories
    all_saved_trajectories.sort_by_key(|t| t.detection_objids.len());
    all_saved_trajectories.dedup_by_key(|t| t.detection_objids.clone());


    // Save results
    println!("Total number of saved trajectories: {}", all_saved_trajectories.len());
    let elapsed = start.elapsed();
    println!("Elapsed time: {:?}", elapsed);

    // // Comment out when not orbit-fitting
    let mut links_path = PathBuf::from(&config.output_path);
    links_path.push(&config.links_filename);
    let file = File::create(&links_path)?;
    let encoder = GzEncoder::new(file, Compression::default());
    let mut writer = BufWriter::new(encoder);
    serde_json::to_writer(&mut writer, &final_links)?;

    println!("Wrote {} links to {:?}", final_links.len(), links_path);

    // Save the deduplicated trajectories
    let mut trajs_path = PathBuf::from(&config.output_path);
    trajs_path.push(&config.trajectories_filename);
    let traj_file = File::create(&trajs_path)?;
    let traj_encoder = GzEncoder::new(traj_file, Compression::default());
    let mut traj_writer = BufWriter::new(traj_encoder);
    serde_json::to_writer(&mut traj_writer, &all_saved_trajectories)?;
    println!("Wrote {} deduplicated trajectories to {:?}", all_saved_trajectories.len(), trajs_path);
    Ok(())
}



