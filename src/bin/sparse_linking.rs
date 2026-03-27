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
use kiddo::immutable::float::kdtree::ImmutableKdTree;
use kiddo::SquaredEuclidean;

use serde_yaml;
use serde::{Serialize, Deserialize};
use pointdexter::initial_condition::InitialCondition;
use pointdexter::detection::Detection;
// use Clustering-and-Linking::detection::Detection;
use pointdexter::sync::sync_detection_to_orbit as construct_orbit;
use pointdexter::sparse_linking_stuff::config::{Config, Cli};
use pointdexter::sparse_linking_stuff::ccd::*;
use pointdexter::sparse_linking_stuff::non_detection_prob::*;
use pointdexter::sparse_linking_stuff::utils::ExposureRow;
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
use spacerocks::Time;

use std::hash::{Hash, Hasher};
use pointdexter::sparse_linking_stuff::grid::{Cell, IcBoundsKey, refine_with_grid_cache_sparse};
use pointdexter::sparse_linking_stuff::trajectory::*;
use pointdexter::sparse_linking_stuff::utils::*;
use pointdexter::sparse_linking_stuff::pipeline_metrics::PipelineMetrics;
use pointdexter::sparse_linking_stuff::cluster_and_link::*;




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
    let epsilon_values = vec![45.0, 15.0, 5.0, 3.0, 1.0]; // CURRENT USE
    let t_bounds = (0.0, config.t_max);

    let mut kernel = SpiceKernel::new();
    kernel.load_spk(&format!("{}/sb441-n16.bsp", config.spice_path))?;
    kernel.load_spk(&format!("{}/de440s.bsp", config.spice_path))?;
    kernel.load_bpc(&format!("{}/earth_1962_240827_2124_combined.bpc", config.spice_path))?;

    let ic_path = config.initial_conditions.clone();
    let det_path = config.detection_catalog.clone();
    
    // Read in my catalogs
    let ics: Vec<InitialCondition> = pointdexter::io::load_ics::load_initial_conditions(&ic_path, "keplerian", "SSB", config.t_ref)?;
    let detections: Vec<Detection> = pointdexter::io::load_detections::load_detections(&det_path, &"J2000", &kernel)?;
    let all_x3_exps: Vec<ExposureRow> = read_exposure_metadata(&config.all_x3_exps_file)?;

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

    let batch_size = config.batch_size; // Size of this can be adjusted based on memory and performance, just don't want to have too many trajectories in memory at once
    let num_batches = (ics.len() + batch_size - 1) / batch_size;
    println!("Processing {} batches of size {}", num_batches, batch_size);

    for (batch_index, ic_batch) in ics.chunks(batch_size).enumerate() {
        println!("Processing batch {} of {}", batch_index + 1, num_batches);
    
        // Shared bucket for trajectories found in this batch
        let shared_traj_bucket: Arc<Mutex<Vec<(Trajectory, InitialCondition)>>> = Arc::new(Mutex::new(Vec::new()));

        // Iterate over initial conditions in parallel, collecting trajectories
        ic_batch.par_iter().progress_count(ic_batch.len() as u64).for_each(|ic| {
            let mut valid_trajs = run_ic_cluster_and_trajectory_search(
                ic, &detections, &det_map, config.epsilon_arcsec, &config, &pipeline_metrics,
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


        //  // Print the IC IDs of all curent trajectories for debugging
        // for traj in &current_trajs {
        //     println!("Trajectory with center ID: {}", traj.1.id);
        // }
           

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
                .filter_map(|d| d.objid.clone())
                .collect();

            let saved_traj = SavedTrajectory {
                detection_objids,
                ic_params: [ic.r, ic.vr, ic.vo, ic.inc],
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
                    .find(|d| d.objid.as_deref() == Some(objid))
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



