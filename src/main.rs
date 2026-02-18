use serde_yaml;
use clap::Parser;
use itertools::Itertools;
use rayon::prelude::*;
use std::num::NonZero;

use spacerocks::SpiceKernel;
use spacerocks::Time;

use pointdexter::cli::{Config, Cli};
use pointdexter::detection::Detection;
use pointdexter::io::{load_detections, load_initial_conditions};
use pointdexter::sync::*;
use pointdexter::InitialCondition;
use pointdexter::gauss::gauss;

use std::io::Write;

use plotly::{Plot, Scatter};
use plotly::common::Mode;
use indicatif::ParallelProgressIterator;
use indicatif::ProgressIterator;

use kiddo::{KdTree, SquaredEuclidean};
use kiddo::immutable::float::kdtree::ImmutableKdTree;

use plotly::common::Marker;
use plotly::layout::{Axis, Layout};

fn squared_euclid_to_angle_rad(d2: f64) -> f64 {
    // For unit vectors: dot = 1 - d^2/2
    let dot = 1.0 - 0.5 * d2;
    dot.clamp(-1.0, 1.0).acos()
}

pub fn sync_detections_to_orbit(detections: &Vec<Detection>, ic: &InitialCondition) -> Vec<Option<[f64; 3]>> {
    let mut synced_points: Vec<Option<[f64; 3]>> = vec![None; detections.len()];
    for (j, det) in detections.iter().enumerate() {
        if let Some(synced_pos) = sync_detection_to_orbit(det, ic) {
            synced_points[j] = Some(synced_pos);       
        }
    }
    synced_points
}

fn main() -> Result<(), Box<dyn std::error::Error>> {

    type Tree = ImmutableKdTree<f64, u32, 3, 32>; // item is u32 index into `points

    let args = Cli::try_parse()?;
    let f = std::fs::File::open(args.config)?;
    let config: Config = serde_yaml::from_reader(f)?;

    let mut kernel = SpiceKernel::new();
    kernel.load_spk(format!("{}/sb441-n16.bsp", config.spice_path).as_str())?;
    kernel.load_spk(format!("{}/de440s.bsp", config.spice_path).as_str())?;
    kernel.load_bpc(format!("{}/earth_1962_240827_2124_combined.bpc", config.spice_path).as_str())?;

    // Load detections from catalog. They will automatically be transformed to the desired reference plane. 
    // The observer positions and velocities will be rotated accordingly.
    println!("Loading detections...");
    let mut detections = load_detections(&config.detection_catalog, &config.orbit_reference_plane, &kernel)?;
    println!("Loaded {} detections.", detections.len());

    // Load initial conditions from file.
    println!("Loading initial conditions...");
    let mut ics = load_initial_conditions(&config.initial_conditions_file, &config.ic_type, &config.ic_origin, config.reference_epoch)?;
    println!("Loaded {} initial conditions.", ics.len());

    const ARCSEC_PER_RAD: f64 = (180.0 * 3600.0) / std::f64::consts::PI;
    let eps: f64 = config.epsilon_arcsec / ARCSEC_PER_RAD;
    let cluster_radius_squared = 2.0 * (1.0 - eps.cos());

    // open a file to write cluster info to
    let mut cluster_file = std::fs::File::create("clusters5-oopsallfakes.txt")?;
    // wrap in a mutex for thread safety
    let cluster_file = std::sync::Mutex::new(&mut cluster_file);


    // Sync detections to orbits
    let mut idx = 0;
    for ic in ics[..1000].iter() {
    // ics.par_iter().progress_count(ics.len() as u64).for_each(|ic| {


        // timing: 
        let start = std::time::Instant::now();
        let synced_points = sync_detections_to_orbit(&detections, ic);
        let duration = start.elapsed();
        println!("Synced points for IC {}: {}/{} in {:.2?}", ic.id, synced_points.iter().filter(|p| p.is_some()).count(), synced_points.len(), duration);
        // println!("Synced points for IC {}: {}/{}", ic.id, synced_points.iter().filter(|p| p.is_some()).count(), synced_points.len());

        // let mut points: Vec<[f64; 3]> = Vec::new();
        // let mut ids: Vec<u64> = Vec::new();
        // for (i, p) in synced_points.iter().enumerate() {
        //     if let Some(v) = p {
        //         points.push(*v);
        //         ids.push(i as u64);
        //     }
        // }
        // let kdtree: Tree = ImmutableKdTree::new_from_slice(&points);

        // let mut clusters: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();
        // let mut cluster_count = 0;

        // for (idx, point) in synced_points.iter().enumerate() {

        //     if point.is_none() {
        //         continue;
        //     }
        //     let v = point.unwrap();

        //     let mut cluster = kdtree.within_unsorted::<SquaredEuclidean>(&v, cluster_radius_squared);
        //     for nn in cluster.iter_mut() {
        //         nn.item = ids[nn.item as usize] as u32;
        //     }

            
        //     if cluster.len() < config.min_detections {
        //         continue;
        //     }

        //     let central_detection = &detections[idx];
        //     let central_detection_epoch = central_detection.epoch;
            
        //     let indices: Vec<usize> = cluster.iter().map(|nn| nn.item as usize).collect();
        //     let angles: Vec<f64> = cluster.iter().map(|nn| squared_euclid_to_angle_rad(nn.distance)).collect();
        //     let dts: Vec<f64> = indices.iter().map(|&i| (detections[i].epoch - central_detection_epoch).abs()).collect();

        //     // if angle > (dt / max_dt) * eps, or dt > max_dt, then discard
        //     let max_dt = 15.0;
        //     let mut filtered_indices: Vec<usize> = indices.iter()
        //         .zip(angles.iter())
        //         .zip(dts.iter())
        //         .filter(|&((_i, &angle), &dt)| {
        //             if dt > max_dt {
        //                 return false;
        //             }
        //             // add a 1 arcsec buffer to angle threshold
        //             let angle_thresh = ((dt / max_dt) * eps) + (1.0 / ARCSEC_PER_RAD);
        //             angle <= angle_thresh
        //         })
        //         .map(|((&i, &_angle), &_dt)| i)
        //         .collect();

        //     if filtered_indices.len() < config.min_detections {
        //         continue;
        //     }

        //     filtered_indices.sort_unstable();
        //     if clusters.contains(&filtered_indices) {
        //         continue;
        //     }
        //     clusters.insert(filtered_indices.clone());
            

        //     let cluster_detections: Vec<&Detection> = filtered_indices.iter().map(|&i| &detections[i]).collect();

        //     // find the earliest and latest epochs
        //     let min_epoch = cluster_detections.iter().map(|d| d.epoch).fold(f64::INFINITY, |a, b| a.min(b));
        //     let max_epoch = cluster_detections.iter().map(|d| d.epoch).fold(f64::NEG_INFINITY, |a, b| a.max(b));
        //     let duration = max_epoch - min_epoch;
        //     if duration < config.min_duration {
        //         continue;
        //     }

        //     let unique_nights: std::collections::HashSet<i32> = cluster_detections.iter().map(|d| d.epoch.round() as i32).collect();
        //     if unique_nights.len() < config.min_nites {
        //         continue;
        //     }

        //     cluster_count += 1;            

            // // try to optimize the cluster
            // let problem = ClusterProblem {
            //     detections: Arc::new(cluster_detections.iter().map(|&d| d.clone()).collect()),
            //     ic0: ic.clone(),
            // };
            // // initial parameter vector
            // let p0: [f64; 4] = [ic.r, ic.vr, ic.vo, ic.inc];

            // let step = [0.05 * ic.r, 0.1 * ic.vr.abs(), 0.1 * ic.vo.abs(), 1e-4];
            // let (opt_ic, cost) = match optimize_ic(
            //     Arc::new(cluster_detections.iter().map(|&d| d.clone()).collect()),
            //     ic.clone(),
            //     [ic.r, ic.vr, ic.vo, ic.inc],
            //     step,
            // ) {
            //     Ok(res) => res,
            //     Err(e) => {
            //         eprintln!("Optimization error for IC {}: {}", ic.id, e);
            //         continue;
            //     }
            // };
            // println!("Optimized IC {}: r = {:.3}, vr = {:.3}, vo = {:.3}, inc = {:.3} deg, cost = {:.6}", 
            //     ic.id, opt_ic.r, opt_ic.vr, opt_ic.vo, opt_ic.inc.to_degrees(), cost * ARCSEC_PER_RAD);

            // // println!("Initial cluster radius: {}", cost);
            // // write to file:

            // acquire the mutex lock to write to the file
            // let mut file = cluster_file.lock().unwrap();
            // // write the ic.id, and then a list of detection ids
            // write!(file, "{},{}\n", ic.id, filtered_indices.len()).unwrap();
            // for &j in filtered_indices.iter() {
            //     let detid = detections[j].detid.clone().unwrap();
            //     write!(file, "{},", detid).unwrap();
            // }
            // writeln!(file).unwrap();


            // let initial_synced_points = sync_detections_to_orbit(&cluster_detections.iter().map(|&d| d.clone()).collect(), ic);
            // let mut initial_cluster_vectors: Vec<[f64; 3]> = Vec::new();
            // for p in initial_synced_points.iter() {
            //     if let Some(v) = p {
            //         initial_cluster_vectors.push(*v);
            //     }
            // }

            // let initial_phis = initial_cluster_vectors.iter().map(|v| v[1].atan2(v[0])).collect_vec();
            // let initial_thetas = initial_cluster_vectors.iter().map(|v| (v[2].clamp(-1.0, 1.0)).asin()).collect_vec();

            // let final_synced_points = sync_detections_to_orbit(&cluster_detections.iter().map(|&d| d.clone()).collect(), &opt_ic);
            // let mut cluster_vectors: Vec<[f64; 3]> = Vec::new();
            // for p in final_synced_points.iter() {
            //     if let Some(v) = p {
            //         cluster_vectors.push(*v);
            //     }
            // }
            // let phis = cluster_vectors.iter().map(|v| v[1].atan2(v[0])).collect_vec();
            // let thetas = cluster_vectors.iter().map(|v| (v[2].clamp(-1.0, 1.0)).asin()).collect_vec();

            // let mut plot = Plot::new();

            // let scatter = Scatter::new(phis.clone(), thetas)
            //     .name("Cluster Points")
            //     .mode(Mode::Markers)
            //     .marker(Marker::new().color("#000000"));

            // let scatter2 = Scatter::new(initial_phis, initial_thetas)
            //     .name("Initial Points")
            //     .mode(Mode::Markers)
            //     .marker(Marker::new().color("#ff0000"));



            // plot.add_trace(scatter);
            // plot.add_trace(scatter2);
            
            // let layout = Layout::new()
            //     // Anchor the y-axis to the "x" axis
            //     .y_axis(Axis::new().scale_anchor("x"))
            //     // Optionally, set the ratio between them (default is 1)
            //     // .yaxis(Axis::new().scaleanchor("x").scaleratio(1.0)) 
            //     // You might also want to set an explicit size for a truly "square" plot area
            //     .width(500)
            //     .height(500);

            // plot.set_layout(layout);



            // let filename = format!(
            //     "/Users/kjnapier/Desktop/clusters/cluster_{}_{}.png",
            //     ic.id, cluster_count
            // );

            // plot.write_image(&filename, plotly::ImageFormat::PNG, 600, 600, 1.0);


        }   
        
//    });
    // }
    Ok(())
}

