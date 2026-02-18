use std::f64::consts::PI;
use std::fs::File;
use std::time::Instant;
use std::collections::{HashMap, HashSet};

use serde_yaml;
use clap::Parser;

use cdshealpix::nested::map::skymap::CountMap;
use cdshealpix::nested::moc::HasMaxDepth;
use cdshealpix::nested::moc::op::or;
use cdshealpix::nested::{get, Layer};
use cdshealpix::{TRANSITION_LATITUDE};
use cdshealpix::{nside};
use cdshealpix::nested;

use rayon::prelude::*;

use spacerocks::{SpiceKernel};

use kiddo::SquaredEuclidean;
use kiddo::immutable::float::kdtree::ImmutableKdTree;

use pointdexter::cli::{Config, Cli};
use pointdexter::utils::{unit_vector_to_lonlat, lonlat_to_unit_vector, ecliptic_to_equatorial, ecliptic_to_equatorial_vector, cull_cluster, gather_cluster_data, iterative_reject, find_non_strict_peaks, sum_hp_map_neighbors_v2};
use pointdexter::detection::{Detection, ReferencePlane};
use pointdexter::sync::sync_detection_to_orbit;
use pointdexter::load_detections;
use pointdexter::load_ics::load_initial_conditions;
use pointdexter::InitialCondition;

const MAX_DT: f64 = 4000.0; // days
const MIN_UNIQUE_TIMES: usize = 12;
const MIN_POINTS_PER_CELL: usize = 4;
const MIN_NIGHTS: usize = 4;

const ARCSEC_PER_RAD: f64 = 3600.0 * 180.0 / PI;

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

    type Tree = ImmutableKdTree<f64, u32, 3, 32>; // item is u32 index into `points` array, 3D points, leaf capacity 32

    let depth = 15_u8;
    let nside = nside(depth) as u64;
    let npix = 12 * nside * nside;
    println!("Total HEALPix cells at depth {}: {}", depth, npix);
    println!("Nside at depth {}: {} {}", depth, nside, (4.0*PI/(npix as f64)).sqrt()*(180.0/PI)*60.0*60.0);
    let nested: &Layer = get(depth);

    // Calculate cluster angles.
    const EPS: f64 = 7.0 / ARCSEC_PER_RAD; 
    let cluster_r2 = 2.0 * (1.0 - EPS.cos());
    let obliq = 84381.448/3600.0;

    println!("Cluster radius (radians): {}", EPS*ARCSEC_PER_RAD);
    
    // We need to make this more flexible.
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

   
    let start_time = Instant::now();
    let n_ics = ics.len();

    let mut ii = 0;

    ics[..n_ics].par_iter().enumerate().for_each(|(ii, ic)| {
        type Tree = ImmutableKdTree<f64, i32, 3, 32>; // item is u32 index into `points` array, 3D points, leaf capacity 32

        let mut out = String::new();
        
        let nested: &Layer = get(depth);
        
        out.push_str(&format!("# Processing IC {} of {}\n", ii + 1, n_ics));
            
        let mut points: Vec<[f64; 3]> = Vec::with_capacity(detections.len());

        let mut hp_map: HashMap<u64, Vec<i32>> = HashMap::new();
        let mut hp_map_counts: HashMap<u64, usize> = HashMap::new();

            
        let synced_points = sync_detections_to_orbit(&detections, ic);
        let mut points: Vec<[f64; 3]> = Vec::new();
        let mut intids: Vec<u64> = Vec::new();
        for (i, p) in synced_points.iter().enumerate() {
            if let Some(v) = p {
                let (lon, lat) = unit_vector_to_lonlat(v[0], v[1], v[2]);
                let idx = nested.hash(lon, lat);
                *hp_map_counts.entry(idx as u64).or_insert(0) += 1;
                points.push(*v);
                intids.push(i as u64);
            }
        }

        // let nside = nside(depth) as u64;
        // let npix = 12 * nside * nside;

        let radius = 2.0*(4.0*PI/(npix as f64)).sqrt();
        let mut hp_sum_map_v2 = sum_hp_map_neighbors_v2(&hp_map_counts, nested, radius);
        let peaks = find_non_strict_peaks(&hp_sum_map_v2, depth, MIN_UNIQUE_TIMES);

        // Keep track of the mapping from intid → index in points array for efficient lookup after neighbor queries.
        let mut intid_idx_hash: HashMap<i32, usize> = HashMap::new();
        for (i, intid) in intids.iter().enumerate() {
            intid_idx_hash.insert(*intid as i32, i);
        }
        
        // Build a global 3D Kd-tree for all points, to enable fast neighbor queries.
        let tree: Tree = ImmutableKdTree::new_from_slice(&points);
            
        let mut seen: HashSet<Vec<u64>> = HashSet::new();
        let mut not_skipped = 0;
        for peak in peaks.iter() {
            let cell = peak.0;
            let cell_count = peak.1;
            if cell_count < MIN_UNIQUE_TIMES {
                continue;
            }
            not_skipped += 1;
                    
            // Get the center of the cell in (lon, lat), convert to (x, y, z) unit vector, 
            // and query the global tree for neighbors within cluster_r.
            let (lon, lat) = nested.center(cell);
            let (cell_x, cell_y, cell_z) = lonlat_to_unit_vector(lon, lat);
            let cell_point = [cell_x, cell_y, cell_z];

            let mut cluster = tree.within::<SquaredEuclidean>(&cell_point, cluster_r2);

            let res = nested.kth_neighbourhood(cell, 4u32);

            // Hash the sorted intids and skip if we've already seen this cluster.
            let mut intids_key: Vec<u64> = cluster
                .iter()
                .map(|nn| intids[nn.item as usize])
                .collect();
            intids_key.sort_unstable();

            if !seen.insert(intids_key.clone()) {
                continue;
            }

            for nn in cluster.iter_mut() {
                nn.item = intids[nn.item as usize] as i32;
            }

            if cluster.len() < MIN_UNIQUE_TIMES {
                continue;
            }
                                
            // Cull the cluster to remove outliers and construct the Cluster struct with local coordinates.
            let clust = cull_cluster(&detections, &intid_idx_hash, &points, cell_point, &cluster, ic.epoch, MAX_DT, MIN_UNIQUE_TIMES, MIN_NIGHTS);
            if clust.is_none() {
                continue;
            }
            let clust = clust.unwrap();

            let (t, theta_x, theta_y, ast_ucty, mxe, mye, xe, ye) = gather_cluster_data(&clust, &detections, ic.epoch);
            
            let fit_result = iterative_reject(
                &t, &mxe, &mye, &theta_x, &theta_y, &ast_ucty, &ast_ucty,
                5.0,   // sigma threshold
                MIN_POINTS_PER_CELL,  // stop if fewer than a threhold number of points remain
                MIN_UNIQUE_TIMES,  // stop if fewer than a threhold number of unique observation times remain
                MIN_NIGHTS, // stop if fewer than a threshold number of unique nights remain
            );
            
            let fit_final = match fit_result {
                Ok(fit_final) => fit_final,
                Err(e) => {
                    //println!("Fit failed: {}", e);
                    continue;
                }
            };
            let (fit, kept_indices, rejected_indices) = fit_final;
            let params = fit.params;

            out.push_str(&format!("# Fit successful for cell {}: {:?}\n", cell, params));
                            
            let rejected_intids: Vec<usize> = rejected_indices
                .iter()
                .map(|&idx| clust.local_intids[idx] as usize)
                .collect();
            
            out.push_str("# intid     detid         t(days)     theta_x(\")  theta_y(\")  x_model(\")  y_model(\")    sig_x     sig_y\n");
            let mut kept_count = 0;
            for i in 0..t.len() {
                if rejected_indices.contains(&i) {
                    continue;
                }
                kept_count += 1;
                let x_model = params[0] + params[1] * t[i] + params[2] * mxe[i];
                let y_model = params[3] + params[4] * t[i] + params[2] * mye[i];
                let rx = theta_x[i] - x_model;
                let pull_x = rx / ast_ucty[i];
                let ry = theta_y[i] - y_model;
                let pull_y = ry / ast_ucty[i];
                let pull = pull_x.abs().max(pull_y.abs());
                let intid = clust.local_intids[i];
                let detid = &detections[intid as usize].detid.clone().unwrap_or("".to_string());
                out.push_str(&format!("{:8} {:12} {:12.6}   {:8.4}     {:8.4}    {:8.4}   {:8.4}.  {:8.4}.  {:8.4}\n", clust.local_intids[i], detid, t[i], theta_x[i]*206265., theta_y[i]*206265., x_model*206265., y_model*206265., pull_x, pull_y));
            }
            out.push_str(&format!("# Kept {} detections in fit, IC {}\n", kept_count, ii+1));
            out.push_str("# detid         epoch        RA(deg)   sig_x(\")   Dec(deg)  sig_y(\")  obs_x        obs_y           obs_z      obscode  mag\n");
            for i in 0..t.len() {
                if rejected_indices.contains(&i) {
                    continue;
                }
                
                let intid = clust.local_intids[i];
                let detid = &detections[intid as usize].detid.clone().unwrap_or("".to_string());
                let objid = &detections[intid as usize].objid.clone().unwrap_or("".to_string());
                let filt = &detections[intid as usize].filter.clone().unwrap_or("".to_string());
                let rho_hat = detections[clust.local_intids[i] as usize].rho_hat; 
                let lam = rho_hat[1].atan2(rho_hat[0]);
                let beta = rho_hat[2].asin();
                // These are ecliptic latitude and longitude at this point.
                // Convert back to equatorial for output.
                let (ra, dec) = ecliptic_to_equatorial(lam, beta, obliq * std::f64::consts::PI / 180.0);

                let ra_deg = ra * 180.0 / PI;
                let dec_deg: f64 = dec * 180.0 / PI;
                let epoch = detections[clust.local_intids[i] as usize].epoch;
                let obscode = "X05"; // TODO.
                let mag = detections[clust.local_intids[i] as usize].mag.unwrap_or(0.0);
                let mag_sig = detections[clust.local_intids[i] as usize].mag_ucty.unwrap_or(0.0);
                let obs_pos = &detections[clust.local_intids[i] as usize].observer_position;
                let obs_pos2 = ecliptic_to_equatorial_vector(*obs_pos, obliq * std::f64::consts::PI / 180.0);

                out.push_str(&format!("ic_00 {:12} {:12.6} {:10.6}   {:6.3}   {:10.6}   {:6.3}  {:13.10}  {:13.10}  {:13.10}   {}   {:.2} {:.2} {} {} \n", 
                    detid, epoch, ra_deg, ast_ucty[i]*206265., dec_deg, ast_ucty[i]*206265., obs_pos2[0], obs_pos2[1], obs_pos2[2], obscode, mag, mag_sig, filt, objid));
            }
                

        }
        out.push_str(&format!("# Finished processing IC {} of {}: {} cells, {} of which had > {} points.\n", ii+1, n_ics, peaks.len(), not_skipped, MIN_UNIQUE_TIMES));
        print!("{}", out);
    });
                
    
    Ok(())
}












// fn run_ic_clustering(detections: &Vec<Detection>, ic: &InitialCondition, depth: u8, cluster_r2: f64, ii: usize, n_ics: usize)-> String {

//     type Tree = ImmutableKdTree<f64, i32, 3, 32>; // item is u32 index into `points` array, 3D points, leaf capacity 32

//     let mut out = String::new();
    
//     let nested: &Layer = get(depth);
    
//     out.push_str(&format!("# Processing IC {} of {}\n", ii + 1, n_ics));
        
//     let mut points: Vec<[f64; 3]> = Vec::with_capacity(detections.len());

//     let mut hp_map: HashMap<u64, Vec<i32>> = HashMap::new();
//     let mut hp_map_counts: HashMap<u64, usize> = HashMap::new();

        
//     let synced_points = sync_detections_to_orbit(&detections, ic);

//     let mut points: Vec<[f64; 3]> = Vec::new();
//     let mut intids: Vec<u64> = Vec::new();
//     for (i, p) in synced_points.iter().enumerate() {
//         if let Some(v) = p {
//             let (lon, lat) = unit_vector_to_lonlat(v[0], v[1], v[2]);
//             let idx = nested.hash(lon, lat);
//             *hp_map_counts.entry(idx as u64).or_insert(0) += 1;
//             points.push(*v);
//             intids.push(i as u64);
//         }
//     }

//     let nside = nside(depth) as u64;
//     let npix = 12 * nside * nside;

//     let radius = 2.0*(4.0*PI/(npix as f64)).sqrt();
//     let mut hp_sum_map_v2 = sum_hp_map_neighbors_v2(&hp_map_counts, nested, radius);
//     let peaks = find_non_strict_peaks(&hp_sum_map_v2, depth, MIN_UNIQUE_TIMES);

//     // Keep track of the mapping from intid → index in points array for efficient lookup after neighbor queries.
//     let mut intid_idx_hash: HashMap<i32, usize> = HashMap::new();
//     for (i, intid) in intids.iter().enumerate() {
//         intid_idx_hash.insert(*intid, i);
//     }
    
//     // Build a global 3D Kd-tree for all points, to enable fast neighbor queries.
//     let tree: Tree = ImmutableKdTree::new_from_slice(&points);
        
//     let mut seen: HashSet<Vec<i32>> = HashSet::new();
//     let mut not_skipped = 0;
//     for peak in peaks.iter() {
//         let cell = peak.0;
//         let cell_count = peak.1;
//         if cell_count < MIN_UNIQUE_TIMES {
//             continue;
//         }
//         not_skipped += 1;
                
//         // Get the center of the cell in (lon, lat), convert to (x, y, z) unit vector, 
//         // and query the global tree for neighbors within cluster_r.
//         let (lon, lat) = nested.center(cell);
//         let (cell_x, cell_y, cell_z) = lonlat_to_unit_vector(lon, lat);
//         let cell_point = [cell_x, cell_y, cell_z];

//         let mut cluster = tree.within::<SquaredEuclidean>(&cell_point, cluster_r2);

//         let res = nested.kth_neighbourhood(cell, 4u32);

//         // Hash the sorted intids and skip if we've already seen this cluster.
//         let mut intids_key: Vec<i32> = cluster
//             .iter()
//             .map(|nn| intids[nn.item as usize])
//             .collect();
//         intids_key.sort_unstable();

//         if !seen.insert(intids_key.clone()) {
//             continue;
//         }

//         for nn in cluster.iter_mut() {
//             nn.item = intids[nn.item as usize] as i32;
//         }

//         if cluster.len() < MIN_UNIQUE_TIMES {
//             continue;
//         }
                            
//         // Cull the cluster to remove outliers and construct the Cluster struct with local coordinates.
//         let clust = cull_cluster(&detections, &intid_idx_hash, &points, cell_point, &cluster, ic.epoch);
//         if clust.is_none() {
//             continue;
//         }
//         let clust = clust.unwrap();

//         let (t, theta_x, theta_y, ast_ucty, mxe, mye, xe, ye) = gather_cluster_data(&clust, &detections, ic.epoch);
        
//         let fit_result = iterative_reject(
//             &t, &mxe, &mye, &theta_x, &theta_y, &ast_ucty, &ast_ucty,
//             5.0,   // sigma threshold
//             MIN_UNIQUE_TIMES  // stop if fewer than a threhold number of points remain
//         );
        
//         let fit_final = match fit_result {
//             Ok(fit_final) => fit_final,
//             Err(e) => {
//                 //println!("Fit failed: {}", e);
//                 continue;
//             }
//         };
//         let (fit, kept_indices, rejected_indices) = fit_final;
//         let params = fit.params;

//         out.push_str(&format!("# Fit successful for cell {}: {:?}\n", cell, params));
                        
//         let rejected_intids: Vec<usize> = rejected_indices
//             .iter()
//             .map(|&idx| clust.local_intids[idx] as usize)
//             .collect();
        
//         out.push_str("# intid     detid         t(days)     theta_x(\")  theta_y(\")  x_model(\")  y_model(\")    sig_x     sig_y\n");
//         let mut kept_count = 0;
//         for i in 0..t.len() {
//             if rejected_indices.contains(&i) {
//                 continue;
//             }
//             kept_count += 1;
//             let x_model = params[0] + params[1] * t[i] + params[2] * mxe[i];
//             let y_model = params[3] + params[4] * t[i] + params[2] * mye[i];
//             let rx = theta_x[i] - x_model;
//             let pull_x = rx / ast_ucty[i];
//             let ry = theta_y[i] - y_model;
//             let pull_y = ry / ast_ucty[i];
//             let pull = pull_x.abs().max(pull_y.abs());
//             let intid = clust.local_intids[i];
//             let detid = &detections[intid as usize].detid;
//             out.push_str(&format!("{:8} {:12} {:12.6}   {:8.4}     {:8.4}    {:8.4}   {:8.4}.  {:8.4}.  {:8.4}\n", clust.local_intids[i], detid, t[i], theta_x[i]*206265., theta_y[i]*206265., x_model*206265., y_model*206265., pull_x, pull_y));
//         }
//         out.push_str(&format!("# Kept {} detections in fit, IC {}\n", kept_count, ii+1));
//         out.push_str("# detid         epoch        RA(deg)   sig_x(\")   Dec(deg)  sig_y(\")  obs_x        obs_y           obs_z      obscode  mag\n");
//         for i in 0..t.len() {
//             if rejected_indices.contains(&i) {
//                 continue;
//             }
            
//             let intid = clust.local_intids[i];
//             let detid = &detections[intid as usize].detid;
//             let objID = &detections[intid as usize].objID;
//             let filt = &detections[intid as usize].filt;
//             let rho_hat = detections[clust.local_intids[i] as usize].rho_hat; 
//             let lam = rho_hat[1].atan2(rho_hat[0]);
//             let beta = rho_hat[2].asin();
//             // These are ecliptic latitude and longitude at this point.
//             // Convert back to equatorial for output.
//             let (ra, dec) = ecliptic_to_equatorial(lam, beta, obliq * std::f64::consts::PI / 180.0);

//             let ra_deg = ra * 180.0 / PI;
//             let dec_deg: f64 = dec * 180.0 / PI;
//             let epoch = detections[clust.local_intids[i] as usize].epoch;
//             let obscode = obscode;
//             let mag = detections[clust.local_intids[i] as usize].mag;
//             let mag_sig = detections[clust.local_intids[i] as usize].mag_sig;
//             let obs_pos = &detections[clust.local_intids[i] as usize].observer_position;
//             let obs_pos2 = ecliptic_to_equatorial_vector(*obs_pos, obliq * std::f64::consts::PI / 180.0);

//             out.push_str(&format!("ic_00 {:12} {:12.6} {:10.6}   {:6.3}   {:10.6}   {:6.3}  {:13.10}  {:13.10}  {:13.10}   {}   {:.2} {:.2} {} {}\n", 
//                 detid, epoch, ra_deg, ast_ucty[i]*206265., dec_deg, ast_ucty[i]*206265., obs_pos2[0], obs_pos2[1], obs_pos2[2], obscode, mag, mag_sig, filt, objID));
//         }
            

//     }
//     out.push_str(&format!("# Finished processing IC {} of {}: {} cells, {} of which had > {} points.\n", ii+1, n_ics, peaks.len(), not_skipped, MIN_UNIQUE_TIMES));
//     out
// } 