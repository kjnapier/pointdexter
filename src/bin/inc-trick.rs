// use serde_yaml;
// use clap::Parser;
// use rayon::prelude::*;

// use spacerocks::SpiceKernel;
// use spacerocks::Time;
// use spacerocks::constants::SPEED_OF_LIGHT;
// use pointdexter::cli::{Config, Cli};

// use pointdexter::detection::Detection;
// use pointdexter::io::{load_detections, load_initial_conditions};
// use pointdexter::load_initial_conditions3d;
// use pointdexter::sync::*;

// use nalgebra::Vector3;

// use indicatif::ParallelProgressIterator;
// use indicatif::ProgressIterator;

// fn main() -> Result<(), Box<dyn std::error::Error>> {


//     let args = Cli::try_parse()?;
//     let f = std::fs::File::open(args.config)?;
//     let config: Config = serde_yaml::from_reader(f)?;

//     let mut kernel = SpiceKernel::new();
//     kernel.load_spk(format!("{}/sb441-n16.bsp", config.spice_path).as_str())?;
//     kernel.load_spk(format!("{}/de440s.bsp", config.spice_path).as_str())?;
//     kernel.load_bpc(format!("{}/earth_1962_240827_2124_combined.bpc", config.spice_path).as_str())?;

//     // Load detections from catalog. They will automatically be transformed to the desired reference plane. 
//     // The observer positions and velocities will be rotated accordingly.
//     println!("Loading detections...");
//     let mut detections = load_detections(&config.detection_catalog, &config.orbit_reference_plane, &kernel)?;
//     println!("Loaded {} detections.", detections.len());

//     // Load initial conditions from file.
//     println!("Loading initial conditions...");
//     let mut ics = load_initial_conditions3d(&config.initial_conditions_file, &config.ic_origin, config.reference_epoch)?;
//     println!("Loaded {} initial conditions.", ics.len());

//     const ARCSEC_PER_RAD: f64 = (180.0 * 3600.0) / std::f64::consts::PI;

//     // sample inclinations uniformly in cos(theta)
//     let mut inclinations: Vec<f64> = (0..100).map(|i| {
//         let cos_i = 1.0 - 2.0 * (i as f64 + 0.5) / 100.0; // +0.5 for mid-bin sampling
//         cos_i.acos()
//     }).collect();

//     for ic in ics.iter() {

//         let mut rs: Vec<f64> = Vec::with_capacity(detections.len());
//         let mut vrs: Vec<f64> = Vec::with_capacity(detections.len());
//         let mut vos: Vec<f64> = Vec::with_capacity(detections.len());
//         let mut fs: Vec<f64> = Vec::with_capacity(detections.len());
//         let mut gs: Vec<f64> = Vec::with_capacity(detections.len());
//         let mut rhats: Vec<Vector3<f64>> = Vec::with_capacity(detections.len());
//         let mut ahats: Vec<Vector3<f64>> = Vec::with_capacity(detections.len());
//         let mut dhats: Vec<Vector3<f64>> = Vec::with_capacity(detections.len());
//         let mut cos_thetas: Vec<f64> = Vec::with_capacity(detections.len());

//         for det in &detections {
//             let rho_opt = optimize_rho_v2(det, ic, ic.r - 1.0, 1.0); // Keep in mind that there can be multiple solutions, but not for TNOs.

//             let (rho, r) = rho_opt.unwrap_or_else(|| {
//                 println!("Failed to optimize rho for detection at epoch {}. Using fallback values.", det.epoch);
//                 (ic.r, ic.r) // fallback to something reasonable
//             });
            
//             let r_vec = det.observer_position + rho * det.rho_hat;
//             let light_corrected_epoch = det.epoch - rho / SPEED_OF_LIGHT;
//             let (f, g) = ic.fg_at_epoch(ic.epoch + (ic.epoch - light_corrected_epoch));

//             let (x, y, z) = (r_vec[0], r_vec[1], r_vec[2]);
//             let xy_norm = (x * x + y * y).sqrt();

//             let ahat = Vector3::new(-y / xy_norm, x / xy_norm, 0.0);
//             let dhat = Vector3::new(-z * x / (r * xy_norm), -z * y / (r * xy_norm), xy_norm / r);

//             let cos_theta = xy_norm / r;

//             let vo = ic.h / r;
//             let vsq = (ic.energy + ic.mu / r).max(0.0) * 2.0; // ensure non-negative vsq
//             let vr = (vsq - vo * vo).max(0.0).sqrt(); // ensure non-negative argument for sqrt

//             rs.push(r);
//             vrs.push(vr);
//             vos.push(vo);
//             fs.push(f);
//             gs.push(g);
//             rhats.push(r_vec / r);
//             ahats.push(ahat);
//             dhats.push(dhat);
//             cos_thetas.push(cos_theta);
//         }

//         for inc in inclinations.iter() {

//             let start = std::time::Instant::now();
            
//             let mut synced_points: Vec<Option<[f64; 3]>> = Vec::with_capacity(detections.len());
//             let cos_inc = inc.cos();

//             let mut latitude_threshold = inc.clone();
//             if latitude_threshold > std::f64::consts::PI / 2.0 {
//                 latitude_threshold = std::f64::consts::PI - latitude_threshold;
//             }
//             let sin_latitude_threshold = latitude_threshold.sin();

//             for idx in 0..detections.len() {

//                 let rhat = rhats[idx];
//                 if (rhat[2]).abs() > sin_latitude_threshold {
//                     synced_points.push(None);
//                     continue;
//                 }

//                 let r = rs[idx];
//                 let vr = vrs[idx];
//                 let vo = vos[idx];
//                 let f = fs[idx];
//                 let g = gs[idx];
//                 let ahat = ahats[idx];
//                 let dhat = dhats[idx];
//                 let cos_theta = cos_thetas[idx];

//                 let cos_psi = cos_inc / cos_theta;
//                 let sin_psi = (1.0 - cos_psi * cos_psi).clamp(0.0, 1.0).sqrt();

//                 let v_vec = vr * rhat + vo * (cos_psi * ahat + sin_psi * dhat);
//                 let r_vec = r * rhat;

//                 let new_position = r_vec * f + v_vec * g;
//                 let new_pointing = new_position / new_position.norm();
//                 synced_points.push(Some([new_pointing[0], new_pointing[1], new_pointing[2]]));
//             }

//             if synced_points.iter().all(|p| p.is_none()) {
//                 println!("All points are outside the latitude threshold for inclination {}. Skipping.", inc);
//                 continue;
//             }

//             // // now accumulate the synced points into a 2d histogram
//             // let mut phi_min = std::f64::INFINITY;
//             // let mut phi_max = std::f64::NEG_INFINITY;
//             // let mut theta_min = std::f64::INFINITY;;
//             // let mut theta_max = std::f64::NEG_INFINITY;
//             // let mut epoch_min = std::f64::INFINITY;
//             // let mut epoch_max = std::f64::NEG_INFINITY;
//             // let mut thetas: Vec<f64> = Vec::with_capacity(synced_points.len());
//             // let mut phis: Vec<f64> = Vec::with_capacity(synced_points.len());
//             // let mut epochs: Vec<f64> = Vec::with_capacity(synced_points.len());
//             // for (idx, point) in synced_points.iter().enumerate() {
//             //     if point.is_none() {
//             //         continue;
//             //     }
//             //     let point = point.unwrap();
//             //     let x = point[0];
//             //     let y = point[1];
//             //     let z = point[2];
//             //     let theta = (z).asin();
//             //     let phi = y.atan2(x);
//             //     let epoch = detections[idx].epoch;

//             //     if phi < phi_min {
//             //         phi_min = phi;
//             //     }
//             //     if phi > phi_max {
//             //         phi_max = phi;
//             //     }
//             //     if theta < theta_min {
//             //         theta_min = theta;
//             //     }
//             //     if theta > theta_max {
//             //         theta_max = theta;
//             //     }
//             //     if epoch < epoch_min {
//             //         epoch_min = epoch;
//             //     }
//             //     if epoch > epoch_max {
//             //         epoch_max = epoch;
//             //     }
//             //     thetas.push(theta);
//             //     phis.push(phi);
//             //     epochs.push(epoch);
//             // }

//             // // make a histogram grid with 5 arcsec bins
//             // let bin_size = 5.0 / ARCSEC_PER_RAD;
//             // // make time bins in a third dimension.
//             // let time_bin_size = 15.0; // days
//             // let n_time_bins = ((detections.iter().map(|d| d.epoch).fold(0./0., f64::max) - detections.iter().map(|d| d.epoch).fold(0./0., f64::min)) / time_bin_size).ceil() as usize;
//             // let n_phi_bins = ((phi_max - phi_min) / bin_size).ceil() as usize;
//             // let n_theta_bins = ((theta_max - theta_min) / bin_size).ceil() as usize;
//             // let mut histogram = vec![vec![vec![0usize; n_time_bins]; n_phi_bins]; n_theta_bins];

//             // let mut peaks: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
//             // for idx in 0..phis.len() {
//             //     let theta = thetas.get(idx);
//             //     let phi = phis.get(idx);
//             //     let theta = *theta.unwrap();
//             //     let phi = *phi.unwrap();
//             //     let epoch = epochs.get(idx);
//             //     let epoch = *epoch.unwrap();
                
//             //     let phi_bin = ((phi - phi_min) / bin_size).floor() as usize;
//             //     let theta_bin = ((theta - theta_min) / bin_size).floor() as usize;
//             //     let time_bin = ((epoch - epoch_min) / time_bin_size).floor() as usize;

//             //     histogram[theta_bin][phi_bin][time_bin] += 1;
//             //     if histogram[theta_bin][phi_bin][time_bin] > 12 {
//             //         peaks.insert((theta_bin, phi_bin));
//             //     }
//             // }
//             // println!("Found {} peaks for inclination {}.", peaks.len(), inc);
//             let duration = start.elapsed();
//             println!("Processed IC at epoch {} in {:?}", ic.epoch, duration);

//             // make a filename and save the histogram to a text file
//             // let filename = format!("/Users/kjnapier/Desktop/deh/histogram_ic_epoch_{}_inc_{:.4}.txt", ic.epoch, inc);
//             // let mut file = std::fs::File::create(&filename)?;
//             // use std::io::Write;
//             // for theta_bin in 0..n_theta_bins {
//             //     for phi_bin in 0..n_phi_bins {
//             //         write!(file, "{} ", histogram[theta_bin][phi_bin])?;
//             //     }
//             // }
//         }

        

        

//     }   
//     Ok(())
// }



pub fn main() {
    println!("This is a placeholder for the inc-trick binary. The actual code is in matt-0.rs, which is being used for testing and development. This file will eventually be removed.");
}