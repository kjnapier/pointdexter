use std::f64::consts::PI;
use std::collections::HashMap;

use rayon::prelude::*;

use nalgebra::{DMatrix, DVector, Matrix3, Matrix3xX, Unit, Vector3, matrix};

use kiddo::KdTree;
use kiddo::SquaredEuclidean;
use kiddo::NearestNeighbour;

use std::collections::HashSet;
use ordered_float::OrderedFloat;

use std::fs::File;
use std::io;
use std::io::Write;

use std::time::Instant;

use cdshealpix::nested::{get, Layer};
use cdshealpix::{TRANSITION_LATITUDE};

use cdshealpix::{nside};
use cdshealpix::nested;

use pointdexter::Detection;
use pointdexter::InitialCondition;
use pointdexter::sync::sync_detection_to_orbit_with_time;
use pointdexter::load_initial_conditions;
use pointdexter::load_detections;
use pointdexter::cli::{Config, Cli};

use serde_yaml;
use clap::Parser;

use spacerocks::SpiceKernel;


const ARCSEC_PER_HOUR_TO_RAD_PER_DAY: f64 = (24.0 / 3600.0) * PI / 180.0;
const ARCSEC_PER_RAD: f64 = 3600.0 * 180.0 / PI;


// Function to convert a unit vector (x, y, z) to (lon, lat) in radians
fn unit_vector_to_lonlat(x: f64, y: f64, z: f64) -> (f64, f64) {
    let lon = y.atan2(x); // atan2(y, x) for longitude
    let lat = z.asin();    // asin(z) for latitude
    (lon, lat)
}

fn lonlat_to_unit_vector(lon: f64, lat: f64) -> (f64, f64, f64) {
    let x = lat.cos() * lon.cos();
    let y = lat.cos() * lon.sin();
    let z = lat.sin();
    (x, y, z)
}

pub fn xyz_to_proj_matrix(r_ref: Vector3<f64>) -> Matrix3<f64> {
    let x_ref = r_ref.x;
    let y_ref = r_ref.y;
    let z_ref = r_ref.z;

    let r = (x_ref * x_ref + y_ref * y_ref + z_ref * z_ref).sqrt();
    let lon0 = y_ref.atan2(x_ref);
    let lat0 = (z_ref / r).asin();

    let slon0 = lon0.sin();
    let clon0 = lon0.cos();
    let slat0 = lat0.sin();
    let clat0 = lat0.cos();

    Matrix3::new(
        -slon0,           clon0,            0.0,
        -clon0 * slat0,  -slon0 * slat0,   clat0,
         clon0 * clat0,   slon0 * clat0,   slat0,
    )
}

pub fn find_unique_f64_hashset(data: Vec<f64>) -> Vec<f64> {
    let unique_elements: HashSet<OrderedFloat<f64>> = data.into_iter().map(OrderedFloat).collect();
    unique_elements.into_iter().map(|x| x.0).collect()
}



#[derive(Debug)]
struct FitResult {
    params: [f64; 5],  // a, b, c, d, e
    cov: DMatrix<f64>, // 5×5 covariance matrix
    sigma: f64,
    r2: f64,
    chi2: f64,
    dof: f64,
}

/// Weighted coupled least squares with explicit error handling.
fn coupled_weighted_fit(
    t: &[f64],
    f: &[f64],
    g: &[f64],
    x: &[f64],
    y: &[f64],
    sigma_x: &[f64],
    sigma_y: &[f64],
) -> Result<FitResult, String> {
    let n = t.len();
    if n < 3 {
        return Err("Not enough data points to fit.".into());
    }

    // Build design matrix (2n × 5) and observation vector (2n)
    let mut xmat = DMatrix::<f64>::zeros(2 * n, 5);
    let mut yvec = DVector::<f64>::zeros(2 * n);
    let mut w = DVector::<f64>::zeros(2 * n); // inverse variance weights

    for i in 0..n {
        // x_i = a + b*t + c*f
        xmat[(i, 0)] = 1.0;
        xmat[(i, 1)] = t[i];
        xmat[(i, 2)] = f[i];
        yvec[i] = x[i];
        w[i] = 1.0 / sigma_x[i].powi(2);

        // y_i = d + e*t + c*g
        let j = n + i;
        xmat[(j, 2)] = g[i];
        xmat[(j, 3)] = 1.0;
        xmat[(j, 4)] = t[i];
        yvec[j] = y[i];
        w[j] = 1.0 / sigma_y[i].powi(2);
    }

    let w_diag = DMatrix::from_diagonal(&w);
    let xtwx = xmat.transpose() * &w_diag * &xmat;
    let xtwy = xmat.transpose() * &w_diag * &yvec;

    // Try to invert (XᵀWX)
    let xtwx_inv = xtwx.try_inverse().ok_or_else(|| {
        "Matrix (XᵀWX) is singular or ill-conditioned — cannot perform fit.".to_string()
    })?;

    // Solve for parameters
    let beta = &xtwx_inv * xtwy;
    let params = [beta[0], beta[1], beta[2], beta[3], beta[4]];

    // Compute residuals and statistics
    let y_hat = &xmat * &beta;
    let r = &yvec - y_hat;
    let chi2 = r.iter().zip(w.iter()).map(|(&ri, &wi)| wi * ri.powi(2)).sum::<f64>();
    let dof = (2 * n) as f64 - 5.0;
    if dof <= 0.0 {
        return Err("Not enough degrees of freedom.".into());
    }

    let sigma2 = chi2 / dof;
    let cov = &xtwx_inv * sigma2;

    // Weighted R²
    let mean_y = yvec.iter().zip(w.iter()).map(|(&yi, &wi)| wi * yi).sum::<f64>() / w.sum();
    let ss_tot = yvec
        .iter()
        .zip(w.iter())
        .map(|(&yi, &wi)| wi * (yi - mean_y).powi(2))
        .sum::<f64>();
    let r2 = 1.0 - chi2 / ss_tot;

    Ok(FitResult {
        params,
        cov,
        sigma: sigma2.sqrt(),
        r2,
        chi2,
        dof,
    })
}

/// Iteratively reject the single largest outlier until |pull| < threshold or min_points reached.
fn iterative_reject(
    t: &[f64],
    f: &[f64],
    g: &[f64],
    x: &[f64],
    y: &[f64],
    sigma_x: &[f64],
    sigma_y: &[f64],
    pull_threshold: f64,
    min_points: usize,
) -> Result<(FitResult, Vec<usize>, Vec<usize>), String> {
    assert!(t.len() == f.len() && t.len() == g.len());
    let mut keep: Vec<usize> = (0..t.len()).collect();
    let mut rejected: Vec<usize> = Vec::new();

    loop {
        // Subset data
        let t_f: Vec<f64> = keep.iter().map(|&i| t[i]).collect();

        let mut t_n: Vec<i64> = t_f.iter().map(|&ti| ti.round() as i64).collect();
        t_n.sort_unstable();
        t_n.dedup();

        let f_f: Vec<f64> = keep.iter().map(|&i| f[i]).collect();
        let g_f: Vec<f64> = keep.iter().map(|&i| g[i]).collect();
        let x_f: Vec<f64> = keep.iter().map(|&i| x[i]).collect();
        let y_f: Vec<f64> = keep.iter().map(|&i| y[i]).collect();
        let sx_f: Vec<f64> = keep.iter().map(|&i| sigma_x[i]).collect();
        let sy_f: Vec<f64> = keep.iter().map(|&i| sigma_y[i]).collect();

        if t_n.len() < MIN_NIGHTS {
            return Err(format!(
                "Stopped: not enough unique observation times left ({} < {}).",
                t_n.len(),
                MIN_UNIQUE_TIMES
            ));
        }
        // Perform fit
        let fit = match coupled_weighted_fit(&t_f, &f_f, &g_f, &x_f, &y_f, &sx_f, &sy_f) {
            Ok(fit) => fit,
            Err(e) => {
                eprintln!("Fit failed: {e}");
                return Err(format!("Fit failed after rejecting {} points: {}", rejected.len(), e));
            }
        };

        // Compute pulls
        let mut worst_i_local: Option<usize> = None;
        let mut worst_pull = 0.0;
        for (i, &ti) in t_f.iter().enumerate() {
            let xi_model = fit.params[0] + fit.params[1] * ti + fit.params[2] * f_f[i];
            let yi_model = fit.params[3] + fit.params[4] * ti + fit.params[2] * g_f[i];
            let pull_x = (x_f[i] - xi_model) / sx_f[i];
            let pull_y = (y_f[i] - yi_model) / sy_f[i];
            let pull_mag = pull_x.abs().max(pull_y.abs());
            if pull_mag > worst_pull {
                worst_pull = pull_mag;
                worst_i_local = Some(i);
            }
        }

        // Stop if all points are within threshold
        if worst_pull < pull_threshold {
            return Ok((fit, keep, rejected));
        }

        // Reject the worst point
        if let Some(local_i) = worst_i_local {
            let global_i = keep[local_i];
            /*
            println!(
                "Rejecting point {} with |pull| = {:.2} (remaining: {})",
                global_i,
                worst_pull,
                keep.len() - 1
            );
            */
            
            rejected.push(global_i);
            keep.remove(local_i);
        }

        // Stop if not enough data left
        if keep.len() < min_points {
            return Err(format!(
                "Stopped: not enough points left ({} < {}).",
                keep.len(),
                min_points
            ));
        }
    }
}


// --- Point ---
pub struct Point {
    vec: Vector3<f64>,
    ast_ucty: Option<f64>,
    epoch: f64,
    detid: Option<String>,
    theta_x: f64,
    theta_y: f64,
}

#[derive(Debug)]
pub struct Cluster {
    ref_epoch: f64,
    local_intids: Vec<u64>,
    ref_point: [f64; 3],
    local_points: Vec<[f64; 3]>,
    local_thetas: Vec<[f64; 2]>,
    local_thetas_hash: HashMap<u64, [f64; 2]>,
    local_tree: KdTree<f64, 2>,
    local_obs_positions: Vec<Vector3<f64>>,
    unique_cluster_times: Vec<f64>,
    }
    
impl Cluster {
    pub fn from_intids(
        central_intid: i32,
        local_intids: Vec<u64>,
        detections: &Vec<Detection>,
        central_point: [f64; 3],
        points: &Vec<[f64; 3]>,
    ) -> Self {
        let central_detection_epoch = detections[central_intid as usize].epoch;
        
        let mut local_points: Vec<[f64; 3]> = local_intids
                .iter()
                .map(|&intid| points[intid as usize])
                .collect();

        let full_times: Vec<f64> = local_intids
                .iter()
                .map(|&intid| detections[intid as usize].epoch - central_detection_epoch)
                .collect();

        let mut unique_cluster_times = find_unique_f64_hashset(full_times.clone());
        
        // Transform to local tangent plane coordinates.
        let point_vec = Vector3::new(central_point[0], central_point[1], central_point[2]);
        let mat = xyz_to_proj_matrix(point_vec);

        let mut test_points: Vec<Point> = Vec::new();
        let mut local_thetas: Vec<[f64; 2]> = Vec::new();
        let mut local_thetas_hash: HashMap<u64, [f64; 2]> = HashMap::new();
        let mut local_obs_positions: Vec<Vector3<f64>> = Vec::new();

        for (p, local_intid) in local_points
            .iter()
            .map(|&p| p)
            .zip(local_intids.iter().map(|&i| i as u64))
        {
            let det = &detections[local_intid as usize];
            let detid = det.detid.clone();

            let dt = det.epoch - central_detection_epoch;
            let p_vec = Vector3::new(p[0], p[1], p[2]);
            let longitude = p_vec[0].atan2(p_vec[1]);
            let latitude = p_vec[2].asin();
            let proj_vec = mat * p_vec;
            let theta_x = proj_vec[0] / proj_vec[2];
            let theta_y = proj_vec[1] / proj_vec[2];

            let position = det.observer_position.clone();
            let local_obs_pos = mat * position;
            local_obs_positions.push(local_obs_pos);

            let test_point: Point = Point {
                vec: Vector3::new(p[0], p[1], p[2]),
                ast_ucty: det.ast_ucty,
                epoch: det.epoch,
                detid: detid.clone(),
                theta_x: theta_x,
                theta_y: theta_y,
            };
            test_points.push(test_point);   

            local_thetas.push([theta_x, theta_y]);
            local_thetas_hash.insert(local_intid as u64, [theta_x, theta_y]);
            println!("{:10} {:10.4} {:10.4} {:12.9} {:12.9} {:10.6} {:?}", local_intid, theta_x*206265., theta_y*206265., longitude, latitude, det.epoch, detid);
        }

        let mut local_tree: KdTree<f64, 2> = KdTree::with_capacity(10);
            local_tree.extend(
                local_thetas
                    .iter()
                    .map(|&p| p)
                    .zip(local_intids.iter().map(|&i| i as u64)),
            );

        Cluster {
            //intid: central_intid,
            ref_epoch: central_detection_epoch,
            local_intids,
            ref_point: central_point,
            local_points,
            local_thetas,
            local_thetas_hash,
            local_tree,
            local_obs_positions,
            unique_cluster_times,
        }
    }
}





fn cull_cluster(detections: &Vec<Detection>, indices: &Vec<usize>, points: &Vec<[f64; 3]>, point: [f64; 3], cluster: &Vec<NearestNeighbour<f64, u64>>, central_detection_epoch: f64) -> Option<Cluster> {

    // Prune detections that are not compatible with the parent/central detection.

    let mut cluster_hash = Vec::new();
    let mut min_epoch = f64::MAX;
    let mut max_epoch = f64::MIN;

    for neighbor in cluster.clone().into_iter() {
        
        let det = &detections[neighbor.item as usize];

        let dt = (det.epoch - central_detection_epoch).abs();
        let distance = neighbor.distance;


        if dt > MAX_DT {
            continue;
        }

        if det.epoch < min_epoch {
            min_epoch = det.epoch;
        }

        if det.epoch > max_epoch {
            max_epoch = det.epoch;
        }

        cluster_hash.push(neighbor.item as u64);

    }

    let duration = max_epoch - min_epoch;


    let cluster_times: Vec<f64> = cluster_hash
        .iter()
        .map(|&id| detections[id as usize].epoch - central_detection_epoch)
        .collect();

    let mut cluster_nights: Vec<i32> = cluster_times
        .iter()
        .map(|&x| x.round() as i32)
        .collect();

    cluster_nights.sort_unstable();
    cluster_nights.dedup();

    if cluster_nights.len() < MIN_NIGHTS {
        return None;
    }
        

    let mut unique_cluster_times = find_unique_f64_hashset(cluster_times.clone());

    if unique_cluster_times.len() < MIN_UNIQUE_TIMES {
        return None;
    }

    // Transform to local tangent plane coordinates.
    let point_vec = Vector3::new(point[0], point[1], point[2]);
    let mat = xyz_to_proj_matrix(point_vec);

    cluster_hash.sort();
    let mut local_intids: Vec<u64> = cluster_hash.iter().cloned().collect();

    let mut local_points: Vec<[f64; 3]> = indices
                .iter()
                .map(|&idx| points[idx].clone())
                .collect();

    let mut local_thetas: Vec<[f64; 2]> = Vec::new();
    let mut local_thetas_hash: HashMap<u64, [f64; 2]> = HashMap::new();
    let mut local_obs_positions: Vec<Vector3<f64>> = Vec::new();
    for (p, local_intid) in local_points
        .iter()
        .map(|&p| p)
        .zip(local_intids.iter().map(|&i| i as u64))
    {
        let det = &detections[local_intid as usize];
        let detid = det.detid.clone();

        let dt = det.epoch - central_detection_epoch;
        let p_vec = Vector3::new(p[0], p[1], p[2]);
        let longitude = p_vec[0].atan2(p_vec[1]);
        let latitude = p_vec[2].asin();
        let proj_vec = mat * p_vec;
        let theta_x = proj_vec[0] / proj_vec[2];
        let theta_y = proj_vec[1] / proj_vec[2];

        let position = det.observer_position.clone();
        let local_obs_position = mat * position;
        local_obs_positions.push(local_obs_position);

        local_thetas.push([theta_x, theta_y]);
        local_thetas_hash.insert(local_intid as u64, [theta_x, theta_y]);
    }

    //    println!("");
    let mut local_tree: KdTree<f64, 2> = KdTree::with_capacity(10);
    local_tree.extend(
        local_thetas
            .iter()
            .map(|&p| p)
            .zip(local_intids.iter().map(|&i| i as u64)),
    );

    //println!("large cluster: {:?}", local_intids.len());
    return Some(Cluster {
        //intid: central_intid,
        ref_epoch: central_detection_epoch,
        local_intids: local_intids,
        ref_point: point,
        local_points: local_points,
        local_thetas: local_thetas,
        local_thetas_hash: local_thetas_hash,
        local_tree: local_tree,
        local_obs_positions: local_obs_positions,
        unique_cluster_times: unique_cluster_times,
    })
}



const EPS_SMALL: f64 = 4.0 / ARCSEC_PER_RAD; // 45.0 / ARCSEC_PER_RAD;
const SEP_THRESH: f64 = 4.0;


fn main() -> Result<(), Box<dyn std::error::Error>> {

    let args = Cli::try_parse()?;
    let f = std::fs::File::open(args.config)?;
    let config: Config = serde_yaml::from_reader(f)?;

    const EPS = config.epsilon_arcsec / ARCSEC_PER_RAD;
    const MIN_POINTS_PER_CELL: usize = config.min_detections;
    const MIN_NIGHTS: usize = config.min_nites;
    const MAX_DT: f64 = 4000.0;
    const MIN_UNIQUE_TIMES: usize = 12;


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
    let mut initial_conditions = load_initial_conditions(&config.initial_conditions_file, &config.ic_type, &config.ic_origin, config.reference_epoch)?;
    println!("Loaded {} initial conditions.", initial_conditions.len());

    let depth = 13_u8;
    let nside = nside(depth) as u64;
    let npix = 12 * nside * nside;
    println!("Nside at depth {}: {} {}", depth, nside, (4.0*PI/(npix as f64)).sqrt()*(180.0/PI)*60.0*60.0);
    let nested: &Layer = get(depth);

    let cluster_r2 = 2.0 * (1.0 - EPS.cos());
    let cluster_small_r2 = 2.0 * (1.0 - EPS_SMALL.cos());

    let mut points: Vec<[f64; 3]> = Vec::new();
    // open a file to write the clusters to
    let mut cluster_file = std::fs::File::create("clusters.txt")?;
    // write the header


    let start_time = Instant::now();
    let n_ics = initial_conditions.len();

    let mut cluster_set: HashSet<Vec<u64>> = HashSet::new();

    let mut n_clusters = 0;
    let mut ii = 0;
    //for ic in &initial_conditions[..n_ics] {
    // multi-thread over initial_conditions
    initial_conditions[..10].par_iter().enumerate().for_each(|(ii, ic)| {

        println!("Processing initial condition {}/{}...", ii+1, n_ics);

        let mut points: Vec<[f64; 3]> = Vec::with_capacity(detections.len());
        let mut intids: Vec<i32> = Vec::with_capacity(detections.len());
        let mut det_ids: Vec<Option<String>> = Vec::with_capacity(detections.len());

        let mut hp_map: HashMap<u64, Vec<i32>> = HashMap::new();

        // Transform all detections to points on the unit sphere according to the current orbit.
        let mut light_corrected_epochs: Vec<f64> = Vec::with_capacity(detections.len());
        for (idx, detection) in detections[..].iter().enumerate() {
            if let Some((pointing, light_corrected_epoch)) = sync_detection_to_orbit_with_time(&detection, &ic) {
                let (lon, lat) = unit_vector_to_lonlat(pointing[0], pointing[1], pointing[2]);
                
                let hp_idx = nested.hash(lon, lat);

                hp_map.entry(hp_idx as u64).or_insert_with(Vec::new).push(idx as i32);

                points.push(pointing);
                intids.push(idx as i32);
                light_corrected_epochs.push(light_corrected_epoch);
                det_ids.push(detection.detid.clone());
            }
        }
        println!("Transformed {} detections to points for IC {}.", points.len(), ic.id);

         // Iterate over the HashMap
        // Create a new HashMap to store counts of points per cell
        let mut cell_count_map: HashMap<usize, Vec<u64>> = HashMap::new();
        let mut total_points = 0;
        for (key, values) in &hp_map {
            let count = values.len();
            total_points += count;
            cell_count_map.entry(count).or_insert_with(Vec::new).push(*key);
        }

        let mut sorted_keys: Vec<&usize> = cell_count_map.keys().collect();
        sorted_keys.sort_by(|a, b| b.cmp(a));

        let mut ordered_cells = Vec::new();
        for key in sorted_keys {
            let count_vec = cell_count_map.get(key).unwrap();
            ordered_cells.extend(count_vec);
        }

        let mut tree = KdTree::new();

        // let mut intid_idx_hash: HashMap<u64, usize> = HashMap::new();
        // for (i, intid) in intids.iter().enumerate() {
        //     intid_idx_hash.insert(*intid as u64, i);
        // }
        let mut indices: Vec<usize> = Vec::with_capacity(intids.len());
        for (i, _) in intids.iter().enumerate() {
            indices.push(i);
        }

        tree.extend(
            points
                .iter()
                .map(|&p| p)
                .zip(indices.iter().map(|&i| i as u64)),
        );

        // println!("Built KD-Tree with {} points.", points.len());

        // print how many ordered_cells there are
        // println!("Found {} ordered cells.", ordered_cells.len());

        let mut tracklet_id = 0;
        for cell in ordered_cells.iter() {
            let cell_intids = &hp_map[cell];
            if cell_intids.len() < MIN_POINTS_PER_CELL {
                // println!("Skipping cell {} with only {} points.", cell, cell_intids.len());
                break;
            }
            let (lon, lat) = nested.center(*cell);
            let (cell_x, cell_y, cell_z) = lonlat_to_unit_vector(lon, lat);
            let cell_point = [cell_x, cell_y, cell_z];
            let cluster = tree.within::<SquaredEuclidean>(&cell_point, cluster_r2);
            // println!("Cell {}: found {} points in cluster.", cell, cluster.len());
            let clust = cull_cluster(&detections, &indices, &points, cell_point, &cluster, ic.epoch);
            if clust.is_none() {
                continue;
            }
            let clust = clust.unwrap();

            // let ast_ucty: Vec<f64> = clust.local_intids
            //     .iter()
            //     .map(|id| detections[*id as usize].ast_ucty.clone())
            //     .collect();

            // let ast_ucty: Vec<Option<f64>> = clust.local_intids
            //     .iter()
            //     .map(|id| detections[*id as usize].ast_ucty)
            //     .collect();

            // make ast_ucty a vector of 0.15 arcsec in radians
            let ast_ucty: Vec<f64> = clust.local_intids
                .iter()
                .map(|_id| (0.15 / 3600.0) * (PI / 180.0))
                .collect();

            let theta_x : Vec<f64> = clust.local_thetas
                .iter()
                .map(|xy| xy[0])
                .collect();
            let theta_y : Vec<f64> = clust.local_thetas
                .iter()
                .map(|xy| xy[1])
                .collect();
            let t: Vec<f64> = clust.local_intids
                .iter()
                .map(|id| detections[*id as usize].epoch - ic.epoch)
                .collect();
            let mxe : Vec<f64> = clust.local_obs_positions
                .iter()
                .map(|pos| -pos[0])
                .collect();
            let mye : Vec<f64> = clust.local_obs_positions
                .iter()
                .map(|pos| -pos[1])
                 .collect();

            
            let fit_result = iterative_reject(
                &t, &mxe, &mye, &theta_x, &theta_y, &ast_ucty, &ast_ucty,
                5.0,   // 3σ threshold
                MIN_POINTS_PER_CELL     // stop if fewer than 10 points remain
            );

            // println!("# Fit result: {:?}", fit_result);
            
            let fit_final = match fit_result {
                Ok(fit_final) => fit_final,
                Err(e) => {
                    continue;
                }
            };

            let (fit, kept_indices, rejected_indices) = fit_final;

            let params = fit.params;
            println!("# IC {} of {}: {}", ii + 1, n_ics, ic.id);
            println!("# Fit successful for cell {}: {:?}", cell, params);
            
            let rejected_intids: Vec<usize> = rejected_indices
                .iter()
                .map(|&idx| clust.local_intids[idx] as usize)
                .collect();

            println!("# intid     detid         t(days)     theta_x(\")  theta_y(\")  x_model(\")  y_model(\")    sig_x     sig_y"); 
            for i in 0..t.len() {
                if rejected_indices.contains(&i) {
                    continue;
                }
                let x_model = params[0] + params[1] * t[i] + params[2] * mxe[i];
                let y_model = params[3] + params[4] * t[i] + params[2] * mye[i];
                let rx = theta_x[i] - x_model;
                let pull_x = rx / ast_ucty[i];
                let ry = theta_y[i] - y_model;
                let pull_y = ry / ast_ucty[i];
                let pull = pull_x.abs().max(pull_y.abs());
                let intid = clust.local_intids[i];
                let detid = &detections[intid as usize].detid;
                println!("{:8} {:?} {:12.6}   {:8.4}     {:8.4}    {:8.4}   {:8.4}.  {:8.4}.  {:8.4}", clust.local_intids[i], detid, t[i], theta_x[i]*206265., theta_y[i]*206265., x_model*206265., y_model*206265., pull_x, pull_y);
            }

            println!("# detid         epoch        lon(deg)   sig_x(\")   lat(deg)  sig_y(\")  mag");
            for i in 0..t.len() {
                if rejected_indices.contains(&i) {
                    continue;
                }
                let intid = clust.local_intids[i];
                let detid = detections[intid as usize].detid.clone().unwrap();
                let rho_hat = detections[clust.local_intids[i] as usize].rho_hat; 
                let lam = rho_hat[1].atan2(rho_hat[0]);
                let beta = rho_hat[2].asin();
                let epoch = detections[clust.local_intids[i] as usize].epoch;
                let mag = detections[clust.local_intids[i] as usize].mag.unwrap_or(-1.0);

                // println!("{:12} {:12} {:12.6} {:9.5}   {:5.2}   {:9.5}   {:5.2}  {:13.10}  {:13.10}  {:13.10}   {}   {:.2} {:.2} {}", 
                //     ic.id, detid, epoch, ra_deg, ast_ucty[i]*206265., dec_deg, ast_ucty[i]*206265., obs_pos2[0], obs_pos2[1], obs_pos2[2], obscode, mag, mag_sig, filt);
                println!("{:12} {:12} {:12.6} {:9.5}   {:5.2}   {:9.5}   {:5.2}   {:.2}", 
                    ic.id, detid, epoch, lam, ast_ucty[i]*206265., beta, ast_ucty[i]*206265., mag);

            }
             
  
        }
    });
            
    // };

    let duration = start_time.elapsed();
    println!("Time elapsed: {:?}", duration);
       

    Ok(())
}





// //     Ok(())
// // }


//  // Iterate over the HashMap
//         // Create a new HashMap to store counts of points per cell
//         let mut cell_count_map: HashMap<usize, Vec<u64>> = HashMap::new();
//         let mut total_points = 0;
//         for (key, values) in &hp_map {
//             let count = values.len();
//             total_points += count;
//             cell_count_map.entry(count).or_insert_with(Vec::new).push(*key);
//         }

//         let mut sorted_keys: Vec<&usize> = cell_count_map.keys().collect();
//         sorted_keys.sort_by(|a, b| b.cmp(a));

//         let mut ordered_cells = Vec::new();
//         for key in sorted_keys {
//             let count_vec = cell_count_map.get(key).unwrap();
//             ordered_cells.extend(count_vec);
//         }

//         let mut tree = KdTree::new();

//         // let mut intid_idx_hash: HashMap<u64, usize> = HashMap::new();
//         // for (i, intid) in intids.iter().enumerate() {
//         //     intid_idx_hash.insert(*intid as u64, i);
//         // }
//         let mut indices: Vec<usize> = Vec::with_capacity(intids.len());
//         for (i, _) in intids.iter().enumerate() {
//             indices.push(i);
//         }

//         tree.extend(
//             points
//                 .iter()
//                 .map(|&p| p)
//                 .zip(indices.iter().map(|&i| i as u64)),
//         );

//         // println!("Built KD-Tree with {} points.", points.len());

//         // print how many ordered_cells there are
//         // println!("Found {} ordered cells.", ordered_cells.len());

//         let mut tracklet_id = 0;
//         for cell in ordered_cells.iter() {
//             let cell_intids = &hp_map[cell];
//             if cell_intids.len() < MIN_POINTS_PER_CELL {
//                 // println!("Skipping cell {} with only {} points.", cell, cell_intids.len());
//                 break;
//             }
//             let (lon, lat) = nested.center(*cell);
//             let (cell_x, cell_y, cell_z) = lonlat_to_unit_vector(lon, lat);
//             let cell_point = [cell_x, cell_y, cell_z];
//             let cluster = tree.within::<SquaredEuclidean>(&cell_point, cluster_r2);
//             // println!("Cell {}: found {} points in cluster.", cell, cluster.len());
//             let clust = cull_cluster(&detections, &indices, &points, cell_point, &cluster, ic.epoch);
//             if clust.is_none() {
//                 continue;
//             }
//             let clust = clust.unwrap();

//             // let ast_ucty: Vec<f64> = clust.local_intids
//             //     .iter()
//             //     .map(|id| detections[*id as usize].ast_ucty.clone())
//             //     .collect();

//             // let ast_ucty: Vec<Option<f64>> = clust.local_intids
//             //     .iter()
//             //     .map(|id| detections[*id as usize].ast_ucty)
//             //     .collect();

//             // make ast_ucty a vector of 0.15 arcsec in radians
//             let ast_ucty: Vec<f64> = clust.local_intids
//                 .iter()
//                 .map(|_id| (0.15 / 3600.0) * (PI / 180.0))
//                 .collect();

//             let theta_x : Vec<f64> = clust.local_thetas
//                 .iter()
//                 .map(|xy| xy[0])
//                 .collect();
//             let theta_y : Vec<f64> = clust.local_thetas
//                 .iter()
//                 .map(|xy| xy[1])
//                 .collect();
//             let t: Vec<f64> = clust.local_intids
//                 .iter()
//                 .map(|id| detections[*id as usize].epoch - ic.epoch)
//                 .collect();
//             let mxe : Vec<f64> = clust.local_obs_positions
//                 .iter()
//                 .map(|pos| -pos[0])
//                 .collect();
//             let mye : Vec<f64> = clust.local_obs_positions
//                 .iter()
//                 .map(|pos| -pos[1])
//                  .collect();

            
//             let fit_result = iterative_reject(
//                 &t, &mxe, &mye, &theta_x, &theta_y, &ast_ucty, &ast_ucty,
//                 5.0,   // 3σ threshold
//                 MIN_POINTS_PER_CELL     // stop if fewer than 10 points remain
//             );

//             // println!("# Fit result: {:?}", fit_result);
            
//             let fit_final = match fit_result {
//                 Ok(fit_final) => fit_final,
//                 Err(e) => {
//                     continue;
//                 }
//             };

//             let (fit, kept_indices, rejected_indices) = fit_final;

//             let params = fit.params;
//             println!("# IC {} of {}: {}", ii + 1, n_ics, ic.id);
//             println!("# Fit successful for cell {}: {:?}", cell, params);
            
//             let rejected_intids: Vec<usize> = rejected_indices
//                 .iter()
//                 .map(|&idx| clust.local_intids[idx] as usize)
//                 .collect();

//             println!("# intid     detid         t(days)     theta_x(\")  theta_y(\")  x_model(\")  y_model(\")    sig_x     sig_y"); 
//             for i in 0..t.len() {
//                 if rejected_indices.contains(&i) {
//                     continue;
//                 }
//                 let x_model = params[0] + params[1] * t[i] + params[2] * mxe[i];
//                 let y_model = params[3] + params[4] * t[i] + params[2] * mye[i];
//                 let rx = theta_x[i] - x_model;
//                 let pull_x = rx / ast_ucty[i];
//                 let ry = theta_y[i] - y_model;
//                 let pull_y = ry / ast_ucty[i];
//                 let pull = pull_x.abs().max(pull_y.abs());
//                 let intid = clust.local_intids[i];
//                 let detid = &detections[intid as usize].detid;
//                 println!("{:8} {:?} {:12.6}   {:8.4}     {:8.4}    {:8.4}   {:8.4}.  {:8.4}.  {:8.4}", clust.local_intids[i], detid, t[i], theta_x[i]*206265., theta_y[i]*206265., x_model*206265., y_model*206265., pull_x, pull_y);
//             }

//             println!("# detid         epoch        lon(deg)   sig_x(\")   lat(deg)  sig_y(\")  mag");
//             for i in 0..t.len() {
//                 if rejected_indices.contains(&i) {
//                     continue;
//                 }
//                 let intid = clust.local_intids[i];
//                 let detid = detections[intid as usize].detid.clone().unwrap();
//                 let rho_hat = detections[clust.local_intids[i] as usize].rho_hat; 
//                 let lam = rho_hat[1].atan2(rho_hat[0]);
//                 let beta = rho_hat[2].asin();
//                 let epoch = detections[clust.local_intids[i] as usize].epoch;
//                 let mag = detections[clust.local_intids[i] as usize].mag.unwrap_or(-1.0);

//                 // println!("{:12} {:12} {:12.6} {:9.5}   {:5.2}   {:9.5}   {:5.2}  {:13.10}  {:13.10}  {:13.10}   {}   {:.2} {:.2} {}", 
//                 //     ic.id, detid, epoch, ra_deg, ast_ucty[i]*206265., dec_deg, ast_ucty[i]*206265., obs_pos2[0], obs_pos2[1], obs_pos2[2], obscode, mag, mag_sig, filt);
//                 println!("{:12} {:12} {:12.6} {:9.5}   {:5.2}   {:9.5}   {:5.2}   {:.2}", 
//                     ic.id, detid, epoch, lam, ast_ucty[i]*206265., beta, ast_ucty[i]*206265., mag);

//             }
             
  
//         }
//     });
            
//     //};

//     let duration = start_time.elapsed();
//     println!("Time elapsed: {:?}", duration);


//     println!("Matt test program.");