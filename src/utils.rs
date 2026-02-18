use std::collections::HashSet;
use ordered_float::OrderedFloat;
use std::f64::consts::PI;


use nalgebra::{DMatrix, DVector, Matrix3, Matrix3xX, Unit, Vector3, matrix};
use std::collections::HashMap;
use cdshealpix::nested::{get, Layer};

use kiddo::NearestNeighbour;

use crate::Detection;


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

pub fn proj_to_xyz_matrix(r_ref: Vector3<f64>) -> Matrix3<f64> {
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
        -slon0,            -clon0 * slat0,  clon0 * clat0,
         clon0,            -slon0 * slat0,  slon0 * clat0,
         0.0,               clat0,          slat0,
    )
}

pub fn ecliptic_to_equatorial_vector(v_ecl: Vector3<f64>, epsilon: f64) -> Vector3<f64> {
    // Obliquity of the ecliptic (J2000) in radians
    //let epsilon = 23.4392911 * PI / 180.0;
    
    let cos_e = epsilon.cos();
    let sin_e = epsilon.sin();

    // Rotation matrix around the X-axis
    let rot_x = Matrix3::new(
        1.0, 0.0, 0.0,
        0.0, cos_e, -sin_e,
        0.0, sin_e, cos_e,
    );

    rot_x * v_ecl
}

pub fn find_unique_f64_hashset(data: Vec<f64>) -> Vec<f64> {
    let unique_elements: HashSet<OrderedFloat<f64>> = data.into_iter().map(OrderedFloat).collect();

    unique_elements.into_iter().map(|x| x.0).collect()
}


// Function to convert a unit vector (x, y, z) to (lon, lat) in radians
pub fn unit_vector_to_lonlat(x: f64, y: f64, z: f64) -> (f64, f64) {
    let lon = y.atan2(x); // atan2(y, x) for longitude
    let lat = z.asin();    // asin(z) for latitude
    (lon, lat)
}

pub fn lonlat_to_unit_vector(lon: f64, lat: f64) -> (f64, f64, f64) {
    let x = lat.cos() * lon.cos();
    let y = lat.cos() * lon.sin();
    let z = lat.sin();
    (x, y, z)
}



#[derive(Debug)]
pub struct FitResult {
    pub params: [f64; 5],  // a, b, c, d, e
    pub cov: DMatrix<f64>, // 5×5 covariance matrix
    pub sigma: f64,
    pub r2: f64,
    pub chi2: f64,
    pub dof: f64,
}

/// Weighted coupled least squares with explicit error handling.
pub fn coupled_weighted_fit(
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
pub fn iterative_reject(
    t: &[f64],
    f: &[f64],
    g: &[f64],
    x: &[f64],
    y: &[f64],
    sigma_x: &[f64],
    sigma_y: &[f64],
    pull_threshold: f64,
    min_points: usize,
    min_unique_times: usize,
    min_nights: usize,
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

        if t_n.len() < min_nights {
            return Err(format!(
                "Stopped: not enough unique observation times left ({} < {}).",
                t_n.len(),
                min_unique_times
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





#[derive(Debug)]
pub struct Cluster {
    pub ref_epoch: f64,
    pub local_intids: Vec<i32>,
    pub ref_point: [f64; 3],
    pub local_points: Vec<[f64; 3]>,
    pub local_thetas: Vec<[f64; 2]>,
    pub local_obs_positions: Vec<Vector3<f64>>,
    pub unique_cluster_times: Vec<f64>,
    }
    
impl Cluster {
    pub fn from_intids(
        central_intid: i32,
        local_intids: Vec<i32>,
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

        let mut local_thetas: Vec<[f64; 2]> = Vec::new();
        //let mut local_thetas_hash: HashMap<u64, [f64; 2]> = HashMap::new();
        let mut local_obs_positions: Vec<Vector3<f64>> = Vec::new();

        for (p, local_intid) in local_points
            .iter()
            .map(|&p| p)
            .zip(local_intids.iter().map(|&i| i as usize))
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

            local_thetas.push([theta_x, theta_y]);
            //local_thetas_hash.insert(local_intid as u64, [theta_x, theta_y]);
            println!("{:10} {:10.4} {:10.4} {:12.9} {:12.9} {:10.6} {}", local_intid, theta_x*206265., theta_y*206265., longitude, latitude, det.epoch, detid.unwrap());
        }

        /*
        let mut local_tree: KdTree<f64, 2> = KdTree::with_capacity(10);
            local_tree.extend(
                local_thetas
                    .iter()
                    .map(|&p| p)
                    .zip(local_intids.iter().map(|&i| i as u64)),
            );
        */

        Cluster {
            //intid: central_intid,
            ref_epoch: central_detection_epoch,
            local_intids,
            ref_point: central_point,
            local_points,
            local_thetas,
            //local_thetas_hash,
            //local_tree,
            local_obs_positions,
            unique_cluster_times,
        }
    }
}

pub fn vet_tracklet(detections: &Vec<Detection>, intids: &Vec<usize>) -> bool {
    let dets: Vec<&Detection> = intids
        .iter()
        .map(|&intid| &detections[intid as usize])
        .collect();
    true
}

#[inline]
pub fn wrap_0_2pi(x: f64) -> f64 {
    // wraps angle into [0, 2π)
    let y = x % (2.0 * PI);
    if y < 0.0 { y + 2.0 * PI } else { y }
}
/// Convert equatorial (RA, Dec) → ecliptic (λ, β)
/// All angles in **radians**
pub fn equatorial_to_ecliptic(ra: f64, dec: f64, eps: f64) -> (f64, f64) {
    // Obliquity of the ecliptic (J2000)

    // Compute ecliptic longitude (λ) and latitude (β)
    let lambda = (ra.sin() * eps.cos() + dec.tan() * eps.sin()).atan2(ra.cos());
    let beta = (dec.sin() * eps.cos() - dec.cos() * eps.sin() * ra.sin()).asin();

    (lambda, beta)
}

/// Convert ecliptic (lambda, beta) -> equatorial (ra, dec).
/// Inputs/outputs are in radians.
/// Uses obliquity eps (radians).
pub fn ecliptic_to_equatorial(lam: f64, beta: f64, eps: f64) -> (f64, f64) {
    // sin(dec) = sin(beta)*cos(eps) + cos(beta)*sin(eps)*sin(lam)
    let sin_dec = beta.sin() * eps.cos() + beta.cos() * eps.sin() * lam.sin();
    let dec = sin_dec.asin();

    // ra = atan2( cos(beta)*sin(lam)*cos(eps) - sin(beta)*sin(eps), cos(beta)*cos(lam) )
    let y = beta.cos() * lam.sin() * eps.cos() - beta.sin() * eps.sin();
    let x = beta.cos() * lam.cos();
    let ra = wrap_0_2pi(y.atan2(x));

    (ra, dec)
}

pub fn print_type<T>(_: &T) {
    println!("{}", std::any::type_name::<T>());
}

pub fn cull_cluster(detections: &Vec<Detection>, intid_idx_hash: &HashMap<i32, usize>, points: &Vec<[f64; 3]>, point: [f64; 3], cluster: &Vec<NearestNeighbour<f64, i32>>, central_detection_epoch: f64, max_dt: f64, min_unique_times: usize, min_nights: usize) -> Option<Cluster> {

    // Prune detections that are not compatible with the parent/central detection.

    let mut cluster_hash: Vec<i32> = Vec::new();
    let mut min_epoch = f64::MAX;
    let mut max_epoch = f64::MIN;

    for neighbor in cluster.clone().into_iter() {
        
        let det = &detections[neighbor.item as usize];

        let dt = (det.epoch - central_detection_epoch).abs();
        let distance = neighbor.distance;

        if dt > max_dt {
            continue;
        }

        if det.epoch < min_epoch {
            min_epoch = det.epoch;
        }

        if det.epoch > max_epoch {
            max_epoch = det.epoch;
        }

        cluster_hash.push(neighbor.item as i32);

    }

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
    
    if cluster_nights.len() < min_nights {
        return None;
    }
    
    let mut unique_cluster_times = find_unique_f64_hashset(cluster_times.clone());

    if unique_cluster_times.len() < min_unique_times {
        return None;
    }

    // Transform to local tangent plane coordinates.
    let point_vec = Vector3::new(point[0], point[1], point[2]);
    let mat = xyz_to_proj_matrix(point_vec);

    cluster_hash.sort();
    let mut local_intids: Vec<i32> = cluster_hash.iter().cloned().collect();

    let mut local_points: Vec<[f64; 3]> = local_intids
                .iter()
                .map(|intid| points[intid_idx_hash[intid]].clone())
                .collect();

    let mut local_thetas: Vec<[f64; 2]> = Vec::new();
    let mut local_obs_positions: Vec<Vector3<f64>> = Vec::new();
    for (p, local_intid) in local_points
        .iter()
        .map(|&p| p)
        .zip(local_intids.iter().map(|&i| i))
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
        
    }

    return Some(Cluster {
        //intid: central_intid,
        ref_epoch: central_detection_epoch,
        local_intids: local_intids,
        ref_point: point,
        local_points: local_points,
        local_thetas: local_thetas,
        local_obs_positions: local_obs_positions,
        unique_cluster_times: unique_cluster_times,
    })
}



pub fn sum_hp_map_neighbors_v2(hp_map_counts: &HashMap<u64, usize>, nested: &Layer, radius: f64)-> HashMap<u64, usize> {

    let mut hp_sum_map: HashMap<u64, usize> = HashMap::new();

    for (key, values) in hp_map_counts {
        if *values < 2 {
            continue; // skip low-population cells
        }
        let res = nested.kth_neighbourhood(*key, 1u32);
        //let res: Vec<u64> = vec![*key];

        for n in res.iter() {
            let v = hp_map_counts.get(&n);
            let neighbor_value = hp_map_counts.get(&n);
            if let Some(neighbor_value) = neighbor_value {
                let neighbor_count: usize = *neighbor_value;
                hp_sum_map.entry(*key).and_modify(|e| *e += neighbor_count).or_insert(neighbor_count);
            }
        }
    }
    hp_sum_map
}

/// Find non-strict local maxima (peaks) in a HEALPix map.
/// A tile is considered a peak if its value >= all its neighbors and value >= threshold.
/// Returns list of (tile, value) pairs, sorted by decreasing value.
pub fn find_non_strict_peaks(
    map: &HashMap<u64, usize>,
    depth: u8,
    threshold: usize,
) -> Vec<(u64, usize)> {
    let nested = get(depth);
    let mut peaks = Vec::new();

    for (&tile, &value) in map.iter() {
        if value < threshold {
            continue; // skip low values
        }

        let mut is_peak = true;

        for neighbor in nested.kth_neighbourhood(tile, 1) {
            if neighbor == tile {
                continue; // skip self
            }

            let neighbor_value = map.get(&neighbor).copied().unwrap_or(0);
            if neighbor_value > value {
                is_peak = false;
                break;
            }
        }

        if is_peak {
            peaks.push((tile, value));
        }
    }

    // Sort by decreasing value
    peaks.sort_by(|a, b| b.1.cmp(&a.1));

    peaks
}

pub fn gather_cluster_data(clust: &Cluster, detections: &Vec<Detection>, central_detection_epoch: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    // Gather the data for the detections in this cluster, to be used for model fitting.
    let ast_ucty: Vec<f64> = clust.local_intids
        .iter()
        .map(|id| detections[*id as usize].ast_ucty.unwrap_or(0.2))
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
        .map(|id| detections[*id as usize].epoch - central_detection_epoch)
        .collect();
    let mxe : Vec<f64> = clust.local_obs_positions
        .iter()
        .map(|pos| -pos[0])
        .collect();
    let mye : Vec<f64> = clust.local_obs_positions
        .iter()
        .map(|pos| -pos[1])
            .collect();
    let xe : Vec<f64> = clust.local_obs_positions
        .iter()
        .map(|pos| pos[0])
        .collect();
    let ye : Vec<f64> = clust.local_obs_positions
        .iter()
        .map(|pos| pos[1])
            .collect();
    (t, theta_x, theta_y, ast_ucty, mxe, mye, xe, ye)
}