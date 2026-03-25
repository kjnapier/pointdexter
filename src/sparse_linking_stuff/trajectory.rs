// src/trajectory.rs

use std::collections::{HashMap, HashSet};
use bloomfilter::Bloom;
use serde::{Serialize, Deserialize};

use crate::sparse_linking_stuff::Config;

#[derive(Debug, Serialize, Deserialize)]
pub struct SavedTrajectory {
    pub detection_objids: Vec<String>,
    pub ic_params: [f64; 4]  // use native type
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetaPoint {
    pub delta_phi: f64,
    pub delta_theta: f64,
    pub epoch: f64,
    pub magnitude: f64,
    pub night: i64,
    pub detection_count: usize,
    pub original_detections: Vec<RawDetection>,
    pub expnums: Vec<i64>, // Added to track unique exposures
}

#[derive(Debug, Clone, PartialEq)]
pub struct RawDetection {
    pub id: usize,
    pub delta_phi: f64,
    pub delta_theta: f64,
    pub epoch: f64,
    pub magnitude: f64,
    pub night: i64,
    pub expnum: i64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Trajectory {
    pub detection_ids: Vec<usize>,
    pub center_id: u64, // intid of the cluster center
}

#[derive(Debug)]
pub struct TrajectoryFilter {
    accepted_links: Vec<HashSet<usize>>,
    failed_minimals: Bloom<Vec<u64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FakeTrajectoryRecord {
    pub fakeid: String,
    pub detection_ids: Vec<usize>,
    pub ic_id: usize,  // use native type
}

impl Default for TrajectoryFilter {
    fn default() -> Self {
        Self::new()
    }
}

impl TrajectoryFilter {
    pub fn new() -> Self {
        // 5% false positive rate for up to 1 billion failed trajectories
        let expected_trajectories = 1_000_000_000;
        let fp_rate = 0.05;

        Self {
            accepted_links: Vec::new(),
            failed_minimals: Bloom::new_for_fp_rate(expected_trajectories, fp_rate).expect("Failed to create bloom filter"),
        }
    }

    pub fn add_accepted(&mut self, ids: &HashSet<usize>) {
        self.accepted_links.push(ids.clone());
    }

    pub fn add_failed_minimal(&mut self, ids: &HashSet<usize>) {
        let mut sorted_ids: Vec<u64> = ids.iter().map(|&x| x as u64).collect();
        sorted_ids.sort_unstable();
        self.failed_minimals.set(&sorted_ids);
    }

    pub fn should_skip(&self, candidate: &HashSet<usize>) -> bool {
        // 1. Skip if overlaps with any accepted link by ≥ 15 detections
        if self.accepted_links.iter().any(|link| link.intersection(candidate).count() >= 15) {
            return true;
        }

        // 2. Skip if this was previously failed (exact match, but this shouldn't ever happen since we are deduplicating trajectories, but across batches this can happen)
        let mut sorted_ids: Vec<u64> = candidate.iter().map(|&x| x as u64).collect();
        sorted_ids.sort_unstable();
        if self.failed_minimals.check(&sorted_ids) {
            return true;
        }

        false
    }
}

pub fn compute_vector(p1: &MetaPoint, p2: &MetaPoint) -> [f64; 2] {
    [p2.delta_phi - p1.delta_phi, p2.delta_theta - p1.delta_theta]
}

pub fn compute_speed(p1: &MetaPoint, p2: &MetaPoint) -> f64 {
    let dt = p2.epoch - p1.epoch;
    if dt.abs() < std::f64::EPSILON { return 0.0; }
    let disp = compute_vector(p1, p2);
    let dist = (disp[0].powi(2) + disp[1].powi(2)).sqrt();
    dist / dt
}

// Create what are essentially "tracklets" in d_phi, d_theta space by grouping together detections from the same night that are close in space and time.
// Within the meta-point, keep the average position, epoch, magnitude, and the original detections that were grouped together.
pub fn create_meta_points(detections: &[RawDetection], max_distance: f64, max_time_span: f64) -> Vec<MetaPoint> {
    let mut nights: HashMap<i64, Vec<RawDetection>> = HashMap::new();
    for det in detections.iter().cloned() {
        nights.entry(det.night).or_default().push(det.clone());
    }

    let mut meta_points = Vec::new();
    for (night, dets) in nights {
        let mut processed = HashSet::new();
        for (i, det) in dets.iter().enumerate() {
            if processed.contains(&i) { continue; }


            let mut group = vec![det.clone()];
            processed.insert(i);


            for (j, other) in dets.iter().enumerate() {
                if i == j || processed.contains(&j) { continue; }
                let dist = ((det.delta_phi - other.delta_phi).powi(2)
                    + (det.delta_theta - other.delta_theta).powi(2)).sqrt();
                let time_diff = (det.epoch - other.epoch).abs();
                if dist <= max_distance && time_diff <= max_time_span {
                    group.push(other.clone());
                    processed.insert(j);
                }
            }

            // Collect unique expnums from THIS group only
            let expnums: Vec<i64> = group.iter()
                .map(|d| d.expnum)
                .collect::<HashSet<_>>()  // Deduplicate
                .into_iter()
                .collect();


            let len = group.len() as f64;
            let avg_phi = group.iter().map(|d| d.delta_phi).sum::<f64>() / len;
            let avg_theta = group.iter().map(|d| d.delta_theta).sum::<f64>() / len;
            let avg_epoch = group.iter().map(|d| d.epoch).sum::<f64>() / len;
            let avg_mag = group.iter().map(|d| d.magnitude).sum::<f64>() / len;
            meta_points.push(MetaPoint {
                delta_phi: avg_phi,
                delta_theta: avg_theta,
                epoch: avg_epoch,
                magnitude: avg_mag,
                night,
                detection_count: group.len(),
                original_detections: group,
                expnums: expnums.clone(),

            });
        }
    }
    meta_points
}

pub fn speeds_within_tolerance(speed0: f64, speed1: f64, speed_tol: f64,) -> bool {
    (speed1 - speed0).abs() < speed_tol
}

pub fn angle_between_vectors(v1: [f64; 2], v2: [f64; 2]) -> f64 {
    let dot = v1[0] * v2[0] + v1[1] * v2[1];
    let norm1 = (v1[0].powi(2) + v1[1].powi(2)).sqrt();
    let norm2 = (v2[0].powi(2) + v2[1].powi(2)).sqrt();
    if norm1 == 0.0 || norm2 == 0.0 { return 0.0; }
    let cos_theta = (dot / (norm1 * norm2)).clamp(-1.0, 1.0);
    cos_theta.acos().to_degrees()
}

pub fn is_valid_candidate(p1: &MetaPoint, p2: &MetaPoint, p3: &MetaPoint, rel_speed_tol: f64, max_angle_deg: f64,) -> bool {
    if !(p1.epoch < p2.epoch && p2.epoch < p3.epoch) { return false; }
    let speed0 = compute_speed(p1, p2);
    let speed1 = compute_speed(p2, p3);
    if !speeds_within_tolerance(speed0, speed1, rel_speed_tol) { return false; }
    let angle = angle_between_vectors(
        compute_vector(p1, p2),
        compute_vector(p2, p3),
    );
    angle <= max_angle_deg
}

pub fn build_trajectory_graph(detections: &[MetaPoint]) -> HashMap<usize, Vec<usize>> {
    let mut graph: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..detections.len() {
        for j in 0..detections.len() {
            if detections[j].epoch > detections[i].epoch {
                graph.entry(i).or_default().push(j);
            }
        }
    }
    graph
}

pub fn is_valid_trajectory(traj: &[MetaPoint], config: &Config,) -> bool {
    let total_detections: usize = traj.iter().map(|p| p.detection_count).sum();
    let unique_nights: usize = traj.iter().map(|p| p.night).collect::<HashSet<_>>().len();
    let detection_counts: Vec<usize> = traj.iter().map(|p| p.detection_count).collect();

    // Calculate duration (max epoch - min epoch)
    let epochs: Vec<f64> = traj.iter().map(|p| p.epoch).collect();
    let min_epoch = epochs.iter().copied().fold(f64::INFINITY, f64::min);
    let max_epoch = epochs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let duration = max_epoch - min_epoch;

    // Determine how many unique expnums it has
    let unique_expnums: usize = traj.iter()
        .flat_map(|p| p.expnums.iter())
        .collect::<HashSet<_>>()
        .len();

    // Our has_heavy check is hardcoded for now
    let has_heavy = detection_counts.iter().any(|&d| d >= 12);
    let support_count = detection_counts.iter().filter(|&&d| d > 1).count();

    // Require at least one metapoint with >= 2 detections
    let has_tracklet = detection_counts.iter().any(|&d| d >= 2);

    let valid = total_detections >= config.min_detections 
        && unique_nights >= config.min_nites 
        && has_tracklet 
        && duration >= config.min_duration
        && unique_expnums >= config.min_detections;
    let not_lopsided = !(has_heavy && support_count <= 1);
    valid && not_lopsided
}

pub fn find_all_valid_trajectories(graph: &HashMap<usize, Vec<usize>>, detections: &[MetaPoint], rel_speed_tol: f64, max_angle_deg: f64, center_idx: usize, config: &Config,) -> Vec<Vec<MetaPoint>> {
    let mut results = Vec::new();

    fn dfs(
        current_path: Vec<usize>,
        graph: &HashMap<usize, Vec<usize>>,
        detections: &[MetaPoint],
        results: &mut Vec<Vec<MetaPoint>>,
        rel_speed_tol: f64,
        max_angle_deg: f64,
        center_idx: usize,
        config: &Config,
    ) {
        let contains_center = current_path.contains(&center_idx);
        let traj: Vec<MetaPoint> = current_path.iter().map(|&i| detections[i].clone()).collect();
    
        if contains_center {
            let total_dets: usize = traj.iter().map(|p| p.detection_count).sum();
            if total_dets >= config.min_detections && is_valid_trajectory(&traj, config) {
                results.push(traj.clone());
            }
        }
    
        if let Some(neighbors) = graph.get(current_path.last().unwrap()) {
            for &neighbor in neighbors {
                if current_path.contains(&neighbor) {
                    continue;
                }
    
                if current_path.len() >= 2 {
                    let p1 = &detections[current_path[current_path.len() - 2]];
                    let p2 = &detections[current_path[current_path.len() - 1]];
                    let p3 = &detections[neighbor];
                    if !is_valid_candidate(p1, p2, p3, rel_speed_tol, max_angle_deg) {
                        continue;
                    }
                }
    
                let mut next_path = current_path.clone();
                next_path.push(neighbor);
                dfs(next_path, graph, detections, results, rel_speed_tol, max_angle_deg, center_idx, config);
            }
        }
    }
    

    for i in 0..detections.len() {
        dfs(vec![i], graph, detections, &mut results, rel_speed_tol, max_angle_deg, center_idx, config);
    }

    results
}

pub fn graph_based_trajectory_search(detections: &[MetaPoint], rel_speed_tol: f64, max_angle_deg: f64, config: &Config) -> Vec<Vec<MetaPoint>> {
    let graph = build_trajectory_graph(detections);
    let center_idx = find_center_index(detections)
        .expect("No center MetaPoint found with raw detection at origin.");
    find_all_valid_trajectories(&graph, detections, rel_speed_tol, max_angle_deg, center_idx, &config)
}

pub fn find_center_index(detections: &[MetaPoint]) -> Option<usize> {
    detections.iter().enumerate().find_map(|(i, mp)| {
        if mp.original_detections.iter().any(|d| d.delta_phi == 0.0 && d.delta_theta == 0.0) {
            Some(i)
        } else {
            None
        }
    })
}


/// Will be moved to utils
pub fn angular_separation(ref_vec: [f64; 3], vec: [f64; 3]) -> [f64; 2] {
    const RADS_TO_ARCSEC: f64 = 206265.0;
    let phi0 = ref_vec[1].atan2(ref_vec[0]);
    let theta0 = ref_vec[2].asin();
    let phi1 = vec[1].atan2(vec[0]);
    let theta1 = vec[2].asin();
    let mut dphi = (phi1 - phi0) * RADS_TO_ARCSEC;
    dphi *= theta0.cos();
    let dtheta = (theta1 - theta0) * RADS_TO_ARCSEC;
    [dphi, dtheta]
}