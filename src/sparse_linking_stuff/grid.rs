use rayon::prelude::*;
use std::fmt;
use crate::sparse_linking_stuff::orbit::{Orbit, find_max_separation_cheat};
use crate::initial_condition::{InitialCondition, InputFormat};
//use crate::sparse_linking_stuff::InputFormat;

use ordered_float::OrderedFloat;
use std::collections::HashSet;
use spacerocks::Time;
use std::sync::Arc;
use std::sync::Mutex;
use std::collections::HashMap;
use crate::detection::Detection;
use crate::sync::sync_detection_to_orbit as construct_orbit;
use crate::sparse_linking_stuff::trajectory::Trajectory;
use crate::sparse_linking_stuff::utils::angular_separation;
use spacerocks::SpaceRock;

const ARCSEC_PER_RAD: f64 = 3600.0 * 180.0 / std::f64::consts::PI;
const RADS_TO_ARCSEC: f64 = 206265.0;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GridCacheKey {
    pub ic_key: IcBoundsKey,
    pub epsilon: OrderedFloat<f64>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum IcBoundsKey {
    QEF {
        q_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        e_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        f_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        psi_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        epoch: OrderedFloat<f64>,
    },
    KEP {
        a_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        e_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        f_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        psi_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        epoch: OrderedFloat<f64>,
    },
    KEV {
        r_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        vr_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        vo_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        psi_bounds: (OrderedFloat<f64>, OrderedFloat<f64>),
        epoch: OrderedFloat<f64>,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct Cell {
    pub bounds: [(f64, f64); 4],
    pub midpoints: [f64; 4],
}

impl Cell {
    pub fn new(bounds: [(f64, f64); 4]) -> Self {
        let mut midpoints = [0.0; 4];
        for i in 0..4 {
            midpoints[i] = (bounds[i].0 + bounds[i].1) / 2.0;
        }
        Cell { bounds, midpoints }
    }

    pub fn midpoint(&self) -> &[f64; 4] {
        &self.midpoints
    }

    pub fn split(&self, axis: usize) -> (Self, Self) {
        let mut bounds_1 = self.bounds;
        let mut bounds_2 = self.bounds;

        let (min, max) = self.bounds[axis];
        let mid = (min + max) / 2.0;

        bounds_1[axis] = (min, mid);
        bounds_2[axis] = (mid, max);

        (Cell::new(bounds_1), Cell::new(bounds_2))
    }


    pub fn from_initial_condition(ic: &InitialCondition) -> Option<Self> {

        Some(Cell::new([
            (ic.r_min?, ic.r_max?),
            (ic.vr_min?, ic.vr_max?),
            (ic.vo_min?, ic.vo_max?),
            (ic.inc_min?, ic.inc_max?),
        ]))
    }
}

impl fmt::Display for Cell {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.bounds)
    }
}
fn get_orbit(
    values: &[f64; 4],
    format: InputFormat,
    t: f64,
    mu: f64,
) -> Orbit {
    let epoch = Time::new(t, "utc", "jd").unwrap();

    let ic = match format {
        InputFormat::KEV => InitialCondition::from_spherical(
            "temp".to_string(),
            values[0], // r
            values[1], // vr
            values[2], // vo
            values[3], // inc (used as psi)
            0,
            epoch,
            mu,
            None, None, None, None, None, None, None, None,
        ).unwrap(),
        InputFormat::QEF => InitialCondition::from_elements(
            "temp".to_string(),
            values[0], // q
            values[1], // e
            values[3], // inc
            values[2], // true_anomaly
            0,
            epoch,
            mu,
            None, None, None, None, None, None, None, None,
        ).unwrap(),
        InputFormat::KEP => InitialCondition::from_keplerian(
            "temp".to_string(),
            values[0], // q
            values[1], // e
            values[3], // inc
            values[2], // mean_anomaly
            0,
            epoch,
            mu,
            None, None, None, None, None, None, None, None,
        ).unwrap(),
    };

    Orbit::new(
        "temp".to_string(),
        ic.q,
        ic.e,
        ic.true_anomaly,
        ic.inc,  // psi in Orbit corresponds to inc in InitialCondition
        t,
        mu,
    )
}

fn get_separation(o1: &Orbit, o2: &Orbit, t_bounds: (f64, f64)) -> f64 {
    find_max_separation_cheat(o1, o2, t_bounds) * ARCSEC_PER_RAD
}

pub fn calculate_axis_separation(cell: &Cell, axis: usize, t: f64, t_bounds: (f64, f64), mu: f64, format: InputFormat) -> f64 {
    let combos = vec![
        vec![0, 0, 0], vec![0, 0, 1], vec![0, 1, 0], vec![0, 1, 1],
        vec![1, 0, 0], vec![1, 0, 1], vec![1, 1, 0], vec![1, 1, 1]
    ];
    let axes: HashSet<usize> = (0..4).collect();
    let remaining_axes: Vec<usize> = axes.difference(&HashSet::from([axis])).cloned().collect();

    combos.par_iter()
        .map(|combo| {
            let mut orbits = Vec::new();
            let axis_bound = cell.bounds[axis];
            for &v in &[axis_bound.0, axis_bound.1] {
                let mut values = [0.0; 4];
                for (&idx, &ax) in combo.iter().zip(remaining_axes.iter()) {
                    let bound = cell.bounds[ax];
                    values[ax] = if idx == 0 { bound.0 } else { bound.1 };
                }
                values[axis] = v;
                // orbits.push(get_orbit(values[0], values[1], values[2], values[3], t, mu));
                orbits.push(get_orbit(&values, format, t, mu));

            }
            get_separation(&orbits[0], &orbits[1], t_bounds)
        })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(0.0)
}

fn split_cell(cell: &Cell, t: f64, t_bounds: (f64, f64), mu: f64, epsilon: f64, format: InputFormat) -> Option<Vec<Cell>> {

    let mut max_separation = 0.0;
    let mut max_separation_axis = 0;
    let mut should_split = false;

    for axis in 0..4 {

        let separation = calculate_axis_separation(cell, axis, t, t_bounds, mu, format);

        if separation > epsilon {
            should_split = true;
        }
        if separation > max_separation {
            max_separation = separation;
            max_separation_axis = axis;
        }
    }

    if should_split {
        let (cell1, cell2) = cell.split(max_separation_axis);
        Some(vec![cell1, cell2])
    } else {
        None
    }
}

pub fn process_cells_parallel(
    cells: Vec<Cell>,
    t: f64,
    t_bounds: (f64, f64),
    mu: f64,
    epsilon: f64,
    format: InputFormat, // NEW
) -> Vec<Cell> {
    cells
        .into_par_iter()
        .flat_map(|cell| {
            split_cell(&cell, t, t_bounds, mu, epsilon, format)
                .unwrap_or_else(|| vec![cell])
        })
        .collect()
}

pub fn generate_adaptive_grid(
    initial_bounds: [(f64, f64); 4],
    mu: f64,
    t_bounds: (f64, f64),
    epsilon: f64,
    epoch: Time,
    format: InputFormat, // NEW
) -> Vec<Cell> {

    let epoch_val = epoch.epoch;
    let initial_cell = Cell::new(initial_bounds);
    let mut cells = vec![initial_cell];
    let mut prev_len = 0;

    println!("Initial bounds: {:?}", initial_bounds);

    while prev_len != cells.len() {
        prev_len = cells.len();
        cells = cells
            .into_par_iter()
            .flat_map(|c| split_cell(&c, epoch_val, t_bounds, mu, epsilon, format).unwrap_or_else(|| vec![c]))
            .collect();
    }

    cells
}


/// Keep n pairs


use std::cmp::Ordering;
fn refine_single_trajectory(
    traj: &Trajectory,
    original_ic: &InitialCondition,
    grid_cache: &Arc<Mutex<HashMap<GridCacheKey, Vec<Cell>>>>,
    det_map: &HashMap<u64, Detection>,
    epsilon: f64,
    t_bounds: (f64, f64),
    max_keep: usize,
) -> Vec<(Trajectory, InitialCondition)> {
    let t = (t_bounds.0 + t_bounds.1) / 2.0;

    let ic_key = match original_ic.to_bounds_key() {
        Some(k) => k,
        None => return Vec::new(),
    };
    let grid_key = GridCacheKey {
        ic_key,
        epsilon: OrderedFloat(epsilon),
    };

    // Try to read from cache without holding the lock for long
    let maybe_cells = {
        let cache = grid_cache.lock().unwrap();
        cache.get(&grid_key).cloned()
    };

    let cells = if let Some(cells) = maybe_cells {
        cells
    } else {
       
        let mut computed = vec![match Cell::from_initial_condition(original_ic) {
            Some(cell) => cell,
            None => return Vec::new(),
        }];
        let mut prev_len = 0;
        while prev_len != computed.len() {
            prev_len = computed.len();
            computed = process_cells_parallel(
                computed,
                t,
                t_bounds,
                original_ic.mu,
                epsilon,
                InputFormat::KEV, // Assuming KEV format for grid generation; adjust if needed
            );
        }

        // Insert into cache
        {
            let mut cache = grid_cache.lock().unwrap();
            cache.insert(grid_key.clone(), computed.clone());
        }

        computed
    };

    let center_det = match det_map.get(&(traj.center_id as u64)) {
        Some(d) => d,
        None => return Vec::new(),
    };
    let ids = &traj.detection_ids;

    let mut candidates: Vec<(f64, InitialCondition)> = Vec::new();

    for rc in &cells {
        let mid = rc.midpoint();
        let refined_ic = InitialCondition::from_spherical(
            "temp".to_string(),
            mid[0], // r
            mid[1], // vr
            mid[2], // vo
            mid[3], // inc
            0,
            Time::new(original_ic.epoch, "utc", "jd").unwrap(),
            original_ic.mu,
            Some(rc.bounds[0].0), Some(rc.bounds[0].1),
            Some(rc.bounds[1].0), Some(rc.bounds[1].1),
            Some(rc.bounds[2].0), Some(rc.bounds[2].1),
            Some(rc.bounds[3].0), Some(rc.bounds[3].1),
        ).unwrap();

        // Propagate the center detection
        let center_pt = match construct_orbit(center_det, &refined_ic) {
            Some(v) => v,
            None => continue,
        };

        let mut total_sep = 0.0;
        let mut all_within = true;

        for &id in ids {
            let det = match det_map.get(&(id as u64)) {
                Some(d) => d,
                None => {
                    all_within = false;
                    break;
                }
            };

            let pt = match construct_orbit(det, &refined_ic) {
                Some(v) => v,
                None => {
                    all_within = false;
                    break;
                }
            };

            let [dphi, dtheta] = angular_separation(center_pt, pt);
            let sep = (dphi * dphi + dtheta * dtheta).sqrt();


            // Note: +0.5 margin 
            if sep > epsilon + 0.5 {
                all_within = false;
                break;
            }

            total_sep += sep;
        }

        if all_within {
            let avg_sep = total_sep / ids.len() as f64;
            candidates.push((avg_sep, refined_ic));
        }
    }

    if candidates.is_empty() {
        return Vec::new();
    }

    // Sort by average separation (ascending) and keep up to `max_keep`
    candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Greater));
    candidates
        .into_iter()
        .take(max_keep)
        .map(|(_, ic)| (traj.clone(), ic))
        .collect()
}

pub fn refine_with_grid_cache_sparse(
    trajs: &[(Trajectory, InitialCondition)],
    grid_cache: &Arc<Mutex<HashMap<GridCacheKey, Vec<Cell>>>>,
    det_map: &HashMap<u64, Detection>,
    cluster_radius: f64, // unused but kept for compatibility
    epsilon: f64,
    t_bounds: (f64, f64),
    max_keep: usize,
) -> Vec<(Trajectory, InitialCondition)> {
    println!("Refining trajectories with epsilon: {}", epsilon);

    // For each trajectory, return up to `max_keep` candidate pairs, then flatten
    let lists: Vec<Vec<(Trajectory, InitialCondition)>> = trajs
        .par_iter()
        .map(|(traj, original_ic)| {
            refine_single_trajectory(
                traj,
                original_ic,
                grid_cache,
                det_map,
                epsilon,
                t_bounds,
                max_keep,
            )
        })
        .collect();

    lists.into_iter().flatten().collect()
}


//// NOTE: Below this point is a graveyard of older versions of the refinement code. Explore at your own peril.








































// Keep only the "tightest"


// fn refine_single_trajectory(
//     traj: &Trajectory,
//     original_ic: &InitialCondition,
//     grid_cache: &Arc<Mutex<HashMap<GridCacheKey, Vec<Cell>>>>,
//     det_map: &HashMap<u64, Detection>,
//     epsilon: f64,
//     t_bounds: (f64, f64),
//     // fakes: &[SpaceRock],
// ) -> Option<(Trajectory, InitialCondition)> {
//     let t = (t_bounds.0 + t_bounds.1) / 2.0;

//     let ic_key = original_ic.to_bounds_key()?;
//     let grid_key = GridCacheKey { ic_key, epsilon: OrderedFloat(epsilon) };

//     // Try read from cache without holding lock long
//     let maybe_cells = {
//         let cache = grid_cache.lock().unwrap();
//         cache.get(&grid_key).cloned()
//     };

//     let cells = if let Some(cells) = maybe_cells {
//         cells
//     } else {
//         // Compute expensive cells outside lock
//         let mut computed = vec![Cell::from_initial_condition(original_ic)?];
//         let mut prev_len = 0;
//         while prev_len != computed.len() {
//             prev_len = computed.len();
//             computed = process_cells_parallel(
//                 computed, t, t_bounds, original_ic.mu, epsilon, original_ic.format.clone(),
//             );
//         }

//         // Insert into cache
//         {
//             let mut cache = grid_cache.lock().unwrap();
//             cache.insert(grid_key.clone(), computed.clone());
//         }

//         computed
//     };

//     let center_det = det_map.get(&(traj.center_id as u64))?;
//     let ids = &traj.detection_ids;

//     let mut best_ic = None;
//     let mut min_avg_sep = f64::MAX;

//     for rc in &cells {
//         let mid = rc.midpoint();
//         let refined_ic = InitialCondition::from_params(
//             0, original_ic.format,
//             Some(rc.bounds[0].0), Some(rc.bounds[0].1), mid[0],
//             Some(rc.bounds[1].0), Some(rc.bounds[1].1), mid[1],
//             Some(rc.bounds[2].0), Some(rc.bounds[2].1), mid[2],
//             Some(rc.bounds[3].0), Some(rc.bounds[3].1), mid[3],
//             original_ic.epoch.clone(), original_ic.mu,
//         );

//         let (center_pt, _) = match construct_orbit(center_det, &refined_ic, Some(false)) {
//             Some(v) => v,
//             None => continue,
//         };

//         let mut total_sep = 0.0;
//         let mut all_within = true;

//         for &id in ids {
//             let det = match det_map.get(&(id as u64)) {
//                 Some(d) => d,
//                 None => {
//                     all_within = false;
//                     break;
//                 }
//             };

//             let (pt, _) = match construct_orbit(det, &refined_ic, Some(false)) {
//                 Some(v) => v,
//                 None => {
//                     all_within = false;
//                     break;
//                 }
//             };

//             let [dphi, dtheta] = angular_separation(center_pt, pt);
//             let sep = (dphi.powi(2) + dtheta.powi(2)).sqrt();

//             if sep > epsilon + 0.2 {
//                 all_within = false;
//                 break;
//             }

//             total_sep += sep;
//         }

//         if all_within {
//             let avg_sep = total_sep / ids.len() as f64;
//             if avg_sep < min_avg_sep {
//                 min_avg_sep = avg_sep;
//                 best_ic = Some(refined_ic);
//             }
//         }
//     }

//     best_ic.map(|ic| (traj.clone(), ic))
// }



// pub fn refine_with_grid_cache_sparse(
//     trajs: &[(Trajectory, InitialCondition)],
//     grid_cache: &Arc<Mutex<HashMap<GridCacheKey, Vec<Cell>>>>,
//     det_map: &HashMap<u64, Detection>,
//     cluster_radius: f64,
//     epsilon: f64,
//     t_bounds: (f64, f64),
//     // fakes: &[SpaceRock],
// ) -> Vec<(Trajectory, InitialCondition)> {
//     println!("Refining trajectories with epsilon: {}", epsilon);
//     trajs.par_iter()
//         .filter_map(|(traj, original_ic)| {
//             refine_single_trajectory(traj, original_ic, grid_cache, det_map, epsilon, t_bounds)
//         })
//         .collect()
// }











// Version to keep all valid trajectories that survive

// pub fn refine_single_trajectory(
//     traj: &Trajectory,
//     original_ic: &InitialCondition,
//     grid_cache: &Arc<Mutex<HashMap<GridCacheKey, Vec<Cell>>>>,
//     det_map: &HashMap<u64, Detection>,
//     epsilon: f64,
//     t_bounds: (f64, f64),
// ) -> Option<Vec<(Trajectory, InitialCondition)>> {
//     let t = (t_bounds.0 + t_bounds.1) / 2.0;

//     let ic_key = original_ic.to_bounds_key()?;
//     let grid_key = GridCacheKey {
//         ic_key,
//         epsilon: OrderedFloat(epsilon),
//     };

//     // Try to read from cache
//     let maybe_cells = {
//         let cache = grid_cache.lock().unwrap();
//         cache.get(&grid_key).cloned()
//     };

//     let cells = if let Some(cells) = maybe_cells {
//         cells
//     } else {
//         // Compute cells outside lock
//         let mut computed = vec![Cell::from_initial_condition(original_ic)?];
//         let mut prev_len = 0;
//         while prev_len != computed.len() {
//             prev_len = computed.len();
//             computed = process_cells_parallel(
//                 computed,
//                 t,
//                 t_bounds,
//                 original_ic.mu,
//                 epsilon,
//                 original_ic.format.clone(),
//             );
//         }

//         // Insert into cache
//         {
//             let mut cache = grid_cache.lock().unwrap();
//             cache.insert(grid_key.clone(), computed.clone());
//         }

//         computed
//     };

//     let center_det = det_map.get(&(traj.center_id as u64))?;
//     let ids = &traj.detection_ids;

//     let mut valid_results = Vec::new();

//     for rc in &cells {
//         let mid = rc.midpoint();
//         let refined_ic = InitialCondition::from_params(
//             0,
//             original_ic.format,
//             Some(rc.bounds[0].0),
//             Some(rc.bounds[0].1),
//             mid[0],
//             Some(rc.bounds[1].0),
//             Some(rc.bounds[1].1),
//             mid[1],
//             Some(rc.bounds[2].0),
//             Some(rc.bounds[2].1),
//             mid[2],
//             Some(rc.bounds[3].0),
//             Some(rc.bounds[3].1),
//             mid[3],
//             original_ic.epoch.clone(),
//             original_ic.mu,
//         );

//         let (center_pt, _) = match construct_orbit(center_det, &refined_ic, Some(false)) {
//             Some(v) => v,
//             None => continue,
//         };

//         let mut total_sep = 0.0;
//         let mut all_within = true;

//         for &id in ids {
//             let det = match det_map.get(&(id as u64)) {
//                 Some(d) => d,
//                 None => {
//                     all_within = false;
//                     break;
//                 }
//             };

//             let (pt, _) = match construct_orbit(det, &refined_ic, Some(false)) {
//                 Some(v) => v,
//                 None => {
//                     all_within = false;
//                     break;
//                 }
//             };

//             let [dphi, dtheta] = angular_separation(center_pt, pt);
//             let sep = (dphi.powi(2) + dtheta.powi(2)).sqrt();

//             if sep > epsilon + 0.5 {
//             // if sep > epsilon {
//                 all_within = false;
//                 break;
//             }

//             total_sep += sep;
//         }

//         if all_within {
//             valid_results.push((traj.clone(), refined_ic));
//         }
//     }

//     Some(valid_results)
// }


// pub fn refine_with_grid_cache_sparse(
//     trajs: &[(Trajectory, InitialCondition)],
//     grid_cache: &Arc<Mutex<HashMap<GridCacheKey, Vec<Cell>>>>,
//     det_map: &HashMap<u64, Detection>,
//     cluster_radius: f64,
//     epsilon: f64,
//     t_bounds: (f64, f64),
// ) -> Vec<(Trajectory, InitialCondition)> {
//     println!("Refining trajectories with epsilon: {}", epsilon);

//     trajs
//         .par_iter()
//         .filter_map(|(traj, original_ic)| {
//             refine_single_trajectory(traj, original_ic, grid_cache, det_map, epsilon, t_bounds)
//         })
//         .flatten()
//         .collect()
// }
























// Keep based on spread & linearity

// use rayon::prelude::*;
// use std::cmp::Ordering;
// use std::collections::hash_map::Entry;

// use std::hash::{Hash, Hasher};
// use crate::traj_r_refinement::{optimize_r_and_refine_ic, r2_linear, mean, variance};


// // … existing imports (ConstructOrbit, Detection, Trajectory, Cell, etc.) …

// /// Compute φ and θ from a 3D position vector.
// #[inline]
// fn phi_theta_from_point(pt: &[f64; 3]) -> (f64, f64) {
//     let phi = pt[1].atan2(pt[0]); // y.atan2(x)
//     let r = (pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]).sqrt();
//     let theta = (pt[2] / r).asin();
//     (phi, theta)
// }

// /// Compute spread (variance of φ + variance of θ) and linearity (mean R² of φ and θ vs. epochs)
// /// for a given (traj, ic) pair. Returns None if fewer than 3 valid detections.
// fn compute_raw_metrics(
//     traj: &Trajectory,
//     ic: &InitialCondition,
//     det_map: &HashMap<u64, Detection>,
// ) -> Option<(f64, f64)> {
//     // Gather detections in order of epoch
//     let mut dets: Vec<&Detection> = traj
//         .detection_ids
//         .iter()
//         .filter_map(|id| det_map.get(&(*id as u64)))
//         .collect();
//     dets.sort_by(|a, b| {
//         a.epoch
//             .epoch
//             .partial_cmp(&b.epoch.epoch)
//             .unwrap_or(Ordering::Equal)
//     });

//     let mut phi: Vec<f64> = Vec::with_capacity(dets.len());
//     let mut theta: Vec<f64> = Vec::with_capacity(dets.len());
//     let mut epochs: Vec<f64> = Vec::with_capacity(dets.len());

//     for det in dets.iter() {
//         if let Some((pt, _)) = construct_orbit(det, ic, Some(false)) {
//             let (ph, th) = phi_theta_from_point(&pt);
//             phi.push(ph);
//             theta.push(th);
//             epochs.push(det.epoch.epoch);
//         }
//     }

//     if phi.len() < 3 {
//         return None;
//     }

//     // Spread: sum of variances
//     let spread = variance(&phi) + variance(&theta);

//     // Linearity: mean R² of φ and θ vs. time
//     let r2_phi = r2_linear(&phi, &epochs);
//     let r2_theta = r2_linear(&theta, &epochs);
//     let linearity = 0.5 * (r2_phi + r2_theta);

//     Some((spread, linearity))
// }

// /// Unique key per trajectory to group candidates belonging to the same underlying track.
// #[derive(Clone, Debug, Eq)]
// struct TrajKey {
//     center_id: usize,
//     ids: Vec<usize>,
// }
// impl PartialEq for TrajKey {
//     fn eq(&self, other: &Self) -> bool {
//         self.center_id == other.center_id && self.ids == other.ids
//     }
// }
// impl Hash for TrajKey {
//     fn hash<H: Hasher>(&self, state: &mut H) {
//         self.center_id.hash(state);
//         self.ids.hash(state);
//     }
// }

// /// Normalise an array of values in place to [0,1] (min–max scaling).
// fn norm01(values: &mut [f64]) {
//     if values.is_empty() {
//         return;
//     }
//     let mut lo = f64::INFINITY;
//     let mut hi = f64::NEG_INFINITY;
//     for &x in values.iter() {
//         if x.is_finite() {
//             if x < lo {
//                 lo = x;
//             }
//             if x > hi {
//                 hi = x;
//             }
//         }
//     }
//     if !lo.is_finite() || !hi.is_finite() || hi == lo {
//         for v in values.iter_mut() {
//             *v = 0.0;
//         }
//         return;
//     }
//     let range = hi - lo;
//     for v in values.iter_mut() {
//         *v = (*v - lo) / range;
//     }
// }

// /// Compute spread & linearity for each candidate, normalise them across all candidates,
// /// combine into a single score, and keep the best candidate per trajectory.
// pub fn pick_best_by_metric(
//     candidates: &[(Trajectory, InitialCondition)],
//     det_map: &HashMap<u64, Detection>,
//     w_spread: f64,
//     w_linearity: f64,
// ) -> Vec<(Trajectory, InitialCondition)> {
//     // Compute raw metrics (spread, linearity) for each candidate
//     let raw: Vec<Option<(f64, f64)>> = candidates
//         .par_iter()
//         .map(|(traj, ic)| compute_raw_metrics(traj, ic, det_map))
//         .collect();

//     // Extract spreads and linearities to normalise
//     let mut spreads: Vec<f64> = Vec::with_capacity(candidates.len());
//     let mut linearities: Vec<f64> = Vec::with_capacity(candidates.len());
//     for r in &raw {
//         if let Some((sp, lin)) = r {
//             spreads.push(*sp);
//             linearities.push(*lin);
//         } else {
//             spreads.push(f64::INFINITY); // penalise invalid
//             linearities.push(0.0);
//         }
//     }

//     // Normalise
//     let mut spreads_n = spreads.clone();
//     norm01(&mut spreads_n);
//     let mut linearities_n = linearities.clone();
//     norm01(&mut linearities_n);

//     // Build a mapping: TrajKey -> (index, score)
//     let mut best_map: HashMap<TrajKey, (usize, f64)> = HashMap::new();

//     for (idx, (traj, _ic)) in candidates.iter().enumerate() {
//         // Build the key (centre_id + sorted detection_ids)
//         let mut ids = traj.detection_ids.clone();
//         ids.sort_unstable();
//         let key = TrajKey {
//             center_id: traj.center_id as usize,
//             ids,
//         };

//         // Combined score: lower is better
//         let score = w_spread * spreads_n[idx] + w_linearity * (1.0 - linearities_n[idx]);

//         match best_map.entry(key) {
//             Entry::Vacant(v) => {
//                 v.insert((idx, score));
//             }
//             Entry::Occupied(mut o) => {
//                 let (best_idx, best_score) = *o.get();
//                 if score < best_score {
//                     o.insert((idx, score));
//                 } else if (score - best_score).abs() < 1e-12 {
//                     // tie: prefer smaller p1 range, then lex order of params
//                     let ic_new = &candidates[idx].1;
//                     let ic_old = &candidates[best_idx].1;
//                     if better_ic(ic_new, ic_old) {
//                         o.insert((idx, score));
//                     }
//                 }
//             }
//         }
//     }

//     // Collect the winning indices
//     let mut result: Vec<(Trajectory, InitialCondition)> = Vec::with_capacity(best_map.len());
//     for (_, (idx, _)) in best_map {
//         result.push(candidates[idx].clone());
//     }
//     result
// }

// fn pick_top2_by_metric(
//     candidates: &[(Trajectory, InitialCondition)],
//     det_map: &HashMap<u64, Detection>,
//     w_spread: f64,
//     w_linearity: f64,
// ) -> Vec<(Trajectory, InitialCondition)> {
//     // 0) Compute normalised spread and linearity (as in pick_best_by_metric)
//     let raw: Vec<Option<(f64, f64)>> = candidates
//         .par_iter()
//         .map(|(traj, ic)| compute_raw_metrics(traj, ic, det_map))
//         .collect();

//     let mut spreads: Vec<f64> = Vec::with_capacity(candidates.len());
//     let mut linearities: Vec<f64> = Vec::with_capacity(candidates.len());
//     for r in &raw {
//         if let Some((sp, lin)) = r {
//             spreads.push(*sp);
//             linearities.push(*lin);
//         } else {
//             spreads.push(f64::INFINITY);
//             linearities.push(0.0);
//         }
//     }

//     let mut spreads_n = spreads.clone();
//     norm01(&mut spreads_n);
//     let mut linearities_n = linearities.clone();
//     norm01(&mut linearities_n);

//     // 1) Build a map from trajectory key to a vector of (index, score)
//     let mut top_map: HashMap<TrajKey, Vec<(usize, f64)>> = HashMap::new();
//     for (idx, (traj, _)) in candidates.iter().enumerate() {
//         let mut ids = traj.detection_ids.clone();
//         ids.sort_unstable();
//         let key = TrajKey { center_id: traj.center_id as usize, ids };
//         let score = w_spread * spreads_n[idx] + w_linearity * (1.0 - linearities_n[idx]);

//         let entry = top_map.entry(key).or_default();
//         entry.push((idx, score));

//         // Keep the vector sorted and truncated to 2
//         entry.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
//         if entry.len() > 5 {
//             entry.truncate(5);
//         }
//     }

//     // 2) Flatten the results into (traj, ic) pairs
//     let mut results: Vec<(Trajectory, InitialCondition)> = Vec::new();
//     for entry in top_map.values() {
//         for &(idx, _) in entry {
//             results.push(candidates[idx].clone());
//         }
//     }
//     results
// }

// /// Tie-breaker: smaller p1 width preferred; if equal, lexicographic on (p1, p2, p3, p4).
// fn better_ic(a: &InitialCondition, b: &InitialCondition) -> bool {
//     let wa = p1_width(a);
//     let wb = p1_width(b);
//     match (wa, wb) {
//         (Some(wa), Some(wb)) => {
//             if wa < wb - 1e-12 {
//                 return true;
//             }
//             if wa > wb + 1e-12 {
//                 return false;
//             }
//         }
//         (Some(_), None) => return true,
//         (None, Some(_)) => return false,
//         _ => {}
//     }

//     lex_params(a) < lex_params(b)
// }

// #[inline]
// fn p1_width(ic: &InitialCondition) -> Option<f64> {
//     match (ic.p1_min, ic.p1_max) {
//         (Some(a), Some(b)) => Some((b - a).abs()),
//         _ => None,
//     }
// }

// #[inline]
// fn lex_params(ic: &InitialCondition) -> (
//     OrderedFloat<f64>,
//     OrderedFloat<f64>,
//     OrderedFloat<f64>,
//     OrderedFloat<f64>,
// ) {
//     (
//         OrderedFloat(ic.p1),
//         OrderedFloat(ic.p2),
//         OrderedFloat(ic.p3),
//         OrderedFloat(ic.p4),
//     )
// }

// /// Collect all valid (traj, ic) pairs, score them, and return the best per trajectory.
// /// `w_spread` and `w_linearity` control the relative importance of tightness vs. linearity.
// use rayon::prelude::*;

// // ...

// pub fn refine_with_grid_cache_sparse_metric(
//     trajs: &[(Trajectory, InitialCondition)],
//     grid_cache: &Arc<Mutex<HashMap<GridCacheKey, Vec<Cell>>>>,
//     det_map: &HashMap<u64, Detection>,
//     epsilon: f64,
//     t_bounds: (f64, f64),
//     w_spread: f64,
//     w_linearity: f64,
// ) -> Vec<(Trajectory, InitialCondition)> {
//     // 1) Collect all valid refined pairs in parallel
//     let candidates: Vec<(Trajectory, InitialCondition)> = trajs
//         .par_iter()
//         .filter_map(|(traj, ic)| {
//             refine_single_trajectory_all(traj, ic, grid_cache, det_map, epsilon, t_bounds)
//         })
//         .flat_map(|vec| vec.into_par_iter())  // note: into_par_iter() produces a parallel iterator
//         .collect();

//     // 2) Score them and keep the best per trajectory
//     // pick_best_by_metric(&candidates, det_map, w_spread, w_linearity)
//     pick_top2_by_metric(&candidates, det_map, w_spread, w_linearity)
// }



// /// Version of refine_single_trajectory that returns **all** valid (traj, ic) pairs.
// /// all detections within `epsilon` of the centre.
// pub fn refine_single_trajectory_all(
//     traj: &Trajectory,
//     original_ic: &InitialCondition,
//     grid_cache: &Arc<Mutex<HashMap<GridCacheKey, Vec<Cell>>>>,
//     det_map: &HashMap<u64, Detection>,
//     epsilon: f64,
//     t_bounds: (f64, f64),
// ) -> Option<Vec<(Trajectory, InitialCondition)>> {
//     let t = (t_bounds.0 + t_bounds.1) / 2.0;
//     let ic_key = original_ic.to_bounds_key()?;
//     let grid_key = GridCacheKey {
//         ic_key,
//         epsilon: OrderedFloat(epsilon),
//     };

//     // Get or compute cells
//     let cells = {
//         let maybe_cells = {
//             let cache = grid_cache.lock().unwrap();
//             cache.get(&grid_key).cloned()
//         };
//         if let Some(c) = maybe_cells {
//             c
//         } else {
//             let mut computed = vec![Cell::from_initial_condition(original_ic)?];
//             let mut prev_len = 0;
//             while prev_len != computed.len() {
//                 prev_len = computed.len();
//                 computed = process_cells_parallel(
//                     computed,
//                     t,
//                     t_bounds,
//                     original_ic.mu,
//                     epsilon,
//                     original_ic.format.clone(),
//                 );
//             }
//             {
//                 let mut cache = grid_cache.lock().unwrap();
//                 cache.insert(grid_key.clone(), computed.clone());
//             }
//             computed
//         }
//     };

//     let center_det = det_map.get(&(traj.center_id as u64))?;
//     let ids = &traj.detection_ids;

//     let mut valid_results = Vec::new();

//     for (i, cell) in cells.iter().enumerate() {
//         let mid = cell.midpoint();
//         let refined_ic = InitialCondition::from_params(
//             i,
//             original_ic.format,
//             Some(cell.bounds[0].0),
//             Some(cell.bounds[0].1),
//             mid[0],
//             Some(cell.bounds[1].0),
//             Some(cell.bounds[1].1),
//             mid[1],
//             Some(cell.bounds[2].0),
//             Some(cell.bounds[2].1),
//             mid[2],
//             Some(cell.bounds[3].0),
//             Some(cell.bounds[3].1),
//             mid[3],
//             original_ic.epoch.clone(),
//             original_ic.mu,
//         );

//         // Propagate the centre
//         let (centre_pt, _) = match construct_orbit(center_det, &refined_ic, Some(false)) {
//             Some(v) => v,
//             None => continue,
//         };

//         // Check each detection
//         let mut all_within = true;
//         for &id in ids {
//             let det = match det_map.get(&(id as u64)) {
//                 Some(d) => d,
//                 None => {
//                     all_within = false;
//                     break;
//                 }
//             };
//             let (pt, _) = match construct_orbit(det, &refined_ic, Some(false)) {
//                 Some(v) => v,
//                 None => {
//                     all_within = false;
//                     break;
//                 }
//             };
//             let [dphi, dtheta] = angular_separation(centre_pt, pt);
//             let sep = (dphi * dphi + dtheta * dtheta).sqrt();
//             if sep > epsilon {
//                 all_within = false;
//                 break;
//             }
//         }

//         if all_within {
//             valid_results.push((traj.clone(), refined_ic));
//         }
//     }

//     if valid_results.is_empty() {
//         None
//     } else {
//         Some(valid_results)
//     }
// }
