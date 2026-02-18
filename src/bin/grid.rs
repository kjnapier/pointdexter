
// use spacerocks::transforms::{solve_for_universal_anomaly, stumpff_c, stumpff_s, calc_conic_anomaly_from_true_anomaly, calc_mean_anomaly_from_conic_anomaly};
// use spacerocks::coordinates::Origin;
// use std::collections::VecDeque;
// use std::f64::consts::PI;
// use argmin::prelude::*;
// use argmin::solver::goldensectionsearch::GoldenSectionSearch;
// use argmin::core::ArgminOp;
// use argmin::core::Error;
// use argmin::core::Executor;

// use rayon::prelude::*;
// use indicatif::{ProgressBar, ProgressStyle};
// use std::fs::OpenOptions;
// use std::io::Write;


// // instant
// use std::time::Instant;

// const ARCSEC_PER_RAD: f64 = 3600.0 * 180.0 / PI;

// #[derive(Clone, Debug, PartialEq)]
// struct Cell {
//     bounds: Vec<(f64, f64)>,
// }


// #[derive(Debug, Clone)]
// struct InitialCondition {
//     r: f64,
//     vr: f64,
//     vo: f64,
//     psi: f64,
//     epoch: f64,
//     mu: f64,
//     energy: f64,
// }

// impl InitialCondition {
//     fn new(r: f64, vr: f64, vo: f64, psi: f64, epoch: f64, mu: f64) -> Self {
//         let vsq = vr * vr + vo * vo;
//         let energy = 0.5 * vsq - mu / r;

//         Self {
//             r,
//             vr,
//             vo,
//             psi,
//             epoch,
//             mu,
//             energy,
//         }
//     }

//     fn from_scaled(gamma: f64, gamma_dot: f64, eta_x: f64, eta_y: f64, epoch: f64, mu: f64) -> Self {
//         let r = 1.0 / gamma;
//         let vr = gamma_dot * r;
//         let vo = (eta_x * eta_x + eta_y * eta_y).sqrt() * r;
//         let psi = eta_y.atan2(eta_x);

//         Self::new(r, vr, vo, psi, epoch, mu)
//     }

//     fn f_and_g_and_r(&self, epoch: f64) -> (f64, f64, f64) {
//         let dt = epoch - self.epoch;
//         let alpha = -2.0 * self.energy / self.mu;
//         let s = solve_for_universal_anomaly(self.r, self.vr, alpha, self.mu, dt, 1e-10, 100).unwrap();
//         let f = 1.0 - (s * s / self.r) * stumpff_c(alpha * s * s);
//         let g = dt - (s * s * s / self.mu.sqrt()) * stumpff_s(alpha * s * s);
//         let rsq = f * f * self.r * self.r + 2.0 * f * g * self.r * self.vr + g * g * (self.vo * self.vo + self.vr * self.vr);
//         (f, g, rsq.sqrt())
//     }
// }

// fn nice_acos(x: f64) -> f64 {
//     x.clamp(-1.0, 1.0).acos()
// }

// fn separation(epoch: f64, ic1: &InitialCondition, ic2: &InitialCondition) -> f64 {
//     let (f1, g1, r1) = ic1.f_and_g_and_r(epoch);
//     let (f2, g2, r2) = ic2.f_and_g_and_r(epoch);
//     let term1 = ic1.r * ic2.r * f1 * f2;
//     let term2 = ic1.r * ic2.vr * f1 * g2;
//     let term3 = ic1.vr * ic2.r * g1 * f2;
//     let term4 = g1 * g2 * (ic1.vr * ic2.vr + ic1.vo * ic2.vo * (ic1.psi - ic2.psi).cos());
//     let argument = (term1 + term2 + term3 + term4) / (r1 * r2);
//     nice_acos(argument) * ARCSEC_PER_RAD
// }


// fn maximize_separation(ic1: &InitialCondition, ic2: &InitialCondition, t_bounds: (f64, f64)) -> f64 {
//     // let a = separation(t_bounds.1, ic1, ic2);
//     // let b = separation(t_bounds.0, ic1, ic2);
//     // let mut max_sep = a.max(b);
//     // max_sep
//     separation(t_bounds.1, ic1, ic2)
// }


// fn check_axis(cell: &Cell, axis: usize, t_bounds: (f64, f64), mu: f64, epsilon: f64) -> (bool, f64) {
//     let combos = [0, 1];
//     let mut sep = 0.0;

//     for &i in &combos {
//         for &j in &combos {
//             for &k in &combos {
//                 let mut vertex = [0.0; 4]; //vec![0.0; 4];
//                 let mut idx = 0;
//                 for ax in 0..4 {
//                     if ax == axis {
//                         continue;
//                     }
//                     vertex[ax] = if [i, j, k][idx] == 0 {
//                         cell.bounds[ax].0
//                     } else {
//                         cell.bounds[ax].1
//                     };
//                     idx += 1;
//                 }
//                 let mut orbits = vec![];
//                 for &v in &[cell.bounds[axis].0, cell.bounds[axis].1] {
//                     vertex[axis] = v;
//                     orbits.push(InitialCondition::from_scaled(vertex[0], vertex[1], vertex[2], vertex[3], (t_bounds.0 + t_bounds.1) / 2.0, mu));
//                 }


//                 let max_sep = maximize_separation(&orbits[0], &orbits[1], t_bounds);
//                 if max_sep > epsilon {
//                     return (true, max_sep);
//                 }
//                 if max_sep > sep {
//                     sep = max_sep;
//                 }
//             }
//         }
//     }
//     (false, sep)
// }

// fn split_cell(cell: &Cell, t_bounds: (f64, f64), mu: f64, epsilon: f64) -> Option<Vec<Cell>> {

//     let mut seps = vec![];
//     for axis in 0..4 {
//         let (_, sep) = check_axis(cell, axis, t_bounds, mu, epsilon);
//         seps.push(sep);
//     }

//     let max_sep = seps
//         .iter()
//         .cloned()
//         .max_by(|a, b| a.partial_cmp(b).unwrap())
//         .unwrap();
//     if max_sep < epsilon {
//         return None;
//     }
//     let axis = seps.iter().position(|&x| x == max_sep).unwrap();
//     // println!("Splitting axis {} with separation {:.6}", axis, max_sep);
//     let mut bounds1 = cell.bounds.clone();
//     let mut bounds2 = cell.bounds.clone();
//     let mid = (cell.bounds[axis].0 + cell.bounds[axis].1) / 2.0;
//     bounds1[axis] = (cell.bounds[axis].0, mid);
//     bounds2[axis] = (mid, cell.bounds[axis].1);
//     Some(vec![Cell { bounds: bounds1 }, Cell { bounds: bounds2 }])
// }

// fn log_cell_midpoints(cells: &VecDeque<Cell>, filepath: &str) {


//     let mut file = OpenOptions::new()
//         .create(true)
//         .write(true)
//         .truncate(true)
//         .open(filepath)
//         .expect("Failed to open log file");

//     let line = "gamma,gamma_dot,eta_x,eta_y\n";
//     file.write_all(line.as_bytes()).expect("Write failed");


//     for cell in cells {
//         let midpoint: Vec<String> = cell
//             .bounds
//             .iter()
//             .map(|(lo, hi)| format!("{:.6}", (lo + hi) / 2.0))
//             .collect();
//         let line = format!("{}\n", midpoint.join(","));
//         file.write_all(line.as_bytes()).expect("Write failed");
//     }
// }

// fn main() {
//     let start = Instant::now();

//     let ssb = Origin::ssb();
//     let mu = ssb.mu();
//     let name = "cell_midpoints.csv";
//     let rmin = 60.0;
//     let rmax = 1000.0;
//     let vmax = (2.0 * mu / rmin).sqrt();

//     let bounds = vec![
//         (1.0 / rmax, 1.0 / rmin),
//         (-vmax / rmin, vmax / rmin),
//         (-vmax / rmin, vmax / rmin),
//         (-vmax / rmin, vmax / rmin),
//     ];
   

//     let t_bounds = (0.0, 13.0);
//     let epsilon = 0.15;


//     let mut to_check = VecDeque::new();
//     let mut completed = VecDeque::new();
//     to_check.push_back(Cell { bounds });


//     loop {
//         let old_len = to_check.len();
//         if old_len == 0 {
//             break;
//         }

//         let pb = ProgressBar::new(old_len as u64);
//         pb.set_style(ProgressStyle::default_bar()
//             .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
//             .unwrap()
//             .progress_chars("█▉▊▋▌▍▎▏  "));

//         let new_cells: Vec<_> = to_check
//             .par_iter()
//             .map(|cell| {
//                 let result = split_cell(cell, t_bounds, mu, epsilon);
//                 pb.inc(1);
//                 (cell.clone(), result)
//             })
//             .collect();

//         pb.finish_with_message("Round complete");

//         to_check.clear();
//         for (original_cell, maybe_split) in new_cells {
//             match maybe_split {
//                 Some(splits) => to_check.extend(splits),
//                 None => completed.push_back(original_cell),
//             }
//         }
//     }

//     let r_tol = 2.0;
//     // start splitting cells in gamma until rmax/rmin < r_tol
//     let mut completed = completed; //final_list;
//     loop {
//         let old_len = completed.len();
//         let mut to_check = VecDeque::new();
//         for cell in &completed {
//             let r_min = 1.0 / cell.bounds[0].1;
//             let r_max = 1.0 / cell.bounds[0].0;
//             if r_max / r_min > r_tol {
//                 // make 2 new cells to push to final_to_check
//                 let mut bounds1 = cell.bounds.clone();
//                 let mut bounds2 = cell.bounds.clone();
//                 let mid = (cell.bounds[0].0 + cell.bounds[0].1) / 2.0;
//                 bounds1[0] = (cell.bounds[0].0, mid);
//                 bounds2[0] = (mid, cell.bounds[0].1);
//                 to_check.push_back(Cell { bounds: bounds1 });
//                 to_check.push_back(Cell { bounds: bounds2 });
//                 // print the new r_min and r_max
//             } else {
//                 to_check.push_back(cell.clone());
//             }
//         }
//         if to_check.len() == old_len {
//             break;
//         }
//         completed = to_check;
//     }

//     // keep only cells where at least one point is bound
//     let mut final_cells: Vec<Cell> = vec![];
//     for cell in &completed {
//         let combos = [0, 1];
//         'outer: 
//         for &i in &combos {
//             for &j in &combos {
//                 for &k in &combos {
//                     for &l in &combos {
//                         let vertex = [
//                             if i == 0 { cell.bounds[0].0 } else { cell.bounds[0].1 },
//                             if j == 0 { 0.0 } else { 0.0 },
//                             if k == 0 { cell.bounds[2].0 } else { cell.bounds[2].1 },
//                             if l == 0 { cell.bounds[3].0 } else { cell.bounds[3].1 },
//                         ];
//                         let ic = InitialCondition::from_scaled(vertex[0], vertex[1], vertex[2], vertex[3], (t_bounds.0 + t_bounds.1) / 2.0, mu);
//                         // println!("Vertex: {:?}, Energy: {:.6}", vertex, ic.energy);
//                         if ic.energy < 0.0 {
//                             final_cells.push(cell.clone());
//                             break 'outer;
//                         }
//                     }
//                 }
//             }
//         }
//     } 


//     completed = VecDeque::from(final_cells);
//     println!("{} grid points", completed.len());



//     log_cell_midpoints(&completed, name);

//     println!("\nFinal number of completed cells: {}", completed.len());
//     let duration = start.elapsed();
//     println!("Time taken: {:?}", duration);

//     println!("{} grid points", completed.len());


// }


fn main() {
    println!("This is the grid.rs binary.");
}