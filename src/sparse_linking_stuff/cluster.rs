// use crate::detection::Detection;
// use crate::initial_condition::InitialCondition;
// use crate::construct_orbit::construct_orbit;
// use crate::metrics::{one_minus_dot_product};

// use nalgebra::{DVector, DMatrix};


// #[derive(Debug, Clone)]
// pub struct Cluster<'a> {
//     pub detections: Vec<&'a Detection>,
// }

// impl<'a> Cluster<'a> {

//     pub fn new(detections: Vec<&'a Detection>) -> Self {
//         Self {
//             detections,
//         }
//     }

//     pub fn len(&self) -> usize {
//         self.detections.len()
//     }

//     pub fn timespan(&self) -> f64 {
//         let mut min_time = self.detections[0].epoch.epoch;
//         let mut max_time = self.detections[0].epoch.epoch;
//         for detection in &self.detections {
//             if detection.epoch.epoch < min_time {
//                 min_time = detection.epoch.epoch;
//             }
//             if detection.epoch.epoch > max_time {
//                 max_time = detection.epoch.epoch;
//             }
//         }
//         max_time - min_time
//     }

// }

// pub fn median(x: &Vec<f64>) -> f64 {
//     let mut x = x.clone();
//     x.sort_by(|a, b| a.partial_cmp(b).unwrap());
//     let mid = x.len() / 2;
//     if x.len() % 2 == 0 {
//         (x[mid] + x[mid - 1]) / 2.0
//     } else {
//         x[mid]
//     }
// }

// impl<'a> std::hash::Hash for Cluster<'a> {
//     fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
//         let mut intids = self.detections.iter().map(|d| d.intid).collect::<Vec<_>>();
//         intids.sort();
//         intids.hash(state)
//     }
// }

// impl<'a> std::cmp::PartialEq for Cluster<'a> {
//     fn eq(&self, other: &Self) -> bool {
//         let mut intids = self.detections.iter().map(|d| d.intid).collect::<Vec<_>>();
//         intids.sort();
//         let mut other_intids = other.detections.iter().map(|d| d.intid).collect::<Vec<_>>();
//         other_intids.sort();
//         intids == other_intids
//     }
// }

// impl<'a> std::cmp::Eq for Cluster<'a> {}


