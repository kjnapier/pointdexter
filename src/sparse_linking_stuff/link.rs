use crate::sparse_linking_stuff::orbfit::fitter::FitResult;
//use spacerocks::orbfit::fitter::FitResult;

use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Link {
    pub matching_cluster_id: usize,
    pub objids: Vec<String>,
    pub best_fit: FitResult,
    pub q: f64,
    pub e: f64,
    pub psi: f64,
    pub true_anomaly: f64,
}

impl Link {
    pub fn new(matching_cluster_id: usize, objids: Vec<String>, best_fit: FitResult, q:f64, e:f64, psi:f64, true_anomaly:f64) -> Link {
        Link {
            matching_cluster_id,
            objids,
            best_fit,
            q,
            e,
            psi,
            true_anomaly,
        }
    }
}

// implement Hash for Link
impl std::hash::Hash for Link {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let mut objids = self.objids.clone();
        objids.sort();
        objids.hash(state)
    }
}

// implement ParialEq for Link
impl std::cmp::PartialEq for Link {
    fn eq(&self, other: &Self) -> bool {
        let mut objids = self.objids.clone();
        objids.sort();
        let mut other_objids = other.objids.clone();
        other_objids.sort();
        objids == other_objids
    }
}

impl std::cmp::Eq for Link {}