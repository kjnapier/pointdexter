use rand::seq::IndexedRandom;
use rand::rng;
use crate::sparse_linking_stuff::InitialCondition;

pub fn get_n_random_orbits<'a>(ics: &'a [InitialCondition], n: usize) -> Vec<&'a InitialCondition> {
    // Check if we have enough points to sample
    if ics.len() < n {
        return ics.iter().collect();
    }

    let mut rng = rng();
    // Using choose_multiple with copied references
    ics.choose_multiple(&mut rng, n).collect()
}

pub const ARCSEC_PER_RAD: f64 = 3600.0 * 180.0 / std::f64::consts::PI;