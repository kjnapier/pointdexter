use serde_yaml;
use clap::Parser;
use itertools::Itertools;
use rayon::prelude::*;

use spacerocks::SpiceKernel;
use spacerocks::Time;

use pointdexter::cli::{Config, Cli};
use pointdexter::detection::Detection;
use pointdexter::io::{load_detections, load_initial_conditions};
use pointdexter::sync::*;
use pointdexter::InitialCondition;
use pointdexter::gauss::gauss;
use pointdexter::unit_vector_to_lonlat;
use pointdexter::hpix::vec3_to_nested_ipix;

use std::io::Write;
use std::collections::HashMap;
use std::f64::consts::PI;

use indicatif::ParallelProgressIterator;
use indicatif::ProgressIterator;

use cdshealpix::nested::{get, Layer};

pub fn sync_detections_to_orbit(detections: &Vec<Detection>, ic: &InitialCondition) -> Vec<Option<[f64; 3]>> {
    let mut synced_points: Vec<Option<[f64; 3]>> = vec![None; detections.len()];
    for (j, det) in detections.iter().enumerate() {
        if let Some(synced_pos) = sync_detection_to_orbit(det, ic) {
            synced_points[j] = Some(synced_pos);       
        }
    }
    synced_points
}

pub fn healpix_neighbors(depth: u8, hash: u64) -> Vec<u64> {
    let layer = get(depth);

    let neigh = layer.neighbours(hash, false);

    neigh.values_vec() 
        .into_iter()             // <-- this is the key
        .filter_map(|v| Some(v))     // remove None
        .collect()
}

pub fn count_1(indices: &Vec<u64>) -> Vec<(u64, usize)> {
    if indices.is_empty() {
        return Vec::new();
    }

    let mut sorted = indices.clone();
    sorted.sort_unstable();

    // worst case: all unique
    let mut counts = Vec::with_capacity(sorted.len());

    let mut current = sorted[0];
    let mut count = 1usize;

    for &idx in &sorted[1..] {
        if idx == current {
            count += 1;
        } else {
            counts.push((current, count));
            current = idx;
            count = 1;
        }
    }

    counts.push((current, count));

    // now sort the counts by count
    counts.sort_unstable_by(|a, b| b.1.cmp(&a.1));
    counts
}


pub fn count_2(indices: &Vec<u64>) -> Vec<(u64, usize)> {
    let mut counts: HashMap<u64, usize> = HashMap::new();
    for idx in indices {
        *counts.entry(*idx).or_insert(0) += 1;
    }
    counts.into_iter().collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {

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

    const ARCSEC_PER_RAD: f64 = (180.0 * 3600.0) / std::f64::consts::PI;
    let eps: f64 = config.epsilon_arcsec / ARCSEC_PER_RAD;
    let cluster_radius_squared = 2.0 * (1.0 - eps.cos());

    let mut cluster_file = std::fs::File::create("clusters.txt")?;
    let cluster_file = std::sync::Mutex::new(&mut cluster_file);


    let depth = 15_u8;


    // Sync detections to orbits
    let mut idx = 0;
    for ic in ics[..500].iter() {
    // ics.par_iter().progress_count(ics.len() as u64).for_each(|ic| {
        let start = std::time::Instant::now();

        let mut canvas: HashMap<u64, usize> = HashMap::new();

        let synced_points = sync_detections_to_orbit(&detections, ic);
        
        let mut hpix_indices: Vec<u64> = Vec::with_capacity(synced_points.len() * 8);
        for (i, p) in synced_points.iter().enumerate() {
            if let Some(v) = p {
                let idx = vec3_to_nested_ipix(depth, v[0], v[1], v[2]);
                hpix_indices.push(idx);
            }
        }

        let counts = count_1(&hpix_indices);

        let duration = start.elapsed();
        println!("Synced points for IC {}: {}/{} in {:.2?}", ic.id, synced_points.iter().filter(|p| p.is_some()).count(), synced_points.len(), duration);

        }   
        
//    });
    // }
    Ok(())
}

