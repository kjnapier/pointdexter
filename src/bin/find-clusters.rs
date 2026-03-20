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

use std::io::Write;
// use std::collections::HashMap;
use std::f64::consts::PI;

use indicatif::ParallelProgressIterator;
use indicatif::ProgressIterator;

use hashbrown::HashMap;
use ahash::RandomState;

pub fn sync_detections_to_orbit(detections: &Vec<Detection>, ic: &InitialCondition) -> Vec<Option<[f64; 3]>> {
    let mut synced_points: Vec<Option<[f64; 3]>> = vec![None; detections.len()];
    for (j, det) in detections.iter().enumerate() {
        if let Some(synced_pos) = sync_detection_to_orbit(det, ic) {
            synced_points[j] = Some(synced_pos);       
        }
    }
    synced_points
}

pub type CountsMap = HashMap<u64, u32, RandomState>;
// Return sparse counts + WCS params + bin counts
pub fn build_grid_sparse(alphas: &[f64], betas: &[f64], pixel_scale: f64) -> (CountsMap, f64, f64, f64, usize, usize) {
    debug_assert_eq!(alphas.len(), betas.len());
    let n = alphas.len();
    if n == 0 {
        return (CountsMap::with_hasher(RandomState::new()), 0.0, 0.0, pixel_scale, 0, 0);
    }

    let mut alpha_min = alphas[0];
    let mut alpha_max = alphas[0];
    let mut beta_min  = betas[0];
    let mut beta_max  = betas[0];

    for i in 1..n {
        let a = alphas[i];
        let b = betas[i];
        if a < alpha_min { alpha_min = a; }
        if a > alpha_max { alpha_max = a; }
        if b < beta_min  { beta_min  = b; }
        if b > beta_max  { beta_max  = b; }
    }

    let inv = 1.0 / pixel_scale;
    // let alpha_bins = ((alpha_max - alpha_min) * inv).ceil() as usize;
    // let beta_bins  = ((beta_max  - beta_min)  * inv).ceil() as usize;

    let alpha_bins = (((alpha_max - alpha_min) * inv).ceil() as usize).max(1);
    let beta_bins  = (((beta_max  - beta_min)  * inv).ceil() as usize).max(1);

    
    // let mut counts: HashMap<u64, u32> = HashMap::with_capacity(n * 2);
    let mut counts: CountsMap = CountsMap::with_capacity_and_hasher(n * 2, RandomState::new());


    for (&a, &b) in alphas.iter().zip(betas.iter()) {
        let mut ai = ((a - alpha_min) * inv) as isize; // trunc toward 0; OK since >=0 typically
        let mut bi = ((b - beta_min)  * inv) as isize;

        // Clamp just in case of edge/rounding issues
        if ai < 0 { ai = 0; }
        if bi < 0 { bi = 0; }
        if ai as usize >= alpha_bins { ai = alpha_bins.saturating_sub(1) as isize; }
        if bi as usize >= beta_bins  { bi = beta_bins.saturating_sub(1)  as isize; }

        let key = pack(ai as u32, bi as u32);
        *counts.entry(key).or_insert(0) += 1;
    }

    (counts, alpha_min, beta_min, pixel_scale, alpha_bins, beta_bins)
}


#[inline]
fn pack(i: u32, j: u32) -> u64 {
    ((i as u64) << 32) | (j as u64)
}

#[inline]
fn unpack(key: u64) -> (u32, u32) {
    ((key >> 32) as u32, (key & 0xFFFF_FFFF) as u32)
}

/// Sum counts in the 3x3 neighborhood around each peak (including the center).
pub fn convolve_peaks_sparse(
    counts: &CountsMap,
    peaks: &[(u32, u32)],
    alpha_bins: u32,
    beta_bins: u32,
) -> Vec<u32> {
    let mut out = Vec::with_capacity(peaks.len());

    for &(i, j) in peaks {
        let mut sum: u32 = 0;

        // neighborhood bounds (avoid per-neighbor checks where possible)
        let i0 = i.saturating_sub(1);
        let j0 = j.saturating_sub(1);
        let i1 = (i + 1).min(alpha_bins.saturating_sub(1));
        let j1 = (j + 1).min(beta_bins.saturating_sub(1));

        for ni in i0..=i1 {
            for nj in j0..=j1 {
                if let Some(&v) = counts.get(&pack(ni, nj)) {
                    sum += v;
                }
            }
        }

        out.push(sum);
    }

    out
}

/// Find local maxima among occupied cells, skipping any cell with count <= min_count.
/// A cell is a peak if no neighbor in its 8-neighborhood has a strictly larger count.
/// (Ties are allowed.)
pub fn find_local_maxima_sparse(
    counts: &CountsMap,
    alpha_bins: u32,
    beta_bins: u32,
    min_count: u32, // pass 3 to "skip <= 2"
) -> Vec<(u32, u32)> {
    let mut peaks: Vec<(u32, u32)> = Vec::new();

    for (&key, &count) in counts.iter() {
        if count < min_count {
            continue;
        }

        let (i, j) = unpack(key);

        let i0 = i.saturating_sub(1);
        let j0 = j.saturating_sub(1);
        let i1 = (i + 1).min(alpha_bins.saturating_sub(1));
        let j1 = (j + 1).min(beta_bins.saturating_sub(1));

        let mut is_peak = true;

        'nbrs: for ni in i0..=i1 {
            for nj in j0..=j1 {
                if ni == i && nj == j {
                    continue;
                }
                if let Some(&v) = counts.get(&pack(ni, nj)) {
                    if v > count {
                        is_peak = false;
                        break 'nbrs;
                    }
                }
            }
        }

        if is_peak {
            peaks.push((i, j));
        }
    }

    peaks
}


const ARCSEC_PER_RAD: f64 = (180.0 * 3600.0) / std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let args = Cli::try_parse()?;
    let f = std::fs::File::open(args.config)?;
    let config: Config = serde_yaml::from_reader(f)?;

    let mut kernel = SpiceKernel::new();
    kernel.load_spk(format!("{}/sb441-n16.bsp", config.spice_path).as_str())?;
    kernel.load_spk(format!("{}/de440s.bsp", config.spice_path).as_str())?;
    kernel.load_bpc(format!("{}/earth_1962_240827_2124_combined.bpc", config.spice_path).as_str())?;

    let mut detections = load_detections(&config.detection_catalog, &config.orbit_reference_plane, &kernel)?;
    println!("Loaded {} detections.", detections.len());

    let mut ics = load_initial_conditions(&config.initial_conditions_file, &config.ic_type, &config.ic_origin, config.reference_epoch)?;
    println!("Loaded {} initial conditions.", ics.len());

    
    let eps: f64 = config.epsilon_arcsec / ARCSEC_PER_RAD;
    let pixel_scale = eps / 3.0;

    let mut cluster_file = std::fs::File::create("clusters.txt")?;
    let cluster_file = std::sync::Mutex::new(&mut cluster_file);

    // let alpha_min, beta_min, pixel_scale, alpha_bins, beta_bins = 20;
    let thresh = 10;

    // Sync detections to orbits
    let mut idx = 0;
    // for ic in ics[..].iter() {
    ics.par_iter().progress_count(ics.len() as u64).for_each(|ic| {

        let synced_points = sync_detections_to_orbit(&detections, ic);
        
        let mut phis: Vec<f64> = Vec::with_capacity(synced_points.len());
        let mut thetas: Vec<f64> = Vec::with_capacity(synced_points.len());
        for p_opt in synced_points.iter() {
            if let Some(pos) = p_opt {
                let r = (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]).sqrt();
                let phi = pos[1].atan2(pos[0]);
                let theta = (pos[2] / r).asin();
                phis.push(phi);
                thetas.push(theta);
            }
        }

        let (grid, alpha_min, beta_min, pixel_scale, alpha_bins, beta_bins) = build_grid_sparse(&phis, &thetas, pixel_scale);
        let peaks = find_local_maxima_sparse(&
            grid,
            alpha_bins as u32,
            beta_bins as u32,
            3, // min_count
        );

        let counts = convolve_peaks_sparse(&grid, &peaks, alpha_bins as u32, beta_bins as u32);

        let mut peak_info: Vec<(f64, f64, u32)> = Vec::new();
        for (i, count) in counts.iter().enumerate() {
            if *count >= thresh {
                let (ai, bi) = peaks[i];
                let ai = alpha_min + (ai as f64 + 0.5) * pixel_scale; // center of bin
                let bi = beta_min  + (bi as f64 + 0.5) * pixel_scale;
                peak_info.push((ai, bi, *count));
            }
        }

        // acquire lock and write cluster info to file
        {
            let mut file = cluster_file.lock().unwrap();
            // writeln!(file, "IC {}: Found {} clusters.", ic.id, cluster_count)?;
            for (ai, bi, count) in peak_info.iter() {
                writeln!(file, "{}, {:.10}, {:.10}, {}", ic.id, ai, bi, count).unwrap();
            }
        }


        // }   
        
   });
    // }
    Ok(())
}

