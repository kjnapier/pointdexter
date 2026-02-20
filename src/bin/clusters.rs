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
use std::collections::HashMap;
use std::f64::consts::PI;

use indicatif::ParallelProgressIterator;
use indicatif::ProgressIterator;

pub fn sync_detections_to_orbit(detections: &Vec<Detection>, ic: &InitialCondition) -> Vec<Option<[f64; 3]>> {
    let mut synced_points: Vec<Option<[f64; 3]>> = vec![None; detections.len()];
    for (j, det) in detections.iter().enumerate() {
        if let Some(synced_pos) = sync_detection_to_orbit(det, ic) {
            synced_points[j] = Some(synced_pos);       
        }
    }
    synced_points
}


pub fn build_grid(points: &Vec<Option<[f64; 3]>>, pixel_scale: f64) -> Vec<Vec<usize>> {
    // now accumulate the synced points into a 2d histogram
    let mut phi_min = std::f64::INFINITY;
    let mut phi_max = std::f64::NEG_INFINITY;
    let mut theta_min = std::f64::INFINITY;;
    let mut theta_max = std::f64::NEG_INFINITY;
    let mut thetas: Vec<f64> = Vec::with_capacity(points.len());
    let mut phis: Vec<f64> = Vec::with_capacity(points.len());
    for (idx, point) in points.iter().enumerate() {
        if point.is_none() {
            continue;
        }
        let point = point.unwrap();
        let x = point[0];
        let y = point[1];
        let z = point[2];
        let theta = (z).asin();
        let phi = y.atan2(x);

        if phi < phi_min {
            phi_min = phi;
        }
        if phi > phi_max {
            phi_max = phi;
        }
        if theta < theta_min {
            theta_min = theta;
        }
        if theta > theta_max {
            theta_max = theta;
        }
        thetas.push(theta);
        phis.push(phi);
    }

    // make a histogram grid with 5 arcsec bins
    let bin_size = pixel_scale;
    // make time bins in a third dimension.
    let n_phi_bins = ((phi_max - phi_min) / bin_size).ceil() as usize;
    let n_theta_bins = ((theta_max - theta_min) / bin_size).ceil() as usize;
    // 2d histogram of counts in each bin phi
    let mut histogram = vec![vec![0usize; n_phi_bins]; n_theta_bins];

    let mut peaks: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    for idx in 0..phis.len() {
        let theta = thetas.get(idx);
        let phi = phis.get(idx);
        let theta = *theta.unwrap();
        let phi = *phi.unwrap();
        
        let phi_bin = ((phi - phi_min) / bin_size).floor() as usize;
        let theta_bin = ((theta - theta_min) / bin_size).floor() as usize;
        
        histogram[theta_bin][phi_bin] += 1;
        
    }
    histogram
}

pub fn find_local_maxima(grid: &Vec<Vec<usize>>) -> Vec<(usize, usize)> {
    let mut peaks: Vec<(usize, usize)> = Vec::new();
    for i in 0..grid.len() {
        for j in 0..grid[i].len() {
            let count = grid[i][j];
            if count == 0 {
                continue;
            }
            let mut is_peak = true;
            for di in -1..=1 {
                for dj in -1..=1 {
                    if di == 0 && dj == 0 {
                        continue;
                    }   
                    let ni = i as isize + di;
                    let nj = j as isize + dj;
                    if ni < 0 || ni >= grid.len() as isize || nj < 0 || nj >= grid[i].len() as isize {
                        continue;
                    }
                    if grid[ni as usize][nj as usize] > count {
                        is_peak = false;
                        break;
                    }
                }
                if !is_peak {
                    break;
                }
            }
            if is_peak {
                peaks.push((i, j));
            }
        }
    }
    peaks
}

pub fn convolve_peaks(grid: &Vec<Vec<usize>>, peaks: &Vec<(usize, usize)>) -> Vec<usize> {
    let mut counts: Vec<usize> = Vec::with_capacity(peaks.len());
    for peak in peaks.iter() {
        let i = peak.0;
        let j = peak.1;
        let mut count_sum = 0;
        for di in -1..=1 {
            for dj in -1..=1 {
                let ni = i as isize + di;
                let nj = j as isize + dj;
                if ni < 0 || ni >= grid.len() as isize || nj < 0 || nj >= grid[i].len() as isize {
                    continue;
                }
                count_sum += grid[ni as usize][nj as usize];
            }
        }
        counts.push(count_sum);
    } 
    counts
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


    // Sync detections to orbits
    let mut idx = 0;
    for ic in ics[..].iter() {
    // ics.par_iter().progress_count(ics.len() as u64).for_each(|ic| {
        let start = std::time::Instant::now();

        let synced_points = sync_detections_to_orbit(&detections, ic);
        let grid = build_grid(&synced_points, pixel_scale);

        let peaks = find_local_maxima(&grid);
        let counts = convolve_peaks(&grid, &peaks);
        
        let mut cluster_count = 0;
        for (i, count) in counts.iter().enumerate() {
            if *count >= 12 {
                cluster_count += 1;
                writeln!(cluster_file.lock().unwrap(), "{} {} {} {}", ic.id, peaks[i].0, peaks[i].1, count)?;
            }
        }

        println!("Found {} clusters for IC {}.", cluster_count, ic.id);


        let duration = start.elapsed();
        println!("Synced points for IC {}: {}/{} in {:.2?}", ic.id, synced_points.iter().filter(|p| p.is_some()).count(), synced_points.len(), duration);

        }   
        
//    });
    // }
    Ok(())
}

