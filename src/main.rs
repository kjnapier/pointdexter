use serde_yaml;
use clap::Parser;

use spacerocks::SpiceKernel;

use pointdexter::cli::{Config, Cli};
use pointdexter::detection::Detection;
use pointdexter::io::{load_detections, load_initial_conditions};
use pointdexter::sync::*;
use pointdexter::InitialCondition;

use kiddo::KdTree; 
use kiddo::SquaredEuclidean;

use std::collections::HashSet;
use std::error::Error;
use image::{GrayImage, Luma};

use std::f64::consts::PI;
pub fn save_grid_png(grid: &HitGrid, path: &str) -> Result<(), Box<dyn Error>> {
    let nx = grid.nx as u32;
    let ny = grid.ny as u32;

    // Find max for scaling
    // let max_count = grid.counts.iter().copied().max().unwrap_or(0.0);
    let max_count = grid.counts.iter().fold(0.0_f64, |a, &b| a.max(b));

    // Avoid divide-by-zero; if everything is zero, just write a black image
    let use_log = max_count > 0.0;
    let max_log = if use_log {
        (max_count as f32).ln_1p()
    } else {
        1.0
    };

    let mut img = GrayImage::new(nx, ny);

    for iy in 0..grid.ny {
        for ix in 0..grid.nx {
            let idx = iy * grid.nx + ix;
            let c = grid.counts[idx];

            // Log stretch for dynamic range
            let val: u8 = if use_log {
                let v = (c as f32).ln_1p() / max_log;
                (v.clamp(0.0, 1.0) * 255.0).round() as u8
            } else {
                0
            };

            // Flip vertically so (0,0) is *bottom-left* in the data
            let y_img = ny - 1 - (iy as u32);
            img.put_pixel(ix as u32, y_img, Luma([val]));
        }
    }

    img.save(path)?;
    Ok(())
}




#[derive(Debug)]
pub struct HitGrid {
    /// Row-major counts (iy * nx + ix)
    pub counts: Vec<f64>,
    pub nx: usize,
    pub ny: usize,

    /// Center RA of the grid (true RA, radians)
    pub ra_center: f64,

    /// Minimum RA offset from center (radians) at the left edge of the grid
    /// True RA at left edge = normalize_angle(ra_center + x_min)
    pub x_min: f64,

    /// Minimum Dec of the grid (radians) at the bottom edge of the grid
    pub dec_min: f64,

    /// Pixel scale (radians) in both RA-offset and Dec
    pub pix_scale_rad: f64,
}

impl HitGrid {
    #[inline]
    pub fn index(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    pub fn get(&self, ix: usize, iy: usize) -> Option<f64> {
        if ix < self.nx && iy < self.ny {
            Some(self.counts[self.index(ix, iy)])
        } else {
            None
        }
    }

    /// True RA of left edge (radians), for convenience
    #[inline]
    pub fn ra_min(&self) -> f64 {
        normalize_angle(self.ra_center + self.x_min)
    }

    /// True RA/Dec at the center of pixel (ix, iy), if needed
    pub fn pixel_radec(&self, ix: usize, iy: usize) -> Option<(f64, f64)> {
        if ix >= self.nx || iy >= self.ny {
            return None;
        }
        let x = self.x_min + (ix as f64 + 0.5) * self.pix_scale_rad;
        let ra = normalize_angle(self.ra_center + x);
        let dec = self.dec_min + (iy as f64 + 0.5) * self.pix_scale_rad;
        Some((ra, dec))
    }
}

/// Normalize angle to (-π, π]
#[inline]
fn normalize_angle(angle: f64) -> f64 {
    let two_pi = 2.0 * PI;
    let mut a = angle % two_pi;
    if a <= -PI {
        a += two_pi;
    } else if a > PI {
        a -= two_pi;
    }
    a
}

/// Convert a unit vector (x, y, z) to (ra, dec) in radians
#[inline]
fn vec_to_radec(v: [f64; 3]) -> (f64, f64) {
    let x = v[0];
    let y = v[1];
    let z = v[2];
    let ra = y.atan2(x); // [-π, π)
    let r = (x * x + y * y).sqrt();
    let dec = z.atan2(r); // [-π/2, π/2]
    (ra, dec)
}

/// Build a 2D histogram over unit vectors with given pixel size (arcsec),
/// defining RA in a local coordinate around the circular mean of the data.
pub fn build_hit_grid(points: &[[f64; 3]], pix_scale_arcsec: f64) -> HitGrid {
    assert!(pix_scale_arcsec > 0.0, "pixel scale must be > 0");
    assert!(!points.is_empty(), "no points provided");

    // 1 arcsec in radians
    let arcsec_to_rad = PI / (180.0 * 3600.0);
    let pix_scale_rad = pix_scale_arcsec * arcsec_to_rad;

    let n = points.len();
    let mut ras = Vec::with_capacity(n);
    let mut decs = Vec::with_capacity(n);

    // First pass: convert to (RA, Dec), accumulate for circular mean of RA
    let mut sum_sin = 0.0_f64;
    let mut sum_cos = 0.0_f64;

    for &p in points {
        let (ra, dec) = vec_to_radec(p);
        ras.push(ra);
        decs.push(dec);
        sum_sin += ra.sin();
        sum_cos += ra.cos();
    }

    // Circular mean RA = center of the data on the circle
    let ra_center = sum_sin.atan2(sum_cos);

    // Second pass: compute offsets from center, and min/max of offsets and Dec
    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut dec_min = f64::INFINITY;
    let mut dec_max = f64::NEG_INFINITY;

    // We can reuse ras/decs as we loop, no need to store a second Vec
    let mut offsets: Vec<(f64, f64)> = Vec::with_capacity(n);

    for i in 0..n {
        let ra = ras[i];
        let dec = decs[i];

        // Offset from center, wrapped to (-π, π]
        let x = normalize_angle(ra - ra_center);
        offsets.push((x, dec));

        if x < x_min { x_min = x; }
        if x > x_max { x_max = x; }
        if dec < dec_min { dec_min = dec; }
        if dec > dec_max { dec_max = dec; }
    }

    // Degenerate cases
    if x_max == x_min {
        x_max = x_min + pix_scale_rad;
    }
    if dec_max == dec_min {
        dec_max = dec_min + pix_scale_rad;
    }

    let nx = ((x_max - x_min) / pix_scale_rad).ceil() as usize;
    let ny = ((dec_max - dec_min) / pix_scale_rad).ceil() as usize;

    let mut counts = vec![0.0f64; nx * ny];

    // Third pass: bin each (x, dec)
    for (i, (x, dec)) in offsets.iter().enumerate() {
        let ix = ((x - x_min) / pix_scale_rad).floor() as isize;
        let iy = ((dec - dec_min) / pix_scale_rad).floor() as isize;

        if ix >= 0 && (ix as usize) < nx && iy >= 0 && (iy as usize) < ny {
            let idx = (iy as usize) * nx + (ix as usize);
            counts[idx] += 1.0;
        }
    }

    HitGrid {
        counts,
        nx,
        ny,
        ra_center,
        x_min,
        dec_min,
        pix_scale_rad,
    }
}


pub fn sync_detections_to_orbit(detections: &Vec<Detection>, ic: &InitialCondition) -> Vec<Option<[f64; 3]>> {
    let mut synced_points: Vec<Option<[f64; 3]>> = vec![None; detections.len()];
    for (j, det) in detections.iter().enumerate() {
        if let Some(synced_pos) = sync_detection_to_orbit(det, ic) {
            synced_points[j] = Some(synced_pos);       
        }
    }
    synced_points
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
    let cluster_radius = 2.0 * (1.0 - eps.cos());

    type P = [f64; 3];

    // Sync detections to orbits
    let mut idx = 0;
    for ic in ics[..100].iter() {

        let synced_points: Vec<Option<P>> = sync_detections_to_orbit(&detections, ic);

        // // Collect valid points + payloads contiguously
        // let mut pts = Vec::with_capacity(synced_points.len());
        // for (i, p) in synced_points.iter().enumerate() {
        //     if let Some(v) = p {
        //         pts.push((*v, i as u64));
        //     }
        // }

        let mut points = Vec::with_capacity(synced_points.len());
        for p in synced_points.iter() {
            if let Some(v) = p {
                points.push(*v);
            }
        }

        if points.is_empty() {
            println!("No synced points for IC ID: {}", ic.id);
            continue;
        }

        // // Bulk build (name may be `from_iter`, `from_vec`, `build`, etc.)
        // let tree = KdTree::from_iter(pts.into_iter());

        let grid = build_hit_grid(&points, 5.0);
        let png_path = format!("/Users/kjnapier/Desktop/ims/hitgrid_ic_{:03}.png", idx);
        save_grid_png(&grid, &png_path)?;
        idx += 1;



        

        println!("Initial Condition ID: {}", ic.id);
        // println!("Number of clusters: {}", clusters.len());

    }
    
    Ok(())

}


