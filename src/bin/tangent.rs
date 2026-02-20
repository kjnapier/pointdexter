use pointdexter::cli::{Config, Cli};
use pointdexter::Exposure;
use pointdexter::load_detections;
use pointdexter::xyz_to_proj_matrix;

use spacerocks::{SpiceKernel};
use pointdexter::Detection;
use pointdexter::TangentPlaneExposure;
use nalgebra::Vector3;
use clap::Parser;

use spacerocks::coordinates::Origin;
use spacerocks::time::Time;
use spacerocks::constants::SPEED_OF_LIGHT;

pub fn build_grid(alphas: &Vec<f64>, betas: &Vec<f64>, pixel_scale: f64) -> Vec<Vec<usize>> {
    let mut alpha_min = std::f64::INFINITY;
    let mut alpha_max = std::f64::NEG_INFINITY;
    let mut beta_min = std::f64::INFINITY;;
    let mut beta_max = std::f64::NEG_INFINITY;
    for (alpha, beta) in alphas.iter().zip(betas.iter()) {
        if *alpha < alpha_min {
            alpha_min = *alpha;
        }
        if *alpha > alpha_max {
            alpha_max = *alpha;
        }
        if *beta < beta_min {
            beta_min = *beta;
        }
        if *beta > beta_max {
            beta_max = *beta;
        }
    }
    let alpha_bins = ((alpha_max - alpha_min)/pixel_scale).ceil() as usize;
    let beta_bins = ((beta_max - beta_min)/pixel_scale).ceil() as usize;
    let mut grid = vec![vec![0; beta_bins]; alpha_bins];
    for (alpha, beta) in alphas.iter().zip(betas.iter()) {
        let alpha_idx = ((*alpha - alpha_min)/pixel_scale).floor() as usize;
        let beta_idx = ((*beta - beta_min)/pixel_scale).floor() as usize;
        //println!("alpha: {}, beta: {}, alpha_idx: {}, beta_idx: {}", alpha, beta, alpha_idx, beta_idx);
        grid[alpha_idx][beta_idx] += 1;
    }
    //println!("Grid size: {} x {}", alpha_bins, beta_bins);
    grid
}

pub fn save_grid_png(grid: &Vec<Vec<usize>>, path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use image::{GrayImage, Luma};
    let alpha_bins = grid.len();
    let beta_bins = grid[0].len();
    let mut img = GrayImage::new(beta_bins as u32, alpha_bins as u32);
    for (i, row) in grid.iter().enumerate() {
        for (j, count) in row.iter().enumerate() {
            let intensity = (*count as f64 / 10.0).min(1.0) * 255.0; // scale counts to [0, 255], adjust as needed
            img.put_pixel(j as u32, i as u32, Luma([intensity as u8]));
        }
    }
    img.save(path)?;
    Ok(())
}



pub fn sort_detections_into_exposures(detections: &Vec<Detection>) -> Vec<Exposure> {
    // sort the detections by epoch
    let mut detections = detections.clone();
    detections.sort_by(|a, b| a.epoch.partial_cmp(&b.epoch).unwrap());
    let mut exposure_groups: Vec<Vec<Detection>> = Vec::new();
    let mut current_group: Vec<Detection> = Vec::new();
    let mut current_epoch = detections[0].epoch;
    for detection in detections.iter() {
        if (detection.epoch - current_epoch).abs() < 1e-6 {
            current_group.push(detection.clone());
        } else {
            exposure_groups.push(current_group);
            current_group = vec![detection.clone()];
            current_epoch = detection.epoch;
        }
    }
    if !current_group.is_empty() {
        exposure_groups.push(current_group);
    }

    let mut exposures: Vec<Exposure> = Vec::with_capacity(exposure_groups.len());
    for (idx, group) in exposure_groups.iter().enumerate() {
        let exposure = Exposure {
            id: format!("exposure_{}", idx),
            epoch: group[0].epoch,
            filter: None,
            detections: group.iter().map(|d| d.rho_hat).collect(),
            observer_position: group[0].observer_position,
            observer_velocity: group[0].observer_velocity,
            reference_plane: group[0].reference_plane.clone(),
        };
        exposures.push(exposure);
    }
    exposures
}

pub fn transform_exposures_to_tangent_plane(exposures: &Vec<Exposure>, center: Vector3<f64>) -> Vec<TangentPlaneExposure> {
    // do stuff
    exposures.iter().map(|e| e.transform_to_tangent_plane(center)).collect()
}

pub fn main() ->  Result<(), Box<dyn std::error::Error>> {
    // We need to make this more flexible.
    let args = Cli::try_parse()?;
    let f = std::fs::File::open(args.config)?;
    let config: Config = serde_yaml::from_reader(f)?;

    let mut kernel = SpiceKernel::new();
    kernel.load_spk(format!("{}/sb441-n16.bsp", config.spice_path).as_str())?;
    kernel.load_spk(format!("{}/de440s.bsp", config.spice_path).as_str())?;
    kernel.load_bpc(format!("{}/earth_1962_240827_2124_combined.bpc", config.spice_path).as_str())?;

    let spacerock_origin = Origin::from_str("ssb")?;
    let mu = spacerock_origin.mu();
    println!("Loaded SPICE kernels. Solar system barycenter mu: {}", mu);

    let mut detections = load_detections(&config.detection_catalog, &config.orbit_reference_plane, &kernel)?;
    let exposures = sort_detections_into_exposures(&detections);
    println!("Loaded {} exposures", exposures.len());

    let mut sum = 0;
    for exposure in exposures.iter() {
        println!("Exposure {} has {} detections at epoch {}", exposure.id, exposure.detections.len(), exposure.epoch);
        sum += exposure.detections.len();
    }  
    println!("Total detections: {}", sum);

    //let ref_vec = exposures[0].detections[0];
    // Calculate the reference vector as the mean of the detections in all exposures.
    let mut ref_vec = Vector3::zeros();
    let mut count = 0;
    for exposure in exposures.iter() {
        for det in exposure.detections.iter() {
            ref_vec += det;
            count += 1;
        }
    }
    ref_vec /= count as f64;
    ref_vec = ref_vec.normalize();

    let mut tangent_exposures = Vec::new();

    for exposure in exposures.iter() {
        let tangent_exposure = exposure.transform_to_tangent_plane(ref_vec);
        tangent_exposures.push(tangent_exposure);
    }

    let ref_epoch = config.reference_epoch;
    let gamma = 1.0/40.0; 
    let z0 = 1.0/gamma;
    let gdot = 0.0; 
    let GMtot = mu;

    let mm_sqr = (GMtot * gamma * gamma * gamma);

    // Calculate the mean theta_x and theta_y across all exposures.
    let mut mean_theta_x = 0.0;
    let mut mean_theta_y = 0.0;
    let mut count = 0;
    for exposure in tangent_exposures.iter() {
        for tx in exposure.theta_x.iter() {
            mean_theta_x += tx;
            count += 1;
        }
        for ty in exposure.theta_y.iter() {
            mean_theta_y += ty;
        }
    }
    mean_theta_x /= count as f64;
    mean_theta_y /= count as f64;

    let adot = 24.0 * (0.5/3600.0) * std::f64::consts::PI/180.0; // 0.5 arcsec/day in radians/day 
    let bdot = 24.0 * (0.1/3600.0) * std::f64::consts::PI/180.0; // 0.1 arcsec/day in radians/day
    

    let mut alphas = Vec::with_capacity(count);
    let mut betas = Vec::with_capacity(count);

    println!("Calculating alpha and beta for each detection...");
    let start = std::time::Instant::now();
    
    for exposure in tangent_exposures.iter() {
        let rho2 = (1.0 + mean_theta_x * mean_theta_x + mean_theta_y * mean_theta_y)*(z0 - exposure.xyz_e.z).powi(2);
        let rho = rho2.sqrt();

        // Improve light time correction at some point.
        let dt = rho/SPEED_OF_LIGHT;

        let t = exposure.epoch - ref_epoch;
        let tp = t - dt;
        let f = (1.0 - 0.5 * mm_sqr * tp * tp);
        let g = tp;
        let fac = (1.0 + g/f * gdot - gamma/f * exposure.xyz_e.z);
        for (theta_x, theta_y) in exposure.theta_x.iter().zip(exposure.theta_y.iter()) {
            let phi_x = theta_x * fac + gamma/f * exposure.xyz_e.x;
            let phi_y = theta_y * fac + gamma/f * exposure.xyz_e.y;
            let alpha = phi_x - g/f * adot;
            let beta = phi_y - g/f * bdot;
            //println!("Exposure {}: alpha = {}, beta = {}", exposure.id, alpha, beta);
            alphas.push(alpha);
            betas.push(beta);
        }
    }
    let duration = start.elapsed();
    println!("Finished calculating alpha and beta for {} detections in {:?}.", alphas.len(), duration);

    println!("Calculated alpha and beta for {} detections", alphas.len());

    // Build the grid and save it as a PNG
    // set pixel_scale to be 30 arcseconds per pixel in radians
    let pixel_scale = 30.0 * (std::f64::consts::PI/180.0) / 3600.0;
    
    let grid = build_grid(&alphas, &betas, pixel_scale);
    save_grid_png(&grid, "tangent_plane.png")?;

    // For a set of initial conditions:
        // Find peaks
        // Save the peaks to a file for later analysis.

    
    
   
    Ok(())

    

}