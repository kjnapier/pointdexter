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

    let adot = 0.1;
    let bdot = 0.0;

    for exposure in tangent_exposures.iter() {
        let rho2 = (1.0 + mean_theta_x * mean_theta_x + mean_theta_y * mean_theta_y)*(z0 - exposure.xyz_e.z).powi(2);
        let rho = rho2.sqrt();
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
            println!("Exposure {}: alpha = {}, beta = {}", exposure.id, alpha, beta);
        }
    }
    
   
    Ok(())

    

}