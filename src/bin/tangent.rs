use pointdexter::cli::{Config, Cli};
use pointdexter::Exposure;
use pointdexter::load_detections;

use spacerocks::{SpiceKernel};
use pointdexter::Detection;
use pointdexter::TangentPlaneExposure;
use nalgebra::Vector3;
use clap::Parser;


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

    let mut detections = load_detections(&config.detection_catalog, &config.orbit_reference_plane, &kernel)?;
    let exposures = sort_detections_into_exposures(&detections);
    println!("Loaded {} exposures", exposures.len());



    Ok(())

    

}