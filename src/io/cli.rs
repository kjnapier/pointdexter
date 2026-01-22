use clap::Parser;
use std::path::PathBuf;


/// Search for a pattern in a file and display the lines that contain it.
#[derive(Parser, Debug)]
pub struct Cli {
    pub config: PathBuf,
}

#[derive(serde::Deserialize, Debug)]
pub struct Config {

    pub epsilon_arcsec: f64,
    pub min_detections: usize,
    pub max_detections: usize,
    pub min_duration: f64,
    pub min_nites: usize,

    pub spice_path: String,
    pub detection_catalog: String,

    pub initial_conditions_file: String,
    pub reference_epoch: f64,
    pub orbit_reference_plane: String,
    pub ic_type: String,
    pub ic_origin: String,

    pub output_path: String,
    pub clusters_filename: String,

}