use clap::Parser;

/// Search for a pattern in a file and display the lines that contain it.
#[derive(Parser)]
pub struct Cli {
    pub config: std::path::PathBuf,
}

#[derive(serde::Deserialize, Debug)]
pub struct Config {
    pub spice_path: String,
    pub detection_catalog: String,
    pub all_detections_catalog: String,
    pub all_x3_exps_file: String,
    pub initial_conditions: String,
    pub fake_catalog: String,
    pub output_path: String,
    pub links_filename: String,
    pub clusters_filename: String,
    pub epsilon_arcsec: f64,
    pub min_detections: usize,
    pub max_detections: usize,
    pub min_duration: f64,
    pub min_nites: usize,
    pub t_max: f64,
    pub trajectories_filename: String,
    pub invalid_expnums: Vec<i64>,
    pub invalid_ccds: Vec<i32>,

    pub eta0: f64,
    pub m50: f64,
    pub sigma: f64,
    pub prob_threshold: f64,
    pub non_detection: bool,
    pub t_ref: f64,
    pub batch_size: usize,
}