use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Adding in a struct to see how the data is refined through the pipeline
#[derive(Debug)]
pub struct PipelineMetrics {
    cluster_count: AtomicUsize,
    initial_trajectory_count: AtomicUsize,
    post_epsilon_trajectory_count: AtomicUsize,
    post_orbitfit_trajectory_count: AtomicUsize,
    final_trajectory_count: AtomicUsize,
}

impl PipelineMetrics {
    pub fn new() -> Self {
        Self {
            cluster_count: AtomicUsize::new(0),
            initial_trajectory_count: AtomicUsize::new(0),
            post_epsilon_trajectory_count: AtomicUsize::new(0),
            post_orbitfit_trajectory_count: AtomicUsize::new(0),
            final_trajectory_count: AtomicUsize::new(0),
        }
    }
    
    // Just increment counters - no deduplication
    pub fn track_cluster(&self, _detection_ids: &[usize]) {
        self.cluster_count.fetch_add(1, Ordering::Relaxed);
    }
    
    pub fn track_initial_trajectory(&self, _detection_ids: &[usize]) {
        self.initial_trajectory_count.fetch_add(1, Ordering::Relaxed);
    }
    
    pub fn track_post_epsilon_trajectory(&self, _detection_ids: &[usize]) {
        self.post_epsilon_trajectory_count.fetch_add(1, Ordering::Relaxed);
    }

    pub fn track_post_orbitfit_trajectory(&self, _detection_ids: &HashSet<usize>) {
        self.post_orbitfit_trajectory_count.fetch_add(1, Ordering::Relaxed);
    }
    
    pub fn track_final_trajectory(&self, _detection_ids: &HashSet<usize>) {
        self.final_trajectory_count.fetch_add(1, Ordering::Relaxed);
    }
    
    pub fn print_stage_summary(&self, stage_name: &str) {
        println!("\n{}", stage_name);
        println!("1. Total clusters: {}", self.cluster_count.load(Ordering::Relaxed));
        println!("2. Total initial trajectories: {}", self.initial_trajectory_count.load(Ordering::Relaxed));
        println!("3. Total post-epsilon cascade trajectories (1\"): {}", self.post_epsilon_trajectory_count.load(Ordering::Relaxed));
        println!("4. Total post-orbit-fit trajectories: {}", self.post_orbitfit_trajectory_count.load(Ordering::Relaxed));
        println!("5. Total final trajectories (post orbit-fit & non-det prob): {}\n", self.final_trajectory_count.load(Ordering::Relaxed));
    }
}