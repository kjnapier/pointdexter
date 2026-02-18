use argmin::solver::particleswarm::ParticleSwarm;
use argmin::{
    core::{observers::ObserverMode},
};
use argmin_observer_slog::SlogLogger;

use std::sync::Arc;
use argmin::core::{CostFunction, Error};

use anyhow::Result;
use anyhow;


#[derive(Clone)]
struct ClusterProblem {
    detections: Arc<Vec<Detection>>,
    ic0: InitialCondition,
}

impl ClusterProblem {
    fn ic_from_p(&self, p: &Vec<f64>) -> InitialCondition {
        debug_assert!(p.len() == 4);
        // let mut ic = self.ic0.clone();
        let ic = InitialCondition::from_spherical(
            self.ic0.id.clone(),
            p[0],   // r
            p[1],   // vr
            p[2],   // vo
            p[3],   // inc
            self.ic0.kappa,
            Time::new(self.ic0.epoch, "tdb", "jd").unwrap(),
            self.ic0.mu,
        ).unwrap();
        
        ic
    }
}

impl CostFunction for ClusterProblem {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, p: &Self::Param) -> std::result::Result<Self::Output, Error> {
        let ic = self.ic_from_p(p);
        Ok(robust_cluster_radius(&self.detections, &ic)) // 3 arcsec in radians
    }
}



use argmin::core::Executor;
use argmin::solver::neldermead::NelderMead;

fn optimize_ic(
    detections: Arc<Vec<Detection>>,
    ic_guess: InitialCondition,
    p0: [f64; 4],
    step: [f64; 4],
) -> Result<(InitialCondition, f64)> {
    let problem = ClusterProblem { detections: detections.clone(), ic0: ic_guess.clone() };

    let x0 = vec![p0[0], p0[1], p0[2], p0[3]];

    let simplex: Vec<Vec<f64>> = vec![
        x0.clone(),
        vec![x0[0] + step[0], x0[1],           x0[2],           x0[3]],
        vec![x0[0],           x0[1] + step[1], x0[2],           x0[3]],
        vec![x0[0],           x0[1],           x0[2] + step[2], x0[3]],
        vec![x0[0],           x0[1],           x0[2],           x0[3] + step[3]],
    ];
    
    let solver = NelderMead::new(simplex)
        .with_sd_tolerance(1e-12)?;

    let res = Executor::new(problem, solver)
        .configure(|state| state.max_iters(10_000))
        .add_observer(SlogLogger::term(), ObserverMode::Always)
        .run()
        .map_err(|e| anyhow::anyhow!(e.to_string()))?;

    let best_p = res.state().best_param.as_ref().unwrap();
    let best_cost = res.state().best_cost;
    let mut best_ic = InitialCondition::from_spherical(
        ic_guess.id.clone(),
        best_p[0],
        best_p[1],
        best_p[2],
        best_p[3],
        ic_guess.kappa,
        Time::new(ic_guess.epoch, "tdb", "jd").unwrap(),
        ic_guess.mu,
    ).unwrap();

    Ok((best_ic, best_cost))
}

/// Robust "cluster radius":
/// 1) robust center via trimmed mean on the sphere
/// 2) return median angular distance to that center
pub fn robust_cluster_radius(detections: &Vec<Detection>, ic: &InitialCondition) -> f64 {
    let synced_points = sync_detections_to_orbit(detections, ic);

    let mut points: Vec<[f64; 3]> = Vec::new();
    for p in synced_points.iter() {
        if let Some(v) = p {
            points.push(*v);
        }
    }
    if points.len() < 10 {
        return f64::INFINITY;
    }

    // calculate the barycenter of the cluster, using a median
    let mut xs: Vec<f64> = points.iter().map(|v| v[0]).collect();
    let mut ys: Vec<f64> = points.iter().map(|v| v[1]).collect();
    let mut zs: Vec<f64> = points.iter().map(|v| v[2]).collect();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
    zs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = points.len() / 2;
    let c: [f64; 3] = if points.len() % 2 == 0 {
        [
            (xs[mid - 1] + xs[mid]) / 2.0,
            (ys[mid - 1] + ys[mid]) / 2.0,
            (zs[mid - 1] + zs[mid]) / 2.0,
        ]
    } else {
        [
            xs[mid],
            ys[mid],
            zs[mid],
        ]
    };
    // normalize c
    let norm = (c[0]*c[0] + c[1]*c[1] + c[2]*c[2]).sqrt();
    let c: [f64; 3] = [c[0]/norm, c[1]/norm, c[2]/norm];
    

    // get the standard deviation of angles to the center
    let mut angles: Vec<f64> = Vec::new();
    for v in points.iter() {
        let dx = v[0] - c[0];
        let dy = v[1] - c[1];
        let dz = v[2] - c[2];
        let d2 = dx*dx + dy*dy + dz*dz;
        let angle = squared_euclid_to_angle_rad(d2);
        angles.push(angle);
    }

    // calculate the sigma-clipped mean and stddev
    angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = angles.len();
    let n_sig = 5.0;
    let mut clipped_angles: Vec<f64> = angles.clone();
    for _ in 0..10 {
        let mean_angle = clipped_angles.iter().sum::<f64>() / (clipped_angles.len() as f64);
        let std_angle = (clipped_angles.iter().map(|&a| (a - mean_angle).powi(2)).sum::<f64>() / (clipped_angles.len() as f64)).sqrt();
        let upper_bound = mean_angle + n_sig * std_angle;
        clipped_angles = clipped_angles.into_iter().filter(|&a| a <= upper_bound).collect();
    }

    let robust_mean_angle = clipped_angles.iter().sum::<f64>() / (clipped_angles.len() as f64);
    let angle_diffs: Vec<f64> = clipped_angles.iter().map(|&a| (a - robust_mean_angle).abs()).collect();

    // the astrometric uncertaity is about 0.3 arcseconds. The posterior should be gaussian, centered at the robust mean angle.
    // calculate the likelihood of the robust mean angle given the astrometric uncertainty
    let astrometric_uncertainty = 0.3 / 3600.0 * std::f64::consts::PI / 180.0; // in radians
    let variance = astrometric_uncertainty.powi(2);
    let log_likelihood: f64 = angle_diffs.iter().map(|&a| {
        -0.5 * (a.powi(2) / variance) - 0.5 * (2.0 * std::f64::consts::PI * variance).ln()
    }).sum();
    -log_likelihood
    
}