use crate::chebyshev::{chebyshev_eval, fit_chebyshev_direct};
use spacerocks::transforms::calc_true_anomaly_from_mean_anomaly;
use spacerocks::time::Time;
use spacerocks::transforms::{solve_for_universal_anomaly, stumpff_c, stumpff_s};
use ordered_float::OrderedFloat;

use crate::sparse_linking_stuff::grid::IcBoundsKey;

fn nice_acos(x: f64) -> f64 {
    x.clamp(-1.0, 1.0).acos()
}

#[derive(Debug, Clone, Copy)]
pub enum InputFormat {
    KEV, // (r, vr, vo, psi)
    QEF, // (q, e, f, psi)
    KEP, // (a, e, f, psi)
}


#[derive(Debug, Clone)]
pub struct InitialCondition {
    pub id: String,
    pub r: f64,
    pub vr: f64,
    pub vo: f64,
    pub inc: f64,
    pub kappa: i32,
    pub epoch: f64,
    pub mu: f64,
    pub h: f64,
    pub alpha: f64,
    pub interpolation_bounds: f64,
    pub r_poly: Vec<f64>,
    pub vr_poly: Vec<f64>,
    pub f_poly: Vec<f64>,
    pub g_poly: Vec<f64>,
    pub cos_inc: f64,
    pub sin_latitude_threshold: f64,
    pub energy: f64,
    pub q: f64,
    pub e: f64,
    pub true_anomaly: f64,

    // bounds
    pub r_min: Option<f64>,
    pub r_max: Option<f64>,
    pub vr_min: Option<f64>,
    pub vr_max: Option<f64>,
    pub vo_min: Option<f64>,
    pub vo_max: Option<f64>,
    pub inc_min: Option<f64>,
    pub inc_max: Option<f64>,
}

impl InitialCondition {
    pub fn from_elements(
        id: String,
        q: f64,
        e: f64,
        inc: f64,
        true_anomaly: f64,
        kappa: i32,
        epoch: Time,
        mu: f64,
        r_min: Option<f64>, r_max: Option<f64>,
        vr_min: Option<f64>, vr_max: Option<f64>,
        vo_min: Option<f64>, vo_max: Option<f64>,
        inc_min: Option<f64>, inc_max: Option<f64>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if e > 1.0 {
            if true_anomaly.abs() > nice_acos(-1.0 / e) {
                panic!("True anomaly out of bounds for eccentricity.");
            }
        }

        let p = q * (1.0 + e);
        let h = (p * mu).sqrt();
        let r0 = p / (1.0 + e * true_anomaly.cos());
        let vo = h / r0;
        let vr = (mu / h) * e * true_anomaly.sin();

        Self::from_spherical(
            id, r0, vr, vo, inc, kappa, epoch, mu,
            r_min, r_max, vr_min, vr_max, vo_min, vo_max, inc_min, inc_max,
        )
    }

    pub fn from_keplerian(
        id: String,
        q: f64,
        e: f64,
        inc: f64,
        mean_anomaly: f64,
        kappa: i32,
        epoch: Time,
        mu: f64,
        r_min: Option<f64>, r_max: Option<f64>,
        vr_min: Option<f64>, vr_max: Option<f64>,
        vo_min: Option<f64>, vo_max: Option<f64>,
        inc_min: Option<f64>, inc_max: Option<f64>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let true_anomaly = calc_true_anomaly_from_mean_anomaly(e, mean_anomaly).unwrap();

        let p = q * (1.0 + e);
        let h = (p * mu).sqrt();
        let r0 = p / (1.0 + e * true_anomaly.cos());
        let vo = h / r0;
        let vr = (mu / h) * e * true_anomaly.sin();

        Self::from_spherical(
            id, r0, vr, vo, inc, kappa, epoch, mu,
            r_min, r_max, vr_min, vr_max, vo_min, vo_max, inc_min, inc_max,
        )
    }

    pub fn from_spherical(
        id: String,
        r: f64,
        vr: f64,
        vo: f64,
        inc: f64,
        kappa: i32,
        epoch: Time,
        mu: f64,
        r_min: Option<f64>, r_max: Option<f64>,
        vr_min: Option<f64>, vr_max: Option<f64>,
        vo_min: Option<f64>, vo_max: Option<f64>,
        inc_min: Option<f64>, inc_max: Option<f64>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let epoch_jd = epoch.tdb().jd();
        let h = r * vo;
        let energy = 0.5 * (vr * vr + vo * vo) - mu / r;
        let alpha = -2.0 * energy / mu;

        let radicand = 1.0 + 2.0 * energy * h.powi(2) / mu.powi(2);
        let e = radicand.max(0.0).sqrt();

        let vsq = vr.powi(2) + vo.powi(2);
        let a = 1.0 / (2.0 / r - vsq / mu);
        let q = a * (1.0 - e);

        let cosf = (r * vo.powi(2) / mu - 1.0) / e;
        let mut true_anomaly = nice_acos(cosf);

        let interpolation_bounds = 365.25;

        let dt_values: Vec<f64> = (-(interpolation_bounds as i32)..=(interpolation_bounds as i32))
            .step_by(7)
            .map(|v| v as f64)
            .collect();

        let mut s_values = Vec::new();
        for &dt in &dt_values {
            let s = solve_for_universal_anomaly(r, vr, alpha, mu, dt, 1e-12, 100)?;
            s_values.push(s);
        }

        let stumpff_c_values: Vec<f64> =
            s_values.iter().map(|&s| stumpff_c(alpha * s * s)).collect();
        let stumpff_s_values: Vec<f64> =
            s_values.iter().map(|&s| stumpff_s(alpha * s * s)).collect();

        let f_values: Vec<f64> = s_values
            .iter()
            .zip(stumpff_c_values.iter())
            .map(|(&s, &c)| 1.0 - s * s / r * c)
            .collect();

        let g_values: Vec<f64> = s_values
            .iter()
            .zip(stumpff_s_values.iter())
            .enumerate()
            .map(|(i, (&s, &sval))| dt_values[i] - s.powi(3) / mu.sqrt() * sval)
            .collect();

        let rsq_values: Vec<f64> = f_values
            .iter()
            .zip(g_values.iter())
            .map(|(&f, &g)| {
                f.powi(2) * r.powi(2) + 2.0 * f * g * r * vr + g.powi(2) * (vo.powi(2) + vr.powi(2))
            })
            .collect();

        let r_values: Vec<f64> = rsq_values.iter().map(|&val| val.sqrt()).collect();

        let fdot_values: Vec<f64> = s_values
            .iter()
            .zip(r_values.iter())
            .zip(stumpff_s_values.iter())
            .map(|((&s, &rval), &sval)| s * mu.sqrt() / (r * rval) * (alpha * s * s * sval - 1.0))
            .collect();

        let gdot_values: Vec<f64> = s_values
            .iter()
            .zip(r_values.iter())
            .zip(stumpff_c_values.iter())
            .map(|((&s, &rval), &c)| 1.0 - s * s / rval * c)
            .collect();

        let r_vr_values: Vec<f64> = f_values
            .iter()
            .zip(fdot_values.iter())
            .map(|(&f, &fdot)| f * fdot * r * r)
            .zip(
                g_values
                    .iter()
                    .zip(gdot_values.iter())
                    .map(|(&g, &gdot)| g * gdot * (vo * vo + vr * vr)),
            )
            .zip(
                f_values
                    .iter()
                    .zip(gdot_values.iter())
                    .zip(g_values.iter().zip(fdot_values.iter()))
                    .map(|((&f, &gdot), (&g, &fdot))| r * vr * (f * gdot + g * fdot)),
            )
            .map(|((a, b), c)| a + b + c)
            .collect();

        let vr_values: Vec<f64> = r_vr_values
            .iter()
            .zip(r_values.iter())
            .map(|(&rv, &rval)| rv / rval)
            .collect();

        let scaled_dt_values: Vec<f64> = dt_values
            .iter()
            .map(|&dt| dt / interpolation_bounds)
            .collect();

        let r_poly  = fit_chebyshev_direct(&scaled_dt_values, &r_values,  3);
        let vr_poly = fit_chebyshev_direct(&scaled_dt_values, &vr_values, 3);
        let f_poly  = fit_chebyshev_direct(&scaled_dt_values, &f_values,  3);
        let g_poly  = fit_chebyshev_direct(&scaled_dt_values, &g_values,  3);

        let mut latitude_threshold = inc;
        if latitude_threshold > std::f64::consts::PI / 2.0 {
            latitude_threshold = std::f64::consts::PI - latitude_threshold;
        }
        let sin_latitude_threshold = latitude_threshold.sin();

        Ok(InitialCondition {
            id,
            r,
            vr,
            vo,
            inc,
            kappa,
            epoch: epoch_jd,
            mu,
            h,
            alpha,
            interpolation_bounds,
            r_poly,
            vr_poly,
            f_poly,
            g_poly,
            cos_inc: inc.cos(),
            sin_latitude_threshold,
            energy,
            q,
            e, 
            true_anomaly,
            r_min,  r_max,
            vr_min, vr_max,
            vo_min, vo_max,
            inc_min, inc_max,
        })
    }

    pub fn vr_at_epoch(&self, epoch: f64) -> f64 {
        let scaled_dt = (epoch - self.epoch) / self.interpolation_bounds;
        chebyshev_eval(&self.vr_poly, scaled_dt)
    }

    pub fn r_at_epoch(&self, epoch: f64) -> f64 {
        let scaled_dt = (epoch - self.epoch) / self.interpolation_bounds;
        chebyshev_eval(&self.r_poly, scaled_dt)
    }

    pub fn r_and_vr_at_epoch(&self, epoch: f64) -> (f64, f64) {
        let scaled_dt = (epoch - self.epoch) / self.interpolation_bounds;
        (chebyshev_eval(&self.r_poly, scaled_dt), chebyshev_eval(&self.vr_poly, scaled_dt))
    }

    pub fn fg_at_epoch(&self, epoch: f64) -> (f64, f64) {
        let scaled_dt = (epoch - self.epoch) / self.interpolation_bounds;
        (chebyshev_eval(&self.f_poly, scaled_dt), chebyshev_eval(&self.g_poly, scaled_dt))
    }

    pub fn f_at_epoch(&self, epoch: f64) -> f64 {
        let scaled_dt = (epoch - self.epoch) / self.interpolation_bounds;
        chebyshev_eval(&self.f_poly, scaled_dt)
    }

    pub fn g_at_epoch(&self, epoch: f64) -> f64 {
        let scaled_dt = (epoch - self.epoch) / self.interpolation_bounds;
        chebyshev_eval(&self.g_poly, scaled_dt)
    }



    pub fn to_bounds_key(&self) -> Option<IcBoundsKey> {
        Some(IcBoundsKey::KEV {
            r_bounds:   (OrderedFloat(self.r_min?),   OrderedFloat(self.r_max?)),
            vr_bounds:  (OrderedFloat(self.vr_min?),  OrderedFloat(self.vr_max?)),
            vo_bounds:  (OrderedFloat(self.vo_min?),  OrderedFloat(self.vo_max?)),
            psi_bounds: (OrderedFloat(self.inc_min?), OrderedFloat(self.inc_max?)),
            epoch:      OrderedFloat(self.epoch),
        })
    }
}