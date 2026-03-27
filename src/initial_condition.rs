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
    KEV, // (r, vr, vo, inc)
    KEP, // (q, e, f, inc)
}

#[derive(Debug, Clone)]
pub struct ElementBounds {
    pub format: InputFormat,
    pub bounds: Vec<(Option<f64>, Option<f64>)>,
}

impl ElementBounds {
    pub fn kev(
        r: (Option<f64>, Option<f64>),
        vr: (Option<f64>, Option<f64>),
        vo: (Option<f64>, Option<f64>),
        inc: (Option<f64>, Option<f64>),
    ) -> Self {
        Self {
            format: InputFormat::KEV,
            bounds: vec![r, vr, vo, inc],
        }
    }

    pub fn kep(
        q: (Option<f64>, Option<f64>),
        e: (Option<f64>, Option<f64>),
        f: (Option<f64>, Option<f64>),
        inc: (Option<f64>, Option<f64>),
    ) -> Self {
        Self {
            format: InputFormat::KEP,
            bounds: vec![q, e, f, inc],
        }
    }

    pub fn get(&self, idx: usize) -> (Option<f64>, Option<f64>) {
        self.bounds.get(idx).copied().unwrap_or((None, None))
    }

    pub fn input_format(&self) -> InputFormat {
        self.format
    }
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
    pub bounds: Option<ElementBounds>,

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
        bounds: Option<ElementBounds>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if e > 1.0 && true_anomaly.abs() > nice_acos(-1.0 / e) {
            panic!("True anomaly out of bounds for eccentricity.");
        }

        let p = q * (1.0 + e);
        let h = (p * mu).sqrt();
        let r0 = p / (1.0 + e * true_anomaly.cos());
        let vo = h / r0;
        let vr = (mu / h) * e * true_anomaly.sin();

        Self::from_spherical_inner(id, r0, vr, vo, inc, kappa, epoch, mu, bounds)
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
        bounds: Option<ElementBounds>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let true_anomaly = calc_true_anomaly_from_mean_anomaly(e, mean_anomaly)?;

        let p = q * (1.0 + e);
        let h = (p * mu).sqrt();
        let r0 = p / (1.0 + e * true_anomaly.cos());
        let vo = h / r0;
        let vr = (mu / h) * e * true_anomaly.sin();

        Self::from_spherical_inner(id, r0, vr, vo, inc, kappa, epoch, mu, bounds)
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
        bounds: Option<ElementBounds>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        Self::from_spherical_inner(id, r, vr, vo, inc, kappa, epoch, mu, bounds)
    }

    fn from_spherical_inner(
        id: String,
        r: f64,
        vr: f64,
        vo: f64,
        inc: f64,
        kappa: i32,
        epoch: Time,
        mu: f64,
        bounds: Option<ElementBounds>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let epoch_jd = epoch.tdb().jd();
        let h = r * vo;
        let energy = 0.5 * (vr * vr + vo * vo) - mu / r;
        let alpha = -2.0 * energy / mu;

        let e = (1.0 + 2.0 * energy * h.powi(2) / mu.powi(2)).max(0.0).sqrt();
        let vsq = vr.powi(2) + vo.powi(2);
        let a = 1.0 / (2.0 / r - vsq / mu);
        let q = a * (1.0 - e);
        let true_anomaly = nice_acos((r * vo.powi(2) / mu - 1.0) / e);

        let interpolation_bounds = 365.25;

        let dt_values: Vec<f64> = (-(interpolation_bounds as i32)..=(interpolation_bounds as i32))
            .step_by(7)
            .map(|v| v as f64)
            .collect();

        let s_values: Vec<f64> = dt_values
            .iter()
            .map(|&dt| solve_for_universal_anomaly(r, vr, alpha, mu, dt, 1e-7, 100))
            .collect::<Result<Vec<_>, _>>()?;

        let stumpff_c_values: Vec<f64> = s_values.iter().map(|&s| stumpff_c(alpha * s * s)).collect();
        let stumpff_s_values: Vec<f64> = s_values.iter().map(|&s| stumpff_s(alpha * s * s)).collect();

        let f_values: Vec<f64> = s_values
            .iter()
            .zip(&stumpff_c_values)
            .map(|(&s, &c)| 1.0 - s * s / r * c)
            .collect();

        let g_values: Vec<f64> = dt_values
            .iter()
            .zip(&s_values)
            .zip(&stumpff_s_values)
            .map(|((&dt, &s), &sv)| dt - s.powi(3) / mu.sqrt() * sv)
            .collect();

        let r_values: Vec<f64> = f_values
            .iter()
            .zip(&g_values)
            .map(|(&f, &g)| {
                (f.powi(2) * r.powi(2)
                    + 2.0 * f * g * r * vr
                    + g.powi(2) * (vr.powi(2) + vo.powi(2)))
                .sqrt()
            })
            .collect();

        let fdot_values: Vec<f64> = s_values
            .iter()
            .zip(&r_values)
            .zip(&stumpff_s_values)
            .map(|((&s, &rv), &sv)| s * mu.sqrt() / (r * rv) * (alpha * s * s * sv - 1.0))
            .collect();

        let gdot_values: Vec<f64> = s_values
            .iter()
            .zip(&r_values)
            .zip(&stumpff_c_values)
            .map(|((&s, &rv), &c)| 1.0 - s * s / rv * c)
            .collect();

        let vr_values: Vec<f64> = f_values
            .iter()
            .zip(&fdot_values)
            .zip(g_values.iter().zip(&gdot_values))
            .zip(&r_values)
            .map(|(((&f, &fdot), (&g, &gdot)), &rv)| {
                (f * fdot * r.powi(2)
                    + g * gdot * (vr.powi(2) + vo.powi(2))
                    + r * vr * (f * gdot + g * fdot))
                    / rv
            })
            .collect();

        let scaled_dt_values: Vec<f64> = dt_values
            .iter()
            .map(|&dt| dt / interpolation_bounds)
            .collect();

        let r_poly  = fit_chebyshev_direct(&scaled_dt_values, &r_values,  3);
        let vr_poly = fit_chebyshev_direct(&scaled_dt_values, &vr_values, 3);
        let f_poly  = fit_chebyshev_direct(&scaled_dt_values, &f_values,  3);
        let g_poly  = fit_chebyshev_direct(&scaled_dt_values, &g_values,  3);

        let latitude_threshold = if inc > std::f64::consts::PI / 2.0 {
            std::f64::consts::PI - inc
        } else {
            inc
        };

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
            sin_latitude_threshold: latitude_threshold.sin(),
            energy,
            q,
            e,
            true_anomaly,
            bounds,
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
        let b = self.bounds.as_ref()?;
        match b.format {
            InputFormat::KEV => Some(IcBoundsKey::KEV {
                r_bounds:   (OrderedFloat(b.get(0).0?), OrderedFloat(b.get(0).1?)),
                vr_bounds:  (OrderedFloat(b.get(1).0?), OrderedFloat(b.get(1).1?)),
                vo_bounds:  (OrderedFloat(b.get(2).0?), OrderedFloat(b.get(2).1?)),
                inc_bounds: (OrderedFloat(b.get(3).0?), OrderedFloat(b.get(3).1?)),
                epoch:      OrderedFloat(self.epoch),
            }),
            InputFormat::KEP => Some(IcBoundsKey::KEP {
                q_bounds:   (OrderedFloat(b.get(0).0?), OrderedFloat(b.get(0).1?)),
                e_bounds:   (OrderedFloat(b.get(1).0?), OrderedFloat(b.get(1).1?)),
                f_bounds:   (OrderedFloat(b.get(2).0?), OrderedFloat(b.get(2).1?)),
                inc_bounds: (OrderedFloat(b.get(3).0?), OrderedFloat(b.get(3).1?)),
                epoch:      OrderedFloat(self.epoch),
            }),
        }
    }
}