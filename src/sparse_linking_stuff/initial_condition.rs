// use spacerocks::time::Time;
// use spacerocks::transforms::{solve_for_universal_anomaly, stumpff_s, stumpff_c, calc_conic_anomaly_from_true_anomaly, calc_mean_anomaly_from_conic_anomaly};
// use spacerocks::OrbitType;

// use std::f64::consts::PI;

// use crate::grid::IcBoundsKey;
// use ordered_float::OrderedFloat;

// pub fn nice_acos(x: f64) -> f64 {
//     if x > 1.0 {
//         0.0
//     } else if x < -1.0 {
//         PI
//     } else {
//         x.acos()
//     }
// }


// #[derive(Debug, Clone)]
// pub struct InitialCondition {
//     pub id: usize,
//     pub q_min: Option<f64>,
//     pub q_max: Option<f64>,
//     pub q: f64,
//     pub e_min: Option<f64>,
//     pub e_max: Option<f64>,
//     pub e: f64,
//     pub psi_min: Option<f64>,
//     pub psi_max: Option<f64>,
//     pub psi: f64,
//     pub true_anomaly_min: Option<f64>,
//     pub true_anomaly_max: Option<f64>,
//     pub true_anomaly: f64,
//     pub epoch: Time,

//     pub h: f64,
//     pub energy: f64,
//     pub mean_anomaly: f64,

//     pub r0: f64,
//     pub vr: f64,
//     pub vo: f64,

//     pub mu: f64,
//     pub orbit_type: OrbitType,
// }

// impl InitialCondition {
//     pub fn new(id: usize, q_min: Option<f64>, q_max: Option<f64>, q: f64, 
//                 e_min: Option<f64>, e_max: Option<f64>, e: f64, 
//                 true_anomaly_min: Option<f64>, true_anomaly_max: Option<f64>, true_anomaly: f64, 
//                 psi_min: Option<f64>, psi_max: Option<f64>, psi: f64, epoch: Time, mu: f64) -> Self {
//         let p = q * (1.0 + e); 
//         let h = (mu * p).sqrt();
//         let energy = -0.5 *(mu.powi(2) / h.powi(2)) * (1.0 - e.powi(2));
//         let r0 = p / (1.0 + e * true_anomaly.cos());
//         let vo = h / r0;
//         let vr = (mu / h) * e * true_anomaly.sin();

//         let eccentric_anomaly = calc_conic_anomaly_from_true_anomaly(e, true_anomaly).expect("Failed to calculate eccentric anomaly");
//         let mean_anomaly = calc_mean_anomaly_from_conic_anomaly(e, eccentric_anomaly).expect("Failed to calculate mean anomaly");


//         let orbit_type = OrbitType::from_eccentricity(e, 1e-10).expect("Failed to determine orbit type");
//         Self { 
//             id, 
//             q_min, q_max, q, 
//             e_min, e_max, e, 
//             true_anomaly_min, true_anomaly_max, true_anomaly, 
//             psi_min, psi_max, psi, 
//             epoch, 
//             h, energy, mean_anomaly, 
//             r0, vr, vo, 
//             mu, orbit_type 
//         }
//     }

//     pub fn chi(&self, epoch: Time) -> f64 {
//         let dt = epoch.epoch - self.epoch.epoch;

        
//         let s = solve_for_universal_anomaly(
//             self.r0, 
//             self.vr, 
//             self.alpha(), 
//             self.mu, 
//             dt,
//             1e-10, // was 1e-14
//             100 // Number of iterations 
//         ).unwrap_or_else(|_| panic!(
//             "Failed to solve for universal_anomaly\n\
//              r: {}, vr: {}, alpha: {}, mu: {}, dt: {}", 
//             self.r0, self.vr, self.alpha(), self.mu, dt
//         )); 

//         s
//     }

//     pub fn lagrange_f_and_g(&self, epoch: Time) -> (f64, f64) {
    
//         let dt = epoch.epoch - self.epoch.epoch;
//         let chi = self.chi(epoch);
        
//         // Calculate Stumpff functions
//         let c = stumpff_c(self.alpha() * chi.powi(2));
//         let s = stumpff_s(self.alpha() * chi.powi(2));
        
//         let f = 1.0 - chi.powi(2) * c / self.r0;
//         let g = dt - chi.powi(3) * s / self.mu.sqrt();

//         (f, g)
//     }

//     pub fn calculate_r(&self, epoch: Time) -> f64 {
        
//         // Get Lagrange coefficients
//         let (f,g) = self.lagrange_f_and_g(epoch);


//         // Calculate r**2 from Lagrange coefficients and v0 = vr^2 + vo^2
//         let rsq = f.powi(2) * self.r0.powi(2) + 
//                  2.0 * f * g * self.r0 * self.vr + 
//                  g.powi(2) * (self.vo.powi(2) + self.vr.powi(2));
        
//         rsq.sqrt()
//     }

//     pub fn a(&self) -> Option<f64> {
//         match self.orbit_type {
//             OrbitType::Hyperbolic => Some(self.q / (self.e - 1.0)),
//             OrbitType::Parabolic => None,
//             OrbitType::Elliptical => Some(self.q / (1.0 - self.e)),
//             OrbitType::Circular => Some(self.q),
//             _ => None,
//         }
//     }

//     pub fn n(&self) -> f64 {
//         match self.a() {
//             Some(a) => (self.mu / a.powi(3)).sqrt(),
//             None => (self.mu / self.q.powi(3)).sqrt(),
//         }
//     }

//     pub fn alpha(&self) -> f64 {
//         -2.0 * self.energy / self.mu
//     }


//     pub fn separation(&self, other: &InitialCondition, epoch: Time) -> f64 {

//         let self_r = self.calculate_r(epoch.clone());
//         let (self_f, self_g) = self.lagrange_f_and_g(epoch.clone());

//         let other_r = other.calculate_r(epoch.clone());
//         let (other_f, other_g) = other.lagrange_f_and_g(epoch.clone());

//         // Calculate common products
//         let r_prod = self.r0 * other.r0;
//         let r_new_prod = self_r * other_r;
//         let psi_cos = (self.psi - other.psi).cos();
        
//         // Calculate all terms maintaining original logic
//         let first_term = r_prod * self_f * other_f;
//         let second_term = self.r0 * other.vr * self_f * other_g;
//         let third_term = self.vr * other.r0 * self_g * other_f;
//         let fourth_term = self_g * other_g * (self.vr * other.vr + self.vo * other.vo * psi_cos);

//         // Combine and normalize
//         let argument = (first_term + second_term + third_term + fourth_term) / r_new_prod;
        
//         nice_acos(argument)
//     }
// }

// impl InitialCondition {
//     pub fn to_bounds_key(&self) -> Option<IcBoundsKey> {
//         Some(IcBoundsKey {
//             q_bounds: (OrderedFloat(self.q_min?), OrderedFloat(self.q_max?)),
//             e_bounds: (OrderedFloat(self.e_min?), OrderedFloat(self.e_max?)),
//             f_bounds: (OrderedFloat(self.true_anomaly_min?), OrderedFloat(self.true_anomaly_max?)),
//             psi_bounds: (OrderedFloat(self.psi_min?), OrderedFloat(self.psi_max?)),
//             epoch: OrderedFloat(self.epoch.epoch),
//         })
//     }
// }

// impl std::fmt::Display for InitialCondition {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         // write!(f, "r: {}, vr: {}, vo: {}, psi: {}, orbit_type: {:?}", 
//         //        self.r0, self.vr, self.vo, self.psi, self.orbit_type)
//         write!(f, "q: {}, e: {}, psi: {}, true_anomaly: {}, epoch: {}, 
//                     h: {}, energy: {}, mean_anomaly: {}, r0: {}, vr: {}, 
//                     vo: {}, mu: {}, orbit_type: {:?}", self.q, self.e, self.psi,
//                     self.true_anomaly, self.epoch.epoch, self.h, self.energy, self.mean_anomaly,
//                     self.r0, self.vr, self.vo, self.mu, self.orbit_type)
//     }
// }











///// NEW!!!!!!!!
use spacerocks::time::Time;
use spacerocks::transforms::{solve_for_universal_anomaly, stumpff_s, stumpff_c, calc_conic_anomaly_from_true_anomaly, calc_mean_anomaly_from_conic_anomaly};
use spacerocks::OrbitType;

use std::f64::consts::PI;

use crate::sparse_linking_stuff::grid::IcBoundsKey;
use ordered_float::OrderedFloat;

// use serde::{Serialize, Deserialize};

pub fn nice_acos(x: f64) -> f64 {
    if x > 1.0 {
        0.0
    } else if x < -1.0 {
        PI
    } else {
        x.acos()
    }
}

#[derive(Debug, Clone, Copy)]
pub enum InputFormat {
    KEV, // (r, vr, vo, psi)
    QEF, // (q, e, f, psi)
    KEP, // (a, e, f, psi)
}


#[derive(Debug, Clone)]
pub struct InitialCondition {
    pub id: usize,
    pub format: InputFormat,

    pub p1_min: Option<f64>,
    pub p1_max: Option<f64>,
    pub p1: f64,

    pub p2_min: Option<f64>,
    pub p2_max: Option<f64>,
    pub p2: f64,

    pub p3_min: Option<f64>,
    pub p3_max: Option<f64>,
    pub p3: f64,

    pub p4_min: Option<f64>,
    pub p4_max: Option<f64>,
    pub p4: f64,

    pub epoch: Time,
    pub mu: f64,

    // Derived quantities (standardized form)
    pub q: f64,
    pub e: f64,
    pub true_anomaly: f64,
    pub psi: f64,
    pub h: f64,
    pub energy: f64,
    pub mean_anomaly: f64,
    pub r0: f64,
    pub vr: f64,
    pub vo: f64,
    pub orbit_type: OrbitType,
}


impl InitialCondition {
    pub fn from_params(
        id: usize,
        format: InputFormat,
        p1_min: Option<f64>, p1_max: Option<f64>, p1: f64,
        p2_min: Option<f64>, p2_max: Option<f64>, p2: f64,
        p3_min: Option<f64>, p3_max: Option<f64>, p3: f64,
        p4_min: Option<f64>, p4_max: Option<f64>, p4: f64,
        epoch: Time,
        mu: f64,
    ) -> Self {
        let (q, e, true_anomaly, psi, h, energy, mean_anomaly, r0, vr, vo, orbit_type) =
            match format {
                InputFormat::KEV => {
                    let r = p1;
                    let vr = p2;
                    let vo = p3;
                    let psi = p4;
                    let h = r * vo;
                    let vsq = vr.powi(2) + vo.powi(2);
                    let energy = vsq / 2.0 - mu / r;
                    // Ensure that e is non-negative by taking max with 0.0 before sqrt for floating point
                    let radicand = 1.0 + 2.0 * energy * h.powi(2) / mu.powi(2);
                    let e = radicand.max(0.0).sqrt();
                    //let e = if (e - 1.0).abs() < 1e-8 { 1.0 - 1e-8 } else { e };
                    let a = 1.0 / (2.0 / r - vsq / mu);
                    let q = a * (1.0 - e);
                    let orbit_type = OrbitType::from_eccentricity(e, 1e-10).unwrap();
                    let cosf = (r * vo.powi(2) / mu - 1.0) / e;
                    let mut f = nice_acos(cosf);
                    if vr < 0.0 {
                        f = if orbit_type == OrbitType::Elliptical {
                            2.0 * PI - f
                        } else {
                            -f
                        };
                    }
                    let ecc = calc_conic_anomaly_from_true_anomaly(e, f).unwrap();
                    let m = calc_mean_anomaly_from_conic_anomaly(e, ecc).unwrap();
                    (q, e, f, psi, h, energy, m, r, vr, vo, orbit_type)
                }

                InputFormat::KEP => {
                    let a = p1;
                    let e = p2;
                    let f = p3;
                    let psi = p4;
                    let hsq = a * mu * (1.0 - e.powi(2));
                    let h = hsq.sqrt();
                    let r = hsq / mu * (1.0 / (1.0 + e * f.cos()));
                    let vsq = mu * (2.0 / r - 1.0 / a);
                    let vo = ((e.powi(2) - 1.0) / r.powi(2) * (vsq / mu.powi(2) - 2.0 / (mu * r)).powi(-1)).sqrt();
                    let vr_sq = vsq - vo.powi(2);
                    let vr = if vr_sq < 0.0 { 0.0 } else { vr_sq.sqrt() * if f > PI { -1.0 } else { 1.0 } };
                    let orbit_type = OrbitType::from_eccentricity(e, 1e-10).unwrap();
                    let m = calc_mean_anomaly_from_conic_anomaly(e, calc_conic_anomaly_from_true_anomaly(e, f).unwrap()).unwrap();
                    let q = match orbit_type {
                        OrbitType::Circular => a,
                        _ => a * (1.0 - e),
                    };
                    (q, e, f, psi, h, -0.5 * mu.powi(2) / hsq * (1.0 - e.powi(2)), m, r, vr, vo, orbit_type)
                }

                InputFormat::QEF => {
                    let q = p1;
                    let e = p2;
                    let f = p3;
                    let psi = p4;
                    let p = q * (1.0 + e);
                    let h = (mu * p).sqrt();
                    let r = p / (1.0 + e * f.cos());
                    let vo = h / r;
                    let vr = (mu / h) * e * f.sin();
                    let energy = -0.5 * mu.powi(2) / h.powi(2) * (1.0 - e.powi(2));
                    let m = calc_mean_anomaly_from_conic_anomaly(e, calc_conic_anomaly_from_true_anomaly(e, f).unwrap()).unwrap();
                    let orbit_type = OrbitType::from_eccentricity(e, 1e-10).unwrap();
                    (q, e, f, psi, h, energy, m, r, vr, vo, orbit_type)
                }
            };

        Self {
            id,
            format,
            p1_min,
            p1_max,
            p1,
            p2_min,
            p2_max,
            p2,
            p3_min,
            p3_max,
            p3,
            p4_min,
            p4_max,
            p4,
            epoch,
            mu,
            q,
            e,
            true_anomaly,
            psi,
            h,
            energy,
            mean_anomaly,
            r0,
            vr,
            vo,
            orbit_type,
        }
    }



    pub fn chi(&self, epoch: Time) -> Option<f64> {
        let dt = epoch.epoch - self.epoch.epoch;
        
        match solve_for_universal_anomaly(
            self.r0, 
            self.vr, 
            self.alpha(), 
            self.mu, 
            dt,
            1e-7, 
            100 
        ) {
            Ok(val) => Some(val),
            Err(_) => {
                println!(
                    "Failed to solve for universal_anomaly: r={}, vr={}, alpha={}, mu={}, dt={}",
                    self.r0, self.vr, self.alpha(), self.mu, dt
                );
                None
            }
        }
    }


    pub fn lagrange_f_and_g(&self, epoch: Time) -> Option<(f64, f64)> {
        let dt = epoch.epoch - self.epoch.epoch;
        let chi = self.chi(epoch)?;  // Propagate None with ?
        
        let c = stumpff_c(self.alpha() * chi.powi(2));
        let s = stumpff_s(self.alpha() * chi.powi(2));
        
        let f = 1.0 - chi.powi(2) * c / self.r0;
        let g = dt - chi.powi(3) * s / self.mu.sqrt();

        Some((f, g))
    }

    pub fn calculate_r(&self, epoch: Time) -> Option<f64> {
        let (f, g) = self.lagrange_f_and_g(epoch)?;

        let rsq = f.powi(2) * self.r0.powi(2) + 
                2.0 * f * g * self.r0 * self.vr + 
                g.powi(2) * (self.vo.powi(2) + self.vr.powi(2));
        
        Some(rsq.sqrt())
    }

    pub fn a(&self) -> Option<f64> {
        match self.orbit_type {
            OrbitType::Hyperbolic => Some(self.q / (self.e - 1.0)),
            OrbitType::Parabolic => None,
            OrbitType::Elliptical => Some(self.q / (1.0 - self.e)),
            OrbitType::Circular => Some(self.q),
            _ => None,
        }
    }

    pub fn n(&self) -> f64 {
        match self.a() {
            Some(a) => (self.mu / a.powi(3)).sqrt(),
            None => (self.mu / self.q.powi(3)).sqrt(),
        }
    }

    pub fn alpha(&self) -> f64 {
        -2.0 * self.energy / self.mu
    }


    pub fn separation(&self, other: &InitialCondition, epoch: Time) -> Option<f64> {
        let self_r = self.calculate_r(epoch.clone())?;
        let (self_f, self_g) = self.lagrange_f_and_g(epoch.clone())?;

        let other_r = other.calculate_r(epoch.clone())?;
        let (other_f, other_g) = other.lagrange_f_and_g(epoch.clone())?;

        let r_prod = self.r0 * other.r0;
        let r_new_prod = self_r * other_r;
        let psi_cos = (self.psi - other.psi).cos();
        
        let first_term = r_prod * self_f * other_f;
        let second_term = self.r0 * other.vr * self_f * other_g;
        let third_term = self.vr * other.r0 * self_g * other_f;
        let fourth_term = self_g * other_g * (self.vr * other.vr + self.vo * other.vo * psi_cos);

        let argument = (first_term + second_term + third_term + fourth_term) / r_new_prod;
        
        Some(nice_acos(argument))
    }
}

impl InitialCondition {
    pub fn to_bounds_key(&self) -> Option<IcBoundsKey> {
        let p1_bounds = (OrderedFloat(self.p1_min?), OrderedFloat(self.p1_max?));
        let p2_bounds = (OrderedFloat(self.p2_min?), OrderedFloat(self.p2_max?));
        let p3_bounds = (OrderedFloat(self.p3_min?), OrderedFloat(self.p3_max?));
        let p4_bounds = (OrderedFloat(self.p4_min?), OrderedFloat(self.p4_max?));
        let epoch = OrderedFloat(self.epoch.epoch);

        match self.format {
            InputFormat::QEF => Some(IcBoundsKey::QEF {
                q_bounds: p1_bounds,
                e_bounds: p2_bounds,
                f_bounds: p3_bounds,
                psi_bounds: p4_bounds,
                epoch,
            }),
            InputFormat::KEP => Some(IcBoundsKey::KEP {
                a_bounds: p1_bounds,
                e_bounds: p2_bounds,
                f_bounds: p3_bounds,
                psi_bounds: p4_bounds,
                epoch,
            }),
            InputFormat::KEV => Some(IcBoundsKey::KEV {
                r_bounds: p1_bounds,
                vr_bounds: p2_bounds,
                vo_bounds: p3_bounds,
                psi_bounds: p4_bounds,
                epoch,
            }),
        }
    }
}


impl std::fmt::Display for InitialCondition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // write!(f, "r: {}, vr: {}, vo: {}, psi: {}, orbit_type: {:?}", 
        //        self.r0, self.vr, self.vo, self.psi, self.orbit_type)
        write!(f, "q: {}, e: {}, psi: {}, true_anomaly: {}, epoch: {}, 
                    h: {}, energy: {}, mean_anomaly: {}, r0: {}, vr: {}, 
                    vo: {}, mu: {}, orbit_type: {:?}", self.q, self.e, self.psi,
                    self.true_anomaly, self.epoch.epoch, self.h, self.energy, self.mean_anomaly,
                    self.r0, self.vr, self.vo, self.mu, self.orbit_type)
    }
}