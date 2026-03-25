use spacerocks::time::Time;
use spacerocks::transforms::{solve_for_universal_anomaly, stumpff_s, stumpff_c, calc_conic_anomaly_from_true_anomaly, calc_mean_anomaly_from_conic_anomaly};
use spacerocks::OrbitType;
use std::f64::consts::PI;
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentOpt;

fn nice_acos(x: f64) -> f64 {
    if x > 1.0 {
        0.0
    } else if x < -1.0 {
        PI
    } else {
        x.acos()
    }
}


#[derive(Debug, Clone)]
pub struct Orbit {
    pub id: String,
    pub q: f64,
    pub e: f64,
    pub psi: f64,
    pub true_anomaly: f64,
    pub epoch: f64,

    pub h: f64,
    pub energy: f64,
    pub mean_anomaly: f64,

    pub r0: f64,
    pub vr: f64,
    pub vo: f64,

    pub mu: f64,
    pub orbit_type: OrbitType,
}

impl Orbit {
    pub fn new(id: String, q: f64, e: f64, true_anomaly: f64, psi: f64, epoch: f64, mu: f64) -> Self {
        let id = String::from("id");
        let p = q * (1.0 + e); 
        let h = (mu * p).sqrt();
        let energy = -0.5 *(mu.powi(2) / h.powi(2)) * (1.0 - e.powi(2));
        let r0 = p / (1.0 + e * true_anomaly.cos());
        let vo = h / r0;
        let vr = (mu / h) * e * true_anomaly.sin();

        let eccentric_anomaly = calc_conic_anomaly_from_true_anomaly(e, true_anomaly).expect("Failed to calculate eccentric anomaly");
        let mean_anomaly = calc_mean_anomaly_from_conic_anomaly(e, eccentric_anomaly).expect("Failed to calculate mean anomaly");


        let orbit_type = OrbitType::from_eccentricity(e, 1e-10).expect("Failed to determine orbit type");
        Self { id, q, e, true_anomaly, psi, epoch, h, energy, mean_anomaly, r0, vr, vo, mu, orbit_type }
    }

    pub fn chi(&self, epoch: f64) -> f64 {
        let dt = epoch - self.epoch;
        
        let s = solve_for_universal_anomaly(
            self.r0, 
            self.vr, 
            self.alpha(), 
            self.mu, 
            dt,
            1e-6, // Need to check on this
            1000 // Number of iterations 
        ).unwrap_or_else(|_| panic!(
            "Failed to solve for universal_anomaly within Orbit method\n\
             r: {}, vr: {}, vo: {}, alpha: {}, mu: {}, dt: {}, q: {}, e: {}, psi: {}, true_anomaly: {:?}", 
            self.r0, self.vr, self.vo, self.alpha(), self.mu, dt, self.q, self.e, self.psi, self.true_anomaly
        )); 

        s
    }

    pub fn lagrange_f_and_g(&self, epoch: f64) -> (f64, f64) {
    
        let dt = epoch - self.epoch;
        let chi = self.chi(epoch);
        
        // Calculate Stumpff functions
        let c = stumpff_c(self.alpha() * chi.powi(2));
        let s = stumpff_s(self.alpha() * chi.powi(2));
        
        let f = 1.0 - chi.powi(2) * c / self.r0;
        let g = dt - chi.powi(3) * s / self.mu.sqrt();

        (f, g)
    }

    pub fn calculate_r(&self, epoch: f64) -> f64 {
        
        // Get Universal Anomaly chi
        let chi = self.chi(epoch);
        
        // Get Lagrange coefficients
        let (f,g) = self.lagrange_f_and_g(epoch);


        // Calculate r**2 from Lagrange coefficients and v0 = vr^2 + vo^2
        let rsq = f.powi(2) * self.r0.powi(2) + 
                 2.0 * f * g * self.r0 * self.vr + 
                 g.powi(2) * (self.vo.powi(2) + self.vr.powi(2));
        
        rsq.sqrt()
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


    pub fn separation(&self, other: &Orbit, epoch: f64) -> f64 {

        let self_r = self.calculate_r(epoch);
        let (self_f, self_g) = self.lagrange_f_and_g(epoch);

        let other_r = other.calculate_r(epoch);
        let (other_f, other_g) = other.lagrange_f_and_g(epoch);

        // Calculate common products
        let r_prod = self.r0 * other.r0;
        let r_new_prod = self_r * other_r;
        let psi_cos = (self.psi - other.psi).cos();
        
        // Calculate all terms maintaining original logic
        let first_term = r_prod * self_f * other_f;
        let second_term = self.r0 * other.vr * self_f * other_g;
        let third_term = self.vr * other.r0 * self_g * other_f;
        let fourth_term = self_g * other_g * (self.vr * other.vr + self.vo * other.vo * psi_cos);

        // Combine and normalize
        let argument = (first_term + second_term + third_term + fourth_term) / r_new_prod;
        
        nice_acos(argument)
    }
}

impl std::fmt::Display for Orbit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // write!(f, "r: {}, vr: {}, vo: {}, psi: {}, orbit_type: {:?}", 
        //        self.r0, self.vr, self.vo, self.psi, self.orbit_type)
        write!(f, "q: {}, e: {}, psi: {}, true_anomaly: {}, epoch: {}, 
                    h: {}, energy: {}, mean_anomaly: {}, r0: {}, vr: {}, 
                    vo: {}, mu: {}, orbit_type: {:?}", self.q, self.e, self.psi,
                    self.true_anomaly, self.epoch, self.h, self.energy, self.mean_anomaly,
                    self.r0, self.vr, self.vo, self.mu, self.orbit_type)
    }
}

// Separation.rs module in my grid making architecture
struct OrbitPair<'a> {
    orbit1: &'a Orbit,
    orbit2: &'a Orbit,
}

impl<'a> CostFunction for OrbitPair<'a> {
    type Param = f64;
    type Output = f64;

    fn cost(&self, &t: &Self::Param) -> Result<Self::Output, Error> {
        Ok(-self.orbit1.separation(self.orbit2, t))
    }
}

pub fn find_max_separation(orbit1: &Orbit, orbit2: &Orbit, bounds: (f64, f64)) -> f64 {
    let (tmin, tmax) = bounds;
    
    let cost = OrbitPair { orbit1, orbit2 };
    let solver = BrentOpt::new(tmin, tmax);

    let res = Executor::new(cost, solver)
        .run()
        .unwrap();

    -res.state.best_cost
}

pub fn find_max_separation_cheat(orbit1: &Orbit, orbit2: &Orbit, bounds: (f64, f64)) -> f64 {
    let (tmin, tmax) = bounds;
    
    // Calculate separation at both boundary points
    let sep_at_tmin = orbit1.separation(orbit2, tmin);
    let sep_at_tmax = orbit1.separation(orbit2, tmax);
    
    // Return the maximum of the two separations
    sep_at_tmin.max(sep_at_tmax)
}