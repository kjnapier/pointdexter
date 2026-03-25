use spacerocks::time::Time;
use spacerocks::observing::{Observer, Observation};

use nalgebra::Vector3;

#[derive(Debug, Clone, PartialEq)]
pub struct Detection {
    pub objid: String,
    pub fakeid: Option<String>,
    pub intid: usize,
    // pub linkid: Option<String>,
    pub ra: f64,
    pub dec: f64,
    pub ra_ucty: f64, 
    pub dec_ucty: f64,
    pub pointing: Vector3<f64>,
    pub epoch: Time,
    pub expnum: i64,
    pub nite: i64,
    pub observer: Observer,
    pub observer_dot_pointing: f64,
    pub observer_distance: f64,
    // pub cos_solar_elongation: f64,
    // pub sin_solar_elongation_sq: f64,
    pub mag: Option<String>,
    pub orbit_id: Option<i64>,

    // New elements with tracklets. Removed some elements above that aren't needed right now
    pub ra_rate: Option<f64>,
    pub dec_rate: Option<f64>,
    pub num_dets: i64,
    pub original_intids: Option<Vec<usize>>,
}

impl Detection {
    pub fn new(
        objid: String, 
        fakeid: Option<String>, 
        intid: usize, 
        // linkid: Option<String>, 
        ra: f64, 
        dec: f64, 
        ra_ucty: f64, 
        dec_ucty: f64, 
        epoch: Time, 
        expnum: i64, 
        //nite: i64,
        observer: Observer, 
        mag: Option<String>,
        orbit_id: Option<i64>, 
        ra_rate: Option<f64>, 
        dec_rate: Option<f64>, 
        num_dets: i64, 
        original_intids: Option<Vec<usize>>) -> Self {

        let pointing = Vector3::new(ra.cos() * dec.cos(), ra.sin() * dec.cos(), dec.sin());
        let observer_dot_pointing = observer.position.dot(&pointing);
        let observer_distance = observer.position.norm();
        let o = observer.position / observer_distance;
        let cos_solar_elongation = -o.dot(&pointing);
        let sin_solar_elongation_sq = 1.0 - cos_solar_elongation.powi(2);

        // Sketchy way to get nite from epoch 
        let nite = epoch.epoch.floor() as i64;  

        Self {
            objid,
            fakeid,
            intid,
            // linkid,
            ra,
            dec,
            ra_ucty,
            dec_ucty,
            pointing: pointing,
            epoch,
            expnum,
            nite,
            observer,
            observer_dot_pointing,
            observer_distance,
            // cos_solar_elongation,
            // sin_solar_elongation_sq,
            mag: mag,
            orbit_id: orbit_id,
            ra_rate,
            dec_rate,
            num_dets,
            original_intids,

        }
    }

    pub fn pointing(&self) -> Vector3<f64> {
        let x = self.ra.cos() * self.dec.cos();
        let y = self.ra.sin() * self.dec.cos();
        let z = self.dec.sin();
        Vector3::new(x, y, z)
    }

    // Right now, this only works for astrometric (ra, dec) observations
    pub fn to_observation(&self) -> Observation {
        let covariance = [[self.ra_ucty.powi(2), 0.0], [0.0, self.dec_ucty.powi(2)]];
        let mag: Option<f64> = self.mag.clone().and_then(|s| s.parse().ok());
        Observation::from_astrometry(
            self.epoch.clone(),
            self.ra.clone(),
            self.dec.clone(),
            self.observer.clone(),
            Some(covariance),
            mag,
            None,
        ).expect("Failed to create observation from detection")
    }
}
