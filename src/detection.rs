use spacerocks::{Time, Observer};
use nalgebra::{Vector3, DMatrix, Matrix3};

// equatorial_to_ecliptic_matrix and ecliptic_to_equatorial_matrix can be used for conversions
pub const EQUATORIAL_TO_ECLIPTIC: Matrix3<f64> = Matrix3::new(1.0, 0.0, 0.0,
                                                              0.0, 0.917_482_062_069_181_8, 0.397_777_155_931_913_7,
                                                              0.0, -0.397_777_155_931_913_7, 0.917_482_062_069_181_8);


pub const ECLIPTIC_TO_EQUATORIAL: Matrix3<f64> = Matrix3::new(1.0, 0.0, 0.0,
                                                              0.0, 0.917_482_062_069_181_8, -0.397_777_155_931_913_7,
                                                              0.0, 0.397_777_155_931_913_7, 0.917_482_062_069_181_8);
#[derive(Debug, Clone, PartialEq)]
pub enum ReferencePlane {
    Ecliptic,
    Equatorial,
}

// make the str for, so you can pass it to functions
impl ReferencePlane {
    pub fn as_str(&self) -> &str {
        match self {
            ReferencePlane::Ecliptic => "ECLIPJ2000",
            ReferencePlane::Equatorial => "J2000",
        }
    }
}

#[derive(Debug, Clone)]
pub struct Detection {

    pub rho_hat: Vector3<f64>,
    pub observer_position: Vector3<f64>,
    pub epoch: f64,
    pub reference_plane: ReferencePlane,

    pub detid: Option<String>,
    pub trackid: Option<String>,    
    pub observer_velocity: Option<Vector3<f64>>,
    
    pub mag: Option<f64>,
    pub mag_ucty: Option<f64>,
    pub filter: Option<String>,

    pub ast_ucty: Option<f64>,
    pub ra_ucty: Option<f64>,
    pub dec_ucty: Option<f64>,
    
    pub objid: Option<String>,
    pub obscode: Option<String>,
    
    pub rho_hat_dot_observer_position: f64,
    pub observer_distance_squared: f64,

    // New additions 
    pub expnum: Option<i64>,
    pub nite: Option<i64>,
    pub observer: Option<Observer>,
    pub ra: Option<f64>,
    pub dec: Option<f64>,

}

impl Detection {
    pub fn new(
        ra: f64,
        dec: f64,
        epoch: Time,
        observer_position: Vector3<f64>,
    ) -> Self {
        
        let rho_hat = Vector3::new(ra.cos() * dec.cos(), ra.sin() * dec.cos(), dec.sin());
        let epoch_jd = epoch.tdb().jd();
        let rho_hat_dot_observer_position = rho_hat.dot(&observer_position);
        let observer_distance_squared = observer_position.dot(&observer_position);

        Detection {
            rho_hat: rho_hat,
            observer_position: observer_position,
            epoch: epoch_jd,
            reference_plane: ReferencePlane::Equatorial,
            ast_ucty: None,
            ra_ucty: None,
            dec_ucty: None,
            observer_velocity: None,
            detid: None,
            objid: None,
            trackid: None,
            mag: None,
            mag_ucty: None,
            filter: None,
            obscode: None,
            rho_hat_dot_observer_position: rho_hat_dot_observer_position,
            observer_distance_squared: observer_distance_squared,
            expnum: None,
            nite: None,
            observer: None,
            ra: None, 
            dec: None,
        }
    }

    pub fn to_equatorial(&mut self) {

        if self.reference_plane == ReferencePlane::Equatorial {
            return;
        }
        
        self.rho_hat = ECLIPTIC_TO_EQUATORIAL * self.rho_hat;
        self.observer_position = ECLIPTIC_TO_EQUATORIAL * self.observer_position;
        if let Some(velocity) = &self.observer_velocity {
            self.observer_velocity = Some(ECLIPTIC_TO_EQUATORIAL * velocity);
        }
       self.reference_plane = ReferencePlane::Equatorial;
    }

    pub fn to_ecliptic(&mut self) {

        if self.reference_plane == ReferencePlane::Ecliptic {
            return;
        }
        // rotate rho_hat and observer_position accordingly
        self.rho_hat = EQUATORIAL_TO_ECLIPTIC * self.rho_hat;
        self.observer_position = EQUATORIAL_TO_ECLIPTIC * self.observer_position;
        if let Some(velocity) = &self.observer_velocity {
            self.observer_velocity = Some(EQUATORIAL_TO_ECLIPTIC * velocity);
        }
       self.reference_plane = ReferencePlane::Ecliptic;
    }

    pub fn pointing(&self) -> Vector3<f64> {
        let x = self.ra.unwrap_or(0.0).cos() * self.dec.unwrap_or(0.0).cos();
        let y = self.ra.unwrap_or(0.0).sin() * self.dec.unwrap_or(0.0).cos();
        let z = self.dec.unwrap_or(0.0).sin();
        Vector3::new(x, y, z)
    }

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

    pub fn from_observer(
        ra: f64,
        dec: f64,
        epoch: Time,
        observer: Observer,
    ) -> Self {
        let mut det = Self::new(ra, dec, epoch, observer.position);
        det.observer = Some(observer);
        det
    }







    pub fn set_observer_position(&mut self, position: Vector3<f64>) {
        self.observer_position = position;
    }

    // make setters
    pub fn set_observer_velocity(&mut self, velocity: Vector3<f64>) {
        self.observer_velocity = Some(velocity);
    }

    pub fn set_detid(&mut self, detid: String) {
        self.detid = Some(detid);
    }

    pub fn set_trackid(&mut self, trackid: String) {
        self.trackid = Some(trackid);
    }

    pub fn set_objid(&mut self, objid: String) {
        self.objid = Some(objid);
    }

    pub fn set_mag(&mut self, mag: f64) {
        self.mag = Some(mag);
    }

    pub fn set_ast_ucty(&mut self, ast_ucty: f64) {
        self.ast_ucty = Some(ast_ucty);
    }

    pub fn set_ra_ucty(&mut self, ra_ucty: f64) {
        self.ra_ucty = Some(ra_ucty);
    }

    pub fn set_dec_ucty(&mut self, dec_ucty: f64) {
        self.dec_ucty = Some(dec_ucty);
    }

    pub fn set_filter(&mut self, filter: String) {
        self.filter = Some(filter);
    }

    pub fn set_ra(&mut self, ra: f64) {
        self.ra = Some(ra);
    }

    pub fn set_dec(&mut self, dec: f64) {
        self.dec = Some(dec);
    }

    pub fn set_expnum(&mut self, expnum: i64) {
        self.expnum = Some(expnum);
    }

    pub fn set_nite(&mut self, nite: i64) {
        self.nite = Some(nite);
    }

    pub fn set_observer(&mut self, observer: Observer) {
        self.observer = Some(observer);
    }



    // make getters
    pub fn rho_hat(&self) -> &Vector3<f64> {
        &self.rho_hat
    }
    pub fn observer_position(&self) -> &Vector3<f64> {
        &self.observer_position
    }
    pub fn epoch(&self) -> f64 {
        self.epoch
    }
    pub fn reference_plane(&self) -> &ReferencePlane {
        &self.reference_plane
    }
    pub fn detid(&self) -> &Option<String> {
        &self.detid
    }
    pub fn trackid(&self) -> &Option<String> {
        &self.trackid
    }
    pub fn observer_velocity(&self) -> &Option<Vector3<f64>> {
        &self.observer_velocity
    }
    pub fn mag(&self) -> &Option<f64> {
        &self.mag
    }
    pub fn ast_ucty(&self) -> &Option<f64> {
        &self.ast_ucty
    }
    pub fn filter(&self) -> &Option<String> {
        &self.filter
    }
    pub fn ra_ucty(&self) -> &Option<f64> {
        &self.ra_ucty
    }
    pub fn dec_ucty(&self) -> &Option<f64> {
        &self.dec_ucty
    }

    pub fn ra(&self) -> Option<f64> {
        self.ra
    }

    pub fn dec(&self) -> Option<f64> {
        self.dec
    }

    pub fn expnum(&self) -> Option<i64> {
        self.expnum
    }

    pub fn nite(&self) -> Option<i64> {
        self.nite
    }

    pub fn observer(&self) -> Option<&Observer> {
        self.observer.as_ref()
    }
}