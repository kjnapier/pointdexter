use crate::ReferencePlane;

use nalgebra::{Vector3, DMatrix, Matrix3};

// equatorial_to_ecliptic_matrix and ecliptic_to_equatorial_matrix can be used for conversions
pub const EQUATORIAL_TO_ECLIPTIC: Matrix3<f64> = Matrix3::new(1.0, 0.0, 0.0,
                                                              0.0, 0.917_482_062_069_181_8, 0.397_777_155_931_913_7,
                                                              0.0, -0.397_777_155_931_913_7, 0.917_482_062_069_181_8);


pub const ECLIPTIC_TO_EQUATORIAL: Matrix3<f64> = Matrix3::new(1.0, 0.0, 0.0,
                                                              0.0, 0.917_482_062_069_181_8, -0.397_777_155_931_913_7,
                                                              0.0, 0.397_777_155_931_913_7, 0.917_482_062_069_181_8);



pub struct TangentPlaneExposure {
    pub id: String,
    pub epoch: f64,

    pub xe: f64,
    pub ye: f64,
    pub ze: f64,

    pub theta_x: Vec<f64>,
    pub theta_y: Vec<f64>,
}


pub struct Exposure {
    pub id: String,
    pub epoch: f64,
    pub filter: Option<String>,

    pub detections: Vec<Vector3<f64>>,
    pub observer_position: Vector3<f64>,
    pub observer_velocity: Option<Vector3<f64>>,

    pub reference_plane: ReferencePlane,
}


impl Exposure {

    pub fn transform_to_tangent_plane(&self, center: Vector3<f64>) -> TangentPlaneExposure {
        // Implement logic here.
        let theta_x = self.detections.iter().map(|d| (d - center).dot(&Vector3::x())).collect();
        let theta_y = self.detections.iter().map(|d| (d - center).dot(&Vector3::y())).collect();
        TangentPlaneExposure {
            id: self.id.clone(),
            epoch: self.epoch,
            xe: self.observer_position.x,  // TODO: THIS IS WRONG
            ye: self.observer_position.y, // TODO: THIS IS WRONG
            ze: self.observer_position.z, // TODO: THIS IS WRONG
            theta_x,
            theta_y,
        }
    }
    
    pub fn to_equatorial(&mut self) {

        if self.reference_plane == ReferencePlane::Equatorial {
            return;
        }

        for detection in self.detections.iter_mut() {
            *detection = ECLIPTIC_TO_EQUATORIAL * *detection;
        }
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
        for detection in self.detections.iter_mut() {
            *detection = EQUATORIAL_TO_ECLIPTIC * *detection;
        }
        self.observer_position = EQUATORIAL_TO_ECLIPTIC * self.observer_position;
        if let Some(velocity) = &self.observer_velocity {
            self.observer_velocity = Some(EQUATORIAL_TO_ECLIPTIC * velocity);
        }
       self.reference_plane = ReferencePlane::Ecliptic;
    }
}
