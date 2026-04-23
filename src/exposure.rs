use crate::ReferencePlane;

use nalgebra::{Vector3, DMatrix, Matrix3};

use crate::xyz_to_proj_matrix;

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

    //pub xyz_e: Vector3<f64>,
    
    pub theta_x0: f64,
    pub theta_y0: f64,

    pub theta_x: Vec<f64>,
    pub theta_y: Vec<f64>,
}


pub struct Exposure {
    pub id: String,
    pub epoch: f64,
    pub filter: Option<String>,
    pub central_rho: Vector3<f64>,

    pub detections: Vec<Vector3<f64>>,
    pub observer_position: Vector3<f64>,
    pub observer_velocity: Option<Vector3<f64>>,

    pub reference_plane: ReferencePlane,
}


impl Exposure {

    pub fn transform_to_tangent_plane(&self, ref_vec: Vector3<f64>) -> TangentPlaneExposure {

        let mut rot = xyz_to_proj_matrix(ref_vec);

        let mut theta_x = Vec::new();
        let mut theta_y = Vec::new();

        for det in self.detections.iter() {
            let proj = rot * det;

            let tx = proj[0]/proj[2];
            let ty = proj[1]/proj[2];
            
            theta_x.push(tx);
            theta_y.push(ty);

        }

        // Calculate the central theta_x0 and theta_y0 for this exposure as the projection of the central rho vector.
        let central_proj = rot * self.central_rho;
        let theta_x0 = central_proj[0]/central_proj[2];
        let theta_y0 = central_proj[1]/central_proj[2];
        

        let xyz_e = rot * self.observer_position;
        let xe = xyz_e[0];
        let ye = xyz_e[1];
        let ze = xyz_e[2];

        TangentPlaneExposure {
            id: self.id.clone(),
            epoch: self.epoch,
            xe,
            ye,
            ze,
            theta_x0,
            theta_y0,
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
