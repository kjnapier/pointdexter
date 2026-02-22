

// load spice kernels
// let mut kernel = SpiceKernel::new();
// kernel.load_spk(format!("{}/sb441-n16.bsp", spice_root).as_str())?;
// kernel.load_spk(format!("{}/de440s.bsp", spice_root).as_str())?;
// kernel.load_bpc(format!("{}/earth_1962_240827_2124_combined.bpc", spice_root).as_str())?;

// let mut observatory = Observatory::from_obscode(obscode)?;

// read in the detections from a file
// let file_path = "/Users/kjnapier/Desktop/link/ps1-catalog-culled.csv";
// let mut rdr = csv::Reader::from_path(file_path)?;

// // Read each row
// for result in rdr.records() {
//     let record = result?;
//     let detID: String = record[2].parse()?;
//     let ra: f64 = record[3].parse()?;
//     let dec: f64 = record[4].parse()?;
//     let astrom_sig: f64 = record[5].parse()?;
//     let epoch: f64 = record[12].parse()?;
//     let intid: i32 = record[14].parse()?;

//     let o = observatory.at(&Time::new(epoch, "utc", "jd")?, "J2000", "ssb", &kernel)?;
//     let observer_position = o.position;
//     let observer_velocity = o.velocity.unwrap();

//     let det = Detection::new(
//         ra * std::f64::consts::PI / 180.0,
//         dec * std::f64::consts::PI / 180.0,
//         astrom_sig * (1. / 3600.) * std::f64::consts::PI / 180.0,
//         Time::new(epoch, "utc", "jd")?,
//         observer_position,
//         observer_velocity,
//         detID,
//     );
//     detections.push(det);
// }

// Read the initial conditions from a file
// let ic_file_path = "../ics_r>200_4year_30arcsec.csv";
// let mut ic_rdr = csv::Reader::from_path(ic_file_path)?;

// let mut initial_conditions: Vec<InitialCondition> = Vec::new();
// let mut count = 0;

// let epoch = 2456569.5;

// println!("Reading initial conditions from {}", ic_file_path);

// for (j, result) in ic_rdr.records().enumerate() {
//     let record = result?;
//     let r: f64 = record[0].parse()?;
//     let vr: f64 = record[1].parse()?;
//     let vo: f64 = record[2].parse()?;
//     let inc: f64 = record[3].parse()?;
//     let kappa_f64: f64 = record[4].parse()?;

//     let kappa = kappa_f64 as i32;

//     let orbID = format!("ic_{:06}", j);
//     let ep = Time::new(epoch, "tdb", "jd")?;

//     let mu = ssb.mu();

//     let ic = InitialCondition::from_spherical(orbID, r, vr, vo, inc, kappa, ep, mu);

//     initial_conditions.push(ic?);
//     count += 1;
// }


// initial_conditions[..50].par_iter().for_each(|ic| {
    //     let mut points: Vec<[f64; 3]> = Vec::with_capacity(detections.len());
    //     for detection in &detections[..] {
    //         if let Some(pointing) = sync_detection_to_orbit(&detection, &ic) {
    //             points.push(pointing);
    //         }
    //     }
    // });




    // let filtered_dts: Vec<f64> = dts.iter().cloned().filter(|&dt| dt <= max_dt).collect();
            // if filtered_dts.len() < config.min_detections {
            //     continue;
            // }

            // let mut kept_indices: Vec<usize> = indices.iter()
            //     .zip(dts.iter())
            //     .filter(|&(_i, &dt)| dt <= max_dt)
            //     .map(|(&i, &_dt)| i)
            //     .collect();





             // let mut thresh = ic.inc;
    // if ic.inc > std::f64::consts::PI / 2.0 {
    //     thresh = std::f64::consts::PI - ic.inc;
    // }
    // if (r_vec[2] / r).asin().abs() > thresh {
    //     return None;
    // }

    // for exposure in tangent_exposures.iter() {
    //     let rho2 = (1.0 + mean_theta_x * mean_theta_x + mean_theta_y * mean_theta_y)*(z0 - exposure.xyz_e.z).powi(2);
    //     let rho = rho2.sqrt();

    //     // Improve light time correction at some point.
    //     let dt = rho/SPEED_OF_LIGHT;

    //     let t = exposure.epoch - ref_epoch;
    //     let tp = t - dt;
    //     let f = (1.0 - 0.5 * mm_sqr * tp * tp);
    //     let g = tp;
    //     let fac = (1.0 + g/f * gdot - gamma/f * exposure.xyz_e.z);
    //     for (theta_x, theta_y) in exposure.theta_x.iter().zip(exposure.theta_y.iter()) {
    //         let phi_x = theta_x * fac + gamma/f * exposure.xyz_e.x;
    //         let phi_y = theta_y * fac + gamma/f * exposure.xyz_e.y;
    //         let alpha = phi_x - g/f * adot;
    //         let beta = phi_y - g/f * bdot;
    //         //println!("Exposure {}: alpha = {}, beta = {}", exposure.id, alpha, beta);
    //         alphas.push(alpha);
    //         betas.push(beta);
    //     }
    // }



    // // We can set the maximum value of adot and bdot based on the value of gamma.
            
    // let theta_dot_max = (2.0 * mu * gamma * gamma * gamma).sqrt(); // max angular rate in radians/day
    // // print the max angular rate in arcsec/hour for readability
    // println!("theta_dot_max: {:.3} arcsec/hour", theta_dot_max * 3600.0 * 180.0 / std::f64::consts::PI / 24.0);
    // println!("approx period: {:.2} years", (2.0 * std::f64::consts::PI / theta_dot_max) / 365.25);

    // // Next, we set the step size for our grid in adot and bdot based on the time span
    // // of our observations and the search radius.
    // // The units should be radians/day for both adot and bdot, since they represent angular velocities in the tangent plane.    
    // let dtheta_dot = 3.0 * pixel_scale / time_span; // step size in radians/day, chosen so that the object moves by 3 pixels over the time span of the observations
    // println!("{}", theta_dot_max/dtheta_dot);

    // // Now we can calculate the number of steps in adot and bdot directions.

    // let adot_steps = 2 * (theta_dot_max / dtheta_dot).ceil() as usize;
    // let bdot_steps = 2 * (theta_dot_max / dtheta_dot).ceil() as usize;

    // println!("theta_dot_max: {:.3} arcsec/day, dtheta_dot: {:.3} arcsec/day, adot_steps: {}, bdot_steps: {}", 
    //     theta_dot_max * 3600.0 * 180.0 / std::f64::consts::PI, 
    //     dtheta_dot * 3600.0 * 180.0 / std::f64::consts::PI, 
    //     adot_steps, 
    //     bdot_steps
    // );

    // // Finally, we can calculate the grid in adot and bdot.
    // // We will center the grid around 0, so we will have adot and bdot values ranging from -theta_dot_max to +theta_dot_max. 
    // // We will exclude points such that sqrt(adot^2 + bdot^2) > theta_dot_max to avoid unphysical orbits.
    // let mut adot_bdot_grid: Vec<(f64, f64)> = Vec::new();
    // for i in 0..(adot_steps+1) {
    //     let adot = -theta_dot_max + i as f64 * dtheta_dot;
    //     for j in 0..(bdot_steps+1) {
    //         let bdot = -theta_dot_max + j as f64 * dtheta_dot;
    //         if (adot*adot + bdot*bdot).sqrt() <= theta_dot_max {
    //             adot_bdot_grid.push((adot, bdot));
    //         }
    //     }
    // }
    // println!("Grid size: {} adot steps x {} bdot steps = {} total points", adot_steps, bdot_steps, adot_bdot_grid.len());
    
    