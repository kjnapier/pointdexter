

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