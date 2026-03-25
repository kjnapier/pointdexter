use crate::sparse_linking_stuff::Detection;
use crate::sparse_linking_stuff::InitialCondition;
use crate::sparse_linking_stuff::InputFormat;
use crate::sparse_linking_stuff::utils::ARCSEC_PER_RAD;

use spacerocks::time::Time;
use spacerocks::observing::Observatory;
use spacerocks::constants::MU_BARY;
use spacerocks::SpaceRock;
use spacerocks::transforms::calc_true_anomaly_from_mean_anomaly;
use spacerocks::coordinates::Origin;
use spacerocks::SpiceKernel;

use indicatif::ProgressIterator;
use polars::prelude::*;
use std::collections::{HashMap, HashSet};

use std::fs::File;
use csv::ReaderBuilder;

use std::sync::Arc;

pub fn catalog_to_detections(path: &str, kernel: &SpiceKernel) -> Result<Vec<Detection>, Box<dyn std::error::Error>> {
    let mut catalog = CsvReadOptions::default()
        .with_infer_schema_length(Some(100))
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(path.into()))?
        .finish()?;

    // Check which tracklet columns exist
    let schema = catalog.schema();
    let has_ra_rate = schema.get("ra_rate").is_some();
    let has_dec_rate = schema.get("dec_rate").is_some();
    let has_num_dets = schema.get("num_detections").is_some();
    let has_original_intids = schema.get("original_intids").is_some();

    let ras: Vec<f64> = catalog["ra"].f64()?.into_no_null_iter().collect();
    let decs: Vec<f64> = catalog["dec"].f64()?.into_no_null_iter().collect();

    // Convert ra and dec to radians
    let ras: Vec<f64> = ras.iter().map(|x| x * std::f64::consts::PI / 180.0).collect();
    let decs: Vec<f64> = decs.iter().map(|x| x * std::f64::consts::PI / 180.0).collect();

    let epochs: Vec<f64> = catalog["epoch"].f64()?.into_no_null_iter().collect();
    let obscodes: Vec<&str> = catalog["obscode"].str()?.into_no_null_iter().collect();
    let expnums: Vec<i64> = catalog["expnum"].i64()?.into_no_null_iter().collect();
    let objids: Vec<&str> = catalog["objid"].str()?.into_no_null_iter().collect();
    let fakeids: Vec<Option<&str>> = catalog["fakeid"].str()?.into_iter().collect();
    let magnitudes: Vec<Option<f64>> = catalog["mag"].f64()?.into_iter().collect();
    let orbit_ids: Vec<Option<i64>> = catalog["orbit_id"]
        .cast(&DataType::Int64)?
        .i64()?
        .into_iter()
        .collect();

    // Conditionally read tracklet fields based on column existence
    let ra_rates: Vec<Option<f64>> = if has_ra_rate {
        catalog["ra_rate"].f64()?.into_iter().collect()
    } else {
        vec![None; catalog.height()]
    };

    let dec_rates: Vec<Option<f64>> = if has_dec_rate {
        catalog["dec_rate"].f64()?.into_iter().collect()
    } else {
        vec![None; catalog.height()]
    };

    let num_dets: Vec<i64> = if has_num_dets {
        catalog["num_detections"].i64()?.into_no_null_iter().collect()
    } else {
        vec![1; catalog.height()]
    };

    let original_objids_raw: Vec<&str> = if has_original_intids {
        catalog["original_intids"].str()?.into_no_null_iter().collect()
    } else {
        vec![""; catalog.height()]
    };

    // Make a hash table of observatories
    let mut observatories = HashMap::new();
    for obscode in obscodes.iter().collect::<HashSet<_>>() {
        let observatory = Observatory::from_obscode(obscode)?;
        observatories.insert(obscode, observatory);
    }

    // Make a hash table of observer positions at the epoch of each expnum
    let mut observer_positions = HashMap::new();
    println!("Computing observer positions for each detection.");
    for (idx, expnum) in expnums.iter().progress().enumerate() {
        let epoch = Time::new(epochs[idx], "utc", "jd")
            .map_err(|e| Box::<dyn std::error::Error>::from(format!("Error constructing Time for epoch {}: {}", epochs[idx], e)))?;
        let observatory = observatories.get(&obscodes[idx]).expect("Observatory not found in hash table");
        let mut observer = observatory.at(&epoch.clone(), "J2000", "SSB", kernel)
            .expect("Failed to create observer");
        observer_positions.insert(expnum, observer);
    }

    // Make a vector of integer ids from 0 to catalog.height()
    let integer_ids: Vec<usize> = (0..catalog.height() as usize).collect();

    // make a vector of Detection objects
    let mut detections = Vec::new();

    for idx in 0..catalog.height() {
        let ra = ras[idx];
        let dec = decs[idx];
        let epoch = Time::new(epochs[idx], "utc", "jd")
            .map_err(|e| Box::<dyn std::error::Error>::from(format!("Error constructing Time for epoch {}: {}", epochs[idx], e)))?;
        let expnum = expnums[idx];
        let observer = observer_positions.get(&expnum)
            .expect("Observer position not found for expnum")
            .clone();
        let intid = integer_ids[idx];
        let fakeid = match fakeids[idx] {
            Some(fakeid) => Some(fakeid.to_string()),
            None => None,
        };
        let mag = match magnitudes[idx] {
            Some(mag) => Some(mag.to_string()),
            None => None,
        };

        let orbit_id = match orbit_ids[idx] {
            Some(orbit_id) => Some(orbit_id),
            None => None,
        };

        let ra_rate = ra_rates[idx];
        let dec_rate = dec_rates[idx];
        let num_dets_val = num_dets[idx];

        // Parse pipe-separated original_objids into Vec<usize>
        let original_objids: Option<Vec<usize>> = if has_original_intids && !original_objids_raw[idx].is_empty() {
            let parsed_ids: Vec<usize> = original_objids_raw[idx]
                .split('|')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse::<usize>().ok())
                .collect();
            if parsed_ids.is_empty() {
                None
            } else {
                Some(parsed_ids)
            }
        } else {
            None
        };

        let detection = Detection::new(
            objids[idx].to_string(),
            fakeid,
            intid,
            ra,
            dec,
            0.15 / ARCSEC_PER_RAD,
            0.15 / ARCSEC_PER_RAD,
            epoch,
            expnum,
            observer,
            mag,
            orbit_id,
            ra_rate,
            dec_rate,
            num_dets_val,
            original_objids,
        );
        detections.push(detection);
    }

    Ok(detections)
}


pub fn catalog_to_fakes(path: &str) -> Result<Vec<SpaceRock>, Box<dyn std::error::Error>> {
    let mut catalog = CsvReadOptions::default()
        .with_infer_schema_length(Some(100))
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(path.into()))?
        .finish()?;

    let fakeids: Vec<&str> = catalog["FAKEID"].str()?.into_no_null_iter().collect();
    let a: Vec<f64> = catalog["a"].f64()?.into_no_null_iter().collect();
    let e: Vec<f64> = catalog["e"].f64()?.into_no_null_iter().collect();
    let inc: Vec<f64> = catalog["inc"].f64()?.into_no_null_iter().collect();
    let arg: Vec<f64> = catalog["omega"].f64()?.into_no_null_iter().collect();
    let node: Vec<f64> = catalog["Omega"].f64()?.into_no_null_iter().collect();
    let M: Vec<f64> = catalog["M"].f64()?.into_no_null_iter().collect();
    let epoch: Vec<f64> = catalog["EPOCH"].f64()?.into_no_null_iter().collect();
    let H: Vec<f64> = catalog["H"].f64()?.into_no_null_iter().collect();

    let mut fakes = Vec::new();

    for idx in 0..catalog.height() {
        let fakeid = fakeids[idx];
        let a = a[idx];
        let e = e[idx];
        let inc = inc[idx] * std::f64::consts::PI / 180.0;
        let arg = arg[idx] * std::f64::consts::PI / 180.0;
        let node = node[idx] * std::f64::consts::PI / 180.0;
        let M = M[idx] * std::f64::consts::PI / 180.0;
        let epoch = Time::new(epoch[idx], "utc", "jd").expect("Error in constructing Time");
        let H = H[idx];

        let q = a * (1.0 - e);
        let f = calc_true_anomaly_from_mean_anomaly(e, M).expect("Error in calc_true_anomaly_from_mean_anomaly");

        let mut fake = SpaceRock::from_kepler(fakeid, q, e, inc, arg, node, f, epoch, "ECLIPJ2000", "SSB")?;
        fake.set_absolute_magnitude(H);
        fakes.push(fake);
    }
    Ok(fakes)
}


pub fn read_initial_conditions_bounded(
    path: &str,
    mu: f64,
) -> Result<Vec<InitialCondition>, Box<dyn std::error::Error>> {
    let df = CsvReadOptions::default()
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(path.into()))?
        .finish()?;

    // Determine input format by first column
    let headers = df.get_column_names();
    let format = match headers[0].as_str() {
        "q_min" => InputFormat::QEF,
        "a_min" => InputFormat::KEP,
        "r_min" => InputFormat::KEV,
        _ => return Err(format!("Unrecognized input format from column '{}'", headers[0]).into()),
    };

    // Map format to field name prefixes
    let (p1_name, p2_name, p3_name, p4_name) = match format {
        InputFormat::QEF => ("q", "e", "f", "psi"),
        InputFormat::KEP => ("a", "e", "f", "psi"),
        InputFormat::KEV => ("r", "vr", "vo", "psi"),
    };

    let height = df.height();
    let mut ics = Vec::with_capacity(height);

    for idx in 0..height {
        let id = idx;

        let p1_min = df.column(&format!("{}_min", p1_name))?.f64()?.get(idx);
        let p1_max = df.column(&format!("{}_max", p1_name))?.f64()?.get(idx);
        let p1 = df[p1_name].f64()?.get(idx).unwrap();

        let p2_min = df.column(&format!("{}_min", p2_name))?.f64()?.get(idx);
        let p2_max = df.column(&format!("{}_max", p2_name))?.f64()?.get(idx);
        let p2 = df[p2_name].f64()?.get(idx).unwrap();

        let p3_min = df.column(&format!("{}_min", p3_name))?.f64()?.get(idx);
        let p3_max = df.column(&format!("{}_max", p3_name))?.f64()?.get(idx);
        let p3 = df[p3_name].f64()?.get(idx).unwrap();

        let p4_min = df.column(&format!("{}_min", p4_name))?.f64()?.get(idx);
        let p4_max = df.column(&format!("{}_max", p4_name))?.f64()?.get(idx);
        let p4 = df[p4_name].f64()?.get(idx).unwrap();

        let epoch_val = df["epoch"].f64()?.get(idx).unwrap();
        let epoch = Time::new(epoch_val, "utc", "jd")?;

        let ic = InitialCondition::from_params(
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
        );

        ics.push(ic);
    }

    Ok(ics)
}

#[derive(Debug, Clone)]
pub struct ExposureRow {
    pub expnum: i64,
    pub radeg: f64,
    pub decdeg: f64,
    pub epoch: f64,
}

pub fn read_exposure_metadata(path: &str) -> Result<Vec<ExposureRow>, Box<dyn std::error::Error>> {
    let df = CsvReadOptions::default()
        .with_infer_schema_length(Some(100))
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(path.into()))?
        .finish()?;

    let expnums: Vec<i64> = df["EXPNUM"]
        .i64()?
        .into_no_null_iter()
        .map(|x| x as i64)
        .collect();
    let radeg: Vec<f64> = df["RADEG"].f64()?.into_no_null_iter().collect();
    let decdeg: Vec<f64> = df["DECDEG"].f64()?.into_no_null_iter().collect();
    let epoch: Vec<f64> = df["EPOCH"].f64()?.into_no_null_iter().collect();

    let mut exposures = Vec::with_capacity(df.height());
    for i in 0..df.height() {
        exposures.push(ExposureRow {
            expnum: expnums[i],
            radeg: radeg[i],
            decdeg: decdeg[i],
            epoch: epoch[i],
        });
    }

    Ok(exposures)
}