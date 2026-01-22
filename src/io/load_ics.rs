// use polars::prelude::*;
// use crate::initial_condition::InitialCondition;

// use spacerocks::time::Time;
// use spacerocks::coordinates::Origin;

// pub fn load_initial_conditions(path: &str, method: &str, origin: &str, reference_epoch_jd: f64) -> Result<Vec<InitialCondition>, Box<dyn std::error::Error>> {

//     let mut spacerock_origin = Origin::from_str(origin)?;
//     let mu = spacerock_origin.mu();
    
//     let kind = method.to_lowercase();
//     match kind.as_str() {
//         "keplerian" => read_initial_conditions_kep(path, mu, reference_epoch_jd),
//         "spherical" => read_initial_conditions_sph(path, mu, reference_epoch_jd),
//         _ => Err(Box::new(std::io::Error::new(std::io::ErrorKind::InvalidInput, "Invalid method for reading initial conditions. Use 'keplerian' or 'spherical'."))),
//     }
// }

// pub fn read_initial_conditions_kep(path: &str, mu: f64, reference_epoch_jd: f64) -> Result<Vec<InitialCondition>, Box<dyn std::error::Error>> {

//     let initial_conditions = CsvReadOptions::default()
//         .with_has_header(true)
//         .try_into_reader_with_file_path(Some(path.into()))?
//         .finish()?;

//     let epoch = Time::new(reference_epoch_jd, "utc", "jd")?;

//     let mut ics: Vec<InitialCondition> = Vec::new();
//     for idx in 0..initial_conditions.height() {
//         // id is a string column
//         let id = initial_conditions["id"].utf8()?.get(idx).unwrap();
//         let q = initial_conditions["q"].f64()?.get(idx).unwrap();
//         let e = initial_conditions["e"].f64()?.get(idx).unwrap();
//         let inc = initial_conditions["inc"].f64()?.get(idx).unwrap();
//         let M0 = initial_conditions["M0"].f64()?.get(idx).unwrap();
//         let kappa = initial_conditions["kappa"].f64()?.get(idx).unwrap() as i32;
//         let ic = InitialCondition::from_elements(id.to_string(), q, e, inc, M0, kappa, epoch.clone(), mu)?;
//         ics.push(ic);
//     }
//     Ok(ics)
// }


// pub fn read_initial_conditions_sph(path: &str, mu: f64, reference_epoch_jd: f64) -> Result<Vec<InitialCondition>, Box<dyn std::error::Error>> {
//     let initial_conditions = CsvReadOptions::default()
//         .with_has_header(true)
//         .try_into_reader_with_file_path(Some(path.into()))?
//         .finish()?;

//     let epoch = Time::new(reference_epoch_jd, "utc", "jd")?;

//     let mut ics: Vec<InitialCondition> = Vec::new();
//     for idx in 0..initial_conditions.height() {
//         let id = initial_conditions["id"].utf8()?.get(idx).unwrap();
//         let r = initial_conditions["r"].f64()?.get(idx).unwrap();
//         let vr = initial_conditions["vr"].f64()?.get(idx).unwrap();
//         let vo = initial_conditions["vo"].f64()?.get(idx).unwrap();
//         let inc = initial_conditions["inc"].f64()?.get(idx).unwrap();
//         let kappa = initial_conditions["kappa"].f64()?.get(idx).unwrap() as i32;
//         let ic = InitialCondition::from_spherical(id.to_string(), r, vr, vo, inc, kappa, epoch.clone(), mu)?;
//         ics.push(ic);
// }
//     Ok(ics)
// }


use polars::prelude::*;
use crate::initial_condition::InitialCondition;

use spacerocks::coordinates::Origin;
use spacerocks::time::Time;

pub fn load_initial_conditions(
    path: &str,
    method: &str,
    origin: &str,
    reference_epoch_jd: f64,
) -> Result<Vec<InitialCondition>, Box<dyn std::error::Error>> {
    let spacerock_origin = Origin::from_str(origin)?;
    let mu = spacerock_origin.mu();

    match method.to_lowercase().as_str() {
        "keplerian" => read_initial_conditions_kep(path, mu, reference_epoch_jd),
        "spherical" => read_initial_conditions_sph(path, mu, reference_epoch_jd),
        _ => Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Invalid method for reading initial conditions. Use 'keplerian' or 'spherical'.",
        ))),
    }
}

fn read_csv(path: &str) -> Result<DataFrame, Box<dyn std::error::Error>> {
    Ok(CsvReadOptions::default()
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(path.into()))?
        .finish()?)
}

/// Read an "id" value from a DataFrame row, robust across Polars versions and CSV typing.
/// - String column -> value
/// - UInt8 column -> single char
/// - Numeric column -> stringified number
fn read_id(df: &DataFrame, row: usize) -> Result<String, Box<dyn std::error::Error>> {
    let col = df.column("id")?; // &Column
    let s = col
        .as_series()
        .ok_or_else(|| PolarsError::ComputeError("Column 'id' is not a Series".into()))?;

    // Preferred: actual string column
    if let Ok(ca) = s.str() {
        return ca
            .get(row)
            .map(str::to_owned)
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }

    // Legacy: single-byte "char"
    if let Ok(ca) = s.u8() {
        return ca
            .get(row)
            .map(|b| (b as char).to_string())
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }

    // Common CSV mess: numeric IDs; stringify them.
    if let Ok(ca) = s.i64() {
        return ca
            .get(row)
            .map(|v| v.to_string())
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }
    if let Ok(ca) = s.u64() {
        return ca
            .get(row)
            .map(|v| v.to_string())
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }
    if let Ok(ca) = s.f64() {
        return ca
            .get(row)
            .map(|v| v.to_string())
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }

    Err(format!(
        "Column 'id' has unsupported dtype: {:?} (expected String/UInt8/numeric)",
        s.dtype()
    )
    .into())
}

pub fn read_initial_conditions_kep(
    path: &str,
    mu: f64,
    reference_epoch_jd: f64,
) -> Result<Vec<InitialCondition>, Box<dyn std::error::Error>> {
    let df = read_csv(path)?;
    let epoch = Time::new(reference_epoch_jd, "utc", "jd")?;

    let q = df.column("q")?.as_series().unwrap().f64()?;
    let e = df.column("e")?.as_series().unwrap().f64()?;
    let inc = df.column("inc")?.as_series().unwrap().f64()?;
    let m0 = df.column("M0")?.as_series().unwrap().f64()?;
    let kappa = df.column("kappa")?.as_series().unwrap().f64()?; // keep as f64 then cast per-row

    let mut ics = Vec::with_capacity(df.height());
    for row in 0..df.height() {
        let id = read_id(&df, row)?;

        let ic = InitialCondition::from_elements(
            id,
            q.get(row).ok_or_else(|| format!("Row {row}: 'q' is null"))?,
            e.get(row).ok_or_else(|| format!("Row {row}: 'e' is null"))?,
            inc.get(row).ok_or_else(|| format!("Row {row}: 'inc' is null"))?,
            m0.get(row).ok_or_else(|| format!("Row {row}: 'M0' is null"))?,
            kappa
                .get(row)
                .ok_or_else(|| format!("Row {row}: 'kappa' is null"))? as i32,
            epoch.clone(),
            mu,
        )?;
        ics.push(ic);
    }
    Ok(ics)
}

pub fn read_initial_conditions_sph(
    path: &str,
    mu: f64,
    reference_epoch_jd: f64,
) -> Result<Vec<InitialCondition>, Box<dyn std::error::Error>> {
    let df = read_csv(path)?;
    let epoch = Time::new(reference_epoch_jd, "utc", "jd")?;

    let r = df.column("r")?.as_series().unwrap().f64()?;
    let vr = df.column("vr")?.as_series().unwrap().f64()?;
    let vo = df.column("vo")?.as_series().unwrap().f64()?;
    let inc = df.column("inc")?.as_series().unwrap().f64()?;
    let kappa = df.column("kappa")?.as_series().unwrap().f64()?;

    let mut ics = Vec::with_capacity(df.height());
    for row in 0..df.height() {
        let id = read_id(&df, row)?;

        let ic = InitialCondition::from_spherical(
            id,
            r.get(row).ok_or_else(|| format!("Row {row}: 'r' is null"))?,
            vr.get(row).ok_or_else(|| format!("Row {row}: 'vr' is null"))?,
            vo.get(row).ok_or_else(|| format!("Row {row}: 'vo' is null"))?,
            inc.get(row).ok_or_else(|| format!("Row {row}: 'inc' is null"))?,
            kappa
                .get(row)
                .ok_or_else(|| format!("Row {row}: 'kappa' is null"))? as i32,
            epoch.clone(),
            mu,
        )?;
        ics.push(ic);
    }
    Ok(ics)
}
