use polars::prelude::*;
use crate::initial_condition::{InitialCondition, ElementBounds};
use crate::initial_condition_3d::InitialCondition3D;

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
        "spherical"  => read_initial_conditions_sph(path, mu, reference_epoch_jd),
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

fn read_id(df: &DataFrame, row: usize) -> Result<String, Box<dyn std::error::Error>> {
    let col = df.column("id")?;
    let s = col
        .as_series()
        .ok_or_else(|| PolarsError::ComputeError("Column 'id' is not a Series".into()))?;

    if let Ok(ca) = s.str() {
        return ca.get(row).map(str::to_owned)
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }
    if let Ok(ca) = s.u8() {
        return ca.get(row).map(|b| (b as char).to_string())
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }
    if let Ok(ca) = s.i64() {
        return ca.get(row).map(|v| v.to_string())
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }
    if let Ok(ca) = s.u64() {
        return ca.get(row).map(|v| v.to_string())
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }
    if let Ok(ca) = s.f64() {
        return ca.get(row).map(|v| v.to_string())
            .ok_or_else(|| format!("Row {row}: 'id' is null").into());
    }

    Err(format!(
        "Column 'id' has unsupported dtype: {:?} (expected String/UInt8/numeric)",
        s.dtype()
    ).into())
}

/// Try to read an optional f64 column from a DataFrame, returning None if absent.
fn optional_f64_column<'a>(df: &'a DataFrame, name: &str) -> Option<ChunkedArray<Float64Type>> {
    df.column(name).ok()?.as_series()?.f64().ok().cloned()
}

/// Extract an (Option<f64>, Option<f64>) bound pair for a given row from two optional columns.
fn bound_pair(
    min_col: &Option<ChunkedArray<Float64Type>>,
    max_col: &Option<ChunkedArray<Float64Type>>,
    row: usize,
) -> (Option<f64>, Option<f64>) {
    (
        min_col.as_ref().and_then(|c| c.get(row)),
        max_col.as_ref().and_then(|c| c.get(row)),
    )
}

/// Return Some(bounds) only if at least one bound value is present for this row,
/// otherwise None (no bounds columns present at all).
fn make_kep_bounds(
    q_min:  &Option<ChunkedArray<Float64Type>>,
    q_max:  &Option<ChunkedArray<Float64Type>>,
    e_min:  &Option<ChunkedArray<Float64Type>>,
    e_max:  &Option<ChunkedArray<Float64Type>>,
    f_min:  &Option<ChunkedArray<Float64Type>>,
    f_max:  &Option<ChunkedArray<Float64Type>>,
    inc_min: &Option<ChunkedArray<Float64Type>>,
    inc_max: &Option<ChunkedArray<Float64Type>>,
    row: usize,
) -> Option<ElementBounds> {
    // Only construct bounds if at least one bound column exists
    if [q_min, q_max, e_min, e_max, f_min, f_max, inc_min, inc_max]
        .iter()
        .all(|c| c.is_none())
    {
        return None;
    }
    Some(ElementBounds::kep(
        bound_pair(q_min,   q_max,   row),
        bound_pair(e_min,   e_max,   row),
        bound_pair(f_min,   f_max,   row),
        bound_pair(inc_min, inc_max, row),
    ))
}

fn make_kev_bounds(
    r_min:   &Option<ChunkedArray<Float64Type>>,
    r_max:   &Option<ChunkedArray<Float64Type>>,
    vr_min:  &Option<ChunkedArray<Float64Type>>,
    vr_max:  &Option<ChunkedArray<Float64Type>>,
    vo_min:  &Option<ChunkedArray<Float64Type>>,
    vo_max:  &Option<ChunkedArray<Float64Type>>,
    inc_min: &Option<ChunkedArray<Float64Type>>,
    inc_max: &Option<ChunkedArray<Float64Type>>,
    row: usize,
) -> Option<ElementBounds> {
    if [r_min, r_max, vr_min, vr_max, vo_min, vo_max, inc_min, inc_max]
        .iter()
        .all(|c| c.is_none())
    {
        return None;
    }
    Some(ElementBounds::kev(
        bound_pair(r_min,   r_max,   row),
        bound_pair(vr_min,  vr_max,  row),
        bound_pair(vo_min,  vo_max,  row),
        bound_pair(inc_min, inc_max, row),
    ))
}

pub fn read_initial_conditions_kep(
    path: &str,
    mu: f64,
    reference_epoch_jd: f64,
) -> Result<Vec<InitialCondition>, Box<dyn std::error::Error>> {
    let df = read_csv(path)?;
    let epoch = Time::new(reference_epoch_jd, "utc", "jd")?;

    let q            = df.column("q")?.as_series().unwrap().f64()?.clone();
    let e            = df.column("e")?.as_series().unwrap().f64()?.clone();
    let inc          = df.column("inc")?.as_series().unwrap().f64()?.clone();
    let true_anomaly = df.column("true_anomaly")?.as_series().unwrap().f64()?.clone();
    let kappa        = df.column("kappa")?.as_series().unwrap().f64()?.clone();

    let q_min   = optional_f64_column(&df, "q_min");
    let q_max   = optional_f64_column(&df, "q_max");
    let e_min   = optional_f64_column(&df, "e_min");
    let e_max   = optional_f64_column(&df, "e_max");
    let f_min   = optional_f64_column(&df, "true_anomaly_min");
    let f_max   = optional_f64_column(&df, "true_anomaly_max");
    let inc_min = optional_f64_column(&df, "inc_min");
    let inc_max = optional_f64_column(&df, "inc_max");

    let mut ics = Vec::with_capacity(df.height());
    for row in 0..df.height() {
        let id = read_id(&df, row)?;
        let bounds = make_kep_bounds(
            &q_min, &q_max, &e_min, &e_max,
            &f_min, &f_max, &inc_min, &inc_max,
            row,
        );

        let ic = InitialCondition::from_elements(
            id,
            q.get(row).ok_or_else(|| format!("Row {row}: 'q' is null"))?,
            e.get(row).ok_or_else(|| format!("Row {row}: 'e' is null"))?,
            inc.get(row).ok_or_else(|| format!("Row {row}: 'inc' is null"))?,
            true_anomaly.get(row).ok_or_else(|| format!("Row {row}: 'true_anomaly' is null"))?,
            kappa.get(row).ok_or_else(|| format!("Row {row}: 'kappa' is null"))? as i32,
            epoch.clone(),
            mu,
            bounds,
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

    let r     = df.column("r")?.as_series().unwrap().f64()?.clone();
    let vr    = df.column("vr")?.as_series().unwrap().f64()?.clone();
    let vo    = df.column("vo")?.as_series().unwrap().f64()?.clone();
    let inc   = df.column("inc")?.as_series().unwrap().f64()?.clone();
    let kappa = df.column("kappa")?.as_series().unwrap().f64()?.clone();

    let r_min   = optional_f64_column(&df, "r_min");
    let r_max   = optional_f64_column(&df, "r_max");
    let vr_min  = optional_f64_column(&df, "vr_min");
    let vr_max  = optional_f64_column(&df, "vr_max");
    let vo_min  = optional_f64_column(&df, "vo_min");
    let vo_max  = optional_f64_column(&df, "vo_max");
    let inc_min = optional_f64_column(&df, "inc_min");
    let inc_max = optional_f64_column(&df, "inc_max");

    let mut ics = Vec::with_capacity(df.height());
    for row in 0..df.height() {
        let id = read_id(&df, row)?;
        let bounds = make_kev_bounds(
            &r_min, &r_max, &vr_min, &vr_max,
            &vo_min, &vo_max, &inc_min, &inc_max,
            row,
        );

        let ic = InitialCondition::from_spherical(
            id,
            r.get(row).ok_or_else(|| format!("Row {row}: 'r' is null"))?,
            vr.get(row).ok_or_else(|| format!("Row {row}: 'vr' is null"))?,
            vo.get(row).ok_or_else(|| format!("Row {row}: 'vo' is null"))?,
            inc.get(row).ok_or_else(|| format!("Row {row}: 'inc' is null"))?,
            kappa.get(row).ok_or_else(|| format!("Row {row}: 'kappa' is null"))? as i32,
            epoch.clone(),
            mu,
            bounds,
        )?;
        ics.push(ic);
    }
    Ok(ics)
}