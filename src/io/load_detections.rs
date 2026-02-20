use crate::Detection;
use spacerocks::time::Time;
use spacerocks::observing::Observatory;
use spacerocks::SpiceKernel;


use polars::prelude::*;
use nalgebra::Vector3;

pub fn load_detections(file_path: &str, reference_plane: &str, kernel: &SpiceKernel) -> Result<Vec<Detection>, Box<dyn std::error::Error>> {
    let df = CsvReadOptions::default()
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file_path.into()))?
        .finish()?;

    let n = df.height();

    // ---------- helpers ----------
    #[inline]
    fn has_col(df: &DataFrame, name: &str) -> bool {
        df.get_column_index(name).is_some()
    }

    fn col_series<'a>(df: &'a DataFrame, name: &str) -> PolarsResult<&'a Series> {
        let c = df.column(name)?;
        c.as_series().ok_or_else(|| {
            PolarsError::ComputeError(format!("Column '{name}' is not a Series").into())
        })
    }

    fn req_f64<'a>(df: &'a DataFrame, name: &str) -> PolarsResult<&'a Float64Chunked> {
        let ca = col_series(df, name)?.f64()?;
        if ca.null_count() != 0 {
            return Err(PolarsError::ComputeError(
                format!("Required column '{name}' contains nulls.").into(),
            ));
        }
        Ok(ca)
    }

    fn opt_f64<'a>(df: &'a DataFrame, name: &str) -> PolarsResult<Option<&'a Float64Chunked>> {
        if !has_col(df, name) {
            return Ok(None);
        }
        Ok(Some(col_series(df, name)?.f64()?))
    }

    enum StringCol<'a> {
        Utf8(&'a StringChunked),
        U8(&'a UInt8Chunked),
        I64(&'a Int64Chunked),
        I32(&'a Int32Chunked),
        U64(&'a UInt64Chunked),
        U32(&'a UInt32Chunked),
        // If you ever have categorical, easiest is to cast beforehand (see note below)
    }

    impl<'a> StringCol<'a> {
        fn get(&self, idx: usize) -> Option<String> {
            match self {
                Self::Utf8(ca) => ca.get(idx).map(str::to_owned),
                Self::U8(ca) => ca.get(idx).map(|b| (b as char).to_string()),
                Self::I64(ca) => ca.get(idx).map(|v| v.to_string()),
                Self::I32(ca) => ca.get(idx).map(|v| v.to_string()),
                Self::U64(ca) => ca.get(idx).map(|v| v.to_string()),
                Self::U32(ca) => ca.get(idx).map(|v| v.to_string()),
            }
        }
    }

    fn opt_string_like<'a>(df: &'a DataFrame, name: &str) -> PolarsResult<Option<StringCol<'a>>> {
        if !has_col(df, name) {
            return Ok(None);
        }
        let s = col_series(df, name)?;

        if let Ok(ca) = s.str() {
            return Ok(Some(StringCol::Utf8(ca)));
        }
        if let Ok(ca) = s.u8() {
            return Ok(Some(StringCol::U8(ca)));
        }
        if let Ok(ca) = s.i64() {
            return Ok(Some(StringCol::I64(ca)));
        }
        if let Ok(ca) = s.i32() {
            return Ok(Some(StringCol::I32(ca)));
        }
        if let Ok(ca) = s.u64() {
            return Ok(Some(StringCol::U64(ca)));
        }
        if let Ok(ca) = s.u32() {
            return Ok(Some(StringCol::U32(ca)));
        }

        Err(PolarsError::ComputeError(
            format!(
                "Column '{name}' must be Utf8 or integer type (i64/i32/u64/u32/u8). Got {:?}.",
                s.dtype()
            )
            .into(),
        ))
    }


    // enum StringCol<'a> {
    //     Str(&'a StringChunked),
    //     U8(&'a UInt8Chunked),
    // }
    // impl<'a> StringCol<'a> {
    //     fn get(&self, idx: usize) -> Option<String> {
    //         match self {
    //             Self::Str(ca) => ca.get(idx).map(str::to_owned),
    //             Self::U8(ca) => ca.get(idx).map(|b| (b as char).to_string()),
    //         }
    //     }
    // }

    // fn opt_string_like<'a>(df: &'a DataFrame, name: &str) -> PolarsResult<Option<StringCol<'a>>> {
    //     if !has_col(df, name) {
    //         return Ok(None);
    //     }
    //     let s = col_series(df, name)?;

    //     if let Ok(ca) = s.str() {
    //         return Ok(Some(StringCol::Str(ca)));
    //     }
    //     if let Ok(ca) = s.u8() {
    //         return Ok(Some(StringCol::U8(ca)));
    //     }

    //     Err(PolarsError::ComputeError(
    //         format!("Column '{name}' must be String or UInt8.").into(),
    //     ))
    // }

    // ---------- schema checks ----------
    for &name in &["ra", "dec", "epoch"] {
        if !has_col(&df, name) {
            return Err(format!("Required field '{name}' is missing from detection catalog.").into());
        }
    }

    // ---------- required data ----------
    let ra_deg = req_f64(&df, "ra")?;
    let dec_deg = req_f64(&df, "dec")?;
    let epoch_jd = req_f64(&df, "epoch")?;

    // guard against any odd length mismatch
    for (name, len) in [("ra", ra_deg.len()), ("dec", dec_deg.len()), ("epoch", epoch_jd.len())] {
        if len != n {
            return Err(format!("Column '{name}' length mismatch: expected {n}, got {len}.").into());
        }
    }

    // ---------- observer position inputs (row-level mix allowed) ----------
    let obscode = opt_string_like(&df, "obscode")?;

    // get all unique obscodes if obscode column present
    let mut unique_obscodes = std::collections::HashSet::new();
    if let Some(col) = &obscode {
        for i in 0..n {
            if let Some(code) = col.get(i) {
                unique_obscodes.insert(code);
            }
        }
    }

    // make a hash map of obscodes to observatory objects
    let mut obscode_map = std::collections::HashMap::new();
    for code in unique_obscodes.iter() {
        let obs = Observatory::from_obscode(code)?;
        obscode_map.insert(code.clone(), obs);
    }

    // println!("Unique obscodes in catalog: {:?}", unique_obscodes);

    // xyz columns are optional and may contain nulls per row
    let obs_x = opt_f64(&df, "obs_x")?;
    let obs_y = opt_f64(&df, "obs_y")?;
    let obs_z = opt_f64(&df, "obs_z")?;

    let has_obscode_col = obscode.is_some();
    let has_xyz_cols = obs_x.is_some() && obs_y.is_some() && obs_z.is_some();
    if !has_obscode_col && !has_xyz_cols {
        return Err(
            "Catalog must include either an 'obscode' column or all of 'obs_x','obs_y','obs_z' columns."
                .into(),
        );
    }

    // ---------- optional velocity (row-level nulls allowed) ----------
    let obs_vx = opt_f64(&df, "obs_vx")?;
    let obs_vy = opt_f64(&df, "obs_vy")?;
    let obs_vz = opt_f64(&df, "obs_vz")?;
    let has_v_cols = obs_vx.is_some() && obs_vy.is_some() && obs_vz.is_some();

    // ---------- optional metadata ----------
    let ast_ucty = opt_f64(&df, "ast_ucty")?;
    let magnitude = opt_f64(&df, "mag")?;
    let mag_ucty = opt_f64(&df, "mag_ucty")?;
    let filter = opt_string_like(&df, "filter")?;

    // detid could come in as a string or an int, so we'll just treat it as a string either way.
    let detid = opt_string_like(&df, "detid")?;
    // get the type of the detid column if it exists, and convert to StringCol

    let trackid = opt_string_like(&df, "trackid")?;
    let objid = opt_string_like(&df, "objid")?;
    let obscode = opt_string_like(&df, "obscode")?;
    let ra_ucty = opt_f64(&df, "ra_ucty")?;
    let dec_ucty = opt_f64(&df, "dec_ucty")?;

    // ---------- build ----------
    let mut out = Vec::with_capacity(n);

    for i in 0..n {
        let ra = ra_deg.get(i).unwrap() * std::f64::consts::PI / 180.0;
        let dec = dec_deg.get(i).unwrap() * std::f64::consts::PI / 180.0;
        let epoch = Time::new(epoch_jd.get(i).unwrap(), "utc", "jd")?;

        // Per-row rule:
        // 1) if obscode present on that row -> use it
        // 2) else if obs_x/y/z all present on that row -> use xyz
        // 3) else error
        let observer_position = if let Some(col) = &obscode {
            if let Some(code) = col.get(i) {
                // TODO: replace with real obscode -> (x,y,z) lookup
                // let _ = code;
                // Vector3::new(0.0, 0.0, 0.0)
                let observatory = obscode_map.get(&code).ok_or_else(|| {
                    format!("Row {i}: obscode '{code}' not found in obscode map.")
                })?;
                let observer = observatory.at(&epoch, "J2000", "ssb", &kernel)?;
                observer.position
            } else {
                match (
                    obs_x.and_then(|c| c.get(i)),
                    obs_y.and_then(|c| c.get(i)),
                    obs_z.and_then(|c| c.get(i)),
                ) {
                    (Some(x), Some(y), Some(z)) => Vector3::new(x, y, z),
                    _ => {
                        return Err(format!(
                            "Row {i}: need either non-null 'obscode' or non-null 'obs_x','obs_y','obs_z'."
                        )
                        .into())
                    }
                }
            }
        } else {
            match (
                obs_x.and_then(|c| c.get(i)),
                obs_y.and_then(|c| c.get(i)),
                obs_z.and_then(|c| c.get(i)),
            ) {
                (Some(x), Some(y), Some(z)) => Vector3::new(x, y, z),
                _ => {
                    return Err(format!(
                        "Row {i}: missing observer position; no 'obscode' column and xyz not all present."
                    )
                    .into())
                }
            }
        };

        let mut det = Detection::new(ra, dec, epoch, observer_position);

        // velocity: only set if all 3 columns exist AND all 3 values present on this row
        if has_v_cols {
            if let (Some(vx), Some(vy), Some(vz)) = (
                obs_vx.and_then(|c| c.get(i)),
                obs_vy.and_then(|c| c.get(i)),
                obs_vz.and_then(|c| c.get(i)),
            ) {
                det.observer_velocity = Some(Vector3::new(vx, vy, vz));
            }
        }

        if let Some(ca) = ast_ucty {
            det.ast_ucty = ca.get(i);
        }
        if let Some(ca) = magnitude {
            det.mag = ca.get(i);
        }
        if let Some(col) = &filter {
            det.filter = col.get(i);
        }
        if let Some(col) = &detid {
            det.detid = col.get(i);
        }
        if let Some(col) = &trackid {
            det.trackid = col.get(i);
        }

        if let Some(col) = &mag_ucty {
            det.mag_ucty = col.get(i);
        }

        if let Some(col) = &objid {
            det.objid = col.get(i);
        }

        if let Some(col) = &obscode {
            det.obscode = col.get(i);
        }

        if let Some(col) = &ra_ucty {
            det.ra_ucty = col.get(i);
        }

        if let Some(col) = &dec_ucty {
            det.dec_ucty = col.get(i);
        }

        if reference_plane.to_lowercase() == "ecliptic" {
            det.to_ecliptic();
        } else {
            det.to_equatorial();
        }

        out.push(det);
    }

    Ok(out)
}