use std::f64::consts::PI;

const TRANSITION_Z: f64 = 2.0 / 3.0;
const TRANSITION_Z_INV: f64 = 3.0 / 2.0;

#[inline]
fn nside(depth: u8) -> u32 {
    1u32 << depth
}

/// Z-order (Morton) interleave, exactly like cdshealpix::nested::gpu::ij2z (but widened).
#[inline]
fn ij2z(mut i: u32, mut j: u32) -> u64 {
    i |= j << 16;
    j = (i ^ (i >> 8)) & 0x0000_FF00; i = i ^ j ^ (j << 8);
    j = (i ^ (i >> 4)) & 0x00F0_00F0; i = i ^ j ^ (j << 4);
    j = (i ^ (i >> 2)) & 0x0C0C_0C0C; i = i ^ j ^ (j << 2);
    j = (i ^ (i >> 1)) & 0x2222_2222; i = i ^ j ^ (j << 1);
    i as u64
}

/// Matches cdshealpix::nested::gpu::xpm1_and_q
#[inline]
fn xpm1_and_q(x: f64, y: f64) -> (f64, u8) {
    let x_neg = (x < 0.0) as u8;
    let y_neg = (y < 0.0) as u8;
    let q = (x_neg + y_neg) | (y_neg << 1);

    // lon in [0, pi/2]
    let lon = y.abs().atan2(x.abs());
    let x02 = lon * 4.0 / PI; // in [0,2]

    if x_neg != y_neg {
        (1.0 - x02, q)
    } else {
        (x02 - 1.0, q)
    }
}

/// Matches cdshealpix::nested::gpu::{one_minus_z_pos, one_minus_z_neg}
#[inline]
fn one_minus_z_pos(x: f64, y: f64, z: f64) -> f64 {
    debug_assert!(z > 0.0);
    let d2 = x * x + y * y;
    if d2 < 1e-1 {
        d2 * (0.5 + d2 * (0.125 + d2 * (0.0625 + d2 * (0.0390625 + d2 * 0.02734375))))
    } else {
        1.0 - z
    }
}

#[inline]
fn one_minus_z_neg(x: f64, y: f64, z: f64) -> f64 {
    debug_assert!(z < 0.0);
    let d2 = x * x + y * y;
    if d2 < 1e-1 {
        d2 * (0.5 + d2 * (0.125 + d2 * (0.0625 + d2 * (0.0390625 + d2 * 0.02734375))))
    } else {
        z + 1.0
    }
}

/// Vec (unit) -> NESTED ipix, cdshealpix-gpu-compatible logic, but returns u64 and supports depth=15+
pub fn vec3_to_nested_ipix(depth: u8, x: f64, y: f64, z: f64) -> u64 {
    debug_assert!((-1.0..=1.0).contains(&x));
    debug_assert!((-1.0..=1.0).contains(&y));
    debug_assert!((-1.0..=1.0).contains(&z));
    debug_assert!((x * x + y * y + z * z - 1.0).abs() < 1e-9);

    let nside = nside(depth);
    let half_nside = nside as f64 * 0.5;

    let (x_pm1, q) = xpm1_and_q(x, y);

    let (d0h, x_in_d0c, y_in_d0c) = if z > TRANSITION_Z {
        // North polar cap (Collignon)
        let s = (3.0 * one_minus_z_pos(x, y, z)).sqrt();
        (q, x_pm1 * s, 2.0 - s)
    } else if z < -TRANSITION_Z {
        // South polar cap (Collignon)
        let s = (3.0 * one_minus_z_neg(x, y, z)).sqrt();
        (q + 8, x_pm1 * s, s)
    } else {
        // Equatorial region (cylindrical equal-area)
        let y_pm1 = z * TRANSITION_Z_INV;

        // Determine which “diamond” we are in (exact cdshealpix logic)
        let q01 = (x_pm1 >  y_pm1) as u8;   // 0/1
        let q12 = (x_pm1 >= -y_pm1) as u8;  // 0\1
        let q03 = 1 - q12;                 // 1\0
        let q1  = q01 & q12;               // 1 if q1 else 0

        // Local projected coords inside base cell reference frame
        let x_proj = x_pm1 - ((q01 + q12) as i8 - 1) as f64;
        let y_proj = y_pm1 + (q01 + q03) as f64;

        let d0h = ((q01 + q03) << 2) + ((q + q1) & 3);
        (d0h, x_proj, y_proj)
    };

    // Coords inside the base cell (this transform matters!)
    let xf = half_nside * (y_in_d0c + x_in_d0c);
    let yf = half_nside * (y_in_d0c - x_in_d0c);

    let mut i = xf as i64;
    let mut j = yf as i64;

    // Clamp rare numerical edge cases (cdshealpix does i==nside -> i-=1)
    if i < 0 { i = 0; }
    if j < 0 { j = 0; }
    if i as u32 >= nside { i = (nside - 1) as i64; }
    if j as u32 >= nside { j = (nside - 1) as i64; }

    let i = i as u32;
    let j = j as u32;

    // Assemble nested hash
    ((d0h as u64) << ((depth as u64) * 2)) | ij2z(i, j)
}

