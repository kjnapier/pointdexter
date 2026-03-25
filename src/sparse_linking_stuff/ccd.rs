use std::collections::HashMap;
use once_cell::sync::Lazy;

use std::f64::consts::PI;

pub static CCD_BOUNDS: Lazy<HashMap<&'static str, (f64, f64, f64, f64)>> = Lazy::new(|| {
    let mut m = HashMap::new();
    m.insert("N1", (-1.0811, -0.782681, -0.157306, -0.00750506));
    m.insert("N2", (-0.771362, -0.472493, -0.157385, -0.00749848));
    m.insert("N3", (-0.461205, -0.161464, -0.157448, -0.00749265));
    m.insert("N4", (-0.150127, 0.149894, -0.15747, -0.00749085));
    m.insert("N5", (0.161033, 0.460796, -0.157638, -0.0074294));
    m.insert("N6", (0.472171, 0.771045, -0.157286, -0.00740563));
    m.insert("N7", (0.782398, 1.08083, -0.157141, -0.0074798));
    m.insert("N8", (-0.92615, -0.627492, -0.321782, -0.172004));
    m.insert("N9", (-0.616455, -0.317043, -0.322077, -0.172189));
    m.insert("N10", (-0.305679, -0.00571999, -0.322071, -0.17217));
    m.insert("N11", (0.00565427, 0.305554, -0.322243, -0.172254));
    m.insert("N12", (0.31684, 0.616183, -0.322099, -0.172063));
    m.insert("N13", (0.627264, 0.925858, -0.321792, -0.171887));
    m.insert("N14", (-0.926057, -0.62726, -0.485961, -0.336213));
    m.insert("N15", (-0.616498, -0.317089, -0.486444, -0.336606));
    m.insert("N16", (-0.30558, -0.00578257, -0.486753, -0.336864));
    m.insert("N17", (0.00532179, 0.305123, -0.486814, -0.33687));
    m.insert("N18", (0.316662, 0.616018, -0.486495, -0.336537));
    m.insert("N19", (0.62708, 0.92578, -0.485992, -0.336061));
    m.insert("N20", (-0.770814, -0.471826, -0.650617, -0.500679));
    m.insert("N21", (-0.460777, -0.161224, -0.650817, -0.501097));
    m.insert("N22", (-0.149847, 0.149886, -0.650816, -0.501308));
    m.insert("N23", (0.161001, 0.460566, -0.650946, -0.501263));
    m.insert("N24", (0.47163, 0.770632, -0.650495, -0.500592));
    m.insert("N25", (-0.615548, -0.316352, -0.814774, -0.665052));
    m.insert("N26", (-0.305399, -0.00591217, -0.814862, -0.665489));
    m.insert("N27", (0.00550714, 0.304979, -0.815022, -0.665418));
    m.insert("N28", (0.316126, 0.615276, -0.814707, -0.664908));
    m.insert("N29", (-0.46018, -0.16101, -0.97887, -0.829315));
    m.insert("N31", (0.160884, 0.460147, -0.978775, -0.829426));
    m.insert("S1", (-1.08096, -0.782554, 0.00715956, 0.15689));
    m.insert("S2", (-0.7713, -0.47242, 0.0074194, 0.157269));
    m.insert("S3", (-0.4611, -0.161377, 0.00723009, 0.157192));
    m.insert("S4", (-0.149836, 0.150222, 0.00737069, 0.157441));
    m.insert("S5", (0.161297, 0.461031, 0.0072399, 0.1572));
    m.insert("S6", (0.472537, 0.771441, 0.00728934, 0.157137));
    m.insert("S7", (0.782516, 1.08097, 0.00742809, 0.15709));
    m.insert("S8", (-0.92583, -0.627259, 0.171786, 0.32173));
    m.insert("S9", (-0.616329, -0.31694, 0.171889, 0.321823));
    m.insert("S10", (-0.305695, -0.00579187, 0.172216, 0.322179));
    m.insert("S11", (0.00556739, 0.305472, 0.172237, 0.322278));
    m.insert("S12", (0.316973, 0.61631, 0.172015, 0.322057));
    m.insert("S13", (0.627389, 0.925972, 0.171749, 0.321672));
    m.insert("S14", (-0.925847, -0.627123, 0.335898, 0.48578));
    m.insert("S15", (-0.616201, -0.316839, 0.336498, 0.486438));
    m.insert("S16", (-0.305558, -0.00574858, 0.336904, 0.486749));
    m.insert("S17", (0.00557115, 0.305423, 0.33675, 0.486491));
    m.insert("S18", (0.316635, 0.615931, 0.33649, 0.486573));
    m.insert("S19", (0.627207, 0.925969, 0.336118, 0.485923));
    m.insert("S20", (-0.770675, -0.471718, 0.500411, 0.65042));
    m.insert("S21", (-0.46072, -0.161101, 0.501198, 0.650786));
    m.insert("S22", (-0.149915, 0.14982, 0.501334, 0.650856));
    m.insert("S23", (0.160973, 0.460482, 0.501075, 0.650896));
    m.insert("S24", (0.47167, 0.770647, 0.50045, 0.650441));
    m.insert("S25", (-0.615564, -0.316325, 0.66501, 0.814674));
    m.insert("S26", (-0.30512, -0.0056517, 0.665531, 0.81505));
    m.insert("S27", (0.00560886, 0.305082, 0.665509, 0.815022));
    m.insert("S28", (0.316158, 0.615391, 0.665058, 0.814732));
    m.insert("S29", (-0.46021, -0.160988, 0.829248, 0.978699));
    m.insert("S30", (-0.150043, 0.149464, 0.829007, 0.978648));
    m.insert("S31", (0.160898, 0.460111, 0.82932, 0.978804));
    m
});

pub static CCD_NUM: Lazy<HashMap<&'static str, i32>> = Lazy::new(|| {
    let mut m = HashMap::new();
    m.insert("N1", 32); m.insert("N2", 33); m.insert("N3", 34); m.insert("N4", 35); m.insert("N5", 36); m.insert("N6", 37); m.insert("N7", 38);
    m.insert("N8", 39); m.insert("N9", 40); m.insert("N10", 41); m.insert("N11", 42); m.insert("N12", 43); m.insert("N13", 44);
    m.insert("N14", 45); m.insert("N15", 46); m.insert("N16", 47); m.insert("N17", 48); m.insert("N18", 49); m.insert("N19", 50);
    m.insert("N20", 51); m.insert("N21", 52); m.insert("N22", 53); m.insert("N23", 54); m.insert("N24", 55); m.insert("N25", 56);
    m.insert("N26", 57); m.insert("N27", 58); m.insert("N28", 59); m.insert("N29", 60); m.insert("N30", 61); m.insert("N31", 62);
    m.insert("S1", 25); m.insert("S2", 26); m.insert("S3", 27); m.insert("S4", 28); m.insert("S5", 29); m.insert("S6", 30); m.insert("S7", 31);
    m.insert("S8", 19); m.insert("S9", 20); m.insert("S10", 21); m.insert("S11", 22); m.insert("S12", 23); m.insert("S13", 24);
    m.insert("S14", 13); m.insert("S15", 14); m.insert("S16", 15); m.insert("S17", 16); m.insert("S18", 17); m.insert("S19", 18);
    m.insert("S20", 8); m.insert("S21", 9); m.insert("S22", 10); m.insert("S23", 11); m.insert("S24", 12); m.insert("S25", 4);
    m.insert("S26", 5); m.insert("S27", 6); m.insert("S28", 7); m.insert("S29", 1); m.insert("S30", 2); m.insert("S31", 3);
    m.insert("None", -99);
    m
});

// /// Compute the CCD that contains the point (ra, dec) relative to the exposure center.
// /// All coordinates should be in **radians**.
// pub fn compute_chip(rock_ra: f64, rock_dec: f64, exp_ra: f64, exp_dec: f64) -> (String, i32) {
//     // ΔRA in degrees, wrapped to (-180, 180], scaled by cos(dec)
//     let mut delta_ra = rock_ra - exp_ra;
//     if delta_ra > PI {
//         delta_ra -= 2.0 * PI;
//     }
//     let delta_ra_deg = delta_ra.to_degrees() * exp_dec.cos();

//     // ΔDec in degrees, wrapped to (-180, 180]
//     let mut delta_dec = rock_dec - exp_dec;
//     if delta_dec > PI {
//         delta_dec -= 2.0 * PI;
//     }
//     let delta_dec_deg = delta_dec.to_degrees();

//     for (ccd_name, &(xmin, xmax, ymin, ymax)) in CCD_BOUNDS.iter() {
//         if delta_ra_deg > xmin && delta_ra_deg < xmax && delta_dec_deg > ymin && delta_dec_deg < ymax {
//             let ccd_num = *CCD_NUM.get(ccd_name).unwrap_or(&-99);
//             return (ccd_name.to_string(), ccd_num);
//         }
//     }

//     ("None".to_string(), -99)
// }



const PIXEL_SCALE_RAD: f64 = 0.2637 * PI / (180.0 * 3600.0); // radians per pixel
const EDGE_MARGIN_DEG: f64 = 20.0 * PIXEL_SCALE_RAD * 180.0 / PI; // 20 pixels in degrees

// If a source falls near a the edge of a ccd, treat it as if it is not on a valid chip
pub fn compute_chip(rock_ra: f64, rock_dec: f64, exp_ra: f64, exp_dec: f64) -> (String, i32) {
    // ΔRA in degrees, wrapped to (-180, 180], scaled by cos(dec)
    let mut delta_ra = rock_ra - exp_ra;
    if delta_ra > PI {
        delta_ra -= 2.0 * PI;
    }
    let delta_ra_deg = delta_ra.to_degrees() * exp_dec.cos();

    // ΔDec in degrees, wrapped to (-180, 180]
    let mut delta_dec = rock_dec - exp_dec;
    if delta_dec > PI {
        delta_dec -= 2.0 * PI;
    }
    let delta_dec_deg = delta_dec.to_degrees();

    for (ccd_name, &(xmin, xmax, ymin, ymax)) in CCD_BOUNDS.iter() {
        if delta_ra_deg  > xmin + EDGE_MARGIN_DEG && delta_ra_deg  < xmax - EDGE_MARGIN_DEG
        && delta_dec_deg > ymin + EDGE_MARGIN_DEG && delta_dec_deg < ymax - EDGE_MARGIN_DEG
        {
            let ccd_num = *CCD_NUM.get(ccd_name).unwrap_or(&-99);
            return (ccd_name.to_string(), ccd_num);
        }
    }

    ("None".to_string(), -99)
}