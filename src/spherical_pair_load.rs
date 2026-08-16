//! The adapter between the existing detection loader and Stage 3's [`Obs`].
//!
//! This is the seam where units, frames and time scales get decided, so it is written to make each
//! of those a stated decision rather than an inherited assumption. What was verified about
//! `io::load_detections` before writing it:
//!
//! * ✅ **Frame.** The obscode path builds the observer as
//!   `observatory.at(&epoch, "J2000", "ssb", &kernel)` -- **barycentric**, which is the origin C2
//!   requires. `mu` and the origin move together, and a heliocentric observer against
//!   `Origin::SSB.mu()` fails the grid and ladder by ~0.067%.
//!   🔴 The `obs_x/obs_y/obs_z` fallback path takes the observer straight from the CSV with **no
//!   frame check at all**. A file written in heliocentric coordinates loads without complaint.
//!   [`ObsSet::from_detections`] cannot detect that; it is stated here and re-stated in the config.
//!
//! * ✅ **Time scale.** `Detection::epoch` is `epoch.tdb().jd()`, a TDB Julian date, which is what
//!   [`Obs::epoch`] means. No conversion, and none should be added.
//!
//! * 🔴 **Uncertainty units are ASYMMETRIC in the loader.** `ast_ucty` is arcsec on input and is
//!   converted to radians (with a guard that rejects already-radian input). `ra_ucty` and
//!   `dec_ucty` are stored **exactly as they appear in the file, unconverted and unguarded**.
//!   Reading them as radians is a 206265x error, and it would land in both the gate radius and the
//!   chi2 denominator -- i.e. it would look like a beautifully tight gate.
//!
//! Because of that last point the sigma source is an explicit [`SigmaSource`] rather than a
//! best-effort search: which column is authoritative is a property of the dataset, and a wrong
//! guess is silent. Rubin reports per-source `raErr`/`decErr` and those are what it should use;
//! a pipeline-synthesised `astrom_sig = 0.805/SNR` is **not** Rubin's sigma and should not be
//! mistaken for it.

use std::collections::BTreeMap;

use crate::detection::Detection;
use crate::spherical_pair_anchor::Obs;

const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

/// Where the per-detection astrometric sigma comes from. No default: see the module docs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SigmaSource {
    /// `Detection::ast_ucty`, which the loader has already converted to radians.
    Ast,
    /// `ra_ucty`/`dec_ucty`, which the loader leaves in the file's units -- assumed **arcsec** and
    /// converted here.
    ///
    /// 🔴 `ra_includes_cos_dec` is not guessable from the data. If the file reports a great-circle
    /// RA error the factor is already applied; if it reports `d(RA)` it is not, and at |dec| = 60
    /// that is a factor of 2 on one axis. Whoever writes the config knows; this code must not
    /// pretend to.
    RaDec { ra_includes_cos_dec: bool },
    /// A single value for every detection, in arcsec. For data that reports no uncertainty at all.
    /// ⚠️ A run leaning on this should say so: sigma sets the gate width AND the chi2 denominator.
    Fixed { arcsec: f64 },
}

/// What the adapter did, so a caller can tell "no detections" from "no usable sigmas".
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LoadReport {
    pub offered: usize,
    pub converted: usize,
    /// Rejected because the chosen sigma source was absent or non-positive on that row.
    pub no_sigma: usize,
    /// Rejected because the line of sight was not a finite unit vector.
    pub bad_geometry: usize,
}

/// Detections in the form Stage 3 consumes, with their grouping.
pub struct ObsSet {
    pub obs: Vec<Obs>,
    pub report: LoadReport,
}

impl ObsSet {
    /// Convert loaded detections, taking sigma from the stated source.
    ///
    /// `dets` must all be on the same reference plane. C2's geometry is a set of vectors, so either
    /// plane works -- but only if `rho_hat` and `observer_position` are on the SAME one, which
    /// `Detection::to_ecliptic`/`to_equatorial` guarantee by rotating both together.
    pub fn from_detections(dets: &[Detection], source: SigmaSource) -> Self {
        let mut obs = Vec::with_capacity(dets.len());
        let mut report = LoadReport { offered: dets.len(), ..Default::default() };
        for (i, d) in dets.iter().enumerate() {
            let sigma = match source {
                SigmaSource::Ast => d.ast_ucty,
                SigmaSource::RaDec { ra_includes_cos_dec } => {
                    match (d.ra_ucty, d.dec_ucty) {
                        (Some(sra), Some(sdec)) => {
                            // the loader does NOT convert these -- arcsec in the file, arcsec here
                            let dec = d.rho_hat.z.clamp(-1.0, 1.0).asin();
                            let sra = if ra_includes_cos_dec { sra } else { sra * dec.cos() };
                            // one per-coordinate sigma from two axes
                            Some(((sra * sra + sdec * sdec) * 0.5).sqrt() * ARCSEC)
                        }
                        _ => None,
                    }
                }
                SigmaSource::Fixed { arcsec } => Some(arcsec * ARCSEC),
            };
            let Some(sigma) = sigma.filter(|s| s.is_finite() && *s > 0.0) else {
                report.no_sigma += 1;
                continue;
            };
            let n = d.rho_hat.norm();
            if !n.is_finite() || (n - 1.0).abs() > 1e-6 || !d.observer_position.iter().all(|v| v.is_finite())
            {
                report.bad_geometry += 1;
                continue;
            }
            obs.push(Obs {
                id: i as u32,
                epoch: d.epoch,
                rho_hat: d.rho_hat / n,
                observer: d.observer_position,
                sigma,
            });
        }
        report.converted = obs.len();
        ObsSet { obs, report }
    }

    /// Group detections into visits by epoch.
    ///
    /// 🔴 An exact float epoch is not a key. Detections from one exposure carry the same nominal
    /// time but can differ in the last bits, and separately, a pipeline that rounds differently
    /// produces epochs a few seconds apart for the same exposure. `tol_seconds` is therefore
    /// explicit; the same lesson cost this project a cross-match once, where a plus-or-minus 60 s
    /// window was NARROWER than a systematic 69 s offset between two time scales.
    pub fn visits(&self, tol_seconds: f64) -> Vec<Vec<usize>> {
        let mut order: Vec<usize> = (0..self.obs.len()).collect();
        order.sort_by(|&a, &b| self.obs[a].epoch.total_cmp(&self.obs[b].epoch));
        let tol_days = tol_seconds / 86_400.0;
        let mut out: Vec<Vec<usize>> = Vec::new();
        for i in order {
            match out.last_mut() {
                Some(v) if (self.obs[i].epoch - self.obs[v[0]].epoch).abs() <= tol_days => v.push(i),
                _ => out.push(vec![i]),
            }
        }
        out
    }

    /// Group visits into nights.
    ///
    /// ⚠️ **This is a fallback, not the right answer.** The correct night key is the survey's own
    /// `dayObs`: a floor of the epoch disagrees with it on real data (48 of 924 in one audit), and
    /// a validator that recomputed the key the same wrong way agreed with the bug. Where the input
    /// carries a night label, group on that and do not call this.
    ///
    /// `boundary_days` is the offset of the night boundary from the JD tick; 0.5 puts it at local
    /// noon for a UT-referenced date, which is the convention `floor(mjd - 0.5)` encodes.
    pub fn nights(&self, visits: &[Vec<usize>], boundary_days: f64) -> BTreeMap<i64, Vec<usize>> {
        let mut out: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
        for (vi, v) in visits.iter().enumerate() {
            let key = (self.obs[v[0]].epoch - 2_400_000.5 - boundary_days).floor() as i64;
            out.entry(key).or_default().push(vi);
        }
        out
    }

    /// The population sigma used to size a KD query so one tree serves every pair.
    ///
    /// The cap must admit every pair the exact per-pair gate might keep, and the exact test then
    /// rejects the rest. Under-admission is invisible by construction: nothing downstream can
    /// observe a pair that was never proposed.
    ///
    /// 🔴 A PERCENTILE CANNOT BOUND A QUADRATURE, and this took `sqrt(2) * p97.5` until
    /// 2026-08-16. A pair's sigma is `sqrt(s_i^2 + s_j^2)`; bounding it by `sqrt(2) * s_q`
    /// requires BOTH `s_i` and `s_j` to be at or below `s_q`, and at q = 0.975 that fails for
    /// ~0.06% of pairs -- the pairs with the two widest sigmas, which are exactly the wide,
    /// high-rate, disproportionately chance pairs the gate finds easiest. Dropping them
    /// FLATTERS the gate: `RESULT_c2_kdcap_correction.md` found the same defect on the Python
    /// side (there capping on ONE detection's sigma, 1234 tracklets against 1318) and its
    /// corrected numbers moved UPWARD.
    ///
    /// `sqrt(2) * max` is strict: for any pair, `sqrt(s_i^2 + s_j^2) <= sqrt(2) * max(s)`. The
    /// note measured the cost of being strict on real Rubin sky at 1.1% in radius and 2% in
    /// area, against a sigma distribution tight enough that there was never a speed argument
    /// for the percentile.
    pub fn sigma_cap(&self) -> f64 {
        // NaN-safe: `total_cmp` orders NaN above every real value, so a NaN sigma would become
        // the max and blow the cap open. Filter rather than sort.
        self.obs
            .iter()
            .map(|o| o.sigma)
            .filter(|s| s.is_finite())
            .fold(0.0_f64, f64::max)
            * std::f64::consts::SQRT_2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detection::ReferencePlane;
    use nalgebra::Vector3;

    fn det(epoch: f64, sigma_ast: Option<f64>, ra_ucty: Option<f64>, dec_ucty: Option<f64>) -> Detection {
        Detection {
            rho_hat: Vector3::new(0.2, 0.9, 0.1).normalize(),
            observer_position: Vector3::new(0.4, -0.9, 0.0),
            epoch,
            reference_plane: ReferencePlane::Ecliptic,
            detid: None,
            trackid: None,
            observer_velocity: None,
            mag: None,
            mag_ucty: None,
            filter: None,
            ast_ucty: sigma_ast,
            ra_ucty,
            dec_ucty,
            objid: None,
            obscode: None,
            rho_hat_dot_observer_position: 0.0,
            observer_distance_squared: 0.0,
        }
    }

    #[test]
    fn ast_ucty_is_taken_as_radians_and_radec_as_arcsec() {
        // 🔴 the asymmetry the loader actually has, asserted so a future refactor cannot quietly
        // make the two paths agree in the wrong direction
        let rad = 0.12 * ARCSEC;
        let a = ObsSet::from_detections(&[det(2_460_000.5, Some(rad), None, None)], SigmaSource::Ast);
        assert_eq!(a.obs.len(), 1);
        assert!((a.obs[0].sigma - rad).abs() < 1e-18, "ast_ucty must pass through unscaled");

        // same physical sigma, expressed the way the ra/dec columns express it
        let b = ObsSet::from_detections(
            &[det(2_460_000.5, None, Some(0.12), Some(0.12))],
            SigmaSource::RaDec { ra_includes_cos_dec: true },
        );
        assert_eq!(b.obs.len(), 1);
        assert!(
            (b.obs[0].sigma - rad).abs() / rad < 1e-12,
            "ra/dec must be converted from arcsec: got {} want {rad}",
            b.obs[0].sigma
        );
    }

    #[test]
    fn the_cos_dec_convention_changes_the_answer_and_so_must_be_stated() {
        // dec ~ +5.7 deg for this rho_hat, so cos(dec) is close to 1; use a steeper one
        let mut d = det(2_460_000.5, None, Some(0.20), Some(0.10));
        d.rho_hat = Vector3::new(0.1, 0.4, 0.91).normalize(); // dec ~ 65 deg
        let with = ObsSet::from_detections(&[d.clone()], SigmaSource::RaDec { ra_includes_cos_dec: true });
        let without = ObsSet::from_detections(&[d], SigmaSource::RaDec { ra_includes_cos_dec: false });
        let ratio = with.obs[0].sigma / without.obs[0].sigma;
        assert!(ratio > 1.2, "at high |dec| the convention must matter; ratio {ratio}");
    }

    #[test]
    fn a_detection_without_the_chosen_sigma_is_dropped_and_counted() {
        let set = ObsSet::from_detections(
            &[
                det(2_460_000.5, Some(0.12 * ARCSEC), None, None),
                det(2_460_000.5, None, None, None),
                det(2_460_000.5, Some(0.0), None, None), // zero is not a sigma
            ],
            SigmaSource::Ast,
        );
        assert_eq!(set.obs.len(), 1);
        assert_eq!(set.report.offered, 3);
        assert_eq!(set.report.no_sigma, 2);
        assert_eq!(set.report.converted, 1);
    }

    #[test]
    fn visits_group_by_epoch_within_the_stated_tolerance() {
        let s = 1.0 / 86_400.0;
        let set = ObsSet::from_detections(
            &[
                det(2_460_000.500_00, Some(1e-6), None, None),
                det(2_460_000.500_00 + 0.3 * s, Some(1e-6), None, None), // same exposure
                det(2_460_000.500_00 + 40.0 * s, Some(1e-6), None, None), // a different one
            ],
            SigmaSource::Ast,
        );
        let v = set.visits(5.0);
        assert_eq!(v.len(), 2, "5 s tolerance must merge 0.3 s and split 40 s");
        assert_eq!(v[0].len(), 2);
        assert_eq!(v[1].len(), 1);
        // 🔴 and an exact-equality grouping would have produced three visits
        assert_eq!(set.visits(0.0).len(), 3);
    }

    #[test]
    fn nights_split_at_the_stated_boundary() {
        let set = ObsSet::from_detections(
            &[
                det(2_460_000.3, Some(1e-6), None, None),
                det(2_460_000.7, Some(1e-6), None, None),
                det(2_460_001.7, Some(1e-6), None, None),
            ],
            SigmaSource::Ast,
        );
        // The first two are MJD 59999.8 and 60000.2: 0.4 d apart, straddling the MJD tick.
        // ⭐ That is exactly the case the noon boundary exists for -- `floor(mjd - 0.5)` keeps one
        // night together across midnight UT, where `floor(mjd)` would cut it in half and turn a
        // two-visit night (which can form an anchor) into two single-visit nights (which cannot).
        let v = set.visits(5.0);
        assert_eq!(v.len(), 3);
        let n = set.nights(&v, 0.5);
        assert_eq!(n.len(), 2, "the noon boundary must keep 59999.8 and 60000.2 together");
        assert_eq!(n.values().next().unwrap().len(), 2, "and as one two-visit night");
        // the boundary is a convention, and the wrong one silently destroys anchors
        assert_eq!(set.nights(&v, 0.0).len(), 3, "a midnight boundary splits that night");
    }

    #[test]
    fn the_sigma_cap_strictly_bounds_every_pair_quadrature() {
        let mut ds = Vec::new();
        for i in 0..100 {
            ds.push(det(2_460_000.5, Some((0.05 + 0.001 * i as f64) * ARCSEC), None, None));
        }
        let set = ObsSet::from_detections(&ds, SigmaSource::Ast);
        let cap = set.sigma_cap();
        let smax = (0.05 + 0.001 * 99.0) * ARCSEC;
        assert!((cap / (smax * std::f64::consts::SQRT_2) - 1.0).abs() < 1e-12);

        // 🔴 The property the old p97.5 form did NOT have. Every pair, including the two widest
        // together, must fall inside the cap -- that pair is precisely the one a percentile
        // misses, and it is a wide, high-rate, disproportionately chance pair.
        let s: Vec<f64> = set.obs.iter().map(|o| o.sigma).collect();
        for i in 0..s.len() {
            for j in 0..s.len() {
                assert!(
                    (s[i] * s[i] + s[j] * s[j]).sqrt() <= cap * (1.0 + 1e-12),
                    "pair ({i},{j}) escapes the cap"
                );
            }
        }
        // and the superseded percentile form would have failed exactly that
        let p975 = (0.05 + 0.001 * 97.0) * ARCSEC;
        let old_cap = p975 * std::f64::consts::SQRT_2;
        assert!(
            (smax * smax + smax * smax).sqrt() > old_cap,
            "regression guard: the p97.5 cap must be shown to miss the widest pair"
        );
    }

    #[test]
    fn a_nan_sigma_cannot_blow_the_cap_open() {
        // total_cmp orders NaN above every real value, so a max taken by sorting would return
        // NaN and the cap would admit the whole sky.
        let ds = vec![
            det(2_460_000.5, Some(0.10 * ARCSEC), None, None),
            det(2_460_000.5, Some(f64::NAN), None, None),
            det(2_460_000.5, Some(0.20 * ARCSEC), None, None),
        ];
        let set = ObsSet::from_detections(&ds, SigmaSource::Ast);
        let cap = set.sigma_cap();
        assert!(cap.is_finite(), "a NaN sigma must not reach the cap");
    }
}
