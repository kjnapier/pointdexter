//! Stage 2, second half: barycentric projection and the anchor-pair index.
//!
//! `spherical_pair.rs` solves a pair; `spherical_pair_grid.rs` sizes the `(r, rdot)` ladder.
//! What was missing between them is the step that turns a night's detections into the pairs the
//! solver is asked about: project every detection to a barycentric unit vector at the asserted
//! `r`, index them, and return everything inside the gate.
//!
//! Three design facts this module exists to encode, each measured rather than assumed:
//!
//! 1. ⭐ **De-parallax FIRST.** Gating raw topocentric directions leaves the parallactic reflex at
//!    the centre of the search disc -- at 40 AU opposition that reflex is 92.5"/day against a
//!    19.9"/day gate radius, i.e. 4.7x the radius itself. Projecting through the asserted `r`
//!    removes it exactly, which is what makes the disc small enough to be worth indexing.
//!
//! 2. ⭐ **Spherical, so no tangent plane.** A 3-D index over unit vectors is EXACT: angular
//!    separation `theta` maps to chord `2 sin(theta/2)` with no approximation. A gnomonic index
//!    stretches separations by ~1/cos^2 of the offset from its tangent point, so a fixed L2 radius
//!    keeps only 92.8% of true pairs at 10 deg and 75.9% at 20 deg -- a RECALL loss nothing
//!    downstream can observe. This is a correctness argument, not a preference.
//!
//! 3. 🔴 **The gate carries an astrometric floor.** `spherical_pair::gate_radius` is the pure
//!    dynamical radius, i.e. the `k = 0` case. `RESULT_c2_gate_ksigma.md` measured that the
//!    astrometric term is 59% of the k=3 gate at 40 AU and 83% at 60 AU, and that omitting it puts
//!    the tracklet density `rho` low by 2.6x (40 AU) to 6.3x (60 AU) -- which enters the pair count
//!    SQUARED. `gate_radius_astrometric` composes the two in quadrature; it is a new function
//!    rather than an edit because `gate_radius` is frozen by the `/3` fixture.
//!
//! Nothing here touches an existing pointdexter file.

use kiddo::SquaredEuclidean;
use kiddo::immutable::float::kdtree::ImmutableKdTree;
use nalgebra::Vector3;

use crate::spherical_pair::{gate_radius, range_quadratic};

/// 3-D tree over barycentric unit vectors. Content is the index into `BaryIndex::points`.
pub type BaryTree = ImmutableKdTree<f64, u32, 3, 32>;

/// Chord length subtending `theta` on the unit sphere. Exact, not a small-angle form.
#[inline]
pub fn chord_of_angle(theta: f64) -> f64 {
    2.0 * (0.5 * theta).sin()
}

/// Inverse of [`chord_of_angle`]. Clamped so floating point at the antipode cannot produce NaN.
///
/// 🔴 **Ill-conditioned near `pi`.** `d(asin x)/dx` diverges as `x -> 1`, so recovering an angle
/// from a chord loses precision as the separation approaches antipodal: the round trip holds to
/// ~1e-16 below 1 rad but only ~1.3e-15 at 170 deg, and it degrades as `1/sqrt(1-(c/2)^2)`.
/// Harmless for C2, which gates at arcseconds to a few degrees -- recorded because a caller
/// comparing near-antipodal separations through this function would be measuring the conditioning
/// rather than the sky.
#[inline]
pub fn angle_of_chord(chord: f64) -> f64 {
    2.0 * (0.5 * chord).clamp(-1.0, 1.0).asin()
}

/// Barycentric unit vector of a detection, given the asserted heliocentric distance.
///
/// `None` when the asserted `r` is unreachable along this line of sight, or the root is behind the
/// observer. Both are geometric rejections of the *hypothesis*, not failures of the detection:
/// a detection dropped here is simply not compatible with this shell.
pub fn bary_hat(
    observer: &Vector3<f64>,
    rho_hat: &Vector3<f64>,
    r_assumed: f64,
) -> Option<Vector3<f64>> {
    let p = range_quadratic(observer, rho_hat, r_assumed)?;
    let n = p.norm();
    if !(n > 0.0) {
        return None;
    }
    Some(p / n)
}

/// The gate as measured: the dynamical radius in quadrature with the k-sigma astrometric term.
///
/// `sigma_quad` is the two detections' astrometric sigmas summed in quadrature, in RADIANS, and
/// `k` is the sigma multiple. At `k = 0` this is exactly [`gate_radius`] -- the property the rho
/// measurement leaned on, so that only the gate changed between runs, and it is asserted in the
/// tests below.
///
/// 🔴 Evaluate at the cell's LOWER `r` edge, as [`gate_radius`] documents. The buffer is dominated
/// by the `r`-cell and the lower edge is exact where a centre value needs a 2.8x blind pad.
pub fn gate_radius_astrometric(
    r_lower_edge: f64,
    dt: f64,
    mu: f64,
    sigma_quad: f64,
    k: f64,
) -> f64 {
    let dynamical = gate_radius(r_lower_edge, dt, mu);
    let astrometric = k * sigma_quad;
    (dynamical * dynamical + astrometric * astrometric).sqrt()
}

/// How many detections a projection dropped, and why. Reported rather than swallowed: a shell that
/// silently loses half its detections looks like a shell with no objects in it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ProjectionReport {
    pub offered: usize,
    pub indexed: usize,
    pub unreachable: usize,
}

/// A searchable set of detections projected to one asserted distance.
pub struct BaryIndex {
    tree: BaryTree,
    /// Unit vectors actually indexed, parallel to `ids`.
    pub points: Vec<[f64; 3]>,
    /// Caller's identifier for each indexed point. Not the input position: unreachable detections
    /// are absent, so an index into the *input* would be wrong by exactly the dropped count.
    pub ids: Vec<u32>,
    pub report: ProjectionReport,
}

impl BaryIndex {
    /// Project every detection to the asserted `r` and index the ones that survive.
    ///
    /// `observers` may be per-detection (one entry each) or per-visit shared by the caller; this
    /// takes them per detection so a single index can span visits, which is what an intra-night
    /// tracklet build needs.
    pub fn build(
        observers: &[Vector3<f64>],
        rho_hats: &[Vector3<f64>],
        ids: &[u32],
        r_assumed: f64,
    ) -> Self {
        assert_eq!(observers.len(), rho_hats.len(), "observers/rho_hats length mismatch");
        assert_eq!(ids.len(), rho_hats.len(), "ids/rho_hats length mismatch");
        let mut points = Vec::with_capacity(rho_hats.len());
        let mut kept = Vec::with_capacity(ids.len());
        for i in 0..rho_hats.len() {
            if let Some(b) = bary_hat(&observers[i], &rho_hats[i], r_assumed) {
                points.push([b.x, b.y, b.z]);
                kept.push(ids[i]);
            }
        }
        let report = ProjectionReport {
            offered: rho_hats.len(),
            indexed: points.len(),
            unreachable: rho_hats.len() - points.len(),
        };
        let tree: BaryTree = ImmutableKdTree::new_from_slice(&points);
        BaryIndex { tree, points, ids: kept, report }
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Positions in `points` within `theta` (radians) of `query`.
    ///
    /// The query converts to a chord and kiddo compares SQUARED euclidean distance, so the
    /// threshold is `chord^2`. That composition is exact for unit vectors -- there is no pad and
    /// no tangent point.
    pub fn within_idx(&self, query: &Vector3<f64>, theta: f64) -> Vec<usize> {
        if theta <= 0.0 || self.points.is_empty() {
            return Vec::new();
        }
        // Beyond pi the chord stops increasing, so a naive chord threshold would silently shrink
        // the query. Nothing in C2 gates that wide, but a caller sizing a gate at an inner shell
        // can reach it, and a silently-narrowed query is a recall loss nothing downstream sees.
        if theta >= std::f64::consts::PI {
            return (0..self.points.len()).collect();
        }
        let c = chord_of_angle(theta);
        let q = [query.x, query.y, query.z];
        self.tree
            .within::<SquaredEuclidean>(&q, c * c)
            .into_iter()
            .map(|nn| nn.item as usize)
            .collect()
    }

    /// As [`Self::within_idx`], but returning the caller's ids.
    pub fn within(&self, query: &Vector3<f64>, theta: f64) -> Vec<u32> {
        self.within_idx(query, theta).into_iter().map(|i| self.ids[i]).collect()
    }
}

/// Result of pairing two projected sets inside one gate.
#[derive(Debug, Clone)]
pub struct PairSet {
    /// `(id in a, id in b)`.
    pub pairs: Vec<(u32, u32)>,
    pub n_a: usize,
    pub n_b: usize,
    pub theta: f64,
    /// `n_a * n_b * pi*theta^2 / omega` -- the `N * rho * pi theta^2` law the build plan asks
    /// Stage 2 to check measurement against. `f64::NAN` when no solid angle was supplied.
    pub predicted: f64,
}

impl PairSet {
    /// Measured / predicted. Above 1 means the field is clustered relative to uniform, which is
    /// the normal case for real sky and is exactly why the prediction is a check and not a
    /// substitute for the measurement.
    pub fn excess(&self) -> f64 {
        self.pairs.len() as f64 / self.predicted
    }
}

/// Every pair of (a, b) whose barycentric directions lie within `theta`.
///
/// `omega_sr` is the solid angle the two sets were drawn from, used only for the uniform-field
/// prediction; pass `None` to skip it.
pub fn pair_within_gate(
    a: &BaryIndex,
    b: &BaryIndex,
    theta: f64,
    omega_sr: Option<f64>,
) -> PairSet {
    let mut pairs = Vec::new();
    for (i, p) in a.points.iter().enumerate() {
        let q = Vector3::new(p[0], p[1], p[2]);
        for j in b.within_idx(&q, theta) {
            pairs.push((a.ids[i], b.ids[j]));
        }
    }
    let predicted = match omega_sr {
        Some(w) if w > 0.0 => {
            let disc = std::f64::consts::PI * theta * theta;
            a.len() as f64 * b.len() as f64 * disc / w
        }
        _ => f64::NAN,
    };
    PairSet { pairs, n_a: a.len(), n_b: b.len(), theta, predicted }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spacerocks::coordinates::Origin;

    /// Resolved the way every other C2 module resolves it, never a transcribed literal:
    /// the `/3` frame move is only safe if the port matches what it CALLS.
    fn mu() -> f64 {
        Origin::SSB.mu()
    }

    fn hat(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z).normalize()
    }

    /// Deterministic unit vectors in a cap of half-angle `cap` about `axis`, from a cheap LCG.
    /// A seeded generator rather than `rand`: this repo does not carry one and a test that needs
    /// a new dependency in someone else's Cargo.toml is a test that will not be run.
    fn cap_points(n: usize, axis: Vector3<f64>, cap: f64, seed: u64) -> Vec<Vector3<f64>> {
        let mut s = seed;
        let mut next = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 11) as f64) / ((1u64 << 53) as f64)
        };
        let e1 = if axis.z.abs() < 0.9 {
            axis.cross(&Vector3::z()).normalize()
        } else {
            axis.cross(&Vector3::x()).normalize()
        };
        let e2 = axis.cross(&e1);
        (0..n)
            .map(|_| {
                // uniform in solid angle within the cap
                let u = next();
                let cos_t = 1.0 - u * (1.0 - cap.cos());
                let sin_t = (1.0 - cos_t * cos_t).max(0.0).sqrt();
                let phi = 2.0 * std::f64::consts::PI * next();
                (cos_t * axis + sin_t * (phi.cos() * e1 + phi.sin() * e2)).normalize()
            })
            .collect()
    }

    #[test]
    fn chord_and_angle_invert_exactly_where_c2_operates() {
        // the claim the whole "no tangent plane" argument rests on, asserted over the range C2
        // actually gates: the widest anchor gate in the band is ~2 deg, and 1 rad is 29x that
        for i in 0..2000 {
            let theta = 1.0 * (i as f64) / 2000.0;
            let back = angle_of_chord(chord_of_angle(theta));
            assert!((back - theta).abs() < 1e-15, "theta {theta}: round trip gave {back}");
        }
    }

    #[test]
    fn the_chord_inverse_degrades_near_pi_as_its_conditioning_says_it_must() {
        // 🔴 Not a defect and not something to paper over with a loose tolerance everywhere: the
        // error is bounded by the derivative of asin, which diverges at the antipode. Asserting
        // the BOUND keeps the fact visible; a blanket 1e-14 would have hidden it and would also
        // have passed on a genuinely wrong implementation.
        let eps = f64::EPSILON;
        for i in 1..2000 {
            let theta = std::f64::consts::PI * (i as f64) / 2000.0;
            let c = chord_of_angle(theta);
            let back = angle_of_chord(c);
            // d(theta)/d(c) = 2 / sqrt(4 - c^2); one ulp of c maps to this much theta
            let condition = 2.0 / (4.0 - c * c).max(f64::MIN_POSITIVE).sqrt();
            let bound = 8.0 * eps * condition.max(1.0);
            assert!(
                (back - theta).abs() <= bound,
                "theta {theta}: error {:.3e} exceeds the conditioning bound {bound:.3e}",
                (back - theta).abs()
            );
        }
        // And the conditioning really does bite only at the far end. 🔴 Compared as the WORST
        // error over a range, not at two sampled angles: the round trip is exactly zero at most
        // individual points, so a two-point comparison measures which sample happened to be
        // unlucky. An earlier version of this assertion did exactly that and failed on 0 > 0.
        let worst = |lo: f64, hi: f64| {
            (0..500)
                .map(|i| {
                    let t = lo + (hi - lo) * (i as f64) / 500.0;
                    (angle_of_chord(chord_of_angle(t)) - t).abs()
                })
                .fold(0.0_f64, f64::max)
        };
        let near_pi = worst(3.0, std::f64::consts::PI * 0.9999);
        let small = worst(0.0, 0.1);
        assert!(
            near_pi > 10.0 * small.max(f64::EPSILON),
            "antipodal worst {near_pi:.3e} should dominate small-angle worst {small:.3e}"
        );
    }

    #[test]
    fn chord_matches_the_actual_separation_of_unit_vectors() {
        let a = hat(1.0, 0.0, 0.0);
        for i in 1..500 {
            let theta = std::f64::consts::PI * (i as f64) / 500.0;
            let b = hat(theta.cos(), theta.sin(), 0.0);
            let chord = (b - a).norm();
            assert!(
                (chord - chord_of_angle(theta)).abs() < 1e-15,
                "theta {theta}: chord {chord} vs {}",
                chord_of_angle(theta)
            );
        }
    }

    #[test]
    fn bary_hat_is_a_unit_vector_at_the_asserted_distance() {
        let obs = Vector3::new(0.3, -0.9, 0.0);
        let los = hat(0.2, 0.97, 0.05);
        let r = 42.0;
        let b = bary_hat(&obs, &los, r).expect("42 AU is reachable from 1 AU");
        assert!((b.norm() - 1.0).abs() < 1e-15);
        // and it really is the direction of a point at |p| = r
        let p = range_quadratic(&obs, &los, r).unwrap();
        assert!((p.norm() - r).abs() < 1e-9, "|p| = {} not {r}", p.norm());
        assert!((b - p.normalize()).norm() < 1e-15);
    }

    #[test]
    fn bary_hat_rejects_an_unreachable_distance() {
        // an observer 1 AU out cannot see anything at 0.1 AU along a line pointing outward
        let obs = Vector3::new(1.0, 0.0, 0.0);
        let los = hat(1.0, 0.0, 0.0);
        assert!(bary_hat(&obs, &los, 0.1).is_none());
        assert!(bary_hat(&obs, &los, -5.0).is_none());
    }

    #[test]
    fn the_astrometric_gate_reduces_to_the_dynamical_one_at_k_zero() {
        // 🔴 load-bearing: the rho measurement's "k = 0 reproduces section 3 exactly" depends on
        // this, so only the gate changed between those runs.
        for &r in &[40.0, 60.0, 150.0, 800.0] {
            for &dt in &[0.02, 0.5, 30.0, 365.25] {
                let bare = gate_radius(r, dt, mu());
                let with = gate_radius_astrometric(r, dt, mu(), 1.0e-6, 0.0);
                assert_eq!(bare, with, "r {r} dt {dt}");
            }
        }
    }

    #[test]
    fn the_astrometric_term_dominates_the_far_shells() {
        // The structural fact from RESULT_c2_gate_ksigma.md: the floor stops the gate shrinking
        // as r^-1.5. A test that only asserted "gate > 0" would pass on a missing term.
        let sigma = 0.142 / 206_264.806_247_096_36; // radians, the measured median
        let dt = 0.5;
        let inner = gate_radius_astrometric(40.0, dt, mu(), sigma, 3.0);
        let outer = gate_radius_astrometric(800.0, dt, mu(), sigma, 3.0);
        let inner_dyn = gate_radius(40.0, dt, mu());
        let outer_dyn = gate_radius(800.0, dt, mu());
        let astrom = 3.0 * sigma;

        // ✅ the dynamical radius reproduces the published table: 9.92" at 40 AU, 0.08" at 1000 AU
        let asec = |x: f64| x * 206_264.806_247_096_36;
        assert!((asec(inner_dyn) - 9.92).abs() < 0.01, "40 AU gave {}\"", asec(inner_dyn));
        assert!(
            (asec(gate_radius(1000.0, dt, mu())) - 0.0794).abs() < 0.001,
            "1000 AU gave {}\"",
            asec(gate_radius(1000.0, dt, mu()))
        );

        // 🔴 Assert WHICH TERM WINS, which is the structural claim, rather than a ratio picked by
        // eye. An earlier version of this test demanded outer/outer_dyn > 5.0; the derived value
        // is 3.97, so the threshold was measuring my guess and not the gate.
        assert!(astrom < 0.1 * inner_dyn, "at 40 AU the dynamics must dominate");
        assert!(astrom > 3.0 * outer_dyn, "at 800 AU the astrometry must dominate");
        assert!(inner / inner_dyn < 1.01, "inner: {} vs {}", inner, inner_dyn);
        let outer_ratio = outer / outer_dyn;
        assert!(
            (outer_ratio - 3.97).abs() < 0.05,
            "outer/outer_dyn is derived, not chosen: expected ~3.97, got {outer_ratio}"
        );
        // and the gate does NOT fall as fast as the dynamical radius: the floor is why per-node
        // tracklet cost stops falling with r, so neighbouring far shells select nearly one set
        assert!(inner / outer < inner_dyn / outer_dyn / 3.0);
    }

    #[test]
    fn within_agrees_with_brute_force_angular_separation() {
        // the load-bearing index test: kiddo's squared-chord query must select exactly the set a
        // direct angular comparison selects, with no pad either way
        let axis = hat(0.3, 0.5, 0.81);
        let pts = cap_points(400, axis, 0.05, 12345);
        let obs: Vec<Vector3<f64>> = (0..pts.len()).map(|_| Vector3::new(0.4, -0.85, 0.0)).collect();
        let ids: Vec<u32> = (0..pts.len() as u32).collect();
        let idx = BaryIndex::build(&obs, &pts, &ids, 45.0);
        assert_eq!(idx.len(), pts.len(), "all of these are reachable at 45 AU");

        for &theta in &[1e-5, 1e-4, 1e-3, 0.01, 0.03] {
            for probe in 0..25 {
                let q = Vector3::new(
                    idx.points[probe][0],
                    idx.points[probe][1],
                    idx.points[probe][2],
                );
                let mut got = idx.within_idx(&q, theta);
                got.sort_unstable();
                let mut want: Vec<usize> = (0..idx.len())
                    .filter(|&j| {
                        let p = Vector3::new(idx.points[j][0], idx.points[j][1], idx.points[j][2]);
                        let sep = angle_of_chord((p - q).norm());
                        sep <= theta
                    })
                    .collect();
                want.sort_unstable();
                assert_eq!(got, want, "theta {theta}, probe {probe}");
            }
        }
    }

    #[test]
    fn projection_reports_what_it_dropped() {
        let obs = vec![Vector3::new(1.0, 0.0, 0.0); 3];
        let los = vec![hat(1.0, 0.0, 0.0), hat(-1.0, 0.0, 0.0), hat(0.0, 1.0, 0.0)];
        let ids = vec![7u32, 8, 9];
        // 0.5 AU is behind or unreachable for the outward line of sight from 1 AU
        let idx = BaryIndex::build(&obs, &los, &ids, 0.5);
        assert_eq!(idx.report.offered, 3);
        assert_eq!(idx.report.indexed + idx.report.unreachable, 3);
        assert_eq!(idx.ids.len(), idx.len());
        // 🔴 ids must track the SURVIVORS: an index into the input would be off by the drop count
        for (slot, id) in idx.ids.iter().enumerate() {
            assert!(ids.contains(id), "slot {slot} carries an id that was never offered");
        }
    }

    #[test]
    fn pair_count_follows_n_rho_pi_theta_squared_on_a_uniform_field() {
        // the build plan's Stage 2 gate: measured vs predicted N*rho*pi*theta^2
        let axis = hat(0.0, 0.0, 1.0);
        let cap = 0.10_f64;
        let omega = 2.0 * std::f64::consts::PI * (1.0 - cap.cos()); // solid angle of the cap
        let obs: Vec<Vector3<f64>> = vec![Vector3::new(0.0, -1.0, 0.0); 900];
        let a_pts = cap_points(900, axis, cap, 999);
        let b_pts = cap_points(900, axis, cap, 4242);
        let ids: Vec<u32> = (0..900u32).collect();
        // 🔴 project both at the SAME r: the prediction is about a shared sphere, and two shells
        // would compare sets that do not live on one
        let a = BaryIndex::build(&obs, &a_pts, &ids, 500.0);
        let b = BaryIndex::build(&obs, &b_pts, &ids, 500.0);
        for &theta in &[0.002, 0.005, 0.01] {
            let ps = pair_within_gate(&a, &b, theta, Some(omega));
            let ratio = ps.excess();
            assert!(
                ratio > 0.75 && ratio < 1.35,
                "theta {theta}: measured {} predicted {:.1} (ratio {ratio:.3})",
                ps.pairs.len(),
                ps.predicted
            );
        }
    }

    #[test]
    fn pairs_are_exactly_the_brute_force_set() {
        let axis = hat(0.6, 0.0, 0.8);
        let obs: Vec<Vector3<f64>> = vec![Vector3::new(0.5, -0.8, 0.0); 120];
        let a_pts = cap_points(120, axis, 0.02, 5);
        let b_pts = cap_points(120, axis, 0.02, 6);
        let ids: Vec<u32> = (0..120u32).collect();
        let a = BaryIndex::build(&obs, &a_pts, &ids, 60.0);
        let b = BaryIndex::build(&obs, &b_pts, &ids, 60.0);
        let theta = 0.004;
        let mut got = pair_within_gate(&a, &b, theta, None).pairs;
        got.sort_unstable();
        let mut want = Vec::new();
        for i in 0..a.len() {
            for j in 0..b.len() {
                let p = Vector3::new(a.points[i][0], a.points[i][1], a.points[i][2]);
                let q = Vector3::new(b.points[j][0], b.points[j][1], b.points[j][2]);
                if angle_of_chord((p - q).norm()) <= theta {
                    want.push((a.ids[i], b.ids[j]));
                }
            }
        }
        want.sort_unstable();
        assert_eq!(got, want);
    }

    #[test]
    fn a_tangent_plane_index_would_lose_pairs_that_this_one_keeps() {
        // 🔴 The recall argument, as a test rather than a claim: gnomonic projection about a
        // tangent point stretches separations by ~1/cos^2 of the offset, so a fixed L2 radius in
        // the plane is a SHRINKING angular radius away from the tangent point. This asserts the
        // 3-D index does not have that failure -- the same separation is accepted regardless of
        // where on the sphere the pair sits.
        let theta = 0.01_f64;
        let obs = vec![Vector3::new(0.0, -1.0, 0.0); 2];
        for &off in &[0.0_f64, 0.1, 0.2, 0.35] {
            let a = hat(off.sin(), 0.0, off.cos());
            // b is exactly theta away from a, in the plane containing both and the pole
            let ang = off + theta;
            let b = hat(ang.sin(), 0.0, ang.cos());
            let ids = vec![0u32, 1];
            let idx = BaryIndex::build(&obs, &[a, b], &ids, 300.0);
            assert_eq!(idx.len(), 2);
            let q = Vector3::new(idx.points[0][0], idx.points[0][1], idx.points[0][2]);
            let hits = idx.within_idx(&q, theta * 1.000_001);
            assert_eq!(hits.len(), 2, "offset {off}: the pair must survive at any offset");
        }
    }
}
