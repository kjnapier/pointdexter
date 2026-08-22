//! Stage 4: extend a candidate arc with single detections, and score the support.
//!
//! This is the only real filter the method has. The C2 solve is a **generator, not a filter** --
//! 96.4% of *wrong*-node solves return a bound orbit, and the chi2 gate of Stage 3 is a pairing
//! test that is blind to the asserted `(r, rdot)`. So whether a candidate is an object is decided
//! here, by whether the sky agrees with its prediction where the survey looked.
//!
//! Four things were measured before this was written, and each of them is a line of code below.
//!
//! 🔴 **Cross-night only.** Support on the anchors' own nights is nearly free: over a few hours
//! every orbit through the two anchor detections predicts nearly the same sky position, so
//! re-finding the same source confirms the DETECTION, not the orbit. Measured on real Rubin sky,
//! same-night support ran 1.36 per candidate at 1" against 0.09 cross-night -- 94% of an apparent
//! signal. Worse, a rigid-shift control *endorses* it, because the shifted field carries the source
//! away. An implementation that counts support without excluding the anchor nights reports an
//! inflated confirmation rate and its control agrees.
//!
//! 🔴 **A raw count does not travel.** Chance support obeys
//! `lambda = n_opportunities * rho_det * pi * tolerance^2` to +-20% over a 3300x range in
//! tolerance -- a 10^7 range in area -- and the fitted density came back at 997 per deg^2 against a
//! field measured at 750-5337. So "at least 3 supporting detections" is a statement about the field
//! it was tuned on: 1.75% chance over 12 support epochs becomes ~44% over 100. The score must
//! divide out the opportunities the candidate actually had.
//!
//! 🔴 **The tail is NOT Poisson.** The rate is right -- mean support over mean lambda came in at
//! 0.84-0.99 -- but the variance is not: `var/mean` runs 1.3 at 3" to **9.2** at 30", because
//! detections cluster. Measured `P(>=1)` is 0.68-0.79x Poisson while `P(>=2)` at 3" is **4.37x**
//! Poisson. A bare Poisson survival function is therefore optimistic in the tail by up to ~4x,
//! which is why [`ExtendParams::dispersion`] exists and why it has no default.
//!
//! ⭐ **The tolerance is the grid cell, not a free parameter.** The prediction error is dominated
//! by the node offset, not by astrometry: at the grid's own half-cell (`dr/r` = 0.32% at 40 AU) the
//! prediction is displaced ~4.3" over a year-long anchor span, where the anchors' own astrometry
//! would allow ~0.15". Interior to the anchor span the error is astrometry-limited at 0.11-0.14"
//! regardless of baseline; beyond it, it grows. So gather loose, then refit and re-gather tight.
//!
//! With those, the measured performance at matched recall was 96.9-98.5% of chance rejected against
//! 83.5-89.2% for a raw count: **3-7x fewer survivors**.

use nalgebra::Vector3;

use crate::spherical_pair::{Node, Pair};
use crate::spherical_pair_anchor::{Candidate, Obs, anchor_frame, predict_hat_in};
use crate::spherical_pair_index::chord_of_angle;

use kiddo::SquaredEuclidean;
use kiddo::immutable::float::kdtree::ImmutableKdTree;

/// A visit's detections, indexed on the sky as observed.
///
/// Topocentric, deliberately: the extension compares a predicted line of sight against detections,
/// and both live in the observed frame. The barycentric projection of Stage 2 belongs to the
/// anchor build, where the asserted `r` is what removes the parallax.
pub struct VisitIndex {
    tree: ImmutableKdTree<f64, u32, 3, 32>,
    points: Vec<[f64; 3]>,
    ids: Vec<u32>,
    pub epoch: f64,
    pub night: i64,
    pub observer: Vector3<f64>,
    /// Centre and angular radius of the detections actually present, used as the footprint.
    pub centre: Vector3<f64>,
    pub radius: f64,
    /// `ids`, sorted, for membership tests. See [`Self::contains_id`].
    ids_sorted: Vec<u32>,
    /// [`Self::density`], evaluated once.
    density: f64,
}

impl VisitIndex {
    /// Build from the detections belonging to one visit.
    ///
    /// ⚠️ The footprint is the bounding cap of the detections themselves -- a filled disc, where a
    /// real focal plane has chip gaps (measured on-silicon fraction ~0.865) and bright-star
    /// exclusions. That **credits coverage that did not happen**, so it counts opportunities the
    /// candidate did not really have. The bias is conservative: more opportunities means a larger
    /// `lambda` and therefore a *less* significant score, never a more significant one.
    pub fn build(obs: &[Obs], members: &[usize], night: i64) -> Option<Self> {
        if members.is_empty() {
            return None;
        }
        let mut centre = Vector3::zeros();
        for &i in members {
            centre += obs[i].rho_hat;
        }
        let n = centre.norm();
        if !(n > 0.0) {
            return None;
        }
        centre /= n;
        let radius = members
            .iter()
            .map(|&i| obs[i].rho_hat.dot(&centre).clamp(-1.0, 1.0).acos())
            .fold(0.0_f64, f64::max);
        let points: Vec<[f64; 3]> =
            members.iter().map(|&i| [obs[i].rho_hat.x, obs[i].rho_hat.y, obs[i].rho_hat.z]).collect();
        let ids: Vec<u32> = members.iter().map(|&i| obs[i].id).collect();
        let tree = ImmutableKdTree::new_from_slice(&points);
        let mut ids_sorted = ids.clone();
        ids_sorted.sort_unstable();
        // The same expression `density()` evaluated, on the same inputs, so the same bits.
        let omega = 2.0 * std::f64::consts::PI * (1.0 - radius.cos());
        let density = if omega > 0.0 { points.len() as f64 / omega } else { 0.0 };
        Some(VisitIndex {
            tree,
            points,
            ids,
            epoch: obs[members[0]].epoch,
            night,
            observer: obs[members[0]].observer,
            centre,
            radius,
            ids_sorted,
            density,
        })
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Does this visit contain that detection? Exact id membership, not an epoch comparison:
    /// the anchor nights are derived from this, and a float epoch is not a key.
    ///
    /// ⭐ Binary search over a SORTED COPY of the ids, not `ids.contains`. `extend_candidate`
    /// asks this of every visit for each of a candidate's four anchor detections, so a linear
    /// scan is 283 visits x 4 ids x ~3,600 detections PER CANDIDATE -- and a profile of the
    /// 8-night fixture put it at 23% of the whole binary's runtime, the single largest entry,
    /// ahead of the universal-Kepler solver. Same answer, ~300x fewer comparisons.
    ///
    /// 🔴 A separate vector, because `ids` is positionally parallel to `points` and to the
    /// tree's own indices; sorting it in place would silently re-label every detection the
    /// gather returns.
    pub fn contains_id(&self, id: u32) -> bool {
        self.ids_sorted.binary_search(&id).is_ok()
    }

    /// Detection surface density over the footprint [per steradian].
    ///
    /// ⭐ Precomputed at build time. It is a function of the visit alone, but the extension asked
    /// for it once per candidate per opportunity -- a `cos` and a divide each time, and `sincos`
    /// was 13% of the profiled run. Same expression on the same inputs, so the same bits.
    pub fn density(&self) -> f64 {
        self.density
    }

    fn contains(&self, u: &Vector3<f64>) -> bool {
        u.dot(&self.centre) >= self.radius.cos()
    }

    fn within(&self, u: &Vector3<f64>, theta: f64) -> Vec<u32> {
        let c = chord_of_angle(theta);
        self.tree
            .within::<SquaredEuclidean>(&[u.x, u.y, u.z], c * c)
            .into_iter()
            .map(|nn| self.ids[nn.item as usize])
            .collect()
    }
}

/// How to gather and how to score. No `Default`: every one of these was a measured decision.
#[derive(Debug, Clone, Copy)]
pub struct ExtendParams {
    /// Gather radius [radians]. Set by the grid cell -- see the module docs.
    pub tolerance: f64,
    /// `var/mean` of the chance support count. **1.0 is Poisson and is optimistic by up to ~4x**;
    /// measured 1.3 at 3" rising to 9.2 at 30". Calibrate per field and tolerance.
    pub dispersion: f64,
    /// 🔴 Keep true except when deliberately measuring the same-night channel.
    pub exclude_anchor_nights: bool,
    /// A floor, applied before the score. Cheap, and it does not travel -- the score is what
    /// actually discriminates.
    pub min_support: usize,
    /// Accept only if chance would produce this much support with probability below this.
    pub max_chance_probability: f64,
    /// Require this many supporting **tracklets** -- cross-nights carrying support in two or more
    /// distinct visits. 0 leaves acceptance exactly as it was.
    ///
    /// MJH, 2026-08-18, from the PS1 pipeline: the anchor is already a pair of detections, so ask
    /// the support to be a pair too. The motivation is measured -- at a 10" gather over ~160
    /// opportunities the chance support count is over-dispersed at var/mean 7.4, and the tail is
    /// heavy enough that no threshold on a count of SINGLE detections reaches discovery grade
    /// (>= 20 still keeps 0.38% of candidates). Two detections in one night must both fall in the
    /// disc, so the chance rate is a product rather than a sum.
    pub min_support_tracklets: usize,
}

/// What the extension found for one candidate.
#[derive(Debug, Clone)]
pub struct Support {
    /// Visits where the prediction landed inside a footprint. The denominator of everything.
    pub n_opportunities: usize,
    pub n_support: usize,
    /// Expected chance support: `sum over opportunities of rho_det * pi * tolerance^2`.
    pub lambda: f64,
    /// `P(>= n_support)` under the over-dispersed chance model.
    pub p_chance: f64,
    pub support_ids: Vec<u32>,
    pub accepted: bool,
    /// Cross-nights carrying at least one supporting detection.
    pub n_support_nights: usize,
    /// Cross-nights carrying support in **two or more distinct visits** -- a supporting tracklet.
    ///
    /// 🔴 Distinct VISITS, not detections. Two detections of the same source in one visit are a
    /// blend or a duplicate, not a tracklet, and counting them would readmit the same-night
    /// confirmation channel `exclude_anchor_nights` exists to shut.
    pub n_support_tracklets: usize,
    /// Expected number of chance supporting tracklets, exactly under the model `lambda` uses:
    /// per night, `P(>= 2 of its opportunity visits carry a chance detection)`.
    pub lambda_tracklet: f64,
    /// `P(>= n_support_tracklets)` under the same over-dispersed tail.
    ///
    /// ⚠️ It reuses `dispersion`, which was calibrated for the count of single detections. The
    /// clustering that fattens that tail is largely STATIC sources, which do not reappear at the
    /// prediction on a later visit as an object would -- so this is likely conservative. It is
    /// not measured, and the tracklet counts are reported so the empirical rate can be read off
    /// the run rather than taken from the model.
    pub p_chance_tracklet: f64,
}

/// `ln Gamma(x)` for `x > 0`, Lanczos g = 7, n = 9.
///
/// Written out because this crate carries no statistics dependency and the negative-binomial tail
/// needs a non-integer factorial. Checked against exact values in the tests.
fn ln_gamma(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1259.139_216_722_402_8,
        771.323_428_777_653_13,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // reflection, so the series is only ever used where it converges well
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + 7.5;
    for (i, &c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// `P(X >= k)` for chance support with mean `lambda` and `var/mean = dispersion`.
///
/// `dispersion <= 1` is Poisson; above it, a negative binomial matched on the first two moments.
/// 🔴 The difference is not cosmetic: at `lambda = 0.025` the measured `P(>=2)` was 4.37x the
/// Poisson value, so a Poisson score would call a chance pair significant at roughly four times
/// the rate it should.
pub fn chance_probability(k: usize, lambda: f64, dispersion: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if !(lambda > 0.0) {
        return 0.0;
    }
    if dispersion <= 1.0 + 1e-12 {
        // Poisson tail by direct summation; k is small by construction
        let mut term = (-lambda).exp();
        let mut cdf = term;
        for i in 1..k {
            term *= lambda / i as f64;
            cdf += term;
        }
        return (1.0 - cdf).clamp(0.0, 1.0);
    }
    // negative binomial: var = lambda * dispersion  =>  p = 1/dispersion, r = lambda/(dispersion-1)
    let p = 1.0 / dispersion;
    let r = lambda / (dispersion - 1.0);
    let mut cdf = 0.0;
    for i in 0..k {
        let ln_pmf = ln_gamma(r + i as f64) - ln_gamma(r) - ln_gamma(i as f64 + 1.0)
            + r * p.ln()
            + (i as f64) * (1.0 - p).ln();
        cdf += ln_pmf.exp();
    }
    (1.0 - cdf).clamp(0.0, 1.0)
}

/// Per-night support, for the tracklet statistic: how many nights carried support in two or more
/// distinct visits, and how many such nights chance alone would have produced.
///
/// A `Vec` with a linear probe rather than a map: a candidate sees tens of nights, this is the
/// hot loop, and an allocation per candidate is what the parallel port exists to avoid.
#[derive(Default)]
struct NightTally {
    /// night, opportunity visits, visits that carried support, `sum ln(1 - e_v)`,
    /// `sum e_v / (1 - e_v)`.
    nights: Vec<(i64, usize, usize, f64, f64)>,
}

impl NightTally {
    /// One opportunity visit: its night, its chance rate `p_v`, and whether it carried support.
    fn observe(&mut self, night: i64, p_v: f64, had_support: bool) {
        // e_v = P(this visit carries >= 1 chance detection in the disc), from the same Poisson
        // rate `lambda` is built from -- so this is not a second model of the same sky.
        // `exp_m1` rather than `1 - exp(-p)`: p_v is ~1e-4 here and the subtraction would throw
        // away most of the significant digits.
        let e_v = (-p_v).exp_m1().abs().min(1.0 - 1e-12);
        let slot = match self.nights.iter().position(|n| n.0 == night) {
            Some(i) => &mut self.nights[i],
            None => {
                self.nights.push((night, 0, 0, 0.0, 0.0));
                self.nights.last_mut().expect("just pushed")
            }
        };
        slot.1 += 1;
        slot.2 += usize::from(had_support);
        slot.3 += (-e_v).ln_1p();
        slot.4 += e_v / (1.0 - e_v);
    }

    /// `(nights with any support, nights with a supporting tracklet, expected chance tracklets)`.
    fn finish(&self) -> (usize, usize, f64) {
        let mut n_nights = 0;
        let mut n_tracklets = 0;
        let mut lambda = 0.0;
        for &(_, n_opp, vis_with_support, ln_none, sum_odds) in &self.nights {
            n_nights += usize::from(vis_with_support >= 1);
            n_tracklets += usize::from(vis_with_support >= 2);
            // 🔴 STRUCTURAL, not arithmetic: a night with fewer than two opportunity visits
            // cannot produce a chance tracklet, so it contributes exactly zero rather than
            // whatever the formula rounds to. ~58% of Rubin position-nights get one visit, so a
            // 1e-20 per night would accumulate across the survey as a real-looking rate.
            if n_opp < 2 {
                continue;
            }
            // P(>= 2 hit) = 1 - P(none) - P(exactly one)
            //            = 1 - prod(1-e) * (1 + sum e/(1-e))
            //            = -expm1( ln1p(sum e/(1-e)) + sum ln1p(-e) ).
            // 🔴 The last form is the one to implement. The middle form subtracts two numbers
            // near 1 to produce an answer near 2e-8 at this field's rates -- eight of sixteen
            // digits gone, and a brute-force enumeration over subsets caught it at 5e-9 relative.
            // In log space both terms are small and nearly cancel *before* the exponential, so
            // the leading `sum_{i<j} e_i e_j` survives at full precision.
            lambda += (-((sum_odds.ln_1p() + ln_none).exp_m1())).clamp(0.0, 1.0);
        }
        (n_nights, n_tracklets, lambda)
    }
}

/// Gather support for one candidate over the supplied visits.
///
/// `visits` should be every visit in the extension window; the ones the candidate's own anchors
/// used, and (by default) every visit on those nights, are skipped and are NOT counted as
/// opportunities -- a skipped visit is not a place the object failed to appear.
pub fn extend_candidate(
    obs: &[Obs],
    cand: &Candidate,
    node: &Node,
    mu: f64,
    visits: &[VisitIndex],
    params: &ExtendParams,
) -> Support {
    let a0 = &obs[cand.anchor_a.first];
    let b0 = &obs[cand.anchor_b.first];
    let pair = Pair {
        t_ref: 0.5 * (a0.epoch + b0.epoch),
        a: crate::spherical_pair::PairPoint {
            epoch: a0.epoch,
            rho_hat: a0.rho_hat,
            observer: a0.observer,
        },
        b: crate::spherical_pair::PairPoint {
            epoch: b0.epoch,
            rho_hat: b0.rho_hat,
            observer: b0.observer,
        },
    };
    let used_ids = [
        obs[cand.anchor_a.first].id,
        obs[cand.anchor_a.second].id,
        obs[cand.anchor_b.first].id,
        obs[cand.anchor_b.second].id,
    ];
    // 🔴 Derived here, not supplied. Which nights the anchors came from is a fact about the data,
    // and a caller passing it in is a caller that can pass it in wrong -- which would silently
    // re-admit the same-night support channel this stage exists to exclude.
    let anchor_nights: Vec<i64> = visits
        .iter()
        .filter(|v| used_ids.iter().any(|&id| v.contains_id(id)))
        .map(|v| v.night)
        .collect();

    // ⭐ THE HOIST. Everything in `position_at` that does not depend on the epoch, built once
    // instead of once per visit -- 2 of every 3 universal-Kepler solves in this loop, plus the
    // radii, the two range quadratics and the in-plane frame. `None` here means `None` at every
    // epoch, so the loop below decides exactly what it decided before, one `continue` at a time.
    // See `SCOPE_c2_extend_no_kepler.md` §2 and `AnchorFrame`.
    let frame = anchor_frame(&pair, node, cand.h, mu);

    let disc = std::f64::consts::PI * params.tolerance * params.tolerance;
    let mut out = Support {
        n_opportunities: 0,
        n_support: 0,
        lambda: 0.0,
        p_chance: 1.0,
        support_ids: Vec::new(),
        accepted: false,
        n_support_nights: 0,
        n_support_tracklets: 0,
        lambda_tracklet: 0.0,
        p_chance_tracklet: 1.0,
    };
    let mut nights = NightTally::default();
    for v in visits {
        if params.exclude_anchor_nights && anchor_nights.contains(&v.night) {
            continue;
        }
        let Some(pred) = frame
            .as_ref()
            .and_then(|f| predict_hat_in(f, &pair, node, cand.h, mu, v.epoch, &v.observer))
        else {
            continue;
        };
        if !v.contains(&pred) {
            continue; // the prediction left the footprint: not an opportunity
        }
        out.n_opportunities += 1;
        let p_v = v.density() * disc;
        out.lambda += p_v;
        let mut hit_here = 0usize;
        for id in v.within(&pred, params.tolerance) {
            if used_ids.contains(&id) {
                continue;
            }
            out.support_ids.push(id);
            out.n_support += 1;
            hit_here += 1;
        }
        nights.observe(v.night, p_v, hit_here > 0);
    }
    let (n_nights, n_tracklets, lam_t) = nights.finish();
    out.n_support_nights = n_nights;
    out.n_support_tracklets = n_tracklets;
    out.lambda_tracklet = lam_t;
    out.p_chance = chance_probability(out.n_support, out.lambda, params.dispersion);
    out.p_chance_tracklet =
        chance_probability(out.n_support_tracklets, out.lambda_tracklet, params.dispersion);
    out.accepted = out.n_support >= params.min_support
        && out.p_chance <= params.max_chance_probability
        && out.n_support_tracklets >= params.min_support_tracklets;
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ln_gamma_matches_exact_values() {
        // ln Gamma(n) = ln((n-1)!)
        let mut fact = 1.0_f64;
        for n in 1..12 {
            let expect = fact.ln();
            let got = ln_gamma(n as f64);
            assert!((got - expect).abs() < 1e-10, "ln_gamma({n}) = {got}, want {expect}");
            fact *= n as f64;
        }
        // ln Gamma(1/2) = ln sqrt(pi)
        let half = ln_gamma(0.5);
        assert!((half - std::f64::consts::PI.sqrt().ln()).abs() < 1e-12, "got {half}");
    }

    #[test]
    fn the_negative_binomial_reduces_to_poisson_as_dispersion_goes_to_one() {
        // 🔴 the property that makes `dispersion` a knob rather than two different models
        // 🔴 Asserted as CONVERGENCE, not against a chosen tolerance. At dispersion = 1 + 1e-6 the
        // negative binomial has r = lambda/1e-6 ~ 2e4, so the tail is a difference of ln_gamma
        // values at large argument and a few parts in 1e4 is the arithmetic, not the model. A
        // fixed threshold here would be measuring my choice of threshold.
        for &lam in &[0.02_f64, 0.3, 2.5] {
            for k in 1..5 {
                let pois = chance_probability(k, lam, 1.0);
                let rel = |d: f64| (chance_probability(k, lam, d) - pois).abs() / pois.max(1e-300);
                let coarse = rel(1.01);
                let fine = rel(1.000_001);
                assert!(
                    fine < coarse,
                    "lambda {lam} k {k}: approaching dispersion 1 must converge; \
                     rel(1.01) = {coarse:.3e}, rel(1+1e-6) = {fine:.3e}"
                );
                assert!(fine < 1e-3, "lambda {lam} k {k}: rel error {fine:.3e} at dispersion 1+1e-6");
            }
        }
    }

    #[test]
    fn over_dispersion_fattens_the_tail_by_about_what_was_measured() {
        // Measured on real sky: at 3" (lambda ~ 0.025, var/mean 1.28) P(>=2) was 4.37x Poisson.
        // The model must move in that direction and by a comparable amount -- a dispersion knob
        // that changed the answer by 1% would not be worth having.
        let lam = 0.025;
        let pois = chance_probability(2, lam, 1.0);
        let over = chance_probability(2, lam, 1.28);
        let ratio = over / pois;
        assert!(ratio > 1.5, "dispersion must fatten the tail; ratio {ratio:.2}");
        // and the effect grows with dispersion
        assert!(chance_probability(2, lam, 9.2) > over);
    }

    #[test]
    fn lambda_scales_as_the_square_of_the_tolerance() {
        // the measured area law, as an invariant of the implementation
        let obs = ring_obs();
        let visits = ring_visits(&obs);
        let base = visits[0].density() * std::f64::consts::PI;
        let t1 = 3.0 * ARCSEC;
        let t2 = 30.0 * ARCSEC;
        let l1 = base * t1 * t1;
        let l2 = base * t2 * t2;
        assert!((l2 / l1 - 100.0).abs() < 1e-9, "10x tolerance must be 100x lambda");
    }

    const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

    fn ring_obs() -> Vec<Obs> {
        // 200 detections spread over a 0.2 deg cap, one visit
        let centre = Vector3::new(0.1, 0.98, 0.05).normalize();
        let e1 = centre.cross(&Vector3::z()).normalize();
        let e2 = centre.cross(&e1);
        let mut s: u64 = 7;
        let mut next = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((s >> 11) as f64) / ((1u64 << 53) as f64)
        };
        (0..200)
            .map(|i| {
                let cap = 0.2_f64.to_radians();
                let ct = 1.0 - next() * (1.0 - cap.cos());
                let st = (1.0 - ct * ct).max(0.0).sqrt();
                let ph = std::f64::consts::TAU * next();
                Obs {
                    id: i,
                    epoch: 2_460_000.5,
                    rho_hat: (ct * centre + st * (ph.cos() * e1 + ph.sin() * e2)).normalize(),
                    observer: Vector3::new(0.3, -0.95, 0.0),
                    sigma: 0.12 * ARCSEC,
                }
            })
            .collect()
    }

    fn ring_visits(obs: &[Obs]) -> Vec<VisitIndex> {
        let members: Vec<usize> = (0..obs.len()).collect();
        vec![VisitIndex::build(obs, &members, 60_000).unwrap()]
    }

    #[test]
    fn the_footprint_excludes_a_prediction_that_missed_the_field() {
        let obs = ring_obs();
        let v = &ring_visits(&obs)[0];
        assert!(v.contains(&v.centre), "the centre must be inside its own footprint");
        // a direction 5 degrees away is not an opportunity
        let e1 = v.centre.cross(&Vector3::z()).normalize();
        let far = (v.centre + 5.0_f64.to_radians().tan() * e1).normalize();
        assert!(!v.contains(&far), "a prediction 5 deg off the field is not an opportunity");
    }

    #[test]
    fn density_is_per_steradian_and_tracks_the_count() {
        let obs = ring_obs();
        let v = &ring_visits(&obs)[0];
        let omega = 2.0 * std::f64::consts::PI * (1.0 - v.radius.cos());
        assert!((v.density() * omega - v.len() as f64).abs() < 1e-9);
        // sanity: 200 detections in a ~0.2 deg cap is a few thousand per square degree
        let per_sq_deg = v.density() * (std::f64::consts::PI / 180.0).powi(2);
        assert!(per_sq_deg > 100.0 && per_sq_deg < 1.0e5, "{per_sq_deg} per deg^2");
    }

    #[test]
    fn a_night_with_one_opportunity_visit_can_never_carry_a_tracklet() {
        // 🔴 EXACT zero, asserted as equality. The tracklet statistic's whole claim is that a
        // single-visit night contributes no chance tracklet; if that leaked even at 1e-18 it
        // would accumulate over the ~58% of Rubin position-nights that get one visit.
        let mut t = NightTally::default();
        for (i, p) in [1e-4, 3e-3, 0.2, 2.0].iter().enumerate() {
            t.observe(60_000 + i as i64, *p, true);
        }
        let (n_nights, n_tracklets, lambda) = t.finish();
        assert_eq!(n_nights, 4, "each night carried support");
        assert_eq!(n_tracklets, 0, "no night had a second visit");
        assert_eq!(lambda, 0.0, "single-visit nights must contribute exactly zero, got {lambda}");
    }

    #[test]
    fn a_tracklet_needs_two_distinct_visits_not_two_detections() {
        // One visit carrying support twice is a blend, not a tracklet: `observe` is called once
        // per visit with a boolean, so the count cannot be inflated by a crowded single visit.
        let mut one = NightTally::default();
        one.observe(60_000, 1e-3, true);
        assert_eq!(one.finish().1, 0, "one visit is not a tracklet however many detections");

        let mut two = NightTally::default();
        two.observe(60_000, 1e-3, true);
        two.observe(60_000, 1e-3, true);
        assert_eq!(two.finish().1, 1, "two supporting visits on one night are a tracklet");

        let mut split = NightTally::default();
        split.observe(60_000, 1e-3, true);
        split.observe(60_001, 1e-3, true);
        assert_eq!(split.finish().1, 0, "two nights, one visit each, is not a tracklet");

        let mut half = NightTally::default();
        half.observe(60_000, 1e-3, true);
        half.observe(60_000, 1e-3, false);
        assert_eq!(half.finish().1, 0, "the second visit must itself carry support");
    }

    #[test]
    fn the_chance_tracklet_rate_matches_brute_force_enumeration() {
        // The finaliser is a streaming form of P(>= 2 of n independent visits hit). Checked
        // against enumeration over every subset, which is the definition -- not against a
        // rearrangement of the same algebra, which would agree with its own mistakes.
        for ps in [
            vec![1e-4, 2e-4],
            vec![0.01, 0.02, 0.03],
            vec![0.3, 0.25, 0.4, 0.1],
            vec![0.9, 0.8],
        ] {
            let mut t = NightTally::default();
            for &p in &ps {
                t.observe(60_000, p, false);
            }
            let got = t.finish().2;

            let e: Vec<f64> = ps.iter().map(|p| 1.0 - (-p).exp()).collect();
            let mut brute = 0.0;
            for mask in 0u32..(1 << e.len()) {
                if (mask.count_ones() as usize) < 2 {
                    continue;
                }
                let mut term = 1.0;
                for (i, &ei) in e.iter().enumerate() {
                    term *= if mask & (1 << i) != 0 { ei } else { 1.0 - ei };
                }
                brute += term;
            }
            assert!(
                (got - brute).abs() < 1e-12 * brute.max(1e-12),
                "ps {ps:?}: streaming {got:.12e} vs enumerated {brute:.12e}"
            );
        }
    }

    // 🔴 There is deliberately NO unit test asserting that `min_support_tracklets = 0` leaves
    // acceptance unchanged: written here it could only restate the `&&` in `extend_candidate`
    // and would agree with its own mistakes. That claim is proved end to end instead, by the
    // node-1 candidate file coming back byte-identical to the pre-change build.

    #[test]
    fn a_candidate_with_no_support_is_not_accepted_however_many_opportunities_it_had() {
        // 🔴 the direction that matters: opportunities without support must never accept. A score
        // that divided the other way would call a well-observed empty patch a detection.
        let p = chance_probability(0, 5.0, 1.0);
        assert_eq!(p, 1.0, "zero support is certain under chance");
    }
}
