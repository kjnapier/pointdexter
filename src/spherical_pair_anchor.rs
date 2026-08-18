//! Stage 3: anchors -> anchor pairs -> solve -> vetted candidate states.
//!
//! MJH's architecture (2026-08-15): **each anchor is a pair of detections on one night**, and the
//! support along the arc is single detections (stage 4). So this stage is:
//!
//! ```text
//!   detections --(barycentric gate, intra-night)--> anchors
//!   anchors x anchors --(barycentric gate, cross-night)--> anchor pairs
//!   anchor pair --(solve_h, state_from_solution)--> candidate state
//!   candidate --(chi2 on the two UNUSED detections)--> vetted candidate
//! ```
//!
//! Two facts make the last two steps what they are, and both were measured rather than assumed:
//!
//! 🔴 **The solve is NOT symmetric in its endpoints**, and the two implementations fail
//! DIFFERENTLY -- which is worth stating precisely, because the Python measurement's description
//! does not carry over.
//!
//! * In the Python prototype a reversed pair returned a **retrograde orbit that solved silently**,
//!   bound, raising nothing, with the residual inflated ~20x (median chi2 24.6 against 1.32 on the
//!   same night pair). A purity failure, and a flattering one: it makes the gate look stronger.
//! * In this port a reversed pair is **rejected outright** (`NoBoundRoot`). Stage 1 accumulates a
//!   SIGNED swept angle and returns `None` past `pi` rather than folding, so the backwards sweep
//!   cannot produce a root at all. Verified in the tests below.
//!
//! ⇒ here the risk is **recall, not purity**: an unordered anchor pair does not lie to you, it
//! disappears, and a real object goes with it. The ordering is therefore enforced once, in
//! [`solve_anchor_pair`], so no caller can drop candidates by handing its anchors over backwards.
//!
//! ⭐ **The chi2 has dof 4 and its weights are levered.** An anchor pair consumes 4 numbers (two
//! unit vectors) with `(r, rdot)` asserted, so the two detections NOT used by the solve supply
//! 4 residual coordinates and 4 - 0 = 4 degrees of freedom; the 95% cut is 9.488. The solve pins
//! both anchors exactly, so their astrometric errors ride into the prediction: the near anchor with
//! coefficient ~1 over the short intra-anchor step, the far one levered down by
//! `(anchor baseline / pair baseline)`. Getting that charge wrong does not land the injection on
//! chi2_4, which is what the self-test here checks.
//!
//! Measured performance of the gate, for whoever sets the threshold: it rejects **0.88** of chance
//! pairs at Rubin's typical ~0.6 h anchor baseline and **0.997** at a 6 h baseline -- but the survey
//! supplies a 6 h baseline on only 2.3% of tracklet-forming position-nights, so 0.88 is the
//! operating point, not 0.997.

use nalgebra::Vector3;

use crate::spherical_pair::{
    CanonicalState, GuidedCheck, Node, Pair, PairPoint, ScanCensus, Solution, canonical_step,
    hypot2, range_quadratic, swept_angle,
};
use crate::spherical_pair_index::{BaryIndex, gate_radius_astrometric};

/// 95% of chi2 with 4 degrees of freedom: the measured pre-extension gate.
pub const CHI2_CUT_DOF4: f64 = 9.488;

/// One detection, in the only terms this stage needs.
///
/// Deliberately not `crate::detection::Detection`: that type is tangent-plane oriented and owned by
/// the existing pipeline, and Stage 3 must not acquire a dependency on a schema it does not
/// control. The adapter from whatever the loader produces belongs at the binary's edge.
#[derive(Debug, Clone, Copy)]
pub struct Obs {
    pub id: u32,
    /// TDB Julian date.
    pub epoch: f64,
    /// Topocentric unit vector to the source.
    pub rho_hat: Vector3<f64>,
    /// Barycentric observer position [AU].
    pub observer: Vector3<f64>,
    /// Astrometric uncertainty [radians], per coordinate.
    pub sigma: f64,
}

impl Obs {
    fn point(&self) -> PairPoint {
        PairPoint { epoch: self.epoch, rho_hat: self.rho_hat, observer: self.observer }
    }
}

/// An anchor: two detections on one night, ordered in time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Anchor {
    /// Index into the caller's `&[Obs]`, earlier epoch.
    pub first: usize,
    /// Index into the caller's `&[Obs]`, later epoch.
    pub second: usize,
}

/// Build anchors by pairing two visits through the barycentric gate at the asserted distance.
///
/// `visit_a` and `visit_b` are index lists into `obs`; they are expected to be two visits of one
/// night. ⭐ The gate is applied in BARYCENTRIC unit vectors at `r_lower_edge`, which removes the
/// parallactic reflex before the comparison -- gating raw topocentric rate leaves that reflex at
/// 4.7x the gate radius at 40 AU, so the disc would be centred on the wrong point.
///
/// The KD query uses a population-wide `sigma_cap` so one tree serves every pair; the exact
/// per-pair gate is then re-applied, so the cap only ever admits extra candidates for the exact
/// test to reject. 🔴 That order matters: capping *after* would silently drop true pairs.
/// 🔴 TWO r's, and they are not the same number (split 2026-08-16).
///
/// * `r_assert` is the HYPOTHESIS: the distance at which the detections are placed to make them
///   barycentric. It is the node's own r.
/// * `r_gate_lower` is a BOUND over the node's cell: `gate_radius ~ r^-3/2`, so the pairing
///   radius must be evaluated at the cell's SMALLEST r, whose members move fastest. That is what
///   `gate_radius`'s "the lower edge is exact where a centre value needs a 2.8x blind pad" means.
///
/// One parameter served both until the sized ladder was wired, and with hand-placed nodes that
/// was invisible -- an explicit node has no cell, so the two coincide. The ladder emits
/// `r = 1/gamma` at the cell's LOWER gamma edge, i.e. its LARGEST r, which is the value that
/// makes the gate smallest. At a 56%-wide cell that under-sizes the pairing radius by ~3x and
/// the object's own anchors are never proposed. Under-admission is invisible downstream by
/// construction: nothing can observe a pair that was not offered.
pub fn build_anchors(
    obs: &[Obs],
    visit_a: &[usize],
    visit_b: &[usize],
    r_assert: f64,
    r_gate_lower: f64,
    k: f64,
    mu: f64,
    sigma_cap: f64,
) -> Vec<Anchor> {
    if visit_a.is_empty() || visit_b.is_empty() {
        return Vec::new();
    }
    let dt = (obs[visit_b[0]].epoch - obs[visit_a[0]].epoch).abs();
    if !(dt > 0.0) {
        return Vec::new();
    }
    let pack = |v: &[usize]| {
        let o: Vec<Vector3<f64>> = v.iter().map(|&i| obs[i].observer).collect();
        let u: Vec<Vector3<f64>> = v.iter().map(|&i| obs[i].rho_hat).collect();
        let ids: Vec<u32> = (0..v.len() as u32).collect();
        BaryIndex::build(&o, &u, &ids, r_assert)
    };
    let ia = pack(visit_a);
    let ib = pack(visit_b);
    if ia.is_empty() || ib.is_empty() {
        return Vec::new();
    }
    let cap = gate_radius_astrometric(r_gate_lower, dt, mu, sigma_cap, k);

    let mut out = Vec::new();
    for (slot_a, p) in ia.points.iter().enumerate() {
        let q = Vector3::new(p[0], p[1], p[2]);
        for slot_b in ib.within_idx(&q, cap) {
            let i = visit_a[ia.ids[slot_a] as usize];
            let j = visit_b[ib.ids[slot_b] as usize];
            // the exact per-pair gate: this pair's own sigmas, in quadrature
            let sq = (obs[i].sigma * obs[i].sigma + obs[j].sigma * obs[j].sigma).sqrt();
            let exact = gate_radius_astrometric(r_gate_lower, dt, mu, sq, k);
            let b = Vector3::new(ib.points[slot_b][0], ib.points[slot_b][1], ib.points[slot_b][2]);
            let sep = crate::spherical_pair_index::angle_of_chord((q - b).norm());
            if sep > exact {
                continue;
            }
            // 🔴 ordered in time, always: everything downstream reads the direction of motion off
            // the endpoint order, and an unordered anchor poisons the solve silently
            out.push(if obs[i].epoch <= obs[j].epoch {
                Anchor { first: i, second: j }
            } else {
                Anchor { first: j, second: i }
            });
        }
    }
    out
}

/// A solved, vetted anchor pair.
#[derive(Debug, Clone, Copy)]
pub struct Candidate {
    /// Barycentric state at [`Self::epoch`].
    pub state: [f64; 6],
    /// Epoch of the state: the EARLIER anchor's first detection.
    pub epoch: f64,
    /// Angular momentum the solve converged on.
    pub h: f64,
    /// Residual of the two detections the solve did not use. dof 4.
    pub chi2: f64,
    pub anchor_a: Anchor,
    pub anchor_b: Anchor,
}

/// Why an anchor pair produced no candidate. A bare `None` conflates facts that mean different
/// things -- the same objection Stage 1's [`Solution`] answers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Reject {
    /// The asserted node admits no bound orbit at all.
    UnboundNode,
    /// Bound node, but only a hyperbolic orbit threads both anchors. A real, free rejection.
    NoBoundRoot,
    /// The solve converged but the geometry could not be rebuilt.
    NoState,
    /// A residual epoch could not be predicted.
    NoPrediction,
    /// Solved and predicted, but failed the chi2 gate. Carries the value so a caller can retune
    /// without re-running the solve.
    Chi2 { chi2: f64 },
}

/// Barycentric position at `t`, for an orbit already pinned by `pair` and `h`.
///
/// ⭐ Built from Stage 1's own canonical step rather than from a second propagator: that routine is
/// gated to 1.2e-16 against the frozen oracle, and a parallel implementation is a thing that can
/// drift from it. The plane and the reference direction come from the two anchors, exactly as
/// `state_from_solution` derives them, so the prediction is the same orbit continued.
pub fn position_at(pair: &Pair, node: &Node, h: f64, mu: f64, t: f64) -> Option<Vector3<f64>> {
    let frame = anchor_frame(pair, node, h, mu)?;
    position_at_in(&frame, pair, node, h, mu, t)
}

/// The part of [`position_at`] that does not depend on `t`.
///
/// ⭐ Hoisted because Stage 4 asks for a position at EVERY visit epoch in the window, and two of
/// `position_at`'s three `canonical_step` calls have no `t` in them: over this run's 283-visit
/// window, 566 of the 849 universal-Kepler solves per candidate were recomputing two fixed
/// values. See `SCOPE_c2_extend_no_kepler.md` §2.
///
/// 🔴 The move is BIT-EXACT and is meant to stay that way: same function, same arguments, in the
/// same order, nothing reassociated. Its regression test is a candidate CSV that compares BYTE
/// FOR BYTE with one from before the hoist -- not an approximate tolerance, which would let a
/// real numerical change hide inside an "acceptable" difference.
#[derive(Debug, Clone, Copy)]
pub struct AnchorFrame {
    /// Canonical state at anchor A's epoch. `swept_angle` measures the sweep from here.
    s_a: CanonicalState,
    /// In-plane basis built from the two anchor positions: `a_hat` along anchor A.
    a_hat: Vector3<f64>,
    t_hat: Vector3<f64>,
}

/// Build the `t`-invariant frame once, for reuse at many `t`.
///
/// Returns `None` in exactly the cases [`position_at`] returned `None` for a `t`-invariant
/// reason; the `t`-dependent ones stay in [`position_at_in`]. A caller that gets `None` here
/// gets `None` at every `t`, which is why the extension can hoist the call out of its loop
/// without changing what any single visit decides.
pub fn anchor_frame(pair: &Pair, node: &Node, h: f64, mu: f64) -> Option<AnchorFrame> {
    let s_a: CanonicalState = canonical_step(node, h, pair.a.epoch - pair.t_ref, mu)?;
    let s_b: CanonicalState = canonical_step(node, h, pair.b.epoch - pair.t_ref, mu)?;

    let r_a = hypot2(s_a.pos);
    let r_b = hypot2(s_b.pos);
    if !(r_a > 0.0) || !(r_b > 0.0) {
        return None;
    }
    let p_a = range_quadratic(&pair.a.observer, &pair.a.rho_hat, r_a)?;
    let p_b = range_quadratic(&pair.b.observer, &pair.b.rho_hat, r_b)?;

    let n = p_a.cross(&p_b);
    let nn = n.norm();
    if nn == 0.0 {
        return None;
    }
    let n_hat = n / nn;
    let a_hat = p_a / p_a.norm();
    let t_hat = n_hat.cross(&a_hat);
    Some(AnchorFrame { s_a, a_hat, t_hat })
}

/// [`position_at`] with the `t`-invariant half already built.
pub fn position_at_in(
    frame: &AnchorFrame,
    pair: &Pair,
    node: &Node,
    h: f64,
    mu: f64,
    t: f64,
) -> Option<Vector3<f64>> {
    let s_t: CanonicalState = canonical_step(node, h, t - pair.t_ref, mu)?;
    let r_t = hypot2(s_t.pos);
    if !(r_t > 0.0) {
        return None;
    }
    // sweep from anchor A to t, in the same sense the solve used
    let d_nu = swept_angle(node, h, mu, &frame.s_a, &s_t, t - pair.a.epoch);
    Some(r_t * (d_nu.cos() * frame.a_hat + d_nu.sin() * frame.t_hat))
}

/// Topocentric unit vector predicted at `t`, as seen from `observer`.
pub fn predict_hat(
    pair: &Pair,
    node: &Node,
    h: f64,
    mu: f64,
    t: f64,
    observer: &Vector3<f64>,
) -> Option<Vector3<f64>> {
    let p = position_at(pair, node, h, mu, t)?;
    topocentric_hat(&p, observer)
}

/// [`predict_hat`] with the `t`-invariant half already built.
pub fn predict_hat_in(
    frame: &AnchorFrame,
    pair: &Pair,
    node: &Node,
    h: f64,
    mu: f64,
    t: f64,
    observer: &Vector3<f64>,
) -> Option<Vector3<f64>> {
    let p = position_at_in(frame, pair, node, h, mu, t)?;
    topocentric_hat(&p, observer)
}

/// Reduce a barycentric position to the unit vector an observer sees.
fn topocentric_hat(p: &Vector3<f64>, observer: &Vector3<f64>) -> Option<Vector3<f64>> {
    let geo = p - observer;
    let n = geo.norm();
    if !(n > 0.0) {
        return None;
    }
    Some(geo / n)
}

/// Solve one anchor pair and vet it.
///
/// 🔴 The two anchors are ordered by epoch here and nowhere else, so no caller can get it wrong.
pub fn solve_anchor_pair(
    obs: &[Obs],
    anchor_a: Anchor,
    anchor_b: Anchor,
    node: &Node,
    mu: f64,
    chi2_cut: f64,
) -> Result<Candidate, Reject> {
    // 🔴 THE ORDERING GUARD. See the module docs: reversed endpoints return a retrograde orbit
    // that solves silently with a ~20x residual.
    let (anchor_a, anchor_b) = if obs[anchor_a.first].epoch <= obs[anchor_b.first].epoch {
        (anchor_a, anchor_b)
    } else {
        (anchor_b, anchor_a)
    };
    let a0 = &obs[anchor_a.first];
    let a1 = &obs[anchor_a.second];
    let b0 = &obs[anchor_b.first];
    let b1 = &obs[anchor_b.second];

    let t_ref = 0.5 * (a0.epoch + b0.epoch);
    let pair = Pair { t_ref, a: a0.point(), b: b0.point() };
    let h = match pair.solve_h(node, mu) {
        Solution::UnboundNode => return Err(Reject::UnboundNode),
        Solution::NoBoundRoot { .. } => return Err(Reject::NoBoundRoot),
        Solution::Bound { h, .. } => h,
    };
    let state = pair.state_from_solution(h, node, mu).ok_or(Reject::NoState)?;

    let dt_pair = (b0.epoch - a0.epoch).abs().max(f64::MIN_POSITIVE);
    let mut chi2 = 0.0;
    // The two detections the solve did NOT consume: one per anchor. 4 coordinates, 0 fitted.
    for (res, near, far) in [(a1, a0, b0), (b1, b0, a0)] {
        let pred = predict_hat(&pair, node, h, mu, res.epoch, &res.observer)
            .ok_or(Reject::NoPrediction)?;
        // tangential part of the residual: the radial component is not observed
        let d = res.rho_hat - (res.rho_hat.dot(&pred)) * pred;
        // the near anchor's error rides in with coefficient ~1 over the short intra-anchor step;
        // the far one is levered down by (anchor baseline / pair baseline)
        let lever = (res.epoch - near.epoch).abs() / dt_pair;
        let s2 = res.sigma * res.sigma
            + near.sigma * near.sigma
            + (lever * far.sigma) * (lever * far.sigma);
        if !(s2 > 0.0) {
            return Err(Reject::NoPrediction);
        }
        chi2 += d.dot(&d) / s2;
    }
    if !(chi2 < chi2_cut) {
        return Err(Reject::Chi2 { chi2 });
    }
    Ok(Candidate { state, epoch: a0.epoch, h, chi2, anchor_a, anchor_b })
}

/// Pair every A-anchor against every B-anchor inside the gate, solve, and keep what passes.
///
/// Returns the candidates and a tally of why the rest were rejected -- a stage that reports only
/// its survivors cannot tell "the node is empty" from "the gate ate everything".
pub fn anchor_pairs_and_solve(
    obs: &[Obs],
    anchors_a: &[Anchor],
    anchors_b: &[Anchor],
    node: &Node,
    k: f64,
    mu: f64,
    sigma_cap: f64,
    chi2_cut: f64,
) -> (Vec<Candidate>, RejectTally) {
    let mut tally = RejectTally::default();
    let mut out = Vec::new();
    if anchors_a.is_empty() || anchors_b.is_empty() {
        return (out, tally);
    }
    // gate on the widest anchor separation the two sets can produce: it must not reject a pair for
    // sitting at the far end of its own night
    let dt_max = anchors_a
        .iter()
        .flat_map(|a| anchors_b.iter().map(move |b| (obs[b.first].epoch - obs[a.first].epoch).abs()))
        .fold(0.0_f64, f64::max);
    let cap = gate_radius_astrometric(node.r, dt_max, mu, sigma_cap, k);

    let pack = |set: &[Anchor]| {
        let o: Vec<Vector3<f64>> = set.iter().map(|a| obs[a.first].observer).collect();
        let u: Vec<Vector3<f64>> = set.iter().map(|a| obs[a.first].rho_hat).collect();
        let ids: Vec<u32> = (0..set.len() as u32).collect();
        BaryIndex::build(&o, &u, &ids, node.r)
    };
    let ia = pack(anchors_a);
    let ib = pack(anchors_b);

    for (slot, p) in ia.points.iter().enumerate() {
        let q = Vector3::new(p[0], p[1], p[2]);
        for hit in ib.within_idx(&q, cap) {
            let a = anchors_a[ia.ids[slot] as usize];
            let b = anchors_b[ib.ids[hit] as usize];
            match solve_anchor_pair(obs, a, b, node, mu, chi2_cut) {
                Ok(c) => out.push(c),
                Err(e) => tally.count(e),
            }
        }
    }
    (out, tally)
}

/// Census variant of [`anchor_pairs_and_solve`]: same gate, same pairs, same `Pair`, but walk
/// the whole `F(h)` scan and count brackets instead of returning at the first.
///
/// **DIAGNOSTIC ONLY.** It answers whether an initial-guess or bracket-lookup scheme could
/// safely replace the blind 80-point scan -- which is true only if `F` is single-rooted.
///
/// 🔴 The gate and the endpoint ordering are duplicated from `anchor_pairs_and_solve` /
/// `solve_anchor_pair` rather than shared, because those functions solve and this one must not.
/// The duplication is the risk: a census taken over a DIFFERENT pair set measures nothing about
/// the pairs the search actually solves. Both copies must stay in step, and the ordering guard
/// below is not optional -- reversed endpoints return a retrograde orbit that solves silently.
pub fn anchor_pairs_and_census(
    obs: &[Obs],
    anchors_a: &[Anchor],
    anchors_b: &[Anchor],
    node: &Node,
    k: f64,
    mu: f64,
    sigma_cap: f64,
) -> (Vec<(f64, ScanCensus, GuidedCheck)>, Vec<Pair>) {
    let mut out = Vec::new();
    // 🔴 The pairs where `solve_h` and `solve_h_guided` reached DIFFERENT outcomes, carried out
    // of the SAME loop rather than re-derived by a second pass over a re-built pair set. A
    // duplicated gate would hand back pairs that are not the ones that disagreed.
    let mut disagree = Vec::new();
    if anchors_a.is_empty() || anchors_b.is_empty() {
        return (out, disagree);
    }
    let dt_max = anchors_a
        .iter()
        .flat_map(|a| anchors_b.iter().map(move |b| (obs[b.first].epoch - obs[a.first].epoch).abs()))
        .fold(0.0_f64, f64::max);
    let cap = gate_radius_astrometric(node.r, dt_max, mu, sigma_cap, k);

    let pack = |set: &[Anchor]| {
        let o: Vec<Vector3<f64>> = set.iter().map(|a| obs[a.first].observer).collect();
        let u: Vec<Vector3<f64>> = set.iter().map(|a| obs[a.first].rho_hat).collect();
        let ids: Vec<u32> = (0..set.len() as u32).collect();
        BaryIndex::build(&o, &u, &ids, node.r)
    };
    let ia = pack(anchors_a);
    let ib = pack(anchors_b);

    for (slot, p) in ia.points.iter().enumerate() {
        let q = Vector3::new(p[0], p[1], p[2]);
        for hit in ib.within_idx(&q, cap) {
            let a = anchors_a[ia.ids[slot] as usize];
            let b = anchors_b[ib.ids[hit] as usize];
            // 🔴 THE ORDERING GUARD, as in solve_anchor_pair.
            let (a, b) = if obs[a.first].epoch <= obs[b.first].epoch { (a, b) } else { (b, a) };
            let (a0, b0) = (&obs[a.first], &obs[b.first]);
            let t_ref = 0.5 * (a0.epoch + b0.epoch);
            let pair = Pair { t_ref, a: a0.point(), b: b0.point() };
            if let Some(c) = pair.scan_census(node, mu) {
                let g = pair.guided_check(node, mu);
                if !g.outcome_agrees {
                    disagree.push(pair.clone());
                }
                out.push(((b0.epoch - a0.epoch).abs(), c, g));
            }
        }
    }
    (out, disagree)
}

/// Why anchor pairs did not become candidates.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RejectTally {
    pub unbound_node: usize,
    pub no_bound_root: usize,
    pub no_state: usize,
    pub no_prediction: usize,
    pub chi2: usize,
}

impl RejectTally {
    fn count(&mut self, r: Reject) {
        match r {
            Reject::UnboundNode => self.unbound_node += 1,
            Reject::NoBoundRoot => self.no_bound_root += 1,
            Reject::NoState => self.no_state += 1,
            Reject::NoPrediction => self.no_prediction += 1,
            Reject::Chi2 { .. } => self.chi2 += 1,
        }
    }

    pub fn total(&self) -> usize {
        self.unbound_node + self.no_bound_root + self.no_state + self.no_prediction + self.chi2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spacerocks::coordinates::Origin;

    fn mu() -> f64 {
        Origin::SSB.mu()
    }

    /// A synthetic object built from the canonical formulation, observed from a circular observer.
    ///
    /// ⭐ Generated by the SAME machinery the solve inverts. That is deliberate and it is what the
    /// Python self-test does: it does not validate the physics (the fixture does that), it
    /// validates that pairing -> solve -> chi2 INVERTS its own generation. That is exactly the
    /// property the endpoint-order bug broke, and it is how that bug was found.
    struct Truth {
        node: Node,
        h: f64,
        pair: Pair,
    }

    fn observer_at(t: f64) -> Vector3<f64> {
        // a 1 AU circular observer in the ecliptic; only its motion matters here
        let w = std::f64::consts::TAU / 365.25;
        Vector3::new((w * t).cos(), (w * t).sin(), 0.0)
    }

    /// Build a truth orbit and return an `Obs` at each requested epoch.
    fn synth(epochs: &[f64], r: f64, rdot: f64, sigma: f64) -> (Truth, Vec<Obs>) {
        let node = Node { r, rdot };
        // a bound h comfortably inside the ceiling
        let h_max = crate::spherical_pair::h_max(&node, mu()).expect("node must be bound");
        let h = 0.7 * h_max;
        // 🔴 The node is (r, rdot) AT t_ref, so the generator must use the SAME t_ref the solve
        // will: `solve_anchor_pair` takes the midpoint of the two anchors' first detections.
        // Generating at epochs[0] instead makes "r = 45" mean a different orbit -- rdot carries it
        // 1.5e-3 AU over the 15-day offset -- and the recovered h then misses by ~2.6e-4 relative
        // for reasons that look like a solver defect and are not.
        let t_ref = 0.5 * (epochs[0] + epochs[2]);

        // seed the plane with two arbitrary lines of sight, then read the true positions off the
        // canonical step -- this is the generator, and it needs a Pair only to fix the plane
        let seed = Pair {
            t_ref,
            a: PairPoint {
                epoch: epochs[0],
                rho_hat: Vector3::new(0.0, 1.0, 0.05).normalize(),
                observer: observer_at(epochs[0]),
            },
            b: PairPoint {
                epoch: epochs[1],
                rho_hat: Vector3::new(0.02, 1.0, 0.05).normalize(),
                observer: observer_at(epochs[1]),
            },
        };
        let mut out = Vec::new();
        for (i, &t) in epochs.iter().enumerate() {
            let p = position_at(&seed, &node, h, mu(), t).expect("truth position");
            let o = observer_at(t);
            let geo = p - o;
            out.push(Obs {
                id: i as u32,
                epoch: t,
                rho_hat: geo / geo.norm(),
                observer: o,
                sigma,
            });
        }
        (Truth { node, h, pair: seed }, out)
    }

    /// Two anchors: (0,1) on the first night, (2,3) a while later.
    fn epochs() -> Vec<f64> {
        let n0 = 2_460_000.5_f64;
        vec![n0, n0 + 0.025, n0 + 30.0, n0 + 30.025]
    }

    #[test]
    fn a_noiseless_synthetic_object_solves_to_a_vanishing_chi2() {
        let (truth, obs) = synth(&epochs(), 45.0, -1.0e-4, 0.12 / 206_264.806_247_096_36);
        let a = Anchor { first: 0, second: 1 };
        let b = Anchor { first: 2, second: 3 };
        let c = solve_anchor_pair(&obs, a, b, &truth.node, mu(), CHI2_CUT_DOF4)
            .expect("the true node must solve and pass");
        // with no astrometric noise the two unused detections must be predicted essentially exactly
        assert!(c.chi2 < 1e-6, "chi2 {} should be ~0 on a noiseless object", c.chi2);
        // and the recovered angular momentum is the generator's
        assert!(
            (c.h / truth.h - 1.0).abs() < 1e-9,
            "h {} vs truth {}",
            c.h,
            truth.h
        );
    }

    #[test]
    fn the_endpoint_order_is_enforced_and_it_matters() {
        let (truth, obs) = synth(&epochs(), 45.0, -1.0e-4, 0.12 / 206_264.806_247_096_36);
        let a = Anchor { first: 0, second: 1 };
        let b = Anchor { first: 2, second: 3 };

        // 🔴 the guard: handing the pair over reversed must give the IDENTICAL candidate
        let fwd = solve_anchor_pair(&obs, a, b, &truth.node, mu(), CHI2_CUT_DOF4).unwrap();
        let rev = solve_anchor_pair(&obs, b, a, &truth.node, mu(), CHI2_CUT_DOF4).unwrap();
        assert_eq!(fwd.epoch, rev.epoch, "the state epoch must be the earlier anchor's");
        assert!((fwd.h - rev.h).abs() <= f64::EPSILON * fwd.h.abs() * 8.0);
        assert!((fwd.chi2 - rev.chi2).abs() < 1e-9);

        // 🔴 And the guard is load-bearing -- but NOT the way the Python measurement described.
        // Build the Pair by hand in the wrong order: Stage 1's signed sweep runs past pi and the
        // solve refuses, so a reversed pair is LOST rather than answered wrongly. Without this
        // assertion the equality above would pass on a no-op guard.
        let t_ref = 0.5 * (obs[0].epoch + obs[2].epoch);
        let backwards = Pair { t_ref, a: obs[2].point(), b: obs[0].point() };
        match backwards.solve_h(&truth.node, mu()) {
            Solution::NoBoundRoot { .. } => {}
            Solution::Bound { h, .. } => panic!(
                "a reversed pair solved to h {h}: the pi guard in Stage 1 has stopped working, \
                 and this port has acquired the prototype's retrograde-orbit trap"
            ),
            Solution::UnboundNode => panic!("the node is bound; this cannot be UnboundNode"),
        }
        // the ORDERED pair, by contrast, predicts the unused detection to sub-milliarcsecond
        let fwd_pred = predict_hat(
            &Pair { t_ref, a: obs[0].point(), b: obs[2].point() },
            &truth.node,
            fwd.h,
            mu(),
            obs[1].epoch,
            &obs[1].observer,
        )
        .unwrap();
        let fwd_resid = (obs[1].rho_hat - obs[1].rho_hat.dot(&fwd_pred) * fwd_pred).norm()
            * 206_264.806_247_096_36;
        assert!(fwd_resid < 1e-3, "ordered residual {fwd_resid:.3e}\" should be ~0");
    }

    #[test]
    fn a_wrong_node_is_rejected_by_the_chi2_gate() {
        let (truth, obs) = synth(&epochs(), 45.0, -1.0e-4, 0.12 / 206_264.806_247_096_36);
        let a = Anchor { first: 0, second: 1 };
        let b = Anchor { first: 2, second: 3 };
        // the true node passes
        assert!(solve_anchor_pair(&obs, a, b, &truth.node, mu(), CHI2_CUT_DOF4).is_ok());
        // 🔴 but a wrong r must be REJECTED BY THE RESIDUAL, not by the solve: the solve is a
        // generator and will happily return a bound orbit at the wrong node
        let wrong = Node { r: 45.0 * 1.25, rdot: truth.node.rdot };
        match solve_anchor_pair(&obs, a, b, &wrong, mu(), CHI2_CUT_DOF4) {
            Err(Reject::Chi2 { chi2 }) => assert!(chi2 > CHI2_CUT_DOF4),
            Err(Reject::NoBoundRoot) => { /* also a real rejection, and free */ }
            other => panic!("a 25% wrong r should not pass: {other:?}"),
        }
    }

    #[test]
    fn anchors_are_ordered_in_time_however_the_visits_are_supplied() {
        let (_t, obs) = synth(&epochs(), 45.0, -1.0e-4, 0.12 / 206_264.806_247_096_36);
        let sigma_cap = obs[0].sigma * std::f64::consts::SQRT_2;
        // supply the LATER visit first
        let anchors = build_anchors(&obs, &[1], &[0], 45.0, 45.0, 3.0, mu(), sigma_cap);
        assert!(!anchors.is_empty(), "the object's own two detections must pair");
        for a in &anchors {
            assert!(
                obs[a.first].epoch <= obs[a.second].epoch,
                "anchor {a:?} is not ordered in time"
            );
        }
    }

    #[test]
    fn the_tally_accounts_for_every_pair_that_did_not_survive() {
        let (truth, obs) = synth(&epochs(), 45.0, -1.0e-4, 0.12 / 206_264.806_247_096_36);
        let a = vec![Anchor { first: 0, second: 1 }];
        let b = vec![Anchor { first: 2, second: 3 }];
        let sigma_cap = obs[0].sigma * std::f64::consts::SQRT_2;
        let (cands, tally) =
            anchor_pairs_and_solve(&obs, &a, &b, &truth.node, 3.0, mu(), sigma_cap, CHI2_CUT_DOF4);
        assert_eq!(cands.len() + tally.total(), 1, "every gated pair must be accounted for");
        assert_eq!(cands.len(), 1, "the true node should keep its own object");
    }

    /// 🔴 THE HOIST IS BIT-EXACT, and this is the test that keeps it so.
    ///
    /// `anchor_frame` + `position_at_in` must reproduce `position_at` to the LAST BIT at every
    /// epoch -- compared with `to_bits`, not with a tolerance. The entire argument for hoisting
    /// the `t`-invariant work out of Stage 4's loop is that nothing was recomputed differently
    /// and nothing was reassociated; a tolerance test would let a real numerical change hide
    /// inside an "acceptable" difference, which is exactly the class of silent drift this file
    /// keeps finding. The end-to-end version of this check is a candidate CSV that compares byte
    /// for byte (`_c2_hoist_regress*.yaml`); this is its unit-level guard.
    #[test]
    fn hoisted_frame_reproduces_position_at_bit_for_bit() {
        let node = Node { r: 45.0, rdot: 1.0e-3 };
        let h = 0.7 * crate::spherical_pair::h_max(&node, mu()).expect("node must be bound");
        let (t_a, t_b) = (0.0, 20.0);
        let pair = Pair {
            t_ref: 0.5 * (t_a + t_b),
            a: PairPoint {
                epoch: t_a,
                rho_hat: Vector3::new(0.0, 1.0, 0.05).normalize(),
                observer: observer_at(t_a),
            },
            b: PairPoint {
                epoch: t_b,
                rho_hat: Vector3::new(0.02, 1.0, 0.05).normalize(),
                observer: observer_at(t_b),
            },
        };
        let frame = anchor_frame(&pair, &node, h, mu()).expect("the anchors must build a frame");

        let mut compared = 0usize;
        for k in 0..400 {
            // well outside the anchor span in both directions: the hoist must not be exact only
            // where the two anchors are.
            let t = -60.0 + 0.4 * k as f64;
            let direct = position_at(&pair, &node, h, mu(), t);
            let hoisted = position_at_in(&frame, &pair, &node, h, mu(), t);
            match (direct, hoisted) {
                (Some(d), Some(x)) => {
                    for i in 0..3 {
                        assert_eq!(
                            d[i].to_bits(),
                            x[i].to_bits(),
                            "component {i} at t = {t} differs: {} vs {}",
                            d[i],
                            x[i]
                        );
                    }
                    compared += 1;
                }
                (None, None) => {}
                (d, x) => panic!("one path produced a position and the other did not at t = {t}: \
                                  {d:?} vs {x:?}"),
            }
        }
        // 🔴 A test that compared nothing would pass. Assert it actually ran.
        assert!(compared > 300, "only {compared} epochs produced a position; the fixture is wrong");
    }
}
