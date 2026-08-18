//! The spherical two-point method ("C2") -- Stage 1 of `PLAN_c2_pointdexter_rust.md`.
//!
//! C1 (`matt-0`, `tangent_v2`) GRIDS the velocity. C2 MEASURES it: assert `r` and `rdot` at
//! `t_ref`, take two topocentric unit vectors, and the system is exactly determined --
//! `6 = 2 + 4`. The one leftover scalar is the angular momentum `h`, and it falls out of a
//! **1-D** root find
//!
//! ```text
//!     F(h) = angle(r_vec_0, r_vec_1) - dnu(h) = 0
//! ```
//!
//! not Herget's 2-D `(rho_1, rho_3)`, because `(r, rdot)` is imposed analytically rather than
//! as residuals. `F` is monotonic in `h` with a single root, so a bracketed root find needs no
//! partials at all.
//!
//! This module is deliberately **I/O-free and `InitialCondition`-free**. The radial solve is
//! `spacerocks::solve_for_universal_anomaly` plus the f-and-g functions, so the 105-node
//! Chebyshev table -- and with it the silent `+-365.25 d` extrapolation of
//! `initial_condition.rs` -- is never built in the inner loop. Nothing here touches an
//! existing pointdexter file beyond one `pub mod` line in `lib.rs`.
//!
//! # The frame is BARYCENTRIC, and that is load-bearing
//!
//! `mu` and the origin move together: the constant is `Origin::SSB.mu()` and the observer
//! positions must be barycentric too. The range quadratic solves `|R + rho*rho_hat| = r`, so
//! `R` and the asserted `r` have to share an origin. A heliocentric mix-up is a ~0.13% error
//! in `mu` -- eight orders of magnitude above this module's test tolerances, and it reads as a
//! porting bug rather than a frame mismatch.
//!
//! # What the oracle is
//!
//! Every function here is gated against `fixtures/fh_fixture.json`, schema
//! `ssolink.c2.fh_fixture/3`, generated in ssolink by `_fh_fixture.py` from the Python
//! prototype. See `tests/spherical_pair_fixture.rs`.

use nalgebra::Vector3;
use spacerocks::transforms::solve_for_universal_anomaly;

/// Residual tolerance for `solve_for_universal_anomaly`.
///
/// Its criterion is `|f(chi)| <= tol` on a residual whose units are AU^1.5 and which ends in
/// `sqrt(mu)*dt`, i.e. it is effectively a TIME tolerance `tol/sqrt(mu)` -- independent of
/// distance, unlike the absolute-`s` tolerance that made `universal-kepler.c` fail worse the
/// farther out it ran. At `1e-12` that is ~6e-11 d, so a TNO's position lands to ~1e-13 AU:
/// two orders inside the fixture's `1e-12` budget, and still ~300x above the double-precision
/// noise floor of the residual itself (`sqrt(mu)*dt ~ 12` over a 2-yr baseline).
pub const ANOMALY_TOL: f64 = 1e-12;

/// Iteration cap for the same solver. It is bracketed with a bisection fallback, so this is a
/// backstop, not the convergence path.
pub const ANOMALY_MAX_ITER: usize = 100;

/// Fallback tolerance for the ONE case where `ANOMALY_TOL` is unreachable.
///
/// 🔴 Not a loosening of the solve. `spherical_pair::canonical_step` retries at this tolerance
/// only when the solver returns `Err`, and `polish_anomaly` then restores full precision on the
/// series residual. Measured (`RESULT_c2_nobound_conflation.md`): on near-parabolic nodes
/// (`a ~ 5e3-6e4 AU`) at short baselines, `solve_for_universal_anomaly` **stagnates** against
/// `ANOMALY_TOL` -- it fails identically at 100, 10^4 and **10^6** iterations, so the cap is not
/// the problem -- while succeeding at `1e-10` and `1e-8` and returning an `s` that agrees to 6-7
/// figures across all three. The root is well determined; only the acceptance test fails.
///
/// 🔴 It must NOT be applied globally. `ANOMALY_TOL`'s own note puts 1e-12 at ~1e-13 AU of
/// position, two orders inside the fixture's 1e-12 budget; a global 1e-8 would land far outside
/// it. The retry is reachable only where the code previously returned `None`, so every solve that
/// already succeeded -- and the frozen oracle -- is bit-unchanged.
pub const ANOMALY_TOL_RETRY: f64 = 1e-8;

/// How many times the `ANOMALY_TOL_RETRY` path has fired.
///
/// 🔴 The failure this repairs was silent for the method's whole lifetime because a solver
/// failure was indistinguishable from a physical rejection. A repair that is ALSO silent leaves
/// the next regression just as invisible: if the dependency's stagnating region ever widens, this
/// counter is the only thing that says so. A run should report it.
pub static ANOMALY_RETRY_COUNT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// `(r, rdot, h, dt)` of the FIRST retry, so a run can hand back an exact reproducer rather than a
/// count with no way to re-enter the case. The stagnating set is ~1e-12 wide in `h/h_max`, so a
/// point recorded to anything less than full precision cannot be re-found.
pub static ANOMALY_RETRY_FIRST: std::sync::OnceLock<[f64; 4]> = std::sync::OnceLock::new();

/// Half-cell tolerances of the root find on `h`, relative and absolute.
const H_XTOL: f64 = 1e-14;
const H_RTOL: f64 = 1e-15;
const H_MAX_ITER: usize = 200;

/// Number of geometric samples across the physical bracket `(0, h_max]`.
///
/// Matches the Python prototype's scan so both sides bracket the same sign change. `F` is
/// monotonic, so this is about *finding* the bracket, not about resolving the root.
pub const H_SCAN_N: usize = 80;

/// The ends of the physical bracket, as fractions of `h_max`. `1 - 1e-9` rather than `1.0`
/// because `F` is evaluated *at* the ceiling, where the orbit is marginally bound.
///
/// 🔴 `pub` so a diagnostic that re-walks the range walks the SAME range. A transcribed `1e-6`
/// that drifted from this one would report a root the solver "missed" that is simply outside
/// the solver's domain.
pub const H_SCAN_LO_FRAC: f64 = 1e-6;
pub const H_SCAN_HI_FRAC: f64 = 1.0 - 1e-9;

/// A grid node: the two quantities C2 asserts at `t_ref`.
///
/// `rdot` is the barycentric radial velocity, NOT BK's rate ratio `zdot/z`. Both appear in this
/// codebase (`tangent_v2.rs:430` uses the ratio) and they are not the same quantity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Node {
    /// Barycentric distance at `t_ref`, AU.
    pub r: f64,
    /// Barycentric radial velocity at `t_ref`, AU/day.
    pub rdot: f64,
}

/// One of the two detections a pair is built from.
#[derive(Debug, Clone, PartialEq)]
pub struct PairPoint {
    /// JD, TDB.
    pub epoch: f64,
    /// Topocentric line of sight, unit, ecliptic J2000.
    pub rho_hat: Vector3<f64>,
    /// Observer position, **barycentric** ecliptic J2000, AU.
    pub observer: Vector3<f64>,
}

/// A pair plus the epoch its node is asserted at.
///
/// `t_ref` lives here rather than being passed alongside because it is not separable from the
/// two points: `Node` means nothing without the epoch it is asserted at, and a pair whose
/// `t_ref` drifted from the one the grid was built for is a silent physical error.
#[derive(Debug, Clone, PartialEq)]
pub struct Pair {
    pub t_ref: f64,
    pub a: PairPoint,
    pub b: PairPoint,
}

/// What the `h` solve concluded, and **why**.
///
/// A bare `Option<f64>` was the prototype's first shape and it is not good enough: `None`
/// conflates "the asserted node is itself unbound", "no bound orbit threads both points" and
/// "the solver failed". Section 5 of the design note established that a wrong node ALWAYS has a
/// root -- the solve is a generator, not a filter -- so a bare `None` is never evidence about
/// the node. The rejection that *is* real lives in the second variant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Solution {
    /// `2mu/r - rdot^2 <= 0`: no bound orbit passes through the asserted node at all,
    /// independent of the observations.
    UnboundNode,
    /// The node is bound, but `F` holds its sign across the whole of `(0, h_max]`: the only
    /// orbit through both points at this asserted `r` is hyperbolic. A genuine rejection, and
    /// free -- the bracket scan has already been run. Both ends are carried so a caller (or a
    /// test) can check the conclusion was reached for the right reason.
    NoBoundRoot { h_max: f64, f_lo: f64, f_hi: f64 },
    /// A bound orbit threads both points.
    Bound { h: f64, h_max: f64 },
}

impl Solution {
    /// The angular momentum, if one was found.
    pub fn h(&self) -> Option<f64> {
        match self {
            Solution::Bound { h, .. } => Some(*h),
            _ => None,
        }
    }
}

/// The in-plane state of the canonical orbit: 2-D because it is planar by construction.
///
/// The plane and the in-plane phase are arbitrary -- only the radial solution and the swept
/// angle are used downstream, and both are invariant to that choice. Keeping it 2-D is what
/// makes the swept angle signed for free (see [`swept_angle`]).
#[derive(Debug, Clone, Copy)]
pub(crate) struct CanonicalState {
    pub(crate) pos: [f64; 2],
    pub(crate) vel: [f64; 2],
}

pub(crate) fn hypot2(v: [f64; 2]) -> f64 {
    (v[0] * v[0] + v[1] * v[1]).sqrt()
}

/// Below this the Stumpff series is used instead of the closed forms. At `|z| = 1` the series
/// still converges in ~9 terms, so this is well inside its comfortable range and well outside
/// the range where the closed forms have lost digits.
const STUMPFF_SERIES_Z: f64 = 1.0;

/// Stumpff `C(z) = c2(z)` and `S(z) = c3(z)`, evaluated by series near zero.
///
/// 🔴 `spacerocks::transforms::{stumpff_c, stumpff_s}` compute `(1 - cos sqrt z)/z` and
/// `(sqrt z - sin sqrt z)/z^1.5` unguarded, so their relative error grows like `1/z` as `z`
/// approaches zero -- roughly `4e-16/z` for `C`. That is invisible for a typical TNO node
/// (`z ~ 3e-5`, error `1e-11`, worth 1e-14 AU of position) but it is NOT invisible where C2
/// actually looks: the bracket scan evaluates `F` at `h_max*(1 - 1e-9)`, where the orbit is
/// parabolic to within double precision, `z ~ 8e-14`, and the closed forms carry a 1e-3
/// relative error. Measured against a DOP853 integration, calling them there put `r(t)` 7.8e-7
/// AU wrong -- six orders outside this port's own fixture tolerance.
///
/// The series `C = sum (-z)^k/(2k+2)!`, `S = sum (-z)^k/(2k+3)!` has no such cancellation: at
/// small `z` the first term simply dominates.
fn stumpff_c2_c3(z: f64) -> (f64, f64) {
    if z.abs() < STUMPFF_SERIES_Z {
        let (mut c2, mut c3) = (0.0f64, 0.0f64);
        let (mut term_c2, mut term_c3) = (0.5f64, 1.0f64 / 6.0); // k = 0
        let mut k = 0.0f64;
        loop {
            c2 += term_c2;
            c3 += term_c3;
            if term_c2.abs() < f64::EPSILON * 1e-2 && term_c3.abs() < f64::EPSILON * 1e-2 {
                return (c2, c3);
            }
            // term_{k+1} = -z * term_k / ((2k+3)(2k+4)) for C, /((2k+4)(2k+5)) for S.
            term_c2 *= -z / ((2.0 * k + 3.0) * (2.0 * k + 4.0));
            term_c3 *= -z / ((2.0 * k + 4.0) * (2.0 * k + 5.0));
            k += 1.0;
        }
    } else if z > 0.0 {
        let rz = z.sqrt();
        ((1.0 - rz.cos()) / z, (rz - rz.sin()) / (rz * rz * rz))
    } else {
        // 🔴 `(cosh(sqrt(-z)) - 1) / (-z)`, note the sign of the denominator.
        // `spacerocks::stumpff_c` divides by `z` on this branch, which returns a NEGATIVE c2 for
        // every hyperbolic argument. C2's own scan stays inside `(0, h_max]`, so `z >= 0` and
        // the branch is never taken here -- but anything that reuses those functions for
        // unbound orbits should know.
        let rz = (-z).sqrt();
        ((rz.cosh() - 1.0) / (-z), (rz.sinh() - rz) / (rz * rz * rz))
    }
}

/// Reciprocal semi-major axis, `alpha = -2E/mu`. Positive for a bound orbit.
///
/// Written as a factored difference of squares rather than as `2/r - v^2/mu`:
///
/// ```text
///     alpha = (h_max^2 - h^2) / (mu r^2) = (h_max - h)(h_max + h) / (mu r^2)
/// ```
///
/// The direct form cancels at the top of the bracket scan (`h = h_max*(1 - 1e-9)`, where the two
/// terms agree to ~1e-9 relative): measured there, it returns `alpha` 3.3e-6 relative away from
/// the factored value. This form is exact up to the ceiling and gives `alpha = 0` at `h = h_max`
/// identically, so the bound/unbound boundary is sharp rather than fuzzy at the one `h` where
/// that matters.
///
/// 🔴 Honest caveat: this is DEFENSIVE, not a fix for anything observed. The 7.8e-7 AU error
/// this port originally had at the ceiling came from the Stumpff closed forms
/// ([`stumpff_c2_c3`]), not from here -- `alpha` enters only through `z = alpha*s^2` (where
/// `s^2` is ~1e-3 and `C(z) -> 1/2` regardless) and through `1 - alpha*r0` (~4e-9), so a 3e-6
/// relative error in it moves nothing measurable. Reverting to the direct form survives the
/// mutation suite, which is the correct reading: no test depends on this, and none should be
/// contrived to.
///
/// The subtraction that remains, `2mu/r - rdot^2`, only cancels when the node itself is at
/// escape -- and that case is classified as [`Solution::UnboundNode`] before any of this runs.
fn alpha_of(node: &Node, h: f64, mu: f64) -> f64 {
    let h_max_sq = 2.0 * mu / node.r - node.rdot * node.rdot;
    if h_max_sq > 0.0 {
        let hm = node.r * h_max_sq.sqrt();
        (hm - h) * (hm + h) / (mu * node.r * node.r)
    } else {
        let v_t = h / node.r;
        (h_max_sq - v_t * v_t) / mu
    }
}

/// The bound-orbit ceiling on `h` for an asserted node: `h_max = r*sqrt(2mu/r - rdot^2)`.
///
/// `None` when the node is itself unbound. This is the whole physical bracket -- deliberately
/// NOT a window around a hinted `h`. Production has no hint; it runs the solve precisely
/// because the orbit is unknown, and a hint-centred window returns `None` for every node whose
/// root falls outside it. Measured on the fixture's object, a node at `dr/r = -60%` has its
/// root at `0.159*h_true` and one at `-80%` at `0.039*h_true`: both inside the physics, both
/// outside any plausible window.
pub fn h_max(node: &Node, mu: f64) -> Option<f64> {
    if !(node.r > 0.0) || !node.rdot.is_finite() || !(mu > 0.0) {
        return None;
    }
    let h_max_sq = 2.0 * mu / node.r - node.rdot * node.rdot;
    if !(h_max_sq > 0.0) {
        return None;
    }
    Some(node.r * h_max_sq.sqrt())
}

/// Newton refinement of the universal anomaly on the series-evaluated Kepler residual.
///
/// The residual is the standard one -- `(r0*vr0/sqrt(mu))*s^2*C + (1 - alpha*r0)*s^3*S + r0*s
/// - sqrt(mu)*dt` -- and its derivative is `r(s)`, which is positive, so the root is simple and
/// Newton converges quadratically from the neighbourhood the library solver already lands in.
/// Two steps are ample; the loop exits as soon as a step stops moving `s`.
fn polish_anomaly(mut s: f64, node: &Node, alpha: f64, mu: f64, dt: f64) -> f64 {
    let sqrt_mu = mu.sqrt();
    for _ in 0..2 {
        let z = alpha * s * s;
        let (c, sc) = stumpff_c2_c3(z);
        let f = (node.r * node.rdot / sqrt_mu) * s * s * c
            + (1.0 - alpha * node.r) * s * s * s * sc
            + node.r * s
            - sqrt_mu * dt;
        let df = (node.r * node.rdot / sqrt_mu) * s * (1.0 - z * sc)
            + (1.0 - alpha * node.r) * s * s * c
            + node.r;
        if !(df.abs() > 0.0) || !f.is_finite() {
            break;
        }
        let ds = f / df;
        s -= ds;
        if !s.is_finite() || ds.abs() <= f64::EPSILON * s.abs() {
            break;
        }
    }
    s
}

/// Advance the canonical in-plane state by `dt` with universal variables.
///
/// This is the whole radial solve. `None` on any non-finite input or a solver that fails to
/// bracket -- 🔴 which a caller must CLASSIFY rather than fold into "no solution": the solver's
/// `Err` on bracketing is a different fact from `F` holding its sign.
pub(crate) fn canonical_step(node: &Node, h: f64, dt: f64, mu: f64) -> Option<CanonicalState> {
    if !(node.r > 0.0) || !(h > 0.0) || !node.rdot.is_finite() || !dt.is_finite() || !(mu > 0.0) {
        return None;
    }
    let v_t = h / node.r;
    let alpha = alpha_of(node, h, mu);
    // 🔴 Retry on `Err`, do not discard the root. The dependency STAGNATES against `ANOMALY_TOL`
    // near the parabolic ceiling rather than converging slowly (10^6 iterations fail identically),
    // and the old `.ok()?` turned that into `None` -- which `solve_h` then reported as the
    // PHYSICAL rejection `NoBoundRoot`. See `RESULT_c2_nobound_conflation.md`.
    //
    // Accuracy is not traded away: `polish_anomaly` below is Newton on the series residual, whose
    // derivative is `r(s) > 0`, so the residual is strictly increasing and has EXACTLY ONE root.
    // A looser seed therefore cannot land on a different root -- it can only start further from
    // the same one, and two quadratic steps from ~1e-8 reach machine precision.
    let solve = |tol: f64| {
        solve_for_universal_anomaly(node.r, node.rdot, alpha, mu, dt, tol, ANOMALY_MAX_ITER).ok()
    };
    let s = match solve(ANOMALY_TOL) {
        Some(s) => s,
        None => {
            let s = solve(ANOMALY_TOL_RETRY)?;
            ANOMALY_RETRY_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let _ = ANOMALY_RETRY_FIRST.set([node.r, node.rdot, h, dt]);
            s
        }
    };
    if !s.is_finite() {
        return None;
    }

    // 🔴 Polish `s` on a well-conditioned residual before using it. `solve_for_universal_anomaly`
    // converges on its OWN residual, which is built from the unguarded Stumpff closed forms, so
    // near the parabolic ceiling it converges tightly to a slightly wrong root. Two Newton steps
    // on the same residual evaluated with the series (below) cost nothing and remove it;
    // `df/ds = r > 0` with a single root, so this cannot wander.
    let s = polish_anomaly(s, node, alpha, mu, dt);

    let z = alpha * s * s;
    let (c, sc) = stumpff_c2_c3(z);
    let sqrt_mu = mu.sqrt();

    // f and g in the convention initial_condition.rs already uses, so the two radial solves in
    // this repo cannot silently disagree.
    let f = 1.0 - s * s / node.r * c;
    let g = dt - s * s * s / sqrt_mu * sc;

    // r0_vec = (r, 0), v0_vec = (rdot, h/r): the canonical frame's only content.
    let pos = [f * node.r + g * node.rdot, g * v_t];
    let r = hypot2(pos);
    if !(r > 0.0) || !r.is_finite() {
        return None;
    }

    let f_dot = sqrt_mu / (r * node.r) * (alpha * s * s * s * sc - s);
    let g_dot = 1.0 - s * s / r * c;
    let vel = [f_dot * node.r + g_dot * node.rdot, g_dot * v_t];
    if !vel[0].is_finite() || !vel[1].is_finite() {
        return None;
    }

    Some(CanonicalState { pos, vel })
}

/// `(r(t), rdot(t))` from the asserted node and a trial `h`.
///
/// No `InitialCondition`, no Chebyshev table, so the interpolation-bounds hazard of
/// `initial_condition.rs` does not arise for the solve. (That precompute IS the right trade for
/// the extension step -- one fixed orbit, many epochs -- but not here, where every node in the
/// grid is a different orbit.)
pub fn radial_at(node: &Node, h: f64, dt: f64, mu: f64) -> Option<(f64, f64)> {
    let st = canonical_step(node, h, dt, mu)?;
    let r = hypot2(st.pos);
    let rdot = (st.pos[0] * st.vel[0] + st.pos[1] * st.vel[1]) / r;
    Some((r, rdot))
}

/// Barycentric position at distance `r_assumed` along the line of sight from `observer`.
///
/// The algebra of `sync.rs:24-27`, transcribed rather than called: `optimize_rho_v2` takes an
/// `&InitialCondition` and reaches for `r_at_epoch()`, which drags in the Chebyshev precompute
/// this module exists to avoid.
///
/// `None` when the asserted distance is unreachable along this line of sight (negative
/// discriminant) or when the root is behind the observer. Both are real geometric rejections,
/// not failures.
pub fn range_quadratic(observer: &Vector3<f64>, rho_hat: &Vector3<f64>, r_assumed: f64) -> Option<Vector3<f64>> {
    if !(r_assumed > 0.0) {
        return None;
    }
    let b = observer.dot(rho_hat);
    let disc = b * b + r_assumed * r_assumed - observer.dot(observer);
    if disc < 0.0 {
        return None;
    }
    let rho = -b + disc.sqrt();
    if rho <= 0.0 {
        return None;
    }
    Some(observer + rho * rho_hat)
}

fn wrap_0_2pi(x: f64) -> f64 {
    let t = x % std::f64::consts::TAU;
    if t < 0.0 {
        t + std::f64::consts::TAU
    } else {
        t
    }
}

/// The swept true anomaly between two canonical states, **accumulated signed**.
///
/// 🔴 This is the one place the port deliberately improves on the Python prototype. The
/// prototype takes `arctan2(|cross|, dot)` of the two propagated position vectors, which is
/// unsigned and therefore FOLDS at 180 deg -- a Trojan reaches that in ~6 yr (design note
/// sec 9). Here the canonical frame is 2-D, so each position has an honest in-plane angle and
/// the sweep is `wrap(phi1 - phi0)` plus whole revolutions.
///
/// The revolution count is exact, not estimated: a full period sweeps exactly `2*pi`, so the
/// leftover always sweeps less than `2*pi`, and `floor(dt/T)` is the count. An unbound orbit
/// sweeps less than `2*pi` in total, so there `T` does not exist and the count is zero.
pub(crate) fn swept_angle(node: &Node, h: f64, mu: f64, s0: &CanonicalState, s1: &CanonicalState, dt: f64) -> f64 {
    let phi0 = s0.pos[1].atan2(s0.pos[0]);
    let phi1 = s1.pos[1].atan2(s1.pos[0]);
    // h > 0 by construction, so motion is counter-clockwise and the sweep is positive.
    let partial = wrap_0_2pi(phi1 - phi0);

    let alpha = alpha_of(node, h, mu);
    if alpha > 0.0 && dt > 0.0 {
        let period = std::f64::consts::TAU / (mu.sqrt() * alpha * alpha.sqrt());
        if period.is_finite() && period > 0.0 {
            return std::f64::consts::TAU * (dt / period).floor() + partial;
        }
    }
    partial
}

/// First half-width of the guided bracket, as a multiplicative factor either side of the guess.
///
/// Chosen from measurement, not taste (`RESULT_c2_signcensus.md`): `h_guess/h_root` sits inside
/// +-0.3% at `rdot = 0` and inside +-5% at p99 near the bound limit, so one step of 1.2 brackets
/// the overwhelming majority on the first try.
const H_GUESS_KAPPA: f64 = 1.2;

/// How many geometric expansions before the guess is declared uninformative and the full scan
/// takes over. `1.2^12 ~ 8.9`, which covers the measured 5.4x tail with room to spare.
const H_GUESS_MAX_EXPAND: usize = 12;

/// The scan grid `solve_h` walks, geometric in `h/h_max`.
///
/// 🔴 Extracted so that a diagnostic which re-walks the scan CANNOT silently drift from the scan
/// the solver actually uses. A census taken on a grid that is not this one measures nothing.
pub fn h_scan_grid(h_max: f64) -> impl Iterator<Item = f64> {
    let ln_lo = H_SCAN_LO_FRAC.ln();
    let ln_hi = H_SCAN_HI_FRAC.ln();
    (0..H_SCAN_N)
        .map(move |i| h_max * (ln_lo + (ln_hi - ln_lo) * (i as f64) / ((H_SCAN_N - 1) as f64)).exp())
}

/// What one anchor pair's `F(h)` scan looks like. **DIAGNOSTIC ONLY** -- nothing in the search
/// consumes this, and `solve_h` is unchanged by its existence.
///
/// It exists to answer one question: is `F` single-rooted over the scan? `solve_h` returns the
/// first bracket it finds, so if there are several, a lookup grid of initial guesses would land
/// on a different root and quietly change which orbit each pair resolves to -- candidates would
/// still solve, residuals would still look right, and every internal check would agree.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanCensus {
    /// Sign changes of `F` between ADJACENT finite scan points: how many brackets existed.
    pub sign_changes: usize,
    /// Scan points where `F` was finite. `H_SCAN_N` minus this is domain gap.
    pub finite: usize,
    /// Maximal runs of finite points. >1 means the domain is fragmented; a sign change spanning
    /// a gap is not counted here, and `solve_h` does not count one either.
    pub segments: usize,
    /// `h/h_max` at the lower edge of the FIRST bracket -- the root `solve_h` returns, and the
    /// one any initial guess must reproduce. NaN when there is no bracket.
    pub first_bracket_frac: f64,
    /// `h/h_max` at the UPPER edge of the first bracket. One geometric grid step above the
    /// lower edge, ~1.19x. Recorded rather than derived so a consumer cannot assume the step.
    pub first_bracket_hi_frac: f64,
    /// `h/h_max` of the CONVERGED root, from the same Brent `solve_h` uses. NaN if no bracket.
    ///
    /// 🔴 Scoring a guess against the bracket's centre instead of this measures the GRID, not
    /// the guess: the bracket is ~1.19x wide, so a perfect guess still scatters over
    /// +-1.0914x -- and that is exactly the spread the first run of this census produced.
    /// A guess error below one half grid step is invisible without converging first.
    pub root_frac: f64,
    /// `h/h_max` of the closed-form guess from [`Pair::h_guess`]. NaN when the geometry cannot
    /// supply one. Carried here so the guess can be scored against the bracket that contains the
    /// true root, on the same pairs, without a second pass.
    pub h_guess_frac: f64,
    /// `h/h_max` at the lower edge of the LAST bracket. Equals `first_bracket_frac` when
    /// single-rooted; the spread between them is how far a guess could be wrong.
    pub last_bracket_frac: f64,
}

/// `solve_h` against `solve_h_guided` on one pair. **DIAGNOSTIC ONLY.**
///
/// 🔴 Both solvers are run and compared here rather than asserted equal, because the guided one
/// is NOT bit-identical by construction: Brent from a different bracket converges to the same
/// root, not to the same last bits. The acceptance criterion is a relative tolerance plus zero
/// tolerance on a changed OUTCOME.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GuidedCheck {
    /// Both returned `Bound`, or neither did.
    pub outcome_agrees: bool,
    /// `|h_guided - h_scan| / h_scan` when both are `Bound`; NaN otherwise.
    pub rel_diff: f64,
    /// Whether the scan-based solver found a root at all -- the population split that decides
    /// whether the guided path is a win at this node.
    pub scan_bound: bool,
    pub ns_scan: u64,
    pub ns_guided: u64,
}

impl Pair {
    /// Run both solvers on this pair and compare. See [`GuidedCheck`].
    pub fn guided_check(&self, node: &Node, mu: f64) -> GuidedCheck {
        let t0 = std::time::Instant::now();
        let a = self.solve_h(node, mu);
        let ns_scan = t0.elapsed().as_nanos() as u64;
        let t1 = std::time::Instant::now();
        let b = self.solve_h_guided(node, mu);
        let ns_guided = t1.elapsed().as_nanos() as u64;

        let ha = if let Solution::Bound { h, .. } = a { Some(h) } else { None };
        let hb = if let Solution::Bound { h, .. } = b { Some(h) } else { None };
        GuidedCheck {
            outcome_agrees: ha.is_some() == hb.is_some(),
            rel_diff: match (ha, hb) {
                (Some(x), Some(y)) if x != 0.0 => ((y - x) / x).abs(),
                _ => f64::NAN,
            },
            scan_bound: ha.is_some(),
            ns_scan,
            ns_guided,
        }
    }

    /// Walk `solve_h`'s scan and count the brackets instead of returning at the first.
    ///
    /// Mirrors `solve_h`'s gap handling exactly: a non-finite `F` resets the run, so a sign
    /// change across a domain gap is not a bracket here or there.
    pub fn scan_census(&self, node: &Node, mu: f64) -> Option<ScanCensus> {
        let h_max = h_max(node, mu)?;
        let mut out = ScanCensus {
            sign_changes: 0,
            finite: 0,
            segments: 0,
            first_bracket_frac: f64::NAN,
            first_bracket_hi_frac: f64::NAN,
            root_frac: f64::NAN,
            last_bracket_frac: f64::NAN,
            h_guess_frac: self.h_guess(node).map(|h| h / h_max).unwrap_or(f64::NAN),
        };
        let f = |h: f64| self.f_of_h(h, node, mu);
        let mut prev: Option<(f64, f64)> = None;
        for h in h_scan_grid(h_max) {
            let Some(v) = f(h) else {
                prev = None;
                continue;
            };
            out.finite += 1;
            if prev.is_none() {
                out.segments += 1;
            }
            if let Some((h_prev, v_prev)) = prev {
                if v_prev * v < 0.0 {
                    out.sign_changes += 1;
                    out.last_bracket_frac = h_prev / h_max;
                    if out.sign_changes == 1 {
                        out.first_bracket_frac = h_prev / h_max;
                        out.first_bracket_hi_frac = h / h_max;
                        // The same Brent `solve_h` runs, on the same bracket, so `root_frac` is
                        // the root the solver would have returned -- not an approximation to it.
                        out.root_frac =
                            brent(&f, h_prev, v_prev, h, v).map(|r| r / h_max).unwrap_or(f64::NAN);
                    }
                }
            }
            prev = Some((h, v));
        }
        Some(out)
    }

    /// `solve_h`, started from [`Pair::h_guess`] instead of from a blind 80-point scan.
    ///
    /// Bracket the guess, widen geometrically until `F` changes sign across the bracket, then
    /// hand that bracket to the same Brent. Typical cost is ~2 evaluations to bracket plus the
    /// Brent, against 160 evaluations for the scan.
    ///
    /// 🔴 **It expands until it OBSERVES a sign change; it never assumes one.** That is what
    /// makes it safe to start from an approximation: a guess that is wrong merely costs extra
    /// evaluations, and a guess that is useless costs a fallback to the full scan. It cannot
    /// silently converge to the wrong thing.
    ///
    /// 🔴 It is only equivalent to `solve_h` because `F` was measured single-rooted over
    /// 4.7M gated pairs (`RESULT_c2_signcensus.md`). With two roots, "the root nearest the
    /// guess" and "the smallest-h root" are different answers and this would quietly return the
    /// other one. Re-run that census before trusting this at baselines beyond ~41 d.
    ///
    /// ⚠️ For a pair with NO root this is *slower* than `solve_h`: it pays the expansion and
    /// then the full scan anyway. That population is ~1% at `rdot = 0` but over 90% near the
    /// bound limit, so the net gain is node-dependent -- measured, not assumed, in the census.
    pub fn solve_h_guided(&self, node: &Node, mu: f64) -> Solution {
        let Some(h_max) = h_max(node, mu) else {
            return Solution::UnboundNode;
        };
        let Some(h0) = self.h_guess(node) else {
            return self.solve_h(node, mu);
        };
        let f = |h: f64| self.f_of_h(h, node, mu);
        let h_lo = h_max * H_SCAN_LO_FRAC;
        let h_hi = h_max * H_SCAN_HI_FRAC;
        // The guess is geometry, not a bounded quantity: it can land outside the physical range.
        let h0 = h0.clamp(h_lo, h_hi);

        let mut kappa = H_GUESS_KAPPA;
        for _ in 0..H_GUESS_MAX_EXPAND {
            let lo = (h0 / kappa).max(h_lo);
            let hi = (h0 * kappa).min(h_hi);
            if let (Some(a), Some(b)) = (f(lo), f(hi)) {
                // A single root inside [lo, hi] is exactly what a sign change at the ENDS
                // detects; widening is what moves the ends past it when it is outside.
                if a * b < 0.0 {
                    if let Some(root) = brent(&f, lo, a, hi, b) {
                        return Solution::Bound { h: root, h_max };
                    }
                }
            }
            if lo <= h_lo && hi >= h_hi {
                break; // the whole physical range, and still no sign change at the ends
            }
            kappa *= H_GUESS_KAPPA;
        }
        // 🔴 The scan is the authority, not this. Never report NoBoundRoot on the strength of
        // the expansion alone -- that would convert a measured property into an assumption.
        self.solve_h(node, mu)
    }

    /// A closed-form initial guess for the root of `F`, from the pair's geometry alone.
    ///
    /// `r^2 dnu/dt = h`, so matching the observed barycentric opening angle across the baseline
    /// gives `h ~ theta * r^2 / dt`. Two `range_quadratic` calls and an `atan2`: **no
    /// propagation, no Kepler solve, and nothing the search does not already have.**
    ///
    /// 🔴 Approximate by construction. It evaluates both endpoints at the node's asserted `r`,
    /// where the solve lets `r_a` and `r_b` float with `h`, and it takes the sweep as linear in
    /// time. It is a starting point for a bracket, **never an answer** -- and it is only safe as
    /// a starting point because `F` was measured single-rooted (RESULT_c2_signcensus.md).
    pub fn h_guess(&self, node: &Node) -> Option<f64> {
        let dt = self.b.epoch - self.a.epoch;
        if !(dt > 0.0) {
            return None;
        }
        let p0 = range_quadratic(&self.a.observer, &self.a.rho_hat, node.r)?;
        let p1 = range_quadratic(&self.b.observer, &self.b.rho_hat, node.r)?;
        let theta = p0.cross(&p1).norm().atan2(p0.dot(&p1));
        let h = theta * node.r * node.r / dt;
        if h.is_finite() && h > 0.0 { Some(h) } else { None }
    }

    /// The residual: geometric opening angle minus the dynamical true-anomaly sweep.
    ///
    /// `None` is a classified non-answer, not zero: either the trial `h` cannot be propagated,
    /// or a line of sight cannot reach the implied distance, or -- see below -- the sweep has
    /// grown past what the geometric side can represent.
    ///
    /// 🔴 The geometric term is an angle between two position vectors, so it lives in
    /// `[0, pi]` and cannot represent a sweep beyond half a turn. Rather than fold silently
    /// (the prototype's behaviour), a sweep at or past `pi` returns `None`. C2's target regime
    /// is TNOs over two or more oppositions, where a 45 AU object sweeps a couple of degrees a
    /// year, so this never fires there -- but the Centaur/Trojan regime of sec 9 reaches it,
    /// and that is exactly where the prototype would have returned a plausible wrong number.
    pub fn f_of_h(&self, h: f64, node: &Node, mu: f64) -> Option<f64> {
        if !(h > 0.0) {
            return None;
        }
        let s0 = canonical_step(node, h, self.a.epoch - self.t_ref, mu)?;
        let s1 = canonical_step(node, h, self.b.epoch - self.t_ref, mu)?;

        let d_nu = swept_angle(node, h, mu, &s0, &s1, self.b.epoch - self.a.epoch);
        if !(d_nu < std::f64::consts::PI) {
            return None;
        }

        let p0 = range_quadratic(&self.a.observer, &self.a.rho_hat, hypot2(s0.pos))?;
        let p1 = range_quadratic(&self.b.observer, &self.b.rho_hat, hypot2(s1.pos))?;
        let d_theta = p0.cross(&p1).norm().atan2(p0.dot(&p1));

        let f = d_theta - d_nu;
        if f.is_finite() { Some(f) } else { None }
    }

    /// Bracket and solve `F(h) = 0` over the physical range `(0, h_max]`.
    ///
    /// No hint, because production has none. The scan is geometric in `h/h_max`, so the deep
    /// roots -- down to `0.04*h_max` on the fixture's own cases -- are sampled as densely in
    /// log space as the shallow ones.
    ///
    /// 🔴 It returns at the FIRST sign change, so "the" solution is the **smallest-h** root in
    /// the bracket. Any future initial guess or bracket lookup converges to whichever root sits
    /// nearest the guess, which is the same answer ONLY where `F` is single-rooted. Measure that
    /// with [`Pair::scan_census`] before adopting one -- it is not visible in any output here.
    pub fn solve_h(&self, node: &Node, mu: f64) -> Solution {
        let Some(h_max) = h_max(node, mu) else {
            return Solution::UnboundNode;
        };

        let f = |h: f64| self.f_of_h(h, node, mu);
        let h_lo = h_max * H_SCAN_LO_FRAC;
        let h_hi = h_max * H_SCAN_HI_FRAC;

        let mut prev: Option<(f64, f64)> = None;
        for h in h_scan_grid(h_max) {
            let Some(v) = f(h) else {
                // A non-finite stretch is a gap in the domain, not a sign change across it.
                prev = None;
                continue;
            };
            if let Some((h_prev, v_prev)) = prev {
                if v_prev * v < 0.0 {
                    if let Some(root) = brent(&f, h_prev, v_prev, h, v) {
                        return Solution::Bound { h: root, h_max };
                    }
                }
            }
            prev = Some((h, v));
        }

        Solution::NoBoundRoot {
            h_max,
            f_lo: f(h_lo).unwrap_or(f64::NAN),
            f_hi: f(h_hi).unwrap_or(f64::NAN),
        }
    }

    /// Rebuild the full 3-D barycentric state at the FIRST point's epoch from a converged `h`.
    ///
    /// The orbit plane comes from the two geometric positions, and the transverse direction
    /// from `n_hat x p0_hat` -- which is what fixes the sign of the velocity. Getting that
    /// backwards leaves `|r|`, `rdot` and `|r x v|` all correct, which is why the fixture pins
    /// the line-of-sight residual as well.
    pub fn state_from_solution(&self, h: f64, node: &Node, mu: f64) -> Option<[f64; 6]> {
        let (r0, rdot0) = radial_at(node, h, self.a.epoch - self.t_ref, mu)?;
        let (r1, _) = radial_at(node, h, self.b.epoch - self.t_ref, mu)?;

        let p0 = range_quadratic(&self.a.observer, &self.a.rho_hat, r0)?;
        let p1 = range_quadratic(&self.b.observer, &self.b.rho_hat, r1)?;

        let n = p0.cross(&p1);
        let nn = n.norm();
        if nn == 0.0 {
            return None;
        }
        let n_hat = n / nn;
        let p0_hat = p0 / p0.norm();
        let t_hat = n_hat.cross(&p0_hat); // in-plane, along the motion
        let v0 = rdot0 * p0_hat + (h / r0) * t_hat;

        Some([p0[0], p0[1], p0[2], v0[0], v0[1], v0[2]])
    }
}

/// MJH's pairing gate: the widest barycentric direction change reachable from an asserted `r`.
///
/// `theta_max = sqrt(2mu/r^3)*dt`, from taking `r(t)` linear and mapping to barycentric unit
/// vectors, so the maximum rate is the escape speed over `r`. `rdot` does not appear on
/// purpose: the bound is the whole speed budget, and any nonzero `rdot` spends part of it
/// radially, so this bounds every node at that `r`.
///
/// 🔴 Evaluate it at the grid cell's LOWER `r` edge, not its centre. The buffer is dominated
/// by the `r`-cell, and the lower edge is exact where a centre value needs a 2.8x blind pad.
///
/// 🔴 It bounds the BARYCENTRIC set. Against raw topocentric sky positions it is short by the
/// observer's own parallax -- measured 5.5x at 30 d and 3.0x at 120 d -- which cancels only at
/// an integer year. A production gate that applies this to sky positions at a non-integer-year
/// baseline silently drops true pairs.
pub fn gate_radius(r_lower_edge: f64, dt: f64, mu: f64) -> f64 {
    (2.0 * mu / (r_lower_edge * r_lower_edge * r_lower_edge)).sqrt() * dt
}

/// Brent's method on a bracket whose end values are already known.
///
/// Written out rather than pulled in: the only dependency that would supply it is a numerics
/// crate this repo does not carry, and the routine is short enough that a transcription is
/// cheaper than a new dependency in someone else's Cargo.toml.
///
/// 🔴 `pub` so a diagnostic can trace THIS routine by wrapping `f`, rather than re-implementing
/// it. A replica that drifted by one acceptance condition would take a different path through
/// the bracket and would not reproduce the failure it is being used to explain.
pub fn brent<F>(f: &F, mut a: f64, mut fa: f64, mut b: f64, mut fb: f64) -> Option<f64>
where
    F: Fn(f64) -> Option<f64>,
{
    if fa * fb > 0.0 {
        return None;
    }
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }

    let mut c = a;
    let mut fc = fa;
    let mut d = 0.0;
    let mut used_bisection = true;

    for _ in 0..H_MAX_ITER {
        if fb == 0.0 {
            return Some(b);
        }
        let tol = 2.0 * H_RTOL * b.abs() + 0.5 * H_XTOL;
        if (b - a).abs() < tol {
            return Some(b);
        }

        // Inverse quadratic interpolation when three distinct points are available, secant
        // otherwise.
        let mut s = if fa != fc && fb != fc {
            a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            b - fb * (b - a) / (fb - fa)
        };

        // The standard acceptance conditions; any failure falls back to bisection, which is
        // what makes the method as safe as bisection and usually much faster.
        let lo = (3.0 * a + b) / 4.0;
        let (bound_lo, bound_hi) = if lo < b { (lo, b) } else { (b, lo) };
        let reject = !(s > bound_lo && s < bound_hi)
            || (used_bisection && (s - b).abs() >= 0.5 * (b - c).abs())
            || (!used_bisection && (s - b).abs() >= 0.5 * (c - d).abs())
            || !s.is_finite();
        if reject {
            s = 0.5 * (a + b);
            used_bisection = true;
        } else {
            used_bisection = false;
        }

        let fs = f(s)?;
        d = c;
        c = b;
        fc = fb;
        if fa * fs < 0.0 {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut fa, &mut fb);
        }
    }
    Some(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use spacerocks::coordinates::Origin;

    /// The constant this module is meant to run on, resolved the way both existing binaries
    /// resolve it (`matt-0.rs` via config `ic_origin`, `tangent_v2.rs` hardcoded at :759).
    fn mu() -> f64 {
        Origin::SSB.mu()
    }

    fn node() -> Node {
        Node { r: 45.0, rdot: -1.0e-4 }
    }

    #[test]
    fn the_barycentric_constant_is_the_one_the_binaries_resolve() {
        // 🔴 Not a tautology: the two wrong answers sit 0.134% away and both are reachable by
        // habit -- Gauss's k^2 (which spacerocks carries, but only to normalise its MASSES
        // table) and Origin::SUN.
        assert_eq!(mu(), 2.9630927493968080e-04);
        assert!((mu() / Origin::SUN.mu() - 1.0).abs() > 1e-3);
        assert!((mu() / 0.01720209895_f64.powi(2) - 1.0).abs() > 1e-3);
    }

    #[test]
    fn the_radial_solve_returns_the_node_at_dt_zero() {
        // The identity that catches a sign or units error in the universal-variable step.
        // `r` comes back bit-exact; `rdot` is recovered as `(pos . vel)/|pos|`, so it makes one
        // multiply-divide round trip through `r` and lands within an ulp -- a property of the
        // recovery, not slack in the solve.
        let n = node();
        let h = 0.5 * h_max(&n, mu()).unwrap();
        let (r, rdot) = radial_at(&n, h, 0.0, mu()).unwrap();
        assert_eq!(r, n.r);
        assert!((rdot / n.rdot - 1.0).abs() <= f64::EPSILON, "rdot {rdot} != {}", n.rdot);
    }

    #[test]
    fn the_stumpff_series_agrees_with_the_closed_forms_where_those_are_trustworthy() {
        // Above |z| ~ 1e-2 the closed forms have not lost meaningful precision, so the series
        // must reproduce them -- this is what says the series itself is right.
        for z in [-4.0, -1.0, -0.5, 0.05, 0.5, 1.0, 4.0, 25.0] {
            let (c2, c3) = stumpff_c2_c3(z);
            let (want_c2, want_c3) = if z > 0.0 {
                let rz: f64 = z.sqrt();
                ((1.0 - rz.cos()) / z, (rz - rz.sin()) / (rz * rz * rz))
            } else {
                let rz: f64 = (-z).sqrt();
                ((rz.cosh() - 1.0) / (-z), (rz.sinh() - rz) / (rz * rz * rz))
            };
            assert!((c2 / want_c2 - 1.0).abs() < 1e-14, "z={z}: c2 {c2} != {want_c2}");
            assert!((c3 / want_c3 - 1.0).abs() < 1e-14, "z={z}: c3 {c3} != {want_c3}");
        }
    }

    #[test]
    fn the_stumpff_series_is_exact_where_the_closed_forms_collapse() {
        // 🔴 The regression for the near-parabolic defect. At z = 1e-13 the truncated series is
        // the truth to 1e-30, so this is a hard bound, not a comparison.
        let z = 1e-13;
        let (c2, c3) = stumpff_c2_c3(z);
        assert!((c2 - (0.5 - z / 24.0)).abs() < 1e-18, "c2 = {c2}");
        assert!((c3 - (1.0 / 6.0 - z / 120.0)).abs() < 1e-18, "c3 = {c3}");

        // ...and the unguarded closed forms are visibly not usable there. Measured, then
        // asserted well inside the measurement, so this documents the hazard rather than
        // depending on a particular libm's rounding.
        let rz: f64 = z.sqrt();
        let unguarded_c2 = (1.0 - rz.cos()) / z;
        let unguarded_c3 = (rz - rz.sin()) / (rz * rz * rz);
        let err = (unguarded_c2 - 0.5).abs().max((unguarded_c3 - 1.0 / 6.0).abs() * 3.0);
        assert!(err > 1e-6, "the closed forms did not visibly fail: c2={unguarded_c2} c3={unguarded_c3}");
    }

    #[test]
    fn the_radial_solve_is_right_at_the_marginally_bound_ceiling() {
        // 🔴 The end-to-end regression for the same defect, at the point the bracket scan
        // actually visits: h = h_max*(1 - 1e-9), where alpha ~ 5e-11 and the orbit is parabolic
        // to within double precision.
        //
        // The reference values are NOT from this code. They are the fixture's dr_frac = +0.50
        // node, propagated by a DOP853 integration at rtol 1e-13 (scipy), which
        // universal_kepler.py reproduces bit-for-bit. Before the Stumpff series and the Newton
        // polish, this port landed 7.8e-7 AU away -- six orders outside the fixture's tolerance.
        let n = Node { r: 78.62213858937609, rdot: 0.00010575694221156902 };
        let hm = h_max(&n, mu()).unwrap();
        assert!((hm - 0.2156937672846777).abs() < 1e-15, "h_max = {hm}");
        let h = hm * (1.0 - 1e-9);
        for (dt, want) in [(-182.5, 78.60363410650568), (182.5, 78.6422348690123)] {
            let (r, _) = radial_at(&n, h, dt, mu()).unwrap();
            assert!((r - want).abs() < 1e-12, "dt={dt}: r = {r}, DOP853 says {want}");
        }
    }

    #[test]
    fn the_canonical_step_conserves_angular_momentum() {
        // The invariant that catches a wrong f-and-g pairing: fdot/gdot errors leave |r| and
        // rdot plausible but break r x v.
        let n = node();
        let h = 0.6 * h_max(&n, mu()).unwrap();
        for dt in [-1460.0, -90.0, 7.0, 365.25, 3650.0] {
            let st = canonical_step(&n, h, dt, mu()).unwrap();
            let h_out = st.pos[0] * st.vel[1] - st.pos[1] * st.vel[0];
            assert!((h_out / h - 1.0).abs() < 1e-12, "dt={dt}: h {h_out} != {h}");
        }
    }

    #[test]
    fn h_max_is_the_bound_orbit_ceiling() {
        let n = node();
        let hm = h_max(&n, mu()).unwrap();
        // At h_max the orbit is marginally bound: alpha crosses zero.
        assert!(alpha_of(&n, hm * (1.0 - 1e-12), mu()) > 0.0);
        assert!(alpha_of(&n, hm * (1.0 + 1e-9), mu()) < 0.0);
        // ...and it crosses EXACTLY there, which is what the factored form buys: the direct
        // `2/r - v^2/mu` leaves ~1e-19 of cancellation noise at the ceiling, so the boundary
        // between "bound" and "unbound" would be fuzzy at the one h where it has to be sharp.
        assert_eq!(alpha_of(&n, hm, mu()), 0.0);
        // A node whose radial speed alone exceeds escape has no bound orbit at all.
        assert!(h_max(&Node { r: 45.0, rdot: 1.0 }, mu()).is_none());
    }

    #[test]
    fn the_range_quadratic_lands_at_the_asserted_distance() {
        let observer = Vector3::new(0.4, -0.9, 1e-5);
        let rho_hat = Vector3::new(0.2, -0.97, 0.1).normalize();
        let p = range_quadratic(&observer, &rho_hat, 45.0).unwrap();
        assert!((p.norm() / 45.0 - 1.0).abs() < 1e-14);
        // ...and in front of the observer, on the line of sight.
        let geo = p - observer;
        assert!(geo.dot(&rho_hat) > 0.0);
        assert!((geo.normalize() - rho_hat).norm() < 1e-14);
    }

    #[test]
    fn an_unreachable_distance_is_rejected_not_clamped() {
        let observer = Vector3::new(1.0, 0.0, 0.0);
        // Looking straight along +x from 1 AU, nothing on this ray is at 0.5 AU from the origin.
        let rho_hat = Vector3::new(1.0, 0.0, 0.0);
        assert!(range_quadratic(&observer, &rho_hat, 0.5).is_none());
    }

    /// A Jupiter-like circular orbit: `h = sqrt(mu*r)` at `rdot = 0`, period ~11.9 yr. The
    /// eccentric alternative (`0.9*h_max`) does NOT reach half a turn in 6 yr -- it sweeps
    /// 2.05 rad -- so the demonstration needs the near-circular case and a longer baseline.
    fn trojan_like() -> (Node, f64, f64) {
        let n = Node { r: 5.2, rdot: 0.0 };
        let h = (mu() * n.r).sqrt();
        (n, h, 8.0 * 365.25)
    }

    #[test]
    fn the_swept_angle_does_not_fold_past_half_a_turn() {
        // 🔴 The hazard the prototype carried: unsigned arctan2 folds at pi, and the folded value
        // is a plausible-looking wrong number rather than a failure.
        let (n, h, dt) = trojan_like();
        let s0 = canonical_step(&n, h, 0.0, mu()).unwrap();
        let s1 = canonical_step(&n, h, dt, mu()).unwrap();
        let sweep = swept_angle(&n, h, mu(), &s0, &s1, dt);

        let folded = {
            // The prototype's form, in the canonical plane.
            let dot = s0.pos[0] * s1.pos[0] + s0.pos[1] * s1.pos[1];
            let cross = (s0.pos[0] * s1.pos[1] - s0.pos[1] * s1.pos[0]).abs();
            cross.atan2(dot)
        };
        assert!(sweep > std::f64::consts::PI, "sweep {sweep} did not exceed pi");
        assert!(folded <= std::f64::consts::PI);
        assert!(sweep - folded > 0.1, "the fold is not being demonstrated");
    }

    #[test]
    fn a_pair_spanning_more_than_one_revolution_accumulates_the_revolutions() {
        // 🔴 The 8-yr case above still lies inside ONE revolution, so it exercises the guard but
        // not the revolution count -- dropping the count entirely would still pass it. Past a
        // full period the wrapped angle comes back UNDER pi, so a port without the count reports
        // a small, plausible sweep for a pair that has been round the Sun and back: the worst
        // possible failure, because nothing downstream looks wrong.
        let (n, h, _) = trojan_like();
        let period = std::f64::consts::TAU * (n.r * n.r * n.r / mu()).sqrt(); // circular: ~11.9 yr
        let dt = 1.26 * period;
        let s0 = canonical_step(&n, h, 0.0, mu()).unwrap();
        let s1 = canonical_step(&n, h, dt, mu()).unwrap();

        let sweep = swept_angle(&n, h, mu(), &s0, &s1, dt);
        assert!((sweep / (std::f64::consts::TAU * 1.26) - 1.0).abs() < 1e-6,
                "sweep {sweep} is not 1.26 revolutions");
        assert!(sweep > std::f64::consts::TAU, "the revolution is not being accumulated");

        // Without the count the leftover is 0.26 of a turn -- under pi, so it would sail
        // straight past the fold guard.
        let partial = wrap_0_2pi(s1.pos[1].atan2(s1.pos[0]) - s0.pos[1].atan2(s0.pos[0]));
        assert!(partial < std::f64::consts::PI, "the trap is not being demonstrated");
    }

    #[test]
    fn f_of_h_refuses_a_folded_sweep_instead_of_returning_a_number() {
        let (n, h, dt) = trojan_like();
        let t_ref = 2460676.5;
        let pair = Pair {
            t_ref,
            a: PairPoint {
                epoch: t_ref,
                rho_hat: Vector3::new(1.0, 0.0, 0.0),
                observer: Vector3::new(0.0, 1.0, 0.0),
            },
            b: PairPoint {
                epoch: t_ref + dt,
                rho_hat: Vector3::new(-1.0, 0.0, 0.0),
                observer: Vector3::new(0.0, 1.0, 0.0),
            },
        };
        assert!(pair.f_of_h(h, &n, mu()).is_none());
        // ...and the same geometry over a baseline inside half a turn still answers, so the guard
        // is the sweep and not something about this synthetic pair.
        let short = Pair { b: PairPoint { epoch: t_ref + 365.25, ..pair.b.clone() }, ..pair.clone() };
        assert!(short.f_of_h(h, &n, mu()).is_some());
    }

    /// 🔴 REGRESSION, `RESULT_c2_nobound_conflation.md`. Captured by `ANOMALY_RETRY_FIRST` from a
    /// real census run; full precision is not decoration, the stagnating set is ~1e-12 wide in
    /// `h/h_max` and a point rounded to fewer digits does not re-enter it.
    ///
    /// Before the retry, `canonical_step` returned `None` here, `f_of_h` was non-finite at this
    /// exact `h`, `brent` aborted on it, and `solve_h` reported the PHYSICAL rejection
    /// `NoBoundRoot`. The pair was a real bound orbit with a real sign change on the scan grid.
    #[test]
    fn canonical_step_survives_the_dependency_stagnating_near_parabolic() {
        let node = Node { r: 3.40000000000000000e1, rdot: 4.00000000000000008e-3 };
        let h = 4.05837124060365595e-2;
        let dt = 1.03887515938840806e1;
        let mu = mu();

        // 1. The tight tolerance still fails. If a dependency upgrade ever fixes this, THIS
        //    assertion fires and tells us the retry has become dead code -- rather than the retry
        //    quietly covering for something that no longer happens.
        let alpha = alpha_of(&node, h, mu);
        assert!(
            solve_for_universal_anomaly(node.r, node.rdot, alpha, mu, dt, ANOMALY_TOL, ANOMALY_MAX_ITER)
                .is_err(),
            "the dependency now converges at ANOMALY_TOL here; the retry may be removable"
        );

        // 2. The step nonetheless produces a state.
        let st = canonical_step(&node, h, dt, mu).expect("retry path must produce a state");

        assert!(hypot2(st.pos).is_finite());

        // 3. 🔴 And it is ACCURATE -- the whole question, since a retry that merely returned A
        //    NUMBER would be worse than the `None` it replaces.
        //
        //    ⚠️ NOT tested by energy conservation, which is what this test first tried. The f-and-g
        //    functions place the state on the SAME CONIC for any `s` whatever, so specific energy
        //    is conserved identically (measured: relative error exactly 0.0) even when `s` does
        //    not solve Kepler's equation at all. It tests the orbit and never the TIME, i.e. it is
        //    a control the construction already guarantees.
        //
        //    The quantity that does move is the Kepler residual itself: `s` is defined by
        //    `f(s) = 0`, and the looser seed is accepted at ~`ANOMALY_TOL_RETRY`. This asserts
        //    that `polish_anomaly` closes that gap rather than the retry banking a sloppy root.
        let sqrt_mu = mu.sqrt();
        let residual = |s: f64| {
            let z = alpha * s * s;
            let (c, sc) = stumpff_c2_c3(z);
            (node.r * node.rdot / sqrt_mu) * s * s * c
                + (1.0 - alpha * node.r) * s * s * s * sc
                + node.r * s
                - sqrt_mu * dt
        };
        let s_seed =
            solve_for_universal_anomaly(node.r, node.rdot, alpha, mu, dt, ANOMALY_TOL_RETRY, ANOMALY_MAX_ITER)
                .expect("the relaxed tolerance is what makes the retry possible");
        let s_polished = polish_anomaly(s_seed, &node, alpha, mu, dt);
        let (before, after) = (residual(s_seed).abs(), residual(s_polished).abs());
        assert!(
            after < 1e-15,
            "polish left the Kepler residual at {after:.3e} (seed was {before:.3e})"
        );
        assert!(after < before, "polish did not improve the seed: {before:.3e} -> {after:.3e}");
    }

    #[test]
    fn brent_finds_a_root_it_is_handed() {
        let f = |x: f64| Some(x * x * x - 2.0 * x - 5.0); // classic, root at 2.0945514815...
        let root = brent(&f, 2.0, f(2.0).unwrap(), 3.0, f(3.0).unwrap()).unwrap();
        assert!((root - 2.0945514815423265).abs() < 1e-12, "root was {root}");
    }

    #[test]
    fn brent_refuses_a_bracket_with_no_sign_change() {
        let f = |x: f64| Some(x * x + 1.0);
        assert!(brent(&f, 1.0, f(1.0).unwrap(), 2.0, f(2.0).unwrap()).is_none());
    }

    #[test]
    fn the_gate_is_the_documented_closed_form_and_scales_as_expected() {
        let mu = mu();
        assert!((gate_radius(45.0, 365.25, mu) - (2.0 * mu / 45.0_f64.powi(3)).sqrt() * 365.25).abs() < 1e-18);
        // Linear in dt, and r^-1.5 in distance: a farther cell needs a tighter search radius.
        assert!((gate_radius(45.0, 240.0, mu) / gate_radius(45.0, 120.0, mu) - 2.0).abs() < 1e-12);
        let ratio = gate_radius(100.0, 30.0, mu) / gate_radius(400.0, 30.0, mu);
        assert!((ratio - 4.0_f64.powf(1.5)).abs() < 1e-12, "ratio {ratio}");
    }
}
