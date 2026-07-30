//! The C2 `(r, rdot)` grid -- Stage 2 of `PLAN_c2_pointdexter_rust.md`.
//!
//! C2 asserts `r` and `rdot` at `t_ref` and measures everything else, so the only thing it has
//! to grid is that two-dimensional node set. The cell is **not a fitted power law**: it is a
//! closed form evaluated from the arc's own observer geometry
//! (`RESULT_c2_analytic_spacing.md`), which predicts every half-cell the offset scans measured
//! to a median 1.01 (`r`) / 1.00 (`rdot`):
//!
//! ```text
//!     parallax   dtheta/dgamma = -E_perp(t)        =>  dgamma = eps / W
//!     zoom       dtheta/dgdot  = -theta_vec(t)*t   =>  drdot  = r * eps / Z
//! ```
//!
//! with `gamma = 1/r`, `gdot = rdot/r`, and `W` / `Z` the **anchored residuals** of those two
//! signatures over the extension window.
//!
//! This module is a transcription of ssolink's `kernels/spherical-pair/grid.py`, gated on the
//! `grid_cases` / `ladder_case` blocks of `fixtures/fh_fixture.json` (schema
//! `ssolink.c2.fh_fixture/3`). It is deliberately **ephemeris-free and orbit-free**: everything
//! consumes a [`Geometry`] -- epochs, observer positions, lines of sight -- and returns numbers.
//! That is what makes it testable with no ephemeris and no network, and it is why the port is a
//! transcription rather than a reimplementation: the binary already holds exactly these arrays.
//!
//! # Three things a reader coming from `tangent_v2.rs` / `grid.rs` will want to "improve"
//!
//! All three are load-bearing, and all three are wrong to change:
//!
//! * 🔴 **The projector is the affine function through the two ANCHORS**, not the least-squares
//!   `{1, t}` detrend of `tangent_v2.rs:515-520`. C2's solve pins the direction *exactly* at
//!   `t0` and `t1` by re-solving `h`; C1's fit does not. Sizing C2's cells with C1's projector
//!   sizes them for the wrong method. The single assertion that separates a correct port from
//!   one that copied C1 is that [`anchored_residual`] vanishes **exactly** at both anchors.
//! * 🔴 **Parallax and mean motion are different KINDS of tolerance.** `dgamma = eps/W` is
//!   *absolute* in gamma (=> uniform in `1/r`); the mean-motion branch reduces to
//!   `dgamma/gamma = eps/(3A)`, which is *relative* (=> uniform in `log r`, exactly what
//!   `grid.rs:251-278` implements as `r_max/r_min < r_tol`). Parallax binds by 8-22,000x over
//!   the TNO range and the two cross only inside ~10-20 AU, so **do not reuse the ratio
//!   criterion for C2**. [`binding_branch`] returns both so the choice can be checked rather
//!   than assumed.
//! * 🔴 **The ladder is re-sized at every shell it lands on.** `W` carries no `r`, but it drifts
//!   ~1.14x over 80-1600 AU and moves up to **2.5x with ecliptic latitude**. A fixed step would
//!   assume it is exactly constant; re-evaluating costs one lever evaluation per shell and keeps
//!   the grid from quietly under-sampling.
//!
//! # One grid is one sky region
//!
//! 🔴 [`build_ladder`] emits a ladder for the geometry its callback hands back. Building one
//! global grid at the ecliptic silently loses high-latitude objects; building one globally at
//! high latitude pays for that everywhere. And the bill is not the lever ratio: `W` and the
//! shell count both go ~2.5x from the ecliptic to `beta = -57 deg`, but **nodes go 4.9x**,
//! because shells are uniform in gamma so a tighter step adds them preferentially at small `r`,
//! where the `rdot` axis has not yet collapsed to one node (`RESULT_c2_grid_ladder.md`).
//!
//! # `mu` is barycentric
//!
//! As everywhere in C2: `Origin::SSB.mu()`, matching the frame the `Geometry` is expressed in.
//! `mu` enters only through [`rdot_span`] and the mean-motion branch, so a heliocentric slip
//! here is a quiet 0.07% in the `rdot` node count rather than a loud failure.

use nalgebra::Vector3;

/// What went wrong, kept as a type rather than an `Option` for one reason: a caller that
/// silently accepted an infinite cell would emit a one-node grid and report full coverage.
#[derive(Debug, Clone, PartialEq)]
pub enum GridError {
    /// The window carries no leverage on an axis, so that axis has no finite cell. The usual
    /// cause is a window whose epochs all coincide with the anchors, where the anchored residual
    /// vanishes identically by construction.
    DegenerateGeometry(String),
    /// The inputs are not a self-consistent geometry (mismatched lengths, empty window).
    Invalid(String),
    /// The ladder ran away. `eps_rad` is probably too small, or the geometry degenerate.
    LadderOverflow { max_nodes: usize, r_min: f64 },
}

impl std::fmt::Display for GridError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GridError::DegenerateGeometry(m) => write!(f, "degenerate geometry: {m}"),
            GridError::Invalid(m) => write!(f, "invalid geometry: {m}"),
            GridError::LadderOverflow { max_nodes, r_min } => write!(
                f,
                "ladder exceeded {max_nodes} nodes before reaching r_min={r_min}; \
                 eps_rad is probably too small or the geometry degenerate"
            ),
        }
    }
}

impl std::error::Error for GridError {}

/// Which statistic collapses the per-epoch residuals into a lever.
///
/// The median is what every recorded number was measured with. `P90` is kept because the choice
/// is a policy question -- how much of the window the cell must hold to `eps` -- not a fact, and
/// a port that hardcoded the median would make that policy invisible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stat {
    Median,
    P90,
}

/// The arc geometry both levers are computed from. Angles in radians, lengths in AU.
///
/// `t0`/`t1` are the two **anchor** epochs -- the pair the C2 solve consumes -- and `ts` are the
/// epochs of the **extension** window, the exposures the candidate orbit will be tested against.
///
/// 🔴 Sizing on the window rather than on one epoch is not optional. A wrong `r` is a parallax
/// error that nearly cancels at opposition, so a single opposition epoch is a near-perfect null
/// for it and reads ~100x too forgiving (`RESULT_c2_strict_metric.md`). The two axes notch at
/// *opposite* epochs -- the `rdot` error is an along-track drift that peaks exactly where the
/// parallax error vanishes -- so no blanket factor can rescue a single-epoch metric on both.
///
/// `u*` are unit lines of sight and `E*` the observer positions, in the same (barycentric) frame.
#[derive(Debug, Clone, PartialEq)]
pub struct Geometry {
    pub t0: f64,
    pub t1: f64,
    pub e0: Vector3<f64>,
    pub e1: Vector3<f64>,
    pub u0: Vector3<f64>,
    pub u1: Vector3<f64>,
    pub ts: Vec<f64>,
    pub es: Vec<Vector3<f64>>,
    pub us: Vec<Vector3<f64>>,
}

impl Geometry {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        t0: f64,
        t1: f64,
        e0: Vector3<f64>,
        e1: Vector3<f64>,
        u0: Vector3<f64>,
        u1: Vector3<f64>,
        ts: Vec<f64>,
        es: Vec<Vector3<f64>>,
        us: Vec<Vector3<f64>>,
    ) -> Result<Self, GridError> {
        if ts.is_empty() {
            return Err(GridError::DegenerateGeometry("extension window is empty".into()));
        }
        if es.len() != ts.len() || us.len() != ts.len() {
            return Err(GridError::Invalid(format!(
                "es/us must have {} entries; got {} / {}",
                ts.len(),
                es.len(),
                us.len()
            )));
        }
        if t1 == t0 {
            return Err(GridError::DegenerateGeometry("the two anchor epochs coincide".into()));
        }
        Ok(Geometry { t0, t1, e0, e1, u0, u1, ts, es, us })
    }

    pub fn baseline_days(&self) -> f64 {
        (self.t1 - self.t0).abs()
    }
}

// ------------------------------------------------------------------------------- the projector

/// Per-epoch norm of `vals(t)` after removing the affine function through the two anchors.
///
/// ⭐ This is the C2 projector, and the reason a signature that is affine in `t` costs nothing:
/// the solve reproduces such a signature exactly by re-solving `h`, so only the *unabsorbed*
/// part sets the cell. By construction the result is **identically zero at `t0` and `t1`** --
/// exactly zero, not small -- which is the assertion that catches a port that reached for the
/// least-squares detrend next door in `tangent_v2.rs`.
///
/// 🔴 One deliberate divergence from the Python prototype, and the reason it is here: the
/// prototype writes the interpolation as `v0 + s*(v1 - v0)`, which is exact at `t0` (where
/// `s = 0`) but **not** at `t1` -- `v0 + 1.0*(v1 - v0)` misses `v1` by an ulp, so the residual
/// there came out at 1.1e-16 rather than 0. That is why `test_c2_grid.py` asserts the anchor
/// property with `abs=1e-14` instead of exact equality. The endpoint form `(1-s)*v0 + s*v1` is
/// algebraically identical, costs the same, and is exact at *both* ends, so this port can gate
/// the invariant as an equality. The two forms differ by an ulp, four orders inside the
/// fixture's `lever_rel = 1e-12`, so nothing the oracle records moves.
pub fn anchored_residual(
    vals: &[Vector3<f64>],
    ts: &[f64],
    t0: f64,
    t1: f64,
    v0: Vector3<f64>,
    v1: Vector3<f64>,
) -> Vec<f64> {
    ts.iter()
        .zip(vals.iter())
        .map(|(t, v)| {
            let s = (t - t0) / (t1 - t0);
            (v - ((1.0 - s) * v0 + s * v1)).norm()
        })
        .collect()
}

/// Scalar form of [`anchored_residual`], for the 1-D mean-motion signature.
fn anchored_residual_scalar(vals: &[f64], ts: &[f64], t0: f64, t1: f64, v0: f64, v1: f64) -> Vec<f64> {
    ts.iter()
        .zip(vals.iter())
        .map(|(t, v)| {
            let s = (t - t0) / (t1 - t0);
            (v - ((1.0 - s) * v0 + s * v1)).abs()
        })
        .collect()
}

/// 🔴 Matches numpy's conventions exactly, because the fixture was generated through them: the
/// median of an even-length sample is the **mean of the two middle values** (not the lower one),
/// and `p90` is numpy's default *linear* interpolation between order statistics. The extension
/// windows in the fixture have 24 epochs, so the even case is the one that is actually exercised.
fn stat_of(values: &[f64], stat: Stat) -> f64 {
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).expect("NaN in residuals"));
    let n = v.len();
    match stat {
        Stat::Median => {
            if n % 2 == 1 {
                v[n / 2]
            } else {
                0.5 * (v[n / 2 - 1] + v[n / 2])
            }
        }
        Stat::P90 => {
            let pos = 0.90 * (n as f64 - 1.0);
            let lo = pos.floor() as usize;
            let hi = pos.ceil() as usize;
            v[lo] + (pos - lo as f64) * (v[hi] - v[lo])
        }
    }
}

/// Observer position perpendicular to the line of sight -- the spherical `(xe, ye)`.
fn perp(e: Vector3<f64>, u: Vector3<f64>) -> Vector3<f64> {
    e - e.dot(&u) * u
}

// ------------------------------------------------------------------------------ the two levers

/// `W` [AU]: the unabsorbed parallax signature. `dgamma = eps / W`.
///
/// ⭐ Note what is *not* in here: `r`. `W` is pure observer geometry, which is why uniform-in-
/// `1/r` is **derived rather than fitted** -- asserting `r` *is* asserting `gamma`, and a gamma
/// error shows up as `dtheta ~ dgamma * W` with `W` a lever independent of the orbit and of the
/// reach. It does depend strongly on *direction*: 0.60 AU at the ecliptic against 1.49 AU at
/// `beta = -57 deg`. A longer parallax lever means a TIGHTER cell, so high-latitude fields are
/// easier to measure and more expensive to cover.
pub fn parallax_lever(geom: &Geometry, stat: Stat) -> f64 {
    let ws: Vec<Vector3<f64>> = geom
        .es
        .iter()
        .zip(geom.us.iter())
        .map(|(e, u)| perp(*e, *u))
        .collect();
    let d = anchored_residual(
        &ws,
        &geom.ts,
        geom.t0,
        geom.t1,
        perp(geom.e0, geom.u0),
        perp(geom.e1, geom.u1),
    );
    stat_of(&d, stat)
}

/// `Z` [rad day]: the unabsorbed zoom signature. `dgdot = eps / Z`, `drdot = r * dgdot`.
///
/// A wrong `rdot` rescales the apparent size of the sky track about the anchor direction and
/// grows linearly in `t`, so its error peaks at the window *end* -- exactly where the parallax
/// error notches.
///
/// 🔴 `Z` *does* depend on `r` (through the sky track), because it mixes the parallactic reflex
/// with the object's own motion in a reach-dependent proportion. That mixture is what the fitted
/// `r^1.90` / `r^2.48` exponents were describing. **Measure it per shell; do not hardcode an
/// exponent** -- they are `Z`'s r-dependence at two particular reaches, not a law.
pub fn zoom_lever(geom: &Geometry, stat: Stat) -> f64 {
    let u0 = geom.u0;
    let t0 = geom.t0;
    // Small-angle offset from the anchor direction, scaled by elapsed time.
    let sig = |u: Vector3<f64>, t: f64| (u - u.dot(&u0) * u0) * (t - t0);

    let zs: Vec<Vector3<f64>> = geom
        .ts
        .iter()
        .zip(geom.us.iter())
        .map(|(t, u)| sig(*u, *t))
        .collect();
    let d = anchored_residual(&zs, &geom.ts, t0, geom.t1, sig(u0, t0), sig(geom.u1, geom.t1));
    stat_of(&d, stat)
}

/// `T2` [day^2]: the unabsorbed `t^2` signature, for the sky-curvature branch.
///
/// Kept so the two branches can be *compared* rather than assumed. The swept angle is
/// `A = 0.5*mu*gamma^3*t^2`, so this branch's tolerance is **relative** --
/// `dgamma/gamma = eps/(3A)` -- which is the `log r` grid, not the `1/r` one. See
/// [`binding_branch`].
pub fn meanmotion_lever(geom: &Geometry, stat: Stat) -> f64 {
    let t0 = geom.t0;
    let t1 = geom.t1;
    let vals: Vec<f64> = geom.ts.iter().map(|t| (t - t0) * (t - t0)).collect();
    let d = anchored_residual_scalar(&vals, &geom.ts, t0, t1, 0.0, (t1 - t0) * (t1 - t0));
    stat_of(&d, stat)
}

// ------------------------------------------------------------------------------------- cells

/// Half-cell in `gamma = 1/r` [AU^-1] for an astrometric budget `eps_rad`.
pub fn gamma_cell(w: f64, eps_rad: f64) -> Result<f64, GridError> {
    if !(w > 0.0) {
        return Err(GridError::DegenerateGeometry(
            "parallax lever W is zero: the window has no parallax leverage".into(),
        ));
    }
    Ok(eps_rad / w)
}

/// Half-cell in `rdot` [AU/day] at distance `r`.
pub fn rdot_cell(r: f64, z: f64, eps_rad: f64) -> Result<f64, GridError> {
    if !(z > 0.0) {
        return Err(GridError::DegenerateGeometry(
            "zoom lever Z is zero: the window has no along-track leverage".into(),
        ));
    }
    Ok(r * eps_rad / z)
}

/// Full width of the physically reachable `rdot` at `r`: bound orbits only.
///
/// `+-sqrt(2mu/r)` is escape, so this is the widest the axis could ever need to be. A real
/// population is narrower, which makes a node count built on it an **upper bound** -- the safe
/// direction for a grid-sizing claim.
pub fn rdot_span(r: f64, mu: f64) -> f64 {
    2.0 * (2.0 * mu / r).sqrt()
}

/// How many `rdot` nodes a shell needs. At least one: the shell always exists.
pub fn n_rdot_nodes(r: f64, drdot: f64, mu: f64) -> Result<usize, GridError> {
    let n = (rdot_span(r, mu) / (2.0 * drdot)).ceil();
    if !n.is_finite() {
        return Err(GridError::DegenerateGeometry(format!(
            "rdot cell {drdot} gives a non-finite node count at r={r}"
        )));
    }
    Ok((n as usize).max(1))
}

/// Which mechanism sets the `r` cell, and by how much.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Binds {
    Parallax,
    MeanMotion,
}

impl Binds {
    /// The spelling the fixture uses.
    pub fn as_str(&self) -> &'static str {
        match self {
            Binds::Parallax => "parallax",
            Binds::MeanMotion => "meanmotion",
        }
    }
}

/// Both `r`-cell branches as *relative* cells, so they are comparable at all, plus their ratio.
///
/// `ratio > 1` means the mean-motion tolerance is looser, i.e. **parallax binds** and the grid
/// variable is `1/r`. Over the TNO range that ratio runs 8 to 22,000; the two cross only around
/// 10-20 AU, inside the Centaur regime where the whole C2 strategy inverts anyway (short
/// baselines, not long ones -- design note §9).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BindingBranch {
    pub r: f64,
    pub parallax_rel: f64,
    pub meanmotion_rel: f64,
    pub ratio: f64,
    pub binds: Binds,
}

pub fn binding_branch(
    geom: &Geometry,
    r: f64,
    eps_rad: f64,
    mu: f64,
    stat: Stat,
) -> Result<BindingBranch, GridError> {
    let gamma = 1.0 / r;
    let par_rel = gamma_cell(parallax_lever(geom, stat), eps_rad)? / gamma;
    let t2 = meanmotion_lever(geom, stat);
    let a = 0.5 * mu * gamma * gamma * gamma * t2;
    if !(a > 0.0) {
        return Err(GridError::DegenerateGeometry("mean-motion lever is zero".into()));
    }
    let mm_rel = eps_rad / (3.0 * a);
    Ok(BindingBranch {
        r,
        parallax_rel: par_rel,
        meanmotion_rel: mm_rel,
        ratio: mm_rel / par_rel,
        binds: if mm_rel > par_rel { Binds::Parallax } else { Binds::MeanMotion },
    })
}

// ------------------------------------------------------------------------------------ the ladder

/// One `(r, rdot)` grid node, carrying the cell it was sized with.
///
/// Named `GridNode` and not `Node` because [`crate::spherical_pair::Node`] is the *asserted*
/// pair `(r, rdot)` the solve consumes: this is that plus the cell it stands for, and the two
/// are re-exported into the same namespace by `lib.rs`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridNode {
    pub r: f64,
    pub rdot: f64,
    /// Half-cell in gamma [AU^-1].
    pub dgamma: f64,
    /// Half-cell in rdot [AU/day].
    pub drdot: f64,
    /// Parallax lever [AU].
    pub w: f64,
    /// Zoom lever [rad day].
    pub z: f64,
}

impl GridNode {
    /// The node as the solve wants it.
    pub fn node(&self) -> crate::spherical_pair::Node {
        crate::spherical_pair::Node { r: self.r, rdot: self.rdot }
    }
}

/// Emit the `(r, rdot)` grid over `[r_min, r_max]`.
///
/// `geometry_at(r)` must return the [`Geometry`] for a trial object at distance `r` -- the
/// caller owns the ephemeris and the sky position. ⭐ Because both levers are ephemeris
/// quantities, a run sizes its grid from its **own** cadence and elongation coverage rather than
/// from the synthetic geometry the closed form was measured on.
///
/// Shells march outward in `gamma` from `1/r_max`, stepping by the full cell `2*dgamma`
/// evaluated **at the shell just placed** -- see the module docs for why the step is not fixed.
///
/// 🔴 The result is a ladder for **one sky region**. Call this per band and concatenate.
pub fn build_ladder<F>(
    geometry_at: F,
    r_min: f64,
    r_max: f64,
    eps_rad: f64,
    mu: f64,
    stat: Stat,
    max_nodes: usize,
) -> Result<Vec<GridNode>, GridError>
where
    F: Fn(f64) -> Result<Geometry, GridError>,
{
    if !(r_min < r_max) {
        return Err(GridError::Invalid(format!("need r_min < r_max; got {r_min} < {r_max}")));
    }

    let mut nodes: Vec<GridNode> = Vec::new();
    let mut gamma = 1.0 / r_max;
    let gamma_stop = 1.0 / r_min;
    while gamma < gamma_stop {
        let r = 1.0 / gamma;
        let geom = geometry_at(r)?;
        let w = parallax_lever(&geom, stat);
        let z = zoom_lever(&geom, stat);
        let dgamma = gamma_cell(w, eps_rad)?;
        let drdot = rdot_cell(r, z, eps_rad)?;

        let n = n_rdot_nodes(r, drdot, mu)?;
        // 🔴 A deliberate divergence from the Python prototype, in the error path only: it checks
        // `max_nodes` *after* emitting a whole shell, so a too-small `eps_rad` makes it try to
        // allocate one shell of ~1e10 nodes and die by OOM rather than by its own guard. Checking
        // the shell's width first cannot change any successful ladder.
        if n > max_nodes {
            return Err(GridError::LadderOverflow { max_nodes, r_min });
        }
        let span = rdot_span(r, mu);
        // Nodes at cell centres across the bound span, so the extremes are covered by half a
        // cell rather than sitting exactly on an edge.
        for k in 0..n {
            let rdot = -0.5 * span + span * (k as f64 + 0.5) / n as f64;
            nodes.push(GridNode { r, rdot, dgamma, drdot, w, z });
        }

        gamma += 2.0 * dgamma;
        if nodes.len() > max_nodes {
            return Err(GridError::LadderOverflow { max_nodes, r_min });
        }
    }
    Ok(nodes)
}

/// Shape of a ladder: shells, total nodes, and where the `rdot` axis collapses to one.
///
/// ⭐ The collapse distance is the useful one to quote -- beyond it the entire *bound* `rdot`
/// span fits inside a single cell, so the generator emits one `rdot = 0` node and stops. It is a
/// **span** criterion, so unlike the older "the object's own rdot" figures it assumes nothing
/// about a typical object.
#[derive(Debug, Clone, PartialEq)]
pub struct LadderSummary {
    pub shells: usize,
    pub nodes: usize,
    pub r_min: Option<f64>,
    pub r_max: Option<f64>,
    pub max_rdot_nodes: usize,
    pub single_rdot_shells: usize,
    pub collapse_au: Option<f64>,
}

pub fn summarize(nodes: &[GridNode]) -> LadderSummary {
    if nodes.is_empty() {
        return LadderSummary {
            shells: 0,
            nodes: 0,
            r_min: None,
            r_max: None,
            max_rdot_nodes: 0,
            single_rdot_shells: 0,
            collapse_au: None,
        };
    }

    // Shells in emission order (gamma strictly increasing => r strictly decreasing), counted by
    // exact `r` equality, which is what makes them the same shell in the first place.
    let mut shells: Vec<(f64, usize)> = Vec::new();
    for nd in nodes {
        match shells.iter_mut().find(|(r, _)| *r == nd.r) {
            Some((_, count)) => *count += 1,
            None => shells.push((nd.r, 1)),
        }
    }
    let mut ascending = shells.clone();
    ascending.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("NaN shell radius"));

    // The collapse is where the axis becomes single AND stays single outward: walk in from the
    // largest `r` while the count is one.
    let mut collapse = None;
    for (r, count) in ascending.iter().rev() {
        if *count != 1 {
            break;
        }
        collapse = Some(*r);
    }

    LadderSummary {
        shells: shells.len(),
        nodes: nodes.len(),
        r_min: Some(ascending.first().unwrap().0),
        r_max: Some(ascending.last().unwrap().0),
        max_rdot_nodes: shells.iter().map(|(_, c)| *c).max().unwrap(),
        single_rdot_shells: shells.iter().filter(|(_, c)| *c == 1).count(),
        collapse_au: collapse,
    }
}

// ------------------------------------------------------------------------------------- tests

#[cfg(test)]
mod tests {
    use super::*;
    use spacerocks::coordinates::Origin;

    fn mu() -> f64 {
        Origin::SSB.mu()
    }

    /// A hand-built window: a circular observer, a fixed line of sight, epochs that are NOT
    /// symmetric about the anchors (so a least-squares detrend and the anchored one differ).
    fn geom() -> Geometry {
        let t0 = 0.0;
        let t1 = 365.25;
        let obs = |t: f64| {
            let th = 2.0 * std::f64::consts::PI * t / 365.25;
            Vector3::new(th.cos(), th.sin(), 0.0)
        };
        let los = |t: f64| {
            // A slow drift, so `Z` is non-degenerate.
            let a = 0.20 + 1.0e-4 * t;
            Vector3::new(a.cos(), a.sin(), 0.05).normalize()
        };
        let ts: Vec<f64> = (0..24).map(|k| 500.0 + 40.0 * k as f64).collect();
        let es = ts.iter().map(|t| obs(*t)).collect();
        let us = ts.iter().map(|t| los(*t)).collect();
        Geometry::new(t0, t1, obs(t0), obs(t1), los(t0), los(t1), ts, es, us).unwrap()
    }

    #[test]
    fn the_projector_vanishes_exactly_at_both_anchors() {
        // 🔴 THE assertion of this module. Exactly zero, not "small": the anchored projector is
        // an interpolation through the two anchor values, so a port that reached for the
        // least-squares `{1, t}` detrend of `tangent_v2.rs:515-520` fails here and nowhere else
        // in an obvious way -- it would merely size every cell slightly wrong, forever.
        let g = geom();
        let v0 = Vector3::new(0.3, -0.7, 0.1);
        let v1 = Vector3::new(-1.1, 0.4, 2.0);
        let ts = vec![g.t0, g.t1, 123.0];
        let vals = vec![v0, v1, Vector3::new(5.0, 5.0, 5.0)];
        let d = anchored_residual(&vals, &ts, g.t0, g.t1, v0, v1);
        assert_eq!(d[0], 0.0);
        // 🔴 The one at `t1` is what forced the endpoint form of the interpolation: with the
        // prototype's `v0 + s*(v1 - v0)` this lands at 1.1e-16, and the Python suite has to
        // assert the invariant with a tolerance instead of an equality.
        assert_eq!(d[1], 0.0);
        assert!(d[2] > 0.0, "an off-anchor epoch must not be absorbed");
    }

    #[test]
    fn an_affine_signature_is_absorbed_entirely() {
        // The physical content of the projector: anything affine in `t` costs the cell nothing,
        // because the solve reproduces it by re-solving `h`.
        let t0 = 10.0;
        let t1 = 400.0;
        let a = Vector3::new(1.0, -2.0, 0.5);
        let b = Vector3::new(0.01, 0.02, -0.03);
        let ts: Vec<f64> = (0..7).map(|k| 50.0 + 37.0 * k as f64).collect();
        let vals: Vec<Vector3<f64>> = ts.iter().map(|t| a + b * *t).collect();
        let d = anchored_residual(&vals, &ts, t0, t1, a + b * t0, a + b * t1);
        // Relative to the signal, not absolute: the anchors are pinned exactly, but an off-anchor
        // `a + b*t` is itself only accurate to an ulp of its own magnitude (~10 here), so the
        // absorbed residual floors at ~1e-16 of that and cannot be asserted as a bare zero.
        let scale = vals.iter().map(|v| v.norm()).fold(0.0, f64::max);
        for x in d {
            assert!(x < 1e-14 * scale, "affine part not absorbed: {x}");
        }
    }

    #[test]
    fn the_median_of_an_even_sample_is_the_mean_of_the_two_middles() {
        // Not pedantry: every fixture window has 24 epochs, so this is the branch the whole
        // gate runs through. Taking the lower middle instead would shift every lever.
        let v = [4.0, 1.0, 3.0, 2.0];
        assert_eq!(stat_of(&v, Stat::Median), 2.5);
        assert_eq!(stat_of(&[3.0, 1.0, 2.0], Stat::Median), 2.0);
        // numpy's default linear interpolation, not a nearest-rank percentile.
        assert!((stat_of(&[0.0, 1.0, 2.0, 3.0], Stat::P90) - 2.7).abs() < 1e-12);
    }

    #[test]
    fn the_parallax_lever_carries_no_r_and_the_zoom_lever_does() {
        // `W` is pure observer geometry -- that is *why* the grid is uniform in 1/r. `Z` is not:
        // it mixes the parallactic reflex with the object's own motion, which is why it must be
        // re-measured per shell instead of carrying a fitted exponent.
        let g = geom();
        let w = parallax_lever(&g, Stat::Median);
        let z = zoom_lever(&g, Stat::Median);
        assert!(w > 0.0 && z > 0.0);
        // The gamma cell is set by W alone: same geometry, same cell, whatever r we ask about.
        let c1 = gamma_cell(w, 1e-6).unwrap();
        let c2 = gamma_cell(parallax_lever(&g, Stat::Median), 1e-6).unwrap();
        assert_eq!(c1, c2);
        // The rdot cell scales linearly with r at fixed geometry.
        let a = rdot_cell(40.0, z, 1e-6).unwrap();
        let b = rdot_cell(80.0, z, 1e-6).unwrap();
        assert!((b / a - 2.0).abs() < 1e-12);
    }

    #[test]
    fn a_window_that_is_only_the_anchors_is_rejected_not_infinite() {
        // The failure mode this error type exists for: a caller that accepted an infinite cell
        // would emit a one-node grid and report full coverage of the sky.
        let g = geom();
        let ts = vec![g.t0, g.t1];
        let es = vec![g.e0, g.e1];
        let us = vec![g.u0, g.u1];
        let degenerate = Geometry::new(g.t0, g.t1, g.e0, g.e1, g.u0, g.u1, ts, es, us).unwrap();
        let w = parallax_lever(&degenerate, Stat::Median);
        assert_eq!(w, 0.0);
        assert!(matches!(gamma_cell(w, 1e-6), Err(GridError::DegenerateGeometry(_))));
        assert!(matches!(
            zoom_lever(&degenerate, Stat::Median),
            x if x == 0.0
        ));
        assert!(matches!(rdot_cell(45.0, 0.0, 1e-6), Err(GridError::DegenerateGeometry(_))));
    }

    #[test]
    fn parallax_binds_over_the_tno_range() {
        // The branch decision, checked rather than assumed. `ratio > 1` means the mean-motion
        // tolerance is the looser one, so parallax sets the cell and the grid variable is 1/r.
        let g = geom();
        for r in [40.0, 100.0, 400.0, 1600.0] {
            let b = binding_branch(&g, r, 4.848e-5, mu(), Stat::Median).unwrap();
            assert_eq!(b.binds, Binds::Parallax, "mean motion bound at r={r}");
            assert!(b.ratio > 1.0);
        }
        // ...and it is *increasingly* the binding one farther out: the mean-motion signal dies
        // as gamma^3 while the parallax lever does not move.
        let near = binding_branch(&g, 40.0, 4.848e-5, mu(), Stat::Median).unwrap();
        let far = binding_branch(&g, 1600.0, 4.848e-5, mu(), Stat::Median).unwrap();
        assert!(far.ratio > near.ratio);
    }

    #[test]
    fn the_rdot_axis_collapses_to_one_node_far_out_and_the_ladder_says_where() {
        let g = geom();
        let ladder = build_ladder(|_r| Ok(g.clone()), 200.0, 1600.0, 4.848e-5, mu(), Stat::Median, 100_000)
            .unwrap();
        let s = summarize(&ladder);
        assert!(s.shells > 1 && s.nodes >= s.shells);
        // Shells are placed in gamma, so they come out from far to near, strictly decreasing.
        assert!(ladder.windows(2).all(|w| w[0].r >= w[1].r));
        // The outermost shell is exactly r_max; the ladder stops before r_min.
        assert_eq!(s.r_max.unwrap(), 1600.0);
        assert!(s.r_min.unwrap() >= 200.0);
        // The collapse is a suffix property: every shell at or beyond it has one node.
        let collapse = s.collapse_au.expect("the rdot axis must collapse somewhere out here");
        for nd in &ladder {
            if nd.r > collapse {
                assert_eq!(
                    ladder.iter().filter(|m| m.r == nd.r).count(),
                    1,
                    "shell at r={} beyond the collapse has more than one rdot node",
                    nd.r
                );
            }
        }
    }

    #[test]
    fn a_single_rdot_shell_sits_at_zero_and_a_multi_shell_straddles_it() {
        // Nodes are at cell CENTRES across the bound span, so a collapsed shell lands exactly on
        // rdot = 0 rather than at an edge -- and a shell with an even count straddles it
        // symmetrically. An off-by-one in the centring shows up here and nowhere else.
        let g = geom();
        let ladder =
            build_ladder(|_r| Ok(g.clone()), 200.0, 1600.0, 4.848e-5, mu(), Stat::Median, 100_000).unwrap();
        let outer: Vec<&GridNode> = ladder.iter().filter(|n| n.r == 1600.0).collect();
        assert_eq!(outer.len(), 1);
        assert_eq!(outer[0].rdot, 0.0);

        let inner_r = ladder.last().unwrap().r;
        let inner: Vec<&GridNode> = ladder.iter().filter(|n| n.r == inner_r).collect();
        if inner.len() > 1 {
            let sum: f64 = inner.iter().map(|n| n.rdot).sum();
            assert!(sum.abs() < 1e-15, "rdot nodes are not symmetric about zero: {sum}");
            let span = rdot_span(inner_r, mu());
            assert!(inner.iter().all(|n| n.rdot.abs() <= 0.5 * span));
        }
    }

    #[test]
    fn the_ladder_steps_by_the_cell_it_just_placed() {
        // The re-evaluation that keeps the grid honest where W drifts: consecutive shells are
        // one FULL cell apart in gamma, measured at the inner one.
        let g = geom();
        let ladder =
            build_ladder(|_r| Ok(g.clone()), 200.0, 1600.0, 4.848e-5, mu(), Stat::Median, 100_000).unwrap();
        let mut shells: Vec<f64> = Vec::new();
        for nd in &ladder {
            if shells.last() != Some(&nd.r) {
                shells.push(nd.r);
            }
        }
        for w in shells.windows(2) {
            let dgamma = ladder.iter().find(|n| n.r == w[0]).unwrap().dgamma;
            assert!(((1.0 / w[1] - 1.0 / w[0]) - 2.0 * dgamma).abs() < 1e-15 * (1.0 / w[1]));
        }
    }

    #[test]
    fn a_runaway_ladder_errors_rather_than_returning_a_giant_grid() {
        // Far too many shells for the cap: the guard fires on the accumulated count.
        let g = geom();
        let err = build_ladder(|_r| Ok(g.clone()), 5.0, 1600.0, 4.848e-5, mu(), Stat::Median, 50)
            .unwrap_err();
        assert!(matches!(err, GridError::LadderOverflow { .. }));
    }

    #[test]
    fn one_absurdly_wide_shell_errors_before_it_is_allocated() {
        // The other half of the same guard, and the reason it is checked before the inner loop:
        // a tiny eps_rad makes a SINGLE shell ask for ~1e10 rdot nodes, which the prototype would
        // try to allocate before its own post-hoc check ever ran.
        let g = geom();
        let err = build_ladder(|_r| Ok(g.clone()), 1000.0, 1600.0, 1e-14, mu(), Stat::Median, 10_000)
            .unwrap_err();
        assert!(matches!(err, GridError::LadderOverflow { .. }));
    }

    #[test]
    fn an_empty_ladder_summarizes_to_nothing_rather_than_panicking() {
        let s = summarize(&[]);
        assert_eq!(s.nodes, 0);
        assert_eq!(s.collapse_au, None);
    }
}
