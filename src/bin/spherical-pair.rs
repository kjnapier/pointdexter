//! Stage 0: the C2 spherical-pair binary -- config schema, validation, and a runnable self-check.
//!
//! The build plan's gate for this stage is one sentence: **every C++ hardcode comes from config,
//! never transcribed.** The reference linker carries `t_ref = 2456635.5`, `time_diff = 120`,
//! `ang_thresh = 40`, `chi2_thresh = 2.0` and a `5..1e5` AU band baked into its source; ssolink's
//! own working principle is the same one ("no absolute paths in code; read from `config/`"). So
//! this file defines the schema and refuses to run without it, rather than defaulting its way to a
//! plausible answer.
//!
//! Three choices worth stating, because each is the opposite of what a convenience default does:
//!
//! * 🔴 **No `serde(default)` on a physical parameter.** A missing `k_sigma` must be an error, not
//!   a silent 0 -- and 0 is precisely the value that reproduces the *uncorrected* gate whose
//!   omission put `rho` low by 2.6-6.3x, entering the pair count squared. A default here would
//!   reintroduce a measured bug quietly.
//! * 🔴 **`deny_unknown_fields`.** A misspelled key that is ignored is a config that does not
//!   describe the run, which is how a note ends up citing parameters the run never used.
//! * 🔴 **`mu` is named, not numbered.** The config names an origin and the code calls
//!   `Origin::SSB.mu()`; it never carries the constant. The `/3` frame move showed a transcribed
//!   literal 1.7e-11 away from what the port CALLS still fails the gate by ~0.067%.
//!
//! What this does NOT do yet: load detections or run a search. Stage 3 owns the data layer, and
//! wiring it blind would be the "finished note citing another run's output" failure in code form.
//! `--self-check` exercises the projection, the index and the gate on synthetic geometry, so the
//! stage is runnable and provably wired before any of that exists.

use std::path::PathBuf;

use clap::Parser;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use spacerocks::coordinates::Origin;

use spacerocks::SpiceKernel;
use spacerocks::transforms::solve_for_universal_anomaly;

use pointdexter::io::load_detections::load_detections;
use pointdexter::spherical_pair::{
    H_SCAN_HI_FRAC, H_SCAN_LO_FRAC, Node, Pair, Solution, gate_radius, h_max, h_scan_grid,
    range_quadratic,
};
use pointdexter::spherical_pair_anchor::{Anchor, anchor_pairs_and_solve, build_anchors};
use pointdexter::spherical_pair_extend::{ExtendParams, VisitIndex, extend_candidate};
use pointdexter::spherical_pair_index::{
    BaryIndex, angle_of_chord, chord_of_angle, gate_radius_astrometric, pair_within_gate,
};
use pointdexter::spherical_pair_grid::{
    Geometry, GridError, Stat, build_ladder, summarize,
};
use pointdexter::spherical_pair_load::{ObsSet, SigmaSource};

const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

#[derive(Parser, Debug)]
#[command(name = "spherical-pair", about = "C2 spherical two-point search (stage 0 skeleton)")]
struct Cli {
    /// YAML config. Required: this binary has no working defaults by design.
    #[arg(short, long)]
    config: PathBuf,

    /// Print the resolved configuration and exit.
    #[arg(long)]
    print_config: bool,

    /// Exercise projection, index and gate on synthetic geometry, and report.
    #[arg(long)]
    self_check: bool,

    /// Load the configured detections and run stages 2-3: anchors, anchor pairs, solve, chi2 vet.
    #[arg(long)]
    run: bool,

    /// Size the (r, rdot) ladder from this run's geometry, print its shape, and exit. Sizing is
    /// cheap; the search that follows is not, so the grid a config implies should be inspectable
    /// without paying for a run to find out.
    #[arg(long)]
    ladder_only: bool,

    /// DIAGNOSTIC: over the same gated pairs a `--run` would solve, count how many brackets the
    /// `F(h)` scan contains instead of taking the first. `solve_h` returns the smallest-h root,
    /// so a lookup grid of initial guesses is only safe where that root is the only one. Writes
    /// no candidates and solves nothing.
    #[arg(long)]
    sign_census: bool,
}

/// Which origin the frame and `mu` are taken from. Named, never numeric.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum OriginName {
    Ssb,
    Sun,
}

impl OriginName {
    fn mu(self) -> f64 {
        match self {
            // 🔴 Origin::SSB.mu(), not spacerocks::constants::MU_BARY: the two differ by 1.7e-11
            // and the port must match what it CALLS.
            OriginName::Ssb => Origin::SSB.mu(),
            OriginName::Sun => Origin::SUN.mu(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Frame {
    /// MJH: for TNO work the barycentre is the appropriate origin.
    origin: OriginName,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Search {
    /// Replaces the reference's baked `5 .. 1e5` AU band.
    r_min_au: f64,
    r_max_au: f64,
    /// Grid budget. Visibly a choice, not a derivation: a p90 variant moves node counts 1.6-3x
    /// and changes no scaling law.
    eps_arcsec: f64,
    /// Extension window. Reach is the cheapest knob in the method -- 1 yr needs ~103 nodes where
    /// 3 yr needs ~304 -- so it is configuration, not a constant.
    reach_days: f64,
    max_nodes: usize,
    /// Which statistic of the anchored residual sets the lever: "median" or "p90".
    ///
    /// 🔴 A POLICY, not a fact -- how much of the extension window the cell must hold to `eps`
    /// (`spherical_pair_grid::Stat`). It is in the config precisely because hardcoding it makes
    /// the policy invisible, and the residual is identically ZERO at both anchors: a median over
    /// a window whose epochs cluster near an anchor reports the lever far too small, the cell far
    /// too large, and the grid quietly stops covering.
    lever_stat: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Gate {
    /// Sigma multiple on the astrometric term. 🔴 No default: 0 silently reproduces the
    /// uncorrected gate.
    k_sigma: f64,
    /// Used only where a detection carries no reported uncertainty. 🔴 Rubin DOES report
    /// per-source sigma; this is a fallback, and a run that leans on it should say so.
    fallback_astrom_sigma_arcsec: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct AnchorCfg {
    /// An anchor is a PAIR OF DETECTIONS on one night (MJH, 2026-08-15).
    max_intra_night_hours: f64,
    /// Replaces the reference's `time_diff = 120`.
    max_anchor_baseline_days: f64,
    /// Replaces the reference's `chi2_thresh = 2.0`. 9.488 is 95% of chi2_4, the measured gate.
    chi2_max: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Extension {
    /// Gather tolerance. 🔴 This is NOT free: it is set by the grid cell, because a node offset of
    /// half a cell already displaces the prediction by ~4" over a year-long span. Astrometry alone
    /// would allow ~0.15".
    tolerance_arcsec: f64,
    /// Minimum supporting detections. 🔴 A raw count does not travel between fields -- chance
    /// support scales as n_opportunities * rho * pi * tol^2 -- so this is a floor, and the score
    /// below is what actually discriminates.
    min_support: usize,
    /// Reject a candidate whose support is no better than chance at this probability.
    max_chance_probability: f64,
    /// `var/mean` of the chance support count. 🔴 **1.0 is Poisson and is optimistic in the tail by
    /// up to ~4x** -- measured 1.3 at a 3" gather rising to 9.2 at 30", because detections cluster.
    /// No default: a silent 1.0 would make every score look better than it is.
    chance_dispersion: f64,
    /// 🔴 Support on the anchors' own nights is nearly free (it confirms the detection, not the
    /// orbit) and the shift control endorses it. Keep true unless deliberately measuring that.
    exclude_anchor_nights: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Io {
    detections: PathBuf,
    output: PathBuf,
    /// Directory holding the SPK/BPC kernels the detection loader needs to place the observer
    /// from an `obscode`. Optional: a catalogue that supplies `obs_x/obs_y/obs_z` directly needs
    /// no kernel, and loading one that is not there would fail a run that does not use it.
    /// 🔴 That path is also the one with NO frame check -- see `observer_positions_are_barycentric`.
    spice_path: Option<PathBuf>,
    /// "ecliptic" or "equatorial". Either works -- C2's geometry is vectors -- provided `rho_hat`
    /// and `observer_position` are on the SAME plane, which the loader guarantees by rotating both.
    reference_plane: String,
}

/// Which detection column is authoritative for the astrometric sigma, and how to read it.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum SigmaSourceName {
    /// `ast_ucty`; the loader has already converted it from arcsec to radians.
    Ast,
    /// `ra_ucty`/`dec_ucty`; the loader leaves these RAW, so they are read as arcsec here.
    Radec,
    /// A single configured value for every detection.
    Fixed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Data {
    sigma_source: SigmaSourceName,
    /// Only consulted for `radec`. 🔴 Not guessable from the data, and worth a factor of 2 on one
    /// axis at |dec| = 60: does `ra_ucty` already carry the cos(dec) factor?
    sigma_ra_includes_cos_dec: bool,
    /// Only consulted for `fixed`.
    sigma_fixed_arcsec: f64,
    /// Epoch tolerance for grouping detections into one visit. An exact float epoch is not a key.
    visit_tolerance_seconds: f64,
    /// Night boundary offset from the JD tick; 0.5 is the `floor(mjd - 0.5)` noon convention.
    /// ⚠️ A fallback. Where the input carries the survey's own `dayObs`, group on that instead.
    night_boundary_days: f64,
    // 🔴 `sigma_cap_percentile` was REMOVED 2026-08-16. The shared KD query is now sized on
    // `sqrt(2) * max(sigma)`, a strict bound on every pair's quadrature, so there is no
    // percentile left to choose. The field is deliberately NOT accepted-and-ignored: this
    // struct is `deny_unknown_fields`, so an old config fails loudly rather than describing a
    // run by a parameter it did not use. See ObsSet::sigma_cap.
    /// 🔴 Required acknowledgement. The loader's `obs_x/obs_y/obs_z` path takes the observer
    /// position from the file with NO frame check; C2 needs it BARYCENTRIC, because mu and the
    /// origin move together. Setting this false refuses to run rather than producing a silently
    /// heliocentric search.
    observer_positions_are_barycentric: bool,
    /// OPTIONAL override. Leave empty (the default) and the run sizes its own `(r, rdot)` ladder
    /// from its OWN cadence and elongation coverage via `spherical_pair_grid::build_ladder`.
    ///
    /// 🔴 An explicit list is a DIAGNOSTIC, not a search. Placing nodes by hand around an answer
    /// you already know measures whether the chain can recover an object it was pointed at; it
    /// does not measure whether the chain can find one. Both are useful and they are not the
    /// same claim, so a run that sets this says so loudly in its banner.
    #[serde(default)]
    explicit_nodes: Vec<NodeSpec>,
}

/// The [`Geometry`] the ladder sizes itself from: this run's own cadence and pointing.
///
/// Three choices, each of which changes the cell and so is stated rather than buried:
///
/// * **`t0`/`t1` are the WIDEST anchor baseline the run can form**, not a typical one. The
///   projector is the affine function through the two anchors and is identically zero there, so
///   the unabsorbed residual -- and with it the lever -- grows with the baseline. Sizing on a
///   short baseline would emit a grid too coarse for the long pairs the same run will search.
///   The widest baseline gives the tightest cell, which is the conservative direction.
/// * **The line of sight is the region's mean `rho_hat`**, not a per-epoch object track. `W` is
///   direction-dependent but orbit-free, and the exported region is ~1 deg across, so one
///   direction sizes the whole region. 🔴 This is what makes the ladder per-REGION; a run
///   spanning many degrees of ecliptic latitude must ladder per band and concatenate.
/// * **One entry per VISIT**, not per detection: the extension window is a cadence, and
///   weighting it by detection count would let one crowded visit set the grid.
struct RunCadence {
    ts: Vec<f64>,
    es: Vec<Vector3<f64>>,
    u: Vector3<f64>,
    t0: f64,
    t1: f64,
    e0: Vector3<f64>,
    e1: Vector3<f64>,
}

fn run_cadence(
    set: &ObsSet,
    visits: &[Vec<usize>],
    night_keys: &[i64],
    nights: &std::collections::BTreeMap<i64, Vec<usize>>,
) -> Result<RunCadence, String> {
    if night_keys.len() < 2 {
        return Err("ladder needs at least two nights to define an anchor baseline".into());
    }
    let mut u_sum = Vector3::zeros();
    let (mut ts, mut es) = (Vec::new(), Vec::new());
    for v in visits.iter() {
        let o = &set.obs[v[0]];
        ts.push(o.epoch);
        es.push(o.observer);
        u_sum += v.iter().fold(Vector3::zeros(), |a, &i| a + set.obs[i].rho_hat);
    }
    let first = nights[&night_keys[0]][0];
    let last = *nights[night_keys.last().expect("night_keys non-empty")]
        .last()
        .expect("a night has visits");
    let (a, b) = (&set.obs[visits[first][0]], &set.obs[visits[last][0]]);
    Ok(RunCadence {
        ts,
        es,
        u: u_sum.normalize(),
        t0: a.epoch,
        t1: b.epoch,
        e0: a.observer,
        e1: b.observer,
    })
}

/// The sky track of ONE trial orbit at distance `r`, over this run's cadence.
///
/// `vr`/`vo` are the radial and transverse speeds at `t0`; the transverse direction is `az`
/// radians around the line of sight. Propagation is the crate's own f-and-g
/// ([`InitialCondition3D`]), not a series expansion -- the lever is precisely the part of the
/// signature that SURVIVES removing an affine, so a quadratic-truncated track would approximate
/// the very term being measured.
fn trial_geometry(cad: &RunCadence, r: f64, vr: f64, vo: f64, az: f64, mu: f64)
    -> Result<Geometry, String>
{
    let p0 = range_quadratic(&cad.e0, &cad.u, r)
        .ok_or_else(|| format!("no line-of-sight point at r={r}"))?;
    let rhat = p0.normalize();
    // a transverse basis perpendicular to the heliocentric radius
    let mut a = rhat.cross(&Vector3::z());
    if a.norm() < 1e-8 {
        a = rhat.cross(&Vector3::x());
    }
    let a = a.normalize();
    let b = rhat.cross(&a);
    let that = a * az.cos() + b * az.sin();
    let v0 = rhat * vr + that * vo;

    let ic = pointdexter::initial_condition_3d::InitialCondition3D::from_spherical(
        "trial".into(), r, vr, vo, spacerocks::Time::new(cad.t0, "tdb", "jd")
            .map_err(|e| format!("trial epoch: {e}"))?, mu,
    ).map_err(|e| format!("trial orbit: {e}"))?;

    let mut us = Vec::with_capacity(cad.ts.len());
    for (t, e) in cad.ts.iter().zip(cad.es.iter()) {
        let (f, g) = ic.fg_at_epoch(*t);
        let p = p0 * f + v0 * g;
        let geo = p - e;
        let n = geo.norm();
        if !(n > 0.0) {
            return Err("trial track passed through the observer".into());
        }
        us.push(geo / n);
    }
    let u_at = |t: f64, e: &Vector3<f64>| {
        let (f, g) = ic.fg_at_epoch(t);
        (p0 * f + v0 * g - e).normalize()
    };
    Geometry::new(
        cad.t0, cad.t1, cad.e0, cad.e1,
        u_at(cad.t0, &cad.e0), u_at(cad.t1, &cad.e1),
        cad.ts.clone(), cad.es.clone(), us,
    )
    .map_err(|e: GridError| format!("trial geometry: {e:?}"))
}

/// The trial geometry that sizes shell `r`: the TIGHTEST cell anywhere in the bound domain.
///
/// 🔴 MJH, 2026-08-16: these runs must cover orbits up to the bound limit. `rdot_span` already
/// makes the ṙ *range* the full escape width `2*sqrt(2mu/r)`. This makes the ṙ *cell* match: the
/// zoom lever Z depends on `r` through the sky track and "must be measured per shell, not
/// hardcoded as an exponent" (`spherical_pair_grid::zoom_lever`), and a track is only fixed once
/// the velocity is. So sample the bound velocity domain `vr^2 + vo^2 <= 2mu/r` and keep the
/// sample with the LARGEST Z, i.e. the smallest `drdot`. Sizing on any interior sample would
/// emit a grid too coarse for part of the span it claims to cover.
fn geometry_for_shell(cad: &RunCadence, r: f64, mu: f64) -> Result<Geometry, GridError> {
    let v_esc = (2.0 * mu / r).sqrt();
    let mut best: Option<(f64, Geometry)> = None;
    for &fv in &[0.25_f64, 0.5, 0.75, 1.0] {
        for &fr in &[0.0_f64, 0.5, 0.9] {
            let vr = fr * fv * v_esc;
            let vo = ((fv * v_esc).powi(2) - vr * vr).max(0.0).sqrt();
            for k in 0..4 {
                let az = std::f64::consts::FRAC_PI_2 * k as f64;
                if let Ok(g) = trial_geometry(cad, r, vr, vo, az, mu) {
                    let z = pointdexter::spherical_pair_grid::zoom_lever(&g, Stat::Median).abs();
                    if z.is_finite() && best.as_ref().is_none_or(|(bz, _)| z > *bz) {
                        best = Some((z, g));
                    }
                }
            }
        }
    }
    best.map(|(_, g)| g).ok_or_else(|| {
        GridError::DegenerateGeometry(format!("no usable trial track at r={r}"))
    })
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct NodeSpec {
    r_au: f64,
    rdot_au_per_day: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Config {
    frame: Frame,
    search: Search,
    gate: Gate,
    anchor: AnchorCfg,
    extension: Extension,
    data: Data,
    io: Io,
}

impl Config {
    /// Reject what cannot mean anything, with the offending value named. A validator that only
    /// checks presence lets `r_min > r_max` through to produce an empty ladder that reads as
    /// "no objects here".
    fn validate(&self) -> Result<(), String> {
        let pos = |name: &str, v: f64| -> Result<(), String> {
            if v.is_finite() && v > 0.0 {
                Ok(())
            } else {
                Err(format!("{name} must be finite and positive; got {v}"))
            }
        };
        pos("search.r_min_au", self.search.r_min_au)?;
        pos("search.r_max_au", self.search.r_max_au)?;
        pos("search.eps_arcsec", self.search.eps_arcsec)?;
        pos("search.reach_days", self.search.reach_days)?;
        pos("anchor.max_intra_night_hours", self.anchor.max_intra_night_hours)?;
        pos("anchor.max_anchor_baseline_days", self.anchor.max_anchor_baseline_days)?;
        pos("anchor.chi2_max", self.anchor.chi2_max)?;
        pos("extension.tolerance_arcsec", self.extension.tolerance_arcsec)?;
        pos("gate.fallback_astrom_sigma_arcsec", self.gate.fallback_astrom_sigma_arcsec)?;
        if !(self.search.r_min_au < self.search.r_max_au) {
            return Err(format!(
                "search.r_min_au must be below r_max_au; got {} .. {}",
                self.search.r_min_au, self.search.r_max_au
            ));
        }
        if !(self.gate.k_sigma >= 0.0) || !self.gate.k_sigma.is_finite() {
            return Err(format!("gate.k_sigma must be finite and >= 0; got {}", self.gate.k_sigma));
        }
        if self.search.max_nodes == 0 {
            return Err("search.max_nodes must be nonzero".into());
        }
        let p = self.extension.max_chance_probability;
        if !(p > 0.0 && p < 1.0) {
            return Err(format!("extension.max_chance_probability must be in (0,1); got {p}"));
        }
        if !(self.extension.chance_dispersion >= 1.0) || !self.extension.chance_dispersion.is_finite()
        {
            return Err(format!(
                "extension.chance_dispersion must be finite and >= 1 (1 = Poisson); got {}. \
                 Below 1 would claim chance support is UNDER-dispersed, which nothing measured.",
                self.extension.chance_dispersion
            ));
        }
        pos("data.visit_tolerance_seconds", self.data.visit_tolerance_seconds)?;
        pos("data.sigma_fixed_arcsec", self.data.sigma_fixed_arcsec)?;
        if !self.data.observer_positions_are_barycentric {
            return Err(
                "data.observer_positions_are_barycentric is false: C2 asserts Origin::SSB.mu(), \
                 so a heliocentric observer is a frame error, not a preference. Refusing to run."
                    .into(),
            );
        }
        let plane = self.io.reference_plane.to_lowercase();
        if plane != "ecliptic" && plane != "equatorial" {
            return Err(format!("io.reference_plane must be ecliptic or equatorial; got {plane}"));
        }
        for (i, n) in self.data.explicit_nodes.iter().enumerate() {
            pos(&format!("data.explicit_nodes[{i}].r_au"), n.r_au)?;
            if !n.rdot_au_per_day.is_finite() {
                return Err(format!("data.explicit_nodes[{i}].rdot_au_per_day is not finite"));
            }
            if n.r_au < self.search.r_min_au || n.r_au > self.search.r_max_au {
                return Err(format!(
                    "data.explicit_nodes[{i}].r_au = {} lies outside the configured band {} .. {}",
                    n.r_au, self.search.r_min_au, self.search.r_max_au
                ));
            }
        }
        Ok(())
    }
}

/// Projection + index + gate on synthetic geometry: no ephemeris, no detections, no network.
///
/// The point is not the numbers -- they are synthetic -- but that the stage is wired: a
/// projection that silently dropped everything, or a gate that returned zero, shows here.
fn self_check(cfg: &Config) {
    let mu = cfg.frame.origin.mu();
    println!("\n== self-check (synthetic geometry; not a science result) ==");
    println!("  mu({:?}) = {:.16e}", cfg.frame.origin, mu);

    // an observer 1 AU out, a patch of sky 1 degree across
    let observer = Vector3::new(0.4, -0.91, 0.0);
    let centre = Vector3::new(0.15, 0.98, 0.06).normalize();
    let cap = 1.0_f64.to_radians();
    let omega = 2.0 * std::f64::consts::PI * (1.0 - cap.cos());

    let e1 = centre.cross(&Vector3::z()).normalize();
    let e2 = centre.cross(&e1);
    let n = 4000usize;
    let mut s: u64 = 0x5EED;
    let mut next = || {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((s >> 11) as f64) / ((1u64 << 53) as f64)
    };
    let mut los = Vec::with_capacity(n);
    for _ in 0..n {
        let cos_t = 1.0 - next() * (1.0 - cap.cos());
        let sin_t = (1.0 - cos_t * cos_t).max(0.0).sqrt();
        let phi = 2.0 * std::f64::consts::PI * next();
        los.push((cos_t * centre + sin_t * (phi.cos() * e1 + phi.sin() * e2)).normalize());
    }
    let obs = vec![observer; n];
    let ids: Vec<u32> = (0..n as u32).collect();

    // 🔴 sigma_quad is TWO detections summed in quadrature, not one detection's sigma. Passing a
    // single sigma understates the astrometric term by sqrt(2), and that term is 59% of the gate
    // at 40 AU -- an error that would read as a slightly tight gate rather than as a unit slip.
    let sigma_quad = cfg.gate.fallback_astrom_sigma_arcsec * ARCSEC * std::f64::consts::SQRT_2;

    // 🔴 Two baselines, because they are different questions. `max_intra_night_hours` is the
    // widest span the config ALLOWS; Rubin's actual revisit is ~0.6 h (median intra-night span
    // 0.603 h over 263 nights), and it is the revisit that sets the tracklet gate. Reporting only
    // the maximum makes the astrometric floor look negligible when it dominates in practice.
    let baselines = [("rubin revisit", 0.603 / 24.0), ("config max", cfg.anchor.max_intra_night_hours / 24.0)];

    // 🔴 Pair two DISJOINT halves. Pairing the index against itself counts every point as its own
    // neighbour, which inflated an early version of this report by ~40x and looked like clustering.
    let half = los.len() / 2;
    let (obs_a, obs_b) = obs.split_at(half);
    let (los_a, los_b) = los.split_at(half);
    let ids_a: Vec<u32> = (0..half as u32).collect();
    let ids_b: Vec<u32> = (half as u32..los.len() as u32).collect();
    let omega_half = omega; // both halves are drawn from the same cap

    println!(
        "  {:>14} {:>8} {:>9} {:>12} {:>12} {:>9} {:>9} {:>9}",
        "baseline", "r [AU]", "indexed", "gate k=0 [\"]", "gate k=K [\"]", "astrom %", "pairs", "meas/pred"
    );
    for (label, dt_intra) in baselines {
        for &r in &[cfg.search.r_min_au, cfg.search.r_max_au] {
            let a = BaryIndex::build(obs_a, los_a, &ids_a, r);
            let b = BaryIndex::build(obs_b, los_b, &ids_b, r);
            let bare = gate_radius(r, dt_intra, mu);
            let full = gate_radius_astrometric(r, dt_intra, mu, sigma_quad, cfg.gate.k_sigma);
            let ps = pair_within_gate(&a, &b, full, Some(omega_half));
            println!(
                "  {:>14} {:>8.1} {:>9} {:>12.3} {:>12.3} {:>8.1}% {:>9} {:>9.3}",
                label,
                r,
                a.len() + b.len(),
                bare / ARCSEC,
                full / ARCSEC,
                100.0 * (1.0 - (bare / full).powi(2)),
                ps.pairs.len(),
                ps.excess()
            );
            if a.report.unreachable + b.report.unreachable > 0 {
                println!(
                    "    note: {} of {} lines of sight cannot reach {r} AU",
                    a.report.unreachable + b.report.unreachable,
                    a.report.offered + b.report.offered
                );
            }
        }
    }

    // chord <-> angle is the identity the whole "no tangent plane" argument rests on
    let probe = cfg.extension.tolerance_arcsec * ARCSEC;
    let round = angle_of_chord(chord_of_angle(probe));
    println!(
        "  chord round trip at the extension tolerance: {:.6e} rad -> {:.6e} rad (drift {:.1e})",
        probe,
        round,
        (round - probe).abs()
    );
}

/// Load detections, group them, and run stage 3 for each configured node.
fn run_search(
    cfg: &Config,
    ladder_only: bool,
    sign_census: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let mu = cfg.frame.origin.mu();
    let mut kernel = SpiceKernel::new();
    if let Some(sp) = &cfg.io.spice_path {
        let sp = sp.display();
        kernel.load_spk(&format!("{sp}/sb441-n16.bsp"))?;
        kernel.load_spk(&format!("{sp}/de440s.bsp"))?;
        kernel.load_bpc(&format!("{sp}/earth_1962_240827_2124_combined.bpc"))?;
    } else {
        println!("no io.spice_path: the catalogue must carry obs_x/obs_y/obs_z");
    }

    let path = cfg.io.detections.to_str().ok_or("io.detections is not valid UTF-8")?;
    let dets = load_detections(path, &cfg.io.reference_plane, &kernel)?;

    let source = match cfg.data.sigma_source {
        SigmaSourceName::Ast => SigmaSource::Ast,
        SigmaSourceName::Radec => SigmaSource::RaDec {
            ra_includes_cos_dec: cfg.data.sigma_ra_includes_cos_dec,
        },
        SigmaSourceName::Fixed => SigmaSource::Fixed { arcsec: cfg.data.sigma_fixed_arcsec },
    };
    let set = ObsSet::from_detections(&dets, source);
    let r = set.report;
    println!(
        "loaded {} detections -> {} usable ({:?}); dropped {} without a sigma, {} on geometry",
        r.offered, r.converted, cfg.data.sigma_source, r.no_sigma, r.bad_geometry
    );
    if r.converted == 0 {
        return Err("no detection carried a usable sigma from the configured source".into());
    }
    // 🔴 A silent majority-drop is the failure this reports rather than survives: it looks
    // identical to a sparse field once the counts reach the search.
    if r.no_sigma * 2 > r.offered {
        eprintln!(
            "WARNING: {} of {} detections had no {:?} sigma. Check the column exists and the \
             source is the right one before believing any count below.",
            r.no_sigma, r.offered, cfg.data.sigma_source
        );
    }

    let visits = set.visits(cfg.data.visit_tolerance_seconds);
    let nights = set.nights(&visits, cfg.data.night_boundary_days);
    let cap = set.sigma_cap();
    let anchor_capable = nights.values().filter(|v| v.len() >= 2).count();
    println!(
        "  {} visits over {} nights ({} with >=2 visits, so able to form an anchor); \
         sigma cap {:.4}\" = sqrt(2) * max (strict)",
        visits.len(),
        nights.len(),
        anchor_capable,
        cap / ARCSEC,
    );

    let night_keys: Vec<i64> = nights.keys().copied().collect();
    // one index per visit, for the extension: topocentric, because a prediction is compared with
    // detections as observed
    let mut visit_index: Vec<VisitIndex> = Vec::new();
    for (&nk, vis) in &nights {
        for &vi in vis {
            if let Some(v) = VisitIndex::build(&set.obs, &visits[vi], nk) {
                visit_index.push(v);
            }
        }
    }
    let params = ExtendParams {
        tolerance: cfg.extension.tolerance_arcsec * ARCSEC,
        dispersion: cfg.extension.chance_dispersion,
        exclude_anchor_nights: cfg.extension.exclude_anchor_nights,
        min_support: cfg.extension.min_support,
        max_chance_probability: cfg.extension.max_chance_probability,
    };

    let mut writer = csv::Writer::from_path(&cfg.io.output)?;
    writer.write_record([
        "r_au", "rdot", "night_a", "night_b", "chi2", "n_opp", "n_support", "lambda", "p_chance",
        "accepted", "epoch", "h", "x", "y", "z", "vx", "vy", "vz",
    ])?;

    // ---- the (r, rdot) node set -------------------------------------------------------------
    let specs: Vec<(NodeSpec, f64)> = if cfg.data.explicit_nodes.is_empty() {
        let cad = run_cadence(&set, &visits, &night_keys, &nights)?;
        let stat = match cfg.search.lever_stat.as_str() {
            "median" => Stat::Median,
            "p90" => Stat::P90,
            other => return Err(format!("search.lever_stat must be median|p90; got {other}").into()),
        };
        let eps_rad = cfg.search.eps_arcsec * ARCSEC;
        let nodes = build_ladder(
            |r| geometry_for_shell(&cad, r, mu),
            cfg.search.r_min_au,
            cfg.search.r_max_au,
            eps_rad,
            mu,
            stat,
            cfg.search.max_nodes,
        )
        .map_err(|e| format!("ladder: {e:?}"))?;
        let s = summarize(&nodes);
        println!(
            "  LADDER sized from this run's own geometry: {} shells, {} nodes over {:.1}-{:.1} AU\
             \n    anchor baseline {:.1} d, extension window {} epochs, eps {:.1}\"\
             \n    rdot nodes per shell: max {}, single on {} shells{}",
            s.shells,
            s.nodes,
            s.r_min.unwrap_or(f64::NAN),
            s.r_max.unwrap_or(f64::NAN),
            (cad.t1 - cad.t0).abs(),
            cad.ts.len(),
            cfg.search.eps_arcsec,
            s.max_rdot_nodes,
            s.single_rdot_shells,
            match s.collapse_au {
                Some(r) => format!(", rdot axis collapses beyond {r:.1} AU"),
                None => String::new(),
            }
        );
        if let Some(n0) = nodes.first() {
            println!(
                "    levers ({}): W = {:.4} AU, Z = {:.4} rad.d  =>  half-cell dgamma {:.3e} \
                 (dr/r {:.2}% at {:.1} AU), drdot {:.3e}",
                cfg.search.lever_stat, n0.w, n0.z, n0.dgamma,
                100.0 * n0.dgamma * n0.r, n0.r, n0.drdot,
            );
        }
        if ladder_only {
            return Ok(());
        }
        // 🔴 Carry the cell's LOWER r edge alongside the node. `build_ladder` emits r = 1/gamma
        // at the cell's lower GAMMA edge, i.e. its LARGEST r -- the value that makes the pairing
        // gate smallest, where the cell needs the gate that bounds its fastest (smallest-r)
        // member. See build_anchors.
        nodes
            .iter()
            .map(|n| (NodeSpec { r_au: n.r, rdot_au_per_day: n.rdot },
                      1.0 / (1.0 / n.r + 2.0 * n.dgamma)))
            .collect()
    } else {
        println!(
            "  🔴 EXPLICIT NODES ({} of them): this is a DIAGNOSTIC, not a search. The ladder is \
             bypassed, so recall here means 'recovered when pointed at', not 'found'.",
            cfg.data.explicit_nodes.len()
        );
        // An explicit node has no cell, so assert and bound coincide -- which is exactly why
        // this defect stayed invisible for as long as every run used hand-placed nodes.
        cfg.data.explicit_nodes.iter().map(|n| (*n, n.r_au)).collect()
    };

    for (spec, r_gate_lower) in &specs {
        let node = Node { r: spec.r_au, rdot: spec.rdot_au_per_day };
        // anchors per night: every visit pair inside the configured intra-night span
        let mut per_night: Vec<(i64, Vec<Anchor>)> = Vec::new();
        for &nk in &night_keys {
            let vs = &nights[&nk];
            let mut anchors = Vec::new();
            for a in 0..vs.len() {
                for b in (a + 1)..vs.len() {
                    let (va, vb) = (&visits[vs[a]], &visits[vs[b]]);
                    let dt_h = (set.obs[vb[0]].epoch - set.obs[va[0]].epoch).abs() * 24.0;
                    if dt_h > cfg.anchor.max_intra_night_hours {
                        continue;
                    }
                    anchors.extend(build_anchors(
                        &set.obs, va, vb, spec.r_au, *r_gate_lower, cfg.gate.k_sigma, mu, cap,
                    ));
                }
            }
            if !anchors.is_empty() {
                per_night.push((nk, anchors));
            }
        }
        let n_anchors: usize = per_night.iter().map(|(_, a)| a.len()).sum();

        // --sign-census accumulators; all zero and unread on a normal run.
        let mut census_hist: std::collections::BTreeMap<usize, usize> = Default::default();
        let mut census_pairs = 0usize;
        let mut census_fragmented = 0usize;
        let mut census_multi: Vec<(f64, usize, f64, f64)> = Vec::new();
        let mut census_ratio: Vec<f64> = Vec::new();
        let mut census_no_guess = 0usize;
        let mut guided_rel: Vec<f64> = Vec::new();
        let mut guided_outcome_diff = 0usize;
        let mut guided_disagree_pairs: Vec<pointdexter::spherical_pair::Pair> = Vec::new();
        let (mut guided_ns_scan, mut guided_ns_guided) = (0u64, 0u64);
        let (mut guided_ns_scan_rooted, mut guided_ns_guided_rooted) = (0u64, 0u64);
        let mut guided_rooted = 0usize;

        let mut candidates = 0usize;
        let mut gated = 0usize;
        let mut accepted = 0usize;
        let mut supported = 0usize;
        let mut tally = pointdexter::spherical_pair_anchor::RejectTally::default();
        for i in 0..per_night.len() {
            for j in (i + 1)..per_night.len() {
                let (ka, aa) = (&per_night[i].0, &per_night[i].1);
                let (kb, ab) = (&per_night[j].0, &per_night[j].1);
                let dt_days = (set.obs[ab[0].first].epoch - set.obs[aa[0].first].epoch).abs();
                if dt_days > cfg.anchor.max_anchor_baseline_days {
                    continue;
                }
                if sign_census {
                    // 🔴 Same gate, same pairs, no solve. See anchor_pairs_and_census.
                    let (census_rows, disagreeing) =
                        pointdexter::spherical_pair_anchor::anchor_pairs_and_census(
                            &set.obs, aa, ab, &node, cfg.gate.k_sigma, mu, cap,
                        );
                    guided_disagree_pairs.extend(disagreeing);
                    for (dt_days, c, g) in census_rows {
                        census_pairs += 1;
                        *census_hist.entry(c.sign_changes).or_insert(0usize) += 1;
                        guided_ns_scan += g.ns_scan;
                        guided_ns_guided += g.ns_guided;
                        if !g.outcome_agrees {
                            guided_outcome_diff += 1;
                        }
                        if g.rel_diff.is_finite() {
                            guided_rel.push(g.rel_diff);
                        }
                        if g.scan_bound {
                            guided_ns_scan_rooted += g.ns_scan;
                            guided_ns_guided_rooted += g.ns_guided;
                            guided_rooted += 1;
                        }
                        if c.segments > 1 {
                            census_fragmented += 1;
                        }
                        if c.sign_changes == 1 {
                            // 🔴 Against the CONVERGED root, never the bracket centre -- the
                            // bracket is ~1.19x wide, so the centre measures the grid and
                            // reports +-9% for a perfect guess. Ratio, not difference: h spans
                            // decades across the scan.
                            if c.h_guess_frac.is_finite() && c.root_frac > 0.0 {
                                census_ratio.push(c.h_guess_frac / c.root_frac);
                            } else {
                                census_no_guess += 1;
                            }
                        }
                        if c.sign_changes >= 2 {
                            // Keep the spread: how far apart the roots a guess must choose
                            // between actually are, in the units a guess would be expressed in.
                            census_multi.push((
                                dt_days,
                                c.sign_changes,
                                c.first_bracket_frac,
                                c.last_bracket_frac,
                            ));
                        }
                    }
                    continue;
                }
                let (cands, t) = anchor_pairs_and_solve(
                    &set.obs, aa, ab, &node, cfg.gate.k_sigma, mu, cap, cfg.anchor.chi2_max,
                );
                tally.unbound_node += t.unbound_node;
                tally.no_bound_root += t.no_bound_root;
                tally.no_state += t.no_state;
                tally.no_prediction += t.no_prediction;
                tally.chi2 += t.chi2;
                candidates += cands.len();
                // 🔴 report the GATED count too. Without it "0 candidates, 0 rejected" cannot
                // distinguish "the gate admitted no pairs" from "every pair failed chi2" -- which
                // is the exact confusion the module docs warn about, and this binary had it.
                gated += cands.len() + t.total();
                for c in cands {
                    let sup = extend_candidate(&set.obs, &c, &node, mu, &visit_index, &params);
                    if sup.accepted {
                        accepted += 1;
                    }
                    supported += usize::from(sup.n_support > 0);
                    writer.write_record([
                        format!("{}", spec.r_au),
                        format!("{}", spec.rdot_au_per_day),
                        format!("{ka}"),
                        format!("{kb}"),
                        format!("{:.6}", c.chi2),
                        format!("{}", sup.n_opportunities),
                        format!("{}", sup.n_support),
                        format!("{:.6}", sup.lambda),
                        format!("{:.6e}", sup.p_chance),
                        format!("{}", sup.accepted),
                        format!("{:.9}", c.epoch),
                        format!("{:.12e}", c.h),
                        format!("{:.12e}", c.state[0]),
                        format!("{:.12e}", c.state[1]),
                        format!("{:.12e}", c.state[2]),
                        format!("{:.12e}", c.state[3]),
                        format!("{:.12e}", c.state[4]),
                        format!("{:.12e}", c.state[5]),
                    ])?;
                }
            }
        }
        if sign_census {
            println!(
                "  node r={:.2} rdot={:+.2e}: {n_anchors} anchors over {} nights -> \
                 {census_pairs} gated pairs scanned",
                spec.r_au,
                spec.rdot_au_per_day,
                per_night.len()
            );
            if census_pairs == 0 {
                println!("    (no gated pairs at this node -- nothing to census)");
                continue;
            }
            // 🔴 Report the ZERO-bracket bin too. Without it "single-rooted" cannot be
            // distinguished from "no root at all", and no-bound-root is known to dominate at
            // high |rdot| -- those pairs are not evidence either way.
            for (k, n) in &census_hist {
                println!(
                    "    {k} bracket(s): {n:>10} pairs  ({:.4}%)",
                    100.0 * (*n as f64) / (census_pairs as f64)
                );
            }
            let multi: usize = census_hist.iter().filter(|(k, _)| **k >= 2).map(|(_, n)| *n).sum();
            let rooted: usize = census_hist.iter().filter(|(k, _)| **k >= 1).map(|(_, n)| *n).sum();
            println!(
                "    MULTI-ROOT {multi} of {rooted} pairs that have any root ({:.4}%); \
                 fragmented domain {census_fragmented} ({:.4}% of scanned)",
                if rooted > 0 { 100.0 * (multi as f64) / (rooted as f64) } else { 0.0 },
                100.0 * (census_fragmented as f64) / (census_pairs as f64)
            );
            if !census_ratio.is_empty() {
                census_ratio.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let q = |f: f64| census_ratio[((census_ratio.len() - 1) as f64 * f) as usize];
                println!(
                    "    GUESS h_guess/h_root over {} single-rooted pairs ({} had no guess):",
                    census_ratio.len(),
                    census_no_guess
                );
                println!(
                    "      p01 {:.4}  p05 {:.4}  p50 {:.4}  p95 {:.4}  p99 {:.4}  \
                     min {:.4}  max {:.4}",
                    q(0.01), q(0.05), q(0.50), q(0.95), q(0.99),
                    census_ratio[0],
                    census_ratio[census_ratio.len() - 1]
                );
                // 🔴 The operational number: a bracket [h_guess/k, h_guess*k] must CONTAIN the
                // root, or the solver has to fall back to the full scan. Report the miss rate,
                // because a guess that is usually excellent and occasionally absent is a
                // different engineering problem from one that is uniformly mediocre.
                print!("      inside a factor k of the root:");
                for k in [1.2_f64, 1.5, 2.0, 3.0, 5.0, 10.0] {
                    let n = census_ratio.iter().filter(|r| **r >= 1.0 / k && **r <= k).count();
                    print!("  k={k}: {:.4}%", 100.0 * (n as f64) / (census_ratio.len() as f64));
                }
                println!();
            }
            if !guided_rel.is_empty() || guided_outcome_diff > 0 {
                guided_rel.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let worst = guided_rel.last().copied().unwrap_or(f64::NAN);
                let over_tol = guided_rel.iter().filter(|d| **d > 1e-10).count();
                println!(
                    "    GUIDED vs SCAN: outcome disagreements {guided_outcome_diff}; \
                     |dh|/h over {} both-bound pairs: median {:.3e} max {:.3e}; \
                     over 1e-10: {over_tol}",
                    guided_rel.len(),
                    guided_rel[guided_rel.len() / 2],
                    worst
                );
                // 🔴 Split rooted from unrooted. The guided path pays expansion AND the fallback
                // scan when there is no root, so a single blended speedup hides a real loss on
                // the population that dominates near the bound limit.
                let sp = |a: u64, b: u64| if b > 0 { (a as f64) / (b as f64) } else { f64::NAN };
                println!(
                    "      speedup all pairs {:.2}x ({:.1} vs {:.1} ms); \
                     rooted-only {:.2}x over {guided_rooted} pairs",
                    sp(guided_ns_scan, guided_ns_guided),
                    guided_ns_scan as f64 / 1e6,
                    guided_ns_guided as f64 / 1e6,
                    sp(guided_ns_scan_rooted, guided_ns_guided_rooted),
                );
                report_guided_disagreements(&guided_disagree_pairs, &node, mu);
            }
            // 🔴 Report the retry path even when it never fires -- "0" is the statement that the
            // dependency's stagnating region was not entered, which is exactly what a run after
            // the fix needs to be able to say out loud.
            {
                use pointdexter::spherical_pair::{ANOMALY_RETRY_COUNT, ANOMALY_RETRY_FIRST};
                let n = ANOMALY_RETRY_COUNT.load(std::sync::atomic::Ordering::Relaxed);
                print!("    ANOMALY_TOL_RETRY fired {n} time(s) cumulatively");
                match ANOMALY_RETRY_FIRST.get() {
                    Some([r, rdot, h, dt]) => println!(
                        "; first at r={r:.17e} rdot={rdot:.17e} h={h:.17e} dt={dt:.17e}"
                    ),
                    None => println!(),
                }
            }
            if !census_multi.is_empty() {
                let mut sp: Vec<f64> =
                    census_multi.iter().map(|(_, _, f, l)| (l / f).abs()).collect();
                sp.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let dt: Vec<f64> = census_multi.iter().map(|(d, _, _, _)| *d).collect();
                println!(
                    "    multi-root root spread h_last/h_first: median {:.3}x, max {:.3}x; \
                     baseline dt {:.3}-{:.3} d",
                    sp[sp.len() / 2],
                    sp[sp.len() - 1],
                    dt.iter().cloned().fold(f64::INFINITY, f64::min),
                    dt.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                );
            }
            continue;
        }
        println!(
            "  node r={:.2} rdot={:+.2e}: {n_anchors} anchors over {} nights -> {gated} pairs \
             gated -> {candidates} candidates; rejected {} (unbound {} / no-bound-root {} / \
             no-state {} / no-prediction {} / chi2 {})",
            spec.r_au,
            spec.rdot_au_per_day,
            per_night.len(),
            tally.total(),
            tally.unbound_node,
            tally.no_bound_root,
            tally.no_state,
            tally.no_prediction,
            tally.chi2
        );
        println!(
            "    extension: {supported} of {candidates} candidates found any cross-night support, \
             {accepted} accepted at p <= {:.0e} with >= {} support (dispersion {:.2})",
            cfg.extension.max_chance_probability,
            cfg.extension.min_support,
            cfg.extension.chance_dispersion
        );
    }
    writer.flush()?;
    println!("candidates written to {}", cfg.io.output.display());
    Ok(())
}

fn main() {
    let cli = Cli::parse();
    let text = match std::fs::read_to_string(&cli.config) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("cannot read config {}: {e}", cli.config.display());
            std::process::exit(2);
        }
    };
    let cfg: Config = match serde_yaml::from_str(&text) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("config {} is not valid: {e}", cli.config.display());
            std::process::exit(2);
        }
    };
    if let Err(e) = cfg.validate() {
        eprintln!("config {} is invalid: {e}", cli.config.display());
        std::process::exit(2);
    }
    println!("config {} loaded and validated", cli.config.display());

    if cli.print_config {
        match serde_yaml::to_string(&cfg) {
            Ok(s) => println!("\n{s}"),
            Err(e) => eprintln!("could not re-emit config: {e}"),
        }
    }
    if cli.self_check {
        self_check(&cfg);
    }
    if cli.run || cli.ladder_only || cli.sign_census {
        if let Err(e) = run_search(&cfg, cli.ladder_only, cli.sign_census) {
            eprintln!("run failed: {e}");
            std::process::exit(1);
        }
    }
    if !cli.print_config && !cli.self_check && !cli.run && !cli.ladder_only && !cli.sign_census {
        println!(
            "nothing to do: stage 0 validates configuration and can --self-check.\n\
             the search itself arrives with stage 3 (detections) and stage 4 (extension)."
        );
    }
}

/// 🔴 DIAGNOSTIC. Why did `solve_h_guided` reach a different OUTCOME from `solve_h`?
///
/// `solve_h_guided` falls back to `solve_h` whenever its expansion fails to observe a sign
/// change, so its `Bound` set is a SUPERSET of the scan's and every disagreement is the guided
/// path finding a root the 80-point scan did not. Three mechanisms could do that, and they have
/// very different consequences:
///
///   (a) TWO roots inside ONE coarse cell. An even crossing count is invisible to a scan that
///       only compares adjacent grid points -- and it would mean `RESULT_c2_signcensus.md`'s
///       single-rootedness is a property of the GRID's resolution, not of `F`.
///   (b) A root against a domain gap, where a non-finite `F` sets `prev = None` and the scan
///       cannot see across the reset.
///   (c) `F` grazing zero at the level of its own rounding, in which case the extra "root" is
///       numerical noise and the guided path is manufacturing candidates.
///
/// (c) is separated from (a)/(b) by scale: a real crossing has |F| at the bracket ends
/// comparable to `F`'s range over the scan, a grazing one is orders below it.
fn report_guided_disagreements(pairs: &[Pair], node: &Node, mu: f64) {
    if pairs.is_empty() {
        return;
    }
    // Dense enough that a coarse cell (ratio 1.1911) holds ~250 fine points, so a root PAIR
    // inside one cell is resolved rather than stepped over a second time.
    const N_FINE: usize = 20001;

    println!("    🔴 GUIDED-vs-SCAN DISAGREEMENTS: {} pair(s) at this node", pairs.len());
    for (idx, pair) in pairs.iter().enumerate() {
        let Some(hm) = h_max(node, mu) else { continue };
        let f = |h: f64| pair.f_of_h(h, node, mu);
        let dt = pair.b.epoch - pair.a.epoch;

        // The coarse grid the solver actually walks -- from the solver's own iterator, never a
        // transcription of it.
        let coarse: Vec<(f64, Option<f64>)> = h_scan_grid(hm).map(|h| (h, f(h))).collect();
        let coarse_finite = coarse.iter().filter(|(_, v)| v.is_some()).count();
        let mut coarse_segments = 0usize;
        let mut prev_finite = false;
        for (_, v) in &coarse {
            if v.is_some() && !prev_finite {
                coarse_segments += 1;
            }
            prev_finite = v.is_some();
        }

        // A much finer walk of the SAME range.
        let (ln_lo, ln_hi) = (H_SCAN_LO_FRAC.ln(), H_SCAN_HI_FRAC.ln());
        let fine: Vec<(f64, Option<f64>)> = (0..N_FINE)
            .map(|i| {
                let h = hm * (ln_lo + (ln_hi - ln_lo) * (i as f64) / ((N_FINE - 1) as f64)).exp();
                (h, f(h))
            })
            .collect();
        let f_scale = fine
            .iter()
            .filter_map(|(_, v)| *v)
            .fold(0.0_f64, |m, v| m.max(v.abs()));
        let fine_finite = fine.iter().filter(|(_, v)| v.is_some()).count();

        // Sign changes on the fine grid, mirroring the solver's gap handling (a non-finite
        // point resets the run, so a crossing spanning a gap is not counted).
        let mut roots: Vec<(f64, f64, f64, f64)> = Vec::new(); // (h_lo, h_hi, f_lo, f_hi)
        let mut prev: Option<(f64, f64)> = None;
        for (h, v) in &fine {
            let Some(v) = *v else {
                prev = None;
                continue;
            };
            if let Some((hp, vp)) = prev {
                if vp * v < 0.0 {
                    roots.push((hp, *h, vp, v));
                }
            }
            prev = Some((*h, v));
        }

        let guided = match pair.solve_h_guided(node, mu) {
            Solution::Bound { h, .. } => h / hm,
            _ => f64::NAN,
        };
        println!(
            "      [{idx}] dt={dt:.4} d  h_max={hm:.6e}  |F|max={f_scale:.3e}  \
             coarse finite {coarse_finite}/{} segs {coarse_segments}  fine finite {fine_finite}/{N_FINE}  \
             fine crossings {}  guided h/h_max={guided:.6e}",
            coarse.len(),
            roots.len()
        );

        for (hlo, hhi, flo, fhi) in &roots {
            // Bisect in place: `brent` is private to the library and this only needs to locate
            // the root well enough to say WHICH coarse cell it falls in.
            let (mut a, mut b, mut fa) = (*hlo, *hhi, *flo);
            for _ in 0..80 {
                let m = 0.5 * (a + b);
                match f(m) {
                    Some(fm) if fa * fm < 0.0 => b = m,
                    Some(fm) => {
                        a = m;
                        fa = fm;
                    }
                    None => break,
                }
            }
            let root = 0.5 * (a + b);
            // Which coarse cell contains it, and what the scan saw at that cell's ends.
            let j = coarse.iter().position(|(h, _)| *h > root);
            let cell = match j {
                Some(0) | None => "outside the coarse cells".to_string(),
                Some(j) => {
                    let (h0, v0) = coarse[j - 1];
                    let (h1, v1) = coarse[j];
                    let n_here = roots
                        .iter()
                        .filter(|(rl, rh, _, _)| {
                            let m = 0.5 * (rl + rh);
                            m > h0 && m <= h1
                        })
                        .count();
                    format!(
                        "coarse cell {}..{} [{:.4e},{:.4e}] F=({}, {}) -- {n_here} fine crossing(s) in it",
                        j - 1,
                        j,
                        h0 / hm,
                        h1 / hm,
                        v0.map(|x| format!("{x:+.3e}")).unwrap_or("NONE".into()),
                        v1.map(|x| format!("{x:+.3e}")).unwrap_or("NONE".into()),
                    )
                }
            };
            println!(
                "          root h/h_max={:.6e}  ends F=({:+.3e},{:+.3e})  |F|end/|F|max={:.2e}  {cell}",
                root / hm,
                flo,
                fhi,
                flo.abs().max(fhi.abs()) / f_scale.max(f64::MIN_POSITIVE),
            );

            // 🔴 THE POINT OF THIS DIAGNOSTIC. `solve_h` saw this sign change too -- it walks the
            // same grid -- so the only way it reported NoBoundRoot is that `brent` handed back
            // `None`, and `brent` has exactly two of those: a non-bracketing input (excluded, the
            // ends differ in sign) and `let fs = f(s)?` on an interior trial point.
            //
            // Run the REAL `brent` on the real coarse cell, with `f` wrapped in a recorder, so the
            // trial sequence and the first non-finite point are the solver's own, not a replica's.
            if let Some(j) = j {
                if j > 0 {
                    let (h0, v0) = coarse[j - 1];
                    let (h1, v1) = coarse[j];
                    if let (Some(fa), Some(fb)) = (v0, v1) {
                        let trace = std::cell::RefCell::new(Vec::<(f64, bool)>::new());
                        let traced = |h: f64| {
                            let v = pair.f_of_h(h, node, mu);
                            trace.borrow_mut().push((h, v.is_some()));
                            v
                        };
                        let out = pointdexter::spherical_pair::brent(&traced, h0, fa, h1, fb);
                        let t = trace.borrow();
                        let n_none = t.iter().filter(|(_, ok)| !ok).count();
                        let first_none = t.iter().find(|(_, ok)| !ok).map(|(h, _)| *h / hm);
                        println!(
                            "          BRENT on that cell -> {}  ({} evals, {n_none} non-finite{})",
                            match out {
                                Some(r) => format!("Some({:.6e})", r / hm),
                                None => "🔴 None  => solve_h reports NoBoundRoot".to_string(),
                            },
                            t.len(),
                            match first_none {
                                Some(h) => format!(", first at h/h_max={h:.9e}"),
                                None => String::new(),
                            },
                        );
                        if let Some(hbad) = first_none {
                            // How wide is the non-finite region, and does the scan's own grid
                            // straddle it? A sliver narrower than a fine-grid step is invisible
                            // to every diagnostic that has been run on this so far.
                            // Bisect toward the cell ends, which the scan already evaluated as
                            // finite, so each edge costs 60 evals rather than an unbounded walk.
                            let probe = |x: f64| pair.f_of_h(x * hm, node, mu).is_some();
                            let (mut lo_ok, mut lo) = (h0 / hm, hbad);
                            for _ in 0..60 {
                                let m = 0.5 * (lo_ok + lo);
                                if probe(m) { lo_ok = m } else { lo = m }
                            }
                            let (mut hi_ok, mut hi) = (h1 / hm, hbad);
                            for _ in 0..60 {
                                let m = 0.5 * (hi_ok + hi);
                                if probe(m) { hi_ok = m } else { hi = m }
                            }
                            println!(
                                "          non-finite region spans h/h_max [{lo:.12e}, {hi:.12e}] \
                                 (width {:.3e}); H_SCAN_HI_FRAC={:.12e}",
                                hi - lo,
                                H_SCAN_HI_FRAC,
                            );

                            // WHICH component goes non-finite. `f_of_h` has four ways to return
                            // None; naming the one that fires is the difference between "the
                            // Kepler step failed" and "the geometry was unreachable".
                            let hmid = 0.5 * (lo + hi) * hm;
                            let dta = pair.a.epoch - pair.t_ref;
                            let dtb = pair.b.epoch - pair.t_ref;
                            let ra = pointdexter::spherical_pair::radial_at(node, hmid, dta, mu);
                            let rb = pointdexter::spherical_pair::radial_at(node, hmid, dtb, mu);
                            let v_sq = node.rdot * node.rdot + (hmid / node.r) * (hmid / node.r);
                            let alpha = 2.0 / node.r - v_sq / mu; // 1/a; -> 0 at the parabolic ceiling
                            let geom = match (ra, rb) {
                                (Some((r0, _)), Some((r1, _))) => {
                                    let p0 = range_quadratic(&pair.a.observer, &pair.a.rho_hat, r0);
                                    let p1 = range_quadratic(&pair.b.observer, &pair.b.rho_hat, r1);
                                    format!(
                                        "range_quadratic a={} b={}",
                                        if p0.is_some() { "ok" } else { "🔴 None" },
                                        if p1.is_some() { "ok" } else { "🔴 None" },
                                    )
                                }
                                _ => "not reached".to_string(),
                            };
                            println!(
                                "          at the sliver: canonical_step a={} b={}; {geom}; \
                                 alpha=1/a={alpha:.6e} AU^-1, a={:.4e} AU",
                                if ra.is_some() { "ok" } else { "🔴 None" },
                                if rb.is_some() { "ok" } else { "🔴 None" },
                                1.0 / alpha,
                            );

                            // 🔴 OURS or the DEPENDENCY'S? `canonical_step` can return None from
                            // `solve_for_universal_anomaly(..).ok()?` (spacerocks) or from its own
                            // finite checks after `polish_anomaly` (this repo). Call the
                            // dependency directly on the failing epoch to split the two.
                            let dt_bad = if ra.is_none() { dta } else { dtb };
                            match solve_for_universal_anomaly(
                                node.r,
                                node.rdot,
                                alpha,
                                mu,
                                dt_bad,
                                pointdexter::spherical_pair::ANOMALY_TOL,
                                pointdexter::spherical_pair::ANOMALY_MAX_ITER,
                            ) {
                                Err(e) => println!(
                                    "          => spacerocks solve_for_universal_anomaly ERRORS \
                                     on dt={dt_bad:.6} d: {e:?}  (the dependency, not polish)"
                                ),
                                Ok(s) => println!(
                                    "          => spacerocks returned s={s:.9e} (finite={}); the \
                                     None comes from OUR polish/finite checks, not the dependency",
                                    s.is_finite()
                                ),
                            }
                            // Slow convergence or genuine stagnation? More iterations fixes the
                            // first; only a different residual fixes the second. This decides
                            // whether the repair is a constant or a solver.
                            let variants = [
                                ("tol 1e-12, iter 1e4", 1e-12, 10_000usize),
                                ("tol 1e-12, iter 1e6", 1e-12, 1_000_000usize),
                                ("tol 1e-10, iter 100", 1e-10, 100usize),
                                ("tol 1e-8,  iter 100", 1e-8, 100usize),
                                ("tol 1e-6,  iter 100", 1e-6, 100usize),
                            ];
                            let mut report = Vec::new();
                            for (label, tol, iter) in variants {
                                let r = solve_for_universal_anomaly(
                                    node.r, node.rdot, alpha, mu, dt_bad, tol, iter,
                                );
                                report.push(format!(
                                    "{label}: {}",
                                    match r {
                                        Ok(s) => format!("ok s={s:.6e}"),
                                        Err(_) => "ERR".to_string(),
                                    }
                                ));
                            }
                            println!("          => {}", report.join(" | "));
                        }
                    }
                }
            }
        }
    }
}
