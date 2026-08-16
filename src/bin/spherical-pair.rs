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

use pointdexter::io::load_detections::load_detections;
use pointdexter::spherical_pair::{Node, gate_radius};
use pointdexter::spherical_pair_anchor::{Anchor, anchor_pairs_and_solve, build_anchors};
use pointdexter::spherical_pair_extend::{ExtendParams, VisitIndex, extend_candidate};
use pointdexter::spherical_pair_index::{
    BaryIndex, angle_of_chord, chord_of_angle, gate_radius_astrometric, pair_within_gate,
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
    /// Nodes to search. 🔴 A stopgap: the sized `(r, rdot)` ladder is a closed form evaluated from
    /// the run's OWN cadence and elongation coverage (`spherical_pair_grid::build_ladder`), and
    /// wiring that driver is the next stage. An explicit list is honest about not having it;
    /// emitting a uniform grid here would look like the sized one and would not be.
    explicit_nodes: Vec<NodeSpec>,
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
fn run_search(cfg: &Config) -> Result<(), Box<dyn std::error::Error>> {
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

    for spec in &cfg.data.explicit_nodes {
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
                        &set.obs, va, vb, spec.r_au, cfg.gate.k_sigma, mu, cap,
                    ));
                }
            }
            if !anchors.is_empty() {
                per_night.push((nk, anchors));
            }
        }
        let n_anchors: usize = per_night.iter().map(|(_, a)| a.len()).sum();

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
    if cli.run {
        if let Err(e) = run_search(&cfg) {
            eprintln!("run failed: {e}");
            std::process::exit(1);
        }
    }
    if !cli.print_config && !cli.self_check && !cli.run {
        println!(
            "nothing to do: stage 0 validates configuration and can --self-check.\n\
             the search itself arrives with stage 3 (detections) and stage 4 (extension)."
        );
    }
}
