//! Cross-language gate on `spherical_pair` -- Stage 1 of `PLAN_c2_pointdexter_rust.md`.
//!
//! The oracle is `fixtures/fh_fixture.json`, schema `ssolink.c2.fh_fixture/3`, generated in
//! ssolink by `_fh_fixture.py` from the Python prototype. Its observer positions come from
//! astropy's builtin ephemeris, so freezing them is what lets this test run with **no ephemeris,
//! no network and no Python** -- and makes the comparison exact rather than approximate.
//!
//! 🔴 The rule this file exists to satisfy: **a fixture only one side checks is not an oracle.**
//! ssolink's `tests/test_fh_fixture.py` asserts a set of properties on the same JSON; every one
//! of them that Stage 1 covers is asserted here too, and where Rust can do better than Python it
//! does -- above all `mu`, which is compared to `Origin::SSB.mu()` bit-for-bit rather than to a
//! transcribed literal.
//!
//! Coverage: all six blocks. `cases`, `solve_h_cases`, `state_cases` and `gate_cases` gate
//! Stage 1 (`spherical_pair.rs`); `grid_cases` and `ladder_case` gate Stage 2
//! (`spherical_pair_grid.rs`, the port of ssolink's `kernels/spherical-pair/grid.py`).
//!
//! ⭐ The Stage 2 half is a *stricter* gate than the Stage 1 half, and the fixture says so: the
//! levers and the ladder are pure arithmetic on frozen arrays -- no solver, no ephemeris -- so
//! `lever_rel = 1e-12` there is a real requirement rather than a budget for a different
//! universal-variable implementation.

use nalgebra::Vector3;
use pointdexter::spherical_pair::{gate_radius, h_max, radial_at, range_quadratic, Node, Pair, PairPoint, Solution};
use pointdexter::spherical_pair_grid::{
    binding_branch, build_ladder, gamma_cell, n_rdot_nodes, parallax_lever, rdot_cell, rdot_span,
    summarize, zoom_lever, Geometry, GridNode, Stat,
};
use serde_json::Value;
use spacerocks::coordinates::Origin;

// --------------------------------------------------------------------------------- fixture load

fn doc() -> Value {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/fixtures/fh_fixture.json");
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    serde_json::from_str(&text).expect("fixture is not valid JSON")
}

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

fn vec3(v: &Value) -> Vector3<f64> {
    let a = v.as_array().expect("not an array");
    assert_eq!(a.len(), 3, "expected a 3-vector, got {v}");
    Vector3::new(f(&a[0]), f(&a[1]), f(&a[2]))
}

fn point(v: &Value) -> PairPoint {
    PairPoint {
        epoch: f(&v["epoch"]),
        rho_hat: vec3(&v["rho_hat"]),
        observer: vec3(&v["observer"]),
    }
}

/// A `Pair` from any fixture block carrying `t_ref` and a two-element `points` array.
fn pair_of(case: &Value, t_ref: f64) -> Pair {
    let pts = case["points"].as_array().expect("no points");
    assert_eq!(pts.len(), 2);
    Pair { t_ref, a: point(&pts[0]), b: point(&pts[1]) }
}

fn tol(doc: &Value, key: &str) -> f64 {
    f(&doc["tolerances"][key])
}

fn rel_err(got: f64, want: f64) -> f64 {
    if want == 0.0 { got.abs() } else { (got - want).abs() / want.abs() }
}

// ------------------------------------------------------------------------------ schema and frame

#[test]
fn schema_and_units_are_the_ones_this_port_was_written_against() {
    let d = doc();
    assert_eq!(d["schema"], "ssolink.c2.fh_fixture/3");
    assert_eq!(d["frame"], "barycentric ecliptic J2000");
    assert_eq!(d["units"]["length"], "AU");
    assert_eq!(d["units"]["time"], "day");
    assert_eq!(d["units"]["angle"], "radian");
    assert_eq!(d["units"]["epoch"], "JD TDB");
    assert_eq!(d["units"]["h"], "AU^2/day");
}

#[test]
fn the_fixture_mu_is_bit_for_bit_what_spacerocks_resolves_for_the_ssb() {
    // 🔴 The check only the Rust side can make, and the one that closes a live hazard: a port
    // written next to matt-0.rs/tangent_v2.rs reaches for its neighbours' mu, and if the fixture
    // were heliocentric the gate and ladder would miss by ~0.067% -- eight orders above the
    // tolerances below, reading as a porting bug rather than a frame mismatch.
    //
    // Exact equality, not approx: spacerocks carries a SECOND barycentric constant
    // (constants::MU_BARY, 1.7e-11 away, commented out on the same line as this one). Matching
    // what the code CALLS is the point.
    let d = doc();
    assert_eq!(f(&d["mu"]), Origin::SSB.mu());
    assert_ne!(f(&d["mu"]), Origin::SUN.mu());
}

#[test]
fn the_fixture_covers_the_regime_c2_is_for() {
    // C2's rationale is long arcs; a fixture of short baselines only would gate nothing.
    let d = doc();
    let cases = d["cases"].as_array().unwrap();
    let max_baseline = cases.iter().map(|c| f(&c["baseline_days"])).fold(0.0, f64::max);
    assert!(max_baseline >= 730.0, "no multi-opposition baseline: {max_baseline}");
    let names: std::collections::BTreeSet<&str> =
        cases.iter().map(|c| c["name"].as_str().unwrap()).collect();
    assert!(names.len() >= 3, "needs several eccentricity regimes, got {names:?}");
}

// -------------------------------------------------------------------------------------- /1 cases

/// Everything `f_of_h` and its two intermediates must reproduce, case by case.
#[test]
fn f_of_h_reproduces_the_oracle_at_three_depths() {
    let d = doc();
    let mu = f(&d["mu"]);
    let f_tol = tol(&d, "f_of_h_abs");
    let radial_tol = tol(&d, "radial_rel");
    let pos_tol = tol(&d, "position_abs");

    let mut worst_f = 0.0f64;
    let mut worst_radial = 0.0f64;
    let mut worst_pos = 0.0f64;

    for case in d["cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let bl = f(&case["baseline_days"]);
        let label = format!("{name} {bl}d");
        let t_ref = f(&case["t_ref"]);
        let node = Node { r: f(&case["node"]["r"]), rdot: f(&case["node"]["rdot"]) };
        let h_true = f(&case["h_true"]);
        let pair = pair_of(case, t_ref);

        // depth 1 -- the range quadratic: |position| == r_assumed, and the same position.
        for rq in case["range_quadratic_at_h_true"].as_array().unwrap() {
            let r_assumed = f(&rq["r_assumed"]);
            let want = vec3(&rq["position"]);
            let pt = if f(&rq["epoch"]) == pair.a.epoch { &pair.a } else { &pair.b };
            let got = range_quadratic(&pt.observer, &pt.rho_hat, r_assumed)
                .unwrap_or_else(|| panic!("{label}: range quadratic returned None"));
            assert!(rel_err(got.norm(), r_assumed) < radial_tol, "{label}: |p| != r_assumed");
            worst_pos = worst_pos.max((got - want).norm());
        }

        // depth 2 -- the radial solve, including the dt = 0 identity that catches a sign or
        // units error in the universal-variable step.
        for s in case["radial_at_h_true"].as_array().unwrap() {
            let dt = f(&s["dt"]);
            let (r, rdot) = radial_at(&node, h_true, dt, mu)
                .unwrap_or_else(|| panic!("{label}: radial_at({dt}) returned None"));
            worst_radial = worst_radial.max(rel_err(r, f(&s["r"])));
            // rdot passes through zero, so it gets an absolute floor as well.
            let e = (rdot - f(&s["rdot"])).abs();
            assert!(e < 1e-15 || e / f(&s["rdot"]).abs() < radial_tol, "{label} dt={dt}: rdot");
            if dt == 0.0 {
                assert_eq!(r, node.r, "{label}: radial solve moved the node at dt=0");
                assert_eq!(rdot, node.rdot, "{label}: radial solve moved rdot at dt=0");
            }
        }

        // depth 3 -- the whole residual, at h_true and across the scan.
        let f_at_true = pair.f_of_h(h_true, &node, mu).expect("F(h_true) is None");
        assert!(f_at_true.abs() < f_tol, "{label}: F(h_true) = {f_at_true}");
        worst_f = worst_f.max((f_at_true - f(&case["f_at_h_true"])).abs());

        let scan = &case["f_of_h_scan"];
        let ss: Vec<f64> = scan["h_over_h_true"].as_array().unwrap().iter().map(f).collect();
        let want: Vec<f64> = scan["f"].as_array().unwrap().iter().map(f).collect();
        assert_eq!(ss.len(), want.len());
        assert!((ss[0] - 0.80).abs() < 1e-12 && (ss[ss.len() - 1] - 1.20).abs() < 1e-12,
                "the scan bracket is part of the contract");

        let got: Vec<f64> = ss
            .iter()
            .map(|s| pair.f_of_h(h_true * s, &node, mu).unwrap_or_else(|| panic!("{label}: None at h/h_true={s}")))
            .collect();
        for (i, (g, w)) in got.iter().zip(want.iter()).enumerate() {
            let e = (g - w).abs();
            assert!(e < f_tol, "{label}: scan point {i} (h/h_true={}) off by {e:.3e}", ss[i]);
            worst_f = worst_f.max(e);
        }

        // The properties the scan is FOR -- asserted on the Rust values, not the stored ones,
        // so a port that reproduced the numbers but broke the structure still fails.
        let diffs: Vec<f64> = got.windows(2).map(|w| w[1] - w[0]).collect();
        assert!(diffs.iter().all(|&x| x > 0.0) || diffs.iter().all(|&x| x < 0.0),
                "{label}: F is not monotonic on [0.80, 1.20]*h_true");
        // Monotonic plus opposite-signed ends IS "exactly one root", and unlike counting
        // sign changes it survives an exact zero on the scan -- which this port hits, because
        // F(h_true) comes out at 0.0 rather than the prototype's 1e-17 in some cases.
        assert!(got[0] * got[got.len() - 1] < 0.0, "{label}: no sign change across the scan");
        let i = (0..got.len())
            .find(|&i| got[i] * got[0] <= 0.0)
            .unwrap_or_else(|| panic!("{label}: monotonic but never crosses"));
        assert!(i > 0 && ss[i - 1] <= 1.0 && 1.0 <= ss[i],
                "{label}: the root is in [{}, {}], which does not bracket h_true", ss[i - 1], ss[i]);
    }

    println!(
        "worst |dF| = {worst_f:.3e} rad (tol {f_tol:.0e}); worst radial rel = {worst_radial:.3e} \
         (tol {radial_tol:.0e}); worst |dposition| = {worst_pos:.3e} AU (tol {pos_tol:.0e})"
    );
    assert!(worst_pos < pos_tol, "worst position error {worst_pos:.3e} AU");
    assert!(worst_radial < radial_tol);
}

#[test]
fn conditioning_reproduces_and_is_linear_in_baseline() {
    // The property that makes C2 want LONG arcs. Worth pinning separately: if it broke, the
    // method's whole rationale would be gone and every other test here would still pass.
    //
    // 🔴 The tolerance here is DERIVED, not chosen. The stored value is a central difference
    // `(F+ - F-)/(2*frac)` at `frac = 1e-4`, so whatever absolute disagreement the two `F`
    // implementations have is amplified by `1/(2*frac)` = 5000x. With the fixture's own
    // `f_of_h_abs` budget that is `f_tol/frac` absolute on the sensitivity, and no port can do
    // better without the two sides sharing a propagator. The measured worst is printed below --
    // it lands ~1e5 inside this bound, which is the real statement about agreement.
    let d = doc();
    let mu = f(&d["mu"]);
    let frac = 1e-4;
    let s_tol = tol(&d, "f_of_h_abs") / frac;
    let mut worst_abs = 0.0f64;
    let mut worst_rel = 0.0f64;

    let mut rows: std::collections::BTreeMap<String, Vec<(f64, f64)>> = Default::default();
    for case in d["cases"].as_array().unwrap() {
        let node = Node { r: f(&case["node"]["r"]), rdot: f(&case["node"]["rdot"]) };
        let h = f(&case["h_true"]);
        let pair = pair_of(case, f(&case["t_ref"]));
        let fp = pair.f_of_h(h * (1.0 + frac), &node, mu).unwrap();
        let fm = pair.f_of_h(h * (1.0 - frac), &node, mu).unwrap();
        let s = ((fp - fm) / (2.0 * frac)).abs();

        let want = f(&case["sensitivity_h_dFdh"]);
        assert!(want > 0.0);
        assert!((s - want).abs() < s_tol, "{}: |h dF/dh| {s} != {want}", case["name"]);
        worst_abs = worst_abs.max((s - want).abs());
        worst_rel = worst_rel.max(rel_err(s, want));
        rows.entry(case["name"].as_str().unwrap().to_string())
            .or_default()
            .push((f(&case["baseline_days"]), s));
    }

    println!("worst |h dF/dh| discrepancy: {worst_abs:.3e} abs (bound {s_tol:.0e}), \
              {worst_rel:.3e} rel");

    for (name, mut rs) in rows {
        rs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for w in rs.windows(2) {
            let ((b0, s0), (b1, s1)) = (w[0], w[1]);
            assert!(s1 > s0, "{name}: sensitivity fell from {b0}d to {b1}d");
            // Linear => the sensitivity ratio tracks the baseline ratio. Measured deviation is
            // <= 0.10% over every step here, so 1% bites; 5% would pass on almost anything.
            assert!(rel_err(s1 / s0, b1 / b0) < 0.01,
                    "{name}: {b0}d->{b1}d ratio {:.4} != {:.4}", s1 / s0, b1 / b0);
        }
    }
}

// ------------------------------------------------------------------------------- /2  solve_h

#[test]
fn solve_h_brackets_on_physics_and_classifies_every_outcome() {
    let d = doc();
    let mu = f(&d["mu"]);
    let block = &d["solve_h_cases"];
    let t_ref = f(&block["t_ref"]);
    let h_true = f(&block["h_true"]);
    let pair = pair_of(block, t_ref);
    let root_tol = tol(&d, "solve_h_rel");
    let f_tol = tol(&d, "f_of_h_abs");

    let mut roots: Vec<f64> = Vec::new();
    let mut by_r: Vec<(f64, f64)> = Vec::new();
    let mut statuses: Vec<String> = Vec::new();

    for o in block["offsets"].as_array().unwrap() {
        let dr = f(&o["dr_frac"]);
        let node = Node { r: f(&o["r_assumed"]), rdot: f(&o["rdot_assumed"]) };
        let got = pair.solve_h(&node, mu);
        let want = o["status"].as_str().unwrap();
        statuses.push(want.to_string());

        match (&got, want) {
            (Solution::UnboundNode, "unbound_node") => {
                assert!(o["h_max"].is_null(), "dr={dr}: h_max on an unbound node");
                assert!(h_max(&node, mu).is_none());
            }
            (Solution::NoBoundRoot { h_max: hm, f_lo, f_hi }, "no_bound_root") => {
                assert!(rel_err(*hm, f(&o["h_max"])) < 1e-12, "dr={dr}: h_max");
                // 🔴 The classification must be justified by the bracket ENDS, not by a bare
                // "None": this is what proves both sides rejected for the same reason.
                println!("  dr={dr:+.2} f_lo {f_lo:.12e} vs {:.12e} (d={:.2e}) | \
                          f_hi {f_hi:.12e} vs {:.12e} (d={:.2e})",
                         f(&o["f_at_h_lo"]), (f_lo - f(&o["f_at_h_lo"])).abs(),
                         f(&o["f_at_h_hi"]), (f_hi - f(&o["f_at_h_hi"])).abs());
                assert!((f_lo - f(&o["f_at_h_lo"])).abs() < f_tol, "dr={dr}: F at the low end");
                assert!((f_hi - f(&o["f_at_h_hi"])).abs() < f_tol, "dr={dr}: F at the high end");
                assert!(f_lo * f_hi > 0.0, "dr={dr}: rejection claimed despite a sign change");
            }
            (Solution::Bound { h, h_max: hm }, "bound_root") => {
                let want_root = f(&o["h_root"]);
                assert!(rel_err(*h, want_root) < root_tol, "dr={dr}: root {h} != {want_root}");
                assert!(*h > 0.0 && *h <= *hm, "dr={dr}: root outside (0, h_max]");
                let fr = pair.f_of_h(*h, &node, mu).unwrap();
                assert!(fr.abs() < f_tol, "dr={dr}: F(root) = {fr}");
                assert!(f(&o["f_at_h_lo"]) * f(&o["f_at_h_hi"]) < 0.0,
                        "dr={dr}: root claimed with no sign change across the bracket");
                roots.push(*h / h_true);
                by_r.push((node.r, *h));
                if dr == 0.0 {
                    // The node IS the truth there, so the solve must return the truth.
                    assert!(rel_err(*h, h_true) < root_tol, "the zero offset lost h_true");
                }
            }
            (g, w) => panic!("dr={dr}: got {g:?}, fixture says {w}"),
        }
    }

    // 🔴 The regression this block exists for: a port bracketing h in [0.2, 5]*h_hint passes
    // every /1 case and returns None here. If the fixture stopped reaching that deep, the bug
    // would come back silently -- so assert the coverage, not just the agreement.
    assert!(!roots.is_empty());
    let shallowest = roots.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(shallowest < 0.2, "shallowest root {shallowest:.3} is inside a hinted window");
    assert!(statuses.iter().any(|s| s != "bound_root"), "no rejection case exercised");

    // Monotone in the asserted node: a bigger r needs more angular momentum. A sign or branch
    // error breaks the ordering while leaving individual roots plausible.
    by_r.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let hs: Vec<f64> = by_r.iter().map(|x| x.1).collect();
    assert!(hs.windows(2).all(|w| w[0] <= w[1]), "h_root not monotone in r_assumed: {by_r:?}");
}

// --------------------------------------------------------------------------- /2  state rebuild

#[test]
fn state_from_solution_reproduces_the_state_and_its_invariants() {
    let d = doc();
    let mu = f(&d["mu"]);
    let state_tol = tol(&d, "state_rel");
    assert!(!d["state_cases"].as_array().unwrap().is_empty(), "no state cases");

    for case in d["state_cases"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let node = Node { r: f(&case["node"]["r"]), rdot: f(&case["node"]["rdot"]) };
        let h = f(&case["h"]);
        let pair = pair_of(case, f(&case["t_ref"]));

        // The state is reported at the FIRST observation, not at t_ref.
        assert!((f(&case["epoch"]) - pair.a.epoch).abs() < 1e-9, "{name}: wrong epoch");

        let st = pair.state_from_solution(h, &node, mu).unwrap_or_else(|| panic!("{name}: None"));
        let want: Vec<f64> = case["state"].as_array().unwrap().iter().map(f).collect();
        for (i, (g, w)) in st.iter().zip(want.iter()).enumerate() {
            assert!(g.is_finite());
            // Components pass through zero, so relative alone would be meaningless.
            let e = (g - w).abs();
            assert!(e < 1e-15 || e / w.abs() < state_tol, "{name}: state[{i}] {g} != {w}");
        }

        // The invariants, recomputed here rather than read: comparing six floats alone would let
        // a port that got the orbit-plane SIGN wrong pass on three of them.
        let pos = Vector3::new(st[0], st[1], st[2]);
        let vel = Vector3::new(st[3], st[4], st[5]);
        let inv = &case["invariants"];
        let r0 = pos.norm();
        let rdot0 = pos.dot(&vel) / r0;
        let h0 = pos.cross(&vel).norm();

        assert!(rel_err(r0, f(&inv["r"])) < state_tol, "{name}: |r|");
        assert!(rel_err(h0, h) < state_tol, "{name}: |r x v| = {h0}, h = {h}");
        // ...and against the INDEPENDENT radial solve, which computes r(t0) a different way.
        let (r_exp, rdot_exp) = radial_at(&node, h, pair.a.epoch - pair.t_ref, mu).unwrap();
        assert!(rel_err(r0, r_exp) < state_tol, "{name}: |r| != radial_at");
        assert!(rel_err(r_exp, f(&inv["r_from_radial_solve"])) < state_tol, "{name}: radial r");
        let e = (rdot0 - rdot_exp).abs();
        assert!(e < 1e-15 || e / rdot_exp.abs() < state_tol, "{name}: rdot != radial_at");
        assert!(rel_err(rdot_exp, f(&inv["rdot_from_radial_solve"])) < state_tol, "{name}: radial rdot");

        // The one invariant that ties the state back to the observation it was built from.
        let geo = pos - pair.a.observer;
        let los_resid = (geo / geo.norm() - pair.a.rho_hat).norm();
        assert!(los_resid < 1e-12, "{name}: state is off the line of sight by {los_resid:.3e}");
        assert!(f(&inv["line_of_sight_residual"]) < 1e-12, "{name}: fixture LOS residual");
    }
}

// ------------------------------------------------------------------------------ /2  pairing gate

#[test]
fn the_gate_closed_form_reproduces_and_bounds_the_barycentric_set() {
    let d = doc();
    let mu = f(&d["mu"]);
    let cases = d["gate_cases"]["cases"].as_array().unwrap();
    assert!(!cases.is_empty());

    for g in cases {
        let r = f(&g["r"]);
        let dt = f(&g["dt_days"]);
        let got = gate_radius(r, dt, mu);
        assert!(rel_err(got, f(&g["closed_form_rad"])) < tol(&d, "gate_rel"), "gate r={r} dt={dt}");

        // A gate that reproduced a formula without bounding anything would be useless...
        assert!(got >= f(&g["sampled_barycentric_rad"]), "gate does not bound the set at r={r}");
        assert!(f(&g["barycentric_margin"]) >= 1.0);
        // ...and a hugely over-wide one would bound it while costing the search dearly.
        assert!(f(&g["barycentric_margin"]) < 1.05,
                "gate is {:.2}x the sampled set at r={r}", f(&g["barycentric_margin"]));
    }
}

#[test]
fn the_gate_does_not_bound_sky_positions_except_at_integer_years() {
    // 🔴 Pinned as a POSITIVE property, not a defect: the closed form is barycentric, production
    // sky positions are topocentric, and the observer's own parallax cancels only at an integer
    // year. A port that gates raw sky positions with the bare formula at, say, a 120-day
    // baseline silently drops true pairs.
    let d = doc();
    let cases = d["gate_cases"]["cases"].as_array().unwrap();
    let short: Vec<&Value> = cases.iter().filter(|g| f(&g["dt_days"]) < 200.0).collect();
    let year: Vec<&Value> = cases.iter().filter(|g| (f(&g["dt_days"]) - 365.25).abs() < 1e-9).collect();
    assert!(!short.is_empty() && !year.is_empty());

    for g in &short {
        assert!(f(&g["topocentric_excess"]) > 2.0,
                "parallax excess vanished at dt={}d r={}", f(&g["dt_days"]), f(&g["r"]));
    }
    for g in &year {
        assert!((f(&g["topocentric_excess"]) - 1.0).abs() < 0.05,
                "parallax did not cancel at an integer year: {g}");
    }

    // Farther object => smaller true motion => the observer's motion dominates more.
    for dt in [30.0, 120.0] {
        let mut rows: Vec<(f64, f64)> = cases
            .iter()
            .filter(|g| f(&g["dt_days"]) == dt)
            .map(|g| (f(&g["r"]), f(&g["topocentric_excess"])))
            .collect();
        rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        assert!(rows.windows(2).all(|w| w[0].1 <= w[1].1), "dt={dt}: excess not rising with r");
    }
}

// -------------------------------------------------------------------------- Stage 2's oracle

/// A [`Geometry`] from a fixture `geometry` block. The window arrays are `ts`/`Es`/`us` and the
/// anchors `t0`/`t1`/`E0`/`E1`/`u0`/`u1`, exactly as `grid.py`'s dataclass names them.
fn geometry(v: &Value) -> Geometry {
    let ts: Vec<f64> = v["ts"].as_array().expect("no ts").iter().map(f).collect();
    let es: Vec<Vector3<f64>> = v["Es"].as_array().expect("no Es").iter().map(vec3).collect();
    let us: Vec<Vector3<f64>> = v["us"].as_array().expect("no us").iter().map(vec3).collect();
    Geometry::new(
        f(&v["t0"]),
        f(&v["t1"]),
        vec3(&v["E0"]),
        vec3(&v["E1"]),
        vec3(&v["u0"]),
        vec3(&v["u1"]),
        ts,
        es,
        us,
    )
    .expect("fixture geometry is not self-consistent")
}

#[test]
fn the_grid_levers_and_the_cell_they_set_reproduce_the_oracle() {
    // The whole closed form, case by case: the two levers, the third (mean-motion) lever kept
    // for the branch comparison, and the cell each of them sets.
    let d = doc();
    let mu = f(&d["mu"]);
    let lever_rel = tol(&d, "lever_rel");
    let cases = d["grid_cases"].as_array().unwrap();
    assert!(!cases.is_empty());

    for c in cases {
        let name = c["name"].as_str().unwrap();
        let g = geometry(&c["geometry"]);
        let eps = f(&c["eps_rad"]);

        let w = parallax_lever(&g, Stat::Median);
        let z = zoom_lever(&g, Stat::Median);
        let t2 = pointdexter::spherical_pair_grid::meanmotion_lever(&g, Stat::Median);
        assert!(rel_err(w, f(&c["levers"]["W"])) < lever_rel, "W on {name}");
        assert!(rel_err(z, f(&c["levers"]["Z"])) < lever_rel, "Z on {name}");
        assert!(rel_err(t2, f(&c["levers"]["T2"])) < lever_rel, "T2 on {name}");

        let r = f(&c["cell"]["r"]);
        let dgamma = gamma_cell(w, eps).unwrap();
        let drdot = rdot_cell(r, z, eps).unwrap();
        assert!(rel_err(dgamma, f(&c["cell"]["dgamma"])) < lever_rel, "dgamma on {name}");
        assert!(rel_err(drdot, f(&c["cell"]["drdot"])) < lever_rel, "drdot on {name}");
        assert!(
            rel_err(rdot_span(r, mu), f(&c["cell"]["rdot_span"])) < lever_rel,
            "rdot_span on {name}"
        );
        assert_eq!(
            n_rdot_nodes(r, drdot, mu).unwrap(),
            c["cell"]["n_rdot_nodes"].as_u64().unwrap() as usize,
            "n_rdot_nodes on {name}"
        );
    }
}

#[test]
fn parallax_sets_the_r_cell_in_every_recorded_case() {
    // 🔴 The branch decision, gated rather than assumed. `ratio > 1` means the mean-motion
    // tolerance is the looser one, so parallax binds and the grid variable is 1/r -- NOT the
    // relative `r_max/r_min < r_tol` criterion of `grid.rs:251-278`, which is the mean-motion
    // branch and is looser here by the ratio this test reproduces.
    let d = doc();
    let mu = f(&d["mu"]);
    let lever_rel = tol(&d, "lever_rel");

    for c in d["grid_cases"].as_array().unwrap() {
        let name = c["name"].as_str().unwrap();
        let g = geometry(&c["geometry"]);
        let want = &c["binding_branch"];
        let got = binding_branch(&g, f(&want["r"]), f(&c["eps_rad"]), mu, Stat::Median).unwrap();

        assert!(rel_err(got.parallax_rel, f(&want["parallax_rel"])) < lever_rel, "{name}");
        assert!(rel_err(got.meanmotion_rel, f(&want["meanmotion_rel"])) < lever_rel, "{name}");
        assert!(rel_err(got.ratio, f(&want["ratio"])) < lever_rel, "{name}");
        assert_eq!(got.binds.as_str(), want["binds"].as_str().unwrap(), "{name}");
        assert_eq!(got.binds.as_str(), "parallax", "mean motion bound on {name}");
    }
}

#[test]
fn the_c1_least_squares_projector_would_not_reproduce_these_levers() {
    // 🔴 The porting mistake this module exists to prevent, pinned with a number instead of a
    // comment. `tangent_v2.rs:515-520` detrends the window by least squares on {1, t}; C2's
    // solve pins the direction EXACTLY at both anchors, so its projector is the affine function
    // THROUGH them. The two are not close: the least-squares residual is smaller by up to 6.6x
    // on this fixture -- and smaller W means a WIDER gamma cell, so a port that copied C1 would
    // under-sample the grid by that factor while every test that only checked "is it positive"
    // still passed.
    let d = doc();
    let mut worst = 1.0f64;
    for c in d["grid_cases"].as_array().unwrap() {
        let g = geometry(&c["geometry"]);
        let n = g.ts.len() as f64;

        // Least-squares fit of v(t) = a + b t over the window, by normal equations.
        let (mut st, mut stt) = (0.0, 0.0);
        let (mut sv, mut svt) = (Vector3::zeros(), Vector3::<f64>::zeros());
        let perp: Vec<Vector3<f64>> = g
            .es
            .iter()
            .zip(g.us.iter())
            .map(|(e, u)| e - e.dot(u) * u)
            .collect();
        for (t, v) in g.ts.iter().zip(perp.iter()) {
            st += t;
            stt += t * t;
            sv += *v;
            svt += *v * *t;
        }
        let det = n * stt - st * st;
        let b = (n * svt - st * sv) / det;
        let a = (sv - b * st) / n;
        let mut res: Vec<f64> = g
            .ts
            .iter()
            .zip(perp.iter())
            .map(|(t, v)| (v - (a + b * *t)).norm())
            .collect();
        res.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let k = res.len();
        let lsq_w = 0.5 * (res[k / 2 - 1] + res[k / 2]);

        let ratio = lsq_w / f(&c["levers"]["W"]);
        assert!(ratio < 0.99, "the two projectors agree on {} -- test is not discriminating",
                c["name"].as_str().unwrap());
        worst = worst.min(ratio);
    }
    assert!(worst < 0.3, "the C1 projector's error should reach several-fold; worst was {worst}");
}

#[test]
fn the_ladder_reproduces_the_oracle_node_by_node() {
    let d = doc();
    let mu = f(&d["mu"]);
    let lever_rel = tol(&d, "lever_rel");
    let l = &d["ladder_case"];

    // 🔴 The recorded ladder was built with ONE geometry at every shell. Asserting the flag
    // rather than trusting it keeps a future fixture that varies the geometry per shell -- which
    // is what a real run does -- from being silently checked against a constant closure here.
    assert_eq!(l["constant_geometry"], true);
    let g = geometry(&l["geometry"]);

    let got = build_ladder(
        |_r| Ok(g.clone()),
        f(&l["r_min"]),
        f(&l["r_max"]),
        f(&l["eps_rad"]),
        mu,
        Stat::Median,
        100_000,
    )
    .expect("ladder failed to build");

    let want = l["nodes"].as_array().unwrap();
    assert_eq!(got.len(), want.len(), "node count differs");
    for (i, (a, b)) in got.iter().zip(want.iter()).enumerate() {
        assert!(rel_err(a.r, f(&b["r"])) < lever_rel, "node {i} r");
        assert!(rel_err(a.rdot, f(&b["rdot"])) < lever_rel, "node {i} rdot");
        assert!(rel_err(a.dgamma, f(&b["dgamma"])) < lever_rel, "node {i} dgamma");
        assert!(rel_err(a.drdot, f(&b["drdot"])) < lever_rel, "node {i} drdot");
    }

    // Order is part of the contract: shells march outward in gamma, so `r` comes out strictly
    // decreasing and each shell's rdot nodes are contiguous.
    assert!(got.windows(2).all(|w| w[0].r >= w[1].r));
}

#[test]
fn the_ladder_summary_reproduces_and_the_collapse_is_a_suffix_property() {
    let d = doc();
    let mu = f(&d["mu"]);
    let l = &d["ladder_case"];
    let g = geometry(&l["geometry"]);
    let ladder = build_ladder(
        |_r| Ok(g.clone()),
        f(&l["r_min"]),
        f(&l["r_max"]),
        f(&l["eps_rad"]),
        mu,
        Stat::Median,
        100_000,
    )
    .unwrap();
    let s = summarize(&ladder);
    let want = &l["summary"];

    assert_eq!(s.shells, want["shells"].as_u64().unwrap() as usize);
    assert_eq!(s.nodes, want["nodes"].as_u64().unwrap() as usize);
    assert_eq!(s.max_rdot_nodes, want["max_rdot_nodes"].as_u64().unwrap() as usize);
    assert_eq!(s.single_rdot_shells, want["single_rdot_shells"].as_u64().unwrap() as usize);
    assert!(rel_err(s.r_min.unwrap(), f(&want["r_min"])) < tol(&d, "lever_rel"));
    assert!(rel_err(s.r_max.unwrap(), f(&want["r_max"])) < tol(&d, "lever_rel"));

    // ⭐ The collapse distance is the figure worth quoting -- beyond it the whole BOUND rdot span
    // fits in one cell, so the generator emits a single rdot = 0 node. Reproducing the number is
    // not enough: check it is the property it claims to be.
    let collapse = f(&want["collapse_au"]);
    assert!(rel_err(s.collapse_au.unwrap(), collapse) < tol(&d, "lever_rel"));
    let beyond: Vec<&GridNode> = ladder.iter().filter(|n| n.r >= collapse).collect();
    let shells_beyond: std::collections::BTreeSet<u64> =
        beyond.iter().map(|n| n.r.to_bits()).collect();
    assert_eq!(beyond.len(), shells_beyond.len(), "a shell beyond the collapse has >1 rdot node");
    assert!(beyond.iter().all(|n| n.rdot == 0.0), "a collapsed shell is not centred on rdot = 0");
    // ...and that it is not vacuous: the axis really does open up inside it.
    assert!(ladder.iter().any(|n| n.r < collapse && n.rdot != 0.0));
}
