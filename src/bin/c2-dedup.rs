//! Collapse a C2 candidate file with the tracklet manager, and report the reduction.
//!
//! `NOTE_c2_ps1_escalation_lineage.md` §4 puts the fit tier behind a dedup step measured at
//! 26.7M pairs -> 783,883 accepted -> **35,537 distinct orbits**, a 22x reduction. That figure was
//! obtained by ORBIT-SPACE BINNING (`_c2_dupcollapse_orbitspace.awk`), an explicit PROXY, because
//! the candidate CSV carried no detection ids. `6b232bf` now emits them, so this does the real
//! thing: dedup by DETECTION SET, which is what PS1's `combo_checked` and pytrax's
//! `tracklets_set_` do.
//!
//! ⚠️ A candidate's detection set here is the four anchor ids plus `support_ids` WHEN PRESENT.
//! With `extension.emit_support_ids = false` (the default) the set is the four anchor ids alone,
//! so every set has size 4 and **subset pruning cannot fire** -- a 4-set is never a proper subset
//! of another 4-set. What that configuration measures is the collapse of the SAME anchor pair
//! recurring across ladder nodes. Turn support ids on to exercise the subset half.
//!
//! usage: c2-dedup <candidates.csv> [--all]
//!        --all   include rejected candidates (default: `accepted == true` only, matching the awk)

use std::collections::HashMap;
use std::process::ExitCode;

use pointdexter::tracklet_store::{Insert, Tracklet, TrackletStore};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        eprintln!("usage: c2-dedup <candidates.csv> [--all]");
        return ExitCode::from(2);
    };
    let accepted_only = !args.iter().any(|a| a == "--all");

    let mut rdr = match csv::ReaderBuilder::new().has_headers(true).from_path(path) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("cannot read {path}: {e}");
            return ExitCode::from(1);
        }
    };

    let mut store = TrackletStore::new();
    let (mut rows, mut considered, mut added, mut dup, mut subset, mut malformed) =
        (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    // How many ladder nodes each surviving detection set was produced at -- the multiplicity the
    // dedup is removing. Kept only for the size histogram, so it is bounded by the distinct count.
    let mut node_multiplicity: HashMap<u32, u32> = HashMap::new();
    let mut sizes: HashMap<usize, u64> = HashMap::new();

    let mut rec = csv::StringRecord::new();
    while matches!(rdr.read_record(&mut rec), Ok(true)) {
        rows += 1;
        if rec.len() < 26 {
            malformed += 1;
            continue;
        }
        if accepted_only && rec.get(9) != Some("true") {
            continue;
        }
        let mut ids: Vec<u32> = Vec::with_capacity(8);
        let mut bad = false;
        for i in 22..26 {
            match rec.get(i).and_then(|s| s.trim().parse::<u32>().ok()) {
                Some(v) => ids.push(v),
                None => bad = true,
            }
        }
        if bad {
            malformed += 1;
            continue;
        }
        if let Some(s) = rec.get(26) {
            for tok in s.split_whitespace() {
                if let Ok(v) = tok.parse::<u32>() {
                    ids.push(v);
                }
            }
        }
        let Some(t) = Tracklet::new(ids) else {
            malformed += 1;
            continue;
        };
        considered += 1;
        let n = t.len();
        match store.insert(t) {
            Insert::Added(id) => {
                added += 1;
                *sizes.entry(n).or_insert(0) += 1;
                node_multiplicity.insert(id.0, 1);
            }
            Insert::Duplicate(id) => {
                dup += 1;
                *node_multiplicity.entry(id.0).or_insert(0) += 1;
            }
            Insert::SubsetOf(_) => subset += 1,
        }
    }

    let distinct = store.len() as u64;
    println!("file                {path}");
    println!("rows                {rows}");
    println!(
        "considered          {considered}   ({})",
        if accepted_only { "accepted == true" } else { "all candidates" }
    );
    println!("distinct sets       {distinct}");
    println!("  added             {added}");
    println!("  exact duplicates  {dup}");
    println!("  proper subsets    {subset}");
    if malformed > 0 {
        println!("malformed/skipped   {malformed}");
    }
    if distinct > 0 {
        println!("REDUCTION           {:.1}x", considered as f64 / distinct as f64);
    }

    let mut sz: Vec<_> = sizes.iter().collect();
    sz.sort();
    println!("distinct-set sizes  {:?}", sz);
    if sz.len() == 1 {
        println!(
            "  🔴 every set has the same size, so SUBSET pruning could not fire. Re-run the \\
search with extension.emit_support_ids = true to exercise it."
        );
    }

    let mut mult: Vec<u32> = node_multiplicity.values().copied().collect();
    mult.sort_unstable();
    if !mult.is_empty() {
        let max = *mult.last().unwrap_or(&0);
        let med = mult[mult.len() / 2];
        println!("multiplicity        median {med}  max {max}   (rows per surviving set)");
    }
    ExitCode::SUCCESS
}
