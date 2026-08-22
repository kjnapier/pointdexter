//! A tracklet manager for C2: content-addressed dedup plus reverse-index subset pruning.
//!
//! ## Why this exists
//!
//! `NOTE_c2_ps1_escalation_lineage.md` §4 measures the funnel on node 1 of the full-arc run:
//! **26.7M anchor pairs -> 783,883 candidates with >= 2 supporting tracklets -> 35,537 distinct
//! orbits after dedup.** That last step is a **22x** reduction and the largest single one in the
//! chain. It was not computable until `6b232bf` began emitting the anchor and support detection
//! ids, because a candidate's identity IS its detection set.
//!
//! ## Provenance
//!
//! The design is Michael Lackner's, from his Harvard master's thesis work with M. Holman, as
//! implemented in `pytrax`'s `Discover.hpp`; the reverse-index ancestor is `Chain.hpp` in
//! `ann_1.1.2`. Reviewed against source in `~/ssolink-search/RESULT_c2_tracklet_manager_review.md`;
//! the structures are documented in `import/legacy-search/TRACKLET_MANAGEMENT.md`. **Written
//! fresh** rather than ported: `GAP_ANALYSIS_pytrax_to_pointdexter.md` claims the commented
//! `deh::GroupStore` is "~half-done", but it has only 1 of the 4 structures and its
//! `find_existing_group` is a linear scan.
//!
//! 🔴 Nothing from the Bernstein `ORBFIT` C is used here, so the Numerical-Recipes prohibition in
//! `CLAUDE.md` is not engaged.
//!
//! ## The two operations, and why they are cheap
//!
//! 1. **"Seen this exact set?"** -- content-addressed, O(1). pytrax calls this `tracklets_set_`;
//!    PS1 called it `combo_checked`. 🔴 A linear scan here (what `deh::find_existing_group` does)
//!    is ~2.8e10 set comparisons over this funnel and is not viable.
//! 2. **"Proper subset of something kept?"** -- walk **ONE** member detection's adjacency list.
//!    ⭐ That is exact, and it is the trick both PS1 and the Rust sketch miss: if a candidate is a
//!    subset of some tracklet, then EVERY member lies in that tracklet, so that tracklet appears
//!    in EVERY member's list. Scanning any one list cannot miss it -- no union (PS1 builds one,
//!    costing ~k x too many tests), no intersection. We take the SMALLEST list, which pytrax's
//!    `.back()` does not.
//!
//! ⭐ The two reinforce: because (1) guarantees no duplicate detection sets exist, "proper subset"
//! reduces to a size comparison plus `includes` -- no separate equality test. pytrax's own comment
//! says as much.

use std::collections::HashMap;

/// A tracklet's stable identifier. Monotone, and **never reused** -- a removed id's slot is
/// retired and the counter marches on, so anything holding an id (an orbit fit, a cluster
/// membership) keeps referring to the same thing for the life of a run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct TrackletId(pub u32);

/// A tracklet: a **sorted, duplicate-free** set of detection ids.
///
/// 🔴 The sort is load-bearing, not cosmetic. Subset testing is a merge walk and dedup is by
/// content, and both are wrong on unsorted input. pytrax states this as a precondition
/// (`Discover.hpp` l.1061: *"@pre tracklet detections are sorted in exposure order"*) and enforces
/// it by convention; here the constructor is the only way to build one, so it cannot be violated.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Tracklet(Box<[u32]>);

impl Tracklet {
    /// Sorts and de-duplicates. The only constructor.
    pub fn new(mut ids: Vec<u32>) -> Option<Self> {
        if ids.is_empty() {
            return None;
        }
        ids.sort_unstable();
        ids.dedup();
        Some(Tracklet(ids.into_boxed_slice()))
    }

    pub fn ids(&self) -> &[u32] {
        &self.0
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Is `self` a subset of `other`? Linear merge walk; both sides are sorted by construction.
    fn is_subset_of(&self, other: &Tracklet) -> bool {
        if self.0.len() > other.0.len() {
            return false;
        }
        let mut it = other.0.iter();
        self.0.iter().all(|x| it.any(|y| y == x))
    }
}

/// What an insert did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Insert {
    /// Stored, with a fresh id.
    Added(TrackletId),
    /// An identical detection set already exists.
    Duplicate(TrackletId),
    /// A proper superset already exists, so this carries no information the store lacks.
    SubsetOf(TrackletId),
}

/// The four coupled indices of `TRACKLET_MANAGEMENT.md` §"Data model", kept consistent by
/// `insert` and `remove`.
///
/// ⚠️ Removal retires the id permanently; `next_id` never goes backwards.
#[derive(Debug, Default)]
pub struct TrackletStore {
    by_id: HashMap<TrackletId, Tracklet>,
    /// Content-addressed: the sorted id set -> its tracklet. This is the O(1) "seen it?".
    content: HashMap<Box<[u32]>, TrackletId>,
    /// Reverse index: detection id -> the tracklets containing it. Keeps subset tests local.
    adjacency: HashMap<u32, Vec<TrackletId>>,
    /// Level index: size -> tracklets of that size.
    ///
    /// ⭐ Maintained here at insert, deliberately. PS1's `remove_subsets_v2` recomputed a
    /// length ordering with `std::sort` at the call site **and then iterated a lexicographically
    /// ordered `std::set` instead**, so its longest-first optimisation never actually ran. An
    /// index that is part of the store's invariants cannot be dropped that way.
    by_size: HashMap<usize, Vec<TrackletId>>,
    next_id: u32,
}

impl TrackletStore {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.by_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_id.is_empty()
    }

    pub fn get(&self, id: TrackletId) -> Option<&Tracklet> {
        self.by_id.get(&id)
    }

    /// O(1). The operation that runs most.
    pub fn find_exact(&self, t: &Tracklet) -> Option<TrackletId> {
        self.content.get(&t.0).copied()
    }

    /// Is `t` a PROPER subset of a stored tracklet? Returns the first such superset.
    ///
    /// ⭐ Walks exactly ONE member's adjacency list -- the **smallest** -- which is exact: a
    /// superset of `t` contains every member of `t`, so it appears in every member's list.
    ///
    /// 🔴 "Proper" needs no equality test here: `insert` guarantees no two stored tracklets share
    /// a detection set, so a strictly larger size is sufficient.
    pub fn find_superset(&self, t: &Tracklet) -> Option<TrackletId> {
        let probe = t
            .ids()
            .iter()
            .filter_map(|d| self.adjacency.get(d).map(|v| (v.len(), d)))
            .min()?
            .1;
        // A member with no adjacency entry means nothing stored contains it, so no superset can
        // exist -- `min()` returning None above already covers that case.
        for &cand in self.adjacency.get(probe)? {
            if let Some(other) = self.by_id.get(&cand) {
                if other.len() > t.len() && t.is_subset_of(other) {
                    return Some(cand);
                }
            }
        }
        None
    }

    /// Dedup, then subset-prune, then store. Cheapest test first, as `check_add_tracklet` does.
    pub fn insert(&mut self, t: Tracklet) -> Insert {
        if let Some(id) = self.find_exact(&t) {
            return Insert::Duplicate(id);
        }
        if let Some(id) = self.find_superset(&t) {
            return Insert::SubsetOf(id);
        }
        let id = TrackletId(self.next_id);
        self.next_id += 1;
        for &d in t.ids() {
            self.adjacency.entry(d).or_default().push(id);
        }
        self.by_size.entry(t.len()).or_default().push(id);
        self.content.insert(t.0.clone(), id);
        self.by_id.insert(id, t);
        Insert::Added(id)
    }

    /// Remove a tracklet, keeping all four indices consistent. The id is retired, never reused.
    pub fn remove(&mut self, id: TrackletId) -> Option<Tracklet> {
        let t = self.by_id.remove(&id)?;
        self.content.remove(&t.0);
        if let Some(v) = self.by_size.get_mut(&t.len()) {
            v.retain(|&x| x != id);
            if v.is_empty() {
                self.by_size.remove(&t.len());
            }
        }
        for d in t.ids() {
            if let Some(v) = self.adjacency.get_mut(d) {
                v.retain(|&x| x != id);
                if v.is_empty() {
                    self.adjacency.remove(d);
                }
            }
        }
        Some(t)
    }

    /// Ids of stored tracklets with exactly `n` detections. The level index.
    pub fn ids_by_size(&self, n: usize) -> &[TrackletId] {
        self.by_size.get(&n).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Sizes present, largest first.
    ///
    /// ⭐ Subset pruning wants LONGEST-first, so each subset is retired once by its maximal
    /// superset. ⚠️ Note this is the opposite direction from pytrax's *linking* sweep, which goes
    /// shortest->longest so a track is built by one merge chain and never rediscovered. The two
    /// size orderings are duals; do not copy one where the other belongs.
    pub fn sizes_descending(&self) -> Vec<usize> {
        let mut s: Vec<usize> = self.by_size.keys().copied().collect();
        s.sort_unstable_by(|a, b| b.cmp(a));
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tk(v: &[u32]) -> Tracklet {
        Tracklet::new(v.to_vec()).expect("non-empty")
    }

    /// Brute force: scan EVERY stored tracklet. The reference the shortcut must match.
    fn superset_by_brute_force(store: &TrackletStore, t: &Tracklet) -> bool {
        store
            .by_id
            .values()
            .any(|o| o.len() > t.len() && t.is_subset_of(o))
    }

    #[test]
    fn the_constructor_sorts_and_dedups_so_the_invariant_cannot_be_violated() {
        // 🔴 pytrax states sortedness as a precondition and enforces it by convention
        // (`Discover.hpp` l.1061). Here it is the only constructor, so unsorted input cannot exist.
        assert_eq!(tk(&[9, 2, 5, 2]).ids(), &[2, 5, 9]);
        assert!(Tracklet::new(vec![]).is_none(), "an empty tracklet is not a thing");
    }

    #[test]
    fn exact_dedup_is_by_content_not_by_insertion_order() {
        // ...by CONTENT: the same ids in any order are the same tracklet.
        let mut s = TrackletStore::new();
        let a = match s.insert(tk(&[3, 1, 2])) {
            Insert::Added(id) => id,
            other => panic!("first insert should be Added, got {other:?}"),
        };
        // same set, different order given to the constructor
        assert_eq!(s.insert(tk(&[2, 3, 1])), Insert::Duplicate(a));
        assert_eq!(s.len(), 1, "a duplicate must not create a second entry");
    }

    #[test]
    fn a_proper_subset_is_rejected_but_an_equal_set_is_a_duplicate() {
        let mut s = TrackletStore::new();
        let big = match s.insert(tk(&[1, 2, 3, 4])) {
            Insert::Added(id) => id,
            o => panic!("{o:?}"),
        };
        assert_eq!(s.insert(tk(&[2, 3])), Insert::SubsetOf(big));
        // 🔴 The distinction matters: equality is NOT a proper subset. Because dedup guarantees no
        // two stored sets are equal, `find_superset` can use a bare size test -- this asserts that
        // shortcut does not misclassify the equal case.
        assert_eq!(s.insert(tk(&[1, 2, 3, 4])), Insert::Duplicate(big));
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn a_superset_of_a_stored_tracklet_is_still_added() {
        // The store prunes subsets, not supersets: a longer arc carries new information.
        let mut s = TrackletStore::new();
        s.insert(tk(&[1, 2]));
        assert!(matches!(s.insert(tk(&[1, 2, 3])), Insert::Added(_)));
        assert_eq!(s.len(), 2);
    }

    #[test]
    fn probing_one_adjacency_list_agrees_with_scanning_every_tracklet() {
        // ⭐⭐ THE load-bearing test. `find_superset` walks ONE member's adjacency -- the smallest --
        // on the argument that a superset must appear in EVERY member's list. If that reasoning is
        // wrong the store silently keeps subsets, which is invisible downstream: the funnel just
        // fails to reduce, and nothing looks broken. So it is checked against brute force rather
        // than argued.
        let mut s = TrackletStore::new();
        let mut seed: u64 = 0x9E3779B97F4A7C15;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        // Populate with overlapping sets drawn from a small id pool, so adjacency lists are dense
        // and supersets genuinely occur.
        for _ in 0..400 {
            let n = 2 + (rnd() % 6) as usize;
            let ids: Vec<u32> = (0..n).map(|_| (rnd() % 25) as u32).collect();
            if let Some(t) = Tracklet::new(ids) {
                s.insert(t);
            }
        }
        assert!(s.len() > 20, "store too small to be exercising anything: {}", s.len());

        let mut agreed = 0;
        let mut found = 0;
        for _ in 0..3000 {
            let n = 1 + (rnd() % 5) as usize;
            let ids: Vec<u32> = (0..n).map(|_| (rnd() % 25) as u32).collect();
            let Some(t) = Tracklet::new(ids) else { continue };
            let fast = s.find_superset(&t).is_some();
            let slow = superset_by_brute_force(&s, &t);
            assert_eq!(fast, slow, "one-list probe disagreed with a full scan on {:?}", t.ids());
            agreed += 1;
            found += usize::from(slow);
        }
        assert!(agreed > 2000, "only {agreed} comparisons ran");
        // A test where the answer is always "no" would pass vacuously.
        assert!(found > 50, "only {found} of {agreed} probes found a superset; not exercising the path");
    }

    #[test]
    fn removal_leaves_all_four_indices_consistent() {
        let mut s = TrackletStore::new();
        let a = match s.insert(tk(&[1, 2, 3])) {
            Insert::Added(id) => id,
            o => panic!("{o:?}"),
        };
        s.insert(tk(&[3, 4, 5, 6]));
        let removed = s.remove(a).expect("was present");
        assert_eq!(removed.ids(), &[1, 2, 3]);

        assert!(s.get(a).is_none(), "by_id");
        assert!(s.find_exact(&tk(&[1, 2, 3])).is_none(), "content index still answers");
        assert!(!s.ids_by_size(3).contains(&a), "by_size still lists it");
        for d in [1u32, 2, 3] {
            assert!(
                !s.adjacency.get(&d).map(|v| v.contains(&a)).unwrap_or(false),
                "adjacency for detection {d} still lists the removed tracklet"
            );
        }
        // detection 3 is still in the surviving tracklet, so its entry must persist
        assert!(s.adjacency.contains_key(&3), "shared detection lost its adjacency entry");
        // and the removal must not have stranded a subset that is now insertable again
        assert!(matches!(s.insert(tk(&[1, 2, 3])), Insert::Added(_)));
    }

    #[test]
    fn ids_are_never_reused() {
        // 🔴 Stable ids are what let solutions, clusters and adjacency reference a tracklet for the
        // life of a run. Reuse would silently re-point them.
        let mut s = TrackletStore::new();
        let a = match s.insert(tk(&[1, 2])) {
            Insert::Added(id) => id,
            o => panic!("{o:?}"),
        };
        s.remove(a);
        let b = match s.insert(tk(&[7, 8])) {
            Insert::Added(id) => id,
            o => panic!("{o:?}"),
        };
        assert_ne!(a, b, "a retired id came back");
        assert!(b.0 > a.0);
    }

    #[test]
    fn the_level_index_is_maintained_and_ordered_longest_first() {
        // ⭐ PS1's remove_subsets_v2 sorted longest-first and then iterated a lexicographic set,
        // so the ordering never took effect. Here it is read off a maintained index.
        let mut s = TrackletStore::new();
        s.insert(tk(&[1, 2]));
        s.insert(tk(&[3, 4, 5]));
        s.insert(tk(&[6, 7, 8, 9]));
        assert_eq!(s.sizes_descending(), vec![4, 3, 2]);
        assert_eq!(s.ids_by_size(3).len(), 1);
        assert_eq!(s.ids_by_size(99).len(), 0, "absent size must be empty, not a panic");
    }
}

