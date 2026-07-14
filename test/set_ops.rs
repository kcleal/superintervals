// Set-operation tests for IntervalMap, checked against a brute-force point-set
// reference over a small integer universe.

use std::collections::BTreeSet;
use superintervals::IntervalMap;

const UNIVERSE: i32 = 200;

fn make(intervals: &[(i32, i32)]) -> IntervalMap<i32> {
    let mut m = IntervalMap::new();
    for (i, &(s, e)) in intervals.iter().enumerate() {
        m.add(s, e, i as i32);
    }
    m.build();
    m
}

// The set of integer points covered by a map's raw stored intervals.
fn points_of(m: &IntervalMap<i32>) -> BTreeSet<i32> {
    let mut pts = BTreeSet::new();
    for k in 0..m.size() {
        let it = m.at(k);
        for p in it.start..=it.end {
            pts.insert(p);
        }
    }
    pts
}

// A tiny deterministic LCG so the test needs no external rng crate.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn range(&mut self, n: i32) -> i32 {
        (self.next() % n as u64) as i32
    }
}

#[test]
fn set_ops_match_brute_force() {
    let mut rng = Lcg(0x1234_5678);
    for _ in 0..3000 {
        let na = rng.range(8);
        let nb = rng.range(8);
        let mut av = Vec::new();
        let mut bv = Vec::new();
        for _ in 0..na {
            let s = rng.range(UNIVERSE);
            av.push((s, (s + rng.range(26)).min(UNIVERSE - 1)));
        }
        for _ in 0..nb {
            let s = rng.range(UNIVERSE);
            bv.push((s, (s + rng.range(26)).min(UNIVERSE - 1)));
        }
        let a = make(&av);
        let b = make(&bv);
        let a_pts = points_of(&a);
        let b_pts = points_of(&b);

        // merge_overlaps preserves the point set and is disjoint.
        let m = a.merge_overlaps();
        assert_eq!(points_of(&m), a_pts, "merge point set");
        for k in 1..m.size() {
            assert!(m.starts[k] > m.ends[k - 1], "merged not disjoint");
        }

        // union
        let u = a.union(&b);
        let expect: BTreeSet<i32> = a_pts.union(&b_pts).copied().collect();
        assert_eq!(points_of(&u), expect, "union point set");

        // intersection
        let inter = a.intersection(&b);
        let expect: BTreeSet<i32> = a_pts.intersection(&b_pts).copied().collect();
        assert_eq!(points_of(&inter), expect, "intersection point set");

        // difference
        let d = a.difference(&b);
        let expect: BTreeSet<i32> = a_pts.difference(&b_pts).copied().collect();
        assert_eq!(points_of(&d), expect, "difference point set");

        // symmetric_difference
        let x = a.symmetric_difference(&b);
        let expect: BTreeSet<i32> = a_pts.symmetric_difference(&b_pts).copied().collect();
        assert_eq!(points_of(&x), expect, "symmetric_difference point set");

        // gaps within [0, UNIVERSE-1]
        let g = a.gaps(0, UNIVERSE - 1, -1);
        let expect: BTreeSet<i32> = (0..UNIVERSE).filter(|p| !a_pts.contains(p)).collect();
        assert_eq!(points_of(&g), expect, "gaps point set");

        // span
        match a.span() {
            Some((lo, hi)) => {
                assert_eq!(lo, *a_pts.iter().next().unwrap());
                assert_eq!(hi, *a_pts.iter().next_back().unwrap());
            }
            None => assert!(a_pts.is_empty()),
        }
    }
}

#[test]
fn merge_combiner_runs() {
    // [1,5,10] and [3,8,20] coalesce into [1,8]; combiner sums the data.
    let mut a = IntervalMap::new();
    a.add(1, 5, 10);
    a.add(3, 8, 20);
    a.build();
    let m = a.merge_overlaps_with(|x, y| x + y);
    assert_eq!(m.size(), 1);
    let it = m.at(0);
    assert_eq!((it.start, it.end, it.data), (1, 8, 30));
}

#[test]
fn expand_matches_brute_force() {
    // expand grows each interval; check the resulting point set against a manual one.
    let a = make(&[(10, 20), (100, 110)]);
    let e = a.expand(5, 5, i32::MIN, i32::MAX);
    let pts = points_of(&e);
    let expect: BTreeSet<i32> = (5..=25).chain(95..=115).collect();
    assert_eq!(pts, expect, "expand point set");

    // asymmetric
    let asym = a.expand(3, 7, i32::MIN, i32::MAX);
    let it = asym.at(0);
    assert_eq!((it.start, it.end), (7, 27));

    // clamp lo
    let c = a.expand(20, 20, 0, i32::MAX);
    assert_eq!(c.at(0).start, 0);

    // shrink past self -> dropped
    let s = make(&[(10, 12)]);
    let sh = s.expand(-5, -5, i32::MIN, i32::MAX);
    assert_eq!(sh.size(), 0);
}

#[test]
fn flank_matches_bedtools_semantics() {
    let m = make(&[(10, 20)]);
    let f = m.flank(5, 5, i32::MIN, i32::MAX); // [5,9] and [21,25]
    let pts = points_of(&f);
    let expect: BTreeSet<i32> = (5..=9).chain(21..=25).collect();
    assert_eq!(pts, expect, "flank point set");

    // left only
    let fl = m.flank(3, 0, i32::MIN, i32::MAX);
    assert_eq!(fl.size(), 1);
    let it = fl.at(0);
    assert_eq!((it.start, it.end), (7, 9));

    // clamp at lo
    let m2 = make(&[(2, 20)]);
    let fc = m2.flank(10, 0, 0, i32::MAX); // left would be [-8,1], clamp to [0,1]
    assert_eq!(fc.size(), 1);
    let it = fc.at(0);
    assert_eq!((it.start, it.end), (0, 1));
}

#[test]
fn unique_collapses_exact_duplicates() {
    let mut m = IntervalMap::new();
    m.add(10, 20, 1);
    m.add(10, 20, 2);
    m.add(10, 20, 3);
    m.add(5, 8, 9);
    m.add(5, 8, 4);
    m.build();

    // keep first
    let u = m.unique();
    assert_eq!(u.size(), 2);
    assert_eq!((u.at(0).start, u.at(0).end, u.at(0).data), (5, 8, 9));
    assert_eq!((u.at(1).start, u.at(1).end, u.at(1).data), (10, 20, 1));

    // combine: sum payloads
    let us = m.unique_with(|x, y| x + y);
    assert_eq!(us.size(), 2);
    assert_eq!(us.at(0).data, 13); // 9 + 4
    assert_eq!(us.at(1).data, 6); // 1 + 2 + 3

    // distinct-but-overlapping intervals are preserved
    let ov = make(&[(10, 20), (15, 25), (10, 30)]);
    let uo = ov.unique();
    assert_eq!(uo.size(), 3, "overlapping non-duplicates preserved");
}

#[test]
fn results_are_unindexed_then_buildable() {
    let mut a = IntervalMap::new();
    a.add(1, 10, 0);
    a.add(20, 30, 1);
    a.build();
    let mut b = IntervalMap::new();
    b.add(5, 25, 0);
    b.build();

    let mut inter = a.intersection(&b); // [5,10] and [20,25]
    inter.build();
    assert!(inter.has_overlaps(7, 7));
    assert!(!inter.has_overlaps(15, 15));
}
