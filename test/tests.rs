// Integration tests for counting, sorting, iteration and set operations.
use superintervals::IntervalMap;

#[test]
fn count_matches_brute_force_at_block_boundaries() {
    for n in 0..=129 {
        let mut index = IntervalMap::new();
        let mut refs = Vec::new();
        for i in 0..n {
            let interval = match i % 5 {
                0 => (i32::MIN, i32::MAX),
                1 => (-100, (i % 31) as i32 - 15),
                2 => (0, 0),
                3 => ((i % 29) as i32 - 14, 100),
                _ => ((i % 17) as i32, (i % 17) as i32),
            };
            refs.push(interval);
            index.add(interval.0, interval.1, i);
        }
        index.build();
        for a in [i32::MIN, -101, -100, -15, -1, 0, 1, 15, 100, 101, i32::MAX] {
            for b in [i32::MIN, -101, -100, -15, -1, 0, 1, 15, 100, 101, i32::MAX] {
                if a > b {
                    continue;
                }
                let expected = refs.iter().filter(|&&(s, e)| s <= b && e >= a).count();
                assert_eq!(index.count(a, b), expected, "n={n}, query=[{a},{b}]");
            }
        }
    }
}

#[test]
fn hybrid_sort_preserves_payloads_at_threshold_and_signed_extremes() {
    use std::cell::Cell;
    use std::rc::Rc;
    struct Payload {
        id: usize,
        live: Rc<Cell<usize>>,
        clones: Rc<Cell<usize>>,
    }
    impl Clone for Payload {
        fn clone(&self) -> Self {
            self.live.set(self.live.get() + 1);
            self.clones.set(self.clones.get() + 1);
            Self {
                id: self.id,
                live: self.live.clone(),
                clones: self.clones.clone(),
            }
        }
    }
    impl Drop for Payload {
        fn drop(&mut self) {
            self.live.set(self.live.get() - 1);
        }
    }
    for n in [0, 1, 32, 4095, 4096, 4097, 8192] {
        let live = Rc::new(Cell::new(0));
        let clones = Rc::new(Cell::new(0));
        let mut expected = Vec::new();
        let mut map = IntervalMap::new();
        let mut state = 1729u64;
        for id in 0..n {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let a = state as i32;
            let b = (state >> 32) as i32;
            let (start, end) = match id % 13 {
                0 => (i32::MIN, i32::MAX),
                1 => (0, 8),
                2 => (0, 0),
                3 => (i32::MAX, i32::MAX),
                _ => (a.min(b), a.max(b)),
            };
            expected.push((start, end, id));
            live.set(live.get() + 1);
            map.add(
                start,
                end,
                Payload {
                    id,
                    live: live.clone(),
                    clones: clones.clone(),
                },
            );
        }
        expected.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| b.1.cmp(&a.1)));
        map.build();
        assert_eq!(live.get(), n);
        assert_eq!(clones.get(), 0, "full sorting should move owned payloads");
        let actual: Vec<_> = (0..n)
            .map(|i| (map.starts[i], map.ends[i], map.data[i].id))
            .collect();
        assert_eq!(actual, expected);
        for (a, b) in [(i32::MIN, i32::MIN), (-1, 0), (0, 8), (i32::MAX, i32::MAX)] {
            let want = expected
                .iter()
                .filter(|&&(s, e, _)| s <= b && a <= e)
                .count();
            assert_eq!(map.count(a, b), want);
            assert_eq!(map.search_values_iter(a, b).count(), want);
        }
        drop(map);
        assert_eq!(
            live.get(),
            0,
            "every owned payload must be dropped exactly once"
        );
    }
    // Wide payloads use the comparison fallback, with the same coordinate order.
    let mut map = IntervalMap::new();
    for id in (0..5000u64).rev() {
        map.add(id as i32, id as i32 + 3, [id; 16]);
    }
    map.build();
    for i in 0..5000 {
        assert_eq!(map.starts[i], i as i32);
        assert_eq!(map.data[i], [i as u64; 16]);
    }
}

#[test]
fn build_orders_starts_then_longest_end_and_preserves_duplicates() {
    let intervals = [
        (0, 5, 0u32),
        (i32::MIN, i32::MAX, 1),
        (0, 8, 2),
        (i32::MAX, i32::MAX, 3),
        (0, 8, 4),
        (0, 0, 5),
        (i32::MIN, i32::MIN, 6),
        (-1, 1, 7),
    ];
    // Exercise full sorting and the equal-start-group fast path.
    for order in [vec![0, 1, 2, 3, 4, 5, 6, 7], vec![6, 1, 7, 5, 0, 2, 4, 3]] {
        let mut expected: Vec<_> = order.iter().map(|&i| intervals[i]).collect();
        expected.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| b.1.cmp(&a.1)));
        let mut map = IntervalMap::new();
        for &i in &order {
            let (a, b, id) = intervals[i];
            map.add(a, b, id);
        }
        map.build();
        let actual: Vec<_> = (0..map.size())
            .map(|i| (map.starts[i], map.ends[i], map.data[i]))
            .collect();
        assert_eq!(actual, expected);
        for (a, b) in [(i32::MIN, i32::MIN), (-1, 0), (0, 8), (i32::MAX, i32::MAX)] {
            let mut want: Vec<_> = intervals
                .iter()
                .filter(|&&(s, e, _)| s <= b && a <= e)
                .map(|&(_, _, id)| id)
                .collect();
            let mut got: Vec<_> = map.search_values_iter(a, b).copied().collect();
            want.sort_unstable();
            got.sort_unstable();
            assert_eq!(got, want);
            assert_eq!(map.count(a, b), want.len());
        }
    }
}

#[test]
fn branch_links_match_nearest_prior_end_after_build_and_rebuild() {
    fn check<T: Clone>(map: &IntervalMap<T>) {
        assert_eq!(map.branch.len(), map.ends.len());
        for i in 0..map.ends.len() {
            let want = (0..i)
                .rev()
                .find(|&j| map.ends[j] >= map.ends[i])
                .unwrap_or(usize::MAX);
            assert_eq!(map.branch[i], want, "branch at {i}");
        }
    }
    for n in [0, 1, 2, 31, 32, 33, 128, 257] {
        for pattern in 0..4 {
            let mut map = IntervalMap::new();
            for i in 0..n {
                let end = match pattern {
                    0 => i,
                    1 => 2 * n - i,
                    2 => 2 * n,
                    _ => i + (i * 37 % 101),
                };
                map.add(i as i32, end as i32, i);
            }
            map.build();
            check(&map);
            map.add(-1, 10000, n);
            map.build();
            check(&map);
            map.add(-2, -2, n + 1);
            map.build();
            check(&map);
            map.clear();
            map.build();
            check(&map);
        }
    }
}

use std::cell::Cell;
use std::rc::Rc;

struct Payload {
    id: usize,
    text: String,
    clones: Rc<Cell<usize>>,
}
impl Clone for Payload {
    fn clone(&self) -> Self {
        self.clones.set(self.clones.get() + 1);
        Self {
            id: self.id,
            text: self.text.clone(),
            clones: self.clones.clone(),
        }
    }
}

#[test]
fn borrowed_values_preserve_order_and_identity_without_cloning() {
    let clones = Rc::new(Cell::new(0));
    let mut map = IntervalMap::new();
    for (id, (a, b)) in [(10, 20), (0, 100), (10, 20), (30, 30)]
        .into_iter()
        .enumerate()
    {
        map.add(
            a,
            b,
            Payload {
                id,
                text: format!("payload-{id}"),
                clones: clones.clone(),
            },
        );
    }
    map.build();
    clones.set(0);
    for (a, b) in [(-10, -1), (10, 20), (30, 30), (101, 200), (0, 100)] {
        let expected: Vec<_> = map.search_idxs_iter(a, b).collect();
        let borrowed: Vec<_> = map.search_values_iter(a, b).collect();
        assert_eq!(borrowed.len(), expected.len());
        for (&index, value) in expected.iter().zip(&borrowed) {
            assert!(std::ptr::eq(*value, &map.data[index]));
            assert_eq!(value.text, format!("payload-{}", value.id));
        }
        assert_eq!(clones.get(), 0);
        let owned: Vec<_> = map.search_values_iter(a, b).cloned().collect();
        assert_eq!(
            owned.iter().map(|v| v.id).collect::<Vec<_>>(),
            borrowed.iter().map(|v| v.id).collect::<Vec<_>>()
        );
        assert_eq!(clones.get(), owned.len());
        clones.set(0);
    }
    let mut values = map.search_values_iter(10, 20);
    let first = values.next().unwrap();
    let address = first as *const Payload;
    while values.next().is_some() {}
    assert!(values.next().is_none());
    assert!(values.next().is_none());
    assert_eq!(first as *const Payload, address); // borrowed values survive iterator advancement
    assert_eq!(clones.get(), 0);
    let empty = IntervalMap::<String>::new();
    assert!(empty.search_values_iter(0, 100).next().is_none());
}

// Set-operation tests for IntervalMap, checked against a brute-force point-set
// reference over a small integer universe.

use std::collections::BTreeSet;

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
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
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
