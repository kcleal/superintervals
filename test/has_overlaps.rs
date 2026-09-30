// has_overlaps tests, including a randomized check against search_values.

use superintervals::IntervalMap;

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

fn make(intervals: &[(i32, i32)]) -> IntervalMap<i32> {
    let mut m = IntervalMap::new();
    for (i, &(s, e)) in intervals.iter().enumerate() {
        m.add(s, e, i as i32);
    }
    m.build();
    m
}

#[test]
fn has_overlaps_finds_long_interval_before_nested_one() {
    let mut m = IntervalMap::new();
    m.add(0, 99, "long");
    m.add(10, 19, "short");
    m.build();
    let mut found = Vec::new();
    m.search_values(30, 39, &mut found);
    assert_eq!(found, vec!["long"]);
    assert!(m.has_overlaps(30, 39));
}

#[test]
fn has_overlaps_on_empty_map() {
    let mut m: IntervalMap<i32> = IntervalMap::new();
    assert!(!m.has_overlaps(0, 0));
    m.build();
    assert!(!m.has_overlaps(0, 0));
}

#[test]
fn has_overlaps_boundaries() {
    let single = make(&[(10, 20)]);
    assert!(single.has_overlaps(20, 25), "query starts on the interval end");
    assert!(single.has_overlaps(5, 10), "query ends on the interval start");
    assert!(single.has_overlaps(15, 15));
    assert!(single.has_overlaps(0, 30));
    assert!(!single.has_overlaps(21, 25));
    assert!(!single.has_overlaps(5, 9));

    let nested = make(&[(0, 99), (10, 19)]);
    assert!(nested.has_overlaps(99, 120));
    assert!(nested.has_overlaps(20, 20));
    assert!(nested.has_overlaps(5, 5));
    assert!(!nested.has_overlaps(100, 120));
    assert!(!nested.has_overlaps(-10, -1));
}

#[test]
fn has_overlaps_matches_search_values() {
    let mut rng = Lcg(42);
    let lengths = [5, 50, 500];
    let mut found = Vec::new();
    for _ in 0..500 {
        let n = rng.range(60);
        let mut intervals = Vec::new();
        for _ in 0..n {
            let s = rng.range(1000);
            let max_len = lengths[rng.range(3) as usize];
            intervals.push((s, s + rng.range(max_len)));
        }
        let m = make(&intervals);
        for _ in 0..200 {
            let a = rng.range(1150) - 50;
            let b = a + rng.range(100);
            found.clear();
            m.search_values(a, b, &mut found);
            assert_eq!(m.has_overlaps(a, b), !found.is_empty(), "query ({}, {}) on {:?}", a, b, intervals);
        }
    }
}
