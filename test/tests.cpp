// =============================================================================
//  SuperIntervals - test suite and usage examples
// =============================================================================
//
//  This file doubles as documentation: each section is a small, self-contained
//  example of how to use the library, checked with assert(). Reading it top to
//  bottom is a reasonable way to learn the API.
//
//  Core ideas
//  ----------
//   * si::IntervalMap<S, T> stores intervals [start, end] (END-INCLUSIVE) each
//     carrying a data payload of type T. S is the coordinate type (e.g. int).
//   * Usage is always: add() all intervals, call build() once, then query.
//     If you add more intervals after building, call build() again.
//   * Queries come in families that return different things:
//        search_values(...)  -> the data payloads (T)
//        search_idxs(...)    -> internal indices into the map
//        search_keys(...)    -> (start, end) pairs
//        search_items(...)   -> full Interval<S,T> objects
//        count(...)          -> how many intervals overlap
//        coverage(...)       -> count + total overlapping length
//        has_overlaps(...)   -> bool, is there at least one overlap
//   * Results are returned in DESCENDING start order (a property of the
//     algorithm). Reverse the output if you need ascending order.
//
//  Build and run:
//     make run-tests && ./run-tests
//
// =============================================================================

#include "superintervals.hpp"
#include <iostream>
#include <vector>
#include <cassert>
#include <utility>
#include <tuple>
#include <algorithm>

using Map = si::IntervalMap<int, int>;

// -----------------------------------------------------------------------------
//  Small helpers
// -----------------------------------------------------------------------------

// Collect the raw stored intervals of a map as sorted (start, end, data) tuples,
// so results can be compared regardless of internal ordering. Works on an
// unindexed map (build() is not required to read stored intervals).
static std::vector<std::tuple<int, int, int>> collect(Map& m) {
    std::vector<std::tuple<int, int, int>> v;
    for (size_t i = 0; i < m.size(); ++i) {
        si::Interval<int, int> it;
        m.at(i, it);
        v.emplace_back(it.start, it.end, it.data);
    }
    std::sort(v.begin(), v.end());
    return v;
}


// -----------------------------------------------------------------------------
//  1. Basics: add, build, and the query families
// -----------------------------------------------------------------------------
void test_basics() {
    std::cout << "basics, ";
    Map itv;

    // Querying before adding anything is safe and finds nothing.
    assert(!itv.has_overlaps(0, 1));

    // Add intervals as (start, end, data). Here data is a simple id.
    itv.add(10, 20, 0);
    itv.add(11, 12, -1);
    itv.add(13, 14, -1);
    itv.add(15, 16, -1);
    itv.add(25, 29, 4);
    itv.build();  // index once, before querying

    // has_overlaps: is anything overlapping [17, 30]? (yes: [10,20] and [25,29])
    assert(itv.has_overlaps(17, 30));
    assert(!itv.has_overlaps(1, 3));

    // search_values: the data payloads of overlapping intervals.
    // Note results are in descending start order: [25,29] (data 4) then [10,20] (data 0).
    std::vector<int> values;
    itv.search_values(17, 30, values);
    assert(values.size() == 2);
    assert(values[0] == 4 && values[1] == 0);

    // count: just how many overlap (cheaper than collecting them).
    assert(itv.count(17, 30) == 2);

    // search_idxs: internal indices, useful to look up parallel data yourself.
    std::vector<size_t> idxs;
    itv.search_idxs(17, 30, idxs);
    assert(idxs[0] == 4 && idxs[1] == 0);

    // search_keys: the (start, end) pairs of overlapping intervals.
    std::vector<std::pair<int, int>> keys;
    itv.search_keys(17, 30, keys);
    assert(keys[0].first == 25);    // first result's start
    assert(keys[1].second == 20);   // second result's end

    // coverage: { number of overlaps, total overlapping length } within [10, 29].
    std::pair<size_t, int> cov{0, 0};
    itv.coverage(10, 29, cov);
    assert(cov.second == 17);
}


// -----------------------------------------------------------------------------
//  2. Iterating results with a range-based for loop
// -----------------------------------------------------------------------------
void test_iteration() {
    std::cout << "iteration, ";
    Map itv;
    itv.add(1, 10, 0);
    itv.build();

    // search_idxs(start, end) returns a range you can iterate directly,
    // without allocating a result vector.
    size_t count = 0;
    int last_data = 0;
    for (const size_t i : itv.search_idxs(5, 11)) {
        last_data = itv.data[i];
        ++count;
    }
    assert(count == 1 && last_data == 0);

    // search_items(...) yields full Interval<S,T> objects.
    for (const auto& interval : itv.search_items(5, 11)) {
        assert(interval.start == 1 && interval.end == 10 && interval.data == 0);
    }
}


// -----------------------------------------------------------------------------
//  3. Overlap queries against many, varied intervals
//     (nesting, shared endpoints, large containers, duplicates)
// -----------------------------------------------------------------------------
void test_overlap_queries() {
    std::cout << "overlaps, ";
    std::vector<int> a;

    // A large container [0, 250000] plus many smaller intervals; the big one
    // overlaps almost everything, and should be the last (smallest-start) result.
    {
        Map itv;
        itv.add(0, 250000, 0);
        for (int s : {55, 115, 130, 281, 639, 842, 999, 1094, 1157, 1161,
                      1265, 1532, 1590, 1665, 1945, 2384, 2515}) {
            itv.add(s, s + 1000, -1);
        }
        itv.build();
        itv.search_values(1377, 2377, a);
        assert(a.back() == 0 && a.size() == 12);  // [0,250000] is last
        a.clear();
    }

    // Shared endpoints: an interval ending at 31 and one starting at 31 both
    // overlap the point [31, 32] because intervals are end-inclusive.
    {
        Map itv;
        itv.add(3, 40, 0);
        itv.add(4, 5, 4);
        itv.add(6, 7, 4);
        itv.add(10, 31, 5);
        itv.add(31, 32, 5);
        itv.build();
        itv.search_values(31, 32, a); assert(a.size() == 3); a.clear();  // [3,40],[10,31],[31,32]
        itv.search_values(10, 11, a); assert(a.size() == 2); a.clear();
        itv.search_values(4, 7, a);   assert(a.size() == 3); a.clear();
    }

    // Duplicate intervals are all returned.
    {
        Map itv;
        itv.add(3, 40, 0);
        itv.add(3, 40, 4);
        itv.add(3, 40, 4);
        itv.add(3, 4, 4);
        itv.add(35, 50, 4);
        itv.add(40, 400, 5);
        itv.add(40, 400, 4);
        itv.build();
        itv.search_values(38, 41, a); assert(a.size() == 6); a.clear();
        itv.search_values(41, 42, a); assert(a.size() == 3); a.clear();  // [35,50] and the two [40,400]
    }
}


// -----------------------------------------------------------------------------
//  4. coverage(): overlap count and total covered length
// -----------------------------------------------------------------------------
void test_coverage() {
    std::cout << "coverage, ";
    Map itv;
    itv.add(1, 100, 0);
    itv.add(30, 200, 7);
    itv.add(40, 50, 6);
    itv.add(60, 70, 7);
    itv.build();

    std::vector<int> a;
    itv.search_values(55, 65, a);
    assert(a.size() == 3);  // [1,100], [30,200], [60,70]
    a.clear();

    // Within the window [55, 65]: 3 intervals overlap and their clipped lengths
    // sum to 25.
    std::pair<size_t, int> cov{0, 0};
    itv.coverage(55, 65, cov);
    assert(cov.first == 3 && cov.second == 25);
}


// -----------------------------------------------------------------------------
//  5. Edge cases: empty map, single interval, rebuild after adding
// -----------------------------------------------------------------------------
void test_edge_cases() {
    std::cout << "edges, ";

    // Empty map: every query is well-defined and finds nothing.
    {
        Map itv;
        itv.build();
        assert(itv.count(1, 5) == 0);
        assert(!itv.has_overlaps(1, 5));
    }

    // Single interval.
    {
        Map itv;
        itv.add(1, 10, 1);
        itv.build();
        assert(itv.count(1, 5) == 1);
        assert(itv.count(11, 20) == 0);  // no overlap
    }
}


// =============================================================================
//  Set-like operations
//
//  Each operation returns a NEW IntervalMap that is NOT yet indexed - call
//  build() on it before querying (or before passing it as the 'other' argument
//  to another operation). This avoids paying for an index on intermediate
//  results in a chain of operations.
//
//  Notes on semantics:
//   * Intervals are end-inclusive. Merging coalesces OVERLAPPING intervals only;
//     exactly adjacent integer intervals like [1,5] and [6,10] are NOT merged.
//     (This keeps the operations valid for floating-point coordinates too.)
//   * The data payload of merged intervals is decided by an optional combiner
//     lambda; by default the first interval's data is kept.
//   * gaps / difference / symmetric_difference use +/- 1 internally and so
//     require an integral coordinate type.
// =============================================================================

void test_merge_overlaps() {
    std::cout << "merge, ";
    Map a;
    a.add(1, 5, 0);
    a.add(3, 8, 1);    // overlaps [1,5] -> coalesces into [1,8]
    a.add(20, 30, 2);  // disjoint
    a.build();

    Map merged = a.merge_overlaps();   // result is unindexed
    auto v = collect(merged);
    assert(v.size() == 2);
    assert(std::get<0>(v[0]) == 1  && std::get<1>(v[0]) == 8);    // [1,8]
    assert(std::get<0>(v[1]) == 20 && std::get<1>(v[1]) == 30);   // [20,30]

    // A combiner lambda decides the data of merged intervals.
    Map summed = a.merge_overlaps([](const int& x, const int& y) { return x + y; });
    si::Interval<int, int> it; summed.at(0, it);
    assert(it.data == 0 + 1);  // data of [1,5] and [3,8] combined
}

void test_gaps() {
    std::cout << "gaps, ";
    Map a;
    a.add(10, 20, 0);
    a.add(30, 40, 1);
    a.build();

    Map g = a.gaps(0, 50);  // uncovered regions within [0, 50]
    auto v = collect(g);
    assert(v.size() == 3);
    assert(std::get<0>(v[0]) == 0  && std::get<1>(v[0]) == 9);    // before first
    assert(std::get<0>(v[1]) == 21 && std::get<1>(v[1]) == 29);   // between
    assert(std::get<0>(v[2]) == 41 && std::get<1>(v[2]) == 50);   // after last
}

void test_union() {
    std::cout << "union, ";
    Map a, b;
    a.add(1, 10, 0);
    b.add(5, 25, 1);
    a.build(); b.build();

    Map u = a.union_with(b);  // overlapping -> covers [1, 25]
    auto v = collect(u);
    assert(v.size() == 1);
    assert(std::get<0>(v[0]) == 1 && std::get<1>(v[0]) == 25);
}

void test_intersection() {
    std::cout << "intersection, ";
    Map a, b;
    a.add(1, 10, 0);
    a.add(20, 30, 1);
    b.add(5, 25, 2);
    a.build();
    b.build();  // 'other' (b) is queried by index, so it must be built

    Map inter = a.intersection(b);  // overlapping regions: [5,10] and [20,25]
    auto v = collect(inter);
    assert(v.size() == 2);
    assert(std::get<0>(v[0]) == 5  && std::get<1>(v[0]) == 10);
    assert(std::get<0>(v[1]) == 20 && std::get<1>(v[1]) == 25);

    // The result is unindexed; build it and it becomes fully queryable.
    inter.build();
    assert(inter.has_overlaps(7, 7));
    assert(!inter.has_overlaps(15, 15));
}

void test_difference() {
    std::cout << "difference, ";
    Map a, b;
    a.add(1, 10, 0);
    b.add(4, 6, 1);
    a.build(); b.build();

    Map d = a.difference(b);  // A minus B -> [1,3] and [7,10]
    auto v = collect(d);
    assert(v.size() == 2);
    assert(std::get<0>(v[0]) == 1 && std::get<1>(v[0]) == 3);
    assert(std::get<0>(v[1]) == 7 && std::get<1>(v[1]) == 10);
}

void test_symmetric_difference() {
    std::cout << "symmetric_difference, ";
    Map a, b;
    a.add(1, 10, 0);
    b.add(5, 15, 1);
    a.build(); b.build();

    Map x = a.symmetric_difference(b);  // in exactly one set -> [1,4] and [11,15]
    auto v = collect(x);
    assert(v.size() == 2);
    assert(std::get<0>(v[0]) == 1  && std::get<1>(v[0]) == 4);
    assert(std::get<0>(v[1]) == 11 && std::get<1>(v[1]) == 15);
}

void test_span() {
    std::cout << "span, ";
    Map a;
    a.add(10, 20, 0);
    a.add(5, 8, 1);
    a.add(15, 50, 2);
    a.build();

    std::pair<int, int> sp;
    assert(a.span(sp));                       // smallest enclosing interval
    assert(sp.first == 5 && sp.second == 50);

    Map empty;
    assert(!empty.span(sp));  // returns false on an empty set
}


// -----------------------------------------------------------------------------
//  Runner
// -----------------------------------------------------------------------------
int main() {
    std::cout << "\nSuperIntervals tests\n  ";
    test_basics();
    test_iteration();
    test_overlap_queries();
    test_coverage();
    test_edge_cases();
    std::cout << "\n  All query tests passed\n";

    std::cout << "\nSet operation tests\n  ";
    test_merge_overlaps();
    test_gaps();
    test_union();
    test_intersection();
    test_difference();
    test_symmetric_difference();
    test_span();
    std::cout << "\n  All set operation tests passed\n";

    return 0;
}
