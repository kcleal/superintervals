// Version 1.0.0
/**
 * @file c_superintervals.h
 * @brief C implementation of the SuperIntervals data structure for interval intersection queries.
 *
 * SuperIntervals stores a set of (start, end, data) intervals and, after indexing,
 * answers overlap queries against a query range very quickly. Intervals are kept in
 * a flat start-sorted array alongside a single auxiliary "branch" array that lets a
 * backward walk skip over runs of intervals that cannot overlap the query.
 *
 * Typical usage:
 * @code
 *     cSuperIntervals* si = createSuperIntervals();
 *     addInterval(si, 10, 20, 1);
 *     addInterval(si, 15, 25, 2);
 *     indexSuperIntervals(si);          // must be called before querying
 *
 *     cIndexResult res = createIndexResult();
 *     searchIdxs(si, 8, 18, &res);      // res.data holds matching indices
 *     destroyIndexResult(&res);
 *
 *     destroySuperIntervals(si);
 * @endcode
 *
 * @note Intervals are end-INCLUSIVE: [start, end] covers both endpoints.
 * @note indexSuperIntervals() must be called before any query. If more intervals are
 *       added afterwards, call it again.
 * @note Query results are returned in REVERSE position-sorted order (descending start),
 *       which falls naturally out of the backward walk.
 * @note Coordinates are int32_t and data is int32_t (an integer key/handle). Unlike the
 *       C++ template, the C API does not carry an arbitrary payload type.
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef SUPERINTERVALS_HEADER_INCLUDED
#define SUPERINTERVALS_HEADER_INCLUDED 1

/** Library version string. */
#define SUPERINTERVALS_C_VERSION "1.0.0"

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>

#ifdef __AVX2__
    #include <immintrin.h>
#elif defined __ARM_NEON
    #include <arm_neon.h>
#endif

/** Sentinel index meaning "no interval" / "walked off the front of the array". */
#define SI_NONE ((size_t)-1)

/**
 * @brief A single interval: end-inclusive [start, end] with an associated integer value.
 */
typedef struct {
    int32_t start;
    int32_t end;
    int32_t data;
} Interval;

/**
 * @brief A (start, end) coordinate pair, used by searchKeys and the difference operation.
 */
typedef struct {
    int32_t start;
    int32_t end;
} KeyPair;

/**
 * @brief The SuperIntervals container.
 *
 * Fields are exposed for inspection but should be mutated only through the API
 * functions, which maintain the startSorted / endSorted invariants.
 */
typedef struct {
    int32_t* starts;     /**< interval start coordinates */
    int32_t* ends;       /**< interval end coordinates */
    int32_t* data;       /**< interval data values */
    size_t* branch;      /**< skip-ahead index per interval (built by indexSuperIntervals) */
    size_t size;         /**< number of stored intervals */
    size_t capacity;     /**< allocated capacity */
    size_t idx;          /**< scratch cursor used by upperBound */
    bool startSorted;    /**< true while intervals remain in non-decreasing start order */
    bool endSorted;      /**< true while ties on start remain in non-increasing end order */
} cSuperIntervals;

/**
 * @brief A growable result buffer of int32_t values (used for data values, indices).
 *
 * Reused across queries: call clearResult() to reset length without freeing memory.
 */
typedef struct {
    int32_t* data;       /**< result values */
    size_t size;         /**< number of valid entries */
    size_t capacity;     /**< allocated capacity */
} cIndexResult;

/**
 * @brief A growable result buffer of KeyPair values (used by searchKeys).
 */
typedef struct {
    KeyPair* data;
    size_t size;
    size_t capacity;
} cKeyResult;

/**
 * @brief A growable result buffer of Interval values (used by searchItems).
 */
typedef struct {
    Interval* data;
    size_t size;
    size_t capacity;
} cItemResult;

/**
 * @brief A combiner callback: decides the resulting data value when two intervals merge.
 *
 * Receives the accumulated value and the next value; returns the value to keep.
 * Pass NULL to any set operation to use the default (keep the first / left value).
 */
typedef int32_t (*cCombineFn)(int32_t accumulated, int32_t next);

/* ----------------------------------------------------------------------------
 * Lifecycle
 * ------------------------------------------------------------------------- */

/** @brief Allocate and initialise an empty SuperIntervals. @return Owning pointer; free with destroySuperIntervals. */
cSuperIntervals* createSuperIntervals(void);

/** @brief Free a SuperIntervals and all its internal buffers. */
void destroySuperIntervals(cSuperIntervals* si);

/** @brief Remove all intervals but keep allocated capacity. Re-index before querying. */
void clearSuperIntervals(cSuperIntervals* si);

/** @brief Ensure capacity for at least n intervals (does not change size). */
void reserveSuperIntervals(cSuperIntervals* si, size_t n);

/**
 * @brief Append an interval. Intervals are end-inclusive.
 * @param start Interval start.
 * @param end   Interval end (inclusive).
 * @param value Associated integer value.
 * @note Re-index with indexSuperIntervals() before querying.
 */
void addInterval(cSuperIntervals* si, int32_t start, int32_t end, int32_t value);

/** @brief Number of intervals currently stored. */
size_t sizeSuperIntervals(const cSuperIntervals* si);

/* ----------------------------------------------------------------------------
 * Indexing
 * ------------------------------------------------------------------------- */

/** @brief Sort intervals into the canonical (start asc, end desc) order. Called by indexSuperIntervals. */
void sortIntervals(cSuperIntervals* si);

/**
 * @brief Build the branch index. MUST be called after adding intervals and before any query.
 *        Re-call if more intervals are added afterwards.
 */
void indexSuperIntervals(cSuperIntervals* si);

/* ----------------------------------------------------------------------------
 * Element access
 * ------------------------------------------------------------------------- */

/** @brief Copy the interval at the given index into *out. @return true on success, false if out of range. */
bool intervalAt(const cSuperIntervals* si, size_t index, Interval* out);

/** @brief Start coordinate of the interval at index (after indexing this is sorted order). */
int32_t startAt(const cSuperIntervals* si, size_t index);

/** @brief End coordinate of the interval at index. */
int32_t endAt(const cSuperIntervals* si, size_t index);

/** @brief Data value of the interval at index. */
int32_t dataAt(const cSuperIntervals* si, size_t index);

/* ----------------------------------------------------------------------------
 * Internal cursor helper (exposed for advanced use)
 * ------------------------------------------------------------------------- */

/**
 * @brief Position the internal cursor at the last interval whose start is <= value.
 * @return The found index, or SI_NONE if no such interval exists. Also stored in si->idx.
 */
size_t upperBound(cSuperIntervals* si, int32_t value);

/* ----------------------------------------------------------------------------
 * Queries
 * ------------------------------------------------------------------------- */

/** @brief True if any stored interval overlaps [start, end]. */
bool anyOverlaps(cSuperIntervals* si, int32_t start, int32_t end);

/** @brief Count intervals overlapping [start, end]. */
size_t countOverlaps(cSuperIntervals* si, int32_t start, int32_t end);

/**
 * @brief Collect the DATA values of intervals overlapping [start, end].
 * @param found Result buffer (appended to; clear it first if you want fresh results).
 */
void searchValues(cSuperIntervals* si, int32_t start, int32_t end, cIndexResult* found);

/** @brief Collect the INDICES of intervals overlapping [start, end]. */
void searchIdxs(cSuperIntervals* si, int32_t start, int32_t end, cIndexResult* found);

/** @brief Collect the (start, end) KEY pairs of intervals overlapping [start, end]. */
void searchKeys(cSuperIntervals* si, int32_t start, int32_t end, cKeyResult* found);

/** @brief Collect the full Interval records overlapping [start, end]. */
void searchItems(cSuperIntervals* si, int32_t start, int32_t end, cItemResult* found);

/** @brief Collect data values of all intervals containing the single point. */
void searchPoint(cSuperIntervals* si, int32_t point, cIndexResult* found);

/**
 * @brief Count overlaps of [start, end] and sum the overlapping (clipped) lengths.
 * @param count_out  Receives the number of overlapping intervals.
 * @param coverage_out Receives the total covered length within [start, end].
 */
void coverage(cSuperIntervals* si, int32_t start, int32_t end,
              size_t* count_out, int32_t* coverage_out);

/**
 * @brief Legacy alias for searchValues into a raw buffer.
 * @param found      Caller-allocated array large enough to hold the results.
 * @param found_size Receives the number of values written.
 * @deprecated Prefer searchValues with a cIndexResult, which grows automatically.
 */
void findOverlaps(cSuperIntervals* si, int32_t start, int32_t end,
                  int32_t* found, size_t* found_size);

/* ----------------------------------------------------------------------------
 * Set operations
 *
 * Each operation returns a NEW, NOT-yet-indexed cSuperIntervals (owning pointer;
 * free with destroySuperIntervals). Call indexSuperIntervals() on the result before
 * querying it, or before passing it as 'other' to another set operation that queries
 * an index (intersection / difference / symmetric_difference).
 *
 * Intervals are end-inclusive. Two intervals overlap when a.start <= b.end and
 * b.start <= a.end. Merging coalesces overlapping intervals only; exactly adjacent
 * integer intervals (e.g. [1,5] and [6,10]) are NOT merged. Operations that use +/-1
 * (gaps / difference / symmetric_difference) assume integer coordinates.
 *
 * The 'combine' callback resolves the data value when intervals fold together; pass
 * NULL to keep the first/left value.
 * ------------------------------------------------------------------------- */

/**
 * @brief Coalesce this set's intervals into a disjoint, non-overlapping set.
 * @param combine Folds data of merged intervals, or NULL to keep the first.
 * @return New unindexed cSuperIntervals (caller owns).
 */
cSuperIntervals* mergeOverlaps(const cSuperIntervals* si, cCombineFn combine);

/**
 * @brief Compute the uncovered gaps within [lo, hi].
 * @param fill Data value assigned to each gap interval.
 * @return New unindexed cSuperIntervals (caller owns). Requires integer coordinates.
 */
cSuperIntervals* intervalGaps(const cSuperIntervals* si, int32_t lo, int32_t hi, int32_t fill);

/**
 * @brief Union of this set with other (all covered regions, coalesced).
 * @param combine Folds data of merged intervals, or NULL to keep the first.
 * @return New unindexed cSuperIntervals (caller owns).
 */
cSuperIntervals* unionWith(const cSuperIntervals* si, const cSuperIntervals* other, cCombineFn combine);

/**
 * @brief Intersection of this set with other (regions covered by both).
 * @param other  Must already be indexed (it is queried via its branch index).
 * @param combine Produces output data from each overlapping pair (a_data, b_data), or NULL to keep a_data.
 * @return New unindexed cSuperIntervals (caller owns). Pieces are NOT coalesced.
 */
cSuperIntervals* intersection(const cSuperIntervals* si, cSuperIntervals* other, cCombineFn combine);

/**
 * @brief Difference A \ B: regions of this set not covered by other.
 * @param other Must already be indexed. @return New unindexed cSuperIntervals (caller owns). Integer coordinates only.
 */
cSuperIntervals* difference(const cSuperIntervals* si, cSuperIntervals* other);

/**
 * @brief Symmetric difference: regions in exactly one of the two sets.
 * @param other Must already be indexed. @return New unindexed cSuperIntervals (caller owns). Integer coordinates only.
 */
cSuperIntervals* symmetricDifference(cSuperIntervals* si, cSuperIntervals* other);

/**
 * @brief Smallest interval enclosing every stored interval (the span / convex hull).
 * @param lo_out Receives the minimum start, @param hi_out the maximum end.
 * @return true if the set is non-empty (outputs written), false otherwise.
 */
bool intervalSpan(const cSuperIntervals* si, int32_t* lo_out, int32_t* hi_out);

/**
 * @brief Grow (or shrink) every interval, like `bedtools slop`.
 * @param left  Subtracted from each start (extends left); negative shrinks.
 * @param right Added to each end (extends right); negative shrinks.
 * @param lo,hi Clamp the resulting coordinates. Intervals that shrink past themselves are dropped.
 * @return New unindexed cSuperIntervals (caller owns). Data is carried over unchanged.
 */
cSuperIntervals* expandIntervals(const cSuperIntervals* si, int32_t left, int32_t right,
                                 int32_t lo, int32_t hi);

/**
 * @brief Create flanking intervals beside each interval, like `bedtools flank`.
 * @param left  Width of the left flank [start-left, start-1] (0 = none).
 * @param right Width of the right flank [end+1, end+right] (0 = none).
 * @param lo,hi Clamp the flanks. Emits the flanks only (NOT the originals); each inherits its source data.
 * @return New unindexed cSuperIntervals (caller owns).
 */
cSuperIntervals* flankIntervals(const cSuperIntervals* si, int32_t left, int32_t right,
                                int32_t lo, int32_t hi);

/**
 * @brief Keep only one interval per distinct (start, end) coordinate pair.
 * @param combine Folds data of exact duplicates, or NULL to keep the first.
 * @return New unindexed cSuperIntervals (caller owns). Distinct-but-overlapping intervals are preserved.
 */
cSuperIntervals* uniqueIntervals(const cSuperIntervals* si, cCombineFn combine);

/* ----------------------------------------------------------------------------
 * Result buffer helpers
 * ------------------------------------------------------------------------- */

/** @brief Initialise an empty cIndexResult. */
cIndexResult createIndexResult(void);
/** @brief Reset a cIndexResult to empty without freeing its buffer. */
void clearIndexResult(cIndexResult* r);
/** @brief Free a cIndexResult's buffer. */
void destroyIndexResult(cIndexResult* r);

/** @brief Initialise an empty cKeyResult. */
cKeyResult createKeyResult(void);
/** @brief Reset a cKeyResult to empty without freeing its buffer. */
void clearKeyResult(cKeyResult* r);
/** @brief Free a cKeyResult's buffer. */
void destroyKeyResult(cKeyResult* r);

/** @brief Initialise an empty cItemResult. */
cItemResult createItemResult(void);
/** @brief Reset a cItemResult to empty without freeing its buffer. */
void clearItemResult(cItemResult* r);
/** @brief Free a cItemResult's buffer. */
void destroyItemResult(cItemResult* r);

/* ============================================================================
 * Implementation
 * ========================================================================= */

cSuperIntervals* createSuperIntervals(void) {
    cSuperIntervals* si = (cSuperIntervals*)malloc(sizeof(cSuperIntervals));
    si->starts = NULL;
    si->ends = NULL;
    si->data = NULL;
    si->branch = NULL;
    si->size = 0;
    si->capacity = 0;
    si->idx = 0;
    si->startSorted = true;
    si->endSorted = true;
    return si;
}

void destroySuperIntervals(cSuperIntervals* si) {
    if (!si) return;
    free(si->starts);
    free(si->ends);
    free(si->data);
    free(si->branch);
    free(si);
}

void clearSuperIntervals(cSuperIntervals* si) {
    si->size = 0;
    si->idx = 0;
    si->startSorted = true;
    si->endSorted = true;
}

void reserveSuperIntervals(cSuperIntervals* si, size_t n) {
    if (n > si->capacity) {
        si->capacity = n;
        si->starts = (int32_t*)realloc(si->starts, n * sizeof(int32_t));
        si->ends = (int32_t*)realloc(si->ends, n * sizeof(int32_t));
        si->data = (int32_t*)realloc(si->data, n * sizeof(int32_t));
    }
}

void addInterval(cSuperIntervals* si, int32_t start, int32_t end, int32_t value) {
    if (si->size >= si->capacity) {
        size_t new_capacity = si->capacity == 0 ? 1 : si->capacity * 2;
        reserveSuperIntervals(si, new_capacity);
    }

    if (si->startSorted && si->size > 0) {
        si->startSorted = (start < si->starts[si->size - 1]) ? false : true;
        if (si->startSorted && start == si->starts[si->size - 1] && end > si->ends[si->size - 1]) {
            si->endSorted = false;
        }
    }

    si->starts[si->size] = start;
    si->ends[si->size] = end;
    si->data[si->size] = value;
    si->size++;
}

size_t sizeSuperIntervals(const cSuperIntervals* si) {
    return si->size;
}

static int compareIntervalsStart(const void* a, const void* b) {
    const Interval* ia = (const Interval*)a;
    const Interval* ib = (const Interval*)b;
    if (ia->start != ib->start) {
        return (ia->start < ib->start) ? -1 : 1;
    }
    return (ib->end < ia->end) ? -1 : (ib->end > ia->end);  // end descending
}

static int compareIntervalsEnd(const void* a, const void* b) {
    const Interval* ia = (const Interval*)a;
    const Interval* ib = (const Interval*)b;
    return (ib->end < ia->end) ? -1 : (ib->end > ia->end);  // end descending
}

static void sortBlock(cSuperIntervals* si, size_t start_i, size_t end_i,
        int (*compare)(const void*, const void*)) {
    size_t range_size = end_i - start_i;
    Interval* tmp = (Interval*)malloc(range_size * sizeof(Interval));

    for (size_t i = 0; i < range_size; ++i) {
        tmp[i].start = si->starts[start_i + i];
        tmp[i].end = si->ends[start_i + i];
        tmp[i].data = si->data[start_i + i];
    }

    qsort(tmp, range_size, sizeof(Interval), compare);

    for (size_t i = 0; i < range_size; ++i) {
        si->starts[start_i + i] = tmp[i].start;
        si->ends[start_i + i] = tmp[i].end;
        si->data[start_i + i] = tmp[i].data;
    }

    free(tmp);
}

void sortIntervals(cSuperIntervals* si) {
    if (!si->startSorted) {
        sortBlock(si, 0, si->size, compareIntervalsStart);
        si->startSorted = true;
        si->endSorted = true;
    } else if (!si->endSorted) {
        size_t it_start = 0;
        while (it_start < si->size) {
            size_t block_end = it_start + 1;
            bool needs_sort = false;
            while (block_end < si->size && si->starts[block_end] == si->starts[it_start]) {
                if (block_end > it_start && si->ends[block_end] > si->ends[block_end - 1]) {
                    needs_sort = true;
                }
                ++block_end;
            }
            if (needs_sort) {
                sortBlock(si, it_start, block_end, compareIntervalsEnd);
            }
            it_start = block_end;
        }
        si->endSorted = true;
    }
}

void indexSuperIntervals(cSuperIntervals* si) {
    if (si->size == 0) {
        return;
    }

    sortIntervals(si);

    si->branch = (size_t*)realloc(si->branch, si->size * sizeof(size_t));
    memset(si->branch, -1, si->size * sizeof(size_t));

    // branch[i] = nearest j < i with ends[j] >= ends[i], via a monotonic stack.
    size_t* br = (size_t*)malloc(si->size * sizeof(size_t));
    int32_t* br_ends = (int32_t*)malloc(si->size * sizeof(int32_t));
    size_t br_size = 0;

    br[br_size] = 0;
    br_ends[br_size] = si->ends[0];
    br_size++;

    for (size_t i = 1; i < si->size; ++i) {
        while (br_size > 0 && br_ends[br_size - 1] < si->ends[i]) {
            br_size--;
        }
        if (br_size > 0) {
            si->branch[i] = br[br_size - 1];
        }
        br[br_size] = i;
        br_ends[br_size] = si->ends[i];
        br_size++;
    }

    free(br);
    free(br_ends);
    si->idx = 0;
}

bool intervalAt(const cSuperIntervals* si, size_t index, Interval* out) {
    if (index >= si->size) {
        return false;
    }
    out->start = si->starts[index];
    out->end = si->ends[index];
    out->data = si->data[index];
    return true;
}

int32_t startAt(const cSuperIntervals* si, size_t index) { return si->starts[index]; }
int32_t endAt(const cSuperIntervals* si, size_t index)   { return si->ends[index]; }
int32_t dataAt(const cSuperIntervals* si, size_t index)  { return si->data[index]; }

size_t upperBound(cSuperIntervals* si, int32_t value) {
    if (si->size == 0) {
        si->idx = SI_NONE;
        return SI_NONE;
    }
    size_t length = si->size;
    size_t idx = 0;
    const int entries_per_256KB = 256 * 1024 / sizeof(int32_t);
    const int num_per_cache_line = 64 / sizeof(int32_t) > 1 ? 64 / sizeof(int32_t) : 1;

    if (length >= (size_t)entries_per_256KB) {
        while (length >= (size_t)(3 * num_per_cache_line)) {
            size_t half = length / 2;
            idx += (si->starts[idx + half] <= value) * (length - half);
            length = half;
        }
    }
    while (length > 0) {
        size_t half = length / 2;
        idx += (si->starts[idx + half] <= value) * (length - half);
        length = half;
    }

    // idx is now the count of starts <= value; convert to the last such index.
    if (idx == 0) {
        si->idx = SI_NONE;   // no interval starts at or before value
        return SI_NONE;
    }
    idx -= 1;
    si->idx = idx;
    return idx;
}

bool anyOverlaps(cSuperIntervals* si, int32_t start, int32_t end) {
    size_t idx = upperBound(si, end);
    return idx != SI_NONE && start <= si->ends[idx];
}

void searchValues(cSuperIntervals* si, int32_t start, int32_t end, cIndexResult* found) {
    if (si->size == 0) {
        return;
    }
    size_t idx = upperBound(si, end);
    if (idx == SI_NONE) {
        return;
    }
    size_t i = idx;
    while (i != SI_NONE && start <= si->ends[i]) {
        if (found->size == found->capacity) {
            found->capacity = found->capacity ? found->capacity * 2 : 16;
            found->data = (int32_t*)realloc(found->data, found->capacity * sizeof(int32_t));
        }
        found->data[found->size++] = si->data[i];
        --i;
    }
    if (i == SI_NONE) {
        return;
    }
    i = si->branch[i];
    while (i != SI_NONE) {
        if (start <= si->ends[i]) {
            if (found->size == found->capacity) {
                found->capacity = found->capacity ? found->capacity * 2 : 16;
                found->data = (int32_t*)realloc(found->data, found->capacity * sizeof(int32_t));
            }
            found->data[found->size++] = si->data[i];
            --i;
        } else {
            i = si->branch[i];
        }
    }
}

void searchIdxs(cSuperIntervals* si, int32_t start, int32_t end, cIndexResult* found) {
    if (si->size == 0) {
        return;
    }
    size_t idx = upperBound(si, end);
    if (idx == SI_NONE) {
        return;
    }
    size_t i = idx;
    while (i != SI_NONE && start <= si->ends[i]) {
        if (found->size == found->capacity) {
            found->capacity = found->capacity ? found->capacity * 2 : 16;
            found->data = (int32_t*)realloc(found->data, found->capacity * sizeof(int32_t));
        }
        found->data[found->size++] = (int32_t)i;
        --i;
    }
    if (i == SI_NONE) {
        return;
    }
    i = si->branch[i];
    while (i != SI_NONE) {
        if (start <= si->ends[i]) {
            if (found->size == found->capacity) {
                found->capacity = found->capacity ? found->capacity * 2 : 16;
                found->data = (int32_t*)realloc(found->data, found->capacity * sizeof(int32_t));
            }
            found->data[found->size++] = (int32_t)i;
            --i;
        } else {
            i = si->branch[i];
        }
    }
}

void searchKeys(cSuperIntervals* si, int32_t start, int32_t end, cKeyResult* found) {
    if (si->size == 0) {
        return;
    }
    size_t idx = upperBound(si, end);
    if (idx == SI_NONE) {
        return;
    }
    size_t i = idx;
    while (i != SI_NONE && start <= si->ends[i]) {
        if (found->size == found->capacity) {
            found->capacity = found->capacity ? found->capacity * 2 : 16;
            found->data = (KeyPair*)realloc(found->data, found->capacity * sizeof(KeyPair));
        }
        found->data[found->size].start = si->starts[i];
        found->data[found->size].end = si->ends[i];
        found->size++;
        --i;
    }
    if (i == SI_NONE) {
        return;
    }
    i = si->branch[i];
    while (i != SI_NONE) {
        if (start <= si->ends[i]) {
            if (found->size == found->capacity) {
                found->capacity = found->capacity ? found->capacity * 2 : 16;
                found->data = (KeyPair*)realloc(found->data, found->capacity * sizeof(KeyPair));
            }
            found->data[found->size].start = si->starts[i];
            found->data[found->size].end = si->ends[i];
            found->size++;
            --i;
        } else {
            i = si->branch[i];
        }
    }
}

void searchItems(cSuperIntervals* si, int32_t start, int32_t end, cItemResult* found) {
    if (si->size == 0) {
        return;
    }
    size_t idx = upperBound(si, end);
    if (idx == SI_NONE) {
        return;
    }
    size_t i = idx;
    while (i != SI_NONE && start <= si->ends[i]) {
        if (found->size == found->capacity) {
            found->capacity = found->capacity ? found->capacity * 2 : 16;
            found->data = (Interval*)realloc(found->data, found->capacity * sizeof(Interval));
        }
        found->data[found->size].start = si->starts[i];
        found->data[found->size].end = si->ends[i];
        found->data[found->size].data = si->data[i];
        found->size++;
        --i;
    }
    if (i == SI_NONE) {
        return;
    }
    i = si->branch[i];
    while (i != SI_NONE) {
        if (start <= si->ends[i]) {
            if (found->size == found->capacity) {
                found->capacity = found->capacity ? found->capacity * 2 : 16;
                found->data = (Interval*)realloc(found->data, found->capacity * sizeof(Interval));
            }
            found->data[found->size].start = si->starts[i];
            found->data[found->size].end = si->ends[i];
            found->data[found->size].data = si->data[i];
            found->size++;
            --i;
        } else {
            i = si->branch[i];
        }
    }
}

void searchPoint(cSuperIntervals* si, int32_t point, cIndexResult* found) {
    searchValues(si, point, point, found);
}

size_t countOverlaps(cSuperIntervals* si, int32_t start, int32_t end) {
    if (si->size == 0) {
        return 0;
    }
    size_t idx = upperBound(si, end);
    if (idx == SI_NONE) {
        return 0;
    }
    size_t count = 0;
    size_t i = idx;
    while (i != SI_NONE && start <= si->ends[i]) {
        ++count;
        --i;
    }
    if (i == SI_NONE) {
        return count;
    }
    i = si->branch[i];
    while (i != SI_NONE) {
        if (start <= si->ends[i]) {
            ++count;
            --i;
        } else {
            i = si->branch[i];
        }
    }
    return count;
}

void coverage(cSuperIntervals* si, int32_t start, int32_t end,
              size_t* count_out, int32_t* coverage_out) {
    *count_out = 0;
    *coverage_out = 0;
    if (si->size == 0) {
        return;
    }
    size_t idx = upperBound(si, end);
    if (idx == SI_NONE) {
        return;
    }
    size_t i = idx;
    while (i != SI_NONE && start <= si->ends[i]) {
        int32_t hi = si->ends[i] < end ? si->ends[i] : end;
        int32_t lo = si->starts[i] > start ? si->starts[i] : start;
        ++(*count_out);
        *coverage_out += hi - lo;
        --i;
    }
    if (i == SI_NONE) {
        return;
    }
    i = si->branch[i];
    while (i != SI_NONE) {
        if (start <= si->ends[i]) {
            int32_t hi = si->ends[i] < end ? si->ends[i] : end;
            int32_t lo = si->starts[i] > start ? si->starts[i] : start;
            ++(*count_out);
            *coverage_out += hi - lo;
            --i;
        } else {
            i = si->branch[i];
        }
    }
}

void findOverlaps(cSuperIntervals* si, int32_t start, int32_t end,
                  int32_t* found, size_t* found_size) {
    *found_size = 0;
    if (si->size == 0) {
        return;
    }
    size_t idx = upperBound(si, end);
    if (idx == SI_NONE) {
        return;
    }
    size_t i = idx;
    while (i != SI_NONE && start <= si->ends[i]) {
        found[(*found_size)++] = si->data[i];
        --i;
    }
    if (i == SI_NONE) {
        return;
    }
    i = si->branch[i];
    while (i != SI_NONE) {
        if (start <= si->ends[i]) {
            found[(*found_size)++] = si->data[i];
            --i;
        } else {
            i = si->branch[i];
        }
    }
}

/* ----------------------------------------------------------------------------
 * Set operations
 * ------------------------------------------------------------------------- */

static int32_t si_keep_first(int32_t a, int32_t b) { (void)b; return a; }

// Returns indices into starts/ends/data in ascending (start, end-desc) order.
// Cheap when already sorted; otherwise sorts a copied index list.
static size_t* si_sorted_order(const cSuperIntervals* si) {
    size_t n = si->size;
    size_t* order = (size_t*)malloc(n * sizeof(size_t));
    // startSorted/endSorted are maintained incrementally by addInterval()/indexSuperIntervals().
    // When both hold the data is already in (start, end-desc) order, so skip the qsort entirely.
    if (si->startSorted && si->endSorted) {
        for (size_t k = 0; k < n; ++k) order[k] = k;
        return order;
    }
    Interval* tmp = (Interval*)malloc(n * sizeof(Interval));
    for (size_t k = 0; k < n; ++k) {
        tmp[k].start = si->starts[k];
        tmp[k].end = si->ends[k];
        tmp[k].data = (int32_t)k;  // remember original index in data slot
    }
    qsort(tmp, n, sizeof(Interval), compareIntervalsStart);
    for (size_t k = 0; k < n; ++k) {
        order[k] = (size_t)tmp[k].data;
    }
    free(tmp);
    return order;
}

cSuperIntervals* mergeOverlaps(const cSuperIntervals* si, cCombineFn combine) {
    if (!combine) combine = si_keep_first;
    cSuperIntervals* out = createSuperIntervals();
    if (si->size == 0) {
        return out;
    }
    size_t* order = si_sorted_order(si);
    int32_t cur_start = si->starts[order[0]];
    int32_t cur_end = si->ends[order[0]];
    int32_t cur_data = si->data[order[0]];
    for (size_t k = 1; k < si->size; ++k) {
        size_t idx = order[k];
        if (si->starts[idx] <= cur_end) {  // overlap -> extend
            if (si->ends[idx] > cur_end) {
                cur_end = si->ends[idx];
            }
            cur_data = combine(cur_data, si->data[idx]);
        } else {  // disjoint -> flush
            addInterval(out, cur_start, cur_end, cur_data);
            cur_start = si->starts[idx];
            cur_end = si->ends[idx];
            cur_data = si->data[idx];
        }
    }
    addInterval(out, cur_start, cur_end, cur_data);
    free(order);
    return out;
}

cSuperIntervals* intervalGaps(const cSuperIntervals* si, int32_t lo, int32_t hi, int32_t fill) {
    cSuperIntervals* out = createSuperIntervals();
    cSuperIntervals* merged = mergeOverlaps(si, NULL);  // start-sorted, disjoint
    int32_t cursor = lo;
    for (size_t k = 0; k < merged->size; ++k) {
        int32_t s = merged->starts[k];
        int32_t e = merged->ends[k];
        if (e < lo || s > hi) {
            continue;  // entirely outside the span
        }
        if (s > cursor) {
            addInterval(out, cursor, s - 1, fill);  // gap before this interval
        }
        if (e + 1 > cursor) {
            cursor = e + 1;
        }
    }
    if (cursor <= hi) {
        addInterval(out, cursor, hi, fill);  // trailing gap
    }
    destroySuperIntervals(merged);
    return out;
}

cSuperIntervals* unionWith(const cSuperIntervals* si, const cSuperIntervals* other, cCombineFn combine) {
    cSuperIntervals* combined = createSuperIntervals();
    reserveSuperIntervals(combined, si->size + other->size);
    for (size_t k = 0; k < si->size; ++k) {
        addInterval(combined, si->starts[k], si->ends[k], si->data[k]);
    }
    for (size_t k = 0; k < other->size; ++k) {
        addInterval(combined, other->starts[k], other->ends[k], other->data[k]);
    }
    cSuperIntervals* out = mergeOverlaps(combined, combine);
    destroySuperIntervals(combined);
    return out;
}

cSuperIntervals* intersection(const cSuperIntervals* si, cSuperIntervals* other, cCombineFn combine) {
    if (!combine) combine = si_keep_first;
    cSuperIntervals* out = createSuperIntervals();
    cIndexResult hits = createIndexResult();
    for (size_t k = 0; k < si->size; ++k) {
        clearIndexResult(&hits);
        searchIdxs(other, si->starts[k], si->ends[k], &hits);
        for (size_t h = 0; h < hits.size; ++h) {
            size_t oi = (size_t)hits.data[h];
            int32_t s = si->starts[k] > other->starts[oi] ? si->starts[k] : other->starts[oi];
            int32_t e = si->ends[k] < other->ends[oi] ? si->ends[k] : other->ends[oi];
            if (s <= e) {
                addInterval(out, s, e, combine(si->data[k], other->data[oi]));
            }
        }
    }
    destroyIndexResult(&hits);
    return out;
}

static int compareKeyPair(const void* a, const void* b) {
    const KeyPair* ka = (const KeyPair*)a;
    const KeyPair* kb = (const KeyPair*)b;
    if (ka->start != kb->start) {
        return (ka->start < kb->start) ? -1 : 1;
    }
    return (ka->end < kb->end) ? -1 : (ka->end > kb->end);
}

cSuperIntervals* difference(const cSuperIntervals* si, cSuperIntervals* other) {
    cSuperIntervals* out = createSuperIntervals();
    cKeyResult cover = createKeyResult();
    for (size_t k = 0; k < si->size; ++k) {
        clearKeyResult(&cover);
        searchKeys(other, si->starts[k], si->ends[k], &cover);
        qsort(cover.data, cover.size, sizeof(KeyPair), compareKeyPair);
        int32_t cursor = si->starts[k];
        for (size_t c = 0; c < cover.size; ++c) {
            int32_t cs = cover.data[c].start > si->starts[k] ? cover.data[c].start : si->starts[k];
            int32_t ce = cover.data[c].end < si->ends[k] ? cover.data[c].end : si->ends[k];
            if (cs > cursor) {
                addInterval(out, cursor, cs - 1, si->data[k]);  // uncovered piece
            }
            if (ce + 1 > cursor) {
                cursor = ce + 1;
            }
        }
        if (cursor <= si->ends[k]) {
            addInterval(out, cursor, si->ends[k], si->data[k]);  // trailing uncovered piece
        }
    }
    destroyKeyResult(&cover);
    return out;
}

cSuperIntervals* symmetricDifference(cSuperIntervals* si, cSuperIntervals* other) {
    cSuperIntervals* a_minus_b = difference(si, other);
    cSuperIntervals* b_minus_a = difference(other, si);
    cSuperIntervals* out = unionWith(a_minus_b, b_minus_a, NULL);
    destroySuperIntervals(a_minus_b);
    destroySuperIntervals(b_minus_a);
    return out;
}

bool intervalSpan(const cSuperIntervals* si, int32_t* lo_out, int32_t* hi_out) {
    if (si->size == 0) {
        return false;
    }
    int32_t lo = si->starts[0];
    int32_t hi = si->ends[0];
    for (size_t k = 1; k < si->size; ++k) {
        if (si->starts[k] < lo) lo = si->starts[k];
        if (si->ends[k] > hi) hi = si->ends[k];
    }
    *lo_out = lo;
    *hi_out = hi;
    return true;
}

cSuperIntervals* expandIntervals(const cSuperIntervals* si, int32_t left, int32_t right,
                                 int32_t lo, int32_t hi) {
    cSuperIntervals* out = createSuperIntervals();
    reserveSuperIntervals(out, si->size);
    for (size_t k = 0; k < si->size; ++k) {
        // Compute in 64-bit then clamp to [lo, hi]; avoids overflow when lo/hi are
        // the type extremes. left/right may be negative (shrink).
        int64_t s64 = (int64_t)si->starts[k] - left;
        int64_t e64 = (int64_t)si->ends[k] + right;
        if (s64 < lo) s64 = lo;
        if (e64 > hi) e64 = hi;
        if (s64 <= e64) {
            addInterval(out, (int32_t)s64, (int32_t)e64, si->data[k]);
        }
    }
    return out;
}

cSuperIntervals* flankIntervals(const cSuperIntervals* si, int32_t left, int32_t right,
                                int32_t lo, int32_t hi) {
    cSuperIntervals* out = createSuperIntervals();
    for (size_t k = 0; k < si->size; ++k) {
        if (left > 0 && si->starts[k] > lo) {
            int32_t le = si->starts[k] - 1;
            int64_t ls64 = (int64_t)si->starts[k] - left;  // 64-bit avoids underflow
            if (ls64 < lo) ls64 = lo;
            if (ls64 <= le && le <= hi) {
                addInterval(out, (int32_t)ls64, le, si->data[k]);
            }
        }
        if (right > 0 && si->ends[k] < hi) {
            int32_t rs = si->ends[k] + 1;
            int64_t re64 = (int64_t)si->ends[k] + right;  // 64-bit avoids overflow
            if (re64 > hi) re64 = hi;
            if (rs <= re64 && rs >= lo) {
                addInterval(out, rs, (int32_t)re64, si->data[k]);
            }
        }
    }
    return out;
}

cSuperIntervals* uniqueIntervals(const cSuperIntervals* si, cCombineFn combine) {
    if (!combine) combine = si_keep_first;
    cSuperIntervals* out = createSuperIntervals();
    if (si->size == 0) {
        return out;
    }
    size_t* order = si_sorted_order(si);  // (start, end) grouped
    size_t run = order[0];
    int32_t acc = si->data[order[0]];
    for (size_t k = 1; k < si->size; ++k) {
        size_t idx = order[k];
        if (si->starts[idx] == si->starts[run] && si->ends[idx] == si->ends[run]) {
            acc = combine(acc, si->data[idx]);  // same coordinates -> fold
        } else {
            addInterval(out, si->starts[run], si->ends[run], acc);
            run = idx;
            acc = si->data[idx];
        }
    }
    addInterval(out, si->starts[run], si->ends[run], acc);
    free(order);
    return out;
}

/* ----------------------------------------------------------------------------
 * Result buffer helpers
 * ------------------------------------------------------------------------- */

cIndexResult createIndexResult(void) {
    cIndexResult r = { NULL, 0, 0 };
    return r;
}
void clearIndexResult(cIndexResult* r) { r->size = 0; }
void destroyIndexResult(cIndexResult* r) { free(r->data); r->data = NULL; r->size = 0; r->capacity = 0; }

cKeyResult createKeyResult(void) {
    cKeyResult r = { NULL, 0, 0 };
    return r;
}
void clearKeyResult(cKeyResult* r) { r->size = 0; }
void destroyKeyResult(cKeyResult* r) { free(r->data); r->data = NULL; r->size = 0; r->capacity = 0; }

cItemResult createItemResult(void) {
    cItemResult r = { NULL, 0, 0 };
    return r;
}
void clearItemResult(cItemResult* r) { r->size = 0; }
void destroyItemResult(cItemResult* r) { free(r->data); r->data = NULL; r->size = 0; r->capacity = 0; }

#endif  // SUPERINTERVALS_HEADER_INCLUDED

#ifdef __cplusplus
}  // extern "C"
#endif
