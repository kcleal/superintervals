#ifdef __cplusplus
extern "C" {
#endif

#ifndef SUPERINTERVALS_HEADER_INCLUDED
#define SUPERINTERVALS_HEADER_INCLUDED 1
#endif
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>

#ifdef __AVX2__
    #include <immintrin.h>
#elif defined __ARM_NEON
    #include <arm_neon.h>
#endif

// Define the Interval struct
typedef struct {
    int32_t start;
    int32_t end;
    int32_t data;
} Interval;

// Define the SuperIntervals struct
typedef struct {
    int32_t* starts;
    int32_t* ends;
    int32_t* data;
    size_t* branch;
    int32_t* extent;
    size_t size;
    size_t capacity;
    size_t idx;
    bool startSorted;
    bool endSorted;
} cSuperIntervals;

// Function prototypes
cSuperIntervals* createSuperIntervals();
void destroySuperIntervals(cSuperIntervals* si);
void clearSuperIntervals(cSuperIntervals* si);
void reserveSuperIntervals(cSuperIntervals* si, size_t n);
void addInterval(cSuperIntervals* si, int32_t start, int32_t end, int32_t value);
void sortIntervals(cSuperIntervals* si);
void indexSuperIntervals(cSuperIntervals* si, bool use_linear);
void upperBound(cSuperIntervals* si, int32_t value);
bool anyOverlaps(cSuperIntervals* si, int32_t start, int32_t end);
void findOverlaps(cSuperIntervals* si, int32_t start, int32_t end, int32_t* found, size_t* found_size);
size_t countOverlaps(cSuperIntervals* si, int32_t start, int32_t end);

// Helper function prototypes
static void sortBlock(cSuperIntervals* si, size_t start_i, size_t end_i, int (*compare)(const void*, const void*));
static size_t eytzinger_helper(int32_t* arr, size_t n, size_t i, size_t k, int32_t* eytz, size_t* eytz_index);
static int eytzinger(int32_t* arr, size_t n, int32_t* eytz, size_t* eytz_index);

cSuperIntervals* createSuperIntervals() {
    cSuperIntervals* si = (cSuperIntervals*)malloc(sizeof(cSuperIntervals));
    si->starts = NULL;
    si->ends = NULL;
    si->data = NULL;
    si->branch = NULL;
    si->extent = NULL;
    si->size = 0;
    si->capacity = 0;
    si->idx = 0;
    si->startSorted = true;
    si->endSorted = true;
    return si;
}

void destroySuperIntervals(cSuperIntervals* si) {
    free(si->starts);
    free(si->ends);
    free(si->data);
    free(si->branch);
    free(si->extent);
    free(si);
}

void clearSuperIntervals(cSuperIntervals* si) {
    si->size = 0;
    si->idx = 0;
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

static int compareIntervalsStart(const void* a, const void* b) {
    const Interval* ia = (const Interval*)a;
    const Interval* ib = (const Interval*)b;
    if (ia->start != ib->start) {
        return ia->start - ib->start;
    }
    return ib->end - ia->end;  // Sort start, and end in descending order
}

static int compareIntervalsEnd(const void* a, const void* b) {
    const Interval* ia = (const Interval*)a;
    const Interval* ib = (const Interval*)b;
    return ib->end - ia->end;  // Sort by end in descending order
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

static size_t eytzinger_helper(int32_t* arr, size_t n, size_t i, size_t k, int32_t* eytz, size_t* eytz_index) {
    if (k < n) {
        i = eytzinger_helper(arr, n, i, 2*k+1, eytz, eytz_index);
        eytz[k] = arr[i];
        eytz_index[k] = i;
        ++i;
        i = eytzinger_helper(arr, n, i, 2*k + 2, eytz, eytz_index);
    }
    return i;
}

static int eytzinger(int32_t* arr, size_t n, int32_t* eytz, size_t* eytz_index) {
    return eytzinger_helper(arr, n, 0, 0, eytz, eytz_index);
}

void indexSuperIntervals(cSuperIntervals* si, bool use_linear) {
    if (si->size == 0) {
        return;
    }

    sortIntervals(si);

    int32_t* eytz = (int32_t*)malloc((si->size + 1) * sizeof(int32_t));
    size_t* eytz_index = (size_t*)malloc((si->size + 1) * sizeof(size_t));
    eytzinger(si->starts, si->size, eytz, eytz_index);

    si->branch = (size_t*)realloc(si->branch, si->size * sizeof(size_t));
    memset(si->branch, -1, si->size * sizeof(size_t));

    if (!use_linear) {
        si->extent = (int32_t*)malloc(si->size * sizeof(int32_t));
        memcpy(si->extent, si->ends, si->size * sizeof(int32_t));

        for (size_t i = 0; i < si->size - 1; ++i) {
            int32_t e = si->ends[i];
            for (size_t j = i + 1; j < si->size; ++j) {
                if (si->ends[j] >= si->ends[i]) {
                    break;
                }
                si->branch[j] = i;
                if (e > si->extent[j]) {
                    si->extent[j] = e;
                }
            }
        }
    } else {
        // Linear branch implementation
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
    }

    free(eytz);
    free(eytz_index);
    si->idx = 0;
}


void upperBound(cSuperIntervals* si, int32_t value) {
    size_t length = si->size - 1;
    si->idx = 0;
    const int entries_per_256KB = 256 * 1024 / sizeof(int32_t);
    const int num_per_cache_line = 64 / sizeof(int32_t) > 1 ? 64 / sizeof(int32_t) : 1;

    if (length >= entries_per_256KB) {
        while (length >= 3 * num_per_cache_line) {
            size_t half = length / 2;
//            __builtin_prefetch(&si->starts[si->idx + half / 2]);
            size_t first_half1 = si->idx + (length - half);
//            __builtin_prefetch(&si->starts[first_half1 + half / 2]);
            si->idx += (si->starts[si->idx + half] <= value) * (length - half);
            length = half;
        }
    }

    while (length > 0) {
        size_t half = length / 2;
        si->idx += (si->starts[si->idx + half] <= value) * (length - half);
        length = half;
    }

    if (si->idx > 0 && (si->idx == si->size - 1 || si->starts[si->idx] > value)) {
        --si->idx;
    }
}

bool anyOverlaps(cSuperIntervals* si, int32_t start, int32_t end) {
    upperBound(si, end);
    return start <= si->ends[si->idx];
}

void findOverlaps(cSuperIntervals* si, int32_t start, int32_t end, int32_t* found, size_t* found_size) {
    if (si->size == 0) {
        *found_size = 0;
        return;
    }
    upperBound(si, end);
    size_t i = si->idx;
    *found_size = 0;
    while (i > 0) {
        if (start <= si->ends[i]) {
            found[(*found_size)++] = si->data[i];
            --i;

        } else {
            if (si->branch[i] >= i) {  // segfaults here
                break;
            }
            i = si->branch[i];
        }
    }
    if (i == 0 && start <= si->ends[0] && si->starts[0] <= end) {
        found[(*found_size)++] = si->data[0];
    }
}


#ifdef __cplusplus
}  // extern C
#endif