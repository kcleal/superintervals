/*
SuperIntervals - a static data structure for finding interval intersections

Notes
-----
 - intervals are considered end-inclusive
 - the index() function must be called before any queries. If more
   intervals are added, call index() again.

*/
#pragma once

#include <algorithm>
#include <vector>
#include <iostream>
#ifdef __AVX2__
    #include <immintrin.h>
#elif defined __ARM_NEON
    #include <arm_neon.h>
#endif

// S for scalar for start, end. T for data type
template<typename S, typename T>
class SuperIntervals {
    public:
    struct Interval {
        S start, end;
        size_t index;
    };
    alignas(alignof(std::vector<S>)) std::vector<S> starts;
    alignas(alignof(std::vector<S>)) std::vector<S> ends;
    alignas(alignof(size_t)) std::vector<size_t> branch;
    std::vector<Interval> intervals;
    std::vector<T> data;
    size_t idx, n_intervals;
    bool startSorted, endSorted;
    S it_low, it_high;
    SuperIntervals()
        : idx(0)
        , n_intervals(0)
        , startSorted(true)
        , endSorted(true)
        , it_low(0)
        , it_high(0)
        {}

    ~SuperIntervals() = default;

    void clear() {
        intervals.clear(); data.clear(); starts.clear(); ends.clear(); branch.clear(); idx = 0;
    }

    void reserve(size_t n) {
        intervals.reserve(n); data.reserve(n); starts.reserve(n); ends.reserve(n);
    }

    size_t size() {
        return intervals.size();
    }

    void add(S start, S end, const T& value) {
        if (startSorted && !intervals.empty()) {
            startSorted = (start < intervals.back().start) ? false : true;
            if (startSorted && start == intervals.back().start && end > intervals.back().end) {
                endSorted = false;
            }
        }
        intervals.emplace_back() = {start, end, n_intervals};
        data.emplace_back(value);
        ++n_intervals;
    }

    void sortIntervals() {
        size_t i;
        if (intervals.size() <= 1) {
            startSorted = true;
            endSorted = true;
            return;
        }
        if (!startSorted) {
            std::sort(intervals.begin(), intervals.end(),
            [](const Interval& a, const Interval& b) { return (a.start < b.start || (a.start == b.start && a.end > b.end)); });
            std::vector<T> data_copy = data;
            i = 0;
            for (const auto& itv : intervals) {
                starts[i] = itv.start;
                ends[i] = itv.end;
                data[i++] = data_copy[itv.index];
            }
            startSorted = true;
            endSorted = true;
        } else if (!endSorted) {  // only sort parts that need sorting - ends in descending order
            std::vector<T> data_copy = data;
            auto it_start = intervals.begin();
            while (it_start != intervals.end()) {
                auto block_end = it_start + 1;
                bool srt = false;
                while (block_end != intervals.end() && block_end->start == it_start->start) {
                    if (block_end->end > (block_end - 1)->end) {
                        srt = true;
                    }
                    ++block_end;
                }
                if (srt) {
                    std::sort(it_start, block_end, [](const Interval& a, const Interval& b) { return a.end > b.end; });
                    i = std::distance(intervals.begin(), it_start);
                    while (it_start != intervals.end()) {
                        data[i++] = data_copy[it_start->index];
                        if (it_start == block_end) {
                            break;
                        }
                        ++it_start;
                    }
                }
                it_start = block_end;
            }
            i = 0;
            for (const auto& itv : intervals) {
                starts[i] = itv.start;
                ends[i++] = itv.end;
            }
            endSorted = true;
        } else {
            i = 0;
            for (const auto& itv : intervals) {
                starts[i] = itv.start;
                ends[i++] = itv.end;
            }
        }
    }

    void index() {
        n_intervals = intervals.size();
        if (intervals.empty()) {
            return;
        }
        starts.resize(n_intervals);
        ends.resize(n_intervals);

        sortIntervals();

//        std::vector<int> counts(intervals.size(), 0);
//        int max = 0;
//        branch.resize(intervals.size(), SIZE_MAX);
//        for (size_t i=0; i < ends.size() - 1; ++i) {
//            for (size_t j=i + 1; j < ends.size(); ++j) {
//                if (ends[j] >= ends[i]) {
//                    break;
//                }
//                branch[j] = i;
////                counts[j] += 1;
////                if (counts[j] > max) {
////                    max = counts[j];
////                }
//            }
//        }

        // Linear construction, but slower for most cases!
        branch.resize(intervals.size(), SIZE_MAX);
        std::vector<std::pair<S, size_t>> br;
        br.reserve(1000);
        br.emplace_back() = {ends[0], 0};
        for (size_t i=1; i < ends.size(); ++i) {
            while (!br.empty() && br.back().first < ends[i]) {
                br.pop_back();
            }
            if (!br.empty()) {
                branch[i] = br.back().second;
            }
            br.emplace_back() = {ends[i], i};
        }
        idx = 0;
    }

    struct IntervalItem {
        S start, end;
        T data;
    };

    IntervalItem at(size_t index) {
        return IntervalItem(starts[index], ends[index], data[index]);
    }

    void at(size_t index, IntervalItem& itv) {
        itv.start = starts[index];
        itv.end = ends[index];
        itv.data = data[index];
    }

    class Iterator {
    public:
        size_t it_index;

        Iterator(const SuperIntervals* list, size_t index) : super(list) {
            _start = list->it_low;
            _end = list->it_high;
            it_index = index;
        }
        typename SuperIntervals::IntervalItem operator*() const {
            return typename SuperIntervals<S, T>::IntervalItem{super->starts[it_index], super->ends[it_index], super->data[it_index]};
        }
        Iterator& operator++() {
            if (it_index == 0) {
                it_index = SIZE_MAX;
                return *this;
            }
            if (it_index > 0) {
                if (_start <= super->ends[it_index]) {
                    --it_index;
//                    if (_start > super->ends[it_index]) {
//                        it_index = SIZE_MAX;
//                    }
                } else {
                    if (super->branch[it_index] >= it_index) {
                        it_index = SIZE_MAX;
                        return *this;
                    }
                    it_index = super->branch[it_index];
                    if (_start <= super->ends[it_index]) {
                        --it_index;
//                        if (_start > super->ends[it_index]) {
//                            it_index = SIZE_MAX;
//                        }
                    } else {
                        it_index = SIZE_MAX;
                        return *this;
                    }
                }
            }
            return *this;
        }
        bool operator!=(const Iterator& other) const {
            return it_index != other.it_index;
        }
        bool operator==(const Iterator& other) const {
            return it_index == other.it_index;
        }
        Iterator begin() const { return Iterator(super, super->idx); }
        Iterator end() const { return Iterator(super, SIZE_MAX); }
    private:
        S _start, _end;
        const SuperIntervals<S, T>* super;

    };

    Iterator begin() const { return Iterator(this, idx); }
    Iterator end() const { return Iterator(this, 0); }

    void searchInterval(S start, S end) {
        it_low = start; it_high = end;
        upperBound(end);
        if (start > ends[idx] || starts[0] > end) {
            idx = SIZE_MAX;
        }
    }

    inline void upperBound(S value) {  // https://github.com/mh-dm/sb_lower_bound/blob/master/sbpm_lower_bound.h
        size_t length = n_intervals;
        idx = 0;
        constexpr int entries_per_256KB = 256 * 1024 / sizeof(S);
        if (length >= entries_per_256KB) {
            constexpr int num_per_cache_line = std::max(64 / int(sizeof(S)), 1);
            while (length >= 3 * num_per_cache_line) {
                size_t half = length / 2;
                __builtin_prefetch(&starts[idx + half / 2]);
                size_t first_half1 = idx + (length - half);
                __builtin_prefetch(&starts[first_half1 + half / 2]);
                idx += (starts[idx + half] <= value) * (length - half);
                length = half;
            }
        }
        while (length > 0) {
            size_t half = length / 2;
            idx += (starts[idx + half] <= value) * (length - half);
            length = half;
        }
        if (idx > 0 && (idx == n_intervals || starts[idx] > value)) {
            --idx;
        }
    }

    bool anyOverlaps(S start, S end) {
        upperBound(end);
        return start <= ends[idx];
    }

    void findOverlaps(S start, S end, std::vector<T>& found) {
        if (!n_intervals) {
            return;
        }
        if (!found.empty()) {
            found.clear();
        }
        upperBound(end);
        size_t i = idx;
        while (i > 0) {
            if (start <= ends[i--]) {
                found.push_back(i+1);
            } else {
                if (++i; branch[i] >= i) {
                    break;
                }
                i = branch[i];
                __builtin_prefetch(&ends[i]);
            }
        }
        if (i==0 && start <= ends[0] && starts[0] <= end) {
            found.push_back(0);
        }
    }

    size_t countOverlaps(S start, S end) {
        if (!n_intervals) {
            return 0;
        }
        upperBound(end);
        size_t found = 0;
        size_t i = idx;

#ifdef __AVX2__
        __m256i start_vec = _mm256_set1_epi32(start);
        constexpr size_t simd_width = 256 / (sizeof(S) * 8);
        constexpr size_t block = simd_width * 4;
#elif defined __ARM_NEON
        int32x4_t start_vec = vdupq_n_s32(start);
        constexpr size_t simd_width = 128 / (sizeof(S) * 8);
        uint32x4_t ones = vdupq_n_u32(1);
        constexpr size_t block = simd_width * 4;
#else
        constexpr size_t block = 16;
#endif

        while (i > 0) {
            if (start <= ends[i]) {
                ++found;
                --i;

#ifdef __AVX2__
                while (i > block) {
                    size_t count = 0;
                    for (size_t j = i; j > i - block; j -= simd_width) {
                        __m256i ends_vec = _mm256_load_si256((__m256i*)(&ends[j - simd_width + 1]));
                        __m256i cmp_mask = _mm256_cmpgt_epi32(start_vec, ends_vec);
                        int mask = _mm256_movemask_epi8(~cmp_mask);
                        count += _mm_popcnt_u32(mask) / 4;  // Each comparison result is 4 bits
                    }
                    found += count;
                    i -= block;
                    if (count < block) {  // check for a branch
                        break;
                    }
                }
#elif defined __ARM_NEON
                while (i > block) {
                    size_t count = 0;
                    for (size_t j = i; j > i - block; j -= simd_width) { // Neon processes 4 int32 at a time
                        int32x4_t ends_vec = vld1q_s32(&ends[j - simd_width + 1]);
                        uint32x4_t mask = vcleq_s32(start_vec, ends_vec);
                        uint32x4_t bool_mask = vandq_u32(mask, ones); // Convert -1 to 1 for true elements
                        count += vaddvq_u32(bool_mask);
                    }
                    found += count;
                    i -= block;
                    if (count < block) {  // check for a branch
                        break;
                    }
                }
#else
                while (i > block) {
                    size_t count = 0;
                    for (size_t j = i; j > i - block; --j) {
                        count += (start <= ends[j]) ? 1 : 0;
                    }
                    found += count;
                    i -= block;
                    if (count < block) {  // check for a branch
                        break;
                    }
                }
#endif
            } else {
                if (branch[i] >= i) {
                    break;
                }
                i = branch[i];
            }
        }
        if (i==0 && start <= ends[0] && starts[0] <= end) {
            ++found;
        }
        return found;
    }
};
