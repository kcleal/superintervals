
#pragma once

#include <algorithm>
#include <vector>
#include <climits>
#include <iostream>
#include <limits>
#if defined(__AVX2__)
    #include <immintrin.h>
#elif defined(__ARM_NEON__)
    #include <arm_neon.h>
#endif

/**
 * @file SuperIntervals.hpp
 * @brief A static data structure for finding interval intersections
 *
 * SuperIntervals is a template class that provides efficient interval intersection operations.
 * It supports adding intervals, indexing them for fast queries, and performing various
 * intersection operations.
 *
 * @note Intervals are considered end-inclusive
 * @note The index() function must be called before any queries. If more intervals are added, call index() again.
 *
 * @tparam S The scalar type for interval start and end points (e.g., int, float)
 * @tparam T The data type associated with each interval
 */
template<typename S, typename T>
class SuperIntervals {
    public:

    struct Interval {
        S start, end;
        T data;
    };

    alignas(alignof(std::vector<S>)) std::vector<S> starts;
    alignas(alignof(std::vector<S>)) std::vector<S> ends;
    alignas(alignof(size_t)) std::vector<size_t> branch;
    std::vector<T> data;
    size_t idx;
    bool startSorted, endSorted;

    SuperIntervals()
        : idx(0)
        , startSorted(true)
        , endSorted(true)
        , it_low(0)
        , it_high(0)
        {}

    ~SuperIntervals() = default;

    /**
     * @brief Clears all intervals and resets the data structure
     */
    void clear() {
        data.clear(); starts.clear(); ends.clear(); branch.clear(); idx = 0;
    }

    /**
     * @brief Reserves memory for a specified number of intervals
     * @param n Number of intervals to reserve space for
     */
    void reserve(size_t n) {
        data.reserve(n); starts.reserve(n); ends.reserve(n);
    }

    /**
    * @brief Returns the number of intervals in the data structure
    * @return Number of intervals
    */
    size_t size() {
        return starts.size();
    }

    /**
     * @brief Adds a new interval to the data structure
     * @param start Start point of the interval
     * @param end End point of the interval
     * @param value Data associated with the interval
     */
    void add(S start, S end, const T& value) {
        if (startSorted && !starts.empty()) {
            startSorted = (start < starts.back()) ? false : true;
            if (startSorted && start == starts.back() && end > ends.back()) {
                endSorted = false;
            }
        }
        starts.push_back(start);
        ends.push_back(end);
        data.emplace_back(value);
    }

    /**
     * @brief Indexes the intervals.
     *
     * This function must be called after adding intervals and before performing any queries.
     * If more intervals are added after indexing, this function should be called again.
     */
    virtual void index() {
        if (starts.size() == 0) {
            return;
        }
        starts.shrink_to_fit();
        ends.shrink_to_fit();
        data.shrink_to_fit();
        sortIntervals();
        branch.resize(starts.size(), SIZE_MAX);
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

    /**
     * @brief Retrieves an interval at a specific index
     * @param index The index of the interval to retrieve
     * @return The Interval at the specified index
     */
    Interval at(size_t index) {
        return Interval(starts[index], ends[index], data[index]);
    }

    void at(size_t index, Interval& itv) {
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
        typename SuperIntervals::Interval operator*() const {
            return typename SuperIntervals<S, T>::Interval{super->starts[it_index], super->ends[it_index], super->data[it_index]};
        }
        Iterator& operator++() {
            if (it_index == 0) {
                it_index = SIZE_MAX;
                return *this;
            }
            if (it_index > 0) {
                if (_start <= super->ends[it_index]) {
                    --it_index;
                } else {
                    if (super->branch[it_index] >= it_index) {
                        it_index = SIZE_MAX;
                        return *this;
                    }
                    it_index = super->branch[it_index];
                    if (_start <= super->ends[it_index]) {
                        --it_index;
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
    Iterator end() const { return Iterator(this, SIZE_MAX); }

    /**
     * @brief Sets the search interval. Must be called before using the iterator.
     * @param start Start point of the search range
     * @param end End point of the search range
     */
    void searchInterval(const S start, const S end) noexcept {
        if (starts.empty()) {
            return;
        }
        it_low = start; it_high = end;
        upperBound(end);
        if (start > ends[idx] || starts[0] > end) {
            idx = SIZE_MAX;
        }
    }

    virtual inline void upperBound(const S value) noexcept {  // https://github.com/mh-dm/sb_lower_bound/blob/master/sbpm_lower_bound.h
        size_t length = starts.size() - 1;
        idx = 0;
//        constexpr int num_per_cache_line = 3 * hardware_constructive_interference_size;
//        while (length >= num_per_cache_line) {
//            size_t half = length / 2;
////                __builtin_prefetch(&starts[idx + half / 2]);
////            size_t first_half1 = idx + (length - half);
////                __builtin_prefetch(&starts[first_half1 + half / 2]);
//            idx += (starts[idx + half] <= value) * (length - half);
//            length = half;
//        }

        while (length > 0) {
            size_t half = length / 2;
            idx += (starts[idx + half] <= value) * (length - half);
            length = half;
        }
        if (idx > 0 && (idx == length || starts[idx] > value)) {
            --idx;
        }
    }

    bool anyOverlaps(const S start, const S end) noexcept {
        if (starts.empty()) {
            return false;
        }
        upperBound(end);
        size_t i = idx;
        while (i > 0) {
            if (start <= ends[i]) {
                return true;
                --i;
            } else {
                if (branch[i] >= i) {
                    break;
                }
                i = branch[i];
            }
        }
        if (i==0 && start <= ends[0] && starts[0] <= end) {
            return true;
        }
        return false;
    }

    void findOverlaps(const S start, const S end, std::vector<T>& found) {
        if (starts.empty()) {
            return;
        }
        upperBound(end);
        size_t i = idx;
        while (i > 0) {
            if (start <= ends[i]) {
                found.push_back(data[i]);
                --i;
            } else {
                if (branch[i] >= i) {
                    break;
                }
                i = branch[i];
            }
        }
        if (i==0 && start <= ends[0] && starts[0] <= end) {
            found.push_back(data[0]);
        }
    }

    size_t countOverlaps(const S start, const S end) noexcept {
        if (starts.empty()) {
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
                        count += _mm_popcnt_u32(mask);
                    }
                    found += count / 4;  // Each comparison result is 4 bits
                    i -= block;
                    if (count < block) {
                        break;
                    }
                }
#elif defined(__ARM_NEON__)
                while (i > block) {
                    size_t count = 0;
                    uint32x4_t mask, bool_mask;
                    for (size_t j = i; j > i - block; j -= simd_width) { // Neon processes 4 int32 at a time
                        int32x4_t ends_vec = vld1q_s32(&ends[j - simd_width + 1]);
                        mask = vcleq_s32(start_vec, ends_vec);  // True (0xFFFFFFFF) for elements where start_vec <= ends_vec
                        bool_mask = vandq_u32(mask, ones);
                        count += vaddvq_u32(bool_mask);
                    }
                    found += count;
                    i -= block;
//                    if (count < block && vgetq_lane_u32(mask, 0) == 0) {  // check for overlap again, before checking for branch?
                    if (count < block) {  // check for overlap again, before checking for branch?
                        break;
                    }
                }
#else  // Rely on compiler auto vectorize
                while (i > block) {
                    size_t count = 0;
                    for (size_t j = i; j > i - block; --j) {
                        count += (start <= ends[j]) ? 1 : 0;
                    }
                    found += count;
                    i -= block;
                    if (count < block && start > ends[i + 1]) {  // check for a branch
//                    if (count < block) {  // check for a branch
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

    void coverage(const S start, const S end, std::pair<size_t, S> &cov_result) {
        if (starts.empty()) {
            cov_result.first = 0;
            cov_result.second = 0;
            return;
        }
        upperBound(end);
        size_t i = idx;
        size_t cnt = 0;
        S cov = 0;
        while (i > 0) {
            if (start <= ends[i]) {
                ++cnt;
                cov += std::min(ends[i], end) - std::max(starts[i], start);
                --i;
            } else {
                if (branch[i] >= i) {
                    break;
                }
                i = branch[i];
            }
        }
        if (i==0 && start <= ends[0] && starts[0] <= end) {
            cov += std::min(ends[i], end) - std::max(starts[i], start);
            ++cnt;
        }
        cov_result.first = cnt;
        cov_result.second = cov;
    }

    void findStabbed(const S point, std::vector<T>& found) {
        if (starts.empty()) {
            return;
        }
        upperBound(point);
        size_t i = idx;
        while (i > 0) {
            if (point <= ends[i]) {
                found.push_back(data[i]);
                --i;
            } else {
                if (branch[i] >= i) {
                    break;
                }
                i = branch[i];
            }
        }
        if (i==0 && point <= ends[0] && starts[0] <= point) {
            found.push_back(data[0]);
        }
    }

    size_t countStabbed(const S point) noexcept {
        if (starts.empty()) {
            return 0;
        }
        upperBound(point);
        size_t found = 0;
        size_t i = idx;

#ifdef __AVX2__
        __m256i start_vec = _mm256_set1_epi32(point);
        constexpr size_t simd_width = 256 / (sizeof(S) * 8);
        constexpr size_t block = simd_width * 4;
#elif defined __ARM_NEON
        int32x4_t start_vec = vdupq_n_s32(point);
        constexpr size_t simd_width = 128 / (sizeof(S) * 8);
        uint32x4_t ones = vdupq_n_u32(1);
        constexpr size_t block = simd_width * 4;
#else
        constexpr size_t block = 16;
#endif

        while (i > 0) {
            if (point <= ends[i]) {
                ++found;
                --i;

#ifdef __AVX2__
                while (i > block) {
                    size_t count = 0;
                    for (size_t j = i; j > i - block; j -= simd_width) {
                        __m256i ends_vec = _mm256_load_si256((__m256i*)(&ends[j - simd_width + 1]));
                        __m256i cmp_mask = _mm256_cmpgt_epi32(start_vec, ends_vec);
                        int mask = _mm256_movemask_epi8(~cmp_mask);
                        count += _mm_popcnt_u32(mask);
                    }
                    found += count / 4;  // Each comparison result is 4 bits
                    i -= block;
                    if (count < block) {
                        break;
                    }
                }
#elif defined __ARM_NEON
                while (i > block) {
                    size_t count = 0;
                    uint32x4_t bool_mask;
                    for (size_t j = i; j > i - block; j -= simd_width) { // Neon processes 4 int32 at a time
                        int32x4_t ends_vec = vld1q_s32(&ends[j - simd_width + 1]);
                        uint32x4_t mask = vcleq_s32(start_vec, ends_vec);
                        bool_mask = vandq_u32(mask, ones); // Convert -1 to 1 for true elements
                        count += vaddvq_u32(bool_mask);
                    }
                    found += count;
                    i -= block;
                    if (count < block) {  // go check for a branch
                        break;
                    }
                }
#else
                while (i > block) {
                    size_t count = 0;
                    for (size_t j = i; j > i - block; --j) {
                        count += (point <= ends[j]) ? 1 : 0;
                    }
                    found += count;
                    i -= block;
                    if (count < block) {  // go check for a branch
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
        if (i==0 && point <= ends[0] && starts[0] <= point) {
            ++found;
        }
        return found;
    }

    protected:

    S it_low, it_high;
    std::vector<Interval> tmp;

    template<typename CompareFunc>
    void sortBlock(size_t start_i, size_t end_i, CompareFunc compare) {
        size_t range_size = end_i - start_i;
        tmp.resize(range_size);
        for (size_t i = 0; i < range_size; ++i) {
            tmp[i].start = starts[start_i + i];
            tmp[i].end = ends[start_i + i];
            tmp[i].data = data[start_i + i];
        }
        std::sort(tmp.begin(), tmp.end(), compare);
        for (size_t i = 0; i < range_size; ++i) {
            starts[start_i + i] = tmp[i].start;
            ends[start_i + i] = tmp[i].end;
            data[start_i + i] = tmp[i].data;
        }
    }

    void sortIntervals() {
        if (!startSorted) {
            sortBlock(0, starts.size(),
                [](const Interval& a, const Interval& b) { return (a.start < b.start || (a.start == b.start && a.end > b.end)); });
            startSorted = true;
            endSorted = true;
        } else if (!endSorted) {  // only sort parts that need sorting - ends in descending order
            size_t it_start = 0;
            while (it_start < starts.size()) {
                size_t block_end = it_start + 1;
                bool needs_sort = false;
                while (block_end < starts.size() && starts[block_end] == starts[it_start]) {
                    if (block_end > it_start && ends[block_end] > ends[block_end - 1]) {
                        needs_sort = true;
                    }
                    ++block_end;
                }
                if (needs_sort) {
                    sortBlock(it_start, block_end, [](const Interval& a, const Interval& b) { return a.end > b.end; });
                }
                it_start = block_end;
            }
            endSorted = true;
        }
    }

#ifdef __cpp_lib_hardware_interference_size
    static constexpr std::size_t hardware_constructive_interference_size = std::hardware_constructive_interference_size;  // Corresponds to cache line size
#else
    static constexpr std::size_t hardware_constructive_interference_size = 64;
#endif

};


template<typename S, typename T>
class SuperIntervalsEytz : public SuperIntervals<S, T> {
public:

//    alignas(alignof(std::vector<S>)) std::vector<S> extent;

    void index() override {
        if (this->starts.size() == 0) {
            return;
        }
        this->starts.shrink_to_fit();
        this->ends.shrink_to_fit();
        this->data.shrink_to_fit();
        this->sortIntervals();

        eytz.resize(this->starts.size() + 1);
        eytz_index.resize(this->starts.size() + 1);
        eytzinger(&this->starts[0], this->starts.size());

        this->branch.resize(this->starts.size(), SIZE_MAX);
        std::vector<std::pair<S, size_t>> br;
        br.reserve(1000);
        br.emplace_back() = {this->ends[0], 0};
        for (size_t i=1; i < this->ends.size(); ++i) {
            while (!br.empty() && br.back().first < this->ends[i]) {
                br.pop_back();
            }
            if (!br.empty()) {
                this->branch[i] = br.back().second;
            }
            br.emplace_back() = {this->ends[i], i};
        }
        this->idx = 0;
    }

    inline void upperBound(const S x) noexcept override {
         size_t i = 0;
         const size_t n_intervals = this->starts.size();
         while (i < n_intervals) {
             if (eytz[i] > x) {
                 i = 2 * i + 1;
             } else {
                 i = 2 * i + 2;
             }
         }
         int shift = __builtin_ffs(~(i + 1));
         size_t best_idx = (i >> shift) - ((shift > 1) ? 1 : 0);
         this->idx = (best_idx < n_intervals) ? eytz_index[best_idx] : n_intervals - 1;
         if (this->idx > 0 && this->starts[this->idx] > x) {
             --this->idx;
         }
    }

private:
    std::vector<S> eytz;
    std::vector<size_t> eytz_index;

    size_t eytzinger_helper(S* arr, size_t n, size_t i, size_t k) {
        if (k < n) {
            i = eytzinger_helper(arr, n, i, 2*k+1);
            eytz[k] = this->starts[i];
            eytz_index[k] = i;
            ++i;
            i = eytzinger_helper(arr, n, i, 2*k + 2);
        }
        return i;
    }

    int eytzinger(S* arr, size_t n) {
        return eytzinger_helper(arr, n, 0, 0);
    }
};
