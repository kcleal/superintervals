
#pragma once

#include <vector>
#include <iostream>
#ifdef __AVX2__
    #include <immintrin.h>
#elif defined __ARM_NEON
    #include <arm_neon.h>
#endif

// S for scalar for start, end. T for data type
template<typename S, typename T>
class MatryList {
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
    MatryList()
        : idx(0)
        , n_intervals(0)
        , startSorted(true)
        , endSorted(true)
        , it_low(0)
        , it_high(0)
        {}

    ~MatryList() = default;

    struct IntervalItem {
        S start, end;
        T data;
    };
    class Iterator {
    public:
        Iterator(const MatryList* list, size_t index) : matry(list) {
            _start = list->it_low;
            _end = list->it_high;
            it_index = index;
        }
        typename MatryList::IntervalItem operator*() const {
            return typename MatryList<S, T>::IntervalItem{matry->starts[it_index], matry->ends[it_index], matry->data[it_index]};
        }
        Iterator& operator++() {
            if (it_index == 0) {
                it_index = static_cast<size_t>(-1);
                return *this;
            }
            if (it_index > 0) {
                if (_start <= matry->ends[it_index]) {
                    --it_index;
                } else {
                    if (matry->branch[it_index] >= it_index) {
                        it_index = 0;
                        return *this;
                    }
                    it_index = matry->branch[it_index];
                    if (_start <= matry->ends[it_index]) {
                        --it_index;
                    } else {
                        it_index = 0;
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
        Iterator begin() const { return Iterator(matry, matry->idx); }
        Iterator end() const { return Iterator(matry, 0); }
    private:
        S _start, _end;
        const MatryList<S, T>* matry;
        size_t it_index;
    };

    Iterator begin() const { return Iterator(this, idx); }
    Iterator end() const { return Iterator(this, 0); }

    void searchInterval(S start, S end) {
        it_low = start; it_high = end;
        upperBound(end);
    }

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

    inline void upperBound(S value) {  // https://academy.realm.io/posts/how-we-beat-cpp-stl-binary-search/
        size_t size = n_intervals;
        idx = 0;
        while (size > 0) {
            size_t half = size / 2;
            size_t other_half = size - half;
            size_t probe = idx + half;
            size_t other_low = idx + other_half;
            S v = starts[probe];
            size = half;
            idx = v <= value ? other_low : idx;
        }
        if (idx > 0 && (idx == n_intervals || starts[idx] > value)) {
            --idx;
        }
    }

    void sortIntervals() {
        size_t i;
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
                auto block_start = it_start;
                auto block_end = it_start + 1;
                bool srt = false;
                while (block_end != intervals.end() && block_end->start == block_start->start) {
                    if (block_end->end > (block_end - 1)->end) {
                        srt = true;
                    }
                    ++block_end;
                }
                if (srt) {
                    std::sort(block_start, block_end, [](const Interval& a, const Interval& b) { return a.end > b.end; });
                    i = std::distance(intervals.begin(), block_start);
                    while (true) {
                        data[i++] = data_copy[block_start->index];
                        if (block_start == block_end) {
                            break;
                        }
                        ++block_start;
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
        if (n_intervals < 2) {
            return;
        }
        starts.resize(n_intervals);
        ends.resize(n_intervals);
        if (!branch.empty()) {
            branch.clear();
        }

        branch.resize(intervals.size(), 0);

        sortIntervals();
        S last_end = ends[0];
        size_t j = 1;
        for (size_t i=0; i < ends.size() - 1; ++i) {
            if (ends[j] < ends[i]) {
                branch[j] = i;
                last_end = (ends[i] > last_end) ? ends[i] : last_end;
            } else if (ends[j] < last_end) {
                branch[j] = branch[i];
            }
            ++j;
        }
        idx = 0;
    }

    IntervalItem at(size_t index) {
        return IntervalItem(starts[index], ends[index], data[index]);
    }

    void at(size_t index, IntervalItem& itv) {
        itv.start = starts[index];
        itv.end = ends[index];
        itv.data = data[index];
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

        while (i != 0) {
            if (start <= ends[i]) {
                found.push_back(data[i]);
                --i;
            } else {
                i = branch[i];
            }
        }
        if (start <= ends[0] && starts[0] <= end) {
            found.push_back(data[0]);
        }
    }

    size_t countOverlaps(S start, S end) {
        if (!n_intervals) {
            return 0;
        }
        upperBound(end);
        size_t found = 0;
        size_t i = idx;
        constexpr size_t block = 32;

#ifdef __AVX2__
        __m256i start_vec = _mm256_set1_epi32(start);
        constexpr size_t simd_width = 256 / (sizeof(S) * 8);
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
            if (count < block) {
                break;
            }
        }
#elif defined __ARM_NEON
        int32x4_t start_vec = vdupq_n_s32(start);
        constexpr size_t simd_width = 128 / (sizeof(S) * 8);
        uint32x4_t ones = vdupq_n_u32(1);
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
            if (count < block) {
                break;
            }
        }
#else // rely on auto-vectorization
        while (i > block) {
            size_t count = 0;
            for (size_t j = i; j > i - block; --j) {
                count += (start <= ends[j]) ? 1 : 0;
            }
            found += count;
            i -= block;
            if (count < block) {
                break;
            }
        }
#endif
        // process remaining
        while (i != 0) {
            if (start <= ends[i]) {
                ++found;
                --i;
            } else {
                i = branch[i];
            }
        }
        if (start <= ends[0] && starts[0] <= end) {
            ++found;
        }
        return found;
    }
};
