
#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#ifdef __AVX2__
    #include <immintrin.h>
#endif

// S for scalar for start, end. T for data type
template<typename S, typename T>
class MatryList {
    public:
    struct Interval {
        S start, end;
    };
    alignas(sizeof(S)) std::vector<S> starts;
    alignas(sizeof(S)) std::vector<S> ends;
    alignas(sizeof(size_t)) std::vector<size_t> branch;
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
        upper_bound(end);
    }

    void clear() {
        idx = 0;
        intervals.clear(); data.clear(); starts.clear(); branch.clear();
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
        intervals.emplace_back() = {start, end};
        data.emplace_back(value);
    }

    inline void upper_bound(S value) {  // https://academy.realm.io/posts/how-we-beat-cpp-stl-binary-search/
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

    struct tmpItem {
        S start, end;
        T data;
    };

    void index() {
        size_t i;
        n_intervals = intervals.size();
        if (n_intervals < 2) {
            return;
        }
        starts.resize(n_intervals);
        ends.resize(n_intervals);
        if (!startSorted || !endSorted) {  // sort data by start, sort data by end and keep the indexes
            std::vector<tmpItem> tmp(intervals.size());
            i = 0;
            for (const auto &itv : intervals) {
                tmp[i].start = itv.start;
                tmp[i].end = itv.end;
                tmp[i].data = data[i];
                ++i;
            }
            if (!startSorted) {
                std::sort(tmp.begin(), tmp.end(),
                [](const tmpItem& a, const tmpItem& b) {
                    return (a.start < b.start || (a.start == b.start && a.end > b.end)); });
            } else {  // only sort parts that need sorting - ends in descending order
                auto it_start = tmp.begin();
                while (it_start != tmp.end()) {
                    auto block_start = it_start;
                    auto block_end = it_start + 1;
                    bool srt = false;
                    while (block_end != tmp.end() && block_end->start == block_start->start) {
                        if (block_end->end > (block_end - 1)->end) {
                            srt = true;
                        }
                        ++block_end;
                    }
                    if (srt) {
                        std::sort(block_start, block_end, [](const tmpItem& a, const tmpItem& b) { return a.end > b.end; });
                    }
                    it_start = block_end;
                }
            }
            i = 0;
            for (const auto& itv : tmp) {
                intervals[i].start = itv.start;
                intervals[i].end = itv.end;
                data[i] = itv.data;
                ends[i] = itv.end;
                ++i;
            }
            startSorted = true;
            endSorted = true;
        }
        i = 0;
        for (const auto& item: intervals) {
            starts[i] = item.start;
            ends[i] = item.end;
            ++i;
        }
        if (!branch.empty()) {
            branch.clear();
        }
        branch.resize(intervals.size(), SIZE_MAX);
//        for (i=0; i < ends.size() - 1; ++i) {
//            for (size_t j=i + 1; j < ends.size(); ++j) {
//                if (ends[j] >= ends[i]) {
//                    break;
//                }
//                branch[j] = i;
//            }
//        }

        for (i=0; i < ends.size() - 1; ++i) {
            if (ends[i+1] < ends[i]) {
                branch[i+1] = i;
            } else {
                if (branch[i] != SIZE_MAX && ends[i+1] < ends[branch[i]]) {
                    branch[i+1] = branch[i];
                }
            }
        }


// Alternative construction, but slower
//        std::vector<std::pair<S, size_t>> br;
//        br.reserve(1000);
//        br.push_back({ends[0], 0});
//        for (i=1; i < ends.size(); ++i) {
//            while (!br.empty() && br.back().first < ends[i]) {
//                br.pop_back();
//            }
//            if (!br.empty()) {
//                branch[i] = br.back().second;
//            }
//            br.push_back({ends[i], i});
//        }
        idx = 0;
    }

    void findOverlaps(S start, S end, std::vector<size_t>& found) {
        if (!n_intervals) {
            return;
        }
        if (!found.empty()) {
            found.clear();
        }
        upper_bound(end);

        size_t i = idx;
        while (i > 0) {
            if (start <= ends[i--]) {
                found.push_back(i+1);
            } else {
                if (++i; branch[i] >= i) {
                    break;
                }
                i = branch[i];
            }
        }
        if (!i && start <= ends[0]) {
            found.push_back(0);
        }
    }

    size_t countOverlaps(S start, S end) {
        if (!n_intervals) {
            return 0;
        }
        upper_bound(end);
        size_t found = 0;
        size_t i = idx;
        constexpr size_t block = 32;

#ifdef __AVX2__
        __m256i start_vec = _mm256_set1_epi32(start);
        constexpr size_t simd_width = 256 / (sizeof(S) * 8);
#endif

        while (i > 0) {
            if (start > ends[i--]) {
                if (++i; branch[i] >= i) {
                    break;
                }
                i = branch[i];
            } else {
                ++found;

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
            }
        }
        if (i==0 && start <= ends[0]) {
            ++found;
        }
        return found;
    }

};
