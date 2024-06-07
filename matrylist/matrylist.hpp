
#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>

//#include <arm_neon.h>

// S for scalar for start, end. T for data type
template<typename S, typename T>
class MatryList {
    public:
    struct Interval {
        S start, end;
    };
    alignas(64) std::vector<S> starts;
    alignas(64) std::vector<S> ends;
    alignas(64) std::vector<size_t> branch;
    std::vector<Interval> intervals;
    std::vector<T> data;
    size_t idx, n_intervals;
    bool startSorted, endSorted;
    MatryList()
        : idx(0)
        , n_intervals(0)
        , startSorted(true)
        , endSorted(true)
        {}

    ~MatryList() = default;

//    struct IntervalItem {
//        S start, end;
//        T data;
//    };
//    class Iterator {
//    public:
//        Iterator(size_t idx)
//            : it_index(idx) {}
//
//        IntervalItem operator*() const {
//            return IntervalItem{intervals[it_index].start, intervals[it_index].end, data[it_index]};
//        }
//
//        Iterator& operator--() {
//            if (idx > 0) {
//                --idx;
//                if (start <= ends[idx]) {
//                    return *this;
//                } else {
//                    if (branch[idx] >= idx) {
//                        return *this;
//                    }
//                    idx = branch[idx];
//                }
//            }
//
//            --idx;
//
//            while (idx > 0) {
//            if (start <= ends[i--]) {
//                found.push_back(i+1);
//            } else {
//                if (++i; branch[i] >= i) {
//                    break;
//                }
//                i = branch[i];
//            }
//        }
//        if (!i && start <= ends[0]) {
//            found.push_back(0);
//        }
//
//            return *this;
//        }
//
//        bool operator!=(const Iterator& other) const {
//            return it1 != other.it1 || it2 != other.it2 || it3 != other.it3;
//        }
//    private:
//        size_t it_index;
//    };
//
//    Iterator begin() { return Iterator(intervals.data()); }
//    Iterator end() { return Iterator(intervals.data() + intervals.size()); }


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

    // https://academy.realm.io/posts/how-we-beat-cpp-stl-binary-search/
    inline void upper_bound(S value) {
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
        for (i=0; i < ends.size() - 1; ++i) {
            for (size_t j=i + 1; j < ends.size(); ++j) {
                if (ends[j] >= ends[i]) {
                    break;
                }
                branch[j] = i;
            }
        }
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
        const size_t block = 64;
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
        while (i > 0) {
            if (start > ends[i--]) {
                if (++i; branch[i] >= i) {
                    break;
                }
                i = branch[i];
            } else {
                ++found;
            }
        }
        if (i==0 && start <= ends[0]) {
            ++found;
        }
        return found;
    }

//    size_t countOverlaps(S start, S end) {
//
//    }
};
