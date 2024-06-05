
#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>


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
    size_t idx, jdx, n_intervals;
    bool startSorted, endSorted;
    S last_start, last_end;
    MatryList()
        : idx(0)
        , jdx(0)
        , n_intervals(0)
        , startSorted(true)
        , endSorted(true)
        , last_start(std::numeric_limits<S>::min())
        , last_end(std::numeric_limits<S>::min())
        {}

    ~MatryList() = default;

    void clear() {
        idx = 0;
        intervals.clear(); data.clear(); starts.clear(); ends.clear();
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

    std::string intervalStr(Interval& v) {
        std::string s = "(" + std::to_string(v.start) + ", " + std::to_string(v.end) + ")";
        return s;
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
        // sort data by start, sort data by end and keep the indexes
        if (!startSorted || !endSorted) {
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
        for (i=0; i < intervals.size() - 1; ++i) {
            for (size_t j=i + 1; j < intervals.size(); ++j) {
                if (intervals[j].end >= intervals[i].end) {
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
        if (start < last_start || end < last_end || starts[jdx] <= end) {
            upper_bound(end);
        }
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
        last_start = start;
        last_end = end;
    }

    size_t countOverlaps(S start, S end) {
        if (!n_intervals) {
            return 0;
        }
        // try and re-use upper bound
        if (start < last_start || end < last_end || starts[jdx] <= end) {
            upper_bound(end);
        }
        size_t found = 0;
        size_t i = idx;
        while (i > 0) {
            if (start <= ends[i--]) {
                ++found;
            } else {
                if (++i; branch[i] >= i) {
                    break;
                }
                i = branch[i];
            }
        }
        if (i==0 && start <= ends[0]) {
            ++found;
        }
        last_start = start;
        last_end = end;
        return found;
    }
};


// S for scalar for start, end. T for data type
template<typename S, typename T>
class MatryListE {
    public:
    struct Interval {
        S start, end;
    };
    alignas(64) std::vector<S> starts;
    alignas(64) std::vector<S> ends;
    alignas(64) std::vector<S> eytz;
    alignas(64) std::vector<size_t> eytz_index;
    alignas(64) std::vector<size_t> branch;
    std::vector<Interval> intervals;
    std::vector<T> data;
    size_t idx, jdx, n_intervals;
    bool startSorted, endSorted;
    S last_start, last_end;
    MatryListE()
        : idx(0)
        , jdx(0)
        , n_intervals(0)
        , startSorted(true)
        , endSorted(true)
        , last_start(std::numeric_limits<S>::min())
        , last_end(std::numeric_limits<S>::min())
        {}

    ~MatryListE() = default;

    void clear() {
        idx = 0;
        intervals.clear(); data.clear(); starts.clear(); ends.clear(); eytz.clear(); eytz_index.clear();
    }

    void reserve(size_t n) {
        intervals.reserve(n); data.reserve(n); starts.reserve(n); ends.reserve(n); eytz.reserve(n); eytz_index.reserve(n);
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

    std::string intervalStr(Interval& v) {
        std::string s = "(" + std::to_string(v.start) + ", " + std::to_string(v.end) + ")";
        return s;
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
        // sort data by start, sort data by end and keep the indexes
        if (!startSorted || !endSorted) {
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
        for (i=0; i < intervals.size() - 1; ++i) {
            for (size_t j=i + 1; j < intervals.size(); ++j) {
                if (intervals[j].end >= intervals[i].end) {
                    break;
                }
                branch[j] = i;
            }
        }
        idx = 0;

        eytz.resize(n_intervals);
        eytz_index.resize(n_intervals);
        eytzinger(&starts[0], n_intervals);
    }

    size_t eytzinger_helper(S* arr, size_t n, size_t i, size_t k) {
        if (k < n) {
            i = eytzinger_helper(arr, n, i, 2*k+1);
            eytz[k] = starts[i];
            eytz_index[k] = i;
            ++i;
            i = eytzinger_helper(arr, n, i, 2*k + 2);
        }
        return i;
    }

    size_t eytzinger(S* arr, size_t n) {
      return eytzinger_helper(arr, n, 0, 0);
    }

    void eytzinger_upper_bound(S x) {
        size_t i = 0;
        size_t best_idx = n_intervals;
        idx = n_intervals - 1;
        jdx = idx;
        while (i < n_intervals) {
            if (eytz[i] > x) {
                // weirdly this runs faster when this branch is included, better branch prediction?
                if (best_idx == n_intervals || eytz[i] <= eytz[best_idx]) {
                    best_idx = i;  // best candidate closer to x
                }
                i = 2 * i + 1;
            } else {

                i = 2 * i + 2;
            }
        }
        if (best_idx <= idx) {
            idx = eytz_index[best_idx];
            jdx = idx;
            if (idx > 0) {
                if (starts[idx] > x) {
                    --idx;
                }
            }
        }
    }

    void findOverlaps(S start, S end, std::vector<size_t>& found) {
        if (!n_intervals) {
            return;
        }
        if (!found.empty()) {
            found.clear();
        }
        if (start < last_start ||  end < last_end || starts[jdx] <= end) {
            eytzinger_upper_bound(end);
        }
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
        last_start = start;
        last_end = end;
    }

    size_t countOverlaps(S start, S end) {
        if (!n_intervals) {
            return 0;
        }
        // try and re-use upper bound
        if (start < last_start || end < last_end || starts[jdx] <= end) {
            eytzinger_upper_bound(end);
        }
        size_t found = 0;
        size_t i = idx;
        while (i > 0) {
            if (start <= ends[i--]) {
                ++found;
            } else {
                if (++i; branch[i] >= i) {
                    break;
                }
                i = branch[i];
            }
        }
        if (i==0 && start <= ends[0]) {
            ++found;
        }
        last_start = start;
        last_end = end;
        return found;
    }
};

