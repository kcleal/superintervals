
#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cassert>
#include <queue>
#include <iostream>

#include <bitset>

// S for scalar for start, end. T for data type
template<typename S, typename T>
class MatryList {
    public:
    struct Interval {
        S start, end;
    };

    alignas(64) std::vector<S> starts;
    alignas(64) std::vector<S> ends;
    alignas(64) std::vector<S> eytz;
    alignas(64) std::vector<size_t> eytz_index;
    std::vector<Interval> intervals;
    alignas(64) std::vector<size_t> branch;
    std::vector<T> data;
    S distance_threshold;
    size_t idx, n_intervals;
    bool startSorted, endSorted;
    size_t n = n_intervals;
    size_t best_idx_eytz;
    S max_eytz;
    typename std::vector<S>::iterator it_begin, it_end;
    MatryList()
        : distance_threshold(25000)
        , idx(0)
        , n_intervals(0)
        , startSorted(true)
        , endSorted(true)
        , max_eytz(0)
        {}

    ~MatryList() = default;

    struct tmpItem {
        S start, end;
        T data;
    };

    void clear() {
        idx = 0;
        intervals.clear();
        data.clear();
        starts.clear();
        ends.clear();
        eytz.clear();
        eytz_index.clear();
    }

    void reserve(size_t n) {
        intervals.reserve(n);
        data.reserve(n);
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
        ends.push_back(end);
        data.emplace_back(value);
    }

    inline bool is_overlapping_interval(const S x1, const S x2, const Interval& itv) noexcept {
        return std::max(x1, itv.start) <= std::min(x2, itv.end);
    }

    // https://academy.realm.io/posts/how-we-beat-cpp-stl-binary-search/
    inline void fast_upper_bound(S start_val, S value) {

        if (start_val >= starts[idx] && start_val - starts[idx] < distance_threshold) {
            while (idx < n_intervals && intervals[idx].start <= value) {
                ++idx;
            }
            if (idx) {
                --idx;
            }
//            return;
//            }
//            else if (start_val < starts[idx] && starts[idx] - start_val < distance_threshold) {
//                while (intervals[idx].start > value) {
//                    if (idx) {
//                        --idx;
//                    } else {
//                        break;
//                    }
//                }
        } else {
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
    }

    inline void _binary_search(S pos) noexcept {
//            auto up = last_q_start < pos
//                         ? std::upper_bound(it_begin + idx, it_end, pos)
//                         : std::upper_bound(it_begin, it_begin + idx, pos);
        auto up = std::upper_bound(starts.begin(), starts.end(), pos);
        if (up != it_begin && (up == it_end || *up > pos)) {
            --up;
        }
        idx = std::distance(starts.begin(), up);
    }

    std::string intervalStr(Interval& v) {
        std::string s = "(" + std::to_string(v.start) + ", " + std::to_string(v.end) + ")";
        return s;
    }

    void index() {
        size_t i;
        n_intervals = intervals.size();
        if (n_intervals < 2) {
            return;
        }
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

        starts.resize(n_intervals);
        i = 0;
        for (const auto& item: intervals) {
            starts[i] = item.start;
            ++i;
        }
        if (!branch.empty()) {
            branch.clear();
        }
        branch.resize(intervals.size(), SIZE_MAX);

        for (i=0; i < intervals.size() - 1; ++i) {
            for (size_t j=i + 1; j < intervals.size(); ++j) {
                if (intervals[j].end >= intervals[i].end) {
//                if (intervals[j].end > intervals[i].end) {
                    break;
                }
                branch[j] = i;
//                if (branch[j] == SIZE_MAX || i > branch[j]) {
//                    branch[j] = i;
//                }
            }
        }

        idx = 0;
        eytz.resize(n_intervals);
        eytz_index.resize(n_intervals);
        eytzinger(&starts[0], n_intervals);

//        eytz.resize(n_intervals + 1);
//        eytz_index.resize(n_intervals + 1);
//        eytzinger(&starts[0], n_intervals + 1);


//        std::cout << "\nEytzinger array:\n";
//        for (auto &item : eytz) {
//            std::cout << item << ", ";
//        } std::cout << std::endl;
//        std::cout << " best eytz idx=" << best_idx_eytz << std::endl;

//            size_t jj = 0;
//            for (auto item: branch) {
//                std::cout << item << ", ";
//                jj += 1;
//                if (jj > 500) {
//                    break;
//                }
//            } std::cout << std::endl << std::endl << std::endl;
//            std::terminate();
    }

//    size_t eytzinger_helper(S* arr, size_t n, size_t i, size_t k) {
//        if (k < n) {
//            i = eytzinger_helper(arr, n, i, 2*k);
//            eytz[k] = starts[i];
//            eytz_index[k] = i;
//            ++i;
//            i = eytzinger_helper(arr, n, i, 2*k + 1);
//        }
//        return i;
//    }
//
//    int eytzinger(S* arr, size_t n) {
//      return eytzinger_helper(arr, n, 0, 1);
//    }

    void upper_bound2_new3(S x) {  // need to change construction method for this one to work
        size_t i = 1;
        while (i < n_intervals) {
            i = 2 * i + (eytz[i] < x);
//            if (eytz[i] > x) {
//                i = 2 * i;
//            } else {
//                i = 2 * i + 1;
//            }
        }
        i >>= __builtin_ffs(~i);
        idx = (i == 0) ? n_intervals - 1 : eytz_index[i];
        if (idx > 0 && starts[idx] > x) {
            --idx;
        }
    }

    size_t eytzinger_helper(S* arr, size_t n, size_t i, size_t k) {
        if (k < n) {
            i = eytzinger_helper(arr, n, i, 2*k+1);
            eytz[k] = starts[i];
            eytz_index[k] = i;
//            if (starts[i] > max_eytz) {
//                max_eytz = starts[i];
//                best_idx_eytz = k;
//            }
            ++i;
            i = eytzinger_helper(arr, n, i, 2*k + 2);
        }
        return i;
    }

    int eytzinger(S* arr, size_t n) {
      return eytzinger_helper(arr, n, 0, 0);
    }

    void upper_bound2_new2(S x) {
        size_t i = 0;
//        size_t best_idx = n_intervals;
        while (i < n_intervals) {
            if (eytz[i] > x) {
                i = 2 * i + 1;
            } else {
//                if (best_idx <= n_intervals) {
//                    best_idx = i;
//                }
                i = 2 * i + 2;
            }
        }
//        if (i == n_intervals + 1) {
//            idx = n_intervals;
//        } else {
//            ++i;
//            idx = (i >> __builtin_ffs(~i)) - 1;
//        }

//        idx = (idx < n_intervals) ? eytz_index[idx] : n_intervals - 1;
//        if (idx > 0 && starts[idx] > x) {
//            --idx;
//        }

//        idx = best_idx;
        idx = (i == n_intervals + 1) ? n_intervals : (i >> __builtin_ffs(~i)) - 1;

//        std::cout << idx << " " << best_idx << " " << n_intervals << std::endl;
        idx = (idx < n_intervals) ? eytz_index[idx] : n_intervals - 1;
        if (idx > 0 && starts[idx] > x) {
            --idx;
        }

        // i >>= __builtin_ffs(~i);
//        size_t shift = __builtin_ffs(~(i+1));
//        std::cout << "i=" << std::bitset<8>(i) << " " << std::bitset<8>(i+1) << std::endl;
//        std::cout << n_intervals << " " << i << " shift=" << shift << std::endl;
//        size_t j = (i == n_intervals + 1) ? n_intervals : ((i+1) >> shift) - 1;
//        size_t j = i;
//        j >>= __builtin_ffs(~(i << 1));
//        best_idx = j;
//        std::cout << i <<" best idx=" << best_idx << " "
//             << __builtin_ffs(~i) << " " << j << std::endl;
//        idx = (best_idx < n_intervals) ? eytz_index[best_idx] : n_intervals - 1;
//        if (idx > 0 && starts[idx] > x) {
//            --idx;
//        }
    }

    void upper_bound2(S x) {
        size_t i = 0;
        size_t best_idx = n_intervals;
        idx = n_intervals - 1;
        while (i < n_intervals) {
            if (eytz[i] > x) {
                //  best_idx = i;
                // runs faster when this branch is included
                // probably due to better branch prediction
                if (best_idx == n_intervals || eytz[i] <= eytz[best_idx]) {
                    best_idx = i;  // best candidate closer to x
                }
                i = 2 * i + 1;
            } else {
                i = 2 * i + 2;
            }
        }
//        idx = (best_idx < idx) ? eytz_index[best_idx] : idx;
//        if (idx > 0 && starts[idx] > x) {
//            --idx;
//        }
        if (best_idx < idx) {
            idx = eytz_index[best_idx];
            if (idx > 0) {
                if (starts[idx] > x) {
                    --idx;
                }
            }
        }
    }

    void search_overlap(S start, S end, std::vector<size_t>& found) {
        if (!n_intervals) {
            return;
        }
        if (!found.empty()) {
            found.clear();
        }
//            _binary_search(end);
//        fast_upper_bound(start, end);
        upper_bound2(end);
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
//
//        while (true) {
//            if (is_overlapping_interval(start, end, intervals[i])) {
//                found.push_back(i);
//            } else {
//                if (branch[i] >= i) {
//                    break;
//                }
//                i = branch[i] + 1;
//            }
//            if (i == 0) {
//                break;
//            }
//            --i;
//        }
    }

    size_t countOverlapping(S start, S end) {
        if (!n_intervals) {
            return 0;
        }
        upper_bound2(end);
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
        return found;
    }
};


