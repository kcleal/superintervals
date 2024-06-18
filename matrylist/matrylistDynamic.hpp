
#pragma once

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <set>

// S for scalar for start, end. T for data type
template<typename S, typename T>
class MatryListDynamic {
    public:
    struct Interval {
        S start, end;
        T value;
        typename std::multiset<Interval>::iterator branch;
        Interval() : start(0), end(0), value(T{}), branch() {}
        Interval(S s, S e, T v, typename std::multiset<Interval>::iterator endIter)
            : start(s), end(e), value(v), branch(endIter) {}
    };
    struct CustomComparator {
        bool operator()(const Interval& a, const Interval& b) const {
            if (a.start != b.start) {
                return a.start < b.start; // Compare by start value first
            }
            return a.end > b.end; // If start values are equal, compare by end value
        }
    };
    std::multiset<Interval, CustomComparator> intervals;
    Interval tmpInterval;

    struct IntervalItem {
        S start, end;
        T value;
    };

    MatryListDynamic() {}
    ~MatryListDynamic() = default;

    void clear() {
        intervals.clear();
    }

    size_t size() {
        return intervals.size();
    }

    void add(S start, S end, const T& value) {
        intervals.insert({start, end, value, intervals.end()});
    }

    void index() {
        if (intervals.size() < 2) {
            return;
        }
        auto i = intervals.begin();
        auto j = i;
        ++j;
        S last_end = i->end;
        while (j != intervals.end()) {
            if (j->end < i->end) {
                auto j_handle = intervals.extract(j);
                j_handle.value().branch = i;
                j = intervals.insert(std::move(j_handle));
                last_end = (i->end > last_end) ? i->end : last_end;
            } else if (j->end < last_end) {
                auto j_handle = intervals.extract(j);
                j_handle.value().branch = i->branch;
                j = intervals.insert(std::move(j_handle));
            }
            ++i; ++j;
        }
    }

    bool anyOverlaps(S start, S end) {
        tmpInterval.start = end;
        tmpInterval.end = end;
        auto it = intervals.upper_bound(tmpInterval);
        if (it != intervals.begin() && (it == intervals.end() || it->start > end)) {
            --it;
        }
        return start <= it->end;
    }

    void findOverlaps(S start, S end, std::vector<IntervalItem>& found_intervals) {
        if (intervals.empty()) {
            return;
        }
        if (!found_intervals.empty()) {
            found_intervals.clear();
        }
        tmpInterval.start = end;
        tmpInterval.end = end;
        auto it = intervals.upper_bound(tmpInterval);
        if (it != intervals.begin() && (it == intervals.end() || it->start > end)) {
            --it;
        }

        while (it != intervals.begin()) {
            if (start <= it->end) {
                found_intervals.push_back({it->start, it->end, it->value});
                --it;
            } else {
                if (it->branch == intervals.end() || it->branch->start > it->start) {
                    break;
                }
                it = it->branch;
            }
        }
        if (it == intervals.begin() && it->start < end) {
            found_intervals.push_back({it->start, it->end, it->value});
        }
    }

    size_t countOverlaps(S start, S end) {
        size_t count = 0;
        if (intervals.empty()) {
            return count;
        }
        tmpInterval.start = end;
        tmpInterval.end = end;
        auto it = intervals.upper_bound(tmpInterval);
        if (it != intervals.begin() && (it == intervals.end() || it->start > end)) {
            --it;
        }
        while (it != intervals.begin()) {
            if (start <= it->end) {
                count += 1;
                --it;
            } else {
                if (it->branch == intervals.end() || it->branch->start > it->start) {
                    break;
                }
                it = it->branch;
            }
        }
        if (it == intervals.begin() && it->start < end) {
            count += 1;
        }
        return count;
    }

};
