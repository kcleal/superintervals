#include "intervalstab.hpp"
#include "IITree.hpp"
#include "IntervalTree.h"
extern "C" {
    #include "cgranges.h"
    #include "intervaldb.h"
}
#include "superintervals.hpp"


#include <chrono>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <utility>
#include <unordered_map>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::microseconds;
using std::chrono::milliseconds;

namespace Bench {

size_t found, index;
high_resolution_clock::time_point t0, t1;

// Defines what is stored alongside a test interval
typedef int DType;


struct BedInterval {
    int start;
    int end;
};

size_t uSec(high_resolution_clock::time_point& t0) {
    return duration_cast<microseconds>(high_resolution_clock::now() - t0).count();
}

void branch_factor(std::vector<int>& ends) {
    double avg = 0;
    double max_count = 0;
    std::vector<int> counts(ends.size(), 0);
    for (size_t i=0; i < ends.size() - 1; ++i) {
        for (size_t j=i + 1; j < ends.size(); ++j) {
            if (ends[j] >= ends[i]) {
                break;
            }
            counts[j] += 1;
        }
    }
    double sum_count = 0;
    for (const auto &v : counts) {
        if (v > max_count) {
            max_count = (double)v;
        }
        sum_count += (double)v;
    }
    avg = sum_count / ((double)counts.size() + 0.0000001);
    std::cout << "Avg branching: " << avg << " Max branching: " << max_count << std::endl;
}

void load_intervals(const std::string& intervals_file,
                    const std::string& queries_file,
                    std::vector<BedInterval>& intervals,
                    std::vector<BedInterval>& queries) {
    intervals.clear();
    queries.clear();
    std::ifstream intervals_stream(intervals_file);
    std::ifstream queries_stream(queries_file);
    if (!intervals_stream || !queries_stream) {
        std::cerr << "Failed to open input files\n";
        std::exit(-1);
    }
    std::string line;
    while (std::getline(intervals_stream, line)) {
        std::istringstream iss(line);
        std::string token;
        std::getline(iss, token, '\t');
        if (token != "chr1") continue;
        std::getline(iss, token, '\t');
        int start = std::stoi(token);
        std::getline(iss, token, '\t');
        int end = std::stoi(token);
        intervals.emplace_back(BedInterval{std::min(start, end), std::max(start, end)});
    }
    while (std::getline(queries_stream, line)) {
        std::istringstream iss(line);
        std::string token;
        std::getline(iss, token, '\t');
        if (token != "chr1") continue;
        std::getline(iss, token, '\t');
        int start = std::stoi(token);
        std::getline(iss, token, '\t');
        int end = std::stoi(token);
        queries.emplace_back(BedInterval{std::min(start, end), std::max(start, end)});
    }
}


void run_IITree(std::vector<BedInterval>& intervals, std::vector<BedInterval>& queries) {
    alignas(16) std::vector<DType> a, b;
    a.reserve(10000); b.reserve(10000);
    std::cout << "ImplicitITree-C++,";
    t0 = high_resolution_clock::now();
    IITree<int, int> tree;
    index = 0;
    for (const auto& item : intervals) {
        tree.add(item.start, item.end, index);
        index += 1;
    }
    tree.index();
    std::cerr << uSec(t0) << ",";

    found = 0;
    std::vector<size_t> c;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        tree.overlap(item.start, item.end, c);
//        found += c.size();
        // IITree only returns indexes, need to reference actual data for a fair test
        for (const auto &v: c) {
            a.push_back(tree.data(v));
        }
        found += a.size();
        a.clear();
    }
    std::cout << uSec(t1) << "," << found << std::endl;
}

void run_ITree(std::vector<BedInterval>& intervals, std::vector<BedInterval>& queries) {
    alignas(16) std::vector<DType> a, b;
    a.reserve(10000); b.reserve(10000);
    std::cout << "IntervalTree-C++,";
    std::vector<interval_tree::Interval<int, int>> intervals2;
    index = 0;
    for (const auto& item : intervals) {
        intervals2.push_back(interval_tree::Interval<int, int>(item.start, item.end - 1, index));
        index += 1;
    }
    t0 = high_resolution_clock::now();
    interval_tree::IntervalTree<int, int> tree2(std::move(intervals2));
    std::cerr << uSec(t0) << ",";

    found = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        std::vector<interval_tree::Interval<int, int>> result = tree2.findOverlapping(item.start, item.end - 1);
        found += result.size();
    }
    std::cerr << uSec(t1) << "," << found << std::endl;
}

void run_NCLS(std::vector<BedInterval>& intervals, std::vector<BedInterval>& queries) {
    alignas(16) std::vector<DType> a, b;
    a.reserve(10000); b.reserve(10000);
    std::cout << "NCLS-C,";
    t0 = high_resolution_clock::now();

    int nD = (int)intervals.size();
    int* p_n = (int*)malloc(sizeof(int));
    int* p_nlists = (int*)malloc(sizeof(int));
    int* nhits = (int*)malloc(sizeof(int));
    IntervalMap* im = (IntervalMap*)calloc(1, sizeof(IntervalMap));
    CALLOC(im, nD, IntervalMap);
    index = 0;
    for (const auto& item : intervals) {
        im[index].start  = (int)item.start;
        im[index].end  = (int)item.end;
        im[index].target_id = (int)index;
        im[index].sublist = -1;
        index += 1;
    }
    SublistHeader* sh = build_nested_list(im, nD, p_n, p_nlists);
    std::cerr << uSec(t0) << ",";

    t1 = high_resolution_clock::now();
    IntervalIterator *it;
//    IntervalIterator *it_alloc;
    IntervalIterator *it_alloc = interval_iterator_alloc();
    IntervalMap im_buf[50000];
    found = 0;
    for (const auto& item : queries) {
//        it_alloc = interval_iterator_alloc();
        reset_interval_iterator(it_alloc);
        it = it_alloc;
        while(it){
            find_intervals(it, item.start, item.end, im, *p_n, sh, *p_nlists, im_buf, 50000, nhits, &it);
            found += *nhits;
        }
//        free_interval_iterator(it_alloc);
    }
    std::cerr << uSec(t1) << "," << found << std::endl;
}

void run_SuperIntervals(std::vector<BedInterval>& intervals, std::vector<BedInterval>& queries,
                si::IntervalMap<int, DType> &itv, std::string name
                ) {
    alignas(16) std::vector<DType> a, b;
    a.reserve(10000); b.reserve(10000);

    std::cout << name << ",";
    index = 0;
    t0 = high_resolution_clock::now();
    for (const auto& item : intervals) {
        itv.add(item.start, item.end - 1, index);
        index += 1;
    }
    itv.build();
//    branch_factor(itv.ends);
    std::cout << uSec(t0) << ",";  // construct

    found = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        itv.search_values(item.start, item.end - 1, a);
        found += a.size();
        a.clear();
    }
    std::cerr << uSec(t1) << "," << found << ",";  // find all overlapping

    found = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        itv.search_values_large(item.start, item.end - 1, a);
        found += a.size();
        a.clear();
    }
    std::cerr << uSec(t1) << "," << found << ",";  // find all overlapping large

    found = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        found += itv.count(item.start, item.end - 1);
    }
    std::cerr << uSec(t1) << "," << found << ",";  // count all overlapping

    found = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        found += itv.count_large(item.start, item.end - 1);
    }
    std::cerr << uSec(t1) << "," << found << "\n";  // count all overlapping large

}

void run_tools(std::vector<BedInterval>& intervals, std::vector<BedInterval>& queries) {



    run_ITree(intervals, queries);

    run_IITree(intervals, queries);

    run_NCLS(intervals, queries);

    //auto itv = SuperIntervals<int, size_t>();
    auto itv = si::IntervalMap<int, DType>();
    run_SuperIntervals(intervals, queries, itv, "SuperIntervals-C++");

//    auto itv2 = si::IntervalMapEytz<int, DType>();
//    run_SuperIntervals(intervals, queries, itv2, "SuperIntervalsEytz-C++");

}

} // namespace Bench

int main(int argc, char *argv[]) {
    if (argc < 3) {
		printf("Usage: run-cpp-libs <reference.bed> <query.bed>\n");
		return -1;
	}
    std::vector<Bench::BedInterval> intervals;
    std::vector<Bench::BedInterval> queries;

    Bench::load_intervals(std::string(argv[1]), argv[2], intervals, queries);
    Bench::run_tools(intervals, queries);

    return 0;
}
