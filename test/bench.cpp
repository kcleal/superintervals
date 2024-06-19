
#include "IITree.hpp"
#include "IntervalTree.h"
extern "C" {
    #include "intervaldb.h"
}
#include "matrylist.hpp"
#include "matrylistDynamic.hpp"

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
using std::chrono::milliseconds;


struct BedInterval {
    int start;
    int end;
};

size_t uSec(high_resolution_clock::time_point& t0) {
    return duration_cast<milliseconds>(high_resolution_clock::now() - t0).count();
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
        throw std::runtime_error("Failed to open input files");
    }
    std::string line;
    while (std::getline(intervals_stream, line)) {
        std::istringstream iss(line);
        std::string token;
        std::getline(iss, token, '\t');  // Skip the first token (chromosome)
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
        std::getline(iss, token, '\t');
        int start = std::stoi(token);
        std::getline(iss, token, '\t');
        int end = std::stoi(token);
        queries.emplace_back(BedInterval{std::min(start, end), std::max(start, end)});
    }
    std::cerr << "N ref intervals: " << intervals.size() << ", N queries: " << queries.size() << std::endl;
}


void run_tools(std::vector<BedInterval>& intervals, std::vector<BedInterval>& queries) {
    size_t found, index;
    high_resolution_clock::time_point t0, t1;
    std::vector<size_t> a, b;
//    std::unordered_map<size_t, size_t> found_indexes, found_indexes2;

    std::cout << "MatryList\t";
    MatryList<int, size_t> itv;
    t0 = high_resolution_clock::now();
    itv = MatryList<int, size_t>();
    index = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : intervals) {
        itv.add(item.start, item.end, index);
        index += 1;
    }
    itv.index();
    std::cout << uSec(t0) << "\t";  // construct
    found = 0;
//    index = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        itv.findOverlaps(item.start, item.end, a);
//        found_indexes[index] = a.size(); index += 1;
        found += a.size();
    }
    std::cerr << uSec(t1) << "\t" << found << "\t";  // find all overlapping

    found = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        found += itv.countOverlaps(item.start, item.end);
    }
    std::cerr << uSec(t1) << "\t" << found << std::endl;  // count all overlapping

    std::cout << "ImplicitITree\t";
    t0 = high_resolution_clock::now();
    IITree<int, int> tree;
    index = 0;
    for (const auto& item : intervals) {
        tree.add(item.start, item.end, index);
        index += 1;
    }
    tree.index();
    std::cerr << uSec(t0) << "\t";

    found = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        tree.overlap(item.start, item.end, b);
        // Only returns indexes need to enumerate actual data for a fair test
        for (const auto &v: b) {
            index = tree.data(v); // make sure this step is not elided
        }
        found += b.size();
    }
    std::cout << uSec(t1) << "\t" << found << std::endl;


    std::cout << "IntervalTree\t";
    std::vector<interval_tree::Interval<int, int>> intervals2;
    index = 0;
    for (const auto& item : intervals) {
        intervals2.push_back(interval_tree::Interval<int, int>(item.start, item.end, index));
        index += 1;
    }
    t0 = high_resolution_clock::now();
    interval_tree::IntervalTree<int, int> tree2(std::move(intervals2));
    std::cerr << uSec(t0) << "\t";

    found = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        std::vector<interval_tree::Interval<int, int>> result = tree2.findOverlapping(item.start, item.end);
//        found_indexes2[index] = result.size(); index += 1;
        found += result.size();
    }
    std::cerr << uSec(t1) << "\t" << found << std::endl;

//    for (const auto& kv : found_indexes) {
//        if (kv.second != found_indexes2[kv.first])
//        std::cout << kv.first << " " << kv.second << ", " << found_indexes2[kv.second] << std::endl;
//    }

    std::cout << "NestedContList\t";
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
    std::cerr << uSec(t0) << "\t";

    t1 = high_resolution_clock::now();
    IntervalIterator *it;
    IntervalIterator *it_alloc;
    IntervalMap im_buf[1024];
    found = 0;
    for (const auto& item : queries) {
        it_alloc = interval_iterator_alloc();
        it = it_alloc;
        while(it){
            find_intervals(it, item.start, item.end, im, *p_n, sh, *p_nlists, im_buf, 1024, nhits, &it);
            found += *nhits;
        }
        free_interval_iterator(it_alloc);
    }
    std::cerr << uSec(t1) << "\t" << found << std::endl;

}

int main(int argc, char *argv[]) {
    if (argc < 3) {
		printf("Usage: run-cpp-libs <loaded.bed> <streamed.bed>\n");
		return 0;
	}
    std::vector<BedInterval> intervals;
    std::vector<BedInterval> queries;

    load_intervals(std::string(argv[1]), std::string(argv[2]), intervals, queries);
    run_tools(intervals, queries);
}
