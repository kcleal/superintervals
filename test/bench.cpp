
#include "IITree.hpp"
#include "IntervalTree.h"
#include "matrylist.hpp"

#include <chrono>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <utility>
#include <unordered_set>

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
    high_resolution_clock::time_point t0 = high_resolution_clock::now();
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
    std::cout << uSec(t0) << "ms, N ref intervals " << intervals.size() << " N queries " << queries.size() << std::endl;
}


void run_tools(std::vector<BedInterval>& intervals, std::vector<BedInterval>& queries) {
    size_t found;
    int index;
    high_resolution_clock::time_point t0, t1;
    std::vector<size_t> a, b;

    std::cout << "\n MatryList \n";
    MatryList<int, int> itv;
    t0 = high_resolution_clock::now();
    itv = MatryList<int, int>();
    index = 0;
    found = 0;
    t1 = high_resolution_clock::now();
    for (const auto& item : intervals) {
        itv.add(item.start, item.end, index);
        index += 1;
    }
    itv.index();
    std::cout << uSec(t0) << " ms construct" << std::endl;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
//        itv.findOverlaps(item.start, item.end, a);
//        found += a.size();
        found += itv.countOverlaps(item.start, item.end);

    }
    std::cout << uSec(t1) << " ms query" << std::endl;
    std::cout << uSec(t0) << " ms total" << std::endl;
    std::cout << "Found: " << found << std::endl;

    std::cout << "\n IITree \n";
    t0 = high_resolution_clock::now();
    IITree<int, int> tree;
    index = 0;
    found = 0;
    for (const auto& item : intervals) {
        tree.add(item.start, item.end, index);
        index += 1;
    }
    tree.index();
    std::cout << uSec(t0) << " ms construct" << std::endl;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        tree.overlap(item.start, item.end, b);
        found += b.size();
    }
    std::cout << uSec(t1) << " ms query" << std::endl;
    std::cout << uSec(t0) << " ms total" << std::endl;
    std::cout << "Found: " << found << std::endl;

    std::cout << "\n IntervalTree \n";
    t0 = high_resolution_clock::now();
    std::vector<Interval<int, int>> intervals2;
    index = 0;
    found = 0;
    for (const auto& item : intervals) {
        intervals2.push_back(Interval<int, int>(item.start, item.end, index));
        index += 1;
    }
    IntervalTree<int, int> tree2(std::move(intervals2));
    std::cout << uSec(t0) << " ms construct" << std::endl;
    std::vector<Interval<int, int>> results;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        auto result = tree2.findOverlapping(item.start, item.end);
        results.insert(results.end(), result.begin(), result.end());
        found += result.size();
    }
    std::cout << uSec(t1) << " ms query" << std::endl;
    std::cout << uSec(t0) << " ms total" << std::endl;
    std::cout << "Found: " << found << std::endl;

}

int main(int argc, char *argv[]) {
    if (argc < 3) {
		printf("Usage: run-cpp-libs <loaded.bed> <streamed.bed>\n");
		return 0;
	}
    std::vector<BedInterval> intervals;
    std::vector<BedInterval> queries;

    std::cout << "\n****** Reads+genes2 ******\n";
    load_intervals(std::string(argv[1]), std::string(argv[2]), intervals, queries); //, shuffle);
    run_tools(intervals, queries);

}