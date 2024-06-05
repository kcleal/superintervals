
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



void load_intervals(const std::string& intervals_file,
                    const std::string& queries_file,
                    std::vector<BedInterval>& intervals,
                    std::vector<BedInterval>& queries,
                    bool shuffle) {
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

    if (shuffle) {
//        std::random_device rd;
        std::mt19937 g(12345);
        std::shuffle(queries.begin(), queries.end(), g);
    } else {
        std::sort(queries.begin(), queries.end(), [](const BedInterval& a, const BedInterval& b) {
            return (a.start < b.start || (a.start == b.start && a.end > b.end));
        });
    }

//    std::sort(intervals.begin(), intervals.end(), [](const BedInterval& a, const BedInterval& b) {
//        return (a.start < b.start || (a.start == b.start && a.end > b.end));
//    });
    std::cout << " N ref intervals " << intervals.size() << " N queries " << queries.size() << std::endl;

}


void print_vec(std::vector<size_t>& a, MatryList<int, int>& itv) {
    std::cout << " Found:\n";
    for (auto item : a) {
        std::cout << item << " - " << itv.intervals[item].start << " " << itv.intervals[item].end << std::endl;
    }
}

size_t uSec(high_resolution_clock::time_point& t0) {
    return duration_cast<milliseconds>(high_resolution_clock::now() - t0).count();
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

//    itv.add(0, 250000000, -1);
    for (const auto& item : intervals) {
        itv.add(item.start, item.end, index);
        index += 1;
    }
    itv.index();
    std::cout << uSec(t0) << " construct ms" << std::endl;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
//        itv.search_overlap(item.start, item.end, a);
        //found += a.size();
        found += itv.countOverlapping(item.start, item.end);

    }
    std::cout << uSec(t1) << " query ms" << std::endl;
    std::cout << uSec(t0) << " total ms" << std::endl;
    std::cout << "Found: " << found << std::endl;
//    return;

    std::cout << "\n IITree \n";
    t0 = high_resolution_clock::now();
    IITree<int, int> tree;
    index = 0;
    found = 0;
//    tree.add(0, 250000000, -1);
    for (const auto& item : intervals) {
        tree.add(item.start, item.end, index);
        index += 1;
    }
    tree.index();
    std::cout << uSec(t0) << " construct ms" << std::endl;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        tree.overlap(item.start, item.end, b);
        found += b.size();
    }
    std::cout << uSec(t1) << " query ms" << std::endl;
    std::cout << uSec(t0) << " total ms" << std::endl;
    std::cout << "Found: " << found << std::endl;
//    return;

    std::cout << "\n IntervalTree \n";
    t0 = high_resolution_clock::now();
    std::vector<Interval<int, int>> intervals2;
    index = 0;
    found = 0;
//    intervals2.push_back(Interval<int, int>(0, 250000000, -1));
    for (const auto& item : intervals) {
        intervals2.push_back(Interval<int, int>(item.start, item.end, index));
        index += 1;
    }
    IntervalTree<int, int> tree2(std::move(intervals2));
    std::cout << uSec(t0) << " construct ms" << std::endl;
    std::vector<Interval<int, int>> results;
    t1 = high_resolution_clock::now();
    for (const auto& item : queries) {
        auto result = tree2.findOverlapping(item.start, item.end);
        results.insert(results.end(), result.begin(), result.end());
        found += result.size();
    }
    std::cout << uSec(t1) << " query ms" << std::endl;
    std::cout << uSec(t0) << " total ms" << std::endl;
    std::cout << "Found: " << found << std::endl;

//    size_t mim = 1000000;
//    for (const auto& item : queries) {
//
//        auto result = tree2.findOverlapping(item.start, item.end);
//        results.insert(results.end(), result.begin(), result.end());
//
//        auto last_q_start = itv.last_q_start;
//        itv.search_overlap(item.start, item.end, a);
//
//        if (result.size() != a.size()) {
//            if (result.size()) {
//
//                for (auto itemr : a) {
//                    std::cerr << itv.intervals[itemr].start << " " << itv.intervals[itemr].end << std::endl;
//                }
//                std::cerr << result.size() << " " << a.size() << " query was "
//                   << item.start << " " << item.end << " lqs " << last_q_start << std::endl;
//                for (auto itemr : result) {
//                    std::cerr << itemr.start << " " << itemr.stop << " " << itemr.value << std::endl;
//                }
//                break;
//            }
//            mim = std::min(result.size() , mim);
//        }
//    }
//    std::cerr << "mim " << mim << std::endl;
}

int main(int argc, char *argv[]) {

    if (argc < 3) {
		printf("Usage: run-cpp-libs <loaded.bed> <streamed.bed>\n");
		return 0;
	}

    std::vector<BedInterval> intervals;
    std::vector<BedInterval> queries;
    bool shuffle = false;
//    bool shuffle = true;

    std::cout << "\n****** Reads+genes2 ******\n";
    load_intervals(std::string(argv[1]), std::string(argv[2]), intervals, queries, shuffle);
    run_tools(intervals, queries);

//    run_tools(queries, intervals);

//    std::cout << "\n\n***** Random2 ******\n";
//    load_intervals("a.bed", "b.bed", intervals, queries, shuffle);
//    run_tools(intervals, queries);
//    run_tools(queries, intervals);


}