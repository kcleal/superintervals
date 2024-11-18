
#include "superintervals.hpp"
#include <iostream>
#include <vector>
#include <cassert>
#include <utility>


void print_vec(std::vector<size_t>& a, SuperIntervals<int, int>& itv) {
    std::cout << " Found:\n";
    for (auto i : a) {
        std::cout << i << " - " << itv.starts[i] << " " << itv.ends[i] << std::endl;
    }
}


void superTests(SuperIntervals<int, int> &itv, std::string name) {
    std::vector<int> a;
    std::pair<size_t, int> cov_res;
    std::cout << "\n" << name << " tests \n";

    std::cout << "0, ";
    itv.add(10, 20, 0);
    itv.add(11, 12, -1);
    itv.add(13, 14, -1);
    itv.add(15, 16, -1);
    itv.add(25, 29, 4);
    itv.index();
    itv.findOverlaps(17, 30, a);
    assert (a[0] == 4); assert (a[1] == 0);
    itv.coverage(10, 29, cov_res);
    assert (cov_res.second == 17);
    itv.clear(); a.clear();

    std::cout << "1, ";
    itv.add(1, 2, 0);
    itv.add(3, 8, -1);
    itv.add(5, 7, -1);
    itv.add(7, 20, 3); //
    itv.add(9, 10, -1);
    itv.add(13, 15, -1);
    itv.add(15, 16, -1);
    itv.add(19, 30, 7);  //
    itv.add(22, 24, -1);
    itv.add(24, 25, -1);
    itv.add(26, 28, -1);
    itv.add(32, 39, -1);
    itv.add(34, 36, -1);
    itv.add(38, 40, -1);
    itv.index();
    itv.findOverlaps(17, 21, a);
    assert (a[0] == 7); assert (a[1] == 3);
    itv.coverage(17, 21, cov_res);
    assert (cov_res.second == 5);
    itv.clear(); a.clear();

    std::cout << "2, ";
    itv.add(0, 250000, 0);
    itv.add(55, 1055, -1);
    itv.add(115, 1115, -1);
    itv.add(130, 1130, -1);
    itv.add(281, 1281, -1);
    itv.add(639, 1639, -1);
    itv.add(842, 1842, -1);
    itv.add(999, 1999, -1);
    itv.add(1094, 2094, -1);
    itv.add(1157, 2157, -1);
    itv.add(1161, 2161, -1);
    itv.add(1265, 2265, -1);
    itv.add(1532, 2532, -1);
    itv.add(1590, 2590, -1);
    itv.add(1665, 2665, -1);
    itv.add(1945, 2945, -1);
    itv.add(2384, 3384, -1);
    itv.add(2515, 3515, -1);
    itv.index();
    itv.findOverlaps(1377, 2377, a);
    assert (a.back() == 0 && a.size() == 12);
    itv.clear(); a.clear();

    std::cout << "3, ";
    itv.add(0, 400, 0);
    itv.add(2, 10, 0);
    itv.add(4, 6, 0);
    itv.add(6, 7, 0);
    itv.add(9, 20, 1);
    itv.add(15, 70, 2);
    itv.add(19, 30, 3);
    itv.add(29, 40, 4);
    itv.add(39, 50, 5);
    itv.add(49, 60, 6);
    itv.add(58, 59, 7);
    itv.index();
    itv.findOverlaps(1, 5, a);
    assert (a.back() == 0 && a.size() == 3);
    itv.clear(); a.clear();

    std::cout << "4, ";
    itv.add(1, 6100000, 0);
    itv.add(4, 5, 6);
    itv.add(6, 7, 7);
    itv.add(9, 10, 7);
    itv.add(11, 12, 7);
    itv.index();
    itv.findOverlaps(2, 25, a);
    assert (a.back() == 0 && a.size() == 5);
    itv.clear(); a.clear();

    std::cout << "5, ";
    itv.add(1, 100, 0);
    itv.add(30, 200, 7);
    itv.add(40, 50, 6);
    itv.add(60, 70, 7);
    itv.index();
    itv.findOverlaps(55, 65, a);
    assert (a.back() == 0 && a.size() == 3);
    itv.coverage(55, 65, cov_res);
    assert (cov_res.second == 25);
    itv.clear(); a.clear();

    std::cout << "6, ";
    itv.add(10, 1001, 0);
    itv.add(30, 400, 1);
    itv.add(60, 700, 2);
    itv.add(65, 80, 3);
    itv.index();
    itv.findOverlaps(100, 200, a);
    assert (a.back() == 0 && a.size() == 3);
    itv.clear(); a.clear();

    std::cout << "7, ";
    itv.add(3, 40, 0);
    itv.add(10, 30, 5);
    itv.add(20, 25, 5);
    itv.add(22, 24, 6);
    itv.index();
    itv.findOverlaps(31, 32, a);
    assert (a.back() == 0 && a.size() == 1);
    itv.clear(); a.clear();

    std::cout << "8, ";
    itv.add(3, 40, 0);
    itv.add(4, 5, 4);
    itv.add(6, 7, 4);
    itv.add(10, 31, 5);
    itv.add(31, 32, 5);
    itv.index();
    itv.findOverlaps(31, 32, a); assert (a.back() == 0 && a.size() == 3); a.clear();
    itv.findOverlaps(10, 11, a); assert (a.back() == 0 && a.size() == 2); a.clear();
    itv.findOverlaps(8, 40, a); assert (a.back() == 0 && a.size() == 3); a.clear();
    itv.findOverlaps(4, 7, a); assert (a.back() == 0 && a.size() == 3); a.clear();
    itv.clear();

    std::cout << "9, ";
    itv.add(3, 40, 0);
    itv.add(3, 40, 4);
    itv.add(3, 40, 4);
    itv.add(3, 4, 4);
    itv.add(35, 50, 4);
    itv.add(40, 400, 5);
    itv.add(40, 400, 4);
    itv.index();
    itv.findOverlaps(38, 41, a); assert (a.back() == 0 && a.size() == 6); a.clear();
    itv.findOverlaps(41, 42, a); assert (a.back() == 4 && a.size() == 3); a.clear();
    itv.findOverlaps(339, 410, a); assert (a.back() == 5 && a.size() == 2); a.clear();
    itv.clear();

    std::cout << "10, ";
    itv.add(0, 100, 0);
    itv.add(10, 110, 2);
    itv.add(10, 20, 1);
    itv.add(30, 40, 3);
    itv.add(35, 135, 4);
    itv.add(45, 55, 5);
    itv.add(60, 160, 6);
    itv.add(70, 80, 7);
    itv.add(90, 190, 8);
    itv.add(110, 120, 9);
    itv.add(130, 140, 10);
    itv.add(150, 250, 11);
    itv.index();
    itv.findOverlaps(95, 105, a); assert (a.front() == 8 && a.size() == 5);
    itv.clear(); a.clear();

    std::cout << "11, ";
    itv.add(1, 100, 7);
    itv.add(10, 110, 8);
    itv.add(15, 16, 9);
    itv.add(30, 130, 11);
    itv.add(50, 60, 4);
    itv.add(100, 200, 11);
    itv.index();
    itv.findOverlaps(20, 90, a); assert (a.front() == 4 && a.size() == 4);
    itv.clear(); a.clear();

    std::cout << "12, ";
    itv.add(1, 10, 1);
    itv.index();
    size_t count = itv.countOverlaps(1, 5);
    assert (count == 1);
    itv.clear(); a.clear();

    std::cout << "13, ";
    itv.index();
    count = itv.countOverlaps(1, 5);
    assert (count == 0);
    itv.clear(); a.clear();

    std::cout << "14, ";
    itv.add(1, 10, 0);
    itv.index();
    itv.searchInterval(5, 11);
    count = 0;
    int last_d{};
    for (const auto& i : itv) {
        last_d = i.data;
        count += 1;
    }
    assert (count == 1); assert (last_d == 0);
    itv.clear(); a.clear();

    std::cout << "15, ";
    itv.add(10, 11, -1);
    itv.add(1, 100, 1);
    itv.add(1, 1000, 2);
    itv.index();
    itv.searchInterval(5, 11);
    count = 0;
    for (const auto& i : itv) {
        last_d = i.data;
        count += 1;
    }
    assert (count == 3);  assert (last_d == 2);
    itv.clear(); a.clear();

    std::cout << "16\n";
    itv.add(1, 10, 0);
    itv.index();
    itv.searchInterval(11, 12);
    count = 0;
    for (const auto& i : itv) {
        last_d = i.data;
        count += 1;
    }
    assert (count == 0);
    itv.clear(); a.clear();

    std::cout << "All tests passed for " << name << "\n";

}


int main(int argc, char *argv[]) {
    auto itv = SuperIntervals<int, int>();
    superTests(itv, "SuperIntervals");
    auto itv2 = SuperIntervalsEytz<int, int>();
    superTests(itv2, "SuperIntervalsEytz");
    return 0;
}
