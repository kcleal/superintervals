
#include "matrylist.hpp"
#include <iostream>
#include <vector>
#include <cassert>


void print_vec(std::vector<size_t>& a, MatryList<int, int>& itv) {
    std::cout << " Found:\n";
    for (auto item : a) {
        std::cout << item << " - " << itv.intervals[item].start << " " << itv.intervals[item].end << std::endl;
    }
}


int main(int argc, char *argv[]) {

    std::vector<size_t> a;
    a.reserve(1000000);

    std::cout << "\n MatryList tests \n";
    MatryList<int, int> itv;
    itv = MatryList<int, int>();

    itv.add(10, 20, -1);
    itv.add(11, 12, -1);
    itv.add(13, 14, -1);
    itv.add(15, 16, -1);
    itv.add(25, 29, -1);
    itv.index();
    itv.findOverlaps(17, 30, a);
    assert (a[0] == 4); assert (a[1] == 0);
    itv.clear();

    itv.add(1, 2, -1);
    itv.add(3, 8, -1);
    itv.add(5, 7, -1);
    itv.add(7, 20, -1);  // 7
    itv.add(9, 10, -1);  // 3
    itv.add(13, 15, -1); // 6
    itv.add(15, 16, -1); // 3
    itv.add(19, 30, -1); // -1
    itv.add(22, 24, -1);
    itv.add(24, 25, -1);
    itv.add(26, 28, -1);
    itv.add(32, 39, -1);
    itv.add(34, 36, -1);
    itv.add(38, 40, -1);
    itv.index();
    itv.findOverlaps(17, 21, a); //print_vec(a, itv);
    assert (a[0] == 7); assert (a[1] == 3);
    itv.clear();
//    return 0;

    itv.add(0, 250000000, -1);
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
    itv.clear();
//    print_vec(a, itv);
    assert (a.back() == 0 && a.size() == 12);
//    return 0;

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
    itv.clear();
    assert (a.back() == 0 && a.size() == 3);

    itv.add(1, 6100000, 6);
    itv.add(4, 5, 6);
    itv.add(6, 7, 7);
    itv.add(9, 10, 7);
    itv.add(11, 12, 7);
    itv.index();
    itv.findOverlaps(2, 25, a);
    itv.clear();
    assert (a.back() == 0 && a.size() == 5);

    itv.add(1, 100, 6);
    itv.add(30, 200, 7);
    itv.add(40, 50, 6);
    itv.add(60, 70, 7);
    itv.index();
    itv.findOverlaps(55, 65, a);
//    int count = itv.countOverlapping(55, 65);
    itv.clear();
//    print_vec(a, itv);
    assert (a.back() == 0 && a.size() == 3);

    itv.add(10, 1001, 0);
    itv.add(30, 400, 1);
    itv.add(60, 700, 2);
    itv.add(65, 80, 3);
    itv.index();
    itv.findOverlaps(100, 200, a);
    itv.clear();
    assert (a.back() == 0 && a.size() == 3);

    itv.add(3, 40, 4);
    itv.add(10, 30, 5);
    itv.add(20, 25, 5);
    itv.add(22, 24, 6);
    itv.index();
    itv.findOverlaps(31, 32, a);
    itv.clear();
    assert (a.back() == 0 && a.size() == 1);
//    print_vec(a, itv);

    itv.add(3, 40, 4);
    itv.add(4, 5, 4);
    itv.add(6, 7, 4);
    itv.add(10, 31, 5);
    itv.add(31, 32, 5);
    itv.index();
    itv.findOverlaps(31, 32, a); assert (a.back() == 0 && a.size() == 3);
    itv.findOverlaps(10, 11, a); assert (a.back() == 0 && a.size() == 2);
    itv.findOverlaps(8, 40, a); assert (a.back() == 0 && a.size() == 3);
    itv.findOverlaps(4, 7, a); assert (a.back() == 0 && a.size() == 3);
    itv.clear();
//    print_vec(a, itv);

    itv.add(3, 40, 4);
    itv.add(3, 40, 4);
    itv.add(3, 40, 4);
    itv.add(3, 4, 4);
    itv.add(35, 50, 4);
    itv.add(40, 400, 4);
    itv.add(40, 400, 4);
    itv.index();
    itv.findOverlaps(38, 41, a); assert (a.back() == 0 && a.size() == 6);
    itv.findOverlaps(41, 42, a); assert (a.back() == 4 && a.size() == 3);
    itv.findOverlaps(339, 410, a); assert (a.back() == 5 && a.size() == 2);
    itv.clear();
//    print_vec(a, itv);

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
    itv.clear();

    itv.add(1, 100, 7);
    itv.add(10, 110, 8);
    itv.add(15, 16, 9);
    itv.add(30, 130, 11);
    itv.add(50, 60, 11);
    itv.add(100, 200, 11);
    itv.index();
    itv.findOverlaps(20, 90, a); assert (a.front() == 4 && a.size() == 4);
    itv.clear();
//    print_vec(a, itv);

    itv.add(1, 100, 7);
    itv.add(20, 120, 8);
    itv.add(110, 210, 9);
    itv.add(130, 140, 10);
    itv.index();
    itv.findOverlaps(90, 190, a); //assert (a.front() == 5 && a.size() == 4);

//    itv.searchInterval(90, 190);
//    auto iter = itv.begin();
//    for (const auto& item : iter) {
//        std::cout << item.start << " " << item.end << " " << item.data << std::endl;
//    }

    std::cout << "All tests passed\n";
    return 0;

}