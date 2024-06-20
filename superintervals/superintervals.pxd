# distutils: language = c++
from libcpp.vector cimport vector

cdef extern from "superintervals.hpp":

    # cdef struct IntervalItem:
    #     int start, end
    #     int data

    # cdef cppclass Iterator:
    #     Iterator(const SuperIntervals[int, int] *mlist, size_t index)
    #     IntervalItem operator*()
    #     Iterator& operator++()
    #     bint operator!=(const Iterator& other)
    #     bint operator==(const Iterator& other)
    #     Iterator begin()
    #     Iterator end()

    cdef cppclass SuperIntervals[int, int]:
        SuperIntervals() except +
        void add(int start, int end, int value)
        void searchInterval(int start, int end)
        void clear()
        void reserve(size_t n)
        size_t size()
        size_t countOverlaps(int start, int end)
        void findOverlaps(int start, int end, vector[int]& found)

        # Iterator begin()
        # Iterator end()

ctypedef SuperIntervals[int, int] SuperIntervalsCpp
# ctypedef Iterator CppIterator
# ctypedef IntervalItem CppIntervalItem

cdef class IntervalSet:
    cdef SuperIntervalsCpp* thisptr
    cdef vector[int] found
    cpdef add(self, int start, int end, int value)
    cpdef searchInterval(self, int start, int end)
    cpdef clear(self)
    cpdef reserve(self, size_t n)
    cpdef size(self)
    cpdef countOverlaps(self, int start, int end)
    cpdef findOverlaps(self, int start, int end)