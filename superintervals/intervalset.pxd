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

        vector[int] starts, ends, data
        void add(int start, int end, int value)
        void index()
        void searchInterval(int start, int end)
        void clear()
        void reserve(size_t n)
        size_t size()
        size_t countOverlaps(int start, int end)
        void findOverlaps(int start, int end, vector[int]& found)

        # Iterator begin()
        # Iterator end()

# ctypedef Iterator CppIterator
# ctypedef IntervalItem CppIntervalItem

cdef class IntervalSet:
    cdef SuperIntervals* thisptr
    cdef vector[int] found
    cdef bint with_data
    cdef list data
    cdef int n_intervals
    cpdef add(self, int start, int end, value=*)
    cpdef index(self)
    cpdef set_search_interval(self, int start, int end)
    cpdef clear(self)
    cpdef reserve(self, size_t n)
    cpdef size(self)
    cpdef count_overlaps(self, int start, int end)
    cpdef find_overlaps(self, int start, int end)