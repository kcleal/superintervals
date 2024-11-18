# distutils: language = c++
from libcpp.vector cimport vector

cdef extern from "superintervals.hpp":

    struct IntervalItem:
        int start, end
        int data

    cdef cppclass SuperIntervals[int, int]:
        SuperIntervals() except +

        vector[int] starts, ends, data
        size_t idx

        void add(int start, int end, int value)
        void index()
        void searchInterval(int start, int end)
        void clear()
        void reserve(size_t n)
        size_t size()
        bint anyOverlaps(int start, int end)
        size_t countOverlaps(int start, int end)
        void findOverlaps(int start, int end, vector[int]& found)

        cppclass const_iterator
        cppclass Iterator:
            Iterator(const SuperIntervals * list, size_t index)
            IntervalItem operator *() const
            Iterator& operator++()
            bint operator !=(const Iterator& other) const
            bint operator ==(const Iterator& other) const
            Iterator begin() const
            Iterator end() const

            size_t it_index

        Iterator begin() const
        Iterator end() const


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
    cpdef any_overlaps(self, int start, int end)
    cpdef count_overlaps(self, int start, int end)
    cpdef find_overlaps(self, int start, int end)
