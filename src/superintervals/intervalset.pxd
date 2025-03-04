# distutils: language = c++
from libcpp.vector cimport vector

cdef extern from "superintervals.hpp":

    # struct IntervalItem:
    #     int start, end
    #     int data

    cdef cppclass SuperIntervals[int, int]:
        SuperIntervals() except +

        struct Interval:
            int start, end
            int data

        vector[int] starts, ends, data
        size_t idx

        void add(int start, int end, int value)
        void index()
        void at(size_t index, Interval& itv)
        void searchInterval(int start, int end)
        void clear()
        void reserve(size_t n)
        size_t size()
        bint anyOverlaps(int start, int end)
        size_t countOverlaps(int start, int end)
        void findOverlaps(int start, int end, vector[int]& found)
        void findIndexes(int start, int end, vector[size_t]& found_indexes)
        cppclass const_iterator
        cppclass Iterator:
            Iterator(const SuperIntervals * list, size_t index)
            Interval operator *() const
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
    cdef vector[size_t] found_indexes
    cdef list data
    cdef int n_intervals
    cpdef add(self, int start, int end, value=*)
    cpdef add_int_value(self, int start, int end, int value)
    cpdef index(self)
    cpdef at(self, int index)
    cdef void interval_at(self, int index, SuperIntervals.Interval & itv)
    cpdef set_search_interval(self, int start, int end)
    cpdef clear(self)
    cpdef reserve(self, size_t n)
    cpdef size(self)
    cpdef any_overlaps(self, int start, int end)
    cpdef count_overlaps(self, int start, int end)
    cpdef find_overlaps(self, int start, int end)
    cpdef find_indexes(self, int start, int end)
