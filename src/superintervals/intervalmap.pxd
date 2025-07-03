# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.pair cimport pair

cdef extern from "Python.h":
    void Py_INCREF(object)
    void Py_DECREF(object)


cdef extern from "SuperIntervals.hpp" namespace "si":
    cdef cppclass Interval[S, T]:
        S start, end
        T data
        Interval()
        Interval(S s, S e, T d)

    # cdef cppclass IntervalMap[S, T]:
    cdef cppclass CppIntervalMap "si::IntervalMap"[S, T]:
        IntervalMap() except +

        vector[S] starts, ends
        vector[size_t] branch
        vector[T] data
        size_t idx

        void clear()
        void reserve(size_t n)
        size_t size()
        void add(S start, S end, const T& value)
        void build()
        void at(size_t index, Interval[S, T]& itv)
        const Interval[S, T]& at(size_t index) const

        # Search methods
        void upper_bound(const S value) const
        bint has_overlaps(const S start, const S end)
        size_t count_linear(const S start, const S end)
        size_t count(const S start, const S end)
        size_t count_large(const S start, const S end)

        void search_values(const S start, const S end, vector[T]& found)
        void search_values_large(const S start, const S end, vector[T]& found)
        void search_idxs(const S start, const S end, vector[size_t]& found)
        void search_keys(const S start, const S end, vector[pair[S, S]]& found)
        void search_items(const S start, const S end, vector[Interval[S, T]]& found)
        void search_point(const S point, vector[T]& found)
        void coverage(const S start, const S end, pair[size_t, S]& cov_result)

        # Iterator classes not yet implemented!
        # cppclass IndexIterator:
        #     IndexIterator(const IntervalMap * parent, size_t pos)
        #     size_t operator *() const
        #     IndexIterator& operator++()
        #     bint operator !=(const IndexIterator& other) const
        #     bint operator ==(const IndexIterator& other) const
        #
        # cppclass ItemIterator:
        #     ItemIterator(const IntervalMap * parent, size_t pos)
        #     Interval[S, T] operator *() const
        #     ItemIterator& operator++()
        #     bint operator !=(const ItemIterator& other) const
        #     bint operator ==(const ItemIterator& other) const
        #
        # cppclass IndexRange:
        #     IndexRange(const IntervalMap * parent, S start, S end)
        #     IndexIterator begin() const
        #     IndexIterator end() const
        #
        # cppclass ItemRange:
        #     ItemRange(const IntervalMap * parent, S start, S end)
        #     ItemIterator begin() const
        #     ItemIterator end() const
        #
        # IndexRange search_idxs(S start, S end) const
        # ItemRange search_items(S start, S end) const



# Type alias for Python object pointer
ctypedef void * PyObjectPtr

cdef class IntervalMap:
    cdef CppIntervalMap[int, PyObjectPtr] * thisptr
    cdef vector[PyObjectPtr] found_values
    cdef vector[size_t] found_indexes

    cpdef add(self, int start, int end, object value= *)
    cpdef build(self)
    cpdef at(self, int index)
    cpdef starts_at(self, int index)
    cpdef ends_at(self, int index)
    cpdef data_at(self, int index)
    cpdef clear(self)
    cpdef reserve(self, size_t n)
    cpdef size(self)
    cpdef has_overlaps(self, int start, int end)
    cpdef count(self, int start, int end)
    cpdef search_values(self, int start, int end)
    cpdef search_idxs(self, int start, int end)
    cpdef search_keys(self, int start, int end)
    cpdef search_items(self, int start, int end)
    cpdef coverage(self, int start, int end)
