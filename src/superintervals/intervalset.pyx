
from cython.operator cimport dereference, postincrement, preincrement


__all__ = ["IntervalSet"]

cdef class IntervalSet:
    def __cinit__(self, with_data=False):
        self.thisptr = new SuperIntervals()
        self.n_intervals = 0
        self.with_data = with_data
    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    cpdef add(self, int start, int end, value=None):
        self.thisptr.add(start, end, self.n_intervals)
        if self.with_data:
            self.data.append(value)
        self.n_intervals += 1

    cpdef index(self):
        if self.with_data:
            assert self.data == self.n_intervals
        self.thisptr.index()

    cpdef set_search_interval(self, int start, int end):
        self.thisptr.searchInterval(start, end)

    cpdef clear(self):
        self.thisptr.clear()
        self.n_intervals = 0

    cpdef reserve(self, size_t n):
        self.thisptr.reserve(n)

    cpdef size(self):
        return self.thisptr.size()

    cpdef any_overlaps(self, int start, int end):
        return self.thisptr.anyOverlaps(start, end)

    cpdef count_overlaps(self, int start, int end):
        return self.thisptr.countOverlaps(start, end)

    cpdef find_overlaps(self, int start, int end):
        self.found.clear()
        self.thisptr.findOverlaps(start, end, self.found)
        return self.found

    def __iter__(self):
        return IteratorWrapper(self)


cdef class IteratorWrapper:
    cdef SuperIntervals.Iterator * _cpp_iterator
    cdef SuperIntervals * _si

    def __cinit__(self, IntervalSet interval_set):
        self._si = interval_set.thisptr
        self._cpp_iterator = new SuperIntervals.Iterator.Iterator(interval_set.thisptr, interval_set.thisptr.idx)

    def __dealloc__(self):
        del self._cpp_iterator

    def __iter__(self):
        return self

    def __next__(self):
        if self._cpp_iterator[0] == self._cpp_iterator[0].end():
            raise StopIteration

        cdef int start = self._si.starts[self._cpp_iterator.it_index]
        cdef int end = self._si.ends[self._cpp_iterator.it_index]
        cdef int data = self._si.data[self._cpp_iterator.it_index]

        preincrement(self._cpp_iterator[0])
        return start, end, data
