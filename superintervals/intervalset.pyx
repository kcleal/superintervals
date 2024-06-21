

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

    cpdef count_overlaps(self, int start, int end):
        return self.thisptr.countOverlaps(start, end)

    cpdef find_overlaps(self, int start, int end):
        self.found.clear()
        self.thisptr.findOverlaps(start, end, self.found)
        # return [(self.thisptr.starts[i], self.thisptr.ends[i]) for i in self.found]
        return self.found
        # return [self.found[i] for i in range(self.found.size())]

    # def __iter__(self):
    #     return IteratorWrapper(self)

# cdef class IteratorWrapper:
#     cdef CppIterator* iter_ptr
#     cdef MatryList parent
#
#     def __cinit__(self, MatryList parent):
#         self.parent = parent
#         self.iter_ptr = new CppIterator(parent.thisptr.begin())
#
#     def __dealloc__(self):
#         if self.iter_ptr is not NULL:
#             del self.iter_ptr
#
#     def __next__(self):
#         cdef CppIntervalItem item = deref(self.iter_ptr)
#         if self.iter_ptr.operator!=(self.parent.thisptr.end()):
#             inc(self.iter_ptr)
#             return (item.start, item.end, item.data)
#         else:
#             raise StopIteration()