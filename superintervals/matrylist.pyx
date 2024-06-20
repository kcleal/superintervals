
cimport matrylist

from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as inc




cdef class matrylist:
    cdef CppMatryList* thisptr
    cdef vector[size_t] found
    def __cinit__(self):
        self.thisptr = new CppMatryList()

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr

    def add(self, int start, int end, int value):
        self.thisptr.add(start, end, value)

    def searchInterval(self, int start, int end):
        self.thisptr.searchInterval(start, end)

    def clear(self):
        self.thisptr.clear()

    def reserve(self, size_t n):
        self.thisptr.reserve(n)

    def size(self):
        return self.thisptr.size()

    def countOverlaps(self, int start, int end):
        return self.thisptr.countOverlaps(start, end)

    def findOverlaps(self, int start, int end):
        self.found.clear()
        self.thisptr.findOverlaps(start, end, self.found)
        return [self.found[i] for i in range(self.found.size())]

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