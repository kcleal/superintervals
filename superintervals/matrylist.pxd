
from libcpp.vector cimport vector

cdef extern from "superintervals.cpp":

    cdef struct IntervalItem:
        int start, end
        int data

    cdef cppclass Iterator:
        Iterator(const MatryList * mlist, size_t index)
        IntervalItem operator *()
        Iterator& operator++()
        bint operator!=(const Iterator& other)
        bint operator==(const Iterator& other)
        Iterator begin()
        Iterator end()
    cdef cppclass MatryList[int, int]:
        MatryList() except +
        void add(int start, int end, int value)
        void searchInterval(int start, int end)
        void clear()
        void reserve(size_t n)
        size_t size()
        size_t countOverlaps(int start, int end)
        void findOverlaps(int start, int end, vector[size_t]& found)

        Iterator begin()
        Iterator end()


ctypedef matrylist.MatryList[int, int] CppMatryList
ctypedef matrylist.MatryList[int, int].Iterator CppIterator
ctypedef matrylist.MatryList[int, int].IntervalItem CppIntervalItem

cdef class matrylist:
    cdef CppMatryList* thisptr
    cdef vector[size_t] found
    cpdef add(self, int start, int end, int value)
    cpdef searchInterval(self, int start, int end)
    cpdef clear(self)
    cpdef reserve(self, size_t n)
    cpdef size(self)
    cpdef countOverlaps(self, int start, int end)
    cpdef findOverlaps(self, int start, int end)
