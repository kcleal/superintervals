
from cython.operator cimport dereference, postincrement, preincrement


__all__ = ["IntervalSet"]

cdef class IntervalSet:
    """
    Superintervals interval set to manage a collection of intervals, supporting operations such as adding intervals,
    checking overlaps, and querying stored data.
    """
    def __cinit__(self):
        self.thisptr = new SuperIntervals()
        self.n_intervals = 0
        self.data = []
    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    cpdef add(self, int start, int end, value=None):
        """
        Add an interval with an optional associated value.

        Args:
            start (int): The start of the interval (inclusive).
            end (int): The end of the interval (exclusive).
            value (optional): The value to associate with this interval. Defaults to None.

        Updates:
            - Adds the interval to the underlying data structure.
            - Stores the value in the `data` list.
            - Increments `n_intervals`.
        """
        self.thisptr.add(start, end, self.n_intervals)
        self.data.append(value)
        self.n_intervals += 1

    cpdef add_int_value(self, int start, int end, int value):
        """
        Add an interval with an integer value directly.

        Args:
            start (int): The start of the interval (inclusive).
            end (int): The end of the interval (exclusive).
            value (int): The integer value to associate with this interval.

        Updates:
            - Adds the interval and the value to the underlying data structure.
            - Increments `n_intervals`.
        """
        self.thisptr.add(start, end, value)
        self.n_intervals += 1

    cpdef index(self):
        """
        Build an index for the intervals, enabling efficient queries.

        Raises:
            ValueError: If the `data` list length does not match `n_intervals`, 
                        indicating misuse of `add` and `add_int_value`.
        """
        if len(self.data):
            if len(self.data) != self.n_intervals:
                raise ValueError('Value data list != n_intervals. Use only add or add_int_value functions, not both')
        self.thisptr.index()

    cpdef at(self, int index):
        """
        Fetches the interval and data at the given index. Negative indexing is not supported.
        Args:
            index (int): The index of a stored interval.
        Raises:
            IndexError: If the index is out of range.
        Returns:
            tuple: start, end, data
        """
        if self.size() == 0 or index < 0 or index > self.size():
            return IndexError('Index out of range')
        cdef SuperIntervals.Interval itv
        self.thisptr.at(index, itv)
        if len(self.data) == 0:
            return itv.start, itv.end, itv.data
        else:
            return itv.start, itv.end, self.data[itv.data]

    cdef void interval_at(self, int index, SuperIntervals.Interval &itv):
        self.thisptr.at(index, itv)

    cpdef set_search_interval(self, int start, int end):
        """
        Define a search interval for querying overlaps. This is only needed if using the IntervalSet as an iterator.

        Args:
            start (int): The start of the search interval (inclusive).
            end (int): The end of the search interval (exclusive).
        """
        self.thisptr.searchInterval(start, end)

    cpdef clear(self):
        """
        Clear all intervals and associated data.

        Updates:
            - Removes all intervals from the underlying data structure.
            - Resets `n_intervals` to 0.
        """
        self.thisptr.clear()
        self.n_intervals = 0
        self.data = []

    cpdef reserve(self, size_t n):
        """
        Reserve space for a specified number of intervals in the underlying data structure.

        Args:
            n (size_t): The number of intervals to reserve space for.
        """
        self.thisptr.reserve(n)

    cpdef size(self):
        """
        Get the number of intervals in the set.

        Returns:
            int: The number of intervals.
        """
        return self.thisptr.size()

    cpdef any_overlaps(self, int start, int end):
        """
        Check if any intervals overlap with a given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (exclusive).

        Returns:
            bool: True if any intervals overlap with the given range, False otherwise.
        """
        return self.thisptr.anyOverlaps(start, end)

    cpdef count_overlaps(self, int start, int end):
        """
        Count the number of intervals that overlap with a given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (exclusive).

        Returns:
            int: The count of overlapping intervals.
        """
        return self.thisptr.countOverlaps(start, end)

    cpdef find_overlaps(self, int start, int end):
        """
        Find all intervals that overlap with a given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (exclusive).

        Returns:
            list: A list of values associated with the overlapping intervals if `data` is not empty.
                  Otherwise, returns a list of indices of the overlapping intervals.
        """
        self.found.clear()
        self.thisptr.findOverlaps(start, end, self.found)
        cdef int i
        if len(self.data):
            return [self.data[i] for i in self.found]
        return self.found

    cpdef find_indexes(self, int start, int end):
        """
        Find all interval indexes that overlap with a given range. Use 'at()' on an indexes to get interval data

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (exclusive).

        Returns:
            list: Returns a list of indices of the overlapping intervals.
        """
        self.found_indexes.clear()
        self.thisptr.findIndexes(start, end, self.found_indexes)
        return self.found_indexes

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
