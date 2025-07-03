
from libcpp.pair cimport pair

__all__ = ["IntervalMap"]


cdef class IntervalMap:
    """
    SuperIntervals interval map to manage a collection of intervals with associated Python objects,
    supporting operations such as adding intervals, checking overlaps, and querying stored data.
    """

    def __cinit__(self):
        """
        Initialize the IntervalMap.
        """
        self.thisptr = new CppIntervalMap[int, PyObjectPtr]()

    def __dealloc__(self):
        cdef PyObjectPtr obj_ptr
        if self.thisptr:
            for obj_ptr in self.thisptr.data:
                if obj_ptr != NULL:
                    Py_DECREF(<object> obj_ptr)
            del self.thisptr

    def __len__(self):
        return self.size()

    def __getitem__(self, int index):
        return self.at(index)

    cpdef add(self, int start, int end, object value=None):
        """
        Add an interval with an associated Python object.

        Args:
            start (int): The start of the interval (inclusive).
            end (int): The end of the interval (inclusive).
            value (object): The Python object to associate with this interval.

        Updates:
            - Adds the interval to the underlying data structure.
            - Stores a reference to the Python object directly in C++.
        """
        cdef PyObjectPtr obj_ptr = NULL
        if value is not None:
            obj_ptr = <PyObjectPtr> value
            Py_INCREF(value)  # Increment reference count
        self.thisptr.add(start, end, obj_ptr)

    cpdef build(self):
        """
        Builds the superintervals index, must be called before queries are made.
        """
        self.thisptr.build()

    cpdef at(self, int index):
        """
        Fetches the interval and data at the given index. Negative indexing is not supported.

        Args:
            index (int): The index of a stored interval.

        Raises:
            IndexError: If the index is out of range.

        Returns:
            tuple: (start, end, data)
        """
        if self.size() == 0 or index < 0 or index >= self.size():
            raise IndexError('Index out of range')
        if self.thisptr.data[index] != NULL:
            return self.thisptr.starts[index], self.thisptr.ends[index], <object> self.thisptr.data[index]
        else:
            return self.thisptr.starts[index], self.thisptr.ends[index], None

    cpdef starts_at(self, int index):
        """
        Fetches the start position at the given index. Negative indexing is not supported.

        Args:
            index (int): The index of a stored interval.

        Raises:
            IndexError: If the index is out of range.

        Returns:
            tuple: start
        """
        if self.size() == 0 or index < 0 or index >= self.size():
            raise IndexError('Index out of range')
        return self.thisptr.starts[index]

    cpdef ends_at(self, int index):
        """
        Fetches the end position at the given index. Negative indexing is not supported.

        Args:
            index (int): The index of a stored interval.

        Raises:
            IndexError: If the index is out of range.

        Returns:
            tuple: start
        """
        if self.size() == 0 or index < 0 or index >= self.size():
            raise IndexError('Index out of range')
        return self.thisptr.ends[index]

    cpdef data_at(self, int index):
        """
        Fetches the stored data at the given index. Negative indexing is not supported.

        Args:
            index (int): The index of a stored interval.

        Raises:
            IndexError: If the index is out of range.

        Returns:
            tuple: start
        """
        if self.size() == 0 or index < 0 or index >= self.size():
            raise IndexError('Index out of range')
        if self.thisptr.data[index] != NULL:
            return None
        else:
            return self.thisptr.starts[index]

    cpdef clear(self):
        """
        Clear all intervals and associated data.
        """
        cdef PyObjectPtr obj_ptr
        for obj_ptr in self.thisptr.data:
            if obj_ptr != NULL:
                Py_DECREF(<object> obj_ptr)
        self.thisptr.clear()

    cpdef reserve(self, size_t n):
        """
        Reserve space for a specified number of intervals.

        Args:
            n (size_t): The number of intervals to reserve space for.
        """
        self.thisptr.reserve(n)

    cpdef size(self):
        """
        Get the number of intervals in the map.

        Returns:
            int: The number of intervals.
        """
        return self.thisptr.size()

    cpdef has_overlaps(self, int start, int end):
        """
        Check if any intervals overlap with a given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (inclusive).

        Returns:
            bool: True if any intervals overlap with the given range, False otherwise.
        """
        return self.thisptr.has_overlaps(start, end)

    cpdef count(self, int start, int end):
        """
        Count the number of intervals that overlap with a given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (inclusive).

        Returns:
            int: The count of overlapping intervals.
        """
        return self.thisptr.count(start, end)

    cpdef search_values(self, int start, int end):
        """
        Find all Python objects associated with intervals that overlap the given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (inclusive).

        Returns:
            list: A list of Python objects associated with overlapping intervals.
        """
        self.found_values.clear()
        self.thisptr.search_values(start, end, self.found_values)
        cdef list result = [None] * self.found_values.size()
        cdef size_t i
        for i in range(self.found_values.size()):
            if self.found_values[i] != NULL:
                result[i] = <object> self.found_values[i]
        return result

    cpdef search_idxs(self, int start, int end):
        """
        Find indices of all intervals that overlap with a given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (inclusive).

        Returns:
            list: A list of indices of overlapping intervals.
        """
        self.found_indexes.clear()
        self.thisptr.search_idxs(start, end, self.found_indexes)
        return list(self.found_indexes)

    cpdef search_keys(self, int start, int end):
        """
        Find interval keys (start, end pairs) that overlap with a given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (inclusive).

        Returns:
            list: A list of (start, end) tuples for overlapping intervals.
        """
        self.found_indexes.clear()
        self.thisptr.search_idxs(start, end, self.found_indexes)
        cdef list result = [None] * self.found_indexes.size()
        cdef size_t i
        for i in range(self.found_indexes.size()):
            result[i] = (self.thisptr.starts[i], self.thisptr.ends[i])
        return result

    cpdef search_items(self, int start, int end):
        """
        Find complete interval items (start, end, data) that overlap with a given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (inclusive).

        Returns:
            list: A list of (start, end, data) tuples for overlapping intervals.
        """
        self.found_indexes.clear()
        self.thisptr.search_idxs(start, end, self.found_indexes)
        cdef list result = [None] * self.found_indexes.size()
        cdef size_t i
        for i in range(self.found_indexes.size()):
            if self.thisptr.data[i] != NULL:
                result[i] = (self.thisptr.starts[i], self.thisptr.ends[i], <object> self.thisptr.data[i])
            else:
                result[i] = (self.thisptr.starts[i], self.thisptr.ends[i], None)
        return result

    cpdef coverage(self, int start, int end):
        """
        Compute coverage statistics for the given range.

        Args:
            start (int): The start of the range (inclusive).
            end (int): The end of the range (inclusive).

        Returns:
            tuple: (count, total_coverage) where count is number of overlapping intervals
                   and total_coverage is the sum of overlapping lengths.
        """
        cdef pair[size_t, int] cov_result = pair[size_t, int](0, 0)
        self.thisptr.coverage(start, end, cov_result)
        return cov_result.first, cov_result.second
