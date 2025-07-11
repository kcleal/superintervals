
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

    @classmethod
    def from_arrays(cls, starts, ends, values=None):
        """
        Create an IntervalMap from arrays of starts, ends, and optional values.

        This is the most efficient way to create an IntervalMap when you have
        existing data in array format. The resulting IntervalMap is ready to use
        (no need to call build()).

        Args:
            starts: Array-like of start positions (array.array, numpy array, list, etc.)
            ends: Array-like of end positions (array.array, numpy array, list, etc.)
            values: Optional iterable of values to associate with each interval.
                    If None, no values are stored.

        Returns:
            IntervalMap: A new IntervalMap ready for queries

        Examples:
            >>> from array import array
            >>> import numpy as np
            >>>
            >>> # Using array.array
            >>> starts = array('i', [1, 10, 20])
            >>> ends = array('i', [5, 15, 25])
            >>> values = ["gene1", "gene2", "gene3"]
            >>> im = IntervalMap.from_arrays(starts, ends, values)
            >>>
            >>> # Using numpy arrays
            >>> starts = np.array([1, 10, 20], dtype=np.int32)
            >>> ends = np.array([5, 15, 25], dtype=np.int32)
            >>> im = IntervalMap.from_arrays(starts, ends)
            >>>
            >>> # Using lists (will be converted internally)
            >>> im = IntervalMap.from_arrays([1, 10, 20], [5, 15, 25], ["A", "B", "C"])
        """
        cdef IntervalMap instance = cls()
        cdef int[:] start_view
        cdef int[:] end_view
        if hasattr(starts, 'shape'):  # numpy array or already array-like
            start_view = starts
            end_view = ends
        else:  # lists or other iterables
            from array import array
            start_array = array('i', starts)
            end_array = array('i', ends)
            start_view = start_array
            end_view = end_array

        if start_view.shape[0] != end_view.shape[0]:
            raise ValueError("starts and ends must have the same length")

        cdef size_t n = start_view.shape[0]
        cdef size_t i
        cdef object value

        instance.reserve(n)

        if values is None:
            for i in range(n):
                instance.add(start_view[i], end_view[i])
        else:
            if len(values) != n:
                raise ValueError("values length must match starts/ends length")
            for i, value in enumerate(values):
                instance.add(start_view[i], end_view[i], value)

        instance.build()
        return instance

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
        if self.thisptr.data[index] == NULL:
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

    cpdef count_batch(self, int[:] starts, int[:] ends):
        """
        Count overlapping intervals for multiple query ranges.

        Args:
            starts: Memory view of start positions (array.array, numpy array, etc.)
            ends: Memory view of end positions (array.array, numpy array, etc.)

        Returns:
            List: Count of overlapping intervals for each query range

        Example:
            >>> from array import array
            >>> import numpy as np
            >>> im = IntervalMap()
            >>> im.add(1, 10, "A")
            >>> im.add(5, 15, "B")
            >>> im.build()
            >>> 
            >>> # Works with array.array
            >>> starts = array('i', [1, 8, 15])
            >>> ends = array('i', [5, 12, 20])
            >>> counts = im.count_batch(starts, ends)
            >>> 
            >>> # Also works with numpy arrays
            >>> starts_np = np.array([1, 8, 15], dtype=np.int32)
            >>> ends_np = np.array([5, 12, 20], dtype=np.int32)
            >>> counts = im.count_batch(starts_np, ends_np)
        """
        if starts.shape[0] != ends.shape[0]:
            raise ValueError("starts and ends must have the same length")
        cdef size_t n = starts.shape[0]
        cdef list result_array = [0] * n
        cdef size_t i
        for i in range(n):
            result_array[i] = self.thisptr.count(starts[i], ends[i])

        return result_array

    cpdef search_idxs_batch(self, int[:] starts, int[:] ends):
        """
        Find indices of overlapping intervals for multiple query ranges.

        Args:
            starts: Memory view of start positions 
            ends: Memory view of end positions

        Returns:
            list: List of lists, where each sublist contains indices of 
                  overlapping intervals for the corresponding query range

        Example:
            >>> from array import array
            >>> import numpy as np
            >>> im = IntervalMap()
            >>> im.add(1, 10, "A")
            >>> im.add(5, 15, "B")
            >>> im.build()
            >>> 
            >>> # Works with array.array
            >>> starts = array('i', [1, 8])
            >>> ends = array('i', [5, 12])
            >>> results = im.search_idxs_batch(starts, ends)
            >>> 
            >>> # Works with numpy arrays
            >>> starts_np = np.array([1, 8], dtype=np.int32)
            >>> ends_np = np.array([5, 12], dtype=np.int32)
            >>> results = im.search_idxs_batch(starts_np, ends_np)
        """
        if starts.shape[0] != ends.shape[0]:
            raise ValueError("starts and ends must have the same length")

        cdef size_t n = starts.shape[0]
        cdef list results = [[] for _ in range(n)]
        cdef size_t i, j
        for i in range(n):
            self.found_indexes.clear()
            self.thisptr.search_idxs(starts[i], ends[i], self.found_indexes)
            query_result = [0] * self.found_indexes.size()
            for j in range(self.found_indexes.size()):
                query_result[j] = self.found_indexes[j]
            results[i] = query_result

        return results

    cpdef search_values_batch(self, int[:] starts, int[:] ends):
        """
        Find values of overlapping intervals for multiple query ranges.

        Args:
            starts: Memory view of start positions
            ends: Memory view of end positions

        Returns:
            list: List of lists, where each inner list contains values of 
                  overlapping intervals for the corresponding query range

        Example:
            >>> from array import array
            >>> import numpy as np
            >>> im = IntervalMap()
            >>> im.add(1, 10, "A")
            >>> im.add(5, 15, "B")
            >>> im.build()
            >>> 
            >>> # Works with array.array
            >>> starts = array('i', [1, 8])
            >>> ends = array('i', [5, 12])
            >>> results = im.search_values_batch(starts, ends)
            >>> 
            >>> # Works with numpy arrays  
            >>> starts_np = np.array([1, 8], dtype=np.int32)
            >>> ends_np = np.array([5, 12], dtype=np.int32)
            >>> results = im.search_values_batch(starts_np, ends_np)
        """
        if starts.shape[0] != ends.shape[0]:
            raise ValueError("starts and ends must have the same length")

        cdef size_t n = starts.shape[0]
        cdef list results = [[] for _ in range(n)]
        cdef size_t i, j
        cdef list query_result
        for i in range(n):
            self.found_values.clear()
            self.thisptr.search_values(starts[i], ends[i], self.found_values)
            query_result = [None] * self.found_values.size()
            for j in range(self.found_values.size()):
                if self.found_values[j] != NULL:
                    query_result[j] = <object> self.found_values[j]
            results[i] = query_result

        return results
