
from libcpp.pair cimport pair

__all__ = ["IntervalMap"]


# Module-level helpers (Cython does not support closures inside cpdef methods).
def _start_end_key(t):
    # Sort (start, end, value) tuples by coordinates only; values may not be
    # comparable to one another.
    return (t[0], t[1])


def _keep_first(a, b):
    return a


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
            return <object> self.thisptr.data[index]

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
        cdef size_t idx
        for i in range(self.found_indexes.size()):
            idx = self.found_indexes[i]
            result[i] = (self.thisptr.starts[idx], self.thisptr.ends[idx])
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
        cdef size_t idx
        for i in range(self.found_indexes.size()):
            idx = self.found_indexes[i]
            if self.thisptr.data[idx] != NULL:
                result[i] = (self.thisptr.starts[idx], self.thisptr.ends[idx], <object> self.thisptr.data[idx])
            else:
                result[i] = (self.thisptr.starts[idx], self.thisptr.ends[idx], None)
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

    # -------------------------------------------------------------------------
    # Set-like operations
    #
    # Each operation returns a NEW IntervalMap that has ALREADY been built and is
    # ready to query. Intervals are end-inclusive. Two intervals overlap when
    # a_start <= b_end and b_start <= a_end. Merging coalesces overlapping
    # intervals only; exactly adjacent intervals like (1, 5) and (6, 10) are not
    # merged.
    #
    # Data handling: when intervals merge, their Python values are combined into
    # a tuple by default. Pass a `combine` callable taking (a, b) and returning a
    # value to override this (e.g. combine=lambda a, b: a to keep the first).
    # -------------------------------------------------------------------------

    cdef list _sorted_intervals(self):
        # Stored intervals as (start, end, value) tuples, sorted by (start, end).
        cdef size_t n = self.thisptr.size()
        cdef list items = [None] * n
        cdef size_t i
        cdef object v
        for i in range(n):
            if self.thisptr.data[i] != NULL:
                v = <object> self.thisptr.data[i]
            else:
                v = None
            items[i] = (self.thisptr.starts[i], self.thisptr.ends[i], v)
        items.sort(key=_start_end_key)
        return items

    cpdef merge_overlaps(self, combine=None):
        """
        Coalesce overlapping intervals into a disjoint, non-overlapping set.

        Args:
            combine: Optional callable (a, b) -> value to combine the data of
                     intervals that merge together. By default the values are
                     collected into a tuple.

        Returns:
            IntervalMap: A new, built IntervalMap of the merged intervals.
        """
        cdef IntervalMap out = IntervalMap()
        cdef list items = self._sorted_intervals()
        if not items:
            out.build()
            return out
        cur_start, cur_end, cur_data = items[0]
        cdef Py_ssize_t k
        for k in range(1, len(items)):
            s, e, d = items[k]
            if s <= cur_end:  # overlap -> extend
                if e > cur_end:
                    cur_end = e
                cur_data = combine(cur_data, d) if combine is not None else (
                    cur_data + (d,) if isinstance(cur_data, tuple) else (cur_data, d))
            else:  # disjoint -> flush
                out.add(cur_start, cur_end, cur_data)
                cur_start, cur_end, cur_data = s, e, d
        out.add(cur_start, cur_end, cur_data)
        out.build()
        return out

    cpdef union_with(self, IntervalMap other, combine=None):
        """
        Union of this set with another: all covered regions, coalesced.

        Args:
            other (IntervalMap): The other interval map.
            combine: Optional callable (a, b) -> value used when intervals merge.
                     By default merged values are collected into a tuple.

        Returns:
            IntervalMap: A new, built IntervalMap of the merged union.
        """
        cdef IntervalMap combined = IntervalMap()
        cdef list items = self._sorted_intervals() + other._sorted_intervals()
        items.sort(key=_start_end_key)
        for s, e, d in items:
            combined.add(s, e, d)
        return combined.merge_overlaps(combine)

    cpdef intersection(self, IntervalMap other, combine=None):
        """
        Intersection of this set with another: regions covered by both.

        Args:
            other (IntervalMap): The other interval map (must be built).
            combine: Optional callable (a, b) -> value for each overlapping pair,
                     where a is this side's value and b is the other's. By default
                     the two values are paired into a tuple (a, b).

        Returns:
            IntervalMap: A new, built IntervalMap of every overlapping sub-region.
                         The result is not coalesced; call merge_overlaps() on it
                         if a disjoint set is required.
        """
        cdef IntervalMap out = IntervalMap()
        cdef size_t n = self.thisptr.size()
        cdef size_t i, j, h
        cdef int s, e
        cdef object a_val, b_val
        for i in range(n):
            self.found_indexes.clear()
            other.thisptr.search_idxs(self.thisptr.starts[i], self.thisptr.ends[i], other.found_indexes)
            for j in range(other.found_indexes.size()):
                h = other.found_indexes[j]
                s = max(self.thisptr.starts[i], other.thisptr.starts[h])
                e = min(self.thisptr.ends[i], other.thisptr.ends[h])
                if s <= e:
                    a_val = <object> self.thisptr.data[i] if self.thisptr.data[i] != NULL else None
                    b_val = <object> other.thisptr.data[h] if other.thisptr.data[h] != NULL else None
                    out.add(s, e, combine(a_val, b_val) if combine is not None else (a_val, b_val))
        out.build()
        return out

    cpdef difference(self, IntervalMap other):
        """
        Difference: regions of this set not covered by `other` (A minus B).

        Args:
            other (IntervalMap): The other interval map (must be built).

        Returns:
            IntervalMap: A new, built IntervalMap. Output pieces are sub-ranges of
                         this set's intervals and inherit their source value.
        """
        cdef IntervalMap out = IntervalMap()
        cdef size_t n = self.thisptr.size()
        cdef size_t i
        cdef int cursor, cs, ce
        cdef object val
        cdef list cover
        cdef size_t j, h
        for i in range(n):
            other.found_indexes.clear()
            other.thisptr.search_idxs(self.thisptr.starts[i], self.thisptr.ends[i], other.found_indexes)
            cover = []
            for j in range(other.found_indexes.size()):
                h = other.found_indexes[j]
                cover.append((other.thisptr.starts[h], other.thisptr.ends[h]))
            cover.sort()
            val = <object> self.thisptr.data[i] if self.thisptr.data[i] != NULL else None
            cursor = self.thisptr.starts[i]
            for c0, c1 in cover:
                cs = max(c0, self.thisptr.starts[i])
                ce = min(c1, self.thisptr.ends[i])
                if cs > cursor:
                    out.add(cursor, cs - 1, val)
                if ce + 1 > cursor:
                    cursor = ce + 1
            if cursor <= self.thisptr.ends[i]:
                out.add(cursor, self.thisptr.ends[i], val)
        out.build()
        return out

    cpdef symmetric_difference(self, IntervalMap other):
        """
        Symmetric difference: regions in exactly one of the two sets,
        (A minus B) union (B minus A). Both inputs must be built.

        Args:
            other (IntervalMap): The other interval map (must be built).

        Returns:
            IntervalMap: A new, built IntervalMap.
        """
        cdef IntervalMap a_minus_b = self.difference(other)
        cdef IntervalMap b_minus_a = other.difference(self)
        return a_minus_b.union_with(b_minus_a)

    cpdef gaps(self, int lo, int hi, object fill=None):
        """
        Compute the gaps (uncovered regions) within [lo, hi].

        Args:
            lo (int): Lower bound of the span (inclusive).
            hi (int): Upper bound of the span (inclusive).
            fill: Value assigned to gap intervals (which match no input). Default None.

        Returns:
            IntervalMap: A new, built IntervalMap of the complement within [lo, hi].
        """
        cdef IntervalMap out = IntervalMap()
        cdef IntervalMap merged = self.merge_overlaps(_keep_first)  # geometry only
        cdef size_t k
        cdef int cursor = lo
        cdef int s, e
        for k in range(merged.thisptr.size()):
            s = merged.thisptr.starts[k]
            e = merged.thisptr.ends[k]
            if e < lo or s > hi:
                continue
            if s > cursor:
                out.add(cursor, s - 1, fill)
            if e + 1 > cursor:
                cursor = e + 1
        if cursor <= hi:
            out.add(cursor, hi, fill)
        out.build()
        return out

    cpdef span(self):
        """
        Smallest interval enclosing every stored interval (convex hull).

        Returns:
            tuple or None: (min_start, max_end), or None if the map is empty.
        """
        cdef size_t n = self.thisptr.size()
        if n == 0:
            return None
        cdef int lo = self.thisptr.starts[0]
        cdef int hi = self.thisptr.ends[0]
        cdef size_t k
        for k in range(1, n):
            if self.thisptr.starts[k] < lo:
                lo = self.thisptr.starts[k]
            if self.thisptr.ends[k] > hi:
                hi = self.thisptr.ends[k]
        return (lo, hi)

    cpdef expand(self, int left, int right, lo=None, hi=None):
        """
        Grow (or shrink) every interval, like ``bedtools slop``.

        Args:
            left (int): Subtracted from each start (extends left); negative shrinks.
            right (int): Added to each end (extends right); negative shrinks.
            lo: Optional lower clamp; no start falls below this. Default unbounded.
            hi: Optional upper clamp; no end rises above this. Default unbounded.

        Returns:
            IntervalMap: A new, built IntervalMap of the resized intervals. Data is
                         carried over unchanged; intervals that shrink past themselves
                         (start > end) are dropped.
        """
        cdef IntervalMap out = IntervalMap()
        cdef size_t n = self.thisptr.size()
        cdef size_t i
        cdef object s, e, val
        for i in range(n):
            s = self.thisptr.starts[i] - left
            e = self.thisptr.ends[i] + right
            if lo is not None and s < lo:
                s = lo
            if hi is not None and e > hi:
                e = hi
            if s <= e:
                val = <object> self.thisptr.data[i] if self.thisptr.data[i] != NULL else None
                out.add(s, e, val)
        out.build()
        return out

    cpdef flank(self, int left, int right, lo=None, hi=None):
        """
        Create flanking intervals beside each interval, like ``bedtools flank``.

        Emits the left flank ``[start-left, start-1]`` and right flank
        ``[end+1, end+right]`` (NOT the originals). Zero-width or out-of-range
        flanks are skipped; each flank inherits its source interval's value.

        Args:
            left (int): Width of the left flank (0 = none).
            right (int): Width of the right flank (0 = none).
            lo: Optional lower clamp. Default unbounded.
            hi: Optional upper clamp. Default unbounded.

        Returns:
            IntervalMap: A new, built IntervalMap holding the flank intervals only.
        """
        cdef IntervalMap out = IntervalMap()
        cdef size_t n = self.thisptr.size()
        cdef size_t i
        cdef object ls, le, rs, re, val
        for i in range(n):
            val = <object> self.thisptr.data[i] if self.thisptr.data[i] != NULL else None
            if left > 0 and (lo is None or self.thisptr.starts[i] > lo):
                le = self.thisptr.starts[i] - 1
                ls = self.thisptr.starts[i] - left
                if lo is not None and ls < lo:
                    ls = lo
                if ls <= le and (hi is None or le <= hi):
                    out.add(ls, le, val)
            if right > 0 and (hi is None or self.thisptr.ends[i] < hi):
                rs = self.thisptr.ends[i] + 1
                re = self.thisptr.ends[i] + right
                if hi is not None and re > hi:
                    re = hi
                if rs <= re and (lo is None or rs >= lo):
                    out.add(rs, re, val)
        out.build()
        return out

    cpdef unique(self, combine=None):
        """
        Keep only one interval per distinct (start, end) coordinate pair.

        Unlike merge_overlaps(), this collapses only EXACT coordinate duplicates;
        distinct-but-overlapping intervals are all preserved.

        Args:
            combine: Optional callable (a, b) -> value to fold the data of exact
                     duplicates. By default the first interval's value is kept.

        Returns:
            IntervalMap: A new, built IntervalMap with duplicates removed.
        """
        cdef IntervalMap out = IntervalMap()
        cdef list items = self._sorted_intervals()  # (start, end) ascending
        if not items:
            out.build()
            return out
        run_start, run_end, acc = items[0]
        cdef Py_ssize_t k
        for k in range(1, len(items)):
            s, e, d = items[k]
            if s == run_start and e == run_end:
                acc = combine(acc, d) if combine is not None else acc  # fold duplicates
            else:
                out.add(run_start, run_end, acc)
                run_start, run_end, acc = s, e, d
        out.add(run_start, run_end, acc)
        out.build()
        return out
