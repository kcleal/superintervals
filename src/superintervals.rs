
//! This module provides efficient data structures for managing and querying intervals.
//! It includes implementations for standard and Eytzinger layout-based interval storage.

use std::cmp::Ordering;

/// Represents an interval with associated data.
#[derive(Debug, Clone)]
pub struct Interval<T> {
    start: i32,
    end: i32,
    data: T,
}


/// A static data structure for finding interval intersections
///
/// SuperIntervals is a template class that provides efficient interval intersection operations.
/// It supports adding intervals, indexing them for fast queries, and performing various
/// intersection operations.
///
/// Intervals are considered end-inclusive
/// The index() function must be called before any queries. If more intervals are added, call index() again.
pub struct SuperIntervals<T> {
    starts: Vec<i32>,
    ends: Vec<i32>,
    branch: Vec<usize>,
    data: Vec<T>,
    idx: usize,
    start_sorted: bool,
    end_sorted: bool,
}

impl<T> SuperIntervals<T>
where
    T: Clone,
{
    pub fn new() -> Self {
        SuperIntervals {
            starts: Vec::new(),
            ends: Vec::new(),
            branch: Vec::new(),
            data: Vec::new(),
            idx: 0,
            start_sorted: true,
            end_sorted: true,
        }
    }
    /// Clears all intervals from the structure.
    pub fn clear(&mut self) {
        self.starts.clear();
        self.ends.clear();
        self.branch.clear();
        self.data.clear();
        self.idx = 0;
    }
    /// Adds a new interval to the structure.
    ///
    /// # Arguments
    ///
    /// * `start` - The start point of the interval.
    /// * `end` - The end point of the interval.
    /// * `value` - The data associated with the interval.
    pub fn add(&mut self, start: i32, end: i32, value: T) {
        if self.start_sorted && !self.starts.is_empty() {
            self.start_sorted = start >= *self.starts.last().unwrap();
            if self.start_sorted && start == *self.starts.last().unwrap() && end > *self.ends.last().unwrap() {
                self.end_sorted = false;
            }
        }
        self.starts.push(start);
        self.ends.push(end);
        self.data.push(value);
    }
    /// Sorts the intervals by start and -end.
    pub fn sort_intervals(&mut self) {
        if !self.start_sorted {
            self.sort_block(0, self.starts.len(), |a, b| {
                if a.start < b.start {
                    Ordering::Less
                } else if a.start == b.start {
                    if a.end > b.end {
                        Ordering::Less
                    } else {
                        Ordering::Greater
                    }
                } else {
                    Ordering::Greater
                }
            });
            self.start_sorted = true;
            self.end_sorted = true;
        } else if !self.end_sorted {
            let mut it_start = 0;
            while it_start < self.starts.len() {
                let mut block_end = it_start + 1;
                let mut needs_sort = false;
                while block_end < self.starts.len() && self.starts[block_end] == self.starts[it_start] {
                    if block_end > it_start && self.ends[block_end] > self.ends[block_end - 1] {
                        needs_sort = true;
                    }
                    block_end += 1;
                }
                if needs_sort {
                    self.sort_block(it_start, block_end, |a, b| a.end.cmp(&b.end).reverse());
                }
                it_start = block_end;
            }
            self.end_sorted = true;
        }
    }
    /// Indexes the intervals. Must be called before queries are performed.
    pub fn index(&mut self) {
        if self.starts.is_empty() {
            return;
        }
        self.starts.shrink_to_fit();
        self.ends.shrink_to_fit();
        self.data.shrink_to_fit();
        self.sort_intervals();
        self.branch.resize(self.starts.len(), usize::MAX);
        let mut br: Vec<(i32, usize)> = Vec::with_capacity(1000);
        unsafe {
            br.push((*self.ends.get_unchecked(0), 0));
            for i in 1..self.ends.len() {
                while !br.is_empty() && br.last().unwrap().0 < *self.ends.get_unchecked(i) {
                    br.pop();
                }
                if !br.is_empty() {
                    *self.branch.get_unchecked_mut(i) = br.last().unwrap().1;
                }
                br.push((*self.ends.get_unchecked(i), i));
            }
        }
        self.idx = 0;
    }

    pub fn upper_bound(&mut self, value: i32) {
        let mut length = self.starts.len();
        if length == 0 {
            return;
        }
        length -= 1;
        self.idx = 0;
        while length >= 196 {
            let half = length / 2;
            unsafe {
                let _ = self.starts.get_unchecked(self.idx + half / 2);
                let first_half1 = self.idx + (length - half);
                let _ = self.starts.get_unchecked(first_half1 + half / 2);
                if *self.starts.get_unchecked(self.idx + half) <= value {
                    self.idx += length - half;
                }
            }
            length = half;
        }
        while length > 0 {
            let half = length / 2;
            unsafe {
                if *self.starts.get_unchecked(self.idx + half) <= value {
                    self.idx += length - half;
                }
            }
            length = half;
        }
        if self.idx > 0 && (self.idx == self.starts.len() || self.starts[self.idx] > value) {
            self.idx -= 1;
        }
    }
    /// Finds all intervals that overlap with the given range.
    ///
    /// # Arguments
    ///
    /// * `start` - The start of the range to check for overlaps.
    /// * `end` - The end of the range to check for overlaps.
    /// * `found` - A mutable vector to store the overlapping intervals' data.
    pub fn find_overlaps(&mut self, start: i32, end: i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(end);
        let mut i = self.idx;
        unsafe {
            while i > 0 {
                if start <= *self.ends.get_unchecked(i) {
                    found.push(self.data.get_unchecked(i).clone());
                    i -= 1;
                } else {
                    let branch_i = self.branch.get_unchecked(i);
                    if *branch_i >= i {
                        break;
                    }
                    i = *branch_i;
                }
            }
            if i == 0 {
                if start <= *self.ends.get_unchecked(0) && *self.starts.get_unchecked(0) <= end {
                    found.push(self.data.get_unchecked(0).clone());
                }
            }
        }
    }
    /// Counts the number of intervals that overlap with the given range.
    ///
    /// # Arguments
    ///
    /// * `start` - The start of the range to check for overlaps.
    /// * `end` - The end of the range to check for overlaps.
    ///
    /// # Returns
    ///
    /// The number of overlapping intervals.
    pub fn count_overlaps(&mut self, start: i32, end: i32) -> usize {
        if self.starts.is_empty() {
            return 0;
        }
        self.upper_bound(end);
        let mut found: usize = 0;
        let mut i = self.idx;
        unsafe {
            #[cfg(target_arch = "x86_64")]
            {
                use std::arch::x86_64::*;
                let start_vec = _mm256_set1_epi32(start);
                let ones: __m256i = _mm256_set1_epi32(1);
                const SIMD_WIDTH: usize = 8; //usize = 256 / (core::mem::size_of::<i32>() * 8);
                const BLOCK: usize = 32; //SIMD_WIDTH * 8;

                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i > BLOCK {
                            let mut count = 0;
                            for j in (i - BLOCK + 1..=i).rev().step_by(SIMD_WIDTH) {
                                let ends_vec = _mm256_loadu_si256(self.ends.as_ptr().add(j - SIMD_WIDTH + 1) as *const __m256i);
                                // Add one to convert the greater than to greater or equal than
                                let adj_ends_vec = _mm256_add_epi32(ends_vec, ones);
                                let cmp_mask = _mm256_cmpgt_epi32(adj_ends_vec, start_vec);
                                let mask = _mm256_movemask_epi8(cmp_mask);
                                count += mask.count_ones() as usize;
                            }
                            found += count / 4;  // Each comparison result is 4 bits
                            i -= BLOCK;
                            if count < BLOCK {
                                break;
                            }
                        }
                    } else {
                        if *self.branch.get_unchecked(i) >= i {
                            break;
                        }
                        i = *self.branch.get_unchecked(i);
                    }
                }
            }

            #[cfg(target_arch = "aarch64")]
            {
                use std::arch::aarch64::*;
                let start_vec = vdupq_n_s32(start);
                const SIMD_WIDTH: usize = 4; //128 / (core::mem::size_of::<i32>() * 8);
                const BLOCK: usize = 32; // SIMD_WIDTH * 4;
                let ones = vdupq_n_u32(1);
                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i > BLOCK {
                            let mut count = 0;
                            for j in (i - BLOCK + 1..=i).rev().step_by(SIMD_WIDTH) {
                                let ends_vec = vld1q_s32(self.ends.as_ptr().add(j - SIMD_WIDTH + 1) as *const i32);
                                let mask = vcleq_s32(start_vec, ends_vec);
                                let bool_mask = vandq_u32(mask, ones);
                                count += vaddvq_u32(bool_mask) as usize;
                            }
                            found += count;
                            i -= BLOCK;
                            if count < BLOCK {
                                break;
                            }
                        }
                    } else {
                        if *self.branch.get_unchecked(i) >= i {
                            break;
                        }
                        i = *self.branch.get_unchecked(i);
                    }
                }
            }

            #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
            {
                const BLOCK: usize = 16;
                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i > BLOCK {
                            let mut count = 0;
                            for j in (i - BLOCK + 1..=i).rev() {
                                if start <= *self.ends.get_unchecked(j) {
                                    count += 1;
                                }
                            }
                            found += count;
                            i -= BLOCK;
                            if count < BLOCK {
                                break;
                            }
                        }
                    } else {
                        if *self.branch.get_unchecked(i) >= i {
                            break;
                        }
                        i = *self.branch.get_unchecked(i);
                    }
                }
            }
            if i == 0 && start <= *self.ends.get_unchecked(0) && *self.starts.get_unchecked(0) <= end {
                found += 1;
            }
        }
        found
    }

    fn sort_block<F>(&mut self, start_i: usize, end_i: usize, compare: F)
    where
        F: Fn(&Interval<T>, &Interval<T>) -> Ordering,
    {
        unsafe {
            let range_size = end_i - start_i;
            let mut tmp: Vec<Interval<T>> = Vec::with_capacity(range_size);
            for i in 0..range_size {
                tmp.push(Interval {
                    start: *self.starts.get_unchecked(start_i + i),
                    end: *self.ends.get_unchecked(start_i + i),
                    data: (*self.data.get_unchecked(start_i + i)).clone(),
                });
            }
            tmp.sort_by(compare);
            for i in 0..range_size {
                self.starts[start_i + i] = tmp.get_unchecked(i).start;
                self.ends[start_i + i] = tmp.get_unchecked(i).end;
                self.data[start_i + i] = tmp.get_unchecked(i).data.clone();
            }
        }
    }
}


/// A variant of `SuperIntervals` that uses the Eytzinger layout for faster searching.
pub struct SuperIntervalsEytz<T> {
    inner: SuperIntervals<T>,
    eytz: Vec<i32>,
    eytz_index: Vec<usize>,
}

#[inline(always)]
pub fn ffs(x: u32) -> u32 {
    if x == 0 {
        0
    } else {
        x.trailing_zeros() + 1
    }
}

impl<T> SuperIntervalsEytz<T>
where
    T: Clone,
{
    pub fn new() -> Self {
        SuperIntervalsEytz {
            inner: SuperIntervals::new(),
            eytz: Vec::new(),
            eytz_index: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.inner.clear();
        self.eytz.clear();
        self.eytz_index.clear();
    }

    pub fn add(&mut self, start: i32, end: i32, value: T) {
        self.inner.add(start, end, value);
    }

    pub fn sort_intervals(&mut self) {
        self.inner.sort_intervals();
    }

    fn eytzinger_helper(&mut self, mut i: usize, k: usize, n: usize) -> usize {
        unsafe {
            if k < n {
                i = self.eytzinger_helper(i, 2*k+1, n);
                self.eytz[k] = *self.inner.starts.get_unchecked(i);
                self.eytz_index[k] = i;
                i += 1;
                i = self.eytzinger_helper(i, 2*k+2, n);
            }
        }
        i
    }

    pub fn index(&mut self) {
        if self.inner.starts.is_empty() {
            return;
        }
        self.sort_intervals();
        self.eytz.resize(self.inner.starts.len() + 1, 0);
        self.eytz_index.resize(self.inner.starts.len() + 1, 0);
        self.eytzinger_helper(0, 0, self.inner.starts.len());

        self.inner.branch.resize(self.inner.starts.len(), usize::MAX);
        unsafe {
            for i in 0..self.inner.ends.len() - 1 {
                for j in i + 1..self.inner.ends.len() {
                    if self.inner.ends.get_unchecked(j) >= self.inner.ends.get_unchecked(i) {
                        break;
                    }
                    self.inner.branch[j] = i;
                }
            }
        }
        self.inner.idx = 0;
    }

    pub fn upper_bound(&mut self, x: i32) {
        unsafe {
            let mut i: usize = 0;
            let n_intervals: usize = self.inner.starts.len();
            while i < n_intervals {
                if *self.eytz.get_unchecked(i) > x {
                    i = 2 * i + 1;
                } else {
                    i = 2 * i + 2;
                }
            }
            let shift: u32 = ffs(!(i as u32 + 1));
            let best_idx: usize = (i >> shift) - ( if shift > 1 { 1 } else { 0 } );
            self.inner.idx = if best_idx < n_intervals { *self.eytz_index.get_unchecked(best_idx) } else { n_intervals - 1 };
            if self.inner.idx > 0 && *self.inner.starts.get_unchecked(self.inner.idx) > x {
                self.inner.idx -= 1;
            }
        }
    }

    pub fn find_overlaps(&mut self, start: i32, end: i32, found: &mut Vec<T>) {
        self.inner.find_overlaps(start, end, found)
    }

    pub fn count_overlaps(&mut self, start: i32, end: i32) -> usize {
        self.inner.count_overlaps(start, end)
    }
}