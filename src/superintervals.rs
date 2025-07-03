
//! This module provides an associative data structure for performing interval intersection queries.

use std::cmp::Ordering;
use std::cmp::{max, min};
use serde::{Serialize, Deserialize};

/// Represents an interval with associated data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Interval<T> {
    pub start: i32,
    pub end: i32,
    pub data: T,
}

/// A static data structure for finding interval intersections
///
/// IntervalMap is a template class that provides efficient interval intersection operations.
/// It supports adding intervals, indexing them for fast queries, and performing various
/// intersection operations.
///
/// Intervals are considered end-inclusive
/// The build() function must be called before any queries. If more intervals are added, call build() again.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntervalMap<T>
{
    pub starts: Vec<i32>,
    pub ends: Vec<i32>,
    pub branch: Vec<usize>,
    pub data: Vec<T>,
    pub idx: usize,
    pub start_sorted: bool,
    pub end_sorted: bool,
}


impl<T: Clone> IntervalMap<T>
{
    pub fn new() -> Self {
        IntervalMap {
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

    pub fn reserve(&mut self, n: usize) {
        self.starts.reserve(n);
        self.ends.reserve(n);
        self.data.reserve(n);
    }

    pub fn size(&self) -> usize {
        self.starts.len()
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

    pub fn at(&self, index: usize) -> Interval<T> {
        Interval {
                start: self.starts[index],
                end: self.ends[index],
                data: self.data[index].clone(),
        }
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

    /// Build the index. Must be called before queries are performed.
    pub fn build(&mut self) {
        if self.starts.is_empty() {
            return;
        }
        self.sort_intervals();
        self.branch.resize(self.starts.len(), usize::MAX);
        let mut br: Vec<(i32, usize)> = Vec::with_capacity((self.starts.len() / 10) + 1);
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

    #[inline(always)]
    pub fn upper_bound(&mut self, value: i32) {
        let mut length = self.starts.len();
        unsafe {
            self.idx = 0;
            while length > 1 {
                let half = length / 2;
                self.idx += (*self.starts.get_unchecked(self.idx + half) <= value) as usize * (length - half);
                length = half;
            }
//             self.idx = self.idx.wrapping_sub((*self.starts.get_unchecked(self.idx) > value) as usize);
            if *self.starts.get_unchecked(self.idx) > value {
                self.idx = self.idx.wrapping_sub(1);
            }
        }
    }

    #[inline(always)]
    pub fn upper_bound_range(&mut self, value: i32, right: usize) -> usize {
        let mut left = right;
        let mut search_right = right;
        let mut bound = 1;
        // Exponential search to find a smaller range
        unsafe {
            while left > 0 && value < *self.starts.get_unchecked(left) {
                search_right = left;
                left = if bound <= left { left - bound } else { 0 };
                bound *= 2;
            }
            // Binary search in the found range
            let mut length = search_right - left;
            while length > 1 {
                let half = length / 2;
                let condition = *self.starts.get_unchecked(left + half) < value;
                left += if condition { length - half } else { 0 };
                length = half;
            }
            if left == 0 && *self.starts.get_unchecked(left) >= value {
                left = usize::MAX;
            }
        }
        left
    }

    pub fn has_overlaps(&mut self, start: i32, end: i32) -> bool {
        if self.starts.is_empty() {
            return false;
        }
        self.upper_bound(end);
        unsafe {
            self.idx != usize::MAX && start <= *self.ends.get_unchecked(self.idx)
        }
    }

    fn setup_search(&mut self, start: i32, end: i32) -> bool {
        if self.starts.is_empty() {
            return false;
        }
        self.upper_bound(end);
        if self.idx == usize::MAX {
            return false;
        }
        unsafe {
            if start > *self.ends.get_unchecked(self.idx) || *self.starts.get_unchecked(0) > end {
                self.idx = usize::MAX;
                return false;
            }
        }
        true
    }

    // Iterator over indices
    // for idx in tree.search_idxs_iter(10, 20) { ... }
    pub fn search_idxs_iter(&mut self, start: i32, end: i32) -> IndexIterator<T> {
        if !self.setup_search(start, end) {
            return IndexIterator {
                tree: self,
                current_idx: usize::MAX, // Empty iterator
                query_start: start,
            };
        }
        IndexIterator {
            tree: self,
            current_idx: self.idx,
            query_start: start,
        }
    }

    // Iterator over intervals
    pub fn search_items_iter(&mut self, start: i32, end: i32) -> ItemIterator<T> {
        if !self.setup_search(start, end) {
            return ItemIterator {
                tree: self,
                current_idx: usize::MAX, // Empty iterator
                query_start: start,
            };
        }
        ItemIterator {
            tree: self,
            current_idx: self.idx,
            query_start: start,
        }
    }


//     This is the simple algorithm, but suffers performance wise due to the branching
//     pub fn find_overlaps(&mut self, start: i32, end: i32, found: &mut Vec<T>) {
//         if self.starts.is_empty() {
//             return;
//         }
//         self.upper_bound(end);
//         let mut i = self.idx;
//
//         unsafe {
//             while i != usize::MAX {
//                 if start <= *self.ends.get_unchecked(i) {
//                     found.push(self.data.get_unchecked(i).clone());
//                     i = i.wrapping_sub(1);
//                 } else {
//                     i = *self.branch.get_unchecked(i);
//                 }
//             }
//         }
//     }

    /// Finds all intervals that overlap with the given range.
    ///
    /// # Arguments
    ///
    /// * `start` - The start of the range to check for overlaps.
    /// * `end` - The end of the range to check for overlaps.
    /// * `found` - A mutable vector to store the overlapping intervals' data.
    pub fn search_values(&mut self, start: i32, end: i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(end);
        if self.idx == usize::MAX {
            return;
        }
        let mut i = self.idx;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                found.reserve(self.idx + 1);
                found.extend((0..=self.idx).rev().map(|j| {
                        self.data.get_unchecked(j).clone()
                    }));
                return;
            }
            let count = self.idx - i;
            found.reserve(count);
            found.extend(((i + 1)..=self.idx).rev().map(|j| {
                    self.data.get_unchecked(j).clone()
                }));
            i = *self.branch.get_unchecked(i);
            while i != usize::MAX {
                if start <= *self.ends.get_unchecked(i) {
                    found.push(self.data.get_unchecked(i).clone());
                    i = i.wrapping_sub(1);
                } else {
                    i = *self.branch.get_unchecked(i);
                }
            }
        }
    }

    /// Finds all intervals that overlap with the given range. Works best when
    /// query intervals are large compared to stored database intervals. Uses
    /// an exponential search to find overlaps
    ///
    /// # Arguments
    ///
    /// * `start` - The start of the range to check for overlaps.
    /// * `end` - The end of the range to check for overlaps.
    /// * `found` - A mutable vector to store the overlapping intervals' data.
    pub fn search_values_large(&mut self, start: i32, end: i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(end);
        if self.idx == usize::MAX {
            return;
        }
        let mut i = self.upper_bound_range(start, self.idx);
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                found.reserve(self.idx + 1);
                found.extend((0..=self.idx).rev().map(|j| {
                        self.data.get_unchecked(j).clone()
                    }));
                return;
            }
            found.reserve(self.idx - i);
            found.extend(((i + 1)..=self.idx).rev().map(|j| {
                    self.data.get_unchecked(j).clone()
                }));
            i = *self.branch.get_unchecked(i);
            while i != usize::MAX {
                if start <= *self.ends.get_unchecked(i) {
                    found.push(self.data.get_unchecked(i).clone());
                    i = i.wrapping_sub(1);
                } else {
                    i = *self.branch.get_unchecked(i);
                }
            }
        }
    }

    /// Counts all intervals that overlap with the given range.
    ///
    /// # Arguments
    ///
    /// * `start` - The start of the range to check for overlaps.
    /// * `end` - The end of the range to check for overlaps.
    ///
    /// # Returns
    ///
    /// The number of overlapping intervals.
    pub fn count_linear(&mut self, start: i32, end: i32) -> usize {
        if self.starts.is_empty() {
            return 0;
        }
        self.upper_bound(end);
        if self.idx == usize::MAX {
            return 0;
        }
        let mut i = self.idx;
        let mut count: usize = 0;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                return self.idx + 1;
            }
            count += self.idx - i;
            i = *self.branch.get_unchecked(i);
            while i != usize::MAX {
                if start <= *self.ends.get_unchecked(i) {
                    count += 1;
                    i = i.wrapping_sub(1);
                } else {
                    i = *self.branch.get_unchecked(i);
                }
            }
        }
        count
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
    pub fn count(&mut self, start: i32, end: i32) -> usize {
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
                const SIMD_WIDTH: usize = 8;
                const BLOCK: usize = 32;
                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i > BLOCK {
                            let mut count = 0;
                            let mut last_count = 0;
                            for j in (i - BLOCK + 1..=i).rev().step_by(SIMD_WIDTH) {
                                let ends_vec = _mm256_loadu_si256(self.ends.as_ptr().add(j - SIMD_WIDTH + 1) as *const __m256i);
                                // Add one to convert the greater than to greater or equal than
                                let adj_ends_vec = _mm256_add_epi32(ends_vec, ones);
                                let cmp_mask = _mm256_cmpgt_epi32(adj_ends_vec, start_vec);
                                let mask = _mm256_movemask_epi8(cmp_mask);
                                last_count = mask.count_ones() as usize;
                                count += last_count;
                            }
                            found += count / 4;  // Each comparison result is 4 bits
                            i -= BLOCK;
                            if last_count == 0 {
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
                const SIMD_WIDTH: usize = 4;
                const BLOCK: usize = 32;
                let ones = vdupq_n_u32(1);
                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i > BLOCK {
                            let mut count = 0;
                            let mut last_count = 0;
                            for j in (i - BLOCK + 1..=i).rev().step_by(SIMD_WIDTH) {
                                let ends_vec = vld1q_s32(self.ends.as_ptr().add(j - SIMD_WIDTH + 1) as *const i32);
                                let mask = vcleq_s32(start_vec, ends_vec);
                                let bool_mask = vandq_u32(mask, ones);
                                last_count = vaddvq_u32(bool_mask) as usize;
                                count += last_count;
                            }
                            found += count;
                            i -= BLOCK;
                            if last_count == 0 {
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
                while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                    i = i.wrapping_sub(1);
                }
                if i == usize::MAX {
                    return self.idx + 1;
                }
                found += self.idx - i;
                i = *self.branch.get_unchecked(i);
                while i != usize::MAX {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i = i.wrapping_sub(1);
                    } else {
                        i = *self.branch.get_unchecked(i);
                    }
                }
                return found;
            }

            // Last check for 0
            if i == 0 && start <= *self.ends.get_unchecked(0) && *self.starts.get_unchecked(0) <= end {
                found += 1;
            }
        }
        found
    }

    /// Counts all intervals that overlap with the given range. Works best when
    /// query intervals are large compared to stored database intervals. Uses
    /// an exponential search to find overlaps
    ///
    /// # Arguments
    ///
    /// * `start` - The start of the range to check for overlaps.
    /// * `end` - The end of the range to check for overlaps.
    ///
    /// # Returns
    ///
    /// The number of overlapping intervals.
    pub fn count_large(&mut self, start: i32, end: i32) -> usize {
        if self.starts.is_empty() {
            return 0;
        }
        self.upper_bound(end);
        if self.idx == usize::MAX {
            return 0;
        }
        let mut i = self.upper_bound_range(start, self.idx);
        let mut count: usize = 0;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                return self.idx + 1;
            }
            count += self.idx - i;
            i = *self.branch.get_unchecked(i);
            while i != usize::MAX {
                if start <= *self.ends.get_unchecked(i) {
                    count += 1;
                    i = i.wrapping_sub(1);
                } else {
                    i = *self.branch.get_unchecked(i);
                }
            }
        }
        count
    }

    pub fn search_idxs(&mut self, start: i32, end: i32, found: &mut Vec<usize>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(end);
        if self.idx == usize::MAX {
            return;
        }
        let mut i = self.idx;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                found.extend((0..=self.idx).rev());
                return;
            }
            found.extend(((i + 1)..=self.idx).rev());
            i = *self.branch.get_unchecked(i);
            while i != usize::MAX {
                if start <= *self.ends.get_unchecked(i) {
                    found.push(i);
                    i = i.wrapping_sub(1);
                } else {
                    i = *self.branch.get_unchecked(i);
                }
            }
        }
    }

    pub fn search_keys(&mut self, start: i32, end: i32, found: &mut Vec<(i32, i32)>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(end);
        if self.idx == usize::MAX {
            return;
        }
        let mut i = self.idx;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                found.push((*self.starts.get_unchecked(i), *self.ends.get_unchecked(i)));
                i = i.wrapping_sub(1);
            }
            if i != usize::MAX {
                i = *self.branch.get_unchecked(i);
                while i != usize::MAX {
                    if start <= *self.ends.get_unchecked(i) {
                        found.push((*self.starts.get_unchecked(i), *self.ends.get_unchecked(i)));
                        i = i.wrapping_sub(1);
                    } else {
                        i = *self.branch.get_unchecked(i);
                    }
                }
            }
        }
    }

    pub fn search_items(&mut self, start: i32, end: i32, found: &mut Vec<Interval<T>>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(end);
        if self.idx == usize::MAX {
            return;
        }
        let mut i = self.idx;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                found.push(Interval {
                    start: *self.starts.get_unchecked(i),
                    end: *self.ends.get_unchecked(i),
                    data: self.data.get_unchecked(i).clone(),
                });
                i = i.wrapping_sub(1);
            }
            if i != usize::MAX {
                i = *self.branch.get_unchecked(i);
                while i != usize::MAX {
                    if start <= *self.ends.get_unchecked(i) {
                        found.push(Interval {
                            start: *self.starts.get_unchecked(i),
                            end: *self.ends.get_unchecked(i),
                            data: self.data.get_unchecked(i).clone(),
                        });
                        i = i.wrapping_sub(1);
                    } else {
                        i = *self.branch.get_unchecked(i);
                    }
                }
            }
        }
    }

    pub fn search_stabbed(&mut self, point: i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(point);
        if self.idx == usize::MAX {
            return;
        }
        let mut i = self.idx;
        unsafe {
            while i != usize::MAX && point <= *self.ends.get_unchecked(i) {
                found.push(self.data.get_unchecked(i).clone());
                i = i.wrapping_sub(1);
            }
            if i != usize::MAX {
                i = *self.branch.get_unchecked(i);
                while i != usize::MAX {
                    if point <= *self.ends.get_unchecked(i) {
                        found.push(self.data.get_unchecked(i).clone());
                        i = i.wrapping_sub(1);
                    } else {
                        i = *self.branch.get_unchecked(i);
                    }
                }
            }
        }
    }

    /// Counts the total coverage over the query interval.
    ///
    /// # Arguments
    ///
    /// * `start` - The start of the range to check for overlaps.
    /// * `end` - The end of the range to check for overlaps.
    ///
    /// # Returns
    ///
    /// The number of overlapping intervals, plus the coverage.
    pub fn coverage(&mut self, start: i32, end: i32) -> (usize, i32) {
        if self.starts.is_empty() {
            return (0, 0);
        }
        self.upper_bound(end);
        let mut i = self.idx;
        let mut count: usize = 0;
        let mut coverage: i32 = 0;
        unsafe {
            while i != usize::MAX {
                if start <= *self.ends.get_unchecked(i) {
                    count += 1;
                    coverage += min(*self.ends.get_unchecked(i), end) - max(*self.starts.get_unchecked(i), start);
                    i -= 1;
                } else {
                    i = *self.branch.get_unchecked(i);
                }
            }
        }
        (count, coverage)
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

pub struct IndexIterator<'a, T> {
    tree: &'a IntervalMap<T>,
    current_idx: usize,
    query_start: i32,
}

impl<'a, T: Clone> Iterator for IndexIterator<'a, T> {
    type Item = usize;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_idx == usize::MAX {
            return None;
        }
        let result = self.current_idx;
        unsafe {
            if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                self.current_idx = self.current_idx.wrapping_sub(1);
            } else {
                self.current_idx = *self.tree.branch.get_unchecked(self.current_idx);
            }
        }
        Some(result)
    }
}

pub struct ItemIterator<'a, T> {
    tree: &'a IntervalMap<T>,
    current_idx: usize,
    query_start: i32,
}

impl<'a, T: Clone> Iterator for ItemIterator<'a, T> {
    type Item = Interval<T>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_idx == usize::MAX {
            return None;
        }
        let result = Interval {
            start: self.tree.starts[self.current_idx],
            end: self.tree.ends[self.current_idx],
            data: self.tree.data[self.current_idx].clone(),
        };
        unsafe {
            if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                self.current_idx = self.current_idx.wrapping_sub(1);
            } else {
                self.current_idx = *self.tree.branch.get_unchecked(self.current_idx);
            }
        }
        Some(result)
    }
}
