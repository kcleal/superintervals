
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
    }

    #[inline(always)]
    pub fn upper_bound(&self, value: i32) -> usize {
        let mut idx = 0;
        let mut length = self.starts.len();
        unsafe {
            while length > 1 {
                let half = length / 2;
                idx += (*self.starts.get_unchecked(idx + half) <= value) as usize * (length - half);
                length = half;
            }
//             idx = idx.wrapping_sub((*self.starts.get_unchecked(idx) > value) as usize);
            if *self.starts.get_unchecked(idx) > value {
                idx = idx.wrapping_sub(1);
            }
        }
        return idx;
    }

    #[inline(always)]
    pub fn upper_bound_range(&self, value: i32, right: usize) -> usize {
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

    pub fn has_overlaps(&self, start: i32, end: i32) -> bool {
        if self.starts.is_empty() {
            return false;
        }
        let idx = self.upper_bound(end);
        unsafe {
            idx != usize::MAX && start <= *self.ends.get_unchecked(idx)
        }
    }

    fn setup_search(&self, start: i32, end: i32) -> Option<usize> {  // Return Option<usize>
        if self.starts.is_empty() {
            return None;
        }
        let idx = self.upper_bound(end);
        if idx == usize::MAX {
            return None;
        }
        unsafe {
            if start > *self.ends.get_unchecked(idx) || *self.starts.get_unchecked(0) > end {
                return None;
            }
        }
        Some(idx)
    }

    // Iterator over indices
    // for idx in tree.search_idxs_iter(10, 20) { ... }
    pub fn search_idxs_iter(&self, start: i32, end: i32) -> IndexIterator<T> {
        let current_idx = self.setup_search(start, end).unwrap_or(usize::MAX);
        IndexIterator {
            tree: self,
            current_idx,
            query_start: start,
        }
    }

    // Iterator over intervals
    pub fn search_items_iter(&self, start: i32, end: i32) -> ItemIterator<T> {
        let current_idx = self.setup_search(start, end).unwrap_or(usize::MAX);
        ItemIterator {
            tree: self,
            current_idx,
            query_start: start,
        }
    }


//     This is the simple algorithm, but suffers performance wise due to the branching
//     pub fn find_overlaps(&mut self, start: i32, end: i32, found: &mut Vec<T>) {
//         if self.starts.is_empty() {
//             return;
//         }
//         let mut i = self.upper_bound(end);
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
    pub fn search_values(&self, start: i32, end: i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        let idx = self.upper_bound(end);
        if idx == usize::MAX {
            return;
        }
        let mut i = idx;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                found.reserve(idx + 1);
                found.extend((0..=idx).rev().map(|j| {
                        self.data.get_unchecked(j).clone()
                    }));
                return;
            }
            let count = idx - i;
            found.reserve(count);
            found.extend(((i + 1)..=idx).rev().map(|j| {
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
    pub fn search_values_large(&self, start: i32, end: i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        let idx = self.upper_bound(end);
        if idx == usize::MAX {
            return;
        }
        let mut i = self.upper_bound_range(start, idx);
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                found.reserve(idx + 1);
                found.extend((0..=idx).rev().map(|j| {
                        self.data.get_unchecked(j).clone()
                    }));
                return;
            }
            found.reserve(idx - i);
            found.extend(((i + 1)..=idx).rev().map(|j| {
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
    pub fn count_linear(&self, start: i32, end: i32) -> usize {
        if self.starts.is_empty() {
            return 0;
        }
        let idx = self.upper_bound(end);
        if idx == usize::MAX {
            return 0;
        }
        let mut i = idx;
        let mut count: usize = 0;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                return idx + 1;
            }
            count += idx - i;
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
    /// Counts the number of intervals that overlap with the given range.
    pub fn count(&self, start: i32, end: i32) -> usize {
        if self.starts.is_empty() {
            return 0;
        }
        let idx = self.upper_bound(end);
        if idx == usize::MAX {
            return 0;
        }
        let mut i = idx;
        let mut found = 0;

        unsafe {
            #[cfg(target_arch = "x86_64")]
            {
                use std::arch::x86_64::*;
                let start_vec = _mm256_set1_epi32(start);
                const SIMD_WIDTH: usize = 8; // 256 bits / 32 bits = 8 elements
                const BLOCK: usize = SIMD_WIDTH * 4; // 32 elements per block

                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;

                        // Process blocks of data with SIMD
                        while i > BLOCK {
                            let mut count = 0;
                            let mut j = i;
// int mask = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask));
//                             count += _mm_popcnt_u32(mask);
                            // Process SIMD_WIDTH elements at a time
                            while j > i - BLOCK {
                                let end_idx = if j >= SIMD_WIDTH { j - SIMD_WIDTH + 1 } else { 0 };
                                let ends_vec = _mm256_loadu_si256(self.ends.as_ptr().add(end_idx) as *const __m256i);
                                let cmp_mask = _mm256_cmpgt_epi32(start_vec, ends_vec);
//                                 let mask = _mm256_movemask_epi8(cmp_mask);
                                let mask = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask));
                                // Count the number of set bits, each comparison result is 4 bits
//                                 count += (!mask).count_ones() as usize / 4;
                                count += (!mask).count_ones() as usize; // / 4;
                                j = j.saturating_sub(SIMD_WIDTH);
                            }

                            found += count;
                            i -= BLOCK;

                            // Early exit if we didn't find full block worth of matches
                            if count < BLOCK && start > *self.ends.get_unchecked(i + 1) {
                                break;
                            }
                        }
                    } else {
                        if *self.branch.get_unchecked(i) == usize::MAX {
                            return found;
                        }
                        i = *self.branch.get_unchecked(i);
                    }
                }
            }

            #[cfg(target_arch = "aarch64")]
            {
                use std::arch::aarch64::*;
                let start_vec = vdupq_n_s32(start);
                const SIMD_WIDTH: usize = 4; // 128 bits / 32 bits = 4 elements
                const BLOCK: usize = SIMD_WIDTH * 4; // 16 elements per block
                let ones = vdupq_n_u32(1);
                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i > BLOCK {
                            let mut count = 0;
                            let mut j = i;
                            while j > i - BLOCK {
                                let end_idx = if j >= SIMD_WIDTH { j - SIMD_WIDTH + 1 } else { 0 };
                                let ends_vec = vld1q_s32(self.ends.as_ptr().add(end_idx) as *const i32);
//                                 let mask = vcleq_s32(start_vec, ends_vec);
                                let mask = vcgtq_s32(start_vec, ends_vec);
//                                 let bool_mask = vandq_u32(mask, ones);
                                let bool_mask = vaddq_u32(mask, ones);
                                count += vaddvq_u32(bool_mask) as usize;
                                j = j.saturating_sub(SIMD_WIDTH);
                            }
                            found += count;
                            i -= BLOCK;
                            if count < BLOCK {
                                break;
                            }
                        }
                    } else {
                        if *self.branch.get_unchecked(i) == usize::MAX {
                            return found;
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
                            if count < BLOCK && start > *self.ends.get_unchecked(i + 1) {
                                break;
                            }
                        }
                    } else {
                        if *self.branch.get_unchecked(i) == usize::MAX {
                            return found;
                        }
                        i = *self.branch.get_unchecked(i);
                    }
                }
            }
            // Final check for element at index 0
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
    pub fn count_large(&self, start: i32, end: i32) -> usize {
        if self.starts.is_empty() {
            return 0;
        }
        let idx = self.upper_bound(end);
        if idx == usize::MAX {
            return 0;
        }
        let mut i = self.upper_bound_range(start, idx);
        let mut count: usize = 0;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                return idx + 1;
            }
            count += idx - i;
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

    pub fn search_idxs(&self, start: i32, end: i32, found: &mut Vec<usize>) {
        if self.starts.is_empty() {
            return;
        }
        let idx = self.upper_bound(end);
        if idx == usize::MAX {
            return;
        }
        let mut i = idx;
        unsafe {
            while i != usize::MAX && start <= *self.ends.get_unchecked(i) {
                i = i.wrapping_sub(1);
            }
            if i == usize::MAX {
                found.extend((0..=idx).rev());
                return;
            }
            found.extend(((i + 1)..=idx).rev());
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

    pub fn search_keys(&self, start: i32, end: i32, found: &mut Vec<(i32, i32)>) {
        if self.starts.is_empty() {
            return;
        }
        let idx = self.upper_bound(end);
        if idx == usize::MAX {
            return;
        }
        let mut i = idx;
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

    pub fn search_items(&self, start: i32, end: i32, found: &mut Vec<Interval<T>>) {
        if self.starts.is_empty() {
            return;
        }
        let idx = self.upper_bound(end);
        if idx == usize::MAX {
            return;
        }
        let mut i = idx;
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

    pub fn search_stabbed(&self, point: i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        let idx = self.upper_bound(point);
        if idx == usize::MAX {
            return;
        }
        let mut i = idx;
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
    pub fn coverage(&self, start: i32, end: i32) -> (usize, i32) {
        if self.starts.is_empty() {
            return (0, 0);
        }
        let mut i = self.upper_bound(end);
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
