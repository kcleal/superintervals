
//! This module provides an associative data structure for performing interval intersection queries.

use std::cmp::Ordering;
use std::cmp::{max, min};
use serde::{Serialize, Deserialize};
use aligned_vec::{AVec, ConstAlign};

/// Represents an interval with associated data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Interval<T> {
    pub start: i32,
    pub end: i32,
    pub data: T,
}


#[cfg(target_feature = "avx2")]
type AlignedEnds = AVec<i32, ConstAlign<32>>;

#[cfg(all(target_feature = "neon", not(target_feature = "avx2")))]
type AlignedEnds = AVec<i32, ConstAlign<16>>;

#[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
type AlignedEnds = Vec<i32>;

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
//     pub ends: Vec<i32>,
    pub ends: AlignedEnds,
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
//             ends: Vec::new(),
            ends: {
                #[cfg(target_feature = "avx2")]
                { AVec::new(32) }

                #[cfg(all(target_feature = "neon", not(target_feature = "avx2")))]
                { AVec::new(16) }

                #[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
                { Vec::new() }
            },

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
            #[cfg(all(target_arch = "x86_64", target_feature = "avx2", not(feature = "nosimd")))]
            {
                use std::arch::x86_64::*;
                let start_vec = _mm256_set1_epi32(start);
                const SIMD_WIDTH: usize = 8; // 256 bits / 32 bits = 8 elements
                const BLOCK: usize = SIMD_WIDTH * 4; // 32 elements per block, 2 cache lines

                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;

                        // Process blocks of data with SIMD
//                         while i > BLOCK {
//                             let mut count = 0;
//                             let mut j = i;
//                             // Process SIMD_WIDTH elements at a time
//                             while j > i - BLOCK {
//                                 let end_idx = if j >= SIMD_WIDTH { j - SIMD_WIDTH + 1 } else { 0 };
//                                 let ends_vec = _mm256_loadu_si256(self.ends.as_ptr().add(end_idx) as *const __m256i);
//                                 let cmp_mask = _mm256_cmpgt_epi32(start_vec, ends_vec);
//                                 let mask = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask));
//                                 count += 8 - (mask).count_ones() as usize;
//                                 j = j.saturating_sub(SIMD_WIDTH);
//                             }
//
//                             found += count;
//                             i -= BLOCK;
//
//                             // Early exit if we didn't find full block worth of matches
//                             if count < BLOCK && start > *self.ends.get_unchecked(i + 1) {
//                                 break;
//                             }
//                         }
                        while i > BLOCK {
                            let j = i - BLOCK + 1;

                            // Load all 4 vectors
                            let ends_vec0 = _mm256_loadu_si256(self.ends.as_ptr().add(j) as *const __m256i);
                            let ends_vec1 = _mm256_loadu_si256(self.ends.as_ptr().add(j + SIMD_WIDTH) as *const __m256i);
                            let ends_vec2 = _mm256_loadu_si256(self.ends.as_ptr().add(j + 2 * SIMD_WIDTH) as *const __m256i);
                            let ends_vec3 = _mm256_loadu_si256(self.ends.as_ptr().add(j + 3 * SIMD_WIDTH) as *const __m256i);

                            // Compare all vectors
                            let cmp_mask0 = _mm256_cmpgt_epi32(start_vec, ends_vec0);
                            let cmp_mask1 = _mm256_cmpgt_epi32(start_vec, ends_vec1);
                            let cmp_mask2 = _mm256_cmpgt_epi32(start_vec, ends_vec2);
                            let cmp_mask3 = _mm256_cmpgt_epi32(start_vec, ends_vec3);

                            // Extract masks
                            let mask0 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask0));
                            let mask1 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask1));
                            let mask2 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask2));
                            let mask3 = _mm256_movemask_ps(_mm256_castsi256_ps(cmp_mask3));

                            // Count and accumulate
                            let count = (8 - (mask0).count_ones() as usize) + (8 - (mask1).count_ones() as usize) +
                                        (8 - (mask2).count_ones() as usize) + (8 - (mask3).count_ones() as usize);

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

            #[cfg(all(target_arch = "aarch64", target_feature = "neon", not(feature = "nosimd")))]
            {
                use std::arch::aarch64::*;
                let start_vec = vdupq_n_s32(start);
                const SIMD_WIDTH: usize = 4; // 128 bits / 32 bits = 4 elements
                const BLOCK: usize = SIMD_WIDTH * 8; // 2 cache lines
                let ones = vdupq_n_u32(1);
                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
//                         while i > BLOCK {
//                             let mut count = 0;
//                             let mut j = i;
//                             while j > i - BLOCK {
//                                 let end_idx = if j >= SIMD_WIDTH { j - SIMD_WIDTH + 1 } else { 0 };
//                                 let ends_vec = vld1q_s32(self.ends.as_ptr().add(end_idx) as *const i32);
//                                 let mask = vcgtq_s32(start_vec, ends_vec);
//                                 let bool_mask = vaddq_u32(mask, ones);
//                                 count += vaddvq_u32(bool_mask) as usize;
//                                 j = j.saturating_sub(SIMD_WIDTH);
//                             }
//                             found += count;
//                             i -= BLOCK;
//                             if count < BLOCK {
//                                 break;
//                             }
//                         }

                        while i > BLOCK {
                            let j = i - BLOCK + 1;

                            // Load all 8 vectors
                            let ends_vec0 = vld1q_s32(self.ends.as_ptr().add(j) as *const i32);
                            let ends_vec1 = vld1q_s32(self.ends.as_ptr().add(j + SIMD_WIDTH) as *const i32);
                            let ends_vec2 = vld1q_s32(self.ends.as_ptr().add(j + 2 * SIMD_WIDTH) as *const i32);
                            let ends_vec3 = vld1q_s32(self.ends.as_ptr().add(j + 3 * SIMD_WIDTH) as *const i32);
                            let ends_vec4 = vld1q_s32(self.ends.as_ptr().add(j + 4 * SIMD_WIDTH) as *const i32);
                            let ends_vec5 = vld1q_s32(self.ends.as_ptr().add(j + 5 * SIMD_WIDTH) as *const i32);
                            let ends_vec6 = vld1q_s32(self.ends.as_ptr().add(j + 6 * SIMD_WIDTH) as *const i32);
                            let ends_vec7 = vld1q_s32(self.ends.as_ptr().add(j + 7 * SIMD_WIDTH) as *const i32);

                            // Compare all vectors
                            let mask0 = vcgtq_s32(start_vec, ends_vec0);
                            let mask1 = vcgtq_s32(start_vec, ends_vec1);
                            let mask2 = vcgtq_s32(start_vec, ends_vec2);
                            let mask3 = vcgtq_s32(start_vec, ends_vec3);
                            let mask4 = vcgtq_s32(start_vec, ends_vec4);
                            let mask5 = vcgtq_s32(start_vec, ends_vec5);
                            let mask6 = vcgtq_s32(start_vec, ends_vec6);
                            let mask7 = vcgtq_s32(start_vec, ends_vec7);

                            // Convert to boolean masks
                            let bool_mask0 = vaddq_u32(mask0, ones);
                            let bool_mask1 = vaddq_u32(mask1, ones);
                            let bool_mask2 = vaddq_u32(mask2, ones);
                            let bool_mask3 = vaddq_u32(mask3, ones);
                            let bool_mask4 = vaddq_u32(mask4, ones);
                            let bool_mask5 = vaddq_u32(mask5, ones);
                            let bool_mask6 = vaddq_u32(mask6, ones);
                            let bool_mask7 = vaddq_u32(mask7, ones);

                            // Sum all lanes and accumulate
                            let count = vaddvq_u32(bool_mask0) as usize + vaddvq_u32(bool_mask1) as usize +
                                        vaddvq_u32(bool_mask2) as usize + vaddvq_u32(bool_mask3) as usize +
                                        vaddvq_u32(bool_mask4) as usize + vaddvq_u32(bool_mask5) as usize +
                                        vaddvq_u32(bool_mask6) as usize + vaddvq_u32(bool_mask7) as usize;

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

            #[cfg(not(any(
                all(target_arch = "x86_64", target_feature = "avx2", not(feature = "nosimd")),
                all(target_arch = "aarch64", target_feature = "neon", not(feature = "nosimd"))
            )))]
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

    // -------------------------------------------------------------------------
    // Set-like operations.
    //
    // Each operation returns a NEW IntervalMap that is NOT indexed: build() has
    // not been called on it. Call build() on the result before querying it (or
    // before passing it as `other` to another set operation). This keeps
    // intermediate results in a chain from paying for an index they never use.
    // Intervals are end-inclusive throughout.
    //
    // Two intervals overlap when a.start <= b.end && b.start <= a.end. Merging
    // coalesces overlapping intervals only; exactly adjacent intervals such as
    // [1,5] and [6,10] are NOT merged. Callers who want adjacency merging should
    // widen ends by one before adding.
    //
    // Where intervals fold into one output interval, the data payload is decided
    // by a combiner closure. Variants without the `_with` suffix keep the first
    // interval's data.
    // -------------------------------------------------------------------------

    /// Returns indices into starts/ends/data in ascending (start, end) order.
    /// Cheap when already sorted (the common post-build() case).
    fn sorted_order(&self) -> Vec<usize> {
        let mut order: Vec<usize> = (0..self.starts.len()).collect();
        // start_sorted/end_sorted are maintained incrementally by add()/build(). When
        // both hold the data is fully (start,end)-ordered, so skip the verification scan
        // and the sort. (start_sorted alone is insufficient: equal-start intervals may
        // have unordered ends, which would break callers that need exact grouping.)
        let mut sorted = self.start_sorted && self.end_sorted;
        if !sorted {
            sorted = true;
            for k in 1..self.starts.len() {
                if self.starts[k] < self.starts[k - 1]
                    || (self.starts[k] == self.starts[k - 1] && self.ends[k] < self.ends[k - 1])
                {
                    sorted = false;
                    break;
                }
            }
        }
        if !sorted {
            order.sort_by(|&x, &y| {
                self.starts[x]
                    .cmp(&self.starts[y])
                    .then(self.ends[x].cmp(&self.ends[y]))
            });
        }
        order
    }

    /// Coalesces all stored intervals into a disjoint, non-overlapping set,
    /// combining the data of merged intervals with `combine`.
    /// Returns a new, UNINDEXED IntervalMap (call build() before querying).
    pub fn merge_overlaps_with<F>(&self, mut combine: F) -> IntervalMap<T>
    where
        F: FnMut(&T, &T) -> T,
    {
        let mut out = IntervalMap::new();
        if self.starts.is_empty() {
            return out;
        }
        let order = self.sorted_order();
        let mut cur_start = self.starts[order[0]];
        let mut cur_end = self.ends[order[0]];
        let mut cur_data = self.data[order[0]].clone();
        for k in 1..order.len() {
            let idx = order[k];
            if self.starts[idx] <= cur_end {
                if self.ends[idx] > cur_end {
                    cur_end = self.ends[idx];
                }
                cur_data = combine(&cur_data, &self.data[idx]);
            } else {
                out.add(cur_start, cur_end, cur_data);
                cur_start = self.starts[idx];
                cur_end = self.ends[idx];
                cur_data = self.data[idx].clone();
            }
        }
        out.add(cur_start, cur_end, cur_data);
        out
    }

    /// Coalesces overlapping intervals, keeping the first interval's data.
    pub fn merge_overlaps(&self) -> IntervalMap<T> {
        self.merge_overlaps_with(|a, _| a.clone())
    }

    /// Computes the gaps (uncovered regions) within `[lo, hi]`.
    /// Gap intervals match no input, so they carry `fill` as their data.
    /// Returns a new, UNINDEXED IntervalMap (call build() before querying).
    pub fn gaps(&self, lo: i32, hi: i32, fill: T) -> IntervalMap<T> {
        let mut out = IntervalMap::new();
        let merged = self.merge_overlaps();
        let mut cursor = lo;
        for k in 0..merged.starts.len() {
            let s = merged.starts[k];
            let e = merged.ends[k];
            if e < lo || s > hi {
                continue;
            }
            if s > cursor {
                out.add(cursor, s - 1, fill.clone());
            }
            if e + 1 > cursor {
                cursor = e + 1;
            }
        }
        if cursor <= hi {
            out.add(cursor, hi, fill.clone());
        }
        out
    }

    /// Union of this set with another, coalescing overlaps and combining the
    /// data of merged intervals with `combine`.
    /// Returns a new, UNINDEXED IntervalMap (call build() before querying).
    pub fn union_with<F>(&self, other: &IntervalMap<T>, combine: F) -> IntervalMap<T>
    where
        F: FnMut(&T, &T) -> T,
    {
        let mut combined = IntervalMap::new();
        combined.reserve(self.starts.len() + other.starts.len());
        for k in 0..self.starts.len() {
            combined.add(self.starts[k], self.ends[k], self.data[k].clone());
        }
        for k in 0..other.starts.len() {
            combined.add(other.starts[k], other.ends[k], other.data[k].clone());
        }
        combined.merge_overlaps_with(combine)
    }

    /// Union of two sets, keeping the first interval's data on each merge.
    pub fn union(&self, other: &IntervalMap<T>) -> IntervalMap<T> {
        self.union_with(other, |a, _| a.clone())
    }

    /// Intersection of this set with another: every overlapping sub-region.
    /// `other` must already be built (it is queried via its index). The output
    /// data of each piece is produced by `combine(a_data, b_data)`.
    /// The result is NOT coalesced; call merge_overlaps() if a disjoint set is needed.
    /// Returns a new, UNINDEXED IntervalMap (call build() before querying).
    pub fn intersection_with<F>(&self, other: &IntervalMap<T>, mut combine: F) -> IntervalMap<T>
    where
        F: FnMut(&T, &T) -> T,
    {
        let mut out = IntervalMap::new();
        let mut hits: Vec<usize> = Vec::new();
        for k in 0..self.starts.len() {
            hits.clear();
            other.search_idxs(self.starts[k], self.ends[k], &mut hits);
            for &h in &hits {
                let s = max(self.starts[k], other.starts[h]);
                let e = min(self.ends[k], other.ends[h]);
                if s <= e {
                    out.add(s, e, combine(&self.data[k], &other.data[h]));
                }
            }
        }
        out
    }

    /// Intersection of two sets, keeping this side's data for each piece.
    pub fn intersection(&self, other: &IntervalMap<T>) -> IntervalMap<T> {
        self.intersection_with(other, |a, _| a.clone())
    }

    /// Difference: regions of this set not covered by `other` (A \ B).
    /// `other` must already be built. Output pieces are sub-ranges of this set's
    /// intervals and inherit their source data.
    /// Returns a new, UNINDEXED IntervalMap (call build() before querying).
    pub fn difference(&self, other: &IntervalMap<T>) -> IntervalMap<T> {
        let mut out = IntervalMap::new();
        let mut cover: Vec<(i32, i32)> = Vec::new();
        for k in 0..self.starts.len() {
            cover.clear();
            other.search_keys(self.starts[k], self.ends[k], &mut cover);
            cover.sort();
            let mut cursor = self.starts[k];
            for &(cs0, ce0) in &cover {
                let cs = max(cs0, self.starts[k]);
                let ce = min(ce0, self.ends[k]);
                if cs > cursor {
                    out.add(cursor, cs - 1, self.data[k].clone());
                }
                if ce + 1 > cursor {
                    cursor = ce + 1;
                }
            }
            if cursor <= self.ends[k] {
                out.add(cursor, self.ends[k], self.data[k].clone());
            }
        }
        out
    }

    /// Symmetric difference: regions in exactly one of the two sets,
    /// (A \ B) ∪ (B \ A). Both inputs must already be built.
    /// Returns a new, UNINDEXED IntervalMap (call build() before querying).
    pub fn symmetric_difference(&self, other: &IntervalMap<T>) -> IntervalMap<T> {
        let a_minus_b = self.difference(other);
        let b_minus_a = other.difference(self);
        a_minus_b.union(&b_minus_a)
    }

    /// Smallest interval enclosing every stored interval (convex hull / span).
    /// Returns None for an empty set.
    pub fn span(&self) -> Option<(i32, i32)> {
        if self.starts.is_empty() {
            return None;
        }
        let mut lo = self.starts[0];
        let mut hi = self.ends[0];
        for k in 1..self.starts.len() {
            if self.starts[k] < lo {
                lo = self.starts[k];
            }
            if self.ends[k] > hi {
                hi = self.ends[k];
            }
        }
        Some((lo, hi))
    }

    /// Grows (or shrinks) every interval, like `bedtools slop`.
    /// `left` is subtracted from each start, `right` is added to each end;
    /// negative values shrink (intervals that shrink past themselves are dropped).
    /// `lo`/`hi` clamp the resulting coordinates. Data is carried over unchanged.
    /// Returns a new, UNINDEXED IntervalMap (call build() before querying).
    pub fn expand(&self, left: i32, right: i32, lo: i32, hi: i32) -> IntervalMap<T> {
        let mut out = IntervalMap::new();
        out.reserve(self.starts.len());
        for k in 0..self.starts.len() {
            // saturating arithmetic avoids overflow when lo/hi are the type extremes.
            let mut s = self.starts[k].saturating_sub(left);
            if s < lo {
                s = lo;
            }
            let mut e = self.ends[k].saturating_add(right);
            if e > hi {
                e = hi;
            }
            if s <= e {
                out.add(s, e, self.data[k].clone());
            }
        }
        out
    }

    /// Creates flanking intervals beside each interval, like `bedtools flank`.
    /// Emits the left flank `[start-left, start-1]` and right flank `[end+1, end+right]`
    /// (NOT the originals). Zero-width or out-of-range flanks are skipped; each flank
    /// inherits its source interval's data, clamped to `[lo, hi]`.
    /// Returns a new, UNINDEXED IntervalMap (call build() before querying).
    pub fn flank(&self, left: i32, right: i32, lo: i32, hi: i32) -> IntervalMap<T> {
        let mut out = IntervalMap::new();
        for k in 0..self.starts.len() {
            if left > 0 && self.starts[k] > lo {
                let le = self.starts[k] - 1;
                let mut ls = self.starts[k].saturating_sub(left);
                if ls < lo {
                    ls = lo;
                }
                if ls <= le && le <= hi {
                    out.add(ls, le, self.data[k].clone());
                }
            }
            if right > 0 && self.ends[k] < hi {
                let rs = self.ends[k] + 1;
                let mut re = self.ends[k].saturating_add(right);
                if re > hi {
                    re = hi;
                }
                if rs <= re && rs >= lo {
                    out.add(rs, re, self.data[k].clone());
                }
            }
        }
        out
    }

    /// Keeps only one interval per distinct (start, end) coordinate pair, folding
    /// the data of exact duplicates with `combine`. Unlike merge_overlaps(), this
    /// collapses only EXACT coordinate duplicates; distinct-but-overlapping
    /// intervals are all preserved. Output is start-sorted.
    /// Returns a new, UNINDEXED IntervalMap (call build() before querying).
    pub fn unique_with<F>(&self, mut combine: F) -> IntervalMap<T>
    where
        F: FnMut(&T, &T) -> T,
    {
        let mut out = IntervalMap::new();
        if self.starts.is_empty() {
            return out;
        }
        let order = self.sorted_order();
        let mut run = order[0];
        let mut acc = self.data[order[0]].clone();
        for k in 1..order.len() {
            let idx = order[k];
            if self.starts[idx] == self.starts[run] && self.ends[idx] == self.ends[run] {
                acc = combine(&acc, &self.data[idx]);
            } else {
                out.add(self.starts[run], self.ends[run], acc);
                run = idx;
                acc = self.data[idx].clone();
            }
        }
        out.add(self.starts[run], self.ends[run], acc);
        out
    }

    /// Keeps one interval per distinct (start, end) pair, keeping the first data value.
    pub fn unique(&self) -> IntervalMap<T> {
        self.unique_with(|a, _| a.clone())
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

// Iterator interfaces

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
        unsafe {
            // Linear scan backwards
            if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                let value = self.current_idx;
                self.current_idx = self.current_idx.wrapping_sub(1);
                return Some(value);
            }
            // Traverse branch array
            loop {
                self.current_idx = *self.tree.branch.get_unchecked(self.current_idx);
                if self.current_idx == usize::MAX {
                    break;
                }
                if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                    let value = self.current_idx;
                    self.current_idx = self.current_idx.wrapping_sub(1);
                    return Some(value);
                }
            }
        }
        None
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
        unsafe {
            if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                let value = Interval {
                    start: *self.tree.starts.get_unchecked(self.current_idx),
                    end: *self.tree.ends.get_unchecked(self.current_idx),
                    data: self.tree.data.get_unchecked(self.current_idx).clone(),
                };
                self.current_idx = self.current_idx.wrapping_sub(1);
                return Some(value);
            }
            loop {
                self.current_idx = *self.tree.branch.get_unchecked(self.current_idx);
                if self.current_idx == usize::MAX {
                    break;
                }
                if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                    let value = Interval {
                        start: *self.tree.starts.get_unchecked(self.current_idx),
                        end: *self.tree.ends.get_unchecked(self.current_idx),
                        data: self.tree.data.get_unchecked(self.current_idx).clone(),
                    };
                    self.current_idx = self.current_idx.wrapping_sub(1);
                    return Some(value);
                }
            }
        }
        None
    }
}

pub struct KeyIterator<'a, T> {
    tree: &'a IntervalMap<T>,
    current_idx: usize,
    query_start: i32,
}

impl<'a, T: Clone> Iterator for KeyIterator<'a, T> {
    type Item = (i32, i32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_idx == usize::MAX {
            return None;
        }
        unsafe {
            if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                let value = (
                    *self.tree.starts.get_unchecked(self.current_idx),
                    *self.tree.ends.get_unchecked(self.current_idx),
                );
                self.current_idx = self.current_idx.wrapping_sub(1);
                return Some(value);
            }
            loop {
                self.current_idx = *self.tree.branch.get_unchecked(self.current_idx);
                if self.current_idx == usize::MAX {
                    break;
                }
                if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                    let value = (
                        *self.tree.starts.get_unchecked(self.current_idx),
                        *self.tree.ends.get_unchecked(self.current_idx),
                    );
                    self.current_idx = self.current_idx.wrapping_sub(1);
                    return Some(value);
                }
            }
        }
        None
    }
}

pub struct ValueIterator<'a, T> {
    tree: &'a IntervalMap<T>,
    current_idx: usize,
    query_start: i32,
}

impl<'a, T: Clone> Iterator for ValueIterator<'a, T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_idx == usize::MAX {
            return None;
        }
        unsafe {
            if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                let value = self.tree.data.get_unchecked(self.current_idx).clone();
                self.current_idx = self.current_idx.wrapping_sub(1);
                return Some(value);
            }
            loop {
                self.current_idx = *self.tree.branch.get_unchecked(self.current_idx);
                if self.current_idx == usize::MAX {
                    break;
                }
                if self.query_start <= *self.tree.ends.get_unchecked(self.current_idx) {
                    let value = self.tree.data.get_unchecked(self.current_idx).clone();
                    self.current_idx = self.current_idx.wrapping_sub(1);
                    return Some(value);
                }
            }
        }
        None
    }
}

// Updated iterator creation methods for IntervalMap
impl<T: Clone> IntervalMap<T> {
    /// Returns an iterator over indices of intervals that intersect [start, end]
    pub fn search_idxs_iter(&self, start: i32, end: i32) -> IndexIterator<T> {
        let current_idx = if self.starts.is_empty() {
            usize::MAX
        } else {
            self.upper_bound(end)
        };
        IndexIterator {
            tree: self,
            current_idx,
            query_start: start,
        }
    }

    /// Returns an iterator over items (intervals with data) that intersect [start, end]
    pub fn search_items_iter(&self, start: i32, end: i32) -> ItemIterator<T> {
        let current_idx = if self.starts.is_empty() {
            usize::MAX
        } else {
            self.upper_bound(end)
        };
        ItemIterator {
            tree: self,
            current_idx,
            query_start: start,
        }
    }

    /// Returns an iterator over keys (interval start/end pairs) that intersect [start, end]
    pub fn search_keys_iter(&self, start: i32, end: i32) -> KeyIterator<T> {
        let current_idx = if self.starts.is_empty() {
            usize::MAX
        } else {
            self.upper_bound(end)
        };
        KeyIterator {
            tree: self,
            current_idx,
            query_start: start,
        }
    }

    /// Returns an iterator over values (data) of intervals that intersect [start, end]
    pub fn search_values_iter(&self, start: i32, end: i32) -> ValueIterator<T> {
        let current_idx = if self.starts.is_empty() {
            usize::MAX
        } else {
            self.upper_bound(end)
        };
        ValueIterator {
            tree: self,
            current_idx,
            query_start: start,
        }
    }
}
