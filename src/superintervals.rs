use std::cmp::Ordering;

#[derive(Debug, Clone)]
pub struct Interval<T> {
    start: i32,
    end: i32,
    data: T,
}

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

    pub fn clear(&mut self) {
        self.starts.clear();
        self.ends.clear();
        self.branch.clear();
        self.data.clear();
        self.idx = 0;
    }

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

    pub fn index(&mut self) {
        if self.starts.is_empty() {
            return;
        }
        self.sort_intervals();
        self.branch.resize(self.starts.len(), usize::MAX);
        for i in 0..self.ends.len() - 1 {
            for j in i + 1..self.ends.len() {
                if self.ends[j] >= self.ends[i] {
                    break;
                }
                self.branch[j] = i;
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
        const ENTRIES_PER_256KB: usize = 256 * 1024 / std::mem::size_of::<i32>();

        if length >= ENTRIES_PER_256KB {
            let num_per_cache_line: usize = std::cmp::max(64 / std::mem::size_of::<i32>(), 1);
            while length >= 3 * num_per_cache_line {
                let half = length / 2;

                // Using unsafe block to remove bounds checking
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
        }

        while length > 0 {
            let half = length / 2;

            // Using unsafe block to remove bounds checking
            unsafe {
                if *self.starts.get_unchecked(self.idx + half) <= value {
                    self.idx += length - half;
                }
            }

            length = half;
        }

        // Adjust index if necessary
        if self.idx > 0 && (self.idx == self.starts.len() || self.starts[self.idx] > value) {
            self.idx -= 1;
        }
    }

    pub fn find_overlaps(&mut self, start: i32, end: i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(end);
        let mut i = self.idx;
        unsafe {  // remove bounds checking
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
                const SIMD_WIDTH: usize = 8; // 256 bits / 32 bits
                const BLOCK: usize = SIMD_WIDTH * 4;

                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i >= BLOCK {
                            let mut count: usize = 0;
                            for j in (i - BLOCK + 1..=i).rev().step_by(SIMD_WIDTH) {
                                let ends_vec = _mm256_loadu_si256(self.ends.get_unchecked(j - SIMD_WIDTH + 1) as *const i32 as *const __m256i);
                                let cmp_mask = _mm256_cmpgt_epi32(start_vec, ends_vec);
                                let mask = _mm256_movemask_epi8(_mm256_xor_si256(cmp_mask, _mm256_set1_epi32(-1)));
                                count += (_popcnt32(mask) / 4) as usize;
                            }
                            found += count;
                            i -= BLOCK;
                            if count < BLOCK {
                                break;
                            }
                        }
                    } else {
                        let branch_i = *self.branch.get_unchecked(i);
                        if branch_i >= i {
                            break;
                        }
                        i = branch_i;
                    }
                }
            }

            #[cfg(target_arch = "aarch64")]
            {
                use std::arch::aarch64::*;
                let start_vec = vdupq_n_s32(start);
                const SIMD_WIDTH: usize = 4; // 128 bits / 32 bits
                const BLOCK: usize = SIMD_WIDTH * 4;

                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i >= BLOCK {
                            let mut count = 0;
                            for j in (i - BLOCK + 1..=i).step_by(SIMD_WIDTH) {
                                let ends_vec = vld1q_s32(self.ends.get_unchecked(j - SIMD_WIDTH + 1) as *const i32);
                                let mask = vcleq_s32(start_vec, ends_vec);
                                count += vaddvq_u32(vreinterpretq_u32_s32(mask)) as usize;
                            }
                            found += count;
                            i -= BLOCK;
                            if count < BLOCK {
                                break;
                            }
                        }
                    } else {
                        let branch_i = *self.branch.get_unchecked(i);
                        if branch_i >= i {
                            break;
                        }
                        i = branch_i;
                    }
                }
            }

            // Process remaining elements
            while i > 0 {
                if start <= *self.ends.get_unchecked(i) {
                    found += 1;
                    i -= 1;
                } else {
                    let branch_i = *self.branch.get_unchecked(i);
                    if branch_i >= i {
                        break;
                    }
                    i = branch_i;
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
        let range_size = end_i - start_i;
        let mut tmp: Vec<Interval<T>> = Vec::with_capacity(range_size);
        for i in 0..range_size {
            tmp.push(Interval {
                start: self.starts[start_i + i],
                end: self.ends[start_i + i],
                data: self.data[start_i + i].clone(),
            });
        }
        tmp.sort_by(compare);
        for i in 0..range_size {
            self.starts[start_i + i] = tmp[i].start;
            self.ends[start_i + i] = tmp[i].end;
            self.data[start_i + i] = tmp[i].data.clone();
        }
    }
}


// Extended superintervals

pub struct SuperIntervalsE<T> {
    starts: Vec<i32>,
    ends: Vec<i32>,
    branch: Vec<usize>,
    data: Vec<T>,
    eytz: Vec<i32>,
    eytz_index: Vec<usize>,
    idx: usize,
    start_sorted: bool,
    end_sorted: bool,
}

impl<T> SuperIntervalsE<T>
where
    T: Clone,
{
    pub fn new() -> Self {
        SuperIntervalsE {
            starts: Vec::new(),
            ends: Vec::new(),
            branch: Vec::new(),
            data: Vec::new(),
            eytz: Vec::new(),
            eytz_index: Vec::new(),
            idx: 0,
            start_sorted: true,
            end_sorted: true,
        }
    }

    pub fn clear(&mut self) {
        self.starts.clear();
        self.ends.clear();
        self.branch.clear();
        self.data.clear();
        self.idx = 0;
    }

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

    fn eytzinger_helper(&mut self, mut i: usize, k: usize, n: usize) -> usize {
        if k < n {
            i = self.eytzinger_helper(i, 2*k+1, n);
            self.eytz[k] = self.starts[i];
            self.eytz_index[k] = i;
            i += 1;
            i = self.eytzinger_helper(i, 2*k+2, n);
        }
        i
    }

    pub fn index(&mut self) {
        if self.starts.is_empty() {
            return;
        }
        self.sort_intervals();
        self.eytz.resize(self.starts.len() + 1, 0);
        self.eytz_index.resize(self.starts.len() + 1, 0);
        self.eytzinger_helper(0, 0, self.starts.len());

        self.branch.resize(self.starts.len(), usize::MAX);
        for i in 0..self.ends.len() - 1 {
            for j in i + 1..self.ends.len() {
                if self.ends[j] >= self.ends[i] {
                    break;
                }
                self.branch[j] = i;
            }
        }
        self.idx = 0;
    }

    pub fn upper_bound(&mut self, x: i32) {
        let mut i = 0;
        let n_intervals = self.starts.len();
        let mut best_idx = n_intervals;
        unsafe {
            while i < n_intervals {
                if *self.eytz.get_unchecked(i) > x {
                    if best_idx == n_intervals || *self.eytz.get_unchecked(i) <= *self.eytz.get_unchecked(best_idx) {
                        best_idx = i;  // best candidate closer to x
                    }
                    i = 2 * i + 1;
                } else {
                    i = 2 * i + 2;
                }
            }
            self.idx = if best_idx < n_intervals { *self.eytz_index.get_unchecked(best_idx) } else { n_intervals - 1 };
            if self.idx > 0 && *self.starts.get_unchecked(self.idx) > x {
                self.idx -= 1;
            }
        }
    }

    pub fn find_overlaps(&mut self, start: i32, end: i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(end);
        let mut i = self.idx;
        unsafe {  // remove bounds checking
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
                const SIMD_WIDTH: usize = 8; // 256 bits / 32 bits
                const BLOCK: usize = SIMD_WIDTH * 4;

                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i >= BLOCK {
                            let mut count: usize = 0;
                            for j in (i - BLOCK + 1..=i).rev().step_by(SIMD_WIDTH) {
                                let ends_vec = _mm256_loadu_si256(self.ends.get_unchecked(j - SIMD_WIDTH + 1) as *const i32 as *const __m256i);
                                let cmp_mask = _mm256_cmpgt_epi32(start_vec, ends_vec);
                                let mask = _mm256_movemask_epi8(_mm256_xor_si256(cmp_mask, _mm256_set1_epi32(-1)));
                                count += (_popcnt32(mask) / 4) as usize;
                            }
                            found += count;
                            i -= BLOCK;
                            if count < BLOCK {
                                break;
                            }
                        }
                    } else {
                        let branch_i = *self.branch.get_unchecked(i);
                        if branch_i >= i {
                            break;
                        }
                        i = branch_i;
                    }
                }
            }

            #[cfg(target_arch = "aarch64")]
            {
                use std::arch::aarch64::*;
                let start_vec = vdupq_n_s32(start);
                const SIMD_WIDTH: usize = 4; // 128 bits / 32 bits
                const BLOCK: usize = SIMD_WIDTH * 4;

                while i > 0 {
                    if start <= *self.ends.get_unchecked(i) {
                        found += 1;
                        i -= 1;
                        while i >= BLOCK {
                            let mut count = 0;
                            for j in (i - BLOCK + 1..=i).step_by(SIMD_WIDTH) {
                                let ends_vec = vld1q_s32(self.ends.get_unchecked(j - SIMD_WIDTH + 1) as *const i32);
                                let mask = vcleq_s32(start_vec, ends_vec);
                                count += vaddvq_u32(vreinterpretq_u32_s32(mask)) as usize;
                            }
                            found += count;
                            i -= BLOCK;
                            if count < BLOCK {
                                break;
                            }
                        }
                    } else {
                        let branch_i = *self.branch.get_unchecked(i);
                        if branch_i >= i {
                            break;
                        }
                        i = branch_i;
                    }
                }
            }

            // Process remaining elements
            while i > 0 {
                if start <= *self.ends.get_unchecked(i) {
                    found += 1;
                    i -= 1;
                } else {
                    let branch_i = *self.branch.get_unchecked(i);
                    if branch_i >= i {
                        break;
                    }
                    i = branch_i;
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
        let range_size = end_i - start_i;
        let mut tmp: Vec<Interval<T>> = Vec::with_capacity(range_size);
        for i in 0..range_size {
            tmp.push(Interval {
                start: self.starts[start_i + i],
                end: self.ends[start_i + i],
                data: self.data[start_i + i].clone(),
            });
        }
        tmp.sort_by(compare);
        for i in 0..range_size {
            self.starts[start_i + i] = tmp[i].start;
            self.ends[start_i + i] = tmp[i].end;
            self.data[start_i + i] = tmp[i].data.clone();
        }
    }
}