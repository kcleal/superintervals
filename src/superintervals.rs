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
    branch: Vec<Option<usize>>,
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

        self.branch.resize(self.starts.len(), None);
        for i in 0..self.ends.len() - 1 {
            for j in i + 1..self.ends.len() {
                if self.ends[j] >= self.ends[i] {
                    break;
                }
                self.branch[j] = Some(i);
            }
        }

        self.idx = 0;
    }

    fn upper_bound(&mut self, value: &i32) {
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
                    if *self.starts.get_unchecked(self.idx + half) <= *value {
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
                if *self.starts.get_unchecked(self.idx + half) <= *value {
                    self.idx += length - half;
                }
            }

            length = half;
        }

        // Adjust index if necessary
        if self.idx > 0 && (self.idx == self.starts.len() || self.starts[self.idx] > *value) {
            self.idx -= 1;
        }
    }

    pub fn find_overlaps(&mut self, start: &i32, end: &i32, found: &mut Vec<T>) {
        if self.starts.is_empty() {
            return;
        }
        self.upper_bound(end);
        let mut i = self.idx;

        while i > 0 {
            // Using unsafe block to remove bounds checking
            unsafe {
                if *start <= *self.ends.get_unchecked(i) {
                    found.push(self.data.get_unchecked(i).clone());
                    i -= 1;
                } else {
                    if let Some(branch_i) = self.branch.get_unchecked(i) {
                        if *branch_i >= i {
                            break;
                        }
                        i = *branch_i;
                    } else {
                        break;
                    }
                }
            }
        }

        // Check the first element separately
        if i == 0 {
            unsafe {
                if *start <= *self.ends.get_unchecked(0) && *self.starts.get_unchecked(0) <= *end {
                    found.push(self.data.get_unchecked(0).clone());
                }
            }
        }
    }


//     fn upper_bound(&mut self, value: &i32) {
//         let mut length = self.starts.len();
//         if length == 0 {
//             return;
//         }
//         length -= 1;
//         self.idx = 0;
//         const ENTRIES_PER_256KB: usize = 256 * 1024 / std::mem::size_of::<i32>();
//
//         if length >= ENTRIES_PER_256KB {
//             let num_per_cache_line: usize = std::cmp::max(64 / std::mem::size_of::<i32>(), 1);
//             while length >= 3 * num_per_cache_line {
//                 let half = length / 2;
//                 // Simulate prefetching
//                 let _ = &self.starts[self.idx + half / 2];
//                 let first_half1 = self.idx + (length - half);
//                 let _ = &self.starts[first_half1 + half / 2];
//                 if self.starts[self.idx + half] <= *value {
//                     self.idx += length - half;
//                 }
//                 length = half;
//             }
//         }
//
//         while length > 0 {
//             let half = length / 2;
//             if self.starts[self.idx + half] <= *value {
//                 self.idx += length - half;
//             }
//             length = half;
//         }
//
//         if self.idx > 0 && (self.idx == self.starts.len() || self.starts[self.idx] > *value) {
//             self.idx -= 1;
//         }
//     }

//     pub fn find_overlaps(&mut self, start: &i32, end: &i32, found: &mut Vec<T>) {
//         if self.starts.is_empty() {
//             return;
//         }
//         self.upper_bound(end);
//         let mut i = self.idx;
//         while i > 0 {
//             if *start <= self.ends[i] {
//                 found.push(self.data[i].clone());
//                 i -= 1;
//             } else {
//                 if let Some(branch_i) = self.branch[i] {
//                     if branch_i >= i {
//                         break;
//                     }
//                     i = branch_i;
//                 } else {
//                     break;
//                 }
//             }
//         }
//         if i == 0 && *start <= self.ends[0] && self.starts[0] <= *end {
//             found.push(self.data[0].clone());
//         }
//     }

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
