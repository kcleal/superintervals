use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::str;
use std::time::Instant;
extern crate fnv;
use fnv::FnvHashMap;
use clap::Parser;
extern crate libc;

use superintervals::SuperIntervals;
use superintervals::SuperIntervalsEytz;

use bincode;

// Define a trait that all SuperIntervals subclasses implement
pub trait IntervalCollection<T: Clone> {
    fn new() -> Self;
    fn add(&mut self, start: i32, end: i32, value: T);
    fn index(&mut self);
}

// Implement the trait for SuperIntervals
impl<T: Clone> IntervalCollection<T> for SuperIntervals<T> {
    fn new() -> Self {
        SuperIntervals::new()
    }
    fn add(&mut self, start: i32, end: i32, value: T) {
        self.add(start, end, value);
    }
    fn index(&mut self) {
        self.index();
    }
}

// Implement the trait for SuperIntervalsEytz
impl<T: Clone> IntervalCollection<T> for SuperIntervalsEytz<T> {
    fn new() -> Self {
        SuperIntervalsEytz::new()
    }
    fn add(&mut self, start: i32, end: i32, value: T) {
        self.add(start, end, value);
    }
    fn index(&mut self) {
        self.index();
    }
}


type GenericError = Box<dyn Error>;

// Parse a i32 with no checking whatsoever. (e.g. non-number characters will just)
fn i32_from_bytes_uncheckd(s: &[u8]) -> i32 {
    if s.is_empty() {
        0
    } else if s[0] == b'-' {
        -s[1..].iter().fold(0, |a, b| a * 10 + (b & 0x0f) as i32)
    } else {
        s.iter().fold(0, |a, b| a * 10 + (b & 0x0f) as i32)
    }
}

fn parse_bed_line(line: &[u8]) -> (&str, i32, i32) {
    let n = line.len() - 1;
    let mut p = 0;
    for c in &line[p..n] {
        if *c == b'\t' {
            break;
        }
        p += 1;
    }
    let seqname = unsafe { str::from_utf8_unchecked(&line[..p]) };
    p += 1;
    let p0 = p;

    for c in &line[p..n] {
        if *c == b'\t' {
            break;
        }
        p += 1;
    }
    let first = i32_from_bytes_uncheckd(&line[p0..p]);
    p += 1;
    let p0 = p;

    for c in &line[p..n] {
        if *c == b'\t' {
            break;
        }
        p += 1;
    }
    let last = i32_from_bytes_uncheckd(&line[p0..p]) - 1;

    (seqname, first, last)
}


// Read a bed file into a SuperIntervals
fn read_bed_file<I: IntervalCollection<()>>(path: &str) -> Result<FnvHashMap<String, I>, GenericError> {
    let mut nodes = FnvHashMap::default();
    let file = File::open(path)?;
    let mut rdr = BufReader::new(file);
    let mut line = Vec::new();
    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first, last) = parse_bed_line(&line);
        if seqname != "chr1" {
            continue;
        }
        let intervals = nodes.entry(seqname.to_string()).or_insert_with(I::new);
        intervals.add(first, last, ());
        line.clear();
    }
    let now = Instant::now();
    for intervals in nodes.values_mut() {
        intervals.index();
    }
    eprint!("{},", now.elapsed().as_micros());
    std::io::stderr().flush().unwrap();
    Ok(nodes)
}


fn query_bed_files(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
    let file = File::open(filename_b)?;
    let mut rdr = BufReader::new(file);
    let mut ranges: Vec<(i32, i32)> = Vec::new();
    let mut line = Vec::new();
    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (chrom, first, last) = parse_bed_line(&line);
        if chrom != "chr1" {
            continue;
        }
        ranges.push((first, last));
        line.clear();
    }
    //
    eprint!("SuperIntervals-rs,");
    let mut trees: FnvHashMap<String, SuperIntervals<()>> = read_bed_file::<SuperIntervals<()>>(filename_a)?;
    let intervals: &mut SuperIntervals<()> = trees.get_mut("chr1").ok_or("Chromosome intervals not found")?;

    // Verify the deserialized data
    let serialized_size = bincode::serialized_size(&intervals).unwrap();
    assert_ne!(serialized_size, 0);

    // Find overlaps (collecting results)
    let mut total_found = 0;
    let mut results = Vec::new();
    results.reserve(10000);
    let mut now = Instant::now();
    for &(first, last) in &ranges {
        intervals.find_overlaps(first, last, &mut results);
        total_found += results.len();
        results.clear();
    }
    eprint!("{},{},", now.elapsed().as_micros(), total_found);
    std::io::stderr().flush().unwrap();

    // Count overlaps
    let mut n_overlaps = 0;
    now = Instant::now();
    for &(first, last) in &ranges {
        n_overlaps += intervals.count_overlaps(first, last);
    }
    eprint!("{},{}\n", now.elapsed().as_micros(), n_overlaps);
    std::io::stderr().flush().unwrap();

    //
    eprint!("SuperIntervalsEytz-rs,");
    let mut trees2: FnvHashMap<String, SuperIntervalsEytz<()>> = read_bed_file::<SuperIntervalsEytz<()>>(filename_a)?;
    let intervals2: &mut SuperIntervalsEytz<()> = trees2.get_mut("chr1").ok_or("Chromosome intervals not found")?;

    // Find overlaps (collecting results)
    total_found = 0;
    results.clear();
    results.reserve(10000);
    now = Instant::now();
    for &(first, last) in &ranges {
        intervals2.find_overlaps(first, last, &mut results);
        total_found += results.len();
        results.clear();
    }
    eprint!("{},{},", now.elapsed().as_micros(), total_found);
    std::io::stderr().flush().unwrap();

    // Count overlaps
    n_overlaps = 0;
    now = Instant::now();
    for &(first, last) in &ranges {
        n_overlaps += intervals2.count_overlaps(first, last);
    }
    eprint!("{},{}\n", now.elapsed().as_micros(), n_overlaps);
    std::io::stderr().flush().unwrap();

    Ok(())
}


#[derive(Parser, Debug)]
#[command(about = " Find overlaps between two groups of intervals ")]
struct Args {
    /// intervals to index
    #[arg(value_name = "intervals.bed")]
    input1: String,
    /// query intervals
    #[arg(value_name = "queries.bed")]
    input2: String,

    /// compute proportion of queries covered
    #[arg(short = 'c', long)]
    coverage: bool,
}

fn query_bed_files_coverage(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
    let mut trees: FnvHashMap<String, SuperIntervals<()>> = read_bed_file::<SuperIntervals<()>>(filename_a)?;
//     let tree = read_bed_file(filename_a)?;

    let file = File::open(filename_b)?;
    let mut rdr = BufReader::new(file);
    let mut line = Vec::new();

    let mut total_count: usize = 0;
    let now = Instant::now();

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first, last) = parse_bed_line(&line);

        let mut cov: i32 = 0;
        let mut count: usize = 0;

        if let Some(seqname_tree) = trees.get_mut(seqname) {
            let countcov = seqname_tree.coverage(first, last);
            count = countcov.0;
            cov = countcov.1;
        }

        unsafe {
            let linelen = line.len();
            line[linelen - 1] = b'\0';
            libc::printf(
                b"%s\t%u\t%u\n\0".as_ptr() as *const libc::c_char,
                line.as_ptr() as *const libc::c_char,
                count as u32,
                cov,
            );
        }

        total_count += count;

        line.clear();
    }

    eprintln!("overlap: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    eprintln!("With coverage func total overlaps: {}", total_count);

    Ok(())
}

fn main() {
    let matches = Args::parse();
    let input1 = matches.input1.as_str();
    let input2 = matches.input2.as_str();
    let result;
    if matches.coverage {
        result = query_bed_files_coverage(input1, input2);
    } else {
        result = query_bed_files(input1, input2);
    }
    if let Err(err) = result {
        println!("error: {}", err)
    }
}
