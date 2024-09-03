use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str;
use std::time::Instant;

extern crate fnv;
use fnv::FnvHashMap;

use clap::Parser;

extern crate libc;

use superintervals::superintervals::SuperIntervals;


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

type IntervalHashMap = FnvHashMap<String, SuperIntervals< ()>>;

// Read a bed file into a SuperIntervals
fn read_bed_file(path: &str) -> Result<FnvHashMap<String, SuperIntervals< ()>>, GenericError> {
    let mut nodes = IntervalHashMap::default();

    let now = Instant::now();

    let file = File::open(path)?;
    let mut rdr = BufReader::new(file);
    let mut line_count = 0;
    let mut line = Vec::new();

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first, last) = parse_bed_line(&line);

        let intervals = nodes.entry(seqname.to_string()).or_insert_with(SuperIntervals::new);
        intervals.add(first, last, ());

        line_count += 1;
        line.clear();
    }
    eprintln!(
        "reading bed: {}s",
        now.elapsed().as_millis() as f64 / 1000.0
    );
    eprintln!("lines: {}", line_count);
    eprintln!("sequences: {}", nodes.len());

    let now = Instant::now();
    for intervals in nodes.values_mut() {
        intervals.index();
    }
    eprintln!("indexing: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    Ok(nodes)
}

fn query_bed_files(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
    let mut trees = read_bed_file(filename_a)?;

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
    println!("N queries: {}", ranges.len());
    let intervals: &mut SuperIntervals<()> = trees.get_mut("chr1").ok_or("Chromosome intervals not found")?;

    // Find overlaps (collecting results)
    let now = Instant::now();
    let mut total_found = 0;

    let mut results = Vec::new();
    results.reserve(10000);
    for &(first, last) in &ranges {
        intervals.find_overlaps(&first, &last, &mut results);
        total_found += results.len();
        results.clear();
    }

    let find_elapsed = now.elapsed();
    println!("Total found: {}, Find Time taken: {:?}", total_found, find_elapsed);

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

}

fn main() {
    let matches = Args::parse();

    let input1 = matches.input1.as_str();
    let input2 = matches.input2.as_str();

    let result = query_bed_files(input1, input2);

    if let Err(err) = result {
        println!("error: {}", err)
    }
}
