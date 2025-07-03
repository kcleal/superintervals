use coitrees::*;

use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::str;
use std::time::Instant;

extern crate fnv;
use fnv::FnvHashMap;

use clap::Parser;

extern crate libc;

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

type IntervalHashMap = FnvHashMap<String, Vec<Interval<()>>>;

// Read a bed file into a COITree
fn read_bed_file(path: &str, name: &str) -> Result<FnvHashMap<String, COITree<(), u32>>, GenericError> {
    let mut nodes = IntervalHashMap::default();
    let file = File::open(path)?;
    let mut rdr = BufReader::new(file);
    let mut line = Vec::new();
    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first, last) = parse_bed_line(&line);
        if seqname != "chr1" || last < 0 || last < first {
            line.clear();
            continue;
        }
        let node_arr = if let Some(node_arr) = nodes.get_mut(seqname) {
            node_arr
        } else {
            nodes.entry(seqname.to_string()).or_insert(Vec::new())
        };
        node_arr.push(Interval::new(first, last, ()));
        line.clear();
    }
    let now = Instant::now();
    let mut trees = FnvHashMap::<String, COITree<(), u32>>::default();
    for (seqname, seqname_nodes) in nodes {
        trees.insert(seqname, COITree::new(&seqname_nodes));
    }
    eprint!("{},{},", name, now.elapsed().as_micros());
    std::io::stderr().flush().unwrap();
    Ok(trees)
}


fn query_bed_files(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
    let tree = read_bed_file(filename_a, "Coitrees")?;
    let file = File::open(filename_b)?;
    let mut rdr = BufReader::new(file);
    let mut ranges: Vec<(i32, i32)> = Vec::new();
    let mut line = Vec::new();
    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (chrom, first, last) = parse_bed_line(&line);
        if chrom != "chr1" || last < 0 || last < first {
            line.clear();
            continue;
        }
        ranges.push((first, last));
        line.clear();
    }
    let seqname_tree = tree.get("chr1").ok_or("Chromosome tree not found")?;

    // Find overlaps (collecting results)
    let mut total_found = 0;
    let mut results = Vec::new();
    results.reserve(10000);
    let mut now = Instant::now();
    for &(first, last) in &ranges {
        seqname_tree.query(first, last, |node| {
            results.push(node.metadata);
        });
        total_found += results.len();
        results.clear();
    }
    eprint!("{},{},", now.elapsed().as_micros(), total_found);
    std::io::stderr().flush().unwrap();

    // Count overlaps
    now = Instant::now();
    let total_count: usize = ranges.iter()
        .map(|&(first, last)| seqname_tree.query_count(first, last))
        .sum();
    eprint!("{},{}\n", now.elapsed().as_micros(), total_count);
    std::io::stderr().flush().unwrap();

    Ok(())
}


fn query_bed_files_with_sorted_querent(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
    let trees = read_bed_file(filename_a, "Coitrees-s")?;

    let mut querents = FnvHashMap::<String, COITreeSortedQuerent<(), u32>>::default();
    for (seqname, tree) in &trees {
        querents.insert(seqname.clone(), COITreeSortedQuerent::new(tree));
    }

    let file = File::open(filename_b)?;
    let mut rdr = BufReader::new(file);
    let mut ranges: Vec<(i32, i32)> = Vec::new();
    let mut line = Vec::new();
    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (_, first, last) = parse_bed_line(&line);
        ranges.push((first, last));
        line.clear();
    }

    // Use `get_mut` to get a mutable reference
    let seqname_tree = querents.get_mut("chr1").ok_or("Chromosome tree not found")?;

    // Find overlaps (collecting results)
    let mut total_found = 0;
    let mut results = Vec::new();
    results.reserve(10000);
    let now = Instant::now();
    for &(first, last) in &ranges {
        seqname_tree.query(first, last, |node| {
            results.push(node.metadata);
        });
        total_found += results.len();
        results.clear();
    }
    eprint!("{},{}\n", now.elapsed().as_micros(), total_found);
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

    /// use alternative search strategy that's faster if queries are sorted and tend to overlap
    #[arg(short = 's', long = "sorted")]
    use_sorted_querent: bool,

}

fn main() {
    let matches = Args::parse();

    let input1 = matches.input1.as_str();
    let input2 = matches.input2.as_str();

    let result;

    if matches.use_sorted_querent {
        result = query_bed_files_with_sorted_querent(input1, input2);
    } else {
        result = query_bed_files(input1, input2);
    }
    if let Err(err) = result {
        println!("error: {}", err)
    }
}
