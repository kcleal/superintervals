use coitrees::*;

use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
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
fn read_bed_file(path: &str) -> Result<FnvHashMap<String, COITree<(), u32>>, GenericError> {
    let mut nodes = IntervalHashMap::default();

    let now = Instant::now();

    let file = File::open(path)?;
    let mut rdr = BufReader::new(file);
    let mut line_count = 0;
    let mut line = Vec::new();

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first, last) = parse_bed_line(&line);

        let node_arr = if let Some(node_arr) = nodes.get_mut(seqname) {
            node_arr
        } else {
            nodes.entry(seqname.to_string()).or_insert(Vec::new())
        };

        node_arr.push(Interval::new(first, last, ()));

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
    let mut trees = FnvHashMap::<String, COITree<(), u32>>::default();
    for (seqname, seqname_nodes) in nodes {
        trees.insert(seqname, COITree::new(&seqname_nodes));
    }
    eprintln!("veb_order: {}s", now.elapsed().as_millis() as f64 / 1000.0);

    Ok(trees)
}


fn query_bed_files(filename_a: &str, filename_b: &str) -> Result<(), GenericError> {
    let tree = read_bed_file(filename_a)?;

    let file = File::open(filename_b)?;
    let mut rdr = BufReader::new(file);

    let mut ranges: Vec<(i32, i32)> = Vec::new();
    let mut line = Vec::new();

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (_, first, last) = parse_bed_line(&line);
        ranges.push((first, last));
        line.clear();
    }

    let seqname_tree = tree.get("chr1").ok_or("Chromosome tree not found")?;

    // Count overlaps
    let now = Instant::now();
    let total_count: usize = ranges.iter()
        .map(|&(first, last)| seqname_tree.query_count(first, last))
        .sum();
    let count_elapsed = now.elapsed();
    println!("Total count: {}, Count Time taken: {:?}", total_count, count_elapsed);

    // Find overlaps (collecting results)
    let now = Instant::now();
    let mut total_found = 0;

    let mut results = Vec::new();
    results.reserve(10000);
    for &(first, last) in &ranges {
//         let mut found = 0;
        seqname_tree.query(first, last, |node| {
            results.push(node.metadata);
//                 found += 1;
        });
//         results.push(found);
        total_found += results.len();
        results.clear();

    }

    let find_elapsed = now.elapsed();
    println!("Total found: {}, Find Time taken: {:?}", total_found, find_elapsed);

    Ok(())
}

fn query_bed_files_with_sorted_querent(
    filename_a: &str,
    filename_b: &str,
) -> Result<(), GenericError> {
    let trees = read_bed_file(filename_a)?;

    let file = File::open(filename_b)?;
    let mut rdr = BufReader::new(file);
    let mut line = Vec::new();

    let mut total_count: usize = 0;
    let now = Instant::now();

    let mut querents = FnvHashMap::<String, COITreeSortedQuerent<(), u32>>::default();
    for (seqname, tree) in &trees {
        querents.insert(seqname.clone(), COITreeSortedQuerent::new(tree));
    }

    while rdr.read_until(b'\n', &mut line).unwrap() > 0 {
        let (seqname, first, last) = parse_bed_line(&line);

        let mut count: usize = 0;
        if let Some(querent) = querents.get_mut(seqname) {
            querent.query(first, last, |_| count += 1);
        }

        // unfortunately printing in c is quite a bit faster than rust
        unsafe {
            let linelen = line.len();
            line[linelen - 1] = b'\0';
            libc::printf(
                b"%s\t%u\n\0".as_ptr() as *const libc::c_char,
                line.as_ptr() as *const libc::c_char,
                count as u32,
            );
        }

        total_count += count;

        line.clear();
    }

    eprintln!("overlap: {}s", now.elapsed().as_millis() as f64 / 1000.0);
    eprintln!("Sorted func total overlaps: {}", total_count);

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

    /// load both interval sets into memory instead of streaming queries
    #[arg(short = 't')]
    tree_vs_tree: bool,

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
