use coitrees::{COITree, Interval, IntervalTree};
use std::{fs, hint::black_box, ops::ControlFlow, time::Instant};
use superintervals::IntervalMap;
fn load(path: &str) -> Vec<(i32,i32)> {
    fs::read_to_string(path).unwrap().lines().map(|line| {
        let mut cols=line.split('\t'); assert_eq!(cols.next(),Some("chr1"));
        let a=cols.next().unwrap().parse().unwrap();
        let b:i32=cols.next().unwrap().parse().unwrap(); (a,b-1)
    }).collect()
}
// Independent weighted endpoint oracle: sum payload IDs without enumerating hits.
fn expected_checksum(refs: &[(i32,i32)], queries: &[(i32,i32)]) -> u64 {
    let mut starts:Vec<_>=refs.iter().enumerate().map(|(i,&(a,_))|(a,i as u64)).collect();
    let mut ends:Vec<_>=refs.iter().enumerate().map(|(i,&(_,b))|(b,i as u64)).collect();
    starts.sort_unstable(); ends.sort_unstable();
    let prefix=|xs:&[(i32,u64)]| { let mut p=vec![0u64]; for &(_,id) in xs { p.push(p.last().unwrap().wrapping_add(id)); } p };
    let sp=prefix(&starts); let ep=prefix(&ends);
    queries.iter().fold(0u64,|sum,&(a,b)|sum.wrapping_add(
        sp[starts.partition_point(|&(x,_)|x<=b)].wrapping_sub(ep[ends.partition_point(|&(x,_)|x<a)])))
}
fn main() {
    let args:Vec<String>=std::env::args().collect();
    let refs=load(&args[2]); let queries=load(&args[3]);
    let mut results:Vec<u32>=Vec::with_capacity(10000);
    if args[1]=="si" {
        let now=Instant::now(); let mut index=IntervalMap::new();
        for (i,&(a,b)) in refs.iter().enumerate() { index.add(a,b,i as u32); }
        index.build(); let build=now.elapsed().as_micros();
        let now=Instant::now(); let mut hits=0usize;
        for &(a,b) in &queries { index.search_values(a,b,&mut results); hits+=results.len(); black_box(&results); results.clear(); }
        let find=now.elapsed().as_micros();
        let now=Instant::now(); let count:usize=queries.iter().map(|&(a,b)|index.count(a,b)).sum();
        let count_time=now.elapsed().as_micros();
        let now=Instant::now(); let mut stream_hits=0usize; let mut checksum=0u64;
        for &(a,b) in &queries { for id in index.search_values_iter(a,b) { stream_hits+=1; checksum=checksum.wrapping_add(*id as u64); } }
        black_box((stream_hits,checksum)); let stream=now.elapsed().as_micros();
        assert_eq!(checksum,expected_checksum(&refs,&queries),"streamed payload checksum");
        assert_eq!(stream_hits,hits,"streamed hit total");
        println!("SuperIntervals-rs,{build},{find},{hits},{count_time},{count},{stream},{stream_hits}");
    } else if args[1]=="lapper" {
        // Lapper uses unsigned, half-open coordinates; load() returns closed ends.
        let now=Instant::now();
        let nodes:Vec<rust_lapper::Interval<u32,u32>>=refs.iter().enumerate().map(|(i,&(a,b))|
            rust_lapper::Interval { start:a as u32, stop:(b+1) as u32, val:i as u32 }).collect();
        let index=rust_lapper::Lapper::new(nodes);
        let build=now.elapsed().as_micros();
        let now=Instant::now(); let mut hits=0usize;
        for &(a,b) in &queries {
            for node in index.find(a as u32,(b+1) as u32) { results.push(node.val); }
            hits+=results.len(); black_box(&results); results.clear();
        }
        let find=now.elapsed().as_micros();
        let now=Instant::now(); let count:usize=queries.iter().map(|&(a,b)|index.count(a as u32,(b+1) as u32)).sum();
        let count_time=now.elapsed().as_micros();
        let now=Instant::now(); let mut stream_hits=0usize; let mut checksum=0u64;
        for &(a,b) in &queries { for node in index.find(a as u32,(b+1) as u32) { stream_hits+=1; checksum=checksum.wrapping_add(node.val as u64); } }
        black_box((stream_hits,checksum)); let stream=now.elapsed().as_micros();
        assert_eq!(checksum,expected_checksum(&refs,&queries),"streamed payload checksum");
        assert_eq!(stream_hits,hits,"streamed hit total");
        println!("rust-lapper,{build},{find},{hits},{count_time},{count},{stream},{stream_hits}");
    } else if args[1]=="sif" {
        // Use half-open ranges, matching the original BED inputs.
        let now=Instant::now();
        let tree=sif_itree::ITree::<i32,u32>::new(refs.iter().enumerate()
            .map(|(i,&(a,b))| (a..b+1,i as u32)));
        let build=now.elapsed().as_micros();
        let now=Instant::now(); let mut hits=0usize;
        for &(a,b) in &queries {
            let _=tree.query(a..b+1,|(_,id)| {
                results.push(*id);
                ControlFlow::<()>::Continue(())
            });
            hits+=results.len(); black_box(&results); results.clear();
        }
        let find=now.elapsed().as_micros();
        let now=Instant::now(); let mut stream_hits=0usize; let mut checksum=0u64;
        for &(a,b) in &queries {
            let _=tree.query(a..b+1,|(_,id)| {
                stream_hits+=1; checksum=checksum.wrapping_add(*id as u64);
                ControlFlow::<()>::Continue(())
            });
        }
        black_box((stream_hits,checksum)); let stream=now.elapsed().as_micros();
        assert_eq!(checksum,expected_checksum(&refs,&queries),"streamed payload checksum");
        assert_eq!(stream_hits,hits,"streamed hit total");
        // No specialised count API; do not substitute callback enumeration.
        println!("sif-itree,{build},{find},{hits},NA,NA,{stream},{stream_hits}");
    } else {
        assert_eq!(args[1],"coitrees");
        let now=Instant::now();
        let nodes:Vec<Interval<u32>>=refs.iter().enumerate().map(|(i,&(a,b))|Interval::new(a,b,i as u32)).collect();
        let tree=COITree::<u32,u32>::new(&nodes);
        let build=now.elapsed().as_micros();
        let now=Instant::now(); let mut hits=0usize;
        for &(a,b) in &queries {
            tree.query(a,b,|node|results.push(*node.metadata));
            hits+=results.len(); black_box(&results); results.clear();
        }
        let find=now.elapsed().as_micros();
        let now=Instant::now(); let count:usize=queries.iter().map(|&(a,b)|tree.query_count(a,b)).sum();
        let count_time=now.elapsed().as_micros();
        let now=Instant::now(); let mut stream_hits=0usize; let mut checksum=0u64;
        for &(a,b) in &queries { tree.query(a,b,|node| { stream_hits+=1; checksum=checksum.wrapping_add(*node.metadata as u64); }); }
        black_box((stream_hits,checksum)); let stream=now.elapsed().as_micros();
        assert_eq!(checksum,expected_checksum(&refs,&queries),"streamed payload checksum");
        assert_eq!(stream_hits,hits,"streamed hit total");
        println!("Coitrees,{build},{find},{hits},{count_time},{count},{stream},{stream_hits}");
    }
}
