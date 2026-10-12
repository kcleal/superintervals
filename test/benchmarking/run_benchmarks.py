#!/usr/bin/env python3
"""Build and run the manuscript's native interval-library benchmark. Python >=3.10."""
import argparse, bisect, csv, hashlib, json, os, platform, random
import shutil, statistics, subprocess, sys, tarfile, time, urllib.request, html, re, math
from pathlib import Path

URL='https://github.com/kcleal/superintervals/releases/download/v0.2.0/data.tar.gz'
DATA_SHA='c98159506689c412a909505e79576035b42ab50a4930d6c693ef9517756039de'
PAIRS=[('mito','DB53.chrM_reads_as_chr1.bed','DB53.chrM_reads_as_chr1.bed'),
 ('rna-anno','ex-rna.bed','ex-anno.bed'),('anno-rna','ex-anno.bed','ex-rna.bed'),
 ('genes-DB53','ucsc.hg19.genes.bed','chr1.DB53.bed'),
 ('DB53-genes','chr1.DB53.bed','ucsc.hg19.genes.bed'),
 ('ONT-DB53','PAO33946.chr1.bed','chr1.DB53.bed'),
 ('DB53-ONT','chr1.DB53.bed','PAO33946.chr1.bed'),
 ('fixed-width',None,None)]
HERE=Path(__file__).resolve().parent

# Display reciprocal directions together without changing benchmark seeds/order.
DISPLAY_PAIRS=['DB53-ONT','ONT-DB53','DB53-genes','genes-DB53','rna-anno','anno-rna','mito','fixed-width']

def workload_key(workload):
    dataset,order=workload
    return (DISPLAY_PAIRS.index(dataset) if dataset in DISPLAY_PAIRS else len(DISPLAY_PAIRS),
            dataset,order!='sorted')


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def command(args, cwd=None, env=None, log=None, timeout=1800):
    if log is not None: log.write('$ '+ ' '.join(map(str,args))+'\n'); log.flush()
    p=subprocess.run(list(map(str,args)),cwd=cwd,env=env,stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT,timeout=timeout)
    output=p.stdout.decode(errors='replace')
    if log is not None: log.write(output+'\n'); log.flush()
    if p.returncode: raise RuntimeError(f'Command failed ({p.returncode}): {args}\n{output[-5000:]}')
    return output

def prepare_source(repo, out):
    # Capture only benchmark-required source trees; exclude build products and local data.
    dst=out/'source'; dst.mkdir()
    for name in ('src','Cargo.toml'):
        p=repo/name
        if p.is_dir(): shutil.copytree(p,dst/name,ignore=shutil.ignore_patterns('__pycache__','*.so','*.o','*.pyc'))
        else: shutil.copy2(p,dst/name)
    (dst/'test'/'3rd-party').mkdir(parents=True)
    shutil.copy2(repo/'test'/'bench.cpp',dst/'test'/'bench.cpp')
    for p in (repo/'test'/'3rd-party').iterdir():
        if p.is_file(): shutil.copy2(p,dst/'test'/'3rd-party'/p.name)
    # Empty tests/examples satisfy root Cargo manifest; the independent driver uses the library.
    (dst/'test'/'benchmarking').mkdir()
    for name in ('run_benchmarks.py','rust_driver.rs','Cargo.lock','make_site.py'):
        shutil.copy2(HERE/name,dst/'test'/'benchmarking'/name)
    (dst/'examples').mkdir()
    for p in (repo/'examples').glob('*.rs'): shutil.copy2(p,dst/'examples'/p.name)
    for p in (repo/'test').glob('*.rs'): shutil.copy2(p,dst/'test'/p.name)
    cpp=dst/'test'/'bench.cpp'; text=cpp.read_text().replace('high_resolution_clock','steady_clock')
    # Existing driver already converts BED endpoints for end-inclusive libraries.
    # Move the IntervalTree timer ahead of its interval-container preparation.
    old='    std::vector<interval_tree::Interval<int, int>> intervals2;'
    if text.count(old)!=1: raise RuntimeError('Unsupported C++ driver: IntervalTree setup changed')
    text=text.replace(old,'    t0 = steady_clock::now();\n'+old)
    text=text.replace('    t0 = steady_clock::now();\n    interval_tree::IntervalTree<int, int> tree2',
                      '    interval_tree::IntervalTree<int, int> tree2')
    # One native tool per process avoids carry-over allocations between tools.
    old='    Bench::run_tools(intervals, queries);'
    new='''    std::string tool = argc > 3 ? argv[3] : "all";
    if (tool == "iitree") Bench::run_IITree(intervals, queries);
    else if (tool == "intervaltree") Bench::run_ITree(intervals, queries);
    else if (tool == "ncls") Bench::run_NCLS(intervals, queries);
    else if (tool == "si") {
        si::IntervalMap<int, Bench::DType> index;
        Bench::run_SuperIntervals(intervals, queries, index, "SuperIntervals-C++");
    } else return 2;'''
    if text.count(old)!=1: raise RuntimeError('Unsupported C++ driver: main changed')
    text=text.replace(old,new)
    # Remove large-query variants from the captured legacy C++ driver.
    for method in ('search_values_large','count_large'):
        pattern=r'    found = 0;\n    t1 = steady_clock::now\(\);\n    for \(const auto& item : queries\) \{\n(?:(?!    found = 0;).)*?itv\.'+method+r'\(.*?    std::cerr << uSec\(t1\).*?;[^\n]*\n'
        text,n=re.subn(pattern,'',text,flags=re.S)
        if n!=1: raise RuntimeError('Unsupported C++ driver: '+method+' setup changed')
    # Stream actual payloads through the C++ borrowed iterator, matching Rust's
    # checksum consumer. Verify payloads using an independent endpoint oracle.
    count_line = 'std::cerr << uSec(t1) << "," << found << ",";  // count all overlapping'
    stream_code = r'''std::cerr << uSec(t1) << "," << found << ",";  // count all overlapping
    found = 0;
    uint64_t checksum = 0;
    t1 = steady_clock::now();
    for (const auto& item : queries) {
        for (const auto& value : itv.search_values(item.start, item.end - 1)) {
            ++found; checksum += static_cast<uint64_t>(value);
        }
    }
    const auto stream_us = uSec(t1);
    std::vector<std::pair<int, uint64_t>> begins, stops;
    begins.reserve(intervals.size()); stops.reserve(intervals.size());
    for (size_t id = 0; id < intervals.size(); ++id) {
        begins.emplace_back(intervals[id].start, id);
        stops.emplace_back(intervals[id].end, id);
    }
    std::sort(begins.begin(), begins.end()); std::sort(stops.begin(), stops.end());
    std::vector<uint64_t> bp(begins.size()+1, 0), ep(stops.size()+1, 0);
    for (size_t id = 0; id < begins.size(); ++id) {
        bp[id+1] = bp[id] + begins[id].second; ep[id+1] = ep[id] + stops[id].second;
    }
    uint64_t expected = 0;
    for (const auto& item : queries) {
        const size_t lo = std::lower_bound(begins.begin(), begins.end(), std::make_pair(item.end, uint64_t{0})) - begins.begin();
        const size_t hi = std::upper_bound(stops.begin(), stops.end(), std::make_pair(item.start, UINT64_MAX)) - stops.begin();
        expected += bp[lo] - ep[hi];
    }
    if (checksum != expected) throw std::runtime_error("streamed payload checksum mismatch");
    std::cerr << stream_us << "," << found << "\n";
'''
    if text.count(count_line) != 1: raise RuntimeError('Unsupported C++ driver: count output changed')
    text = '#include <stdexcept>\n' + text.replace(count_line, stream_code)
    cpp.write_text(text)
    return dst

def build(src,out,log):
    cc=os.environ.get('CC','cc'); cxx=os.environ.get('CXX','c++')
    flags=['-O3','-march=native']
    objects=[]
    for name in ('intervaldb','cgranges'):
        obj=out/(name+'.o'); objects.append(obj)
        command([cc,*flags,'-c',src/'test'/'3rd-party'/(name+'.c'),'-o',obj],log=log)
    cpp=out/'run-cpp-libs'
    command([cxx,*flags,'-std=c++17','-I'+str(src/'src'),'-I'+str(src/'test'/'3rd-party'),
             src/'test'/'bench.cpp',*objects,'-o',cpp],log=log)
    driver=out/'rust-driver'; (driver/'src').mkdir(parents=True)
    shutil.copy2(HERE/'rust_driver.rs',driver/'src'/'main.rs')
    manifest='''[package]
name = "manuscript-benchmark"
version = "0.1.0"
edition = "2021"
[dependencies]
superintervals = { path = "../source" }
coitrees = "=0.4.0"
rust-lapper = "=1.3.0"
[profile.release]
opt-level = 3
lto = true
codegen-units = 1
'''
    (driver/'Cargo.toml').write_text(manifest)
    shutil.copy2(HERE/'Cargo.lock',driver/'Cargo.lock')
    env=os.environ.copy(); env['RUSTFLAGS']='-C target-cpu=native'
    env['CARGO_TARGET_DIR']=str(out/'rust-target')
    command(['cargo','build','--release','--locked','--manifest-path',driver/'Cargo.toml'],env=env,log=log)
    resolved=json.loads(command(['cargo','metadata','--locked','--format-version','1','--manifest-path',driver/'Cargo.toml'],env=env,log=log))
    (out/'cargo-metadata.json').write_text(json.dumps(resolved,indent=2)+'\n')
    return cpp,out/'rust-target'/'release'/'manuscript-benchmark'

def data(archive,out,log):
    archive=Path(archive) if archive else out/'data.tar.gz'
    if not archive.exists():
        print('Downloading release data',flush=True)
        urllib.request.urlretrieve(URL,archive)
    digest=sha(archive)
    if digest!=DATA_SHA: raise RuntimeError(f'Data archive SHA256 mismatch: {digest}')
    dest=out/'datasets'; dest.mkdir()
    wanted={f'data/{name}' for _,a,b in PAIRS for name in (a,b) if name is not None}
    with tarfile.open(archive,'r:gz') as tar:
        for name in sorted(wanted):
            member=tar.getmember(name)
            if not member.isfile(): raise RuntimeError(f'Unexpected archive member: {name}')
            with tar.extractfile(member) as source, open(dest/Path(name).name,'wb') as target:
                shutil.copyfileobj(source,target)
    return dest,digest

def read_bed(path,limit=None,fraction=0.25,seed=1729,window=0):
    # A filename-specific stream gives reciprocal comparisons the same sample.
    stream_seed=int.from_bytes(hashlib.sha256(f'{seed}:{path.name}'.encode()).digest()[:8],'little')
    sample_rng=random.Random(stream_seed); limit_rng=random.Random(stream_seed ^ 0x9e3779b97f4a7c15)
    rows=[]; source_n=eligible_n=selected_n=0
    with open(path) as f:
        for i,line in enumerate(f,1):
            if line.startswith(('#','track','browser')): continue
            cols=line.rstrip('\n').split('\t')
            if cols[0]!='chr1': continue
            try: a,b=map(int,cols[1:3])
            except (ValueError,TypeError): raise ValueError(f'{path}:{i}: invalid BED record')
            if not 0<=a<b<=2147483647: raise ValueError(f'{path}:{i}: invalid BED range')
            source_n+=1
            if window and a>=window: continue
            eligible_n+=1
            if fraction<1 and sample_rng.random()>=fraction: continue
            selected_n+=1
            if not limit or len(rows)<limit: rows.append((a,b))
            else:
                j=limit_rng.randrange(selected_n)
                if j<limit: rows[j]=(a,b)
    if not rows: raise ValueError(f'No chr1 records selected from {path}')
    return rows,{'source_n':source_n,'eligible_n':eligible_n,'selected_n':selected_n,
                 'retained_n':len(rows),'fraction':fraction,'seed':stream_seed,
                 'start_window_bp':window or None,'reservoir_limit':limit}


def write_bed(path,rows):
    with open(path,'w') as f:
        for a,b in rows: f.write(f'chr1\t{a}\t{b}\n')

def simulate_fixed_width(seed, limit=None, fraction=0.25):
    width=150
    reference_n=min(max(1,int(1_000_000*fraction)), limit or 1_000_000)
    query_n=min(max(1,int(100_000*fraction)), limit or 100_000)
    # Scale the window with the reference count to retain ~10x coverage in smoke runs.
    window=max(width, reference_n*width//10)
    def intervals(n, stream_seed):
        rng=random.Random(stream_seed)
        return [(a,a+width) for a in (rng.randrange(window-width+1) for _ in range(n))]
    parameters={'width_bp':width,'window_bp':window,'reference_n':reference_n,
                'query_n':query_n,'reference_seed':seed,'query_seed':seed+1,
                'distribution':'uniform integer starts, sampled with replacement',
                'reference_coverage':reference_n*width/window,'fraction':fraction}
    return intervals(reference_n,seed),intervals(query_n,seed+1),parameters

def oracle(refs,queries):
    # O((n+q) log n), independent of branch links. BED intervals are nonempty.
    starts=sorted(a for a,b in refs); ends=sorted(b for a,b in refs)
    return sum(bisect.bisect_left(starts,b)-bisect.bisect_right(ends,a) for a,b in queries)

def parse(output):
    rows=[]
    for line in output.strip().splitlines():
        cols=line.split(',')
        if len(cols) not in (4,6,8): raise RuntimeError('Unexpected driver output: '+line)
        tool=cols[0]; vals=list(map(int,cols[1:])); build_us=vals[0]
        ops=['find','count','stream'][:(len(cols)-2)//2]
        for i,op in enumerate(ops): rows.append((tool,op,build_us,vals[1+2*i],vals[2+2*i]))
    return rows

def run_suite(cpp,rust,ref,query,order,seed,log,timeout,on_tool_done=None):
    jobs=[('iitree',[cpp,ref,query,'iitree']),('intervaltree',[cpp,ref,query,'intervaltree']),
          ('ncls',[cpp,ref,query,'ncls']),('si-cpp',[cpp,ref,query,'si']),
          ('si-rust',[rust,'si',ref,query]),('coitrees',[rust,'coitrees',ref,query]),
          ('rust-lapper',[rust,'lapper',ref,query])]
    random.Random(seed).shuffle(jobs)
    result=[]
    for name,cmd in jobs:
        if on_tool_done: print(f'  {name} running',flush=True)
        started=time.perf_counter()
        measured=parse(command(cmd,log=log,timeout=timeout))
        elapsed=time.perf_counter()-started
        result.extend(measured)
        if on_tool_done:
            on_tool_done(measured[0][0],elapsed)
            print(f'  {measured[0][0]}: {elapsed:.2f} s',flush=True)
    return result

def number(value):
    return f'{float(value):,.2f}'.rstrip('0').rstrip('.')


def cpu_fields(meta):
    fields = {x['field'].rstrip(':'): x['data'] for x in meta.get('lscpu', {}).get('lscpu', [])}
    model = fields.get('Model name') or fields.get('Model')
    if not model:
        for line in meta.get('cpu', '').splitlines():
            if line.lower().startswith(('model name', 'hardware')):
                model = line.split(':', 1)[-1].strip()
                break
    return fields, model or meta.get('cpu', 'unavailable').strip().split('\n')[0]


def timing_plot(rows, operation):
    selected=[r for r in rows if r['operation']==operation]
    preferred=['SuperIntervals-rs','SuperIntervals-C++','Coitrees','rust-lapper','ImplicitITree-C++','IntervalTree-C++','NCLS-C']
    tools=[tool for tool in preferred if any(r['tool']==tool for r in selected)]
    workloads=sorted({(r['dataset'],r['order']) for r in selected},key=workload_key)
    left=180;cell=64;top=100;row_height=38;width=left+cell*len(workloads)+15;height=top+row_height*len(tools)+58
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="{operation} time relative to the fastest tool" width="{width}" height="{height}" style="display:block;font:12px sans-serif"><rect width="100%" height="100%" fill="white"/>',
           '<text x="12" y="22" font-size="15" font-weight="bold">Time / fastest · lower is better</text>',
           '<text x="12" y="43" fill="#555">Hover: median, Q1–Q3. Colour limit: 8×.</text>']
    for column,key in enumerate(workloads):
        x=left+column*cell;group=[r for r in selected if (r['dataset'],r['order'])==key]
        fastest=min(float(r['query_median_us']) for r in group)
        reference,query=('fixed width','150 bp') if key[0]=='fixed-width' else (key[0].rsplit("-",1) if "-" in key[0] else (key[0],"self"))
        direction='' if key[0]=='fixed-width' else ' →'
        parts.append(f'<text x="{x+cell/2}" y="60" text-anchor="middle">{html.escape(reference+direction)}</text><text x="{x+cell/2}" y="75" text-anchor="middle">{html.escape(query)}</text><text x="{x+cell/2}" y="91" text-anchor="middle" fill="#555">{key[1]}</text>')
        for row,tool in enumerate(tools):
            y=top+row*row_height;record=next((r for r in group if r['tool']==tool),None)
            label='—';fill='#eee';text='#333';title='Operation not provided'
            if record is not None and fastest>0:
                t=float(record['query_median_us']);ratio=t/fastest
                amount=min(1,max(0,math.log2(ratio)/3))
                rgb=[round(a+(b-a)*amount) for a,b in zip((224,243,239),(177,49,35))]
                fill='#'+''.join(f'{v:02x}' for v in rgb);text='#fff' if amount>.55 else '#222'
                label=f'{ratio:.2f}×' if ratio<10 else f'{ratio:.0f}×'
                title=f"{tool}: {t:g} µs [{record['query_q1_us']}, {record['query_q3_us']}]; {ratio:.2f}× fastest; {key[0]} {key[1]}"
            elif record is not None:
                title='Fastest time is zero at microsecond resolution; ratio unavailable'
            parts.append(f'<g><title>{html.escape(title)}</title><rect x="{x}" y="{y}" width="{cell-2}" height="{row_height-2}" fill="{fill}"/><text x="{x+cell/2-1}" y="{y+24}" text-anchor="middle" fill="{text}">{label}</text></g>')
    for row,tool in enumerate(tools):
        parts.append(f'<text x="12" y="{top+row*row_height+24}" font-size="13">{html.escape(tool)}</text>')
    parts.append(f'<text x="12" y="{height-18}" fill="#555">Absolute timings in expandable tables.</text></svg>')
    return ''.join(parts)


def tool_versions(repo, src, dependencies, revision='HEAD'):
    version=re.search(r'^version\s*=\s*"([^"]+)"', (src/'Cargo.toml').read_text(),re.M).group(1)
    result={}
    for tool,path in [('SuperIntervals-rs','src/superintervals.rs'),
                      ('SuperIntervals-C++','src/superintervals.hpp'),
                      ('ImplicitITree-C++','test/3rd-party/IITree.hpp'),
                      ('IntervalTree-C++','test/3rd-party/IntervalTree.h'),
                      ('NCLS-C','test/3rd-party/intervaldb.c')]:
        date=command(['git','log','-1','--format=%cs',revision,'--',path],cwd=repo).strip()
        captured_sha=sha(src/path)
        committed_sha=hashlib.sha256(command(['git','show',f'{revision}:{path}'],cwd=repo).encode()).hexdigest()
        source_url=f'https://github.com/kcleal/superintervals/blob/{revision}/{path}' if captured_sha==committed_sha else 'https://github.com/kcleal/superintervals'
        result[tool]={'version':version if tool.startswith('SuperIntervals') else None,
                      'source_date':date or None,'sha256':captured_sha,
                      'modified':captured_sha!=committed_sha,'source':source_url}
    for tool,crate in [('Coitrees','coitrees'),('rust-lapper','rust-lapper')]:
        entry=dependencies.get(crate,{})
        if entry.get('version'):
            result[tool]={'version':entry['version'],'source_date':None,
                          'source':f'https://crates.io/crates/{crate}/{entry["version"]}'}
    return result


def generate(directory):
    meta=json.loads((directory/'metadata.json').read_text())
    rows=[]
    if (directory/'summary.csv').exists():
        with (directory/'summary.csv').open(newline='') as f:
            rows=[r for r in csv.DictReader(f) if r['operation'] in ('find','count','stream')]
    # Remove obsolete plots when regenerating an older report.
    for obsolete in directory.glob('rank-*.svg'): obsolete.unlink()
    args=meta['arguments']; fields,model=cpu_fields(meta)
    scope=(f'{100*args["sample_fraction"]:g}% seeded sample' if 'sample_fraction' in args else 'full datasets')
    if args.get('smoke'): scope+=' · up to 1,000 records per input'
    elif args.get('max_records'): scope+=f' · up to {args["max_records"]:,} records per input'
    measured=sorted({r['dataset'] for r in rows},key=lambda label:workload_key((label,'sorted')))
    selected=', '.join(args.get('pairs') or measured or [pair[0] for pair in PAIRS])
    info=[('Status',meta['status']),('CPU',model),('Architecture',meta['machine']),
          ('Platform',meta['platform']),('Scope',scope),('Datasets',selected),
          ('Repetitions / warmups',f'{args["repeats"]} / {args["warmups"]}'),
          ('Seed',args['seed']),('Started UTC',meta['started_utc']),
          ('Source commit',meta.get('git_commit','unavailable')),
          ('Source checkout','modified' if meta.get('git_status') else 'clean'),
          ('Boundary fixture','passed' if meta.get('fixture_passed') else 'not completed')]
    if 'mito_window' in args: info.append(('Mito start window (bp)',args['mito_window'] or 'whole chromosome'))
    if meta.get('archive_sha256'): info.append(('Archive SHA256',meta['archive_sha256']))
    if meta.get('error'): info.append(('Error',meta['error']))
    note='Times are medians in µs, excluding file I/O. Find stores payloads; count returns a total; stream consumes payloads into a checksum.'
    if meta['status']!='complete': note+=' Incomplete run: available results are partial.'
    preferred=['SuperIntervals-rs','SuperIntervals-C++','Coitrees','rust-lapper','ImplicitITree-C++','IntervalTree-C++','NCLS-C']
    workloads=sorted({(r['dataset'],r['order']) for r in rows},key=workload_key)
    lookup={(r['dataset'],r['order'],r['tool'],r['operation']):r for r in rows}
    def table(headers, values):
        return '<div class="table"><table><thead><tr>'+''.join('<th>'+html.escape(str(h))+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(v)+'</td>' for v in row)+'</tr>' for row in values)+'</tbody></table></div>'
    def md_table(headers,values):
        escape=lambda v:str(v).replace('|','\\|').replace('\n',' ')
        return ['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(escape,row))+' |' for row in values]
    md=['# SuperIntervals benchmark report','',f'{model} · {meta["machine"]} · {scope} · {args["repeats"]} repetitions','',note,'']
    body='<h1>SuperIntervals benchmark</h1><p class="run">'+html.escape(f'{model} · {meta["machine"]} · {scope} · {args["repeats"]} repetitions')+'</p><p>'+html.escape(note)+'</p>'
    for op,title in [('find','Find'),('count','Count'),('stream','Stream')]:
        selected_rows=[r for r in rows if r['operation']==op]
        if not selected_rows: continue
        tools=[t for t in preferred if any(r['tool']==t for r in selected_rows)]
        headers=['Dataset','Order',*tools]; values=[]; md_values=[]
        for dataset,order in workloads:
            cells=[]; plain=[]
            for tool in tools:
                r=lookup.get((dataset,order,tool,op))
                if r:
                    detail=f"Q1–Q3: {number(r['query_q1_us'])}–{number(r['query_q3_us'])} µs"
                    cells.append('<span title="'+html.escape(detail,quote=True)+'">'+number(r['query_median_us'])+'</span>')
                    plain.append(f"{number(r['query_median_us'])} [{number(r['query_q1_us'])}, {number(r['query_q3_us'])}]")
                else: cells.append('—'); plain.append('—')
            values.append([html.escape(dataset),html.escape(order),*cells])
            md_values.append([dataset,order,*plain])
        timing=timing_plot(rows,op)
        (directory/f'timing-{op}.svg').write_text(timing)
        body+=f'<section id="{op}"><h2>{title}</h2><div class="chart">'+timing+'</div><details><summary>Absolute query times (µs) — hover for Q1–Q3</summary>'+table(headers,values)+'</details></section>'
        md+=['## '+title,'',f'![Relative query times](timing-{op}.svg)','',*md_table(headers,md_values),'']
    tools=[t for t in preferred if any(r['tool']==t for r in rows)]
    build_values=[]
    for dataset,order in workloads:
        values=[]
        for tool in tools:
            r=next((lookup[dataset,order,tool,op] for op in ('find','count','stream') if (dataset,order,tool,op) in lookup),None)
            values.append(number(r['build_median_us']) if r else '—')
        build_values.append([html.escape(dataset),html.escape(order),*values])
    provenance=meta.get('tool_versions',{})
    version_rows=[]; version_md=[]
    for tool in preferred:
        entry=provenance.get(tool,{})
        if not entry and tool in ('Coitrees','rust-lapper'):
            entry=meta.get('dependencies',{}).get('coitrees' if tool=='Coitrees' else 'rust-lapper',{})
        label=entry.get('version') or (f"Git source date: {entry['source_date']}" if entry.get('source_date') else 'Not recorded')
        if entry.get('modified'): label+=' (modified snapshot)'
        source=entry.get('source','')
        link=('<a href="'+html.escape(source,quote=True)+'">Source</a>') if source.startswith('https://') else '—'
        version_rows.append([html.escape(tool),html.escape(label),link])
        version_md.append([tool,label,source or '—'])
    body+='<h2>Tools used</h2>'+table(['Tool','Version / source date','Source'],version_rows)
    md+=['## Tools used','',*md_table(['Tool','Version / source date','Source'],version_md),'']
    body+='<details class="extra"><summary>Construction times (µs)</summary>'+table(['Dataset','Order',*tools],build_values)+'</details>'
    md+=['## Construction (µs)','',*md_table(['Dataset','Order',*tools],build_values),'']
    inputs=sorted(meta.get('prepared_inputs',{}).items(),key=lambda item:workload_key(tuple(item[0].rsplit('-',1))))
    input_values=[[html.escape(k),f'{v["reference_n"]:,}',f'{v["query_n"]:,}',f'{v["expected_hits"]:,}'] for k,v in inputs]
    if (directory/'runtimes.csv').exists():
        with (directory/'runtimes.csv').open(newline='') as f: runtimes=list(csv.DictReader(f))
        groups={}
        for r in runtimes: groups.setdefault((r['dataset'],r['order'],r['tool']),[]).append(float(r['wall_seconds']))
        runtime_values=[]
        for dataset,order in sorted({(r['dataset'],r['order']) for r in runtimes},key=workload_key):
            cells=[statistics.median(groups[dataset,order,t]) if (dataset,order,t) in groups else None for t in tools]
            total=meta.get('repeat_wall_seconds',{}).get(dataset+'-'+order)
            total=statistics.median(total) if total else sum(c for c in cells if c is not None)
            runtime_values.append([dataset,order,*[f'{c:,.2f}' if c is not None else '—' for c in cells],f'{total:,.2f}'])
        headers=['Dataset','Order',*tools,'Repeat total']
        body+='<details class="extra"><summary>Elapsed time per repeat (seconds, including I/O and validation)</summary>'+table(headers,runtime_values)+'</details>'
        md+=['## Elapsed time per repeat (seconds)','',*md_table(headers,runtime_values),'']
    if meta.get('elapsed_seconds') is not None: info.append(('Total elapsed (s)',number(meta['elapsed_seconds'])))
    if meta.get('build_seconds') is not None: info.append(('Compilation (s)',number(meta['build_seconds'])))
    if meta.get('data_seconds') is not None: info.append(('Dataset extraction (s)',number(meta['data_seconds'])))
    body+='<details class="extra"><summary>Run details and input sizes</summary>'+table(['Field','Value'],[[html.escape(str(k)),html.escape(str(v))] for k,v in info])
    if meta.get('simulation'):
        sim=meta['simulation']
        description=f"Fixed-width simulation: {sim['width_bp']} bp references and queries; uniform starts in {sim['window_bp']:,} bp; {sim['reference_coverage']:g}× reference coverage; seeds {sim['reference_seed']} / {sim['query_seed']}."
        body+='<p>'+html.escape(description)+'</p>'; md+=[description,'']
    body+=table(['Input','References','Queries','Hits'],input_values)
    toolchains={k:v.splitlines()[0] for k,v in meta.get('versions',{}).items()}
    details={'cpu':fields or meta.get('cpu'),'toolchains':toolchains,'flags':{k:meta.get(k) for k in ('CFLAGS','CXXFLAGS','RUSTFLAGS')}}
    body+='<pre>'+html.escape(json.dumps(details,indent=2))+'</pre></details>'
    body+='<p class="files"><a href="raw.csv">Raw CSV</a> · <a href="summary.csv">Summary CSV</a> · <a href="runtimes.csv">Runtimes CSV</a> · <a href="metadata.json">Metadata</a> · <a href="commands.log">Log</a></p>'
    md+=['## Run details','',*md_table(['Field','Value'],info),'',*md_table(['Input','References','Queries','Hits'],input_values),'']
    (directory/'report.md').write_text('\n'.join(md))
    style='body{font:15px/1.5 system-ui,sans-serif;max-width:1250px;margin:32px auto;padding:0 24px;color:#26343a;background:#fafcfc}.run{color:#52656c}h1{font-size:28px;margin-bottom:4px}h2{font-size:21px;margin:30px 0 10px}.chart,.table{overflow:auto}.chart{background:white;border:1px solid #e2e9e9;padding:8px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:7px 10px;border-bottom:1px solid #dde5e5;text-align:right;white-space:nowrap}th{background:#edf2f2}th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left}summary{cursor:pointer;padding:12px 0;color:#245e65}details.extra{margin-top:16px;border-top:1px solid #dde5e5}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}a{color:#245e65}.files{margin-top:24px}'
    (directory/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SuperIntervals benchmark</title><style>'+style+'</style>'+body+'</html>')
    summary=os.environ.get('GITHUB_STEP_SUMMARY')
    if summary:
        with open(summary,'a') as f:
            f.write(f'# Interval benchmark\n\n{model} · {meta["machine"]} · {scope} · {meta["status"]}\n\nHTML report, charts, CSVs and reproducibility files are in the run artifact.\n')
    return len(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=HERE.parents[1])
    p.add_argument('--output',type=Path,required=True,help='New results directory; never overwritten')
    p.add_argument('--report-only',action='store_true',help='Regenerate reports from an existing --output directory')
    p.add_argument('--archive',type=Path,help='Use local release data.tar.gz instead of downloading')
    p.add_argument('--repeats',type=int,default=5); p.add_argument('--seed',type=int,default=1729)
    p.add_argument('--warmups',type=int,default=1); p.add_argument('--smoke',action='store_true',help='Seeded sample capped at 1000 records per input')
    p.add_argument('--sample-fraction',type=float,default=0.25,help='Seeded fraction retained across each whole file (default: 0.25)')
    p.add_argument('--mito-window',type=int,default=1000,help='Keep mitochondrial starts in [0,N); 0 selects the whole chromosome')
    p.add_argument('--pairs',nargs='+',choices=[x[0] for x in PAIRS],help='Select datasets; default includes all seven BED pairs and the fixed-width simulation')
    p.add_argument('--max-records',type=int,help='Cap the seeded sample at N records using reservoir sampling')
    p.add_argument('--timeout',type=int,default=1800,help='Per-tool timeout in seconds')
    a=p.parse_args()
    if a.report_only:
        generate(a.output.expanduser().resolve()); return
    if not 0<a.sample_fraction<=1 or a.mito_window<0: p.error('sample fraction must be in (0,1] and mito window nonnegative')
    if a.repeats<1 or a.warmups<0 or (a.max_records is not None and a.max_records<1): p.error('repeats must be positive and warmups nonnegative')
    repo=a.repo.expanduser().resolve(); out=a.output.expanduser().resolve()
    if not (repo/'test'/'bench.cpp').exists(): p.error('--repo must point to superintervals')
    initial_commit=command(['git','rev-parse','HEAD'],cwd=repo).strip()
    initial_status=command(['git','status','--porcelain'],cwd=repo)
    out.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter()
    meta={'arguments':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
          'platform':platform.platform(),'machine':platform.machine(),'python':sys.version,
          'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'status':'running',
          'CFLAGS':'-O3 -march=native','CXXFLAGS':'-O3 -march=native -std=c++17',
          'RUSTFLAGS':'-C target-cpu=native','coordinate_semantics':'BED half-open, chr1 only',
          'versions':{},'prepared_inputs':{},'source_files':{},
          'runner':{key:os.environ.get(key) for key in ('RUNNER_ARCH','RUNNER_OS','RUNNER_NAME','BENCHMARK_RUNNER_LABEL','GITHUB_REPOSITORY','GITHUB_RUN_ID','GITHUB_RUN_ATTEMPT','GITHUB_SHA')}}
    def save():
        meta['elapsed_seconds']=round(time.perf_counter()-started,3)
        (out/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    save()
    try:
        with open(out/'commands.log','w') as log:
            meta['git_commit']=initial_commit
            meta['git_status']=initial_status
            for name,cmd in {'rustc':['rustc','--version'],'cargo':['cargo','--version'],
                'cc':[os.environ.get('CC','cc'),'--version'],'cxx':[os.environ.get('CXX','c++'),'--version']}.items():
                meta['versions'][name]=command(cmd).strip()
            if platform.system()=='Linux':
                meta['cpu']=Path('/proc/cpuinfo').read_text(); meta['memory']=Path('/proc/meminfo').read_text()
                try: meta['lscpu']=json.loads(command(['lscpu','--json']))
                except (RuntimeError,FileNotFoundError,json.JSONDecodeError) as e: meta['lscpu_error']=str(e)
            elif platform.system()=='Darwin':
                for field,key in [('cpu','machdep.cpu.brand_string'),('memory','hw.memsize')]:
                    try: meta[field]=command(['sysctl','-n',key])
                    except RuntimeError as e: meta[field]='unavailable: '+str(e)
            src=prepare_source(repo,out)

            meta['source_files']={str(f.relative_to(src)):sha(f) for f in sorted(src.rglob('*')) if f.is_file()}
            meta['runner_files']={name:sha(HERE/name) for name in ('run_benchmarks.py','rust_driver.rs','Cargo.lock','make_site.py')}
            stage=time.perf_counter()
            cpp,rust=build(src,out,log)
            resolved=json.loads((out/'cargo-metadata.json').read_text())
            meta['dependencies']={p['name']:{'version':p['version'],'source':p['source']} for p in resolved['packages'] if p['name'] in ('coitrees','rust-lapper')}
            meta['tool_versions']=tool_versions(repo,src,meta['dependencies'],initial_commit)
            meta['build_seconds']=round(time.perf_counter()-stage,3)
            stage=time.perf_counter()
            meta['binary_sha256']={str(f.name):sha(f) for f in (cpp,rust)}
            if a.pairs and set(a.pairs)=={'fixed-width'}:
                dest=None
            else:
                dest,meta['archive_sha256']=data(a.archive,out,log)
            meta['data_seconds']=round(time.perf_counter()-stage,3)
            prepared=out/'prepared'; prepared.mkdir()
            # Boundary, duplicate, singleton and nested intervals, plus a >32-lane block.
            refs=[(0,1),(1,2),(0,10),(0,10),(4,5),(100,200)]+[(i,100-i) for i in range(1,40)]
            queries=[(0,1),(1,2),(2,3),(5,6),(10,11),(99,100),(200,201),(500,501),(0,1000)]
            ref=prepared/'fixture-ref.bed'; query=prepared/'fixture-query.bed'
            write_bed(ref,sorted(refs)); write_bed(query,sorted(queries)); expected=oracle(refs,queries)
            fixture=run_suite(cpp,rust,ref,query,'sorted',a.seed,log,a.timeout)
            if any(row[-1]!=expected for row in fixture): raise RuntimeError(f'Boundary fixture mismatch: expected {expected}, got {fixture}')
            meta['fixture_expected_hits']=expected; meta['fixture_passed']=True
            rows=[]; cache={}
            meta['sampling']={}
            fields=['dataset','order','repeat','tool','operation','reference_n','query_n','build_us','query_us','hits','expected_hits']
            with open(out/'raw.csv','w',newline='') as f, open(out/'runtimes.csv','w',newline='') as runtime_file:
                runtimes=csv.writer(runtime_file)
                runtimes.writerow(['dataset','order','repeat','tool','wall_seconds'])
                writer=csv.DictWriter(f,fieldnames=fields); writer.writeheader()
                for pair_i,(label,rname,qname) in enumerate(PAIRS):
                    if a.pairs and label not in a.pairs: continue
                    limit=1000 if a.smoke else a.max_records
                    if label=='fixed-width':
                        refs,queries,meta['simulation']=simulate_fixed_width(a.seed,limit,a.sample_fraction)
                    else:
                        for name in (rname,qname):
                            if name not in cache:
                                window=a.mito_window if label=='mito' else 0
                                cache[name],meta['sampling'][name]=read_bed(dest/name,limit,a.sample_fraction,a.seed,window)
                        refs=cache[rname]; queries=cache[qname]
                    expected=oracle(refs,queries)
                    print(f'{label}: reference={len(refs):,}, queries={len(queries):,}, expected hits={expected:,}',flush=True)
                    for order in ('sorted','shuffled'):
                        r=sorted(refs); q=sorted(queries)
                        if order=='shuffled':
                            random.Random(a.seed+pair_i*2).shuffle(r); random.Random(a.seed+pair_i*2+1).shuffle(q)
                        ref=prepared/(label+'-'+order+'-ref.bed'); query=prepared/(label+'-'+order+'-query.bed')
                        write_bed(ref,r); write_bed(query,q)
                        meta['prepared_inputs'][label+'-'+order]={'reference_sha256':sha(ref),'query_sha256':sha(query),'reference_n':len(r),'query_n':len(q),'expected_hits':expected}
                        save()
                        for repeat in range(-a.warmups,a.repeats):
                            print(f'{label} {order} '+('warmup' if repeat<0 else f'repeat {repeat+1}/{a.repeats}'),flush=True)
                            def completed(tool,seconds):
                                if repeat>=0:
                                    runtimes.writerow([label,order,repeat+1,tool,f'{seconds:.6f}'])
                                    runtime_file.flush()
                                save()
                            repeat_started=time.perf_counter()
                            result=run_suite(cpp,rust,ref,query,order,a.seed+pair_i*100+repeat,log,a.timeout,completed)
                            seconds=time.perf_counter()-repeat_started
                            print(f'{label} {order}: {seconds:.2f} s per repeat',flush=True)
                            if repeat>=0:
                                meta.setdefault('repeat_wall_seconds',{}).setdefault(label+'-'+order,[]).append(round(seconds,3))
                            for tool,op,bld,elapsed,hits in result:
                                if hits!=expected: raise RuntimeError(f'Hit mismatch {label}/{order}/{tool}/{op}: {hits} != {expected}')
                                if repeat<0: continue
                                row=dict(zip(fields,[label,order,repeat+1,tool,op,len(r),len(q),bld,elapsed,hits,expected]))
                                writer.writerow(row); rows.append(row)
                            f.flush()
                        save()
            groups={}
            for row in rows: groups.setdefault((row['dataset'],row['order'],row['tool'],row['operation']),[]).append(row)
            with open(out/'summary.csv','w',newline='') as f:
                writer=csv.writer(f); writer.writerow(['dataset','order','tool','operation','repeats','build_median_us','query_median_us','query_min_us','query_max_us','query_q1_us','query_q3_us','hits'])
                for key,group in sorted(groups.items(),key=lambda item:(*workload_key(item[0][:2]),item[0][2],item[0][3]!='find')):
                    times=[r['query_us'] for r in group]; qs=statistics.quantiles(times,n=4,method='inclusive') if len(times)>1 else [times[0]]*3
                    writer.writerow([*key,len(times),statistics.median(r['build_us'] for r in group),statistics.median(times),min(times),max(times),qs[0],qs[2],group[0]['hits']])
        meta['status']='complete'; save(); print(f'Results: {out}/summary.csv')
        generate(out)
    except (Exception,KeyboardInterrupt) as e:
        meta['status']='failed'; meta['error']=str(e) or type(e).__name__; save()
        try:
            generate(out)
        except Exception as report_error: print(f'Report generation also failed: {report_error}',file=sys.stderr)
        raise
if __name__=='__main__': main()
