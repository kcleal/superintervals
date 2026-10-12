#!/usr/bin/env python3
"""Combine architecture reports into one static Pages dashboard. No dependencies."""
import argparse, csv, html, json, os, shutil, urllib.error, urllib.request
from pathlib import Path

FILES=('raw.csv','summary.csv','runtimes.csv','metadata.json','commands.log')

def build_site(inputs, output, history_url=None):
    roots=[inputs] if (inputs/'metadata.json').exists() else sorted(p.parent for p in inputs.glob('*/metadata.json'))
    reports={}
    for root in roots:
        meta=json.loads((root/'metadata.json').read_text())
        if meta.get('status')!='complete': raise ValueError(f'Incomplete benchmark: {root}')
        arch={'aarch64':'arm64','arm64':'arm64','x86_64':'x86_64','AMD64':'x86_64'}.get(meta['machine'])
        if not arch or arch in reports: raise ValueError(f'Unknown or duplicate architecture: {meta["machine"]}')
        reports[arch]=(root,meta)
    if not reports: raise ValueError('No completed architecture reports found')
    output.mkdir(parents=True,exist_ok=True)
    previous=[]
    if (output/'history.json').exists(): previous=json.loads((output/'history.json').read_text()).get('runs',[])
    elif history_url:
        try:
            with urllib.request.urlopen(history_url,timeout=10) as response:
                previous=json.loads(response.read(5_000_000)).get('runs',[])
        except (urllib.error.URLError,TimeoutError,ValueError): pass
    if not isinstance(previous,list): previous=[]
    started=min(meta['started_utc'] for _,meta in reports.values())
    record={'run_id':os.environ.get('GITHUB_RUN_ID'),'started_utc':started,
            'commit':next(iter(reports.values()))[1].get('git_commit'),'architectures':{}}
    buttons=[];frames=[]
    for i,arch in enumerate(a for a in ('x86_64','arm64') if a in reports):
        root,meta=reports[arch];dest=output/arch;dest.mkdir(exist_ok=True)
        for name in FILES:
            if (root/name).exists(): shutil.copyfile(root/name,dest/name)
        doc=(root/'report.html').read_text()
        doc=doc.replace('<meta charset="utf-8">',f'<meta charset="utf-8"><base href="{arch}/"><style>body{{margin:0!important;padding:0!important}}h1{{display:none}}</style>',1)
        buttons.append(f'<button role="tab" id="tab-{arch}" aria-controls="panel-{arch}" aria-selected="{str(i==0).lower()}" onclick="selectArch(\'{arch}\')">{arch}</button>')
        frames.append(f'<iframe id="panel-{arch}" role="tabpanel" aria-labelledby="tab-{arch}" title="{arch} benchmark" {"hidden" if i else ""} srcdoc="{html.escape(doc,quote=True)}" onload="wire(this)"></iframe>')
        record['architectures'][arch]={'cpu':meta.get('cpu'),'platform':meta['platform'],
            'arguments':meta['arguments'],'versions':meta['versions'],
            'inputs':meta['prepared_inputs'], 'sampling':meta.get('sampling'),
            'simulation':meta.get('simulation'),'flags':{k:meta.get(k) for k in ('CFLAGS','CXXFLAGS','RUSTFLAGS')},
            'summary':list(csv.DictReader((root/'summary.csv').open()))}
    key=record['run_id'] or started
    history=[record]+[r for r in previous if isinstance(r,dict) and (r.get('run_id') or r.get('started_utc'))!=key]
    history=history[:30]
    (output/'history.json').write_text(json.dumps({'schema':1,'runs':history},separators=(',',':'))+'\n')
    rows=[];repo=os.environ.get('GITHUB_REPOSITORY')
    for run in history:
        commit=html.escape(str(run.get('commit') or 'local')[:8]);run_id=str(run.get('run_id') or '')
        if repo and run_id.isdigit(): commit=f'<a href="https://github.com/{html.escape(repo,quote=True)}/actions/runs/{run_id}">{commit}</a>'
        architectures=run.get('architectures',{})
        args=next(iter(architectures.values()),{}).get('arguments',{})
        scope=f'{100*args["sample_fraction"]:g}%' if 'sample_fraction' in args else 'full'
        window=args.get('mito_window',0);scope+=f' · mito {window} bp' if window else ' · whole mito'
        rows.append('<tr><td>'+html.escape(str(run.get('started_utc','')))+'<td>'+commit+'<td>'+html.escape(', '.join(architectures))+'<td>'+html.escape(scope)+'</tr>')
    style='body{font:15px/1.5 system-ui;max-width:1300px;margin:28px auto;padding:0 24px;color:#26343a;background:#fafcfc}h1{font-size:28px}button{padding:9px 20px;border:1px solid #cbd9d9;background:white;cursor:pointer}button[aria-selected=true]{background:#245e65;color:white}iframe{display:block;width:100%;border:0;height:1600px;margin-top:16px}iframe[hidden]{display:none}a{color:#245e65}table{border-collapse:collapse;font-size:13px}td,th{padding:8px 16px;text-align:left;border-bottom:1px solid #dde5e5}details{margin-top:24px;overflow:auto}summary{cursor:pointer}'
    script='''function resize(f){try{f.style.height=f.contentDocument.documentElement.scrollHeight+'px'}catch(e){}}
function wire(f){resize(f);try{f.contentDocument.addEventListener('toggle',()=>requestAnimationFrame(()=>resize(f)),true)}catch(e){}}
function selectArch(a){document.querySelectorAll('[role=tab]').forEach(b=>b.setAttribute('aria-selected',b.id==='tab-'+a));document.querySelectorAll('iframe').forEach(f=>{f.hidden=f.id!=='panel-'+a;if(!f.hidden)resize(f)})}
window.addEventListener('resize',()=>document.querySelectorAll('iframe:not([hidden])').forEach(resize));
document.querySelector('[role=tablist]').addEventListener('keydown',e=>{if(!['ArrowLeft','ArrowRight'].includes(e.key))return;const b=[...document.querySelectorAll('[role=tab]')];const i=b.indexOf(document.activeElement);if(i<0)return;const next=b[(i+(e.key==='ArrowRight'?1:b.length-1))%b.length];next.click();next.focus();e.preventDefault()});'''
    doc='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SuperIntervals benchmarks</title><style>'+style+'</style><h1>SuperIntervals benchmarks</h1><div role="tablist" aria-label="CPU architecture">'+''.join(buttons)+'</div>'+''.join(frames)+'<details><summary>Recent runs</summary><table><tr><th>UTC<th>Source<th>Architecture<th>Input scope</tr>'+''.join(rows)+'</table><p><a href="history.json">Timing history JSON</a></p></details><script>'+script+'</script></html>'
    (output/'index.html').write_text(doc)
    print(f'Generated {output}/index.html: {", ".join(reports)}; {len(history)} history records')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--history-url');a=p.parse_args();build_site(a.input,a.output,a.history_url)
