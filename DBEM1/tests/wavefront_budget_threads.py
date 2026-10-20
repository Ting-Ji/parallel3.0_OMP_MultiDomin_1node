from pathlib import Path
import os,json,time,sys
sys.path.insert(0,str(Path.cwd()/'DBEM1/tests'))
import run_single_domain_stability as a
root=a.h.ROOT/'validation_runs/wavefront_budget_20261019'
work=root/('threads_'+time.strftime('%Y%m%d_%H%M%S'));work.mkdir()
exe=a.h.ROOT/'validation_runs/wavefront_implementation_20261019/default/base/x64/Release/base.exe'
os.environ['DBEM_AUDIT_DUMP']='1'
records=[]
for t in (1,4,8):
 r,rows=a.run_one(work,exe,dict(name='short_t'+str(t),threads=t,beta=.2,tau=.2))
 assert r['passed_algebra'],r
 records.append(r)
base=Path(records[0]['output']);names=sorted(p.name for p in base.glob('audit_*.bin'))
comparisons=[]
for r in records[1:]:
 other=Path(r['output'])
 bad=[n for n in names if base.joinpath(n).read_bytes()!=other.joinpath(n).read_bytes()]
 same=base.joinpath('validation_state.csv').read_bytes()==other.joinpath('validation_state.csv').read_bytes()
 assert not bad and same
 comparisons.append(dict(threads=[1,r['threads']],matrix_files=len(names),unequal=bad,state_byte_identical=same))
result=dict(runs=records,comparisons=comparisons);(root/'threads_after.json').write_text(json.dumps(result,indent=2));print(json.dumps(comparisons))
