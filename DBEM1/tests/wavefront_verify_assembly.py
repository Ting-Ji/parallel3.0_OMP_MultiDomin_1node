"""Verify deterministic assembly and independent history reconstruction; run in Windows Python."""
from pathlib import Path
import argparse,csv,json,hashlib,math
import numpy as np
from threadpoolctl import threadpool_limits
from wavefront_marching import analyze,matrix
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'validation_runs/wavefront_implementation_20261019'
def latest(pattern):return sorted(BASE.glob(pattern))[-1]
def compare(a,b):
 a,b=Path(a['output']),Path(b['output'])
 names=sorted(p.name for p in a.glob('audit_*.bin'))
 assert names==sorted(p.name for p in b.glob('audit_*.bin'))
 unequal=[n for n in names if a.joinpath(n).read_bytes()!=b.joinpath(n).read_bytes()]
 state='validation_state.csv'
 result=dict(matrix_files=len(names),unequal=unequal,state_byte_identical=a.joinpath(state).read_bytes()==b.joinpath(state).read_bytes())
 assert not unequal and result['state_byte_identical'],result
 return result
def self_error(r):
 out=Path(r['output']);g=matrix(out/'audit_G_0001.bin',True)[:3,:3]
 row=next(csv.DictReader((out/'audit_geometry.csv').open()))
 R=np.array([float(row['t'+str(k)]) for k in range(9)]).reshape(3,3).T
 A=matrix(out/'audit_A.bin')[:3,:3]
 g0=(-A@R.T-g)/4
 n=R[:,2];P=np.outer(n,n);dt=r['dt'];cp=50.;cs=math.sqrt(1250.)
 exact=dt/(2*8)*((np.eye(3)-P)/cs+P/cp)
 # Full-disk box-window U0 and U1 have equal exact integrals.
 return dict(U1_relative_frobenius=float(np.linalg.norm(g-exact)/np.linalg.norm(exact)),
  U1_normal_relative=float((n@g@n)/(n@exact@n)-1),
  U0_normal_relative=float((n@g0@n)/(n@exact@n)-1))
def early(r):
 rows=list(csv.DictReader(Path(r['output']).joinpath('stability_tip.csv').open()))
 e=[float(x['tip_over_unit_static'])-float(x['analytic_nu0']) for x in rows if float(x['tau'])<=2]
 return math.sqrt(sum(x*x for x in e)/len(e))
def main():
 single=json.loads(latest('runs_*_default/summary.json').read_text())
 multi=json.loads(latest('multi_*_default/summary.json').read_text())
 result=dict(single_manifest=str(latest('runs_*_default/summary.json')),multi_manifest=str(latest('multi_*_default/summary.json')),parallel=[],single=[])
 for prefix in ('wf_small','wf_late'):
  rows={r['threads']:r for r in single['runs'] if r['name'].startswith(prefix)}
  for t in (4,8):result['parallel'].append(dict(case=prefix,threads=[1,t],**compare(rows[1],rows[t])))
 rows={r['threads']:r for r in multi['runs'] if r['name'].startswith('same_')}
 for t in (4,8):result['parallel'].append(dict(case='same_two_domain',threads=[1,t],**compare(rows[1],rows[t])))
 with threadpool_limits(limits=1):
  for r in single['runs'][:2]:
   result['single'].append(dict(**analyze(r),self_integrals=self_error(r),early_tip_rmse_tau2=early(r)))
 legacy=list(BASE.glob("runs_*_legacy/summary.json"))
 if legacy:
  result["legacy_manifest"]=str(sorted(legacy)[-1])
  result["legacy"]=[]
  with threadpool_limits(limits=1):
   for r in json.loads(sorted(legacy)[-1].read_text())["runs"][:2]:
    result["legacy"].append(dict(**analyze(r),self_integrals=self_error(r),early_tip_rmse_tau2=early(r),seconds=r["seconds"]))
 target=BASE/'assembly_and_stability.json';target.write_text(json.dumps(result,indent=2))
 print(target,flush=True)
if __name__=='__main__':main()
