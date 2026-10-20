from pathlib import Path
import sys,json,time,os,argparse,hashlib
sys.path.insert(0,str(Path.cwd()/'DBEM1/tests'))
import run_single_domain_stability as a
root=a.h.ROOT/'validation_runs/wavefront_implementation_20261019'
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['default','new','legacy'],default='default');p.add_argument('--quick',action='store_true');args=p.parse_args()
work=root/('runs_'+time.strftime('%Y%m%d_%H%M%S')+'_'+args.mode);work.mkdir()
exe=root/args.mode/'base/x64/Release/base.exe'
original=a.v.run_case
def extended(case,spec,timeout,**kw):return original(case,spec,1800,**kw)
a.v.run_case=extended
os.environ['DBEM_AUDIT_DUMP']='1'
summary=dict(runs=[],mode=args.mode,executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
configs=[dict(name='wf_small_t4',beta=.2,threads=4),
 dict(name='wf_late_t4',mesh=(16,2,2),tau=24.,threads=4),
 dict(name='wf_small_t1',beta=.2,threads=1),dict(name='wf_small_t8',beta=.2,threads=8),
 dict(name='wf_late_t1',mesh=(16,2,2),tau=24.,threads=1),
 dict(name='wf_late_t8',mesh=(16,2,2),tau=24.,threads=8),
 dict(name='wf_nu025',nu=.25,tau=4.,threads=4),dict(name='wf_nu045',nu=.45,tau=4.,threads=4)]
if args.quick:configs=configs[:2]
for config in configs:
 try:r,rows=a.run_one(work,exe,config)
 except Exception as e:r=dict(config,error=repr(e),passed_algebra=False)
 summary['runs'].append(r)
 (work/'summary.json').write_text(json.dumps(summary,indent=2))
 print(json.dumps({k:r.get(k) for k in ('name','seconds','passed_algebra','max_u_over_unit_static','max_relative_residual','error')}),flush=True)
 if not r['passed_algebra']:break
print(work/'summary.json',flush=True)

assert len(summary['runs'])==len(configs) and all(r['passed_algebra'] for r in summary['runs']),summary
