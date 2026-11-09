"""Windows short single/multi-domain, same/mixed material and 1/2-thread checks."""
from pathlib import Path
import argparse,csv,json,math
import run_multidomain_validation as v
p=argparse.ArgumentParser();p.add_argument('--dependency-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();O=a.output
records=[]
specs=[v.CaseSpec('md1',1,3,True),v.CaseSpec('md2_same',2,5,True),v.CaseSpec('md2_mixed',2,5,True,segmented_material=True)]
for spec in specs:
 states=[]
 for threads in (1,2):
  v.DEFAULT_THREADS=threads
  for mode in ('off','on'):
   folder=v.prepare_case(O/'smoke'/f'{mode}_{threads}',a.dependency_root,O/mode/'DBEM1_openmp.exe',spec)
   result=v.run_case(folder,spec,120,require_validation=False)
   assert v.is_success_returncode(result.returncode)
   with result.state_path.open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
   assert rows and all(math.isfinite(float(r[k])) for r in rows for k in v.STATE_FIELDS)
   values=[[float(row[k]) for k in v.STATE_FIELDS] for row in rows]
   states.append(values);records.append({'case':spec.name,'mode':mode,'threads':threads,'rows':len(rows),'returncode':result.returncode})
 assert all(s==states[0] for s in states),'Physical state differs across cache/thread modes'
 print(spec.name,'PASS: cache off/on and 1/2 threads',flush=True)
(O/'smoke.json').write_text(json.dumps({'records':records,'all_states_exactly_equal':True},indent=2))
