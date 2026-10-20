"""Check hard-front time-step perturbations against the full-disk closed U integral."""
from pathlib import Path
import csv,json,math
import numpy as np
ROOT=Path(__file__).resolve().parents[2]/'validation_runs/wavefront_implementation_20261019'
rows=list(csv.DictReader((ROOT/'default/base/wf_blocks.csv').open()))
results=[]
for case in (0,1):
 for kind in (0,1,2):
  group=[r for r in rows if int(r['case'])==case and int(r['kernel'])==kind]
  assert len(group)==5
  base=next(r for r in group if float(r['epsilon'])==0)
  def block(r):return np.array([float(r['b'+str(k)]) for k in range(9)])
  b=block(base);scale=max(np.linalg.norm(b),float(base['scale']),1e-300)
  for r in group:
   change=float(np.linalg.norm(block(r)-b)/scale)
   entry=dict(case=case,kernel=kind,epsilon=float(r['epsilon']),normalized_change=change)
   assert change<=5e-6,entry
   if kind==0:
    dt=float(r['dt']);exact=np.diag([dt/(16*50),dt/(16*math.sqrt(1250)),dt/(16*math.sqrt(1250))]).ravel()
    entry['closed_form_relative_error']=float(np.linalg.norm(block(r)-exact)/np.linalg.norm(exact))
    assert entry['closed_form_relative_error']<=5e-6,entry
   results.append(entry)
(ROOT/'front_perturbation.json').write_text(json.dumps(results,indent=2))
print(json.dumps(dict(cases=len(results),max_normalized_change=max(r['normalized_change'] for r in results),
 max_closed_form_relative_error=max(r.get('closed_form_relative_error',0) for r in results)),indent=2))
