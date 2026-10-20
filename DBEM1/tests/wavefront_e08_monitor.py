"""Extend only our isolated E08 test controller, preserving the running solver."""
from pathlib import Path
import os,sys,json,time,csv,math
import psutil
ROOT=Path(__file__).resolve().parents[2]
work=ROOT/'validation_runs/wavefront_budget_20261019/e08_full_20261019_173023'
record=json.loads((work/'summary.json').read_text());case=Path(record['case'])
solver=psutil.Process(record['pid']);parent=solver.parent()
assert Path(solver.exe()).resolve()==(case/'DBEM1.exe').resolve()
assert parent and any(str(x).endswith('wavefront_e08_full.py') for x in parent.cmdline())
started=solver.create_time()
parent.suspend()
record.update(controller_extended=True,controller_pid=parent.pid)
(work/'summary.json').write_text(json.dumps(record,indent=2))
print(json.dumps(dict(attached=solver.pid,suspended_only_test_controller=parent.pid)),flush=True)
try:
 code=solver.wait(timeout=7200)
finally:
 if not solver.is_running():
  try:parent.terminate()
  except psutil.NoSuchProcess:pass
 else:
  parent.resume()
record.update(returncode=code,seconds=time.time()-started,status='complete' if code==1 else 'failed')
states=list(case.rglob('validation_state.csv'));record['state_files']=[str(p) for p in states]
if states:
 rows=0;finite=True;last=-1
 with states[0].open() as f:
  for row in csv.DictReader(f):
   rows+=1;last=int(row['step'])
   for field in ('ux','uy','uz','tx','ty','tz'):finite=finite and math.isfinite(float(row[field]))
 record.update(state_rows=rows,finite_state=finite,last_state_step=last)
record['passed']=code==1 and record.get('finite_state',False) and record.get('state_rows')==672*501 and record.get('last_state_step')==1000
(work/'summary.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record),flush=True)
assert record['passed'],record
