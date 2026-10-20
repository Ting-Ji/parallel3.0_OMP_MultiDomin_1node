"""Rerun the unmodified E08 D01 input with the repaired local executable."""
from pathlib import Path
import os,sys,subprocess,shutil,time,json,hashlib,csv,math
ROOT=Path(__file__).resolve().parents[2]
name='E08_D01_R001_h050_b060_N1000'
root=ROOT/'validation_runs/wavefront_budget_20261019'
work=root/('e08_full_'+time.strftime('%Y%m%d_%H%M%S'));case=work/'case'
(case/'input').mkdir(parents=True);(case/'output').mkdir();(work/'output').mkdir()
source=ROOT/'DBEM1/input/experiment_08_homogeneous_partition'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
hashes={}
for ext in ('model','bd','DATACARD'):
 original=source/(name+'.'+ext);dest=case/('BEM_DATACARD.DAT' if ext=='DATACARD' else 'input/'+name+'.'+ext)
 shutil.copy2(original,dest);assert sha(original)==sha(dest);hashes[ext]=sha(dest)
exe=root/'after/base/x64/Release/base.exe';shutil.copy2(exe,case/'DBEM1.exe')
for file in ('gaussposition','gaussweight'):shutil.copy2(ROOT/'DBEM1'/file,case/file)
env=os.environ.copy();env.pop('DBEM_AUDIT_KERNEL_ONLY',None)
env.pop('DBEM_AUDIT_DUMP',None)
start=time.perf_counter()
record=dict(case=str(case),executable_sha256=sha(exe),input_sha256=hashes,status='running',steps=1000,threads=64)
with (case/'stdout.log').open('wb') as out,(case/'stderr.log').open('wb') as err:
 p=subprocess.Popen([str(case/'DBEM1.exe')],cwd=case,env=env,stdout=out,stderr=err)
 record['pid']=p.pid;(work/'summary.json').write_text(json.dumps(record,indent=2))
 print(json.dumps(record),flush=True)
 try:code=p.wait(timeout=7200)
 except subprocess.TimeoutExpired:
  p.kill();p.wait();raise
record.update(returncode=code,seconds=time.perf_counter()-start,status='complete' if code==1 else 'failed')
states=list(case.rglob('validation_state.csv'));record['state_files']=[str(p) for p in states]
if states:
 rows=0;finite=True
 with states[0].open() as f:
  for row in csv.DictReader(f):
   rows+=1
   for field in ('ux','uy','uz','tx','ty','tz'):finite=finite and math.isfinite(float(row[field]))
 record.update(state_rows=rows,finite_state=finite)
record['passed']=code==1 and record.get('finite_state',False) and record.get('state_rows')==672*501
(work/'summary.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record),flush=True)
assert record['passed'],record
