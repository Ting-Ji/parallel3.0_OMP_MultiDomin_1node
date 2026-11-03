"""Short end-to-end cases for the integrated wavefront executable."""
from pathlib import Path
import json, shutil, csv, math, sys
import run_multidomain_validation as v
R=Path(__file__).resolve().parents[1];out=R/'tmp/wavefront_integration/smoke';out.mkdir(parents=True,exist_ok=True)
exe=R/'tmp/wavefront_integration/build/DBEM1_openmp.exe'
if '--baseline' in sys.argv:
 exe=R/'tmp/wavefront_integration/baseline/baseline.exe';out=out.parent/'baseline_smoke';out.mkdir(exist_ok=True)
if '--one-thread' in sys.argv:
 out=out.parent/'serial_smoke';out.mkdir(exist_ok=True)
v.DEFAULT_THREADS=2
if '--one-thread' in sys.argv:v.DEFAULT_THREADS=1
specs=[v.CaseSpec('serial_1domain_n3',1,3,False,single_baseline_segments=1),v.CaseSpec('md_1domain_n3',1,3,True),v.CaseSpec('md_2domain_n5',2,5,True),v.CaseSpec('md_2domain_material_n5',2,5,True,segmented_material=True)]
summary=[]
for spec in specs:
 if '--legacy' not in sys.argv and not spec.multi_domain:continue
 if '--md-only' in sys.argv and not spec.multi_domain:continue
 folder=v.prepare_case(out,R,exe,spec)
 if (R/'lib/x64/pthreadVC2.dll').exists():shutil.copy2(R/'lib/x64/pthreadVC2.dll',folder/'pthreadVC2.dll')
 result=v.run_case(folder,spec,120,require_validation=False)
 assert v.is_success_returncode(result.returncode),(spec.name,result.returncode,result.log_text[-4000:])
 assert result.state_path.exists(),(spec.name,'missing state')
 with result.state_path.open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
 for row in rows:
  for name in v.STATE_FIELDS:assert math.isfinite(float(row[name])),(spec.name,row)
 summary.append(dict(name=spec.name,returncode=result.returncode,state_rows=len(rows),state=str(result.state_path)))
 print(spec.name,'PASS',len(rows),flush=True)
(out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
