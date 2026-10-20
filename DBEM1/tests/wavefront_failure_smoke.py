"""Exercise fail-fast propagation in modern and legacy assembly paths."""
from pathlib import Path
import json,shutil,subprocess,os,time
import run_single_domain_stability as a
root=a.h.ROOT/'validation_runs/wavefront_implementation_20261019'
work=root/('failure_'+time.strftime('%Y%m%d_%H%M%S'));work.mkdir();work.joinpath("output").mkdir()
source=sorted(root.glob('runs_*_default/case_wf_small_t4'))[-1]
results=[]
for old in (False,True):
 case=work/('legacy_entry' if old else 'new_framework');(case/'input').mkdir(parents=True);(case/'output').mkdir()
 for name in ('rod_validation.model','rod_validation.bd'):shutil.copy2(source/'input'/name,case/'input'/name)
 for name in ('gaussposition','gaussweight'):shutil.copy2(source/name,case/name)
 shutil.copy2(root/'budget/base/x64/Release/base.exe',case/'DBEM1.exe')
 lines=(source/'BEM_DATACARD.DAT').read_text().splitlines()
 if old:
  lines=lines[:16]
  # The old singular initialization uses shared work arrays: run its one-thread path.
 (case/'BEM_DATACARD.DAT').write_text('\n'.join(lines)+'\n')
 # Thread count is line 2, read from an existing reproducible card.
 if old:
  text=(case/'BEM_DATACARD.DAT').read_text().splitlines();text[1]='1';(case/'BEM_DATACARD.DAT').write_text('\n'.join(text)+'\n')
 env=os.environ.copy();env.pop('DBEM_AUDIT_KERNEL_ONLY',None)
 p=subprocess.run([str(case/'DBEM1.exe')],cwd=case,env=env,capture_output=True,timeout=90)
 log=(p.stdout+p.stderr).decode(errors='replace');(case/'run.log').write_text(log,encoding="utf-8")
 ok=p.returncode in (2,-1,4294967295) and 'budget' in log.lower() and 'wavefront' in log.lower()
 ok=ok and 'MultiDomain GMRES step' not in log
 results.append(dict(legacy_entry=old,returncode=p.returncode,passed=ok,log_tail=log[-1800:]))
(root/'failure_propagation.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
assert all(r['passed'] for r in results)
