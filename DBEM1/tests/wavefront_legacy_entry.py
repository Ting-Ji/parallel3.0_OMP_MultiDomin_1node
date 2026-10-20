"""Compare the retained single-domain entry with a one-domain modern input."""
from pathlib import Path
import json,csv,shutil,time
import run_single_domain_stability as a
root=a.h.ROOT/'validation_runs/wavefront_implementation_20261019'
work=root/('entry_'+time.strftime('%Y%m%d_%H%M%S'));work.mkdir();(work/'output').mkdir()
source=sorted(root.glob('runs_*_default/case_wf_small_t4'))[-1];records=[];states=[]
for old in (False,True):
 case=work/('legacy' if old else 'modern');(case/'input').mkdir(parents=True);(case/'output').mkdir()
 for name in ('rod_validation.model','rod_validation.bd'):shutil.copy2(source/'input'/name,case/'input'/name)
 for name in ('gaussposition','gaussweight'):shutil.copy2(source/name,case/name)
 shutil.copy2(root/'default/base/x64/Release/base.exe',case/'DBEM1.exe')
 lines=(source/'BEM_DATACARD.DAT').read_text().splitlines();lines[1]='1';lines[4]='8'
 if old:lines=lines[:16]
 (case/'BEM_DATACARD.DAT').write_text('\n'.join(lines)+'\n')
 spec=a.v.CaseSpec('entry',1,8,not old,beta_override=.2)
 result=a.v.run_case(case,spec,120,require_validation=True)
 assert result.returncode==1,result.log_text[-2000:]
 state=a.v.read_state(result.state_path);states.append(state)
 records.append(dict(legacy=old,returncode=result.returncode,state_rows=len(state),case=str(case)))
comparison=a.compare(states[0],states[1]);assert comparison['passed'],comparison
(root/'entry_compatibility.json').write_text(json.dumps(dict(runs=records,comparison=comparison),indent=2))
print(json.dumps(comparison))
