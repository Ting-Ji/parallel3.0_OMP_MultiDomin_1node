"""Create a guarded source-only repair bundle for the remote E08 checkout."""
from pathlib import Path
import hashlib,json,shutil,subprocess,zipfile
ROOT=Path(__file__).resolve().parents[2]
dest=ROOT/'validation_runs/wavefront_budget_20261019/repair_bundle';dest.mkdir(exist_ok=True)
names=['wavefront_quadrature.cpp','wavefront_quadrature.h','solver_policy.h','pre_gmres.cpp']
canonical=lambda data:hashlib.sha256(data.replace(b'\r\n',b'\n')).hexdigest()
manifest={}
for name in names:
 before=(ROOT/'validation_runs/wavefront_budget_20261019/before/base/src'/name).read_bytes()
 after=(ROOT/'DBEM1'/name).read_bytes()
 out=dest/'DBEM1'/name;out.parent.mkdir(exist_ok=True);out.write_bytes(after)
 manifest[name]=dict(before=canonical(before),after=canonical(after))
(dest/'manifest.json').write_text(json.dumps(manifest,indent=2))
patch=subprocess.check_output(['git','diff','--ignore-space-at-eol','--','DBEM1/wavefront_quadrature.cpp','DBEM1/wavefront_quadrature.h','DBEM1/solver_policy.h','DBEM1/pre_gmres.cpp'],cwd=ROOT)
(dest/'fix.patch').write_bytes(patch)
apply=r'''from pathlib import Path
import argparse,hashlib,json,shutil,time,os
p=argparse.ArgumentParser(description="Apply the E08 near-zero quadrature fix with source-version checks.")
p.add_argument("--root",type=Path,required=True);args=p.parse_args()
bundle=Path(__file__).resolve().parent;root=args.root.resolve()
manifest=json.loads((bundle/"manifest.json").read_text())
canonical=lambda data:hashlib.sha256(data.replace(b"\r\n",b"\n")).hexdigest()
pending=[]
for name,hashes in manifest.items():
 target=root/"DBEM1"/name
 actual=canonical(target.read_bytes())
 if actual==hashes["after"]:continue
 if actual!=hashes["before"]:raise SystemExit("Different source version; no files changed: "+str(target))
 pending.append(name)
if not pending:print("Already applied");raise SystemExit(0)
backup=root/("wavefront_budget_backup_"+time.strftime("%Y%m%d_%H%M%S"));backup.mkdir()
for name in pending:shutil.copy2(root/"DBEM1"/name,backup/name)
try:
 for name in pending:
  target=root/"DBEM1"/name;temporary=target.with_name(target.name+".wf_budget.tmp")
  temporary.write_bytes((bundle/"DBEM1"/name).read_bytes());os.replace(temporary,target)
except Exception:
 for name in pending:shutil.copy2(backup/name,root/"DBEM1"/name)
 raise
print("Applied; backup:",backup)
print("Rebuild the whole Release/x64 project into separate bin/obj directories before rerunning E08.")
'''
(dest/'apply_fix.py').write_text(apply)
report=ROOT/"wavefront_budget_fix_20261019.md"
if report.exists():shutil.copy2(report,dest/"README.md")
with zipfile.ZipFile(ROOT/'wavefront_budget_fix_20261019.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in dest.rglob('*'):
  if f.is_file():z.write(f,f.relative_to(dest))
print(ROOT/'wavefront_budget_fix_20261019.zip')
