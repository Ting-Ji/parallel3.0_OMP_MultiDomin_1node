from pathlib import Path
import subprocess,sys,json,hashlib,shutil
ROOT=Path(__file__).resolve().parents[2]
base=ROOT/'validation_runs/wavefront_budget_20261019'
bundle=base/'repair_bundle';fixture=base/'apply_fixture';(fixture/'DBEM1').mkdir(parents=True,exist_ok=True)
manifest=json.loads((bundle/'manifest.json').read_text())
for name in manifest:shutil.copy2(base/'before/base/src'/name,fixture/'DBEM1'/name)
for run in range(2):
 p=subprocess.run([sys.executable,str(bundle/'apply_fix.py'),'--root',str(fixture)],capture_output=True)
 assert p.returncode==0,(p.stdout,p.stderr)
canonical=lambda p:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
assert all(canonical(fixture/'DBEM1'/n)==r['after'] for n,r in manifest.items())
bad=fixture/'DBEM1/wavefront_quadrature.cpp';bad.write_bytes(bad.read_bytes()+b'\n// Other user edits\n')
before={n:canonical(fixture/'DBEM1'/n) for n in manifest}
p=subprocess.run([sys.executable,str(bundle/'apply_fix.py'),'--root',str(fixture)],capture_output=True)
assert p.returncode!=0
assert before=={n:canonical(fixture/'DBEM1'/n) for n in manifest}
(base/'bundle_guard_test.json').write_text(json.dumps(dict(initial_apply=True,idempotent=True,conflicting_source_refused_without_changes=True),indent=2))
print('PASS guarded source bundle')
