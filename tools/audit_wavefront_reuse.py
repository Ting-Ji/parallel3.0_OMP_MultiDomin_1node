"""Audit provenance and exact-state comparisons of the completed validation."""
from pathlib import Path
import argparse,csv,hashlib,json,math
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();O=a.output
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifests={mode:json.loads((O/mode/'build_manifest.json').read_text(encoding='utf-8-sig')) for mode in ('off','on')}
assert manifests['off']['source_sha256']==manifests['on']['source_sha256']
for name,expected in manifests['on']['source_sha256'].items():assert sha(R/name.replace('\\','/'))==expected,name
flags=lambda m:[f for f in m['flags'] if not f.startswith('/DDBEM_WAVEFRONT_REUSE_PARTITIONS=')]
assert flags(manifests['off'])==flags(manifests['on'])
records=json.loads((O/'end_to_end.json').read_text());assert len(records)==10
for record in records:
 folder=O/'runs'/record['label']/f'{record["mode"]}_{record["repeat"]}'/'DBEM1'
 assert sha(folder/'DBEM1.exe')==record['exe_sha256']==manifests[record['mode']]['executables']['DBEM1_openmp.exe']
 for name,h in record['input_sha256'].items():assert sha(folder/name.replace('\\','/'))==h
 state=folder/'output'/f'rod_front_case_{record["steps"]}'/'validation_state.csv'
 assert sha(state)==record['state_sha256']
 with state.open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
 assert len(rows)==record['rows']==record['elements']*(record['steps']+1)
 assert set(int(r['step']) for r in rows)==set(range(record['steps']+1))
 assert all(math.isfinite(float(r[k])) for r in rows for k in ('ux','uy','uz','tx','ty','tz'))
 checks=list(csv.DictReader(next(state.parent.rglob('gmres_cpu_validation.csv')).open()))
 assert len(checks)==record['steps']
 assert all(int(r['flag'])==0 and math.isfinite(float(r['relative_residual'])) and float(r['relative_residual'])<=float(r['tolerance']) for r in checks)
for label,n in (('small82',6),('medium328',4)):
 selected=[r for r in records if r['label']==label];assert len(selected)==n
 assert len({r['state_sha256'] for r in selected})==1
 assert all(r['input_sha256']==selected[0]['input_sha256'] for r in selected)
smoke=json.loads((O/'smoke.json').read_text());assert smoke['all_states_exactly_equal'] and len(smoke['records'])==12
for mode in ('off','on'):
 for test in ('wavefront_quadrature_test','wavefront_reuse_test','gmres_openmp_test','rhs_parallel_test'):
  assert b'PASS' in (O/(mode+'_'+test+'.log')).read_bytes()
numeric=lambda mode:[s for s in (O/(mode+'_wavefront_reuse_test.log')).read_bytes().splitlines() if not s.startswith(b'PASS')]
assert numeric('off')==numeric('on')
result={'status':'PASS','source_matches_final_builds':True,'only_build_switch_differs':True,'full_runs':10,'smoke_runs':12,'reference_cases':210,
 'kernel_and_adaptivity_bytes_equal':True,'all_full_state_bytes_equal':True,'full_run_residuals_pass':True,'all_source_and_input_and_exe_hashes_pass':True,
 'timing':json.loads((O/'summary.json').read_text()),'reuse_counters':{mode:(O/(mode+'_wavefront_reuse_test.log')).read_text().splitlines()[-1] for mode in ('off','on')}}
(O/'acceptance.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
