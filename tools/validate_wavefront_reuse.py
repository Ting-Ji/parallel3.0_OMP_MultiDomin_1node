"""Windows native build, exact-result regression and paired end-to-end timing."""
from pathlib import Path
import argparse,csv,hashlib,json,math,os,shutil,statistics,subprocess,sys,time
import run_multidomain_validation as v
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--dependency-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--build',action='store_true');a=p.parse_args();O=a.output;O.mkdir(parents=True,exist_ok=True)
def save(name,data):(O/name).write_text(json.dumps(data,indent=2),encoding='utf-8')
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def call(args,cwd,log):
 with log.open('wb') as f:done=subprocess.run(args,cwd=cwd,stdout=f,stderr=subprocess.STDOUT)
 done.check_returncode()
if a.build:
 for mode,flag in [('off',0),('on',1)]:
  call([sys.executable,str(R/'tools/build_gmres_openmp.py'),'--wavefront-tests','--reuse-tests','--rhs-tests','--verify','--reuse-partitions',str(flag),'--dependency-root',str(a.dependency_root),'--output',str(O/mode)],R,O/(mode+'_build.log'))
ref=R/'tests/wavefront_references.txt'
for mode in ('off','on'):
 for test in ('wavefront_quadrature_test','gmres_openmp_test','rhs_parallel_test','wavefront_reuse_test'):
  call([str(O/mode/(test+'.exe'))]+([str(ref)] if test.startswith('wavefront_') else []),O/mode,O/(mode+'_'+test+'.log'))
# The final PASS summary contains different work counters, so compare numeric lines.
numeric=[]
for mode in ('off','on'):
 numeric.append([line for line in (O/(mode+'_wavefront_reuse_test.log')).read_bytes().splitlines() if not line.startswith(b'PASS')])
assert numeric[0]==numeric[1],'Cached kernel/adaptivity results differ'
results=[]
def run_template(label,template,steps,repeat,mode):
 folder=O/'runs'/label/f'{mode}_{repeat}'/'DBEM1';folder.mkdir(parents=True,exist_ok=True);(folder.parent/'output').mkdir(exist_ok=True)
 exe=O/mode/'DBEM1_openmp.exe'
 expected_ne={'small82':82,'medium328':328}[label]
 if not (folder/'BEM_DATACARD.DAT').exists():
  shutil.copytree(template/'input',folder/'input',dirs_exist_ok=True)
  card=(template/'BEM_DATACARD.DAT').read_text(encoding='utf-8-sig').splitlines();card[4]=str(steps);card[10]='2'
  (folder/'BEM_DATACARD.DAT').write_text('\n'.join(card)+'\n')
  for fn in ('gaussposition','gaussweight'):shutil.copy2(template/fn,folder/fn)
  ne=int((template/'output'/f'rod_front_case_{int((template/"BEM_DATACARD.DAT").read_text().splitlines()[4])}'/'validation_metrics.txt').read_text().split('ElementCount=')[1].splitlines()[0])
  bd=folder/'input/rod_front_case.bd';bd.write_text('\n'.join(bd.read_text().splitlines()[:ne*(steps+2)])+'\n')
 if (folder.parent/'started.json').exists():assert json.loads((folder.parent/'started.json').read_text())['exe_sha256']==digest(exe),'Stored output belongs to another build'
 else:shutil.copy2(exe,folder/'DBEM1.exe')
 state=folder/'output'/f'rod_front_case_{steps}'/'validation_state.csv'
 if not state.exists():
  assert not (folder.parent/'started.json').exists(),'Inspect interrupted run before restarting'
  env=os.environ.copy();env.update(DBEM_GMRES_CPU_MODE='omp',DBEM_RHS_CPU_MODE='parallel',OMP_NUM_THREADS='2',OMP_DYNAMIC='FALSE')
  with (folder.parent/'run.log').open('wb') as log:
   child=subprocess.Popen([str(folder/'DBEM1.exe')],cwd=folder,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
   meta={'label':label,'mode':mode,'repeat':repeat,'pid':child.pid,'start_unix':time.time(),'exe_sha256':digest(exe)}
   (folder.parent/'started.json').write_text(json.dumps(meta));save('live.json',meta);code=child.wait(timeout=1200)
  assert code in (0,1),(label,code)
 out=state.parent;ne=expected_ne
 assert abs(float((folder.parent/'output/dt.txt').read_text())-.0021)<1e-14
 with state.open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
 assert len(rows)==ne*(steps+1) and all(math.isfinite(float(r[k])) for r in rows for k in ('ux','uy','uz','tx','ty','tz'))
 residuals=list(csv.DictReader(next(out.rglob('gmres_cpu_validation.csv')).open()))
 assert len(residuals)==steps and all(int(r['flag'])==0 and float(r['relative_residual'])<=float(r['tolerance']) for r in residuals)
 timing={}
 with (out/'timing.txt').open() as f:
  for row in csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t'):timing[row['stage']]=float(row['cumulative_seconds'])
 record={'label':label,'mode':mode,'repeat':repeat,'elements':ne,'steps':steps,'exe_sha256':digest(exe),'state_sha256':digest(state),'rows':len(rows),'timing':timing,
  'input_sha256':{str(q.relative_to(folder)):digest(q) for q in [folder/'BEM_DATACARD.DAT',*sorted((folder/'input').iterdir())]}}
 results.append(record);save('end_to_end.json',results);print(label,mode,repeat,timing['run.solve'],flush=True)
base=a.dependency_root/'wavefront_topic'
small=base/'area_fraction_assessment/rod_front_candidate/full_rod_beta06_axial/runs/exact/DBEM1'
medium=base/'error_source_20261103/moderate_v2/n328_exact_normal/DBEM1'
for label,template,steps,repeats in [('small82',small,56,3),('medium328',medium,112,2)]:
 for repeat in range(repeats):
  for mode in (('off','on') if repeat%2==0 else ('on','off')):run_template(label,template,steps,repeat,mode)
 records=[r for r in results if r['label']==label]
 assert len({r['state_sha256'] for r in records})==1,'Full physical state differs'
 assert all(r['input_sha256']==records[0]['input_sha256'] for r in records)
summary={}
for label in ('small82','medium328'):
 med={mode:statistics.median(r['timing']['run.solve'] for r in results if r['label']==label and r['mode']==mode) for mode in ('off','on')}
 summary[label]={'median_solve_seconds':med,'speedup':med['off']/med['on'],'identical_state_bytes':True}
save('summary.json',summary);print(json.dumps(summary,indent=2))
