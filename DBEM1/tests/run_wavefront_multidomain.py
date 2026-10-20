from pathlib import Path
import sys,os,json,shutil,time,csv,math,re,argparse,hashlib
sys.path.insert(0,str(Path.cwd()/'DBEM1/tests'))
import run_single_domain_stability as a
v=a.v;root=a.h.ROOT/'validation_runs/wavefront_implementation_20261019'
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['default','new'],default='default');args=p.parse_args()
exe=root/args.mode/'base/x64/Release/base.exe'
work=root/('multi_'+time.strftime('%Y%m%d_%H%M%S')+'_'+args.mode);work.mkdir()
summary=dict(runs=[],mode=args.mode,executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
for label,threads,materials,steps in [
 ('same_t1',1,[(20000,0,8),(20000,0,8)],960),
 ('same_t4',4,[(20000,0,8),(20000,0,8)],960),
 ('same_t8',8,[(20000,0,8),(20000,0,8)],960),
 ('different_t4',4,[(20000,.25,8),(10000,.45,12)],128)]:
 case=work/label;(case/'input').mkdir(parents=True);(case/'output').mkdir();(work/'output').mkdir(exist_ok=True)
 b=v.RodMeshBuilder();nx,ny,nz=8,2,2;p=.1
 for d in range(2):
  x0=d*.5;x1=x0+.5
  faces=[('left',(x0,-p,p),(x0,p,p),(x0,p,-p),(x0,-p,-p),ny,nz),
   ('bottom',(x0,-p,-p),(x1,-p,-p),(x1,-p,p),(x0,-p,p),nx,nz),
   ('top',(x0,p,-p),(x0,p,p),(x1,p,p),(x1,p,-p),nz,nx),
   ('front',(x0,-p,-p),(x0,p,-p),(x1,p,-p),(x1,-p,-p),ny,nx),
   ('back',(x0,-p,p),(x1,-p,p),(x1,p,p),(x0,p,p),nx,ny),
   ('right',(x1,-p,-p),(x1,p,-p),(x1,p,p),(x1,-p,p),ny,nz)]
  for face in faces:b.add_subdivided_face(*face)
 def center(i):return tuple(round(sum(b.points[node-1][k] for node in b.elements[i][:4])/4,12) for k in range(3))
 def normal(i):
  pts=[b.points[j-1] for j in b.elements[i][:3]]
  u=[pts[1][j]-pts[0][j] for j in range(3)];w=[pts[2][j]-pts[0][j] for j in range(3)]
  n=[u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]]
  length=math.sqrt(sum(x*x for x in n));return [x/length for x in n]
 n=len(b.elements);size=n//2
 right=[i for i in range(size) if b.face_names[i]=='right']
 left={center(i):i for i in range(size,n) if b.face_names[i]=='left'}
 pairs=[(1,i+1,2,left[center(i)]+1,-1) for i in right]
 assert len(pairs)==4
 for _,i,_,j,_ in pairs:
  assert center(i-1)==center(j-1) and abs(sum(x*y for x,y in zip(normal(i-1),normal(j-1)))+1)<1e-12
 v.DEFAULT_THREADS=threads;v.DEFAULT_E=materials[0][0];v.DEFAULT_V=materials[0][1];v.DEFAULT_ROU=materials[0][2]
 spec=v.CaseSpec(label,2,steps,True,beta_override=.4)
 v.write_model(case/'input/rod_validation.model',b);v.write_datacard(case/'BEM_DATACARD.DAT','rod_validation',spec,b)
 lines=(case/'BEM_DATACARD.DAT').read_text().splitlines()
 lines=lines[:17]+[str(i+1)+' '+' '.join(map(str,m)) for i,m in enumerate(materials)]+[str(len(pairs))]+[' '.join(map(str,pair)) for pair in pairs]
 (case/'BEM_DATACARD.DAT').write_text('\n'.join(lines)+'\n',encoding='ascii')
 with (case/'input/rod_validation.bd').open('w') as fp:
  for i,face in enumerate(b.face_names):fp.write('123\n' if i<size and face=='left' else '456\n')
  for step in range(steps+1):
   for i,face in enumerate(b.face_names):
    value=1 if i>=size and face=='right' and step>0 else 0
    fp.write(f'0 0 {value}\n')
 for helper in ('gaussposition','gaussweight'):shutil.copy2(a.h.ROOT/'DBEM1'/helper,case/helper)
 shutil.copy2(exe,case/'DBEM1.exe')
 os.environ['DBEM_AUDIT_DUMP']='1'
 start=time.perf_counter();r=v.run_case(case,spec,1800,require_validation=True)
 output=v.result_output_dir(r);rows=v.read_state(r.state_path) if r.state_path.exists() else []
 flags=re.findall(r'MultiDomain GMRES step (\d+): flag=(-?\d+)',r.log_text)
 interfaces=v.check_interface(r,1e-10,1e-8);residual=v.check_residual_probe(r,abs_tol=1e-10,rel_tol=1e-8)
 record=dict(name=label,threads=threads,steps=steps,materials=materials,pairs=pairs,seconds=time.perf_counter()-start,
  output=str(output),case=str(case),returncode=r.returncode,interface_passed=interfaces.passed,
  residual_passed=residual.passed,flags_passed=len(flags)==steps and all(f=='0' for _,f in flags),
  finite_state=bool(rows) and all(math.isfinite(row[f]) for row in rows for f in v.STATE_FIELDS))
 record['passed']=r.returncode==1 and record['interface_passed'] and record['residual_passed'] and record['flags_passed'] and record['finite_state']
 record['max_u_over_static']=max((abs(row[f])*20000 for row in rows for f in ('ux','uy','uz')),default=None)
 summary['runs'].append(record);(work/'summary.json').write_text(json.dumps(summary,indent=2))
 print(json.dumps(record),flush=True)
 if not record['passed']:break
print(work/'summary.json',flush=True)

assert len(summary['runs'])==4 and all(r['passed'] for r in summary['runs']),summary
