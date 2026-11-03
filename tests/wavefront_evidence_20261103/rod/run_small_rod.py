"""Isolated old/new full-rod comparison at beta=.6 (Windows Python)."""
from pathlib import Path
import os,sys,json,re,shutil,subprocess,hashlib,math,csv,time
import xml.etree.ElementTree as ET
R=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(R/'tools'))
import run_multidomain_validation as v
BASE=Path(__file__).resolve().parent
OUT=BASE/('full_rod_beta06_axial' if '--axial' in sys.argv else 'full_rod_beta06');OUT.mkdir(exist_ok=True)
BUILDROOT=BASE/'full_rod_beta06/build' if '--axial' in sys.argv and '--skip-build' in sys.argv else OUT/'build'
def save(name,data):
 p=OUT/name;p.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
state={'status':'building','completed':[],'parameters':{'L':6,'b':.175,'h':.3,'E':20000,'nu':0,'rho':8,'beta':.6,'dt':.0021,'steps':56,'threads':2,'elements':82,'old_subdivisions':3,'gauss_per_subcell_direction':3}}
def update(status):state['status']=status;state['updated_unix']=time.time();save('PROGRESS.json',state)

def build():
 src=OUT/'source';src.mkdir(exist_ok=True)
 for p in (R/'DBEM1').iterdir():
  if p.suffix.lower() in ('.cpp','.h'):shutil.copy2(p,src/p.name)
 p=src/'precompiler.h';s=p.read_text(encoding='utf-8-sig');s,n=re.subn(r'const int NUMPW\s*=\s*6\s*;', 'const int NUMPW = 3;',s);assert n==1;p.write_text(s,encoding='utf-8-sig')
 ns={'m':'http://schemas.microsoft.com/developer/msbuild/2003'}
 sources=[src/x.attrib['Include'] for x in ET.parse(R/'DBEM1/DBEM1.vcxproj').findall('.//m:ClCompile[@Include]',ns)]
 vcvars=Path('C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Auxiliary/Build/vcvars64.bat')
 settings=subprocess.check_output(f'"{vcvars}" >nul && set',shell=True,text=True)
 env={k.upper():val for k,val in os.environ.items()}
 for line in settings.splitlines():
  if '=' in line and not line.startswith('='):k,val=line.split('=',1);env[k.upper()]=val
 cl=shutil.which('cl.exe',path=env['PATH'])
 common=[cl,'/nologo','/O2','/MD','/EHsc','/openmp','/std:c++14','/fp:precise','/DDBEM_POLICY_TECPLOT_OUTPUT=0','/DDBEM_POLICY_VALIDATION_OUTPUT=1','/DDBEM_GMRES_CPU_VERIFY=1','/DDBEM_GMRES_CPU_TESTS=0','/I'+str(R/'include'),'/I'+str(src)]
 for name,flag in [('old3',0),('exact',1)]:
  folder=OUT/'build'/name;folder.mkdir(parents=True,exist_ok=True)
  args=common+[f'/DDBEM_POLICY_WAVEFRONT_INTEGRATION={flag}']+[str(p) for p in sources]+['/Fe:DBEM1.exe','/link','/LIBPATH:'+str(R/'lib/x64'),'pthreadVC2.lib','/STACK:52428800']
  done=subprocess.run(args,cwd=folder,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  (folder/'build.log').write_bytes(done.stdout);done.check_returncode()
  save(f'build_{name}.json',{'flags':common[1:]+[f'/DDBEM_POLICY_WAVEFRONT_INTEGRATION={flag}'],'exe_sha256':sha(folder/'DBEM1.exe'),'source_sha256':{p.name:sha(p) for p in src.iterdir()}})
  print('BUILD',name,'OK',flush=True)

def inputs(name):
 folder=OUT/'runs'/name/'DBEM1';(folder/'input').mkdir(parents=True,exist_ok=True);(folder.parent/'output').mkdir(exist_ok=True)
 v.ROD_HALF_WIDTH=.175/2;v.ROD_LENGTH_PER_DOMAIN=.3;v.DEFAULT_E=20000;v.DEFAULT_V=0.;v.DEFAULT_ROU=8.;v.DEFAULT_THREADS=2;v.DEFAULT_ERROR=1e-10;v.DEFAULT_PRE_LEAF=10
 v.material_for_domain=lambda index,segmented:v.Material(e=20000.,v=0.,rou=8.)
 builder=v.generate_single_boundary_mesh(20);assert len(builder.elements)==82
 spec=v.CaseSpec('rod_front_case',1,56,True,beta_override=.6)
 v.write_model(folder/'input/rod_front_case.model',builder)
 # Boundary files contain local components (tangent1,tangent2,normal).
 # On the right cap the third local direction is global +x.
 with (folder/'input/rod_front_case.bd').open('w') as f:
  for face in builder.face_names:f.write(('123' if face=='left' else '456')+'\n')
  for step in range(57):
   for face in builder.face_names:f.write('0 0 1\n' if face=='right' and step>0 else '0 0 0\n')
 v.write_datacard(folder/'BEM_DATACARD.DAT','rod_front_case',spec,builder)
 for fn in ('gaussposition','gaussweight'):shutil.copy2(R/'DBEM1'/fn,folder/fn)
 shutil.copy2(BUILDROOT/name/'DBEM1.exe',folder/'DBEM1.exe')
 env=os.environ.copy();env.update(DBEM_GMRES_CPU_MODE='omp',DBEM_RHS_CPU_MODE='parallel',OMP_NUM_THREADS='2',OMP_DYNAMIC='FALSE',DBEM_GMRES_CPU_PROFILE='0',DBEM_RHS_CPU_PROFILE='0')
 return folder,spec,env

def analyze(folder,spec):
 output=folder/'output'/f'rod_front_case_{spec.nstep}';path=output/'validation_state.csv'
 with path.open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
 assert len(rows)==82*57,len(rows)
 for row in rows:
  if float(row['x'])>5.999 and int(row['step'])>0:
   assert abs(float(row['tx'])-1)<1e-10 and abs(float(row['ty']))<1e-10 and abs(float(row['tz']))<1e-10,'Applied load is not global axial traction'
 cp=50.;dt=.0021;scale=56*dt/(8*cp);sse=ref2=vector=0.;maximum=lateral=precursor=0;front_sse=front_ref2=0
 probe={}
 for row in rows:
  step=int(row['step']);x=float(row['x']);t=step*dt;tau=t-(6-x)/cp;truth=max(0.,tau)/(8*cp)
  u=[float(row[k]) for k in ['ux','uy','uz']];assert all(math.isfinite(z) for z in u)
  err=u[0]-truth;sse+=err*err;ref2+=truth*truth;vector+=err*err+u[1]**2+u[2]**2;maximum=max(maximum,abs(err));lateral=max(lateral,abs(u[1]),abs(u[2]))
  if tau<0:precursor=max(precursor,abs(u[0]))
  if 0<=tau<=4*dt:front_sse+=err*err;front_ref2+=truth*truth
  if row['surfaceType']=='0' and abs(float(row['y'])+.0875)<1e-8:
   element=int(row['element']);probe.setdefault(element,[]).append({'step':step,'x':x,'t':t,'computed':u[0],'reference':truth})
 checks=list(output.rglob('gmres_cpu_validation.csv'));assert len(checks)==1
 with checks[0].open() as f:res=list(csv.DictReader(f))
 residual_ok=len(res)==56 and all(int(z['flag'])==0 and math.isfinite(float(z['relative_residual'])) and float(z['relative_residual'])<=float(z['tolerance']) for z in res)
 assert residual_ok,'True residual gate failed'
 timings={}
 with (output/'timing.txt').open() as f:
  for row in csv.DictReader((l for l in f if not l.startswith('#')),delimiter='\t'):timings[row['stage']]=float(row['cumulative_seconds'])
 return {'state_file':str(path),'state_rows':len(rows),'axial_relative_L2':math.sqrt(sse/ref2),'vector_relative_L2':math.sqrt(vector/ref2),'max_axial_error_scaled':maximum/scale,'max_lateral_scaled':lateral/scale,'max_precursor_scaled':precursor/scale,'front_4dt_relative_L2':math.sqrt(front_sse/front_ref2),'max_true_residual':max(float(z['relative_residual']) for z in res),'residual_ok':residual_ok,'timings':timings,'probes':probe}

def report():
 lines=['# 小规模长杆旧新积分对照','',f"状态：{state['status']}",'','L=6，截面0.175×0.175，轴向20段，82单元；E=20000，ν=0，ρ=8；beta=0.6，dt=0.0021，56步，2线程。旧法3×3区域、每片3×3高斯；新法精确波前区间。其余源码、输入和求解配置相同。','', '参考：反射前 ux=max(t-(L-x)/50,0)/400，uy=uz=0。时间窗t≤0.1176<L/c=0.12。所有边界中心、所有输出步纳入统计；波前辅助窗预设为到达后4dt。','', '| 方法 | 轴向相对L2 | 向量相对L2 | 最大轴向误差/统一尺度 | 最大横向位移/尺度 | 总耗时/s |','|---|---:|---:|---:|---:|---:|']
 for name in state['completed']:
  r=state['results'][name];lines.append(f"| {name} | {r['axial_relative_L2']:.6g} | {r['vector_relative_L2']:.6g} | {r['max_axial_error_scaled']:.6g} | {r['max_lateral_scaled']:.6g} | {r['timings']['run.total']:.3f} |")
 lines+=['','统一位移尺度为p0*t_end/(rho*cP)，没有用波前刚到达时的近零参考值归一化。原始CSV、每步真实残差和计时完整保留。目标是检验候选，不预设新方法一定优于旧法。','', '载荷审计：.bd采用单元局部坐标，右端写入第三（法向）分量；输出中右端tx=1、ty=tz=0已逐步核验。初次误写切向载荷的记录留在full_rod_beta06目录并明确标为无效轴向对照，本目录只包含正确轴向加载。']
 if state['status']=='completed' and 'exact' in state.get('results',{}):
  old=state['results']['old3']['axial_relative_L2'];new=state['results']['exact']['axial_relative_L2']
  lines+=['','## 当前结论',f'旧法全时程轴向相对L2为{old:.4%}，新法为{new:.4%}。新法误差相对减少{(old-new)/old:.2%}，但新法仍明显失准，不满足“旧法超过10%、新法不超过1%—2%”的预定目标。不能将几何面积误差直接当作整杆响应改善幅度。','', '真实残差均小于1e-10；仍需通过时空加密及初始奇异项核查定位共同误差，尚未据此断言具体误差来源。reference_audit.json中的拟合时间平移与幅值仅用于诊断，没有调整正式解析误差。']
 if state.get('error'):lines+=['','失败原因：'+state['error']]
 (OUT/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')

try:
 update('building')
 if '--skip-build' not in sys.argv:build()
 state['results']={}
 manifests=[]
 for name in ['old3','exact']:
  update('running_'+name);report();folder,spec,env=inputs(name)
  manifest={str(p.relative_to(folder)):sha(p) for p in [folder/'BEM_DATACARD.DAT',*sorted((folder/'input').iterdir())]};manifests.append(manifest);save('input_'+name+'.json',manifest)
  if len(manifests)==2:assert manifests[0]==manifests[1],'Input mismatch'
  if not (folder/'output'/f'rod_front_case_{spec.nstep}'/'validation_state.csv').exists():
   with (folder.parent/'run.log').open('wb') as log:
    child=subprocess.Popen([str(folder/'DBEM1.exe')],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,creationflags=subprocess.CREATE_NO_WINDOW)
    save('running.json',{'variant':name,'pid':child.pid,'start_unix':time.time(),'exe_sha256':sha(folder/'DBEM1.exe')})
    rc=child.wait(timeout=900)
   assert v.is_success_returncode(rc),(name,rc)
  r=analyze(folder,spec);save('result_'+name+'.json',r);state['results'][name]={k:val for k,val in r.items() if k!='probes'};state['completed'].append(name);update('finished_'+name);report()
  print('RESULT',name,json.dumps(state['results'][name],ensure_ascii=True),flush=True)
 update('completed');report()
except BaseException as ex:
 state['error']=str(ex);update('failed');report();raise
