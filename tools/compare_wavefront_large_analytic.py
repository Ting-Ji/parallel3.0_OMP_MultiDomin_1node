"""Compare the three completed old/new runs with the reflected nu=0 exact wave.

Read-only streaming analysis of saved CSV files. No solver is launched.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,time,math,csv
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1]
ROOT=R/'tmp/wavefront_integration/large_comparison'
OUT=ROOT/'analytic_comparison';OUT.mkdir(exist_ok=True)
def save(name,data):(OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def normals(model):
 with model.open() as f:
  npnt=int(f.readline());p=np.loadtxt(f,max_rows=npnt)
  ne=int(f.readline());ids=np.loadtxt(f,dtype=np.int32,max_rows=ne)-1
 e=p[ids[:,1]]-p[ids[:,0]];g=p[ids[:,3]]-p[ids[:,0]];n=np.cross(e,g);n/=np.linalg.norm(n,axis=1)[:,None]
 return n
def reference(L,steps):
 dt=.03;c=5.;pressure=.2;E=200.
 x=np.arange(round(L*8)+1)/8.;t=np.arange(steps+1)*dt
 xx=x[None,:];tt=t[:,None]
 u=np.zeros((steps+1,len(x)));stress=np.zeros_like(u);away=np.ones_like(u,dtype=bool)
 for k in range(int(c*t[-1]/(2*L))+2):
  a=((2*k+1)*L-xx)/c;b=((2*k+1)*L+xx)/c;sign=(-1.)**k
  u+=pressure/(8*c)*sign*(np.maximum(tt-a,0)-np.maximum(tt-b,0))
  ha=np.where(tt>a+1e-12,1.,np.where(tt<a-1e-12,0.,.5));hb=np.where(tt>b+1e-12,1.,np.where(tt<b-1e-12,0.,.5))
  stress+=pressure*sign*(ha+hb)
  away&=(np.abs(tt-a)>2*dt+1e-12)&(np.abs(tt-b)>2*dt+1e-12)
 u[0]=0;stress[0]=0
 return x,t,u,stress,away
def modal_check(L):
 errors=[];c=5.;ratio=.001
 for t in [0.,.3,1.,2.,3.,6.,8.,12.,18.,24.,30.]:
  for fraction in [0.,.125,.25,.5,.75,1.]:
   x=fraction*L;n=np.arange(12000);wave=ratio*sum((-1.)**k*(max(c*t-(2*k+1)*L+x,0)-max(c*t-(2*k+1)*L-x,0)) for k in range(int(c*t/(2*L))+2))
   k=(n+.5)*np.pi/L;modal=ratio*x-2*ratio/L*np.sum((-1.)**n*np.sin(k*x)*np.cos(c*k*t)/(k*k));errors.append(abs(wave-modal))
 return max(errors)
def analyze(result,L,label):
 file=Path(result['state_file']);steps=result['steps'];ne=96000 if L==10 else 192000
 output=file.parent;model=file.parents[2]/'input'/(result['case']+'.model')
 nx=normals(model)[:,0];assert len(nx)==ne
 xs,ts,ref,sigma,away=reference(L,steps);shape=ref.shape;size=ref.size
 sums=np.zeros(size);counts=np.zeros(size);sse=ref2=vsse=0.;maximum=lateral=0.;nrows=0
 traction_sse=traction_ref2=0.;traction_max=0.;stress_kept=0
 load_max=fixed_max=0.;ustat=.001*L
 file_info=file.stat()
 columns=['step','element','x','y','z','ux','uy','uz','tx','ty','tz']
 for idx,df in enumerate(pd.read_csv(file,usecols=columns,chunksize=250000,dtype={'step':np.int16,'element':np.int32})):
  step=df['step'].to_numpy().astype(np.int64);ele=df['element'].to_numpy();xx=df['x'].to_numpy();xi=np.rint((xx+L/2)*8).astype(np.int64)
  assert np.max(np.abs((xx+L/2)*8-xi))<1e-8
  flat=step*shape[1]+xi;truth=ref.ravel()[flat];u=df['ux'].to_numpy();uy=df['uy'].to_numpy();uz=df['uz'].to_numpy()
  assert np.isfinite(df[columns].to_numpy()).all()
  e=u-truth;sse+=float(e@e);ref2+=float(truth@truth);vsse+=float(e@e+uy@uy+uz@uz);maximum=max(maximum,float(np.max(abs(e))));lateral=max(lateral,float(np.max(abs(uy))),float(np.max(abs(uz))))
  sums+=np.bincount(flat,weights=u,minlength=size);counts+=np.bincount(flat,minlength=size)
  predicted=sigma.ravel()[flat]*nx[ele];mask=away.ravel()[flat]
  tx=df['tx'].to_numpy();ty=df['ty'].to_numpy();tz=df['tz'].to_numpy();et=tx-predicted
  traction_sse+=float(et[mask]@et[mask]+ty[mask]@ty[mask]+tz[mask]@tz[mask]);traction_ref2+=float(predicted[mask]@predicted[mask]);stress_kept+=int(mask.sum())
  traction_max=max(traction_max,float(np.max(abs(et[mask]))) if mask.any() else 0)
  loaded=(xi==shape[1]-1)&(step>0)
  if loaded.any():load_max=max(load_max,float(np.max(abs(tx[loaded]-.2))),float(np.max(abs(ty[loaded]))),float(np.max(abs(tz[loaded]))))
  fixed=xi==0
  if fixed.any():fixed_max=max(fixed_max,float(np.max(abs(u[fixed]))),float(np.max(abs(uy[fixed]))),float(np.max(abs(uz[fixed]))))
  nrows+=len(df)
  if idx%40==0:save('CHECKPOINT.json',{'label':label,'rows_processed':nrows,'rows_expected':ne*(steps+1),'updated_unix':time.time()})
 assert nrows==ne*(steps+1);assert file.stat().st_size==file_info.st_size
 means=np.divide(sums,counts,out=np.zeros_like(sums),where=counts>0).reshape(shape);valid=counts.reshape(shape)>0
 end=means[:,-1];endref=ref[:,-1];diff=end-endref
 axis_diff=(means-ref)[valid];axis_ref=ref[valid]
 with (output/'gmres_cpu_validation.csv').open() as f:residuals=list(csv.DictReader(f))
 assert len(residuals)==steps and all(int(z['flag'])==0 and float(z['relative_residual'])<=float(z['tolerance']) for z in residuals)
 time_rows=[{'step':i,'t_ms':float(ts[i]),'end_computed':float(end[i]),'end_reference':float(endref[i]),'end_error_over_ustat':float(diff[i]/ustat)} for i in range(steps+1)]
 with (OUT/f'{label}_end_waveform.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=time_rows[0]);w.writeheader();w.writerows(time_rows)
 with (OUT/f'{label}_axial_means.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['step','t_ms','x_from_fixed','computed_ux','reference_ux'])
  for i in range(steps+1):
   for j in np.flatnonzero(valid[i]):w.writerow([i,ts[i],xs[j],means[i,j],ref[i,j]])
 value={'label':label,'case':result['case'],'steps':steps,'L':L,'t_end_ms':float(ts[-1]),'rows':nrows,'end_relative_L2':float(np.linalg.norm(diff)/np.linalg.norm(endref)),'end_max_error_over_ustat':float(np.max(abs(diff))/ustat),'axial_mean_relative_L2':float(np.linalg.norm(axis_diff)/np.linalg.norm(axis_ref)),'full_field_axial_relative_L2':math.sqrt(sse/ref2),'full_field_vector_relative_L2':math.sqrt(vsse/ref2),'full_field_max_axial_error_over_ustat':maximum/ustat,'max_transverse_displacement_over_ustat':lateral/ustat,'traction_vector_relative_L2_away_jumps':math.sqrt(traction_sse/traction_ref2),'traction_max_axial_error_over_pressure_away_jumps':traction_max/.2,'traction_rows_kept':stress_kept,'load_max_global_error':load_max,'fixed_displacement_max_error':fixed_max,'max_true_residual':max(float(z['relative_residual']) for z in residuals),'state_file':str(file),'state_bytes':file_info.st_size,'state_mtime_ns':file_info.st_mtime_ns,'exe_sha256':result['exe_sha256']}
 assert load_max<1e-10 and fixed_max<1e-10
 save(label+'.json',value);print(json.dumps(value,ensure_ascii=True),flush=True);return value
def report(rows):
 lines=['# 三组已完成算例：旧新积分与解析解对照','', '参考适用性：E=200、nu=0、rho=8、c=5，均匀全局+x阶跃面力0.2，左端固定、侧面自由。nu=0时纯轴向行波严格满足三维方程；即使1000域几何为10×10×10，也不需要细长杆近似。1000域L=10，2000域L=20，dt=0.03ms。','', '参考采用包含全部反射的有限行波和；不将反射前公式延用到长时程。另以12000项模态展开交叉核对。','', '| 算例 | 方法 | 端面均值相对L2 | 端面最大误差/Ustat | 轴向均值相对L2 | 全边界轴向相对L2 |','|---|---|---:|---:|---:|---:|']
 for r in rows:lines.append(f"| {r['case']} / {r['steps']} | {r['label'].split('_')[-1]} | {r['end_relative_L2']:.4%} | {r['end_max_error_over_ustat']:.4%} | {r['axial_mean_relative_L2']:.4%} | {r['full_field_axial_relative_L2']:.4%} |")
 lines+=['', '端面均值为加载端全部1600个单元中心位移的等权均值；轴向均值先按x位置对中心样本等权平均，再对所有x和时刻求相对L2。全边界指标保留界面两侧的边界样本。Ustat=pL/E，1000域0.01、2000域0.02。没有用近零时刻参考值逐点除法。','', '| 算例 | 方法 | 面力向量相对L2（避开跳变±2dt） | 加载面力最大偏差 | 最大真实残差 |','|---|---|---:|---:|---:|']
 for r in rows:lines.append(f"| {r['case']} / {r['steps']} | {r['label'].split('_')[-1]} | {r['traction_vector_relative_L2_away_jumps']:.4%} | {r['load_max_global_error']:.3g} | {r['max_true_residual']:.3g} |")
 lines+=['', '应力/面力在波前有真实跳变，主面力指标排除各位置所有到达/反射事件附近±0.06ms；位移连续，位移指标保留全部时刻。法向从原模型角点计算，参考面力为sigma_xx*n_x，参考横向位移和横向面力为0。','', '本次仅分析已有输出，不重跑求解器、不修改历史文件。不把GMRES残差作为物理精度误差。']
 (OUT/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')

progress=read(ROOT/'PROGRESS.json');assert len(progress['completed'])==3
checks={str(L):modal_check(L) for L in [10.,20.]};save('reference_crosscheck.json',checks)
allrows=[]
for pair in progress['completed']:
 L=10. if 'md1000_' in pair['new']['case'] else 20.
 assert pair['baseline']['manifest']['input_sha256']==pair['new']['manifest']['input_sha256']
 for method,key in [('old','baseline'),('new','new')]:
  result=pair[key];label=f"md{1000 if L==10 else 2000}_n{result['steps']}_{method}"
  cached=OUT/(label+'.json')
  value=read(cached) if cached.exists() else analyze(result,L,label)
  allrows.append(value);save('summary.json',{'status':'running','completed_series':len(allrows),'reference_modal_max_abs_difference':checks,'results':allrows});report(allrows)
save('summary.json',{'status':'completed','completed_series':len(allrows),'reference_modal_max_abs_difference':checks,'results':allrows});report(allrows)
print('ALL_SIX_ANALYTIC_COMPARISONS_COMPLETE',flush=True)
