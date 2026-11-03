"""Re-run the saved 1000/2000-domain cases with integrated wavefront quadrature.

Windows Python; sequential solves, immutable old evidence, per-case checkpoints.
"""
import argparse,json,os,time,subprocess,sys
from pathlib import Path
from datetime import datetime,timezone
from types import SimpleNamespace
import psutil
import validate_gmres_openmp as v

R=v.ROOT
CASES=[('cube_md1000_x10_face16',100,R/'tmp/rhs_parallel/acceptance'),
       ('cube_md2000_x20_face16',100,R/'tmp/rhs_parallel/acceptance'),
       ('cube_md1000_x10_face16',1000,R/'tmp/rhs_parallel/acceptance'),
       ('cube_md2000_x20_face16',800,R/'tmp/rhs_parallel/cube2000_800/gmres_rhs')]
STAGES=['multidomain.assembly.step0','multidomain.assembly.history.integration','multidomain.assembly.total','multidomain.rhs','multidomain.gmres','run.solve','run.total']
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def report(out,state):
 lines=['# 精确波前积分：1000/2000 域旧算例复跑','',f"状态：{state['status']}",f"当前：{state.get('current_case','-')}",'','全部沿用原输入、40线程、GMRES=omp、RHS=parallel，单次运行。旧结果只读。新方法改变积分矩阵，新旧场值差异单独报告，不作为并行等价测试。','', '| 算例 | 步数 | 旧总时长/s | 新总时长/s | 旧/新总耗时 | 位移相对L2 | 面力相对L2 | 新残差最大值 |','|---|---:|---:|---:|---:|---:|---:|---:|']
 for row in state.get('completed',[]):
  a,b=row['baseline'],row['new'];c=row['comparison'];old=a['stage_timings']['run.total'];new=b['stage_timings']['run.total']
  lines.append(f"| {b['case']} | {b['steps']} | {old:.3f} | {new:.3f} | {old/new:.4f} | {c['groups']['displacement']['relative_L2']:.6g} | {c['groups']['traction']['relative_L2']:.6g} | {b['max_relative_residual']:.6g} |")
  lines+=['',f"## {b['case']} / {b['steps']} 步分阶段耗时",'', '| 阶段 | 旧/s | 新/s | 旧/新 |','|---|---:|---:|---:|']
  for stage in STAGES:
   x=a['stage_timings'].get(stage,0);y=b['stage_timings'].get(stage,0)
   lines.append(f'| {stage} | {x:.3f} | {y:.3f} | {x/y if y else 0:.4f} |')
  lines+=['',f"GMRES 内迭代变化：{c['inner_iteration_difference']}；完整状态行数：{c['rows']}。",f"峰值采样 RSS：旧 {a['peak_sampled_rss_GiB']:.3f} GiB，新 {b['peak_sampled_rss_GiB']:.3f} GiB。",'']
 if state.get('error'):lines+=['','停止原因：'+state['error']]
 lines+=['','计时与历史运行跨轮次比较，包含机器状态影响；没有重复取均值。更小求解残差不单独证明物理精度更高。']
 (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=R/'tmp/wavefront_integration/large_comparison');p.add_argument('--exe',type=Path,default=R/'tmp/wavefront_integration/large_build/DBEM1_openmp.exe');p.add_argument('--timeout',type=float,default=43200);p.add_argument('--detach',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 if a.detach:
  if (a.output/'queue.lock').exists():raise RuntimeError('Existing queue lock; inspect before launching')
  with (a.output/'queue.log').open('ab') as log:
   child=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),*[x for x in sys.argv[1:] if x!='--detach']],cwd=R,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.CREATE_NEW_PROCESS_GROUP)
  print(json.dumps({'pid':child.pid,'output':str(a.output)},ensure_ascii=True));return
 lock=a.output/'queue.lock'
 with lock.open('x') as f:json.dump({'pid':os.getpid(),'created':psutil.Process().create_time()},f)
 state={'status':'auditing','started_utc':datetime.now(timezone.utc).isoformat(),'pid':os.getpid(),'exe_sha256':v.sha(a.exe),'completed':[]}
 def save():
  state['updated_utc']=datetime.now(timezone.utc).isoformat();v.save(a.output/'PROGRESS.json',state);report(a.output,state)
 try:
  save();manifest=read(a.exe.parent/'build_manifest.json');assert manifest['executables'][a.exe.name]==state['exe_sha256']
  oldbuild=read(R/'tmp/rhs_parallel/build/build_manifest.json');assert manifest['flags']==oldbuild['flags'],'Build flags differ'
  for c,n,broot in CASES:
   old=read(broot/c/f'omp_{n}_0/result.json');assert old['success'] and old['solved_steps']==n and old['manifest']['threads']==40
   assert Path(old['state_file']).is_file();assert old['exe_sha256']==oldbuild['executables']['DBEM1_openmp.exe']
   for relative,digest in old['manifest']['input_sha256'].items():assert v.sha(broot/c/f'omp_{n}_0/DBEM1'/relative)==digest
  if psutil.disk_usage(str(a.output)).free<100*2**30:raise RuntimeError('Insufficient free storage')
  os.environ.update(DBEM_RHS_CPU_MODE='parallel',DBEM_RHS_CPU_PROFILE='0')
  for c,n,broot in CASES:
   old=read(broot/c/f'omp_{n}_0/result.json');state.update(status='running',current_case=f'{c}/{n}');save()
   while psutil.virtual_memory().available/2**30<old['peak_sampled_rss_GiB']+2:
    state['status']='waiting_for_memory';save();time.sleep(20)
   state['status']='running';save()
   new=v.execute(SimpleNamespace(output=a.output,exe=a.exe,timeout=a.timeout),c,n,'omp',repeat=0,profile=False)
   run=a.output/c/f'omp_{n}_0';log=(run/'run.log').read_bytes().decode('mbcs',errors='replace')
   assert 'RHS CPU mode=parallel threads=40' in log and 'planar wavefront integration=1' in log
   assert old['manifest']['input_sha256']==new['manifest']['input_sha256']
   state['status']='comparing_full_fields';save();check=v.compare(old,new)
   row={'baseline':old,'new':new,'comparison':check}
   v.save(run/'comparison.json',row);state['completed'].append(row);save()
  state['status']='completed';save()
 except BaseException as ex:
  state.update(status='failed',error=str(ex));save();raise
 finally:lock.unlink()
if __name__=='__main__':main()
