"""Grazing-front regression, with independent references and unit/material rescaling."""
from pathlib import Path
import sys,os,subprocess,csv,json,math
import numpy as np
sys.path.insert(0,str(Path.cwd()/'DBEM1/tests'))
import build_discretization_audit as b
from wavefront_reference import plane_reference
b.WORK=b.ROOT/'validation_runs/wavefront_budget_20261019/grazing'
b.KERNEL=(b.ROOT/'DBEM1/tests/wavefront_grazing_probe.cpp').read_text()
b.build('base')
folder=b.WORK/'base';env=os.environ.copy();env['DBEM_AUDIT_KERNEL_ONLY']='1'
p=subprocess.run([str(folder/'x64/Release/base.exe')],cwd=folder,env=env,capture_output=True,timeout=180)
(folder/'probe.log').write_bytes(p.stdout+p.stderr);assert p.returncode==0
rows=list(csv.DictReader((folder/'grazing_cases.csv').open()));checks=[]
uv=np.array([[-1,-1],[1,-1],[1,1],[-1,1],[0,-1],[1,0],[0,1],[-1,0]])
for row in rows:
 geometry=int(row['geometry']);length=float(row['length']);factor=float(row['material']);kind=int(row['kernel'])
 points=np.column_stack([np.full(8,10.),-.25+.25*uv[:,0],.25+.25*uv[:,1]])*length*10
 source=np.array([-3.75,-1 if geometry==0 else -.25,.25])*length*10
 dt=float(row['dt']);lag=66 if kind==2 else 65
 coarse,_=plane_reference(points,source,0.,dt,kind,lag,24)
 reference,scale=plane_reference(points,source,0.,dt,kind,lag,48)
 if kind==0:reference*=10/factor;coarse*=10/factor;scale*=10/factor
 actual=np.array([float(row['b'+str(k)]) for k in range(9)]).reshape(3,3)
 target=1e-6*np.linalg.norm(reference)+1e-10*max(scale,float(row['characteristic']))
 error=float(np.linalg.norm(actual-reference));change=float(np.linalg.norm(reference-coarse))
 ratio=error/max(target,1e-300)
 entry=dict(geometry=geometry,length=length,material=factor,shift=float(row['shift']),kernel=kind,error=error,target=float(target),ratio=ratio,
 reference_change_over_target=change/max(target,1e-300),evaluations=int(row['evaluations']))
 checks.append(entry)
 assert int(row['status'])!=2 and ratio<=1 and entry['reference_change_over_target']<=.01,entry
result=dict(cases=len(checks),maximum_error_over_tolerance=max(r['ratio'] for r in checks),
 maximum_reference_change_over_tolerance=max(r['reference_change_over_target'] for r in checks),maximum_evaluations=max(r['evaluations'] for r in checks),results=checks)
(b.WORK/'independent_checks.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='results'}),flush=True)
