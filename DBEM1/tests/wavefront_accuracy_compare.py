"""Small controlled A/B test of the exact thread-slot integration entry used by multidomain."""
from pathlib import Path
import sys,os,subprocess,json,csv,time,hashlib,math,argparse
import numpy as np
sys.path.insert(0,str(Path.cwd()/"DBEM1/tests"))
import build_discretization_audit as b
from wavefront_reference import plane_reference,bowl_reference,material
parser=argparse.ArgumentParser()
parser.add_argument("--reuse",type=Path,help="Analyze an existing isolated build without rebuilding.")
args=parser.parse_args()
work=args.reuse or (b.ROOT/"validation_runs"/("wavefront_accuracy_"+time.strftime("%Y%m%d_%H%M%S")))
work.mkdir(parents=True,exist_ok=True)
if not args.reuse:
 for mode in ("old","new"):
  b.WORK=work/mode;b.KERNEL=(b.ROOT/"DBEM1/tests/wavefront_accuracy_probe.cpp").read_text()
  b.EXTRA_DEFINES=["DBEM_POLICY_WAVEFRONT_INTEGRATION="+str(int(mode=="new"))]
  b.build("base")
  folder=b.WORK/"base";env=os.environ.copy();env["DBEM_AUDIT_KERNEL_ONLY"]="1"
  result=subprocess.run([str(folder/"x64/Release/base.exe")],cwd=folder,env=env,capture_output=True,timeout=180)
  (folder/"probe.log").write_bytes(result.stdout+result.stderr)
  assert result.returncode==0,(mode,result.returncode,(result.stdout+result.stderr)[-3000:])
  print("PASS actual integration entry: "+mode,flush=True)
geometries={}
with (work/"new/base/accuracy_geometry.csv").open() as f:
 for r in csv.DictReader(f):geometries.setdefault(int(r["case"]),[]).append(r)
assert (work/"old/base/accuracy_geometry.csv").read_bytes()==(work/"new/base/accuracy_geometry.csv").read_bytes()
blocks={}
for mode in ("old","new"):
 with (work/mode/"base/accuracy_blocks.csv").open() as f:
  blocks[mode]={(int(r["case"]),int(r["kernel"])):r for r in csv.DictReader(f)}
checks=[]
for key,row in blocks["new"].items():
 ident,kind=key;nu=float(row["nu"]);dt=float(row["dt"]);lag=int(row["lag"])
 points=np.array([[float(r[k]) for k in ("x","y","z")] for r in geometries[ident]])
 source=np.array([float(geometries[ident][0][k]) for k in ("sx","sy","sz")])
 def reference(q):
  if ident in (0,1):
   cp,cs=material(nu)
   value=np.diag([dt/(16*cp),dt/(16*cs),dt/(16*cs)])
   return value,float(np.linalg.norm(value))
  if ident in (6,7):return bowl_reference(8 if ident==6 else 9,source,nu,dt,kind,lag,q)
  if ident==9:
   val,scale=plane_reference(points*10,source*10,nu,dt,kind,lag,q)
   factor=10 if kind==0 else 1
   return val*factor,scale*factor
  return plane_reference(points,source,nu,dt,kind,lag,q)
 coarse,_=reference(48);ref,scale=reference(72)
 reference_change=float(np.linalg.norm(ref-coarse)/max(np.linalg.norm(ref),scale,1e-300))
 if reference_change>1e-8:
  coarse=ref;ref,scale=reference(112)
  reference_change=float(np.linalg.norm(ref-coarse)/max(np.linalg.norm(ref),scale,1e-300))
 assert reference_change<=1e-8,(key,reference_change)
 norm=float(np.linalg.norm(ref))
 r=dict(case=ident,kernel=("U","T1","T2")[kind],lag=lag,reference_norm=norm,reference_scale=scale,reference_change=reference_change)
 for mode in ("old","new"):
  v=np.array([float(blocks[mode][key]["b"+str(k)]) for k in range(9)]).reshape(3,3)
  assert np.isfinite(v).all()
  err=float(np.linalg.norm(v-ref))
  r[mode+"_absolute_error"]=err
  r[mode+"_relative_error"]=err/max(norm,1e-300)
 r["error_reduction_factor"]=r["old_absolute_error"]/max(r["new_absolute_error"],1e-300)
 r["new_evaluations"]=int(row["evaluations"])
 assert r["new_relative_error"]<=1e-6,r
 checks.append(r);print(json.dumps(r),flush=True)
report=dict(run=str(work),entry="IntDynaUijJudge/IntDynaTijJudge with explicit long thread slot, as multidomain.cpp",
 comparison="identical current source, material, quadrature tables and compiler; only DBEM_POLICY_WAVEFRONT_INTEGRATION switches",
 results=checks,production_files_modified=False)
(work/"local_accuracy.json").write_text(json.dumps(report,indent=2,allow_nan=False))
print("RESULT "+str(work/"local_accuracy.json"),flush=True)
