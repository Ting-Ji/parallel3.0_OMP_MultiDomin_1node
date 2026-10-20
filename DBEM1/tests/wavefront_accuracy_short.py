"""Short DomainCount=1/2 rod A/B tests; no retained legacy solver or long run."""
from pathlib import Path
import os,sys,time,shutil,csv,json,math,hashlib,re,argparse
import numpy as np
sys.path.insert(0,str(Path.cwd()/"DBEM1/tests"))
import run_step0_policy_tests as h
from wavefront_reference import plane_reference
from wavefront_marching import matrix
p=argparse.ArgumentParser();p.add_argument("work",type=Path);args=p.parse_args()
work=args.work.resolve();runs=work/("short_solve_"+time.strftime("%H%M%S"));runs.mkdir()
v=h.v;v.DEFAULT_E=20000.;v.DEFAULT_V=0.;v.DEFAULT_ROU=8.;v.DEFAULT_THREADS=1;v.DEFAULT_ERROR=1e-10
os.environ.pop("DBEM_AUDIT_KERNEL_ONLY",None);os.environ["DBEM_AUDIT_DUMP"]="1"
summary=[]
for domains in (1,2):
 for beta in (.2,.6):
  v.ROD_LENGTH_PER_DOMAIN=1./domains
  builder=v.RodMeshBuilder()
  for domain in range(domains):
   x0=domain/domains;x1=(domain+1)/domains;a=.1;nx=8//domains
   faces=[
    ("left",(x0,-a,a),(x0,a,a),(x0,a,-a),(x0,-a,-a),1,1),
    ("bottom",(x0,-a,-a),(x1,-a,-a),(x1,-a,a),(x0,-a,a),nx,1),
    ("top",(x0,a,-a),(x0,a,a),(x1,a,a),(x1,a,-a),1,nx),
    ("front",(x0,-a,-a),(x0,a,-a),(x1,a,-a),(x1,-a,-a),1,nx),
    ("back",(x0,-a,a),(x1,-a,a),(x1,a,a),(x0,a,a),nx,1),
    ("right",(x1,-a,-a),(x1,a,-a),(x1,a,a),(x1,-a,a),1,1)]
   for face in faces:builder.add_subdivided_face(*face)
  hmin=min(math.dist(builder.points[e[i]-1],builder.points[e[(i+1)%4]-1]) for e in builder.elements for i in range(4))
  assert abs(hmin-.125)<1e-12,hmin
  dt=beta*hmin/50.;steps=math.ceil(2./(50*dt))
  pair=[]
  for mode in ("old","new"):
   case=runs/(f"d{domains}_beta{beta}_{mode}");(case/"input").mkdir(parents=True);(case/"output").mkdir();(runs/"output").mkdir(exist_ok=True)
   exe=work/mode/"base/x64/Release/base.exe";shutil.copy2(exe,case/"DBEM1.exe")
   for helper in ("gaussposition","gaussweight"):shutil.copy2(h.ROOT/"DBEM1"/helper,case/helper)
   spec=v.CaseSpec("accuracy",domains,steps,True,beta_override=beta)
   v.write_model(case/"input/rod_validation.model",builder)
   size=len(builder.elements)//domains
   with (case/"input/rod_validation.bd").open("w") as fp:
    for i,face in enumerate(builder.face_names):fp.write("123\n" if i<size and face=="left" else "456\n")
    for step in range(steps+1):
     for i,face in enumerate(builder.face_names):
      value=1 if i>=size*(domains-1) and face=="right" and step>0 else 0
      fp.write(f"0 0 {value}\n")
   v.write_datacard(case/"BEM_DATACARD.DAT","rod_validation",spec,builder)
   card=case/"BEM_DATACARD.DAT";lines=card.read_text().splitlines()
   for domain in range(domains):lines[17+domain]=f"{domain+1} 20000 0 8"
   card.write_text("\n".join(lines)+"\n")
   begin=time.perf_counter();result=v.run_case(case,spec,120,require_validation=True)
   assert result.returncode==1,result.log_text[-2000:]
   assert f"MultiDomain config enabled, DomainCount = {domains}" in result.log_text
   assert f"wavefront={int(mode=='new')}" in result.log_text
   assert abs(float(result.metrics["Dt"])-dt)<1e-12*dt,result.metrics
   assert all("C1 50 " in result.metrics[f"DomainMaterial[{domain}]"] for domain in range(domains)),result.metrics
   out=v.result_output_dir(result);rows=v.read_state(result.state_path)
   assert all(math.isfinite(r[k]) for r in rows for k in v.STATE_FIELDS)
   residual=list(csv.DictReader((out/"residual_probe.csv").open()))
   assert len(residual)==steps and all(r["solver"]=="multidomain" for r in residual)
   assert all(math.isfinite(float(r[k])) for r in residual for k in ("relativeResidual","maxAbsResidual","l2Residual","l2Rhs"))
   maxres=max(float(r["relativeResidual"]) for r in residual);assert maxres<=1.01e-10
   assert all(abs(r["tx"]-1.)<1e-12 for r in rows if r["step"]>0 and abs(r["x"]-1.)<1e-12),"normal load must be global +x"
   right=[i for i,face in enumerate(builder.face_names) if face=="right" and i>=len(builder.elements)*(domains-1)//domains]
   tip={};reaction={}
   for r in rows:
    if abs(r["x"])<1e-12:reaction.setdefault(int(r["step"]),[]).append(r["tx"])
   for r in rows:
    if int(r["element"]) in right:tip.setdefault(int(r["step"]),[]).append(r["ux"]*20000.)
   curve=[dict(step=k,tau=k*dt*50,numeric=sum(vals)/len(vals),reference=k*dt*50,reaction=sum(reaction[k])/len(reaction[k]),reaction_reference=(-2. if k*dt*50>=1. else 0.)) for k,vals in sorted(tip.items()) if k*dt*50<=2.+1e-12]
   errors=[r["numeric"]-r["reference"] for r in curve]
   reaction_errors=[r["reaction"]-r["reaction_reference"] for r in curve if abs(r["tau"]-1.)>50*dt+1e-12]
   with (case/"tip.csv").open("w",newline="") as fp:
    writer=csv.DictWriter(fp,fieldnames=curve[0]);writer.writeheader();writer.writerows(curve)
   geom=next(csv.DictReader((out/"audit_geometry.csv").open()))
   R=np.array([float(geom["t"+str(k)]) for k in range(9)]).reshape(3,3).T
   g=matrix(out/"audit_G_0001.bin",True)[:3,:3];A=matrix(out/"audit_A.bin")[:3,:3]
   g0=(-A@R.T-g)/4
   a=.1;pts=np.array([[0,-a,a],[0,a,a],[0,a,-a],[0,-a,-a],[0,0,a],[0,a,0],[0,0,-a],[0,-a,0]],float)
   ref1,_=plane_reference(pts,np.zeros(3),0.,dt,0,1,72)
   ref0=np.diag([dt/(16*50),dt/(16*math.sqrt(1250)),dt/(16*math.sqrt(1250))])
   record=dict(mode=mode,domains=domains,beta=beta,elements=len(builder.elements),steps=steps,dt=dt,
    samples=len(curve),traction_rmse_away_front=math.sqrt(sum(e*e for e in reaction_errors)/len(reaction_errors)),rmse_over_static=math.sqrt(sum(e*e for e in errors)/len(errors)),linf_over_static=max(map(abs,errors)),
    maximum_relative_residual=maxres,self_U1_relative_error=float(np.linalg.norm(g-ref1)/np.linalg.norm(ref1)),
    self_U0_relative_error=float(np.linalg.norm(g0-ref0)/np.linalg.norm(ref0)),seconds=time.perf_counter()-begin,
    case=str(case),output=str(out),input_sha256={str(f.relative_to(case)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [case/"BEM_DATACARD.DAT",*sorted((case/"input").iterdir())]})
   if domains==2:
    assert v.check_interface(result,1e-10,1e-8).passed
    record.pop("self_U1_relative_error");record.pop("self_U0_relative_error")
   pair.append(record);print(json.dumps(record),flush=True)
  assert pair[0]["input_sha256"]==pair[1]["input_sha256"]
  summary.append(dict(domains=domains,beta=beta,old=pair[0],new=pair[1],
   rmse_reduction_percent=100*(1-pair[1]["rmse_over_static"]/pair[0]["rmse_over_static"]),
   traction_rmse_reduction_percent=100*(1-pair[1]["traction_rmse_away_front"]/pair[0]["traction_rmse_away_front"])))
  (work/"short_solve_comparison.json").write_text(json.dumps(summary,indent=2))
print(work/"short_solve_comparison.json",flush=True)
