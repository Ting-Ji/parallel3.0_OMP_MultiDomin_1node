"""Original E08 geometry/material/time step, shortened only to two physical steps."""
from pathlib import Path
import os,json,shutil,subprocess,time,csv,math,hashlib,argparse
root=Path.cwd();base=root/"validation_runs/e08_nan_diagnosis_20261019";parser=argparse.ArgumentParser();parser.add_argument("--variant",choices=("before","after"),default="after");args=parser.parse_args();exe=base/args.variant/"base/x64/Release/base.exe"
name="E08_D01_R001_h050_b060_N1000";source=root/"DBEM1/input/experiment_08_homogeneous_partition"
work=base/("e08_new_"+args.variant+"_"+time.strftime("%Y%m%d_%H%M%S"));work.mkdir();records=[];cache_reference=state_reference=None
for threads in (1,4,8,64):
 case=work/("t"+str(threads));(case/"input").mkdir(parents=True);(case/"output").mkdir();(work/"output").mkdir(exist_ok=True)
 for ext in ("model","bd"):shutil.copy2(source/(name+"."+ext),case/"input"/(name+"."+ext))
 card=(source/(name+".DATACARD")).read_text().splitlines();card[1]=str(threads);assert card[4]=="1000";card[4]="2";card += ["1","1 200 0 8","0"];(case/"BEM_DATACARD.DAT").write_text("\n".join(card)+"\n")
 for helper in ("gaussposition","gaussweight"):shutil.copy2(root/"DBEM1"/helper,case/helper)
 env=os.environ.copy();env.pop("DBEM_AUDIT_KERNEL_ONLY",None);env["DBEM_AUDIT_DUMP"]="1";start=time.perf_counter()
 r=subprocess.run([str(exe)],cwd=case,env=env,capture_output=True,timeout=240);(case/"stdout.log").write_bytes(r.stdout);(case/"stderr.log").write_bytes(r.stderr);assert b"MultiDomain config enabled, DomainCount = 1" in r.stdout;assert r.returncode==1,(threads,r.returncode,r.stderr[-2000:])
 state=next((case/"output").rglob("validation_state.csv"));rows=list(csv.DictReader(state.open()));assert len(rows)==1344 and int(rows[-1]["step"])==2
 assert all(math.isfinite(float(x[k])) for x in rows for k in ("ux","uy","uz","tx","ty","tz"))
 residual=list(csv.DictReader(state.with_name("residual_probe.csv").open()));assert len(residual)==2
 assert all(x["solver"]=="multidomain" for x in residual);maxres=max(float(x["relativeResidual"]) for x in residual);assert math.isfinite(maxres) and maxres<=1e-9
 matrix_files=sorted(state.parent.glob("audit_*.bin"));assert len(matrix_files)==6,len(matrix_files);cache=b"".join(f.read_bytes() for f in matrix_files);states=state.read_bytes()
 if cache_reference is None:cache_reference=cache;state_reference=states
 else:assert cache==cache_reference,(threads,"matrix mismatch");assert states==state_reference,(threads,"state mismatch")
 records.append(dict(threads=threads,seconds=time.perf_counter()-start,returncode=r.returncode,state_rows=len(rows),last_step=2,max_relative_residual=maxres,cache_sha256=hashlib.sha256(cache).hexdigest(),state_sha256=hashlib.sha256(states).hexdigest(),case=str(case)))
 print(json.dumps(records[-1]),flush=True)
(work/"summary.json").write_text(json.dumps(dict(framework="new_multi_domain_count_1",variant=args.variant,all_matrices_and_states_byte_identical=True,runs=records),indent=2));print(work/"summary.json")
