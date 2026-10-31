"""Run CPU/GPU comparisons without changing production inputs or binaries."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import time
import numpy as np
from build_windows import ROOT, environment


def main():
    cuda=Path(os.environ.get("CUDA_PATH",r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6"))
    env=environment(cuda)
    build=ROOT/"gpu/build";output=ROOT/"gpu/validation";output.mkdir(exist_ok=True)
    start=time.perf_counter()
    result=subprocess.run([str(build/"benchmark_project.exe")],cwd=build,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
    (output/"benchmark.log").write_bytes(result.stdout)
    text=result.stdout.decode("mbcs",errors="replace")
    if result.returncode:print(text[-5000:]);raise SystemExit(result.returncode)
    rows=[json.loads(line.split("BENCH,",1)[1]) for line in text.splitlines() if "BENCH," in line]
    assert len(rows)==12
    summary=[]
    for n in sorted({r["n"] for r in rows}):
        current=[r for r in rows if r["n"]==n and r["repeat"]>0]
        cpu=statistics.median(r["cpu"] for r in current);gpu=statistics.median(r["gpu"] for r in current)
        summary.append(dict(n=n,dense=current[0]["dense"],cpu_median=cpu,gpu_median=gpu,speedup=cpu/gpu,
                            setup=current[0]["setup"],max_residual=max(r["residual"] for r in current)))
    (output/"benchmark.json").write_text(json.dumps(dict(rows=rows,summary=summary),indent=2))
    # Generate small physical rods in isolated input directories using the
    # project's existing builder, including a mapped artificial interface case.
    spec=importlib.util.spec_from_file_location("rod_builder",ROOT/"tools/generate_multidomain_rod.py")
    rod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=rod;spec.loader.exec_module(rod)
    cases=[]
    for domains in [1,2]:
        outputs={};timings={};run_metadata={}
        for backend in ["cpu","cuda"]:
            run=output/"physical"/f"domains{domains}"/backend/"DBEM1"
            inputs=run/"input";inputs.mkdir(parents=True,exist_ok=True)
            for candidate in [run.parent/"output",run/"output"]:
                owned_output=candidate.resolve()
                assert owned_output.is_relative_to(output.resolve()) and owned_output.name=="output"
                if owned_output.exists():shutil.rmtree(owned_output)
                owned_output.mkdir()
            # Each invocation owns this one named case; repeated verification
            # replaces only its generated files within gpu/validation.
            cfg=rod.Config
            for key,value in dict(domains=domains,case="gpu_rod",out_dir=str(inputs),card_out="",
                                  nx=3 if domains==1 else 2,ny=1,nz=1,length=2.,width=.2,height=.2,
                                  nstep=24,threads=4,E=20000.,nu=.25,rho=8.,beta=.8,
                                  iterations=120,error=1e-9,premaxleafpointnum=4,
                                  load_amplitude=1.,write_single_domain_case=False).items():setattr(cfg,key,value)
            rod.main()
            shutil.copy2(inputs/"gpu_rod.DATACARD",run/"BEM_DATACARD.DAT")
            for asset in ["gaussposition","gaussweight"]:shutil.copy2(ROOT/"DBEM1"/asset,run/asset)
            exe=build/"DBEM1_cuda_optional.exe"
            runenv=dict(env);runenv["DBEM_GMRES_BACKEND"]=backend
            runenv["DBEM_CUDA_GMRES_DLL"]=str(build/"dbem_cuda_gmres.dll") if backend=="cuda" else str(build/"intentionally-missing.dll")
            begin=time.perf_counter()
            process=subprocess.run([str(exe)],cwd=run,env=runenv,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
            timings[backend]=time.perf_counter()-begin
            log=process.stdout.decode("mbcs",errors="replace")
            (run.parent/"run.log").write_bytes(process.stdout)
            files=list(run.parent.rglob("validation_state.csv"))
            if len(files)!=1 or process.returncode not in [0,1]:
                print(log[-6000:]);raise RuntimeError(f"physical run failed: {domains} {backend} rc={process.returncode}")
            with files[0].open(newline="") as stream:outputs[backend]=list(csv.DictReader(stream))
            assert len(outputs[backend])>0
            if backend=="cuda":
                assert log.count("CUDA GMRES setup:")==1,"resident context was not reused"
                assert log.count("CUDA GMRES: flag=0")==24,log[-5000:]
                inner=[int(x) for x in re.findall(r"CUDA GMRES: flag=0 outer=\d+ inner=(\d+)",log)]
                assert max(inner)>1,"coupled iterations were not exercised"
            gmres=re.search(r"GMRESTime=([0-9.eE+-]+)",log)
            run_metadata[backend]={"gmres_seconds":float(gmres.group(1)) if gmres else None,"returncode":process.returncode}
        assert len(outputs["cpu"])==len(outputs["cuda"])
        errors={}
        for prefix in ["u","t"]:
            keys=[prefix+axis for axis in "xyz"]
            a=np.array([[float(r[k]) for k in keys] for r in outputs["cpu"]])
            b=np.array([[float(r[k]) for k in keys] for r in outputs["cuda"]])
            diff=np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-300)
            assert np.max(np.abs(a))>1e-14
            assert np.isfinite(diff) and diff<1e-6,(domains,prefix,diff)
            errors[prefix+"_reference_max"]=float(np.max(np.abs(a)))
            errors[prefix+"_relative_difference"]=float(diff)
            errors[prefix+"_max_absolute_difference"]=float(np.max(np.abs(a-b)))
        cases.append(dict(domains=domains,steps=24,timing=timings,metadata=run_metadata,errors=errors))
    report=dict(passed=True,benchmark=summary,physical=cases,total_seconds=time.perf_counter()-start)
    (output/"validation.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
