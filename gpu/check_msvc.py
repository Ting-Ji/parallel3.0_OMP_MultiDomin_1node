"""Compile the actual GPU bridge users with MSVC SDL and C4996 as error."""
import os,shutil,subprocess,sys
from pathlib import Path
from build_windows import ROOT,environment

def main():
    cuda=Path(os.environ.get("CUDA_PATH",r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6"))
    env=environment(cuda);build=ROOT/"gpu/build";build.mkdir(exist_ok=True)
    cl=shutil.which("cl.exe",path=env["PATH"])
    for source in ["pre_gmres.cpp","multidomain.cpp"]:
        command=[cl,"/nologo","/Zs","/EHsc","/openmp","/std:c++14","/sdl","/we4996",
                 "/I"+os.path.relpath(ROOT/"include",build),"/I"+os.path.relpath(ROOT/"DBEM1",build),
                 os.path.relpath(ROOT/"DBEM1"/source,build)]
        r=subprocess.run(command,cwd=build,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (build/(source+".strict.log")).write_bytes(r.stdout)
        if r.returncode:print(r.stdout.decode("mbcs",errors="replace")[-7000:]);raise SystemExit(r.returncode)
        print("PASS: /sdl /we4996, without _CRT_SECURE_NO_WARNINGS:",source)

if __name__=="__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
