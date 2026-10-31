"""Build an optional CUDA DLL and isolated solver executable with existing tools.

Run using Windows Python, optionally launched from WSL. No SDK installation or
changes to the production executable/project settings are made.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import shutil

ROOT = Path(__file__).resolve().parents[1]


def environment(cuda):
    vswhere = Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")) / "Microsoft Visual Studio/Installer/vswhere.exe"
    vs = subprocess.check_output([str(vswhere), "-latest", "-property", "installationPath"], text=True).strip()
    vcvars = Path(vs) / "VC/Auxiliary/Build/vcvars64.bat"
    output = subprocess.check_output(f'"{vcvars}" >nul && set', shell=True, text=True)
    env = dict(os.environ)
    for line in output.splitlines():
        if "=" in line and not line.startswith("="):
            key, value = line.split("=", 1)
            env[key] = value
    env["PATH"] = str(cuda / "bin") + os.pathsep + env.get("Path", env.get("PATH", ""))
    return env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cuda", type=Path, default=Path(os.environ.get("CUDA_PATH", r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6")))
    parser.add_argument("--arch", default="sm_61", help="CUDA target architecture (local GTX 1050 Ti: sm_61)")
    parser.add_argument("--output", type=Path, default=ROOT / "gpu/build")
    parser.add_argument("--benchmark", action="store_true", help="also build the original CPU versus GPU GMRES harness")
    parser.add_argument("--solver", action="store_true", help="also compile the current full solver to the output directory")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    env = environment(args.cuda)
    cmd = [str(args.cuda / "bin/nvcc.exe"), "-shared", "--use-local-env", "-O3", "-std=c++14", "-arch="+args.arch, "--cudart=shared",
           "-Xcompiler", "/MD,/EHsc", os.path.relpath(ROOT / "gpu/gmres_cuda.cu", args.output), "-lcublas", "-lcusparse",
           "-o", "dbem_cuda_gmres.dll"]
    result = subprocess.run(cmd, cwd=args.output, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (args.output / "cuda_build.log").write_bytes(result.stdout)
    if result.returncode:
        print(result.stdout.decode("mbcs", errors="replace")); raise SystemExit(result.returncode)
    if args.solver:
        import xml.etree.ElementTree as ET
        ns = {"m": "http://schemas.microsoft.com/developer/msbuild/2003"}
        tree = ET.parse(ROOT / "DBEM1/DBEM1.vcxproj")
        sources = [os.path.relpath(ROOT / "DBEM1" / x.attrib["Include"], args.output) for x in tree.findall(".//m:ClCompile[@Include]", ns)]
        cmd = [shutil.which("cl.exe", path=env["PATH"]), "/nologo", "/O2", "/MD", "/EHsc", "/openmp", "/std:c++14", "/sdl", "/we4996",
               "/I"+os.path.relpath(ROOT / "include", args.output), "/I"+os.path.relpath(ROOT / "DBEM1", args.output)] + sources + [
               "/Fe:DBEM1_cuda_optional.exe", "/link", "/LIBPATH:"+os.path.relpath(ROOT / "lib/x64", args.output), "pthreadVC2.lib"]
        result = subprocess.run(cmd, cwd=args.output, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (args.output / "solver_build.log").write_bytes(result.stdout)
        if result.returncode:
            print(result.stdout.decode("mbcs", errors="replace")[-12000:]); raise SystemExit(result.returncode)
    if args.benchmark:
        cl = shutil.which("cl.exe", path=env["PATH"])
        common = [cl, "/nologo", "/O2", "/MD", "/EHsc", "/openmp", "/std:c++14", "/sdl", "/we4996",
                  "/I"+os.path.relpath(ROOT/"include",args.output), "/I"+os.path.relpath(ROOT/"DBEM1",args.output)]
        command = common+["/Dmain=dbem_production_entry","/c",os.path.relpath(ROOT/"DBEM1/DBEM.cpp",args.output),"/Fo:dbem_test_globals.obj"]
        r=subprocess.run(command,cwd=args.output,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (args.output/"benchmark_globals_build.log").write_bytes(r.stdout)
        if r.returncode: print(r.stdout.decode("mbcs",errors="replace"));raise SystemExit(r.returncode)
        import xml.etree.ElementTree as ET
        ns={"m":"http://schemas.microsoft.com/developer/msbuild/2003"}
        tree=ET.parse(ROOT/"DBEM1/DBEM1.vcxproj")
        objects=[Path(x.attrib["Include"]).stem+".obj" for x in tree.findall(".//m:ClCompile[@Include]",ns) if Path(x.attrib["Include"]).stem!="DBEM"]+["dbem_test_globals.obj"]
        command=common+[os.path.relpath(ROOT/"gpu/benchmark_project.cpp",args.output)]+objects+["/Fe:benchmark_project.exe","/link","/MAP:benchmark_project.map","/LIBPATH:"+os.path.relpath(ROOT/"lib/x64",args.output),"pthreadVC2.lib"]
        r=subprocess.run(command,cwd=args.output,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (args.output/"benchmark_build.log").write_bytes(r.stdout)
        if r.returncode:print(r.stdout.decode("mbcs",errors="replace"));raise SystemExit(r.returncode)
    for name in ["gaussposition","gaussweight"]:shutil.copy2(ROOT/"DBEM1"/name,args.output/name)
    print("Built:", args.output)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
