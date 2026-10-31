"""Launch the optional Windows CCSR solver with an explicit CPU/CUDA choice."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend",choices=["cpu","cuda"],default="cpu")
    parser.add_argument("--workdir",type=Path,default=ROOT/"DBEM1",help="directory containing the selected BEM_DATACARD.DAT and Gauss tables")
    args=parser.parse_args()
    if os.name!="nt":parser.error("use Windows Python; it may be launched from WSL")
    for name in ["BEM_DATACARD.DAT","gaussposition","gaussweight"]:
        if not (args.workdir/name).is_file():parser.error(f"missing {args.workdir/name}")
    env=dict(os.environ)
    env["DBEM_GMRES_BACKEND"]=args.backend
    if args.backend=="cuda":
        cuda=Path(env.get("CUDA_PATH",r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6"))
        env["PATH"]=str(cuda/"bin")+os.pathsep+env.get("PATH","")
        env["DBEM_CUDA_GMRES_DLL"]=str(ROOT/"gpu/build/dbem_cuda_gmres.dll")
    # The legacy main opens ../output/log.txt before creating per-case output.
    (args.workdir.resolve().parent/"output").mkdir(exist_ok=True)
    process=subprocess.run([str(ROOT/"gpu/build/DBEM1_cuda_optional.exe")],cwd=args.workdir,env=env)
    # Preserve the project's success code (1) instead of masking solver failures.
    return process.returncode


if __name__=="__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    raise SystemExit(main())
