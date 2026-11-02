"""Build the isolated Windows x64 OpenMP solver and numerical test harness.

Launch with Windows Python from WSL. Uses installed MSVC, no CUDA dependency.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "tmp/gmres_openmp/build")
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    vswhere = Path(os.environ["ProgramFiles(x86)"]) / "Microsoft Visual Studio/Installer/vswhere.exe"
    vs = subprocess.check_output([str(vswhere), "-latest", "-property", "installationPath"], text=True).strip()
    vcvars = Path(vs) / "VC/Auxiliary/Build/vcvars64.bat"
    settings = subprocess.check_output(f'"{vcvars}" >nul && set', shell=True, text=True)
    env = {key.upper(): value for key, value in os.environ.items()}
    for line in settings.splitlines():
        if "=" in line and not line.startswith("="):
            key, value = line.split("=", 1)
            env[key.upper()] = value
    cl = shutil.which("cl.exe", path=env["PATH"])
    common = [cl, "/nologo", "/O2", "/MD", "/EHsc", "/openmp", "/std:c++14", "/fp:precise",
              "/DDBEM_POLICY_TECPLOT_OUTPUT=0", "/DDBEM_GMRES_CPU_VERIFY=1", "/DDBEM_GMRES_CPU_TESTS=1",
              "/I" + str(ROOT / "include"), "/I" + str(ROOT / "DBEM1")]
    ns = {"m": "http://schemas.microsoft.com/developer/msbuild/2003"}
    tree = ET.parse(ROOT / "DBEM1/DBEM1.vcxproj")
    sources = [ROOT / "DBEM1" / x.attrib["Include"] for x in tree.findall(".//m:ClCompile[@Include]", ns)]

    def run(name, command):
        result = subprocess.run(command, cwd=out, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (out / (name + ".log")).write_bytes(result.stdout)
        if result.returncode:
            print(result.stdout.decode("mbcs", errors="replace")[-16000:])
            raise SystemExit(result.returncode)

    link = ["/link", "/LIBPATH:" + str(ROOT / "lib/x64"), "pthreadVC2.lib", "/STACK:52428800"]
    run("solver_build", common + [str(p) for p in sources] + ["/Fe:DBEM1_openmp.exe"] + link)
    run("test_globals_build", common + ["/Dmain=dbem_production_entry", "/c", str(ROOT / "DBEM1/DBEM.cpp"), "/Fo:test_globals.obj"])
    objects = [p.stem + ".obj" for p in sources if p.stem.lower() != "dbem"] + ["test_globals.obj"]
    run("test_build", common + [str(ROOT / "tests/gmres_openmp_test.cpp")] + objects + ["/Fe:gmres_openmp_test.exe"] + link)
    for name in ["gaussposition", "gaussweight"]:
        shutil.copy2(ROOT / "DBEM1" / name, out / name)
    dll = ROOT / "lib/x64/pthreadVC2.dll"
    if dll.exists():
        shutil.copy2(dll, out / dll.name)
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob("*.exe")}
    source_hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted((ROOT / "DBEM1").iterdir()) if p.suffix.lower() in [".h", ".cpp"]}
    (out / "build_manifest.json").write_text(json.dumps({"flags": common[1:], "executables": hashes,
        "source_sha256": source_hashes}, indent=2), encoding="utf-8")
    print(json.dumps(hashes, indent=2))


if __name__ == "__main__":
    main()
