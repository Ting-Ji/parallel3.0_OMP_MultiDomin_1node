"""Build and validate the integrated planar wavefront path (Windows Python)."""
from pathlib import Path
import json, subprocess, sys
R=Path(__file__).resolve().parents[1]
OUT=R/'tmp/wavefront_integration';OUT.mkdir(parents=True,exist_ok=True)
import shutil
shutil.copy2(R/'tests/wavefront_references.txt',OUT/'references.txt')
if '--skip-build' not in sys.argv:
 subprocess.run([sys.executable,str(R/'tools/build_gmres_openmp.py'),'--rhs-tests','--wavefront-tests','--output',str(OUT/'build')],check=True)
for name,args in [('wavefront_quadrature_test',[str(OUT/'references.txt')]),('gmres_openmp_test',[]),('rhs_parallel_test',[])]:
 result=subprocess.run([str(OUT/'build'/f'{name}.exe'),*args],cwd=OUT/'build',stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
 (OUT/f'{name}.log').write_bytes(result.stdout)
 print(result.stdout.decode('utf-8',errors='replace')[-1500:])
 result.check_returncode()
