"""Build/run all local wavefront checks with Windows environment propagation."""
from pathlib import Path
import argparse,os,subprocess,sys
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('--build',action='store_true');args=p.parse_args()
def run(*args,**kwargs):subprocess.run([sys.executable,*map(str,args)],cwd=ROOT,check=True,**kwargs)
if args.build:run(ROOT/'DBEM1/tests/build_wavefront_tests.py','--mode','default')
folder=ROOT/'validation_runs/wavefront_implementation_20261019/default/base'
env=os.environ.copy();env['DBEM_AUDIT_KERNEL_ONLY']='1'
with (folder/'probe.log').open('wb') as log:
 subprocess.run([str(folder/'x64/Release/base.exe')],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
run(ROOT/'DBEM1/tests/wavefront_reference.py','--mode','default')
run(ROOT/'DBEM1/tests/wavefront_perturbation.py')
print('PASS local wavefront tests, geometric partitions, cache policy, and time-step perturbations')
