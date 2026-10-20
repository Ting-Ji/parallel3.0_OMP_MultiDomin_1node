from pathlib import Path
import sys,argparse
sys.path.insert(0,str(Path.cwd()/'DBEM1/tests'))
import build_discretization_audit as b
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['default','new','legacy','budget'],default='new');args=p.parse_args()
b.WORK=b.ROOT/'validation_runs/wavefront_implementation_20261019'/args.mode
b.KERNEL=b.KERNEL.replace('RunDiscretizationKernelProbe','RunIndependentPointKernelProbe')+'\n'+(b.ROOT/'DBEM1/tests/wavefront_probe.cpp').read_text()
b.EXTRA_DEFINES=['DBEM_WAVEFRONT_CACHE_SELF_TEST=1']
if args.mode!='default':b.EXTRA_DEFINES+=['DBEM_POLICY_WAVEFRONT_INTEGRATION='+('0' if args.mode=='legacy' else '1')]
if args.mode=='budget':b.EXTRA_DEFINES+=['DBEM_WAVEFRONT_MAX_EVALUATIONS=8']
b.build('base')
