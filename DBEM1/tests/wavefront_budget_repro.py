"""Isolated reproduction of the E08 wavefront-at-vertex budget failure."""
from pathlib import Path
import argparse,sys,os,subprocess,json
sys.path.insert(0,str(Path.cwd()/'DBEM1/tests'))
import build_discretization_audit as b
p=argparse.ArgumentParser();p.add_argument('--label',default='after');args=p.parse_args()
if args.label=='before':raise SystemExit('before is a frozen pre-fix snapshot; run its saved executable/project instead of overwriting it')
b.WORK=b.ROOT/'validation_runs/wavefront_budget_20261019'/args.label
b.EXTRA_DEFINES=[]
b.KERNEL=r"""
static const int wf_unbuffered=[](){setvbuf(stdout,NULL,_IONBF,0);return 0;}();
static int RunDiscretizationKernelProbe() {
 const double uv[8][2]={{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
 Point points[8];long ids[8]={0,1,2,3,4,5,6,7};
 for(int i=0;i<8;++i){points[i].pt[0]=10.;points[i].pt[1]=-.25+.25*uv[i][0];points[i].pt[2]=.25+.25*uv[i][1];}
 DSquareElement e;e.DomainID=0;e.LocalEleID=22;e.GeoInit(points,ids);
 DynaMat mat=BuildDynaMat(0.,100.,8.,.06);
 Wavefront::Vec source{{-3.75,-1.,.25}};
 auto r=Wavefront::Integrate(e,Wavefront::Kernel::T1,source,65,mat);
 FILE*f=0;fopen_s(&f,"budget_repro.json","w");
 fprintf(f,"{\"status\":%d,\"error\":%.17g,\"scale\":%.17g,\"characteristicScale\":%.17g,\"absoluteScale\":%.17g,\"tolerance\":%.17g,\"evaluations\":%zu,\"depth\":%u,\"block\":[",(int)r.status,r.error,r.scale,r.characteristicScale,r.absoluteScale,r.tolerance,r.statistics.evaluations,r.statistics.maxDepth);
 for(int i=0;i<9;++i)fprintf(f,"%s%.17g",i?",":"",r.block[i]);fprintf(f,"]}\n");fclose(f);
 fprintf(stderr,"%s\n",r.message.c_str());
 return 0;
}
"""
b.build('base')
folder=b.WORK/'base';env=os.environ.copy();env['DBEM_AUDIT_KERNEL_ONLY']='1'
r=subprocess.run([str(folder/'x64/Release/base.exe')],cwd=folder,env=env,capture_output=True,timeout=120)
(folder/'repro.log').write_bytes(r.stdout+r.stderr)
print((folder/'budget_repro.json').read_text(),flush=True)
assert r.returncode==0
if args.label!='before':
 data=json.loads((folder/'budget_repro.json').read_text());assert data['status']==1 and data['evaluations']<50000,data
