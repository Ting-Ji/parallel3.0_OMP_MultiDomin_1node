"""Compare singular initialization against the same serial formula, including 64 workers."""
from pathlib import Path
import argparse,sys,subprocess,os,json
sys.path.insert(0,str(Path.cwd()/"DBEM1/tests"))
import build_discretization_audit as b
p=argparse.ArgumentParser();p.add_argument("--label",required=True);a=p.parse_args()
b.WORK=b.ROOT/"validation_runs/e08_nan_diagnosis_20261019"/a.label
if a.label=="before" and b.WORK.exists():raise SystemExit("before is a frozen pre-fix snapshot; use its saved executable instead of overwriting it")
b.KERNEL=r"""
static int RunDiscretizationKernelProbe(){
 const int count=96;
 DSquareElement::StaticInit(.5,0.,100.,8.,.06,2);
 Point* pts=new Point[8*count];Point* nodes=new Point[count];
 DSquareElement* elements=new DSquareElement[count];
 long** ids=new long*[count];
 const double uv[8][2]={{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
 for(int e=0;e<count;++e){
  ids[e]=new long[1];ids[e][0]=e;
  nodes[e].pt[0]=.37*e;nodes[e].pt[1]=-.21*e;nodes[e].pt[2]=.13*e;
  long geo[8];
  for(int i=0;i<8;++i){geo[i]=8*e+i;for(int k=0;k<3;++k)pts[geo[i]].pt[k]=nodes[e].pt[k];
   pts[geo[i]].pt[(e+1)%3]+=(.17+.003*e)*uv[i][0];
   pts[geo[i]].pt[(e+2)%3]+=(.19+.002*e)*uv[i][1];}
  elements[e].GeoInit(pts,geo);
 }
 DSquareElement::StaticGetBeInfo(count,count,ids);
 std::vector<double> reference(count*27);
 for(int e=0;e<count;++e){elements[e].PthPhyInit(nodes,ids[e],456);
  for(int k=0;k<9;++k){reference[27*e+k]=elements[e].m_OnElementU1ij[0][0][k];reference[27*e+9+k]=elements[e].m_OnElementT1ij[0][0][k];reference[27*e+18+k]=elements[e].m_OnElementT2ij[0][0][k];}}
 FILE*f=0;fopen_s(&f,"singular_race.csv","w");fprintf(f,"threads,round,nonfinite,different,max_abs_difference\n");
 int failures=0;
 for(int threads:{1,4,8,64})for(int round=0;round<4;++round){
  #pragma omp parallel for num_threads(threads) schedule(static)
  for(int e=0;e<count;++e)elements[e].PthPhyInit(nodes,ids[e],456);
  int bad=0,different=0;double difference=0;
  for(int e=0;e<count;++e)for(int k=0;k<27;++k){double value=k<9?elements[e].m_OnElementU1ij[0][0][k]:k<18?elements[e].m_OnElementT1ij[0][0][k-9]:elements[e].m_OnElementT2ij[0][0][k-18];
   if(!std::isfinite(value)){++bad;continue;}
   double delta=std::abs(value-reference[27*e+k]);if(memcmp(&value,&reference[27*e+k],sizeof(double))!=0)++different;if(delta>difference)difference=delta;}
  fprintf(f,"%d,%d,%d,%d,%.17g\n",threads,round,bad,different,difference);fflush(f);
  if(bad||different)++failures;
 }
 fclose(f);
 for(int e=0;e<count;++e)delete[] ids[e];delete[] ids;delete[] elements;delete[] nodes;delete[] pts;DSquareElement::StaticRelease(count);
 return failures?2:0;
}
"""
b.build("base")
folder=b.WORK/"base";env=os.environ.copy();env["DBEM_AUDIT_KERNEL_ONLY"]="1"
r=subprocess.run([str(folder/"x64/Release/base.exe")],cwd=folder,env=env,capture_output=True,timeout=240)
(folder/"run.log").write_bytes(r.stdout+r.stderr)
print((folder/"singular_race.csv").read_text(),flush=True)
(folder/"result.json").write_text(json.dumps(dict(label=a.label,returncode=r.returncode)))
if a.label=="before":assert r.returncode==2,"Race fixture did not reproduce"
else:assert r.returncode==0,"Singular initialization differs from its serial reference"
