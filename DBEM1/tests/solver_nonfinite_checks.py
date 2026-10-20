"""Numerical failures must stop before history propagation or state certification."""
from pathlib import Path
import sys,subprocess,os
sys.path.insert(0,str(Path.cwd()/"DBEM1/tests"))
import build_discretization_audit as b
b.WORK=b.ROOT/"validation_runs/e08_nan_diagnosis_20261019/guards"
b.KERNEL=r"""
static int RunDiscretizationKernelProbe(){
 CCSRMat matrix;matrix.initial(1,1);matrix.row_ptr[0]=0;matrix.row_ptr[1]=1;matrix.col_index[0]=0;
 for(int k=0;k<9;++k)matrix.value[k]=k%4==0?1.:0.;
 PreConditioner pre;pre.m_LeafCount=1;pre.m_LeafPointCount=new long[1]{1};pre.m_LeafBeginID=new long[1]{0};pre.m_RePID=new long[1]{0};pre.IID=new long[1]{0};pre.m_PreM=new Wmatrix[1];pre.m_PreM[0].initial(3,3);pre.m_PreM[0]=0.;for(int k=1;k<=3;++k)pre.m_PreM[0](k,k)=1.;
 Wvector rhs(3),initial(3),solution(3);rhs=0.;initial=0.;solution=0.;rhs[1]=1.;thread_num=1;
 int flag=-1,iter[2]={};
 gmres(nullptr,matrix,rhs,pre,3,1e-9,1,initial,solution,flag,iter);
 int failures=flag!=0||solution[1]!=1.;int caught=0;
 rhs[1]=std::numeric_limits<double>::quiet_NaN();
 try{gmres(nullptr,matrix,rhs,pre,3,1e-9,1,initial,solution,flag,iter);}catch(const std::exception&e){++caught;fprintf(stderr,"EXPECTED %s\n",e.what());}
 rhs[1]=1.;pre.m_PreM[0](1,1)=std::numeric_limits<double>::quiet_NaN();
 try{gmres(nullptr,matrix,rhs,pre,3,1e-9,1,initial,solution,flag,iter);}catch(const std::exception&e){++caught;fprintf(stderr,"EXPECTED %s\n",e.what());}
 matrix.value[0]=std::numeric_limits<double>::quiet_NaN();
 PreConditioner empty;
 try{GMRESPreConditioner(nullptr,empty,1,1,matrix);}catch(const std::exception&e){++caught;fprintf(stderr,"EXPECTED %s\n",e.what());}
 _mkdir("output");Point points[1];points[0].pt[0]=points[0].pt[1]=points[0].pt[2]=0.;
 DSquareElement element;element.m_nodelist=points;element.m_nodeID[0]=0;
 BoundaryValue boundary[1];boundary[0].init(1);boundary[0].lu=0.;boundary[0].lt=0.;boundary[0].lu[1]=std::numeric_limits<double>::quiet_NaN();
 MultiDomainModel model;long geo_ids[8]={};long* geometry[1]={geo_ids};
 const bool wrote=WriteValidationOutputs(&element,points,geometry,nullptr,boundary,1,0,1.,2,model);
 FILE*f=0;fopen_s(&f,DBEMOutputPath("validation_state.csv").c_str(),"r");if(f){fclose(f);++failures;}
 failures+=caught!=3||wrote;
 fprintf(stderr,"checks=5 caught=%d state_rejected=%d failures=%d\n",caught,!wrote,failures);
 return failures?2:0;
}
"""
b.build("base");folder=b.WORK/"base";env=os.environ.copy();env["DBEM_AUDIT_KERNEL_ONLY"]="1"
r=subprocess.run([str(folder/"x64/Release/base.exe")],cwd=folder,env=env,capture_output=True,timeout=120)
(folder/"run.log").write_bytes(r.stdout+r.stderr);print((r.stdout+r.stderr).decode(errors="replace"));assert r.returncode==0 and b"Cannot write non-finite state" in r.stderr
