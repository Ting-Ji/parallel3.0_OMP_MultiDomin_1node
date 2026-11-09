#include "dsquareelement.h"
#include "DBEM.h"
#include "solver_policy.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <thread>

static void require(bool v,const char*s){if(!v)throw std::runtime_error(s);}
static void emit(int row,const char*tag,const Wavefront::Result&r){
 std::cout<<row<<' '<<tag<<' '<<(int)r.status;
 for(double x:r.block)std::cout<<' '<<std::hexfloat<<x;
 std::cout<<' '<<r.error<<' '<<r.scale<<' '<<r.characteristicScale<<' '<<r.absoluteScale<<' '<<r.tolerance
 <<' '<<r.statistics.evaluations<<' '<<r.statistics.leaves<<' '<<r.statistics.maxDepth<<'\n';
}
int main(int argc,char**argv)try{
 require(argc==2,"reference input required");std::ifstream in(argv[1]);int count;in>>count;
 std::size_t builds=0,hits=0,boundsBuilds=0,boundsHits=0;
 for(int row=0;row<count;row++){
  std::string name;in>>name;Point p[8];long ids[8]={0,1,2,3,4,5,6,7};
  for(int i=0;i<4;i++){p[i].pt[0]=0;in>>p[i].pt[1]>>p[i].pt[2];}
  for(int i=4;i<8;i++)for(int j=0;j<3;j++)p[i].pt[j]=(p[i-4].pt[j]+p[(i-3)%4].pt[j])/2;
  double gap,nu,dt,natural;int lag,kind;in>>gap>>nu>>dt>>lag>>kind>>natural;
  for(int j=0;j<9;j++){double ignored;in>>ignored;}require(bool(in),"incomplete references");
  DSquareElement e;e.GeoInit(p,ids);auto mat=BuildDynaMat(nu,20000/(2*(1+nu)),8,dt);
  Wavefront::Vec src{{-gap,0,0}};auto k=(Wavefront::Kernel)kind;
  auto run=[&](const char*tag,Wavefront::Kernel which,long n,const DynaMat&m,const Wavefront::Vec&s,const Wavefront::Options&o=Wavefront::Options()){
   auto r=Wavefront::Integrate(e,which,s,n,m,o);builds+=r.statistics.partitionBuilds;hits+=r.statistics.partitionCacheHits;boundsBuilds+=r.statistics.boundsBuilds;boundsHits+=r.statistics.boundsCacheHits;emit(row,tag,r);return r;
  };
  auto first=run("first",k,lag,mat,src);auto second=run("repeat",k,lag,mat,src);
  require(first.status==second.status&&first.block==second.block&&first.error==second.error&&first.scale==second.scale,"repeat changed integral");
  if(first.statistics.partitionBuilds){
#if DBEM_WAVEFRONT_REUSE_PARTITIONS
   require(second.statistics.partitionCacheHits==1,"same key not reused");
#else
   require(second.statistics.partitionBuilds==1,"disabled reference did not rebuild");
#endif
  }
  long n=kind==2?lag-1:lag;
  run("U",Wavefront::Kernel::U,n,mat,src);
  run("T1",Wavefront::Kernel::T1,n,mat,src);
  auto aligned=run("T2aligned",Wavefront::Kernel::T2,n+1,mat,src);
  if(second.statistics.planar&&aligned.status==Wavefront::Status::Converged){
#if DBEM_WAVEFRONT_REUSE_PARTITIONS
   require(aligned.statistics.partitionCacheHits==1,"aligned T2 window not reused");
#endif
  }
  Wavefront::Options budget;budget.maxEvaluations=1;run("budget",k,lag,mat,src,budget);
  Wavefront::Options tolerance;tolerance.relativeTolerance=1e-5;run("tolerance",k,lag,mat,src,tolerance);
  auto sameSpeed=BuildDynaMat(nu,2*mat.G,2*mat.Rou,dt);
  require(sameSpeed.C1==mat.C1&&sameSpeed.C2==mat.C2,"same-speed material fixture invalid");
  run("sameSpeedMaterial",k,lag,sameSpeed,src);
  auto changed=mat;changed.C1=std::nextafter(mat.C1,INFINITY);run("speed",k,lag,changed,src);
  auto shifted=src;shifted[1]=std::nextafter(0.,1.);auto moved=run("source",k,lag,mat,shifted);
  require(!moved.statistics.partitionCacheHits,"source key was ignored");
  auto time=BuildDynaMat(nu,20000/(2*(1+nu)),8,dt*1.01);auto timed=run("time",k,lag,time,src);
  require(!timed.statistics.partitionCacheHits,"time key was ignored");
  e.GeoInit(p,ids);auto rebuilt=run("geometry",k,lag,mat,src);require(!rebuilt.statistics.partitionCacheHits,"geometry generation was ignored");
  Wavefront::Result parallel[2];std::thread threads[2];
  for(int t=0;t<2;t++)threads[t]=std::thread([&,t]{for(int repeat=0;repeat<2;repeat++)parallel[t]=Wavefront::Integrate(e,k,src,lag,mat);});
  for(auto&t:threads)t.join();
  for(int t=0;t<2;t++){require(parallel[t].block==rebuilt.block&&parallel[t].status==rebuilt.status,"worker result differs");emit(row,t?"thread1":"thread0",parallel[t]);}
  if(row<3){
   p[4].pt[0]=.001;e.GeoInit(p,ids);
   auto curved=run("curved",k,lag,mat,src);run("curvedRepeat",k,lag,mat,src);
   require(!curved.statistics.partitionCacheHits,"curved geometry reused a planar partition");
  }
 }
 std::cerr<<"PASS cases="<<count<<" builds="<<builds<<" hits="<<hits<<" boundsBuilds="<<boundsBuilds<<" boundsHits="<<boundsHits<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 2;}
