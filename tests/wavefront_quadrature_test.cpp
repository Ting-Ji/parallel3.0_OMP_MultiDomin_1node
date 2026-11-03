#include "dsquareelement.h"
#include "DBEM.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <cmath>
#include <algorithm>

static void require(bool ok,const char*message){if(!ok)throw std::runtime_error(message);}
static double norm(const std::array<double,9>&a){double v=0;for(double x:a)v+=x*x;return std::sqrt(v);}
int main(int argc,char**argv)try{
 require(argc==2,"reference input required");
 thread_num=2;AssistUij=new double**[2];AssistTij=new double**[2];
 for(int t=0;t<2;t++){AssistUij[t]=new double*[8];AssistTij[t]=new double*[8];for(int i=0;i<8;i++){AssistUij[t][i]=new double[9]{};AssistTij[t][i]=new double[9]{};}}
 std::ifstream in(argv[1]);require(bool(in),"cannot read references");
 int count;in>>count;double maximum=0;int compared=0;
 for(int row=0;row<count;row++){
  std::string name;in>>name;Point points[8];long ids[8]={0,1,2,3,4,5,6,7};
  for(int i=0;i<4;i++){points[i].pt[0]=0;in>>points[i].pt[1]>>points[i].pt[2];}
  for(int i=4;i<8;i++)for(int k=0;k<3;k++)points[i].pt[k]=(points[i-4].pt[k]+points[(i-3)%4].pt[k])/2;
  double gap,nu,dt,natural;int lag,kind;in>>gap>>nu>>dt>>lag>>kind>>natural;
  std::array<double,9>ref{};for(double&v:ref)in>>v;require(bool(in),"incomplete fixture");
  DSquareElement e;e.GeoInit(points,ids);if(!e.UsesPlanarWavefront()){std::cerr<<name<<" geometry="<<e.m_wavefrontGeometry->reason<<'\n';throw std::runtime_error("planar integration not enabled");}
  auto mat=BuildDynaMat(nu,20000./(2*(1+nu)),8.,dt);Point source;source.pt[0]=-gap;source.pt[1]=source.pt[2]=0;
  Wavefront::Vec src{{-gap,0,0}};auto k=(Wavefront::Kernel)kind;
  auto value=Wavefront::Integrate(e,k,src,lag,mat);
  if(value.status==Wavefront::Status::Failed)throw std::runtime_error(value.message);
  std::array<double,9>delta{};for(int j=0;j<9;j++)delta[j]=value.block[j]-ref[j];
  double allowed=1e-6*norm(ref)+1e-10*natural,ratio=allowed>0?norm(delta)/allowed:0;
  maximum=(std::max)(maximum,ratio);if(ratio>1.1){std::cerr<<name<<" ratio="<<ratio<<'\n';throw std::runtime_error("independent kernel reference mismatch");}
  auto old=e.SetThreadDynaMat(&mat);
  if(kind==0)e.IntDynaUijJudge(source,lag,dt);else e.IntDynaTijJudge(kind,source,lag,dt);
  auto serial=kind==0?DSquareElement::m_AssistUij:DSquareElement::m_AssistTij;
  for(int j=0;j<9;j++)require(serial[0][j]==value.block[j],"serial judge did not use exact integral");
  // First-window on-panel singular terms are supplied by the caller's cache.
  if(kind==1&&(lag>1||gap>0)){
   auto second=Wavefront::Integrate(e,Wavefront::Kernel::T2,src,lag,mat);
   require(second.status!=Wavefront::Status::Failed,"T2 for combined entry failed");
   e.IntDynaTijJudge(source,lag,dt);
   for(int j=0;j<9;j++)require(DSquareElement::m_AssistTij[0][j]==value.block[j]+second.block[j],"combined T entry mismatch");
  }
  e.SetThreadDynaMat(old);
  std::exception_ptr failures[2];std::thread workers[2];
  for(long t=0;t<2;t++)workers[t]=std::thread([&,t](){try{auto previous=e.SetThreadDynaMat(&mat);if(kind==0)e.IntDynaUijJudge(source,lag,dt,t);else e.IntDynaTijJudge(kind,source,lag,dt,t);e.SetThreadDynaMat(previous);}catch(...){failures[t]=std::current_exception();}});
  for(auto&t:workers)t.join();for(auto&ex:failures)if(ex)std::rethrow_exception(ex);
  for(int t=0;t<2;t++)for(int j=0;j<9;j++)require((kind==0?AssistUij:AssistTij)[t][0][j]==value.block[j],"parallel judge mismatch");
  if(row==0){Wavefront::Options tiny;tiny.maxEvaluations=1;auto bad=Wavefront::Integrate(e,k,src,lag,mat,tiny);require(bad.status==Wavefront::Status::Failed&&!bad.message.empty(),"budget failure was not reported");}
  compared++;
 }
 std::cout<<"PASS: "<<compared<<" independent U/T1/T2 blocks; serial/2-thread exact match; budget failure; max mixed ratio="<<maximum<<'\n';
 return 0;
}catch(const std::exception&ex){std::cerr<<"FAIL: "<<ex.what()<<'\n';return 2;}
