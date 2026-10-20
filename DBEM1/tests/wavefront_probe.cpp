#if defined(DBEM_WAVEFRONT_CACHE_SELF_TEST) && DBEM_WAVEFRONT_CACHE_SELF_TEST
extern int WavefrontCachePolicySelfTest();
namespace Wavefront { int PartitionSelfTest(); }
#endif

static int RunDiscretizationKernelProbe() {
 FILE *f=0,*geo=0;
 fopen_s(&f,"wf_blocks.csv","w");fopen_s(&geo,"wf_geometry.csv","w");
 fprintf(f,"case,kernel,lag,epsilon,dt,nu,status,error,scale,evaluations,planar,depth");
 for(int k=0;k<9;++k)fprintf(f,",b%d",k);fprintf(f,"\n");
 fprintf(geo,"case,node,x,y,z,sx,sy,sz\n");
 const double uv[8][2]={{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
 int failures=0;
 for(int id=0;id<13;++id){
  double h=id==1?.05:(((id>=2&&id<=7)||id==12)?.03:.1);
  Point points[8];long ids[8]={0,1,2,3,4,5,6,7};
  Wavefront::Vec source{{0,0,0}};
  if(id==2||id==12)source={{-.01,0,0}};
  if(id==3)source={{-.01,.06,.005}};
  if(id==4)source={{-.00001,0,0}};
  if(id==5)source={{-.01,.03,0}};
  if(id==6)source={{-.01,.03,.03}};
  if(id==7)source={{-.02,0,0}};
  if(id==9)source={{-.01,0,0}};
  for(int i=0;i<8;++i){
   double u=uv[i][0],v=uv[i][1];
   points[i].pt[0]=(id==8||id==9)?(id==8?.01:.03)*(u*u+v*v):0;
   points[i].pt[1]=h*u;points[i].pt[2]=h*v;
   if(id==10)points[i].pt[2]+=.01*u*u; // coplanar curved edges
   if(id==7){points[i].pt[1]=h*u*(1+.2*v);points[i].pt[2]=h*v;}
  }
  // Rotate and translate a non-rectangular case, including its source.
  if(id==11){
   source={{-.01,.02,-.01}};
   for(int i=0;i<8;++i)points[i].pt[1]*=(1+.2*uv[i][1]);
   auto move=[](double*p){double x=p[0],y=p[1];p[0]=.8*x-.6*y+.2;p[1]=.6*x+.8*y-.3;p[2]+=.4;};
   for(auto&p:points)move(p.pt);move(source.data());
  }
  DSquareElement element;element.DomainID=1;element.LocalEleID=id;
  element.GeoInit(points,ids);
  for(int i=0;i<8;++i)fprintf(geo,"%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",id,i,points[i].pt[0],points[i].pt[1],points[i].pt[2],source[0],source[1],source[2]);
  for(int kernel=0;kernel<3;++kernel){
   long lag=(id==4||id==5||id==6)?(kernel==2?1:0):(kernel==2?2:1);
   for(int p=0;p<((id==0||id==1)?5:1);++p){
    const double shifts[5]={0.,-1e-12,1e-12,-1e-9,1e-9};double shift=shifts[p];
    double nu=id==7?.25:(id==11?.45:0.);
    double dt=(id==12?sqrt(.01*.01+2*h*h)/100:.0005)*(1+shift);
    DynaMat mat=BuildDynaMat(nu,20000./(2*(1+nu)),8.,dt);
    auto value=Wavefront::Integrate(element,(Wavefront::Kernel)kernel,source,lag,mat);
    if(value.status==Wavefront::Status::Failed){++failures;fprintf(stderr,"%s\n",value.message.c_str());}
    fprintf(f,"%d,%d,%ld,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%zu,%d,%u",id,kernel,lag,shift,dt,nu,(int)value.status,value.error,value.scale,value.statistics.evaluations,value.statistics.planar,value.statistics.maxDepth);
    for(double b:value.block)fprintf(f,",%.17g",b);fprintf(f,"\n");fflush(f);
   }
  }
 }

#if defined(DBEM_WAVEFRONT_CACHE_SELF_TEST) && DBEM_WAVEFRONT_CACHE_SELF_TEST
 if(Wavefront::PartitionSelfTest()){++failures;fprintf(stderr,"polar partition area self-test failed\n");}
 if (WavefrontCachePolicySelfTest()) { ++failures; fprintf(stderr,"cache policy self-test failed\n"); }
#endif
 int invariantChecks=0;
 auto check=[&](bool ok,const char*name){++invariantChecks;if(!ok){++failures;fprintf(stderr,"invariant failed: %s\n",name);}};
 Point p0[8],p1[8];long ids[8]={0,1,2,3,4,5,6,7};
 for(int i=0;i<8;++i){
  p0[i].pt[0]=p1[i].pt[0]=0;
  p0[i].pt[1]=.03*uv[i][0];p1[i].pt[1]=-.03*uv[i][0];
  p0[i].pt[2]=p1[i].pt[2]=.03*uv[i][1];
 }
 DSquareElement a,b;a.GeoInit(p0,ids);b.GeoInit(p1,ids);
 auto mat=BuildDynaMat(0.,10000.,8.,.0005);
 Wavefront::Vec source{{-.01,.06,.005}};
 for(int k=0;k<3;++k){
  auto x=Wavefront::Integrate(a,(Wavefront::Kernel)k,source,k==2?2:1,mat);
  auto y=Wavefront::Integrate(b,(Wavefront::Kernel)k,source,k==2?2:1,mat);
  double difference=0,length=0;
  for(int j=0;j<9;++j){double v=x.block[j]-(k? -1:1)*y.block[j];difference+=v*v;length+=x.block[j]*x.block[j];}
  check(x.status==Wavefront::Status::Converged&&y.status==Wavefront::Status::Converged&&sqrt(difference)<=2e-6*sqrt(length)+2e-10*x.scale,"normal reversal");
 }
 Wavefront::Options tight;tight.maxEvaluations=8;
 auto budget=Wavefront::Integrate(a,Wavefront::Kernel::U,source,1,mat,tight);
 check(budget.status==Wavefront::Status::Failed&&budget.message.find("budget")!=std::string::npos,"budget is explicit failure");
 auto zero=Wavefront::Integrate(a,Wavefront::Kernel::U,Wavefront::Vec{{10,0,0}},1,mat);
 check(zero.status==Wavefront::Status::Outside&&zero.statistics.evaluations==0,"outside exact zero");
 zero=Wavefront::Integrate(a,Wavefront::Kernel::T2,source,0,mat);
 check(zero.status==Wavefront::Status::Outside&&zero.statistics.evaluations==0,"T2 causal zero");
 auto singular=Wavefront::Integrate(a,Wavefront::Kernel::U,Wavefront::Vec{},0,mat);
 check(singular.status==Wavefront::Status::Failed&&singular.message.find("singular")!=std::string::npos,"singular cache guard");
 for(auto&p:p1)p.pt[2]=0;
 DSquareElement invalid;invalid.GeoInit(p1,ids);
 auto bad=Wavefront::Integrate(invalid,Wavefront::Kernel::U,source,1,mat);
 check(bad.status==Wavefront::Status::Failed,"invalid geometry is failure");
 check(Wavefront::Classify(invalid,source,1,false,mat)==1,"invalid geometry reaches contextual integration failure");
 // The curved saddle extends outside the sphere of its corner nodes.
 Point saddlePoints[8];
 for(int i=0;i<8;++i){
  double u=uv[i][0],v=uv[i][1];
  saddlePoints[i].pt[0]=.1*(u*u-v*v);saddlePoints[i].pt[1]=.03*u;saddlePoints[i].pt[2]=.03*v;
 }
 DSquareElement saddle;saddle.GeoInit(saddlePoints,ids);
 auto saddleMat=BuildDynaMat(0.,10000.,8.,.0004);
 Wavefront::Vec saddleSource{{-.11,0,.03}};
 check(saddle.m_wavefrontGeometry->valid&&saddle.m_wavefrontGeometry->boundingRadius>.1&&
  Wavefront::Classify(saddle,saddleSource,0,false,saddleMat)==1,"cached Bernstein sphere preserves curved mid-edge support");
 auto other=BuildDynaMat(.25,8000.,12.,.001);
 auto previous=DSquareElement::SetThreadDynaMat(&other);
 auto scoped=Wavefront::Integrate(a,Wavefront::Kernel::U,source,1,mat);
 check(&DSquareElement::ActiveDynaMat()==&other&&scoped.status==Wavefront::Status::Converged,"material scope restored");
 DSquareElement::SetThreadDynaMat(previous);

#if DBEM_POLICY_WAVEFRONT_INTEGRATION
 DSquareElement::m_DMat=mat;
 Point localSource;for(int j=0;j<3;++j)localSource.pt[j]=source[j];
 for(int kind=0;kind<3;++kind){
  auto expected=Wavefront::Integrate(a,(Wavefront::Kernel)kind,source,kind==2?2:1,mat);
  int active=kind==0?a.IntDynaUijJudge(localSource,1.,mat.Dt,0):a.IntDynaTijJudge(kind,localSource,kind==2?2.:1.,mat.Dt,0);
  double difference=0;
  double*actual=kind==0?DSquareElement::m_AssistUij[0]:DSquareElement::m_AssistTij[0];
  for(int k=0;k<9;++k)difference+=(actual[k]-expected.block[k])*(actual[k]-expected.block[k]);
  check(active&&sqrt(difference)<=1e-12*expected.scale,"local and thread-slot backend contract");
 }
#endif
 FILE*report=0;fopen_s(&report,"wf_invariants.json","w");
 fprintf(report,"{\"checks\":%d,\"failures\":%d}\n",invariantChecks,failures);fclose(report);
 fclose(f);fclose(geo);if(RunIndependentPointKernelProbe()!=0)++failures;return failures?2:0;
}
