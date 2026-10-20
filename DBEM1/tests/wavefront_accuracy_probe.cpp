static int RunDiscretizationKernelProbe() {
 thread_num=1;
 DSquareElement::StaticInit(.5,0.,10000.,8.,.0005,2);
 AssistUij=new double**[1];AssistTij=new double**[1];
 AssistUij[0]=new double*[8];AssistTij[0]=new double*[8];
 for(int i=0;i<8;++i){AssistUij[0][i]=new double[9];AssistTij[0][i]=new double[9];}
 FILE *f=0,*g=0;fopen_s(&f,"accuracy_blocks.csv","w");fopen_s(&g,"accuracy_geometry.csv","w");
 fprintf(f,"case,kernel,lag,nu,dt,E,active,estimated_error,evaluations");
 for(int k=0;k<9;++k)fprintf(f,",b%d",k);fprintf(f,"\n");
 fprintf(g,"case,node,x,y,z,sx,sy,sz\n");
 const double uv[8][2]={{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
 for(int id=0;id<10;++id){
  double h=(id==0||id==6||id==7)?.1:(id==1?.05:(id==8?.003:(id==9?.25:.03)));
  Point points[8],source;long ids[8]={0,1,2,3,4,5,6,7};
  source.pt[0]=(id<2||id==6)?0.:-.01;source.pt[1]=source.pt[2]=0.;
  if(id==3){source.pt[1]=.06;source.pt[2]=.005;}
  if(id==4)source.pt[0]=-.0001;
  if(id==8)source.pt[0]=-.03;
  if(id==9){source.pt[0]=9.75;source.pt[1]=-.25;source.pt[2]=.75;}
  for(int i=0;i<8;++i){
   double u=uv[i][0],v=uv[i][1];
   points[i].pt[0]=id==6?.01*(u*u+v*v):0.;
   points[i].pt[1]=h*u;points[i].pt[2]=h*v;
   if(id==5)points[i].pt[1]*=1+.2*v;
   if(id==7)points[i].pt[0]=.03*(u*u+v*v);
   if(id==9){points[i].pt[0]=10.;points[i].pt[1]-=.25;points[i].pt[2]+=.25;}
  }
  double nu=id==5?.25:0.;double E=id==9?200.:20000.;double dt=id==9?.06:.0005;
  auto mat=BuildDynaMat(nu,E/(2*(1+nu)),8.,dt);
  auto previous=DSquareElement::SetThreadDynaMat(&mat);
  DSquareElement::m_DMat=mat;
  DSquareElement element;element.DomainID=1;element.LocalEleID=id;element.GeoInit(points,ids);
  for(int i=0;i<8;++i)fprintf(g,"%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",id,i,points[i].pt[0],points[i].pt[1],points[i].pt[2],source.pt[0],source.pt[1],source.pt[2]);
  for(int kind=0;kind<3;++kind){
   if(id<2&&kind)continue;
   long lag=(id==4?(kind==2?1:0):(id==9?(kind==2?3:2):(kind==2?2:1)));
   for(int i=0;i<8;++i)for(int k=0;k<9;++k)AssistUij[0][i][k]=AssistTij[0][i][k]=0.;
   const long slot=0; // Select exactly the overload used in multidomain assembly.
   int active=kind==0?element.IntDynaUijJudge(source,(double)lag,dt,slot):element.IntDynaTijJudge(kind,source,(double)lag,dt,slot);
   double estimated=0.;size_t evaluations=0;
   if(DBEMPolicy::wavefrontIntegration){
    auto check=Wavefront::Integrate(element,(Wavefront::Kernel)kind,Wavefront::Vec{{source.pt[0],source.pt[1],source.pt[2]}},lag,mat);
    if(check.status==Wavefront::Status::Failed)throw std::runtime_error(check.message);
    estimated=check.error;evaluations=check.statistics.evaluations;
    for(int k=0;k<9;++k)if(std::abs(check.block[k]-(kind==0?AssistUij[0][0][k]:AssistTij[0][0][k]))>1e-15)throw std::runtime_error("assembly entry mismatch");
   }
   fprintf(f,"%d,%d,%ld,%.17g,%.17g,%.17g,%d,%.17g,%zu",id,kind,lag,nu,dt,E,active,estimated,evaluations);
   for(int k=0;k<9;++k)fprintf(f,",%.17g",kind==0?AssistUij[0][0][k]:AssistTij[0][0][k]);
   fprintf(f,"\n");fflush(f);
  }
  DSquareElement::SetThreadDynaMat(previous);
 }
 fclose(f);fclose(g);
 for(int i=0;i<8;++i){delete[] AssistUij[0][i];delete[] AssistTij[0][i];}
 delete[] AssistUij[0];delete[] AssistTij[0];delete[] AssistUij;delete[] AssistTij;
 AssistUij=AssistTij=nullptr;
 return 0;
}
