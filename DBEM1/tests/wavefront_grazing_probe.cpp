
static int RunDiscretizationKernelProbe() {
 const double uv[8][2]={{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
 FILE*f=0;fopen_s(&f,"grazing_cases.csv","w");
 fprintf(f,"geometry,length,material,shift,kernel,dt,status,error,scale,characteristic,absolute,tolerance,evaluations");
 for(int k=0;k<9;++k)fprintf(f,",b%d",k);fprintf(f,"\n");
 int failures=0;
 for(int geometry=0;geometry<2;++geometry)
 for(double length:{.001,1.,1000.})for(double multiplier:{.0001,1.,10000.}){
  Point points[8];long ids[8]={0,1,2,3,4,5,6,7};
  for(int i=0;i<8;++i){points[i].pt[0]=10*length;points[i].pt[1]=(-.25+.25*uv[i][0])*length;points[i].pt[2]=(.25+.25*uv[i][1])*length;}
  DSquareElement e;e.DomainID=0;e.LocalEleID=22;e.GeoInit(points,ids);
  Wavefront::Vec source{{-3.75*length,(geometry==0?-1.:-.25)*length,.25*length}};
  for(double shift:{-1e-4,-1e-9,-1e-12,0.,1e-12,1e-9,1e-4})
  for(int kind=0;kind<3;++kind){
   double dt=(geometry==0?.06:13.75/(5.*66))*length*(1+shift);
   DynaMat mat=BuildDynaMat(0.,100.*multiplier,8.*multiplier,dt);
   auto r=Wavefront::Integrate(e,(Wavefront::Kernel)kind,source,kind==2?66:65,mat);
   if(r.status==Wavefront::Status::Failed){++failures;fprintf(stderr,"%s\n",r.message.c_str());}
   fprintf(f,"%d,%.17g,%.17g,%.17g,%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%zu",geometry,length,multiplier,shift,kind,dt,(int)r.status,r.error,r.scale,r.characteristicScale,r.absoluteScale,r.tolerance,r.statistics.evaluations);
   for(double b:r.block)fprintf(f,",%.17g",b);fprintf(f,"\n");fflush(f);
  }
 }
 fclose(f);return failures?2:0;
}
