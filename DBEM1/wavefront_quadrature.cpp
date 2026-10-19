#include "wavefront_quadrature.h"
#include "dsquareelement.h"
#include "solver_policy.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <mutex>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Wavefront {
namespace {
constexpr double pi=3.1415926535897932384626433832795;
constexpr double eps=std::numeric_limits<double>::epsilon();
using Block=std::array<double,9>;
Vec add(const Vec&a,const Vec&b){return {{a[0]+b[0],a[1]+b[1],a[2]+b[2]}};}
Vec sub(const Vec&a,const Vec&b){return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};}
Vec mul(const Vec&a,double s){return {{a[0]*s,a[1]*s,a[2]*s}};}
double dot(const Vec&a,const Vec&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Vec cross(const Vec&a,const Vec&b){return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};}
double norm(const Vec&a){return std::sqrt(dot(a,a));}
double norm(const Block&a){double s=0;for(double x:a)s+=x*x;return std::sqrt(s);}
Vec point(const Geometry&g,double u,double v){
 const double t[8]={1,u,v,u*v,u*u,v*v,u*u*v,u*v*v};Vec p{};
 for(int k=0;k<3;++k)for(int j=0;j<8;++j)p[k]+=g.coefficient[k][j]*t[j];
 return p;
}
Vec du(const Geometry&g,double u,double v){
 Vec r{};for(int k=0;k<3;++k){const double*c=g.coefficient[k];
 r[k]=c[1]+c[3]*v+2*c[4]*u+2*c[6]*u*v+c[7]*v*v;}return r;
}
Vec dv(const Geometry&g,double u,double v){
 Vec r{};for(int k=0;k<3;++k){const double*c=g.coefficient[k];
 r[k]=c[2]+c[3]*u+2*c[5]*v+c[6]*u*u+2*c[7]*u*v;}return r;
}
struct Rule {std::vector<double>x,w;};
Rule rule(unsigned n){
 Rule r;r.x.resize(n);r.w.resize(n);
 for(unsigned i=0;i<(n+1)/2;++i){
  double z=std::cos(pi*(i+.75)/(n+.5)),previous,derivative=0;
  for(unsigned iteration=0;iteration<30;++iteration){
   double p=1,q=0;
   for(unsigned j=1;j<=n;++j){double old=q;q=p;p=((2*j-1)*z*q-(j-1)*old)/j;}
   derivative=n*(z*p-q)/(z*z-1);previous=z;z-=p/derivative;
   if(std::abs(z-previous)<=2*eps)break;
  }
  // Recompute the derivative at the converged node.
  double p=1,q=0;for(unsigned j=1;j<=n;++j){double old=q;q=p;p=((2*j-1)*z*q-(j-1)*old)/j;}
  derivative=n*(z*p-q)/(z*z-1);
  double w=2/((1-z*z)*derivative*derivative);
  r.x[i]=-z;r.x[n-1-i]=z;r.w[i]=w;r.w[n-1-i]=w;
 }return r;
}
struct Rules{Rule low=rule(5),high=rule(9);};
const Rules& rules(){static const Rules r;return r;}
struct MatScope{
 const DynaMat* previous;
 explicit MatScope(const DynaMat& m):previous(DSquareElement::SetThreadDynaMat(&m)){}
 ~MatScope(){DSquareElement::SetThreadDynaMat(previous);}
};
struct Cell{
 double a,b,c,d; // parameter rectangle, or angle and unit radial coordinate
 double radialLow=0,radialHigh=0;
 unsigned depth=0;std::size_t id=0;
 Block value{};double error=0,scale=0,missed=0;bool forced=false;
};
struct Worse{
 bool operator()(const Cell&a,const Cell&b)const{
  if(a.forced!=b.forced)return !a.forced;
  if(a.error!=b.error)return a.error<b.error;
  return a.id>b.id;
 }
};

struct Bounds{double low,high,area;bool regular;};
Bounds bounds(const Geometry&g,const Vec&source,const Cell&c){
 Vec p[3][3],q[3][3],control[3][3];
 for(int i=0;i<3;++i)for(int j=0;j<3;++j)
  p[i][j]=sub(point(g,c.a+(c.b-c.a)*i/2,c.c+(c.d-c.c)*j/2),source);
 for(int j=0;j<3;++j){
  q[0][j]=p[0][j];q[2][j]=p[2][j];
  q[1][j]=sub(mul(p[1][j],2),mul(add(p[0][j],p[2][j]),.5));
 }
 for(int i=0;i<3;++i){
  control[i][0]=q[i][0];control[i][2]=q[i][2];
  control[i][1]=sub(mul(q[i][1],2),mul(add(q[i][0],q[i][2]),.5));
 }
 Vec lo=control[0][0],hi=lo;double magnitude=norm(source)+g.diameter;
 for(auto& row:control)for(auto& v:row)for(int k=0;k<3;++k){lo[k]=(std::min)(lo[k],v[k]);hi[k]=(std::max)(hi[k],v[k]);magnitude=(std::max)(magnitude,std::abs(v[k]));}
 double lower=0,upper=0,pad=64*eps*(magnitude+g.diameter);
 for(int k=0;k<3;++k){
  lo[k]=std::nextafter(lo[k]-pad,-std::numeric_limits<double>::infinity());
  hi[k]=std::nextafter(hi[k]+pad,std::numeric_limits<double>::infinity());
  double near=lo[k]>0?lo[k]:(hi[k]<0?-hi[k]:0);
  double far=(std::max)(std::abs(lo[k]),std::abs(hi[k]));
  lower+=near*near;upper+=far*far;
 }
 double maxU=0,maxV=0;
 Vec uLow{{1e300,1e300,1e300}},uHigh{{-1e300,-1e300,-1e300}},vLow=uLow,vHigh=uHigh;
 for(int i=0;i<2;++i)for(int j=0;j<3;++j){
  Vec value=mul(sub(control[i+1][j],control[i][j]),2);
  for(int k=0;k<3;++k){uLow[k]=(std::min)(uLow[k],value[k]);uHigh[k]=(std::max)(uHigh[k],value[k]);}
 }
 for(int i=0;i<3;++i)for(int j=0;j<2;++j){
  Vec value=mul(sub(control[i][j+1],control[i][j]),2);
  for(int k=0;k<3;++k){vLow[k]=(std::min)(vLow[k],value[k]);vHigh[k]=(std::max)(vHigh[k],value[k]);}
 }
 double derivativePad=128*eps*(magnitude+g.diameter),jacobianLowerSquared=0;
 for(int k=0;k<3;++k){uLow[k]-=derivativePad;uHigh[k]+=derivativePad;vLow[k]-=derivativePad;vHigh[k]+=derivativePad;}
 auto intervalProduct=[](double a,double b,double c,double d){
  std::array<double,4> v{{a*c,a*d,b*c,b*d}};
  return std::array<double,2>{{*(std::min_element)(v.begin(),v.end()),*(std::max_element)(v.begin(),v.end())}};
 };
 for(int k=0;k<3;++k){
  int i=(k+1)%3,j=(k+2)%3;
  auto a=intervalProduct(uLow[i],uHigh[i],vLow[j],vHigh[j]);
  auto b=intervalProduct(uLow[j],uHigh[j],vLow[i],vHigh[i]);
  double lo=a[0]-b[1],hi=a[1]-b[0],near=lo>0?lo:(hi<0?-hi:0);
  jacobianLowerSquared+=near*near;
 }
 bool regular=std::sqrt(jacobianLowerSquared)>1e-14*g.diameter*g.diameter*(c.b-c.a)*(c.d-c.c);
 for(int i=0;i<2;++i)for(int j=0;j<3;++j)maxU=(std::max)(maxU,2*norm(sub(control[i+1][j],control[i][j])));
 for(int i=0;i<3;++i)for(int j=0;j<2;++j)maxV=(std::max)(maxV,2*norm(sub(control[i][j+1],control[i][j])));
 // Bernstein products preserve a tight squared-distance hull without
 // inverting a Vandermonde matrix on very small cells.
 double squared[5][5]{},maximumTerm=0;
 const double choose2[3]={1,2,1},choose4[5]={1,4,6,4,1};
 for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int a=0;a<3;++a)for(int b=0;b<3;++b){
  double term=dot(control[i][j],control[a][b])*
   choose2[i]*choose2[a]/choose4[i+a]*choose2[j]*choose2[b]/choose4[j+b];
  squared[i+a][j+b]+=term;maximumTerm=(std::max)(maximumTerm,std::abs(term));
 }
 double squaredLow=std::numeric_limits<double>::infinity(),squaredHigh=0;
 for(auto& row:squared)for(double v:row){squaredLow=(std::min)(squaredLow,v);squaredHigh=(std::max)(squaredHigh,v);}
 double squaredPad=128*eps*(upper+maximumTerm)+2*pad*std::sqrt(upper);
 lower=(std::max)(lower,(std::max)(0.,squaredLow-squaredPad));
 upper=(std::min)(upper,(std::max)(0.,squaredHigh+squaredPad));
 return {std::sqrt(lower),std::sqrt(upper),maxU*maxV*(1+128*eps),regular};
}
}
namespace {
std::mutex failureMutex;
std::string workerFailure;
}
void RecordWorkerFailure(const char*message) noexcept {
 try {std::lock_guard<std::mutex> lock(failureMutex);if(workerFailure.empty())workerFailure=message;}
 catch (...) {std::fputs("Unable to record integration worker failure\n",stderr);std::abort();}
}
void ThrowWorkerFailure() {
 std::string message;
 {std::lock_guard<std::mutex> lock(failureMutex);message.swap(workerFailure);}
 if(!message.empty())throw std::runtime_error(message);
}
Options::Options():relativeTolerance(DBEM_WAVEFRONT_REL_TOL),
 scaledAbsoluteTolerance(DBEM_WAVEFRONT_ABS_SCALE),maxDepth(DBEM_WAVEFRONT_MAX_DEPTH),
 maxEvaluations(DBEM_WAVEFRONT_MAX_EVALUATIONS){}
const char*Name(Kernel k){return k==Kernel::U?"U":(k==Kernel::T1?"T1":"T2");}

std::shared_ptr<const Geometry> BuildGeometry(DSquareElement&e){
 auto g=std::make_shared<Geometry>();
 for(int k=0;k<3;++k)for(int j=0;j<8;++j)g->coefficient[k][j]=e.m_CoordCoef[k][j];
 const double uv[4][2]={{-1,-1},{1,-1},{1,1},{-1,1}};
 for(int i=0;i<4;++i)g->vertices[i]=point(*g,uv[i][0],uv[i][1]);
 g->origin=point(*g,0,0);
 for(auto&v:g->vertices)g->diameter=(std::max)(g->diameter,2*norm(sub(v,g->origin)));
 if(!(g->diameter>0)||!std::isfinite(g->diameter)){g->reason="zero geometry diameter";return g;}
 for(int k=0;k<3;++k)for(int j=0;j<8;++j)if(!std::isfinite(g->coefficient[k][j])){g->reason="non-finite geometry";return g;}
 Vec n=cross(du(*g,0,0),dv(*g,0,0));double length=norm(n);
 if(!(length>1e-14*g->diameter*g->diameter)){g->reason="degenerate central Jacobian";return g;}
 g->normal=mul(n,1/length);
 for(double u:{-1.,0.,1.})for(double v:{-1.,0.,1.})
  if(!(norm(cross(du(*g,u,v),dv(*g,u,v)))>1e-14*g->diameter*g->diameter)){
   g->reason="degenerate geometric mapping";return g;
  }
 g->valid=true;g->planar=true;
 double tol=1e-12*g->diameter;
 for(int i=0;i<8;++i){
  const Point&p=e.m_pointlist[e.m_pointID[i]];Vec v{{p.pt[0],p.pt[1],p.pt[2]}};
  if(std::abs(dot(sub(v,g->origin),g->normal))>tol)g->planar=false;
  if(i>=4&&norm(sub(v,mul(add(g->vertices[i-4],g->vertices[(i-3)%4]),.5)))>tol)g->planar=false;
 }
 g->axis1=sub(g->vertices[1],g->vertices[0]);
 if(!(norm(g->axis1)>tol)){g->valid=false;g->reason="zero polygon edge";return g;}
 g->axis1=mul(g->axis1,1/norm(g->axis1));g->axis2=cross(g->normal,g->axis1);
 for(int i=0;i<4;++i)
  if(dot(cross(sub(g->vertices[(i+1)%4],g->vertices[i]),sub(g->vertices[(i+2)%4],g->vertices[(i+1)%4])),g->normal)<=tol*g->diameter)
   g->planar=false;
 // Certify a nonzero Jacobian over the complete parameter domain, including
 // between the input nodes. An uncertified fold is rejected, never sampled away.
 std::vector<Cell> pending(1);pending[0].a=pending[0].c=-1;pending[0].b=pending[0].d=1;
 std::size_t checked=0;
 while(!pending.empty()){
  Cell c=pending.back();pending.pop_back();
  if(bounds(*g,Vec{},c).regular)continue;
  if(c.depth>=20||++checked>200000){g->valid=false;g->reason="cannot certify a nondegenerate geometric mapping";return g;}
  double u=(c.a+c.b)/2,v=(c.c+c.d)/2;
  for(int i=0;i<2;++i)for(int j=0;j<2;++j){
   Cell child;child.a=i?u:c.a;child.b=i?c.b:u;child.c=j?v:c.c;child.d=j?c.d:v;child.depth=c.depth+1;pending.push_back(child);
  }
 }
 (void)rules();return g;
}
namespace {

double polynomial(const std::vector<double>&c,double x){
 double value=0;for(auto i=c.rbegin();i!=c.rend();++i)value=value*x+*i;return value;
}
// Isolate real roots on the bounded parameter interval by derivative roots.
// On each resulting monotone interval a sign change has exactly one root.
void roots(std::vector<double>c,double low,double high,std::vector<double>&out){
 while(c.size()>1&&c.back()==0)c.pop_back();
 if(c.size()<2)return;
 double scale=0;for(double x:c)scale=(std::max)(scale,std::abs(x));
 if(!(scale>0))return;for(double&x:c)x/=scale;
 if(c.size()==2){double x=-c[0]/c[1];if(x>low&&x<high)out.push_back(x);return;}
 std::vector<double>derivative(c.size()-1),knots{low};
 for(std::size_t i=1;i<c.size();++i)derivative[i-1]=i*c[i];
 roots(derivative,low,high,knots);knots.push_back(high);std::sort(knots.begin(),knots.end());
 for(std::size_t i=0;i<knots.size();++i){
  double x=knots[i],fx=polynomial(c,x);
  if(i&&i+1<knots.size()&&std::abs(fx)<=64*eps)out.push_back(x);
  if(i==0)continue;
  double a=knots[i-1],b=x,fa=polynomial(c,a),fb=fx;
  if((fa<0)==(fb<0)||fa==0||fb==0)continue;
  for(int iteration=0;iteration<64;++iteration){
   double mid=(a+b)/2,fm=polynomial(c,mid);
   if(mid==a||mid==b)break;
   if((fa<0)==(fm<0)){a=mid;fa=fm;}else{b=mid;fb=fm;}
  }
  double root=(a+b)/2;if(root>low&&root<high)out.push_back(root);
 }
}

using Polynomial=std::vector<double>;
Polynomial product(const Polynomial&a,const Polynomial&b){
 Polynomial r(a.size()+b.size()-1,0);
 for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<b.size();++j)r[i+j]+=a[i]*b[j];
 return r;
}
void term(Polynomial&r,double scale,std::initializer_list<std::pair<const Polynomial*,unsigned>>factors){
 Polynomial value{scale};
 for(auto& factor:factors)for(unsigned j=0;j<factor.second;++j)value=product(value,*factor.first);
 if(r.size()<value.size())r.resize(value.size());
 for(std::size_t i=0;i<value.size();++i)r[i]+=value[i];
}
Polynomial dotPolynomial(const std::vector<Vec>&a,const std::vector<Vec>&b){
 Polynomial r(a.size()+b.size()-1,0);
 for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<b.size();++j)r[i+j]+=dot(a[i],b[j]);
 return r;
}
void trim(Polynomial&p){
 double m=0;for(double v:p)m=(std::max)(m,std::abs(v));
 while(p.size()>1&&std::abs(p.back())<=256*eps*m)p.pop_back();
}
std::vector<double> criticalParameters(const Geometry&g,const Vec&source,const std::vector<double>&fronts){
 // The discriminant of |X(x,y)-source|^2-R^2 in y vanishes at
 // every interior y-tangency. Split outer x integration there as well as at
 // intersections with cell edges, so two Gauss rules cannot miss the same cap.
 std::vector<Vec>A(3),B(3),C(2);
 for(int k=0;k<3;++k){
  A[0][k]=g.coefficient[k][0]-source[k];A[1][k]=g.coefficient[k][1];A[2][k]=g.coefficient[k][4];
  B[0][k]=g.coefficient[k][2];B[1][k]=g.coefficient[k][3];B[2][k]=g.coefficient[k][6];
  C[0][k]=g.coefficient[k][5];C[1][k]=g.coefficient[k][7];
 }
 auto a=dotPolynomial(C,C),b=dotPolynomial(B,C),c=dotPolynomial(B,B),d=dotPolynomial(A,B);
 auto ac=dotPolynomial(A,C);for(std::size_t i=0;i<ac.size();++i)c[i]+=2*ac[i];
 for(double&v:b)v*=2;for(double&v:d)v*=2;
 auto e0=dotPolynomial(A,A);std::vector<double>out;
 bool quartic=false,cubic=false;for(double v:a)quartic|=v!=0;for(double v:b)cubic|=v!=0;
 for(double front:fronts){
  auto e=e0;e[0]-=front*front;Polynomial D;
  if(quartic){
   term(D,256,{{&a,3},{&e,3}});term(D,-192,{{&a,2},{&b,1},{&d,1},{&e,2}});
   term(D,-128,{{&a,2},{&c,2},{&e,2}});term(D,144,{{&a,2},{&c,1},{&d,2},{&e,1}});
   term(D,-27,{{&a,2},{&d,4}});term(D,144,{{&a,1},{&b,2},{&c,1},{&e,2}});
   term(D,-6,{{&a,1},{&b,2},{&d,2},{&e,1}});term(D,-80,{{&a,1},{&b,1},{&c,2},{&d,1},{&e,1}});
   term(D,18,{{&a,1},{&b,1},{&c,1},{&d,3}});term(D,16,{{&a,1},{&c,4},{&e,1}});
   term(D,-4,{{&a,1},{&c,3},{&d,2}});term(D,-27,{{&b,4},{&e,2}});
   term(D,18,{{&b,3},{&c,1},{&d,1},{&e,1}});term(D,-4,{{&b,3},{&d,3}});
   term(D,-4,{{&b,2},{&c,3},{&e,1}});term(D,1,{{&b,2},{&c,2},{&d,2}});
  }else if(cubic){
   term(D,1,{{&c,2},{&d,2}});term(D,-4,{{&b,1},{&d,3}});
   term(D,-4,{{&c,3},{&e,1}});term(D,-27,{{&b,2},{&e,2}});
   term(D,18,{{&b,1},{&c,1},{&d,1},{&e,1}});
  }else{term(D,1,{{&d,2}});term(D,-4,{{&c,1},{&e,1}});}
  trim(D);roots(D,-1.,1.,out);
 }
 std::sort(out.begin(),out.end());out.erase(std::unique(out.begin(),out.end(),[](double a,double b){return std::abs(a-b)<=16*eps;}),out.end());
 return out;
}
struct Context{
 DSquareElement&e;Kernel kind;long lag;const DynaMat&mat;const Options&opt;
 const Geometry&g;Vec source;Result&result;
 double lowerTime,upperTime,inner,outer,signedGap=0,nearDistance=0;
 std::vector<double> fronts,critical;
 std::array<std::array<double,2>,4> polygon{};
 bool planar;
 Context(DSquareElement&ee,Kernel k,long l,const DynaMat&m,const Options&o,const Geometry&gg,const Vec&s,Result&r)
 :e(ee),kind(k),lag(l),mat(m),opt(o),g(gg),source(s),result(r),planar(gg.planar){
  lowerTime=(k==Kernel::T2?l-1:l)*m.Dt;upperTime=lowerTime+m.Dt;
  inner=m.C2*(std::max)(0.,lowerTime);outer=m.C1*upperTime;
  for(double c:{m.C1,m.C2})for(double t:{lowerTime,upperTime})if(t>0)fronts.push_back(c*t);
  std::sort(fronts.begin(),fronts.end());fronts.erase(std::unique(fronts.begin(),fronts.end()),fronts.end());
 }
 bool outside(const Bounds&b)const{return b.high<inner||b.low>=outer;}
 Block kernel(const Vec&v,const Vec&normal){
  if(result.statistics.evaluations>=opt.maxEvaluations)throw std::runtime_error("kernel evaluation budget exceeded");
  ++result.statistics.evaluations;
  double r=norm(v);
  if(!(r>0))throw std::runtime_error("singular point requires the cached singular branch");
  double R[6]={0,r,r*r,r*r*r,r*r*r*r,r*r*r*r*r};
  double RI[4]={0,v[0],v[1],v[2]},n[3]={normal[0],normal[1],normal[2]};
  Block block{};
  if(kind==Kernel::U)e.GetDynaUij(block.data(),(double)lag,mat.Dt,R,RI);
  else e.GetDynaTij(kind==Kernel::T1?1:2,block.data(),(double)lag,mat.Dt,R,RI,n);
  for(double value:block)if(!std::isfinite(value))throw std::runtime_error("non-finite kernel");
  return block;
 }
 bool ray(double theta,double&low,double&high)const{
  double dx=std::cos(theta),dy=std::sin(theta);low=0;high=std::numeric_limits<double>::infinity();
  for(int i=0;i<4;++i){
   const auto&a=polygon[i];const auto&b=polygon[(i+1)%4];
   double ex=b[0]-a[0],ey=b[1]-a[1],den=ex*dy-ey*dx,constant=ex*a[1]-ey*a[0];
   if(std::abs(den)<=8*eps*std::hypot(ex,ey)){if(constant>0)return false;continue;}
   double t=constant/den;
   if(den>0)low=(std::max)(low,t);else high=(std::min)(high,t);
  }
  return high>low&&high>0&&std::isfinite(high);
 }
 Block sample(const Cell&c,double x,double y,double&jac){
  if(planar){
   double lo,hi;if(!ray(x,lo,hi)){jac=0;return {};}
   lo=(std::max)(lo,c.radialLow);hi=(std::min)(hi,c.radialHigh);
   if(!(hi>lo)){jac=0;return {};}
   double rho=lo+(hi-lo)*y;
   Vec v=add(mul(g.normal,-signedGap),add(mul(g.axis1,rho*std::cos(x)),mul(g.axis2,rho*std::sin(x))));
   jac=rho*(hi-lo);return kernel(v,g.normal);
  }
  Vec v=sub(point(g,x,y),source),normal=cross(du(g,x,y),dv(g,x,y));jac=norm(normal);
  if(!(jac>1e-14*g.diameter*g.diameter))throw std::runtime_error("invalid sampled Jacobian");
  return kernel(v,mul(normal,1/jac));
 }
 Block integrateRule(const Cell&c,const Rule&r,double&scale){
  Block sum{};scale=0;
  if(!planar){
   // Retain parameter-cell adaptivity, but split the inner coordinate at all
   // actual distance-level roots. Pure rectangular refinement is too slow for
   // a 1e-6 target on a discontinuous curved-panel integrand.
   std::vector<double>outerCuts{c.a,c.b};
   for(double x:critical)if(x>c.a&&x<c.b)outerCuts.push_back(x);
   for(double y:{c.c,c.d}){
    Vec a=sub(point(g,0,y),source),b=du(g,0,y),d{};
    for(int k=0;k<3;++k)d[k]=g.coefficient[k][4]+g.coefficient[k][6]*y;
    for(double front:fronts){
     std::vector<double>poly{dot(a,a)-front*front,2*dot(a,b),dot(b,b)+2*dot(a,d),2*dot(b,d),dot(d,d)};
     roots(poly,c.a,c.b,outerCuts);
    }
   }
   std::sort(outerCuts.begin(),outerCuts.end());
   outerCuts.erase(std::unique(outerCuts.begin(),outerCuts.end(),[](double a,double b){return std::abs(a-b)<=8*eps;}),outerCuts.end());
   for(std::size_t outer=1;outer<outerCuts.size();++outer)for(std::size_t i=0;i<r.x.size();++i){
    double width=outerCuts[outer]-outerCuts[outer-1];
    double x=(outerCuts[outer]+outerCuts[outer-1])/2+width*r.x[i]/2;
    Vec a=sub(point(g,x,0),source),b=dv(g,x,0),d{};
    for(int k=0;k<3;++k)d[k]=g.coefficient[k][5]+g.coefficient[k][7]*x;
    std::vector<double>cuts{c.c,c.d};
    for(double front:fronts){
     std::vector<double>poly{dot(a,a)-front*front,2*dot(a,b),dot(b,b)+2*dot(a,d),2*dot(b,d),dot(d,d)};
     roots(poly,c.c,c.d,cuts);
    }
    std::sort(cuts.begin(),cuts.end());
    cuts.erase(std::unique(cuts.begin(),cuts.end(),[](double a,double b){return std::abs(a-b)<=8*eps;}),cuts.end());
    for(std::size_t t=1;t<cuts.size();++t)for(std::size_t j=0;j<r.x.size();++j){
     double y=(cuts[t]+cuts[t-1])/2+(cuts[t]-cuts[t-1])*r.x[j]/2,jac;
     Block value=sample(c,x,y,jac);
     double w=r.w[i]*r.w[j]*width*(cuts[t]-cuts[t-1])*jac/4;
     for(int k=0;k<9;++k)sum[k]+=w*value[k];
     scale+=std::abs(w)*norm(value);
    }
   }return sum;
  }
  for(std::size_t i=0;i<r.x.size();++i)for(std::size_t j=0;j<r.x.size();++j){
   double x=(c.a+c.b)/2+(c.b-c.a)*r.x[i]/2,y=(c.c+c.d)/2+(c.d-c.c)*r.x[j]/2,jac;
   Block value=sample(c,x,y,jac);
   double w=r.w[i]*r.w[j]*(c.b-c.a)*(c.d-c.c)*jac/4;
   for(int k=0;k<9;++k)sum[k]+=w*value[k];
   scale+=std::abs(w)*norm(value);
  }return sum;
 }
 void evaluate(Cell&c){
  result.statistics.maxDepth=(std::max)(result.statistics.maxDepth,c.depth);
  bool crossing=false;Bounds b{};
  if(!planar){
   b=bounds(g,source,c);
   if(outside(b)){++result.statistics.skipped;return;}
   for(double front:fronts)if(b.low<front&&front<b.high)crossing=true;
  }
  double lowScale=0;Block low=integrateRule(c,rules().low,lowScale);
  c.value=integrateRule(c,rules().high,c.scale);
  Block difference{};for(int k=0;k<9;++k)difference[k]=c.value[k]-low[k];
  c.error=norm(difference);c.forced=!planar&&((crossing&&c.depth<2)||!b.regular);
  // A thin support region may be missed by both Gauss rules. Its possible
  // contribution must be bounded rather than treating two zeros as convergence.
  if(!planar&&c.scale==0&&lowScale==0){
   double r=(std::max)(b.low,inner);
   if(!(r>0))c.forced=true;
   else{
    double bound;
    if(kind==Kernel::U)bound=5/(4*pi*mat.Rou*mat.C2*mat.C2*r);
    else{
     double lambda=2*mat.G*mat.v/(1-2*mat.v);
     // U_weighted=A(r)I+B(r)vv^T. For a unit linear time hat,
     // ||grad U|| <= [40/(cs^2 r^2)+4/(cs^3 r dt)]/(4*pi*rho).
     // Contraction with the traction operator adds sqrt(3)|lambda|+2G.
     bound=(std::sqrt(3.)*std::abs(lambda)+2*mat.G)/(4*pi*mat.Rou)*
      (40/(mat.C2*mat.C2*r*r)+4/(mat.C2*mat.C2*mat.C2*r*mat.Dt));
    }
    c.missed=b.area*bound;c.error=(std::max)(c.error,c.missed);
   }
  }
  if(!std::isfinite(c.error)||!std::isfinite(c.scale))throw std::runtime_error("non-finite quadrature estimate");
 }
};

std::vector<Cell> planeCells(Context&ctx){
 const Geometry&g=ctx.g;ctx.signedGap=dot(sub(ctx.source,g.origin),g.normal);
 double gap=std::abs(ctx.signedGap);
 std::vector<double> radii{0},angles{0,2*pi};
 bool inside=true;double closest=std::numeric_limits<double>::infinity(),maximum=0;
 for(int i=0;i<4;++i){
  Vec v=sub(g.vertices[i],ctx.source);
  ctx.polygon[i]={{dot(v,g.axis1),dot(v,g.axis2)}};
  maximum=(std::max)(maximum,std::hypot(ctx.polygon[i][0],ctx.polygon[i][1]));
 }
 for(int i=0;i<4;++i){
  auto a=ctx.polygon[i],b=ctx.polygon[(i+1)%4];
  double ex=b[0]-a[0],ey=b[1]-a[1],len2=ex*ex+ey*ey;
  if(ex*(-a[1])-ey*(-a[0])<0)inside=false;
  double t=-(a[0]*ex+a[1]*ey)/len2;t=(std::max)(0.,(std::min)(1.,t));
  closest=(std::min)(closest,std::hypot(a[0]+t*ex,a[1]+t*ey));
 }
 if(inside)closest=0;
 ctx.nearDistance=std::hypot(gap,closest);
 if(ctx.lowerTime==0&&inside&&gap<=64*eps*g.diameter)
  throw std::runtime_error("self singular integral must use the existing initialization cache");
 if(gap>=ctx.outer)return {};
 for(double front:ctx.fronts)if(front>gap){
  double rho=std::sqrt((front-gap)*(front+gap));radii.push_back(rho);
 }
 double outer=std::sqrt((ctx.outer-gap)*(ctx.outer+gap));
 radii.push_back(outer);
 // Resolve near-singular regular blocks instead of letting both rules miss a
 // peak narrower than their first radial node.
 if(ctx.lowerTime==0&&ctx.nearDistance>0)
  for(double rho=ctx.nearDistance;rho<outer&&std::isfinite(rho);rho*=2)radii.push_back(rho);
 std::sort(radii.begin(),radii.end());radii.erase(std::unique(radii.begin(),radii.end()),radii.end());
 auto angle=[&](double x,double y){double a=std::atan2(y,x);if(a<0)a+=2*pi;angles.push_back(a);};
 for(auto&v:ctx.polygon)angle(v[0],v[1]);
 // Topology changes only at vertices and circle/edge intersections.
 for(int i=0;i<4;++i){
  auto a=ctx.polygon[i],b=ctx.polygon[(i+1)%4];
  double ex=b[0]-a[0],ey=b[1]-a[1],A=ex*ex+ey*ey,B=a[0]*ex+a[1]*ey;
  for(double radius:radii)if(radius>0&&radius<=maximum*(1+64*eps)){
   double C=a[0]*a[0]+a[1]*a[1]-radius*radius,D=B*B-A*C;
   double tol=64*eps*(B*B+std::abs(A*C)+A*radius*radius);
   if(D< -tol)continue;D=std::sqrt((std::max)(0.,D));
   for(double t:{(-B-D)/A,(-B+D)/A})
    if(t>=0&&t<=1)angle(a[0]+t*ex,a[1]+t*ey);
  }
 }
 std::sort(angles.begin(),angles.end());
 angles.erase(std::unique(angles.begin(),angles.end(),[](double a,double b){return std::abs(a-b)<=16*eps;}),angles.end());
 if(angles.back()!=2*pi)angles.push_back(2*pi);
 std::vector<Cell> cells;
 double inner=ctx.inner>gap?std::sqrt((ctx.inner-gap)*(ctx.inner+gap)):0;
 for(std::size_t i=1;i<angles.size();++i){
  if(!(angles[i]>angles[i-1]))continue;
  double low,high;if(!ctx.ray((angles[i]+angles[i-1])/2,low,high))continue;
  for(std::size_t j=1;j<radii.size();++j){
   double a=(std::max)(radii[j-1],inner),b=(std::min)(radii[j],outer);
   if(!(b>a)||high<=a||low>=b)continue;
   Cell cell;cell.a=angles[i-1];cell.b=angles[i];cell.c=0;cell.d=1;cell.radialLow=a;cell.radialHigh=b;
   cells.push_back(cell);
  }
 }return cells;
}
}

#if defined(DBEM_WAVEFRONT_CACHE_SELF_TEST) && DBEM_WAVEFRONT_CACHE_SELF_TEST
// Validation-only constant-integrand check of the actual polar partition.
int PartitionSelfTest() {
 const double uv[8][2]={{-1,-1},{1,-1},{1,1},{-1,1},{0,-1},{1,0},{0,1},{-1,0}};
 int failures=0;const Rule quadrature=rule(9);
 for(int scenario=0;scenario<6;++scenario){
  Point points[8];long ids[8]={0,1,2,3,4,5,6,7};
  for(int i=0;i<8;++i){
   points[i].pt[0]=0;points[i].pt[1]=.03*uv[i][0];points[i].pt[2]=.03*uv[i][1];
   if(scenario==5)points[i].pt[1]*=1+.2*uv[i][1];
  }
  DSquareElement element;element.GeoInit(points,ids);
  Vec source{{-.01,0,0}};
  if(scenario==1)source[1]=.06;
  if(scenario==2)source[1]=.03;
  if(scenario==3){source[1]=.03;source[2]=.03;}
  if(scenario==4){source[1]=-.04;source[2]=.01;}
  DynaMat mat=BuildDynaMat(0.,10000.,8.,scenario==0||scenario==5?.001:.002);
  Options options;Result result;
  Context ctx(element,Kernel::U,0,mat,options,*element.m_wavefrontGeometry,source,result);
  auto cells=planeCells(ctx);long double area=0;
  for(const auto&cell:cells)for(int part=0;part<64;++part){
   double a=cell.a+(cell.b-cell.a)*part/64,b=cell.a+(cell.b-cell.a)*(part+1)/64;
   for(std::size_t j=0;j<quadrature.x.size();++j){
    double theta=(a+b)/2+(b-a)*quadrature.x[j]/2,lo,hi;
    if(!ctx.ray(theta,lo,hi))continue;
    lo=(std::max)(lo,cell.radialLow);hi=(std::min)(hi,cell.radialHigh);
    if(hi>lo)area+=(long double)quadrature.w[j]*(b-a)*(hi*hi-lo*lo)/4;
   }
  }
  if(std::abs((double)area-.0036)>1e-10*.0036)++failures;
 }
 return failures;
}
#endif
int Classify(DSquareElement&e,const Vec&source,long lag,bool second,const DynaMat&mat){
 if(lag<0||!(mat.Dt>0)||!(mat.C1>=mat.C2)||!(mat.C2>0)||!std::isfinite(mat.C1))
  return 1; // Let the integrating entry report the complete block context.
 if(second&&lag==0)return 0;
 auto geometry=e.m_wavefrontGeometry?e.m_wavefrontGeometry:BuildGeometry(e);
 // Invalid geometry or source must reach Integrate, which returns a failed
 // result with domain, element, kernel, lag and source instead of being culled.
 if(!geometry->valid)return 1;
 for(double x:source)if(!std::isfinite(x))return 1;
 Cell cell;cell.a=cell.c=-1;cell.b=cell.d=1;
 auto b=bounds(*geometry,source,cell);double low=(second?lag-1:lag)*mat.Dt;
 return b.high<mat.C2*low||b.low>=mat.C1*(low+mat.Dt)?0:1;
}
Result Integrate(DSquareElement&e,Kernel kind,const Vec&source,long lag,const DynaMat&material,const Options&opt){
 Result result;
 try{
  if((kind!=Kernel::U&&kind!=Kernel::T1&&kind!=Kernel::T2)||lag<0||!(material.Dt>0)||!(material.C2>0)||!(material.C1>=material.C2)||
    !std::isfinite(material.C1)||!(material.Rou>0)||!std::isfinite(material.Rou)||!(material.G>0)||!std::isfinite(material.G)||!(material.v> -1&&material.v<.5)||!std::isfinite(material.Dt))
   throw std::runtime_error("invalid material or lag");
  for(double x:source)if(!std::isfinite(x))throw std::runtime_error("non-finite source");
  if(!(opt.relativeTolerance>0)||!std::isfinite(opt.relativeTolerance)||!(opt.scaledAbsoluteTolerance>=0)||!std::isfinite(opt.scaledAbsoluteTolerance)||!opt.maxEvaluations)
   throw std::runtime_error("invalid quadrature options");
  if(kind==Kernel::T2&&lag==0){result.status=Status::Outside;return result;}
  // Geometry is normally constructed once by GeoInit, before worker threads.
  // A locally constructed cache supports standalone element tests as well.
  auto geometry=e.m_wavefrontGeometry?e.m_wavefrontGeometry:BuildGeometry(e);
  if(!geometry->valid)throw std::runtime_error(geometry->reason);
  MatScope scope(material);
  Context ctx(e,kind,lag,material,opt,*geometry,source,result);
  Cell whole;whole.a=whole.c=-1;whole.b=whole.d=1;
  Bounds all=bounds(*geometry,source,whole);
  if(ctx.outside(all)){result.status=Status::Outside;return result;}
  std::vector<Cell> seeds;
  bool anyFront=false;
  for(double front:ctx.fronts)if(all.low<front&&front<all.high)anyFront=true;
  if(geometry->planar&&anyFront)seeds=planeCells(ctx);
  else{
   ctx.planar=false;
   unsigned count=geometry->planar?1:NUMPW;
   for(unsigned i=0;i<count;++i)for(unsigned j=0;j<count;++j){
    Cell c;c.a=-1.+2.*i/count;c.b=-1.+2.*(i+1)/count;
    c.c=-1.+2.*j/count;c.d=-1.+2.*(j+1)/count;seeds.push_back(c);
   }
  }
  if(!ctx.planar&&anyFront)ctx.critical=criticalParameters(*geometry,source,ctx.fronts);
  result.statistics.planar=ctx.planar;
  std::priority_queue<Cell,std::vector<Cell>,Worse> heap;
  Block total{};double error=0,scale=0,missed=0;std::size_t forced=0,id=0,frozen=0;
  auto accumulate=[&](const Cell&c,double sign){
   for(int k=0;k<9;++k)total[k]+=sign*c.value[k];
   error+=sign*c.error;scale+=sign*c.scale;missed+=sign*c.missed;
  };
  auto insert=[&](Cell c){
   c.id=id++;ctx.evaluate(c);
   if(c.error==0&&c.scale==0&&!c.forced)return;
   accumulate(c,1);if(c.forced)++forced;heap.push(c);
  };
  for(auto c:seeds)insert(c);
  while(true){
   error=(std::max)(0.,error);scale=(std::max)(0.,scale);
   result.error=error;result.scale=scale;result.statistics.missedSupportError=missed;
   double target=opt.relativeTolerance*norm(total)+opt.scaledAbsoluteTolerance*scale;
   if(forced==0&&error<=target)break;
   if(heap.empty())throw std::runtime_error("maximum adaptive depth exhausted before global convergence");
   Cell c=heap.top();heap.pop();
   if(c.depth>=opt.maxDepth){++frozen;continue;}
   accumulate(c,-1);if(c.forced)--forced;
   double midx=(c.a+c.b)/2,midy=(c.c+c.d)/2;
   if(midx==c.a||midx==c.b||midy==c.c||midy==c.d)throw std::runtime_error("parameter resolution exhausted");
   for(int i=0;i<2;++i)for(int j=0;j<2;++j){
    Cell child;child.a=i?midx:c.a;child.b=i?c.b:midx;
    child.c=j?midy:c.c;child.d=j?c.d:midy;child.depth=c.depth+1;
    child.radialLow=c.radialLow;child.radialHigh=c.radialHigh;insert(child);
   }
  }
  if(kind==Kernel::U){total[3]=total[1];total[6]=total[2];total[7]=total[5];}
  result.block=total;result.statistics.leaves=heap.size()+frozen;result.status=Status::Converged;
 }catch(const std::exception&ex){
  std::ostringstream text;
  text<<"wavefront integration "<<Name(kind)<<" failed: "<<ex.what()
   <<"; domain="<<e.DomainID<<" element="<<e.LocalEleID
   <<" source=("<<source[0]<<","<<source[1]<<","<<source[2]<<") lag="<<lag
   <<" dt="<<material.Dt<<" error="<<result.error<<" scale="<<result.scale
   <<" missedSupportBound="<<result.statistics.missedSupportError<<" evaluations="<<result.statistics.evaluations<<" depth="<<result.statistics.maxDepth;
  result.status=Status::Failed;result.message=text.str();result.block.fill(0);
 }
 return result;
}
}
