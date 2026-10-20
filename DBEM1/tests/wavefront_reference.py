"""Independent reference kernels and high-order wavefront-split quadrature."""
from pathlib import Path
import csv,json,math,argparse
import numpy as np
from numpy.polynomial.legendre import leggauss

def kernel(v,normal,nu,dt,kind,lag):
 mu=20000/(2*(1+nu));lam=2*mu*nu/(1-2*nu);rho=8.
 cp=np.sqrt(2*(1-nu)*mu/((1-2*nu)*rho));cs=np.sqrt(mu/rho)
 r=np.linalg.norm(v,axis=-1)
 if kind==0:lo,hi,a,b=lag*dt,(lag+1)*dt,1.,0.
 elif kind==1:lo,hi,a,b=lag*dt,(lag+1)*dt,lag+1.,-1/dt
 else:lo,hi,a,b=(lag-1)*dt,lag*dt,1.-lag,1/dt
 tp=r/cp;ts=r/cs
 mp=(tp>=lo)&(tp<hi);ms=(ts>=lo)&(ts<hi)
 pp=np.where(mp,a+b*tp,0.);ss=np.where(ms,a+b*ts,0.)
 low=np.maximum(lo,tp);high=np.minimum(hi,ts)
 F=np.where(high>low,a*(high*high-low*low)/2+b*(high**3-low**3)/3,0.)
 A=ss/(cs*cs*r)-F/r**3
 B=pp/(cp*cp*r**3)-ss/(cs*cs*r**3)+3*F/r**5
 vv=v[..., :,None]*v[...,None,:];eye=np.eye(3)
 if kind==0:return (A[...,None,None]*eye+B[...,None,None]*vv)/(4*np.pi*rho)
 Fprime=r*(ss/cs**2-pp/cp**2)
 Aprime=np.where(ms,b/cs,0.)/(cs*cs*r)-ss/(cs*cs*r*r)-Fprime/r**3+3*F/r**4
 Bprime=np.where(mp,b/cp,0.)/(cp*cp*r**3)-3*pp/(cp*cp*r**4)-np.where(ms,b/cs,0.)/(cs*cs*r**3)+3*ss/(cs*cs*r**4)+3*Fprime/r**5-15*F/r**6
 vn=np.sum(v*normal,axis=-1);nv=normal[..., :,None]*v[...,None,:];vnn=v[..., :,None]*normal[...,None,:]
 result=lam*(Aprime/r+Bprime*r+4*B)[...,None,None]*nv
 result+=mu*((vn*(Aprime/r+B))[...,None,None]*eye+(2*Bprime*vn/r)[...,None,None]*vv+
             (2*B)[...,None,None]*nv+(B+Aprime/r)[...,None,None]*vnn)
 return result/(4*np.pi*rho)

def material(nu):return np.sqrt(20000*(1-nu)/((1+nu)*(1-2*nu)*8)),np.sqrt(20000/(2*(1+nu)*8))
def fronts(nu,dt,kind,lag):
 cp,cs=material(nu);low=(lag-1 if kind==2 else lag)*dt
 return sorted(set(c*t for c in (cp,cs) for t in (low,low+dt) if t>0))

def plane_reference(points,source,nu,dt,kind,lag,q):
 vertices=points[:4];n=np.cross(vertices[1]-vertices[0],vertices[3]-vertices[0]);n/=np.linalg.norm(n)
 e1=vertices[1]-vertices[0];e1/=np.linalg.norm(e1);e2=np.cross(n,e1)
 polygon=np.stack([(vertices-source)@e1,(vertices-source)@e2],axis=1)
 gap=float((source-vertices[0])@n)
 radii=[0.]+[np.sqrt((r-abs(gap))*(r+abs(gap))) for r in fronts(nu,dt,kind,lag) if r>abs(gap)]
 if len(radii)==1:return np.zeros((3,3)),0.
 outer=max(radii);near=float('inf')
 inside=True
 for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
  edge=b-a
  inside=inside and np.cross(edge,-a)>=0
  t=np.clip(-a@edge/(edge@edge),0,1);near=min(near,np.linalg.norm(a+t*edge))
 near=abs(gap) if inside else np.hypot(gap,near)
 if (lag-1 if kind==2 else lag)==0 and near>0:
  x=near
  while x<outer:radii.append(x);x*=2
 radii=sorted(set(radii));angles=[0.,2*np.pi]
 angles+=list(np.mod(np.arctan2(polygon[:,1],polygon[:,0]),2*np.pi))
 for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
  edge=b-a
  for radius in radii:
   roots=np.roots([edge@edge,2*(a@edge),a@a-radius*radius])
   for t in roots:
    if abs(t.imag)<1e-12 and 0<=t.real<=1:
     v=a+t.real*edge;angles.append(float(np.mod(np.arctan2(v[1],v[0]),2*np.pi)))
 angles=sorted(set(angles));gx,gw=leggauss(q)
 vectors=[];weights=[]
 for low,high in zip(angles[:-1],angles[1:]):
  theta=(low+high)/2+(high-low)/2*gx
  for t,wt in zip(theta,gw*(high-low)/2):
   direction=np.array([np.cos(t),np.sin(t)])
   entry,exit=0.,float('inf')
   for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
    edge=b-a;den=float(np.cross(edge,direction));rhs=float(np.cross(edge,a))
    if abs(den)<1e-15*np.linalg.norm(edge):
     if rhs>0:exit=-1.;break
    elif den>0:entry=max(entry,rhs/den)
    else:exit=min(exit,rhs/den)
   if exit<=entry:continue
   for a,b in zip(radii[:-1],radii[1:]):
    a=max(a,entry);b=min(b,exit)
    if b<=a:continue
    rho=(a+b)/2+(b-a)/2*gx
    vectors.append(-gap*n+rho[:,None]*(e1*np.cos(t)+e2*np.sin(t)))
    weights.append(gw*(b-a)/2*rho*wt)
 if not vectors:return np.zeros((3,3)),0.
 v=np.concatenate(vectors);w=np.concatenate(weights);normal=np.broadcast_to(n,v.shape)
 values=kernel(v,normal,nu,dt,kind,lag)
 return np.einsum('kab,k->ab',values,w),float(w@np.linalg.norm(values,axis=(1,2)))

def bowl_reference(ident,source,nu,dt,kind,lag,q):
 h=.1;k=.01 if ident==8 else .03;gap=-source[0]
 fs=fronts(nu,dt,kind,lag);radii=[0.]
 for front in fs:
  # k^2 rho^4 + (h^2+2 k gap) rho^2 + gap^2 = front^2
  B=h*h+2*k*gap;C=gap*gap-front*front
  if C<0:radii.append(np.sqrt(-2*C/(B+np.sqrt(B*B-4*k*k*C))))
 angles=[i*np.pi/4 for i in range(9)]
 for r in radii:
  if 1<r<np.sqrt(2):
   t=np.arccos(1/r)
   for i in range(4):angles.extend([i*np.pi/2+t,i*np.pi/2+np.pi/2-t])
 angles=sorted(set(angles));gx,gw=leggauss(q);vectors=[];normals=[];weights=[]
 for a,b in zip(angles[:-1],angles[1:]):
  for theta,wt in zip((a+b)/2+(b-a)/2*gx,gw*(b-a)/2):
   maximum=1/max(abs(np.cos(theta)),abs(np.sin(theta)))
   cuts=sorted(set([0.,maximum]+[r for r in radii if 0<r<maximum]))
   for low,high in zip(cuts[:-1],cuts[1:]):
    rho=(low+high)/2+(high-low)/2*gx;u=rho*np.cos(theta);v=rho*np.sin(theta)
    vectors.append(np.stack([gap+k*rho*rho,h*u,h*v],axis=1))
    normal=np.stack([np.full(q,h*h),-2*k*h*u,-2*k*h*v],axis=1)
    jac=np.linalg.norm(normal,axis=1);normals.append(normal/jac[:,None])
    weights.append(gw*(high-low)/2*wt*rho*jac)
 v=np.concatenate(vectors);n=np.concatenate(normals);w=np.concatenate(weights);values=kernel(v,n,nu,dt,kind,lag)
 return np.einsum('kab,k->ab',values,w),float(w@np.linalg.norm(values,axis=(1,2)))

def curved_plane_reference(source,nu,dt,kind,lag,q):
 h=.1;k=.01;fs=fronts(nu,dt,kind,lag);outer=[-1.,1.]
 for r in fs:
  if r/h<1:outer.extend([-r/h,r/h])
  for sign in (-1,1):
   for u2 in np.roots([k*k,h*h+2*sign*h*k,h*h-r*r]):
    if abs(u2.imag)<1e-12 and 0<u2.real<1:outer.extend([-np.sqrt(u2.real),np.sqrt(u2.real)])
 outer=sorted(set(outer));gx,gw=leggauss(q);vectors=[];weights=[]
 for a,b in zip(outer[:-1],outer[1:]):
  for u,wu in zip((a+b)/2+(b-a)/2*gx,gw*(b-a)/2):
   cuts=[-1.,1.]
   for r in fs:
    D=r*r-(h*u)**2
    if D>=0:
     cuts.extend(v for v in ((-k*u*u-np.sqrt(D))/h,(-k*u*u+np.sqrt(D))/h) if -1<v<1)
   cuts=sorted(set(cuts))
   for low,high in zip(cuts[:-1],cuts[1:]):
    v=(low+high)/2+(high-low)/2*gx
    vectors.append(np.stack([np.zeros(q),np.full(q,h*u),h*v+k*u*u],axis=1)-source)
    weights.append(gw*(high-low)/2*wu*h*h)
 v=np.concatenate(vectors);w=np.concatenate(weights);n=np.zeros_like(v);n[:,0]=1
 values=kernel(v,n,nu,dt,kind,lag)
 return np.einsum('kab,k->ab',values,w),float(w@np.linalg.norm(values,axis=(1,2)))

def main():
 p=argparse.ArgumentParser();p.add_argument('--mode',default='default',choices=['default','new','legacy']);args=p.parse_args()
 root=Path(__file__).resolve().parents[2]/'validation_runs/wavefront_implementation_20261019'/args.mode/'base'
 groups={}
 for r in csv.DictReader((root/'kernel_probe.csv').open()):
  groups.setdefault((int(r['case']),int(r['n'])),[]).append(r)
 point_errors=[]
 for rs in groups.values():
  r=rs[0];v=np.array([[float(r[k]) for k in ('rx','ry','rz')]]);normal=np.array([[.3,.4,np.sqrt(.75)]])
  for kind,name in ((0,'U'),(1,'T1'),(2,'T2')):
   expected=kernel(v,normal,float(r['nu']),float(r['dt']),kind,int(r['n']))[0]
   actual=np.array([float(x[name]) for x in rs]).reshape(3,3)
   point_errors.append(float(np.linalg.norm(expected-actual)/max(np.linalg.norm(actual),1e-30)))
 assert max(point_errors)<1e-10
 (root.parent.parent/'independent_point_kernel_check.json').write_text(json.dumps(dict(max_relative_error=max(point_errors),samples=len(point_errors)),indent=2))
 geometries={}
 for r in csv.DictReader((root/'wf_geometry.csv').open()):
  geometries.setdefault(int(r['case']),[]).append(r)
 rows=list(csv.DictReader((root/'wf_blocks.csv').open()));summary=[]
 for r in rows:
  ident=int(r['case']);kind=int(r['kernel']);lag=int(r['lag']);shift=float(r['epsilon']);nu=float(r['nu']);dt=float(r['dt'])
  if shift!=0:continue
  points=np.array([[float(x[k]) for k in ('x','y','z')] for x in geometries[ident]])
  source=np.array([float(geometries[ident][0][k]) for k in ('sx','sy','sz')])
  def reference(q):
   if ident in (8,9):return bowl_reference(ident,source,nu,dt,kind,lag,q)
   if ident==10:
    value,scale=curved_plane_reference(source,nu,dt,kind,lag,q)
    # The complete support disk is inside this planar curved-edge panel.
    # Use its closed-form U integral; a fixed outer rule converges slowly at
    # the internal square-root tangencies and is not a valid reference here.
    if kind==0:
     cp,cs=material(nu);value=np.diag([dt/(16*cp),dt/(16*cs),dt/(16*cs)])
    return value,scale
   return plane_reference(points,source,nu,dt,kind,lag,q)
  coarse,_=reference(48);ref,scale=reference(64)
  reference_change=np.linalg.norm(ref-coarse)/(max(np.linalg.norm(ref),scale*1e-3,1e-300))
  actual=np.array([float(r['b'+str(k)]) for k in range(9)]).reshape(3,3)
  error=float(np.linalg.norm(actual-ref));allowed=1e-6*np.linalg.norm(ref)+1e-10*scale
  entry=dict(case=ident,kernel=kind,status=int(r['status']),error=error,allowed=float(allowed),ratio=error/max(allowed,1e-300),
   reference_change=float(reference_change),scale=scale,evaluations=int(r['evaluations']),planar=int(r['planar']))
  summary.append(entry);print(json.dumps(entry),flush=True)
 (root.parent.parent/'local_reference_comparison.json').write_text(json.dumps(summary,indent=2))
 assert all(r['status']!=2 and r['ratio']<=1 and r['reference_change']<=1e-8 for r in summary),[r for r in summary if r['ratio']>1 or r['status']==2 or r['reference_change']>1e-8]
if __name__=='__main__':main()
