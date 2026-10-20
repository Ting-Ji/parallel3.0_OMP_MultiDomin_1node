"""Read actual assembled matrices and independently reconstruct the marching operator."""
from pathlib import Path
import csv,json,struct,os
import numpy as np
from scipy.linalg import eig
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parent
def matrix(path,sym=False):
    data=path if isinstance(path,bytes) else Path(path).read_bytes();n,nz=struct.unpack_from('<ii',data)
    row=np.frombuffer(data,dtype='<i4',count=n+1,offset=8)
    cols=np.frombuffer(data,dtype='<i4',count=nz,offset=8+4*(n+1))
    width=6 if sym else 9
    values=np.frombuffer(data,dtype='<f8',count=nz*width,offset=8+4*(n+1+nz)).reshape(nz,width)
    out=np.zeros((3*n,3*n))
    for r in range(n):
        for k in range(row[r],row[r+1]):
            v=values[k]
            block=np.array([[v[0],v[1],v[2]],[v[1],v[3],v[4]],[v[2],v[4],v[5]]]) if sym else v.reshape(3,3)
            out[3*r:3*r+3,3*cols[k]:3*cols[k]+3]=block
    return out
def analyze(record):
    out=Path(record['output']);A=matrix(out/'audit_A.bin');K=matrix(out/'audit_K.bin');dim=A.shape[0]
    geometry=list(csv.DictReader((out/'audit_geometry.csv').open()))
    rotate=np.zeros((dim,dim));scale=np.ones(dim);bcids=[]
    for row in geometry:
        node=int(row['node']);bc=int(row['bcid']);bcids.append(bc)
        rotate[3*node:3*node+3,3*node:3*node+3]=np.array([float(row['t'+str(k)]) for k in range(9)]).reshape(3,3).T
        if bc==456:scale[3*node:3*node+3]=20000.
    steps=record['steps'];X=np.zeros((dim,steps));rhs=np.zeros_like(X)
    for row in csv.DictReader((out/'rhs_breakdown.csv').open()):
        i=3*int(row['rowElement'])+'xyz'.index(row['component']);n=int(row['step'])-1
        X[i,n]=float(row['solution']);rhs[i,n]=float(row['totalRhs'])
    bd=(Path(record['case'])/'input/rod_validation.bd').read_text().splitlines()
    y=np.array([float(x) for line in bd[2*len(geometry):3*len(geometry)] for x in line.split()])
    q=rotate@y
    traction_mask=np.repeat(np.array(bcids)==123,3)
    Q={};forcing={}
    for item in csv.DictReader((out/'audit_history.csv').open()):
        lag=int(item['lag'])
        T=matrix(out/f'audit_T_{lag:04d}.bin');G=matrix(out/f'audit_G_{lag:04d}.bin',bool(int(item['symmetric'])))
        if not np.any(T) and not np.any(G):continue
        unknown=np.where(traction_mask[np.newaxis,:],G,-T)@rotate
        Q[lag]=unknown
        forcing[lag]=G@q
    reconstructed=np.repeat((K@y)[:,None],steps,axis=1)
    for lag,operator in Q.items():
        for delay,coefficient in ((lag-1,1),(lag,2),(lag+1,1)):
            if delay==0 or delay>=steps:continue
            reconstructed[:,delay:]+=coefficient*(operator@X[:,:steps-delay]+forcing[lag][:,None])
    delta=reconstructed-rhs
    rhs_error=np.max(np.abs(delta))/max(np.max(np.abs(rhs)),1e-300)
    # Delay coefficients: C1=2Q1+Q2; Cj=Q(j-1)+2Qj+Q(j+1), j>=2.
    zero=np.zeros_like(A);last=max(Q)+1
    C=[2*Q.get(1,zero)+Q.get(2,zero)]
    for lag in range(2,last+1):C.append(Q.get(lag-1,zero)+2*Q.get(lag,zero)+Q.get(lag+1,zero))
    As=A/scale[np.newaxis,:];Cs=[c/scale[np.newaxis,:] for c in C]
    normalized=scale[:,None]*X
    start=max(0,steps-256)
    left=normalized[:,start:-1];right=normalized[:,start+1:]
    U,S,Vh=np.linalg.svd(left,full_matrices=False)
    rank=min(20,int(sum(S>S[0]*1e-10)))
    Ur=U[:,:rank];Vr=Vh[:rank].T
    dmd=Ur.T@right@Vr/S[:rank]
    lambdas,W=eig(dmd)
    modes=(right@Vr/S[:rank])@W
    modes_result=[]
    for idx in np.argsort(-np.abs(lambdas))[:8]:
        lam=lambdas[idx];mode=modes[:,idx]
        terms=[c@mode/(lam**(j+1)) for j,c in enumerate(Cs)]
        base=As@mode;defect=base-sum(terms)
        relative=float(np.linalg.norm(defect)/(np.linalg.norm(base)+sum(np.linalg.norm(t) for t in terms)))
        modes_result.append(dict(real=float(lam.real),imag=float(lam.imag),modulus=float(abs(lam)),
          recurrence_mode_defect=relative))
    idx=int(np.argmax(np.abs(lambdas)));lam=lambdas[idx];mode=modes[:,idx].copy()
    pivot=int(np.argmax(np.abs(mode)));mode/=mode[pivot]
    converged=False
    for iteration in range(12):
        P=As.astype(complex);derivative=np.zeros_like(P)
        for j,c in enumerate(Cs,1):
            P-=c*lam**(-j);derivative+=j*c*lam**(-j-1)
        residual=P@mode
        bordered=np.zeros((dim+1,dim+1),complex)
        bordered[:dim,:dim]=P;bordered[:dim,dim]=derivative@mode;bordered[dim,pivot]=1
        correction=np.linalg.solve(bordered,np.r_[-residual,0j])
        mode+=correction[:dim];lam+=correction[dim]
        if abs(correction[dim])<1e-13:
            converged=True;break
    base=As@mode;terms=[c@mode/lam**j for j,c in enumerate(Cs,1)]
    defect=float(np.linalg.norm(base-sum(terms))/(np.linalg.norm(base)+sum(np.linalg.norm(t) for t in terms)))
    refined=dict(real=float(lam.real),imag=float(lam.imag),modulus=float(abs(lam)),
                 converged=converged,iterations=iteration+1,recurrence_mode_defect=defect)
    result=dict(refined_growth_pole=refined,name=record['name'],last_nonzero_kernel_lag=max(Q),delay_order=last,dmd_rank=rank,
      independent_history_rhs_relative_error=float(rhs_error),growth_modes=modes_result)
    print(json.dumps(result),flush=True)
    return result
