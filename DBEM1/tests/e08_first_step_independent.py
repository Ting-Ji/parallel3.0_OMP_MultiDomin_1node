"""Solve the corrected first-step matrix independently of the C++ GMRES code."""
from pathlib import Path
import struct,csv,json
import numpy as np
from scipy.sparse import bsr_matrix
from scipy.sparse.linalg import splu
base=Path("validation_runs/e08_nan_diagnosis_20261019")
work=next(base.glob("e08_short_*/summary.json")).parent
case=work/"t1";cache=case/"output/gmres_ccsr_cache.bin"
with cache.open("rb") as f:
 h=f.read(144);assert struct.unpack_from("<i",h,84)[0]==3
 count=struct.unpack_from("<i",h,20)[0];entries=[struct.unpack("<iiqq",f.read(24)) for _ in range(count)]
 e=next(e for e in entries if e[0]==1 and e[1]==0);f.seek(e[2]);d=f.read(e[3])
n,nz=struct.unpack_from("<ii",d);rp=np.frombuffer(d,"<i4",n+1,8);ci=np.frombuffer(d,"<i4",nz,8+4*(n+1));v=np.frombuffer(d,"<f8",9*nz,8+4*(n+1+nz)).reshape(-1,3,3)
assert np.isfinite(v).all();A=bsr_matrix((v,ci,rp),shape=(n*3,n*3)).tocsc()
rows=[]
with next((case/"output").rglob("rhs_breakdown.csv")).open() as f:
 for row in csv.DictReader(f):
  if int(row["step"])!=1:break
  rows.append(row)
b=np.array([float(r["totalRhs"]) for r in rows]);gmres=np.array([float(r["solution"]) for r in rows]);x=splu(A).solve(b)
result=dict(matrix_dimension=A.shape[0],matrix_nonfinite=0,direct_relative_residual=float(np.linalg.norm(A@x-b)/np.linalg.norm(b)),gmres_relative_residual=float(np.linalg.norm(A@gmres-b)/np.linalg.norm(b)),gmres_vs_direct_relative_difference=float(np.linalg.norm(gmres-x)/np.linalg.norm(x)),max_abs_difference=float(np.max(np.abs(gmres-x))))
assert result["direct_relative_residual"]<1e-12 and result["gmres_relative_residual"]<1e-9 and result["gmres_vs_direct_relative_difference"]<1e-7,result
(base/"independent_first_step.json").write_text(json.dumps(result,indent=2));print(json.dumps(result))
