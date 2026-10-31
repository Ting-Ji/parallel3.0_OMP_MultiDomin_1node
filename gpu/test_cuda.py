"""Numerical and failure tests of the real optional CUDA library (Windows Python)."""
import ctypes as ct
import json
import os
from pathlib import Path
import sys
import numpy as np
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]


class Stats(ct.Structure):
    _fields_=[("flag",ct.c_int),("outer",ct.c_int),("inner",ct.c_int),("iterations",ct.c_int),
              ("relative_residual",ct.c_double),("seconds",ct.c_double),("device_bytes",ct.c_uint64)]


def load():
    cuda=Path(os.environ.get("CUDA_PATH",r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6"))
    directory=os.add_dll_directory(str(cuda/"bin"))
    dll=ct.CDLL(str(ROOT/"gpu/build/dbem_cuda_gmres.dll"))
    ip=ct.POINTER(ct.c_int); dp=ct.POINTER(ct.c_double)
    dll.dbem_gpu_create.argtypes=[ct.c_int,ct.c_int64,ip,ip,dp,ct.c_int64,ip,ip,dp,ct.POINTER(ct.c_void_p),ct.c_char_p,ct.c_int]
    dll.dbem_gpu_solve.argtypes=[ct.c_void_p,dp,dp,ct.c_int,ct.c_int,ct.c_double,dp,ct.POINTER(Stats),ct.c_char_p,ct.c_int]
    dll.dbem_gpu_destroy.argtypes=[ct.c_void_p]
    dll.dbem_gpu_apply.argtypes=[ct.c_void_p,ct.c_int,dp,dp,ct.c_char_p,ct.c_int]
    return dll,directory


def arrays(a):
    a=a.tocsr();a.sort_indices()
    return np.asarray(a.indptr,dtype=np.int32),np.asarray(a.indices,dtype=np.int32),np.asarray(a.data,dtype=np.float64)


def ptr(a):
    return a.ctypes.data_as(ct.POINTER(ct.c_double if a.dtype==np.float64 else ct.c_int))


def main():
    dll,directory=load(); results=[]
    for n,restart,cycles in [(12,12,1),(300,40,4),(93,3,60),(192,40,4)]:
        a=sparse.diags([np.full(n-1,-1.1),np.full(n,4.),np.full(n-1,-.7)],[-1,0,1],format="csr")
        p=sparse.eye(n,format="csr")/4
        if n==192:
            block=np.array([[.1,.2,.05],[-.15,.07,.1],[.03,-.12,.09]])
            a=sparse.kron(sparse.diags([np.ones(63),np.ones(64),np.ones(63)],[-1,0,1]),block,format="csr")+4*sparse.eye(n,format="csr")
            p=sparse.diags(1/(4+.1*(np.arange(n)%3)),format="csr")
        ar,ac,av=arrays(a);pr,pc,pv=arrays(p)
        handle=ct.c_void_p();error=ct.create_string_buffer(1024)
        assert dll.dbem_gpu_create(n,len(av),ptr(ar),ptr(ac),ptr(av),len(pv),ptr(pr),ptr(pc),ptr(pv),ct.byref(handle),error,1024)==0,error.value
        try:
            for operation,reference in [(0,a),(1,p),(2,p@a)]:
                inp=np.sin(np.arange(n)*.013);out=np.zeros(n)
                assert dll.dbem_gpu_apply(handle,operation,ptr(inp),ptr(out),error,1024)==0,error.value
                assert np.allclose(out,reference@inp,rtol=1e-13,atol=1e-13)
            for kind in ["manufactured","new_rhs","zero_rhs","exact_initial","nonzero_initial"]:
                known=np.sin(np.arange(n)*.13)+(0.1 if kind=="new_rhs" else 0)
                b=np.asarray(a@known);initial=np.zeros(n)
                if kind=="zero_rhs":b[:]=0;known[:]=0
                if kind=="exact_initial":initial=known.copy()
                if kind=="nonzero_initial":initial=.3*known
                out=np.zeros(n);stats=Stats()
                assert dll.dbem_gpu_solve(handle,ptr(b),ptr(initial),restart,cycles,1e-10,ptr(out),ct.byref(stats),error,1024)==0,error.value
                residual=np.linalg.norm(b-a@out)/max(np.linalg.norm(b),1)
                assert stats.flag==0 and residual<=1.01e-10,(n,kind,stats.flag,residual)
                assert np.linalg.norm(out-known)/max(np.linalg.norm(known),1)<1e-9
                results.append(dict(n=n,kind=kind,residual=float(residual),iterations=stats.iterations,seconds=stats.seconds))
            b[0]=np.nan
            assert dll.dbem_gpu_solve(handle,ptr(b),ptr(initial),restart,cycles,1e-10,ptr(out),ct.byref(stats),error,1024)!=0
            assert b"non-finite" in error.value
        finally:dll.dbem_gpu_destroy(handle)
    # A zero operator must not be reported as a converged solve.
    a=sparse.csr_matrix((6,6));p=sparse.eye(6,format="csr")
    ar,ac,av=arrays(a);pr,pc,pv=arrays(p);handle=ct.c_void_p()
    assert dll.dbem_gpu_create(6,0,ptr(ar),ptr(ac),ptr(av),6,ptr(pr),ptr(pc),ptr(pv),ct.byref(handle),error,1024)==0,error.value
    try:
        b=np.ones(6);x0=np.zeros(6);out=np.zeros(6);stats=Stats()
        assert dll.dbem_gpu_solve(handle,ptr(b),ptr(x0),6,1,1e-10,ptr(out),ct.byref(stats),error,1024)==0,error.value
        assert stats.flag!=0
    finally:dll.dbem_gpu_destroy(handle)
    # Invalid CSR and non-finite coefficients are rejected before upload.
    av=np.ones(6);ac=np.arange(6,dtype=np.int32);ar=np.arange(7,dtype=np.int32);av[0]=np.inf
    assert dll.dbem_gpu_create(6,6,ptr(ar),ptr(ac),ptr(av),6,ptr(pr),ptr(pc),ptr(pv),ct.byref(handle),error,1024)!=0
    output=ROOT/"gpu/validation";output.mkdir(exist_ok=True)
    (output/"library_tests.json").write_text(json.dumps({"passed":True,"cases":results,"nonfinite_and_singular_checks":True},indent=2))
    print(json.dumps({"passed":True,"solve_cases":len(results),"max_residual":max(r["residual"] for r in results)},indent=2))


if __name__=="__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
