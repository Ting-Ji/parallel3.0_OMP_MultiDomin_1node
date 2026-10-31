#define DBEM_GPU_BUILD
#include "../DBEM1/gpu_gmres_api.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(cudaError_t s) { if(s != cudaSuccess) throw std::runtime_error(cudaGetErrorString(s)); }
void blas(cublasStatus_t s) { if(s != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("cuBLAS status " + std::to_string(int(s))); }
void sparse(cusparseStatus_t s) { if(s != CUSPARSE_STATUS_SUCCESS) throw std::runtime_error("cuSPARSE status " + std::to_string(int(s))); }
void message(char* out, int size, const char* msg) { if(out && size > 0) std::snprintf(out, size, "%s", msg); }

template<class T> struct Buffer {
    T* p = nullptr; size_t count = 0;
    ~Buffer() { if(p) cudaFree(p); }
    void alloc(size_t n) {
        if(n == count) return;
        if(n > SIZE_MAX / sizeof(T)) throw std::runtime_error("GPU allocation size overflow");
        if(p) { check(cudaFree(p)); p = nullptr; count = 0; }
        if(n) check(cudaMalloc(reinterpret_cast<void**>(&p), n * sizeof(T)));
        count = n;
    }
    void upload(const T* in, size_t n) {
        alloc(n); if(n) check(cudaMemcpy(p, in, n*sizeof(T), cudaMemcpyHostToDevice));
    }
};

struct Matrix {
    Buffer<int> row, col; Buffer<double> value;
    cusparseSpMatDescr_t desc = nullptr;
    ~Matrix() { if(desc) cusparseDestroySpMat(desc); }
    void init(int n, int64_t nnz, const int* rp, const int* ci, const double* val) {
        if(nnz < 0 || nnz > INT32_MAX || !rp || (nnz && (!ci || !val)) || rp[0] != 0 || rp[n] != nnz)
            throw std::runtime_error("invalid CSR dimensions");
        for(int i=0;i<n;++i) if(rp[i]<0 || rp[i+1]<rp[i]) throw std::runtime_error("invalid CSR row pointer");
        for(int64_t k=0;k<nnz;++k) if(ci[k]<0 || ci[k]>=n || !std::isfinite(val[k]))
            throw std::runtime_error("invalid CSR column/value");
        row.upload(rp, size_t(n)+1); col.upload(ci, size_t(nnz)); value.upload(val, size_t(nnz));
        sparse(cusparseCreateCsr(&desc,n,n,nnz,row.p,col.p,value.p,CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO,CUDA_R_64F));
    }
};

struct Solver {
    int n = 0, m = 0;
    cublasHandle_t bh = nullptr; cusparseHandle_t sh = nullptr;
    Matrix a, p;
    Buffer<double> vectors, basis, coefficients;
    Buffer<unsigned char> scratchA,scratchP;
    cusparseDnVecDescr_t vx = nullptr, vy = nullptr;
    uint64_t matrix_bytes = 0;
    double *b=nullptr,*x=nullptr,*trial=nullptr,*r=nullptr,*z=nullptr,*w=nullptr,*av=nullptr;
    ~Solver() {
        if(vx) cusparseDestroyDnVec(vx); if(vy) cusparseDestroyDnVec(vy);
        if(bh) cublasDestroy(bh); if(sh) cusparseDestroy(sh);
    }
    void init(int size,int64_t nnz,const int* rp,const int* ci,const double* val,
        int64_t pnnz,const int* pr,const int* pc,const double* pv) {
        if(size<=0 || nnz<0 || pnnz<0) throw std::runtime_error("invalid solver size");
        n=size;
        matrix_bytes=(uint64_t(n)+1)*8 + (uint64_t(nnz)+uint64_t(pnnz))*12;
        size_t free=0,total=0; check(cudaMemGetInfo(&free,&total));
        if(matrix_bytes > free*0.8) throw std::runtime_error("matrices exceed GPU memory budget");
        blas(cublasCreate(&bh)); sparse(cusparseCreate(&sh));
        a.init(n,nnz,rp,ci,val); p.init(n,pnnz,pr,pc,pv);
    }
    void workspace(int restart) {
        restart=(std::min)(restart,n);
        if(m==restart) return;
        size_t free=0,total=0; check(cudaMemGetInfo(&free,&total));
        const uint64_t needed=(uint64_t(n)*(uint64_t(restart)+8)+restart)*sizeof(double);
        const uint64_t old=(vectors.count+basis.count+coefficients.count)*sizeof(double);
        if(needed>(free+old)*0.8) throw std::runtime_error("Krylov basis exceeds GPU memory budget; reduce restart");
        vectors.alloc(size_t(n)*7); basis.alloc(size_t(n)*(size_t(restart)+1)); coefficients.alloc(restart);
        b=vectors.p; x=b+n; trial=x+n; r=trial+n; z=r+n; w=z+n; av=w+n;
        if(!vx) sparse(cusparseCreateDnVec(&vx,n,x,CUDA_R_64F));
        if(!vy) sparse(cusparseCreateDnVec(&vy,n,r,CUDA_R_64F));
        sparse(cusparseDnVecSetValues(vx,x)); sparse(cusparseDnVecSetValues(vy,r));
        size_t sa=0,sp=0; double one=1,zero=0;
        sparse(cusparseSpMV_bufferSize(sh,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,a.desc,vx,&zero,vy,
            CUDA_R_64F,CUSPARSE_SPMV_CSR_ALG2,&sa));
        sparse(cusparseSpMV_bufferSize(sh,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,p.desc,vx,&zero,vy,
            CUDA_R_64F,CUSPARSE_SPMV_CSR_ALG2,&sp));
        // Each operator owns its scratch; do not reuse A scratch for P.
        // In particular a zero-byte P buffer must remain nullptr.
        scratchA.alloc(sa);scratchP.alloc(sp);m=restart;
    }
    void mul(Matrix& mat,double* in,double* out) {
        double one=1,zero=0;
        sparse(cusparseDnVecSetValues(vx,in)); sparse(cusparseDnVecSetValues(vy,out));
        sparse(cusparseSpMV(sh,CUSPARSE_OPERATION_NON_TRANSPOSE,&one,mat.desc,vx,&zero,vy,
            CUDA_R_64F,CUSPARSE_SPMV_CSR_ALG2,(&mat==&a?scratchA.p:scratchP.p)));
    }
    void copy(double* from,double* to) { blas(cublasDcopy(bh,n,from,1,to,1)); }
    void axpy(double alpha,double* in,double* out) { blas(cublasDaxpy(bh,n,&alpha,in,1,out,1)); }
    double norm(double* in) {
        double v=0; blas(cublasDnrm2(bh,n,in,1,&v));
        if(!std::isfinite(v)) throw std::runtime_error("non-finite GPU vector norm"); return v;
    }
    double residual(double* candidate,double bnorm) {
        mul(a,candidate,av); copy(b,r); axpy(-1,av,r); return norm(r)/bnorm;
    }
    void solve(const double* rhs,const double* initial,int restart,int cycles,double tol,double* out,DbemGpuStats& stats) {
        if(!rhs || !initial || !out || restart<=0 || cycles<=0 || !(tol>0) || !std::isfinite(tol))
            throw std::runtime_error("invalid solve arguments");
        for(int i=0;i<n;++i) if(!std::isfinite(rhs[i]) || !std::isfinite(initial[i])) throw std::runtime_error("non-finite solve input");
        workspace(restart); stats={}; stats.flag=1;
        stats.device_bytes=matrix_bytes+(vectors.count+basis.count+coefficients.count)*sizeof(double)+scratchA.count+scratchP.count;
        check(cudaMemcpy(b,rhs,n*sizeof(double),cudaMemcpyHostToDevice));
        check(cudaMemcpy(x,initial,n*sizeof(double),cudaMemcpyHostToDevice));
        double bn=norm(b);
        if(bn==0) { check(cudaMemset(x,0,n*sizeof(double))); stats.flag=0; }
        else if((stats.relative_residual=residual(x,bn))<=tol) stats.flag=0;
        if(stats.flag!=0) {
            std::vector<double> h(size_t(m+1)*m),cs(m),sn(m),g(m+1),y(m);
            for(int outer=0;outer<cycles && stats.flag!=0;++outer) {
                stats.outer=outer+1;
                mul(p,r,z); double beta=norm(z);
                if(beta==0) { stats.flag=3; break; }
                copy(z,basis.p); double inv=1/beta; blas(cublasDscal(bh,n,&inv,basis.p,1));
                std::fill(h.begin(),h.end(),0); std::fill(g.begin(),g.end(),0); g[0]=beta;
                for(int j=0;j<m;++j) {
                    ++stats.iterations; stats.inner=j+1;
                    mul(a,basis.p+size_t(j)*n,av); mul(p,av,w);
                    const double before=norm(w);
                    for(int k=0;k<=j;++k) {
                        double dot=0; blas(cublasDdot(bh,n,basis.p+size_t(k)*n,1,w,1,&dot));
                        if(!std::isfinite(dot)) throw std::runtime_error("non-finite Arnoldi dot product");
                        h[size_t(j)*(m+1)+k]=dot; axpy(-dot,basis.p+size_t(k)*n,w);
                    }
                    double hn=norm(w); h[size_t(j)*(m+1)+j+1]=hn;
                    bool breakdown=hn<=std::numeric_limits<double>::epsilon()*(std::max)(before,1e-300);
                    if(!breakdown) {
                        copy(w,basis.p+size_t(j+1)*n); inv=1/hn;
                        blas(cublasDscal(bh,n,&inv,basis.p+size_t(j+1)*n,1));
                    }
                    for(int k=0;k<j;++k) {
                        double& u=h[size_t(j)*(m+1)+k]; double& v=h[size_t(j)*(m+1)+k+1];
                        double t=cs[k]*u+sn[k]*v; v=-sn[k]*u+cs[k]*v; u=t;
                    }
                    double& d=h[size_t(j)*(m+1)+j]; double sub=h[size_t(j)*(m+1)+j+1];
                    double rho=std::hypot(d,sub);
                    if(rho==0 || !std::isfinite(rho)) { stats.flag=3; break; }
                    cs[j]=d/rho; sn[j]=sub/rho; d=rho; h[size_t(j)*(m+1)+j+1]=0;
                    g[j+1]=-sn[j]*g[j]; g[j]*=cs[j];
                    for(int k=j;k>=0;--k) {
                        double v=g[k]; for(int l=k+1;l<=j;++l) v-=h[size_t(l)*(m+1)+k]*y[l];
                        y[k]=v/h[size_t(k)*(m+1)+k];
                    }
                    check(cudaMemcpy(coefficients.p,y.data(),(j+1)*sizeof(double),cudaMemcpyHostToDevice));
                    copy(x,trial); double one=1;
                    blas(cublasDgemv(bh,CUBLAS_OP_N,n,j+1,&one,basis.p,n,coefficients.p,1,&one,trial,1));
                    stats.relative_residual=residual(trial,bn);
                    if(stats.relative_residual<=tol) { copy(trial,x); stats.flag=0; break; }
                    if(breakdown) { copy(trial,x); stats.flag=3; break; }
                    if(j==m-1) copy(trial,x);
                }
                if(stats.flag==3) break;
            }
        }
        if(bn>0) stats.relative_residual=residual(x,bn);
        check(cudaMemcpy(out,x,n*sizeof(double),cudaMemcpyDeviceToHost));
    }
};
}

extern "C" int dbem_gpu_abi_version() { return 1; }
extern "C" int dbem_gpu_create(int n,int64_t nnz,const int* r,const int* c,const double* v,
    int64_t pn,const int* pr,const int* pc,const double* pv,void** handle,char* err,int size) {
    if(!handle) { message(err,size,"null output handle"); return 1; } *handle=nullptr;
    try { std::unique_ptr<Solver> s(new Solver); s->init(n,nnz,r,c,v,pn,pr,pc,pv); *handle=s.release(); message(err,size,""); return 0; }
    catch(const std::exception& e) { message(err,size,e.what()); return 1; }
}
extern "C" int dbem_gpu_solve(void* handle,const double* rhs,const double* initial,int restart,int cycles,
    double tol,double* out,DbemGpuStats* stats,char* err,int size) {
    try {
        if(!handle || !stats) throw std::runtime_error("null solver/stats");
        auto start=std::chrono::steady_clock::now();
        static_cast<Solver*>(handle)->solve(rhs,initial,restart,cycles,tol,out,*stats);
        stats->seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        message(err,size,""); return 0;
    } catch(const std::exception& e) { message(err,size,e.what()); return 1; }
}
extern "C" void dbem_gpu_destroy(void* handle) { delete static_cast<Solver*>(handle); }
extern "C" int dbem_gpu_apply(void* handle,int preconditioner,const double* input,double* output,char* err,int size) {
    try {
        if(!handle || !input || !output || preconditioner<0 || preconditioner>2) throw std::runtime_error("invalid apply arguments");
        Solver& s=*static_cast<Solver*>(handle);s.workspace(1);
        check(cudaMemcpy(s.x,input,s.n*sizeof(double),cudaMemcpyHostToDevice));
        if(preconditioner==2) {s.mul(s.a,s.x,s.av);s.mul(s.p,s.av,s.w);check(cudaMemcpy(output,s.w,s.n*sizeof(double),cudaMemcpyDeviceToHost));message(err,size,"");return 0;}
        s.mul(preconditioner?s.p:s.a,s.x,s.av);
        check(cudaMemcpy(output,s.av,s.n*sizeof(double),cudaMemcpyDeviceToHost));
        message(err,size,"");return 0;
    } catch(const std::exception& e) {message(err,size,e.what());return 1;}
}
