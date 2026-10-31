#ifndef DBEM_GPU_GMRES_BRIDGE_H
#define DBEM_GPU_GMRES_BRIDGE_H
// Included after pre_gmres.h. The default CPU build has no CUDA link dependency.
#include "gpu_gmres_api.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#endif

namespace DbemGpu {
inline bool Requested() {
    std::string backend;
#ifdef _WIN32
    char* value=nullptr; size_t length=0;
    if(_dupenv_s(&value,&length,"DBEM_GMRES_BACKEND")!=0)
        throw std::runtime_error("cannot read DBEM_GMRES_BACKEND");
    std::unique_ptr<char,decltype(&std::free)> owned(value,&std::free);
    if(owned) backend=owned.get();
#else
    const char* value=std::getenv("DBEM_GMRES_BACKEND");
    if(value) backend=value;
#endif
    if(backend.empty() || backend=="cpu") return false;
    if(backend=="cuda") return true;
    throw std::runtime_error("DBEM_GMRES_BACKEND must be cpu or cuda");
}
struct Csr {
    int n=0; std::vector<int> row,col; std::vector<double> val;
    void append(std::vector<std::pair<int,double> >& entries) {
        std::stable_sort(entries.begin(),entries.end(),[](const std::pair<int,double>& a,const std::pair<int,double>& b){return a.first<b.first;});
        for(size_t i=0;i<entries.size();) {
            int c=entries[i].first; double v=0;
            do {v+=entries[i++].second;} while(i<entries.size() && entries[i].first==c);
            if(c<0 || c>=n || !std::isfinite(v)) throw std::runtime_error("invalid GPU CSR entry");
            if(v!=0) {if(val.size()>=INT_MAX) throw std::runtime_error("GPU CSR exceeds 32-bit indexing"); col.push_back(c);val.push_back(v);}
        }
        row.push_back(int(val.size()));
    }
};
inline Csr PackMatrix(const CCSRMat& a) {
    if(a.n<=0 || a.n>INT_MAX/3 || !a.row_ptr || !a.col_index || !a.value || a.row_ptr[0]!=0)
        throw std::runtime_error("invalid CCSR GPU input");
    Csr out; out.n=int(a.n*3); out.row.push_back(0);
    for(int i=0;i<a.n;++i) {
        if(a.row_ptr[i]<0 || a.row_ptr[i+1]<a.row_ptr[i]) throw std::runtime_error("invalid CCSR row range");
        for(int c=0;c<3;++c) {
            std::vector<std::pair<int,double> > entries;
            for(int k=a.row_ptr[i];k<a.row_ptr[i+1];++k) {
                if(a.col_index[k]<0 || a.col_index[k]>=a.n) throw std::runtime_error("CCSR is not square");
                for(int j=0;j<3;++j) entries.push_back({3*a.col_index[k]+j,a.value[size_t(k)*9+c*3+j]});
            }
            out.append(entries);
        }
    }
    return out;
}
inline Csr PackPreconditioner(const PreConditioner& p,int n) {
    if(p.m_LeafCount<0 || !p.m_LeafPointCount || !p.m_LeafBeginID || !p.m_RePID || !p.IID || !p.m_PreM)
        throw std::runtime_error("invalid GPU preconditioner");
    Csr out; out.n=n;
    std::vector<std::vector<std::pair<int,double> > > rows(n);
    for(long leaf=0;leaf<p.m_LeafCount;++leaf) {
        long count=p.m_LeafPointCount[leaf],begin=p.m_LeafBeginID[leaf];
        if(count<0 || count>n/3 || begin<0 || begin>n/3-count) throw std::runtime_error("invalid preconditioner leaf range");
        for(long j=0;j<count;++j) {
            long outID=p.m_OutputPID?p.m_OutputPID[begin+j]:p.m_RePID[begin+j];
            if(outID<0 || outID>=n/3) throw std::runtime_error("invalid preconditioner output map");
            long target=p.IID[outID],localRow=p.IID[j];
            if(target<0 || target>n-3 || localRow<0 || localRow>p.m_PreM[leaf].m-3) throw std::runtime_error("invalid preconditioner output offset");
            for(long k=0;k<count;++k) {
                long inID=p.m_InputPID?p.m_InputPID[begin+k]:p.m_RePID[begin+k];
                if(inID<0 || inID>=n/3) throw std::runtime_error("invalid preconditioner input map");
                long source=p.IID[inID],localCol=p.IID[k];
                if(source<0 || source>n-3 || localCol<0 || localCol>p.m_PreM[leaf].n-3) throw std::runtime_error("invalid preconditioner input offset");
                for(int r=0;r<3;++r) for(int c=0;c<3;++c)
                    rows[size_t(target+r)].push_back({int(source+c),p.m_PreM[leaf].a[localRow+r][localCol+c]});
            }
        }
    }
    out.row.push_back(0); for(auto& entries:rows) out.append(entries);
    return out;
}

class Session {
    CCSRMat* source; const PreConditioner* pre;
    Csr a,p;
    void* handle=nullptr;
    typedef int (*Create)(int,int64_t,const int*,const int*,const double*,int64_t,const int*,const int*,const double*,void**,char*,int);
    typedef int (*SolveFn)(void*,const double*,const double*,int,int,double,double*,DbemGpuStats*,char*,int);
    typedef void (*Destroy)(void*);
    SolveFn solveFn=nullptr; Destroy destroyFn=nullptr;
#ifdef _WIN32
    HMODULE library=nullptr;
#endif
public:
    Session(CCSRMat& matrix,PreConditioner& conditioner):source(&matrix),pre(&conditioner) {
        auto start=std::chrono::steady_clock::now();
        a=PackMatrix(matrix); p=PackPreconditioner(conditioner,a.n);
#ifdef _WIN32
        std::wstring path;
        wchar_t* rawCustom=nullptr; size_t customLength=0;
        if(_wdupenv_s(&rawCustom,&customLength,L"DBEM_CUDA_GMRES_DLL")!=0)
            throw std::runtime_error("cannot read DBEM_CUDA_GMRES_DLL");
        std::unique_ptr<wchar_t,decltype(&std::free)> customOwned(rawCustom,&std::free);
        const wchar_t* custom=customOwned.get();
        if(custom && *custom) {
            path=custom;
        } else {
            std::vector<wchar_t> buf(32768); DWORD len=GetModuleFileNameW(nullptr,buf.data(),DWORD(buf.size()));
            if(!len || len>=buf.size()) throw std::runtime_error("cannot locate GPU executable directory");
            path.assign(buf.data(),len); path=path.substr(0,path.find_last_of(L"\\/")+1)+L"dbem_cuda_gmres.dll";
        }
        std::vector<wchar_t> absolute(32768);
        DWORD full=GetFullPathNameW(path.c_str(),DWORD(absolute.size()),absolute.data(),nullptr);
        if(!full || full>=absolute.size()) throw std::runtime_error("invalid CUDA GMRES DLL path");
        path.assign(absolute.data(),full);
        library=LoadLibraryExW(path.c_str(),nullptr,LOAD_WITH_ALTERED_SEARCH_PATH);
        if(!library) throw std::runtime_error("cannot load dbem_cuda_gmres.dll / CUDA dependencies; Windows error "+std::to_string(GetLastError()));
        auto version=reinterpret_cast<int(*)()>(GetProcAddress(library,"dbem_gpu_abi_version"));
        auto create=reinterpret_cast<Create>(GetProcAddress(library,"dbem_gpu_create"));
        solveFn=reinterpret_cast<SolveFn>(GetProcAddress(library,"dbem_gpu_solve"));
        destroyFn=reinterpret_cast<Destroy>(GetProcAddress(library,"dbem_gpu_destroy"));
        char error[1024]={};
        if(!version || version()!=1 || !create || !solveFn || !destroyFn) {
            FreeLibrary(library); library=nullptr; throw std::runtime_error("incompatible CUDA GMRES plugin ABI");
        }
        if(create(a.n,a.val.size(),a.row.data(),a.col.data(),a.val.data(),p.val.size(),p.row.data(),p.col.data(),p.val.data(),&handle,error,sizeof(error))) {
            FreeLibrary(library); library=nullptr; throw std::runtime_error(std::string("CUDA GMRES setup: ")+error);
        }
#else
        throw std::runtime_error("this optional CUDA plugin currently supports Windows builds");
#endif
        printf("CUDA GMRES setup: n=%d nnz=%zu pre_nnz=%zu seconds=%.6f; FP64, original preconditioner\n",a.n,a.val.size(),p.val.size(),std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
    }
    Session(const Session&)=delete;
    Session& operator=(const Session&)=delete;
    ~Session() {
        if(handle && destroyFn) destroyFn(handle);
#ifdef _WIN32
        if(library) FreeLibrary(library);
#endif
    }
    bool matches(CCSRMat& matrix,PreConditioner& conditioner) const {return source==&matrix && pre==&conditioner;}
    int solve(Wvector& rhs,Wvector& initial,Wvector& out,long restart,long cycles,double tol,int& flag,int* iter) {
        if(rhs.n!=a.n || initial.n!=a.n || out.n!=a.n || restart<=0 || restart>INT_MAX || cycles<=0 || cycles>INT_MAX)
            throw std::runtime_error("CUDA GMRES vector/restart mismatch");
        char error[1024]={}; DbemGpuStats stats={};
        if(solveFn(handle,rhs.b,initial.b,int(restart),int(cycles),tol,out.b,&stats,error,sizeof(error)))
            throw std::runtime_error(std::string("CUDA GMRES solve: ")+error);
        // Verify the original, unpreconditioned residual on the host as well.
        Wvector ax(a.n); SparseMul(*source,out,ax);
        long double residual2=0,rhs2=0;
        for(int row=0;row<a.n;++row) {
            long double d=rhs.b[row]-ax.b[row]; residual2+=d*d; rhs2+=(long double)rhs.b[row]*rhs.b[row];
        }
        double actual=std::sqrt(double(rhs2>0?residual2/rhs2:residual2));
        if(!std::isfinite(actual)) throw std::runtime_error("non-finite CUDA GMRES host residual");
        flag=stats.flag; if(flag==0 && actual>tol) flag=1;
        iter[0]=stats.outer;iter[1]=stats.inner;
        printf("CUDA GMRES: flag=%d outer=%d inner=%d residual_gpu=%.9e residual_host=%.9e solve_seconds=%.6f device_MiB=%.3f\n",flag,iter[0],iter[1],stats.relative_residual,actual,stats.seconds,double(stats.device_bytes)/(1024*1024));
        return flag==3?-1:1;
    }
};
inline Session*& Current() {static thread_local Session* value=nullptr; return value;}
class Scope {
    Session* previous; std::unique_ptr<Session> owned;
public:
    Scope(CCSRMat& a,PreConditioner& p):previous(Current()) {
        if(Requested()) {owned.reset(new Session(a,p)); Current()=owned.get();}
    }
    ~Scope() {Current()=previous;}
};
inline int Solve(CCSRMat& a,PreConditioner& p,Wvector& b,Wvector& x0,Wvector& x,long restart,double tol,long cycles,int& flag,int* iter) {
    if(Current()) {
        if(!Current()->matches(a,p)) throw std::runtime_error("CUDA GMRES scope/matrix mismatch");
        return Current()->solve(b,x0,x,restart,cycles,tol,flag,iter);
    }
    // Legacy CCSR entry points get an ephemeral context; no stale pointer cache.
    Session temporary(a,p); return temporary.solve(b,x0,x,restart,cycles,tol,flag,iter);
}
}
#endif
