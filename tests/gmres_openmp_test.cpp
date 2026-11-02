#include "pre_gmres.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <vector>

static void require(bool ok, const char* what) { if (!ok) throw std::runtime_error(what); }
static double difference(Wvector& a, Wvector& b) {
    double d = 0.0, s = 0.0;
    for (long i = 1; i <= a.n; ++i) { d += (a[i]-b[i])*(a[i]-b[i]); s += b[i]*b[i]; }
    return std::sqrt(d / (s > 1e-60 ? s : 1e-60));
}
struct Fixture {
    CCSRMat a; PreConditioner p; int blocks;
    Fixture(int count, bool empty = false, bool overlap = false) : blocks(count) {
        a.initial(count, 3*count-2); a.row_ptr[0] = 0; int at = 0;
        for (int i = 0; i < count; ++i) {
            if (!(empty && i == count/2)) for (int j = (i ? i-1 : 0); j <= i+1 && j < count; ++j) {
                a.col_index[at] = j;
                for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c)
                    a.value[size_t(at)*9+3*r+c] = i==j ? (r==c ? 4.0+r*.2 : .03) : (r==c ? -.8 : .02);
                ++at;
            }
            a.row_ptr[i+1] = at;
        }
        p.m_LeafCount = count; p.m_LeafPointCount = new long[count]; p.m_LeafBeginID = new long[count];
        p.m_RePID = new long[count]; p.IID = new long[count]; p.m_InputPID = new long[count];
        p.m_OutputPID = new long[count]; p.m_PreM = new Wmatrix[count];
        for (int i = 0; i < count; ++i) {
            p.m_LeafPointCount[i] = 1; p.m_LeafBeginID[i] = i; p.m_RePID[i] = i; p.IID[i] = 3*i;
            p.m_InputPID[i] = count-1-i; p.m_OutputPID[i] = overlap ? 0 : count-1-i;
            p.m_PreM[i].initial(3,3); p.m_PreM[i] = 0.0;
            for (int c = 1; c <= 3; ++c) p.m_PreM[i](c,c) = 1.0/(4.0+.2*(c-1));
        }
    }
};
static void operators(bool empty, bool overlap, int blocks = 53) {
    Fixture f(blocks, empty, overlap); const int n = 3*f.blocks;
    Wvector x(n), ax(n), px(n), referenceA(n), referenceP(n);
    for (int i = 1; i <= n; ++i) x[i] = std::sin(.13*i);
    SparseMul(f.a, x, referenceA); GMRES_M_X(f.p, x, referenceP);
    double d, norm; GmresCpuTestOperators(f.a, f.p, x, ax, px, d, norm);
    require(difference(ax, referenceA) <= 1e-12, "SpMV operator");
    require(difference(px, referenceP) <= 1e-12, "preconditioner operator");
    require(std::fabs(norm-vector_norm2(x,n)) <= 1e-12*norm, "norm");
    require(std::fabs(d-d_dot(x,x,n)) <= 1e-12*d, "dot");
}
static void solves() {
    Fixture f(53); const int n = 3*f.blocks;
    Wvector exact(n), b(n), x0(n,0), x(n), check(n); GmresCpuWorkspace workspace;
    for (int i = 1; i <= n; ++i) exact[i] = std::sin(.13*i);
    SparseMul(f.a, exact, b); int flag, iter[2];
    _putenv_s("DBEM_GMRES_CPU_MODE", "omp");
    // Small restart forces the formerly unsafe final-column branch and restarts.
    gmres(nullptr,f.a,b,f.p,3,1e-10,30,x0,x,flag,iter,workspace);
    require(flag == 0 && iter[0] > 1, "restarted convergence");
    require(difference(x,exact) <= 1e-8, "known solution");
    const std::vector<double> first(x.b,x.b+n);
    x0 = 0.0; gmres(nullptr,f.a,b,f.p,3,1e-10,30,x0,x,flag,iter,workspace);
    require(flag == 0, "repeated solve");
    for (int i = 0; i < n; ++i) require(x.b[i] == first[i], "deterministic repeat");
    x0 = exact; gmres(nullptr,f.a,b,f.p,3,1e-10,30,x0,x,flag,iter,workspace);
    require(flag == 0 && iter[0] == 0, "exact initial guess");
    for (int i = 1; i <= n; ++i) x0[i] = exact[i]*.3;
    gmres(nullptr,f.a,b,f.p,3,1e-10,30,x0,x,flag,iter,workspace); require(flag == 0, "nonzero initial guess");
    x0 = 0.0; gmres(nullptr,f.a,b,f.p,1,1e-14,1,x0,x,flag,iter,workspace);
    require(flag == 1 && iter[0] == 1 && iter[1] == 1, "iteration limit");
    b = 0.0; x0 = 1.0; gmres(nullptr,f.a,b,f.p,3,1e-10,3,x0,x,flag,iter,workspace);
    require(flag == 0 && vector_norm2(x,n) == 0.0, "zero RHS");
    b[1] = std::numeric_limits<double>::quiet_NaN();
    gmres(nullptr,f.a,b,f.p,3,1e-10,3,x0,x,flag,iter,workspace); require(flag != 0, "NaN rejection");
    b = 1.0; x0 = 0.0;
    for (int at = 0; at < f.a.row_ptr[f.a.n]*9; ++at) f.a.value[at] = 0.0;
    gmres(nullptr,f.a,b,f.p,3,1e-10,3,x0,x,flag,iter,workspace); require(flag != 0, "singular breakdown");
    for (int i = 0; i < f.blocks; ++i) for (int j = f.a.row_ptr[i]; j < f.a.row_ptr[i+1]; ++j)
        if (f.a.col_index[j] == i) for (int r = 0; r < 3; ++r) f.a.value[size_t(j)*9+4*r] = 4.0+.2*r;
    x0 = 0.0; gmres(nullptr,f.a,b,f.p,3,1e-10,3,x0,x,flag,iter,workspace);
    require(flag == 0 && iter[1] == 1, "happy breakdown");
    _putenv_s("DBEM_GMRES_CPU_MODE", "auto"); x0 = 0.0;
    gmres(nullptr,f.a,b,f.p,10,1e-10,3,x0,x,flag,iter); require(flag == 0, "auto small legacy");
}
int main() {
    thread_num = 40;
    try { operators(false,false); operators(true,false); operators(false,true);
        operators(false,false,21847); solves();
        printf("PASS: GMRES operators, mappings, fixed reductions, restart, repeat, initial guesses, limit, zero, nonfinite, breakdown; threads=40\n"); return 0; }
    catch (const std::exception& e) { fprintf(stderr,"FAIL: %s\n",e.what()); return 2; }
}
