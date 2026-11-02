// Implementation included only by pre_gmres.cpp, after the legacy CCSR solver.
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

void GmresCpuWorkspace::Prepare(long size, long basis)
{
    if (n == size && restart == basis) return;
    const size_t count = static_cast<size_t>(size);
    if (count > (std::numeric_limits<size_t>::max)() / sizeof(double) / (basis + 1))
        throw std::length_error("GMRES workspace size overflow");
    V.resize(count * (basis + 1)); W.resize(count * basis);
    r.resize(count); u.resize(count); ax.resize(count);
    partial.resize((count + 1023) / 1024);
    h.resize(basis + 1); q.resize(basis); f.resize(basis);
    QT.finish(); R.finish();
    QT.initial(basis + 1, basis + 1); R.initial(basis, basis);
    n = size; restart = basis;
}

namespace DbemCpu {
enum Stage { Spmv, Preconditioner, Orthogonalization, Norm, VectorUpdate, SmallSystem, Setup };
static const char* stageNames[] = { "spmv", "preconditioner", "orthogonalization",
    "norm", "vector_update", "small_system", "workspace_setup" };

// All methods below are collective: every member of the existing team calls them.
struct Team {
    GmresCpuWorkspace& w;
    const bool profile;
    double clock = 0.0;
    void Begin() {
        if (profile) {
#pragma omp master
            clock = DBEMWallTime();
#pragma omp barrier
        }
    }
    void End(int stage) {
        if (profile) {
#pragma omp barrier
#pragma omp master
            w.seconds[stage] += DBEMWallTime() - clock;
#pragma omp barrier
        }
    }
    double Dot(const double* a, const double* b) {
        // Small long-time TSV models amplify even a changed summation order.
        // Preserve their legacy sum; large systems use fixed parallel blocks.
        if (w.n < 65536) {
#pragma omp single
            {
                w.reduced = 0.0;
                for (long i = 0; i < w.n; ++i) w.reduced += a[i] * b[i];
            }
            const double value = w.reduced;
#pragma omp barrier
            return value;
        }
        const long blocks = static_cast<long>(w.partial.size());
#pragma omp for schedule(static)
        for (long block = 0; block < blocks; ++block) {
            double sum = 0.0;
            const long end = (std::min)(w.n, (block + 1) * 1024);
            for (long i = block * 1024; i < end; ++i) sum += a[i] * b[i];
            w.partial[block] = sum;
        }
#pragma omp single
        {
            w.reduced = 0.0;
            for (long block = 0; block < blocks; ++block) w.reduced += w.partial[block];
        }
        // No thread may overwrite reduced until every thread has copied it.
        const double value = w.reduced;
#pragma omp barrier
        return value;
    }
    double Norm2(const double* a) { return std::sqrt(Dot(a, a)); }
    void Multiply(const CCSRMat& a, const double* x, double* y) {
        const int tid = omp_get_thread_num();
        for (long i = w.rowCuts[tid]; i < w.rowCuts[tid + 1]; ++i) {
            double s0 = 0.0, s1 = 0.0, s2 = 0.0;
            for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
                const double* v = a.value + size_t(j) * 9;
                const double* z = x + size_t(a.col_index[j]) * 3;
                s0 += v[0]*z[0] + v[1]*z[1] + v[2]*z[2];
                s1 += v[3]*z[0] + v[4]*z[1] + v[5]*z[2];
                s2 += v[6]*z[0] + v[7]*z[1] + v[8]*z[2];
            }
            y[3*i] = s0; y[3*i+1] = s1; y[3*i+2] = s2;
        }
#pragma omp barrier
    }
    static void Leaf(PreConditioner& p, long leaf, const double* x, double* y) {
        const long begin = p.m_LeafBeginID[leaf];
        for (long j = 0; j < p.m_LeafPointCount[leaf]; ++j) {
            const long out = p.m_OutputPID ? p.m_OutputPID[begin+j] : p.m_RePID[begin+j];
            for (long k = 0; k < p.m_LeafPointCount[leaf]; ++k) {
                const long in = p.m_InputPID ? p.m_InputPID[begin+k] : p.m_RePID[begin+k];
                for (int m = 0; m < 3; ++m)
                    for (int n = 0; n < 3; ++n)
                        y[p.IID[out]+m] += p.m_PreM[leaf](p.IID[j]+m+1, p.IID[k]+n+1) * x[p.IID[in]+n];
            }
        }
    }
    void Apply(PreConditioner& p, const double* x, double* y) {
#pragma omp for schedule(static)
        for (long i = 0; i < w.n; ++i) y[i] = 0.0;
        if (w.uniqueOutput) {
#pragma omp for schedule(static)
            for (long leaf = 0; leaf < p.m_LeafCount; ++leaf) Leaf(p, leaf, x, y);
        } else {
#pragma omp single
            for (long leaf = 0; leaf < p.m_LeafCount; ++leaf) Leaf(p, leaf, x, y);
        }
    }
};

// Re-evaluated for the owned session, not cached by object address.
static bool UniqueOutput(const PreConditioner& p, long points) {
    std::vector<unsigned char> seen(points, 0);
    bool unique = true;
    for (long leaf = 0; leaf < p.m_LeafCount; ++leaf)
        for (long j = 0; j < p.m_LeafPointCount[leaf]; ++j) {
            const long at = p.m_LeafBeginID[leaf] + j;
            const long out = p.m_OutputPID ? p.m_OutputPID[at] : p.m_RePID[at];
            const long in = p.m_InputPID ? p.m_InputPID[at] : p.m_RePID[at];
            if (out < 0 || out >= points || in < 0 || in >= points)
                throw std::invalid_argument("GMRES preconditioner mapping outside vector");
            if (seen[out]) unique = false;
            seen[out] = 1;
        }
    return unique;
}

static int Solve(CCSRMat& a, Wvector& b, PreConditioner& p, long restart,
    double tol, long maxit, Wvector& x0, Wvector& x, int& flag, int* iter,
    GmresCpuWorkspace& w, long step)
{
    flag = 3; iter[0] = iter[1] = 0;
    const long n = 3 * a.n;
    if (n <= 0 || restart <= 0 || maxit <= 0 || !std::isfinite(tol) || tol <= 0.0 ||
        b.n != n || x.n != n || x0.n != n) return -1;
    const bool profile = []() { const char* v = std::getenv("DBEM_GMRES_CPU_PROFILE"); return v && !strcmp(v, "1"); }();
    std::fill(w.seconds, w.seconds + 7, 0.0);
    const double setupClock = DBEMWallTime();
    try { w.Prepare(n, restart); w.uniqueOutput = UniqueOutput(p, a.n); }
    catch (const std::exception& e) { fprintf(stderr, "GMRES CPU setup: %s\n", e.what()); return -1; }
    w.seconds[Setup] = DBEMWallTime() - setupClock;
    const int requested = thread_num > 0 ? static_cast<int>(thread_num) : 1;
    w.rowCuts.resize(requested + 1);
    double n2b = 0.0, tolb = 0.0, normr = 0.0, phibar = 0.0;
    bool done = false, breakdown = false;
    Team team{w, profile};
    flag = 1;
#pragma omp parallel num_threads(requested) shared(done, breakdown, n2b, tolb, normr, phibar)
    {
#pragma omp single
        {
            w.actualThreads = omp_get_num_threads();
            const long long total = static_cast<long long>(a.row_ptr[a.n]) + a.n;
            w.rowCuts[0] = 0; w.rowCuts[w.actualThreads] = a.n;
            long row = 0;
            for (int t = 1; t < w.actualThreads; ++t) {
                const long long target = total * t / w.actualThreads;
                while (row < a.n && static_cast<long long>(a.row_ptr[row]) + row < target) ++row;
                w.rowCuts[t] = row;
            }
        }
        team.Begin();
#pragma omp for schedule(static)
        for (long k = 0; k < n; ++k) x.b[k] = x0.b[k];
        team.End(VectorUpdate);
        team.Begin();
        const double bn = team.Norm2(b.b);
#pragma omp single
        { n2b = bn; tolb = tol * n2b; if (!std::isfinite(bn)) { flag = 3; done = true; } }
        team.End(Norm);
        if (!done && n2b == 0.0) {
#pragma omp for schedule(static)
            for (long k = 0; k < n; ++k) x.b[k] = 0.0;
#pragma omp single
            { flag = 0; done = true; }
        }
        if (!done) {
            team.Begin(); team.Multiply(a, x.b, w.ax.data()); team.End(Spmv);
            team.Begin();
#pragma omp for schedule(static)
            for (long k = 0; k < n; ++k) w.r[k] = b.b[k] - w.ax[k];
            team.End(VectorUpdate);
            team.Begin(); const double rn = team.Norm2(w.r.data());
#pragma omp single
            { normr = rn; if (!std::isfinite(rn)) { flag = 3; done = true; }
              else if (rn <= tolb) { flag = 0; done = true; } }
            team.End(Norm);
        }
        for (long outer = 1; outer <= maxit && !done; ++outer) {
            team.Begin();
#pragma omp single
            { w.QT = 0.0; w.R = 0.0; }
            team.End(SmallSystem);
            team.Begin(); team.Apply(p, w.r.data(), w.u.data()); team.End(Preconditioner);
            team.Begin(); const double hn = team.Norm2(w.u.data());
#pragma omp single
            { w.h[0] = hn; phibar = hn; w.QT(1,1) = 1.0;
              if (!std::isfinite(hn) || hn == 0.0) { flag = 3; done = true; } }
            team.End(Norm);
            if (done) break;
            team.Begin();
#pragma omp for schedule(static)
            for (long k = 0; k < n; ++k) w.V[k] = w.u[k] / w.h[0];
            team.End(VectorUpdate);
            for (long j = 1; j <= restart && !done; ++j) {
                team.Begin(); team.Multiply(a, w.V.data() + size_t(j-1)*n, w.ax.data()); team.End(Spmv);
                team.Begin(); team.Apply(p, w.ax.data(), w.u.data()); team.End(Preconditioner);
                team.Begin();
                for (long k = 0; k < j; ++k) {
                    const double dot = team.Dot(w.V.data() + size_t(k)*n, w.u.data());
#pragma omp single
                    w.h[k] = dot;
#pragma omp for schedule(static)
                    for (long l = 0; l < n; ++l) w.u[l] -= w.h[k] * w.V[size_t(k)*n+l];
                }
                team.End(Orthogonalization);
                team.Begin(); const double hn2 = team.Norm2(w.u.data());
#pragma omp single
                { w.h[j] = hn2; breakdown = hn2 <= 1.0e-12;
                  if (!std::isfinite(hn2)) { flag = 3; done = true; } }
                team.End(Norm);
                if (done) break;
                team.Begin();
#pragma omp for schedule(static)
                for (long k = 0; k < n; ++k) w.V[size_t(j)*n+k] = breakdown ? 0.0 : w.u[k]/w.h[j];
                team.End(VectorUpdate);
                team.Begin();
#pragma omp single
                {
                    for (long k = 1; k <= j; ++k) {
                        w.R(k,j) = 0.0;
                        for (long l = 1; l <= j; ++l) w.R(k,j) += w.QT(k,l) * w.h[l-1];
                    }
                    const double rt = w.R(j,j), next = w.h[j];
                    double c, s;
                    if (next == 0.0) { c = 1.0; s = 0.0; }
                    else if (std::fabs(next) > std::fabs(rt)) {
                        const double t = rt/next; s = 1.0/std::sqrt(1.0+t*t); c = -t*s;
                    } else { const double t = next/rt; c = 1.0/std::sqrt(1.0+t*t); s = -t*c; }
                    w.R(j,j) = c*rt - s*next;
                    for (long k = 1; k <= j; ++k) w.q[k-1] = w.QT(j,k);
                    for (long k = 1; k <= j; ++k) { w.QT(j,k) = c*w.q[k-1]; w.QT(j+1,k) = s*w.q[k-1]; }
                    w.QT(j,j+1) = -s; w.QT(j+1,j+1) = c;
                    w.f[j-1] = c*phibar; phibar *= s;
                    iter[0] = static_cast<int>(outer); iter[1] = static_cast<int>(j);
                    if (!std::isfinite(w.R(j,j)) || w.R(j,j) == 0.0) { flag = 3; done = true; }
                    if (!done && j == restart) {
                        if (inv_mat(w.R, j) < 0) { flag = 3; done = true; }
                        else for (long k = 1; k <= j; ++k) {
                            w.q[k-1] = 0.0;
                            for (long l = 1; l <= j; ++l) w.q[k-1] += w.R(k,l)*w.f[l-1];
                        }
                    }
                }
                team.End(SmallSystem);
                if (done) break;
                team.Begin();
#pragma omp for schedule(static)
                for (long k = 0; k < n; ++k) {
                    if (j < restart) {
                        double& value = w.W[size_t(k)*restart+j-1];
                        value = w.V[size_t(j-1)*n+k]/w.R(j,j);
                        for (long l = 1; l < j; ++l) value -= w.W[size_t(k)*restart+l-1]*w.R(l,j)/w.R(j,j);
                        x.b[k] += w.f[j-1]*value;
                    } else {
                        double value = 0.0;
                        for (long l = 0; l < j; ++l) value += w.V[size_t(l)*n+k]*w.q[l];
                        x.b[k] = value + x0.b[k];
                    }
                }
                team.End(VectorUpdate);
                team.Begin(); team.Multiply(a, x.b, w.ax.data()); team.End(Spmv);
                team.Begin();
#pragma omp for schedule(static)
                for (long k = 0; k < n; ++k) w.r[k] = b.b[k] - w.ax[k];
                team.End(VectorUpdate);
                team.Begin(); const double rn = team.Norm2(w.r.data());
#pragma omp single
                {
                    normr = rn;
                    if (!std::isfinite(rn)) { flag = 3; done = true; }
                    else if (rn <= tolb) { flag = 0; done = true; }
                    else if (breakdown || (j < restart && w.f[j-1] == 0.0)) { flag = 3; done = true; }
                }
                team.End(Norm);
            }
            if (!done) {
                team.Begin();
#pragma omp for schedule(static)
                for (long k = 0; k < n; ++k) x0.b[k] = x.b[k];
                team.End(VectorUpdate);
            }
        }
    }
    if (profile) {
        for (int stage = 0; stage < 7; ++stage) {
            std::string name = std::string("gmres.cpu.") + stageNames[stage];
            DBEMWriteTiming(name.c_str(), w.seconds[stage], step);
        }
        printf("GMRES_CPU step=%ld n=%ld block_nnz=%d threads=%d outer=%d inner=%d pre_serial=%d relative_residual=%.17g\n",
            step, n, a.row_ptr[a.n], w.actualThreads, iter[0], iter[1], !w.uniqueOutput,
            n2b == 0.0 ? 0.0 : normr/n2b);
    }
    return flag == 3 ? -1 : 1;
}
} // namespace DbemCpu

#ifdef DBEM_GMRES_CPU_TESTS
void GmresCpuTestOperators(CCSRMat& a, PreConditioner& p, Wvector& x,
    Wvector& ax, Wvector& px, double& dot, double& norm)
{
    GmresCpuWorkspace w; w.Prepare(x.n, 2);
    w.uniqueOutput = DbemCpu::UniqueOutput(p, a.n);
    const int requested = thread_num > 0 ? static_cast<int>(thread_num) : 1;
    w.rowCuts.resize(requested + 1);
    DbemCpu::Team team{w, false};
#pragma omp parallel num_threads(requested)
    {
#pragma omp single
        { for (int t = 0; t <= omp_get_num_threads(); ++t) w.rowCuts[t] = static_cast<long long>(a.n)*t/omp_get_num_threads(); }
        team.Multiply(a, x.b, ax.b); team.Apply(p, x.b, px.b);
        const double d = team.Dot(x.b, x.b), n = team.Norm2(x.b);
#pragma omp single
        { dot = d; norm = n; }
    }
}
#endif

int gmres(DSquareElement* elements, CCSRMat& a, Wvector& b, PreConditioner& p,
    long restart, double tol, long maxit, Wvector& x0, Wvector& x, int& flag,
    int* iter, GmresCpuWorkspace& workspace, long step)
{
    const char* mode = std::getenv("DBEM_GMRES_CPU_MODE");
    // Opt-in until physical and performance gates have passed.
    if (!mode) mode = "legacy";
    if (strcmp(mode, "legacy") && strcmp(mode, "omp") && strcmp(mode, "auto")) {
        fprintf(stderr, "Invalid DBEM_GMRES_CPU_MODE: %s\n", mode);
        flag = 3; iter[0] = iter[1] = 0; return -1;
    }
    const bool legacy = !strcmp(mode, "legacy") || (!strcmp(mode, "auto") && 3LL*a.n < 65536);
    int result = legacy ? gmres_ccsr_legacy(elements, a, b, p, restart, tol, maxit, x0, x, flag, iter)
        : DbemCpu::Solve(a, b, p, restart, tol, maxit, x0, x, flag, iter, workspace, step);
#if DBEM_GMRES_CPU_VERIFY
    // Independent validation build only: use original serial block SpMV and norm.
    Wvector check(b.n); SparseMul(a, x, check);
    double residual2 = 0.0, b2 = 0.0;
    for (long i = 1; i <= b.n; ++i) { const double d = b[i] - check[i]; residual2 += d*d; b2 += b[i]*b[i]; }
    const double residual = b2 == 0.0 ? std::sqrt(residual2) : std::sqrt(residual2/b2);
    if (flag == 0 && (!std::isfinite(residual) || residual > tol)) { flag = 3; result = -1; }
    FILE* fp = 0;
    const std::string path = DBEMOutputPath("gmres_cpu_validation.csv");
    const bool exists = [] (const std::string& name) { FILE* f = 0; fopen_s(&f, name.c_str(), "r"); if (f) { fclose(f); return true; } return false; }(path);
    fopen_s(&fp, path.c_str(), "a");
    if (fp) { if (!exists) fprintf(fp, "step,mode,flag,outer,inner,relative_residual,tolerance\n");
        fprintf(fp, "%ld,%s,%d,%d,%d,%.17g,%.17g\n", step, legacy ? "legacy" : "omp", flag, iter[0], iter[1], residual, tol); fclose(fp); }
#endif
    return result;
}

int gmres(DSquareElement* elements, CCSRMat& a, Wvector& b, PreConditioner& p,
    long restart, double tol, long maxit, Wvector& x0, Wvector& x, int& flag, int* iter)
{
    GmresCpuWorkspace workspace;
    return gmres(elements, a, b, p, restart, tol, maxit, x0, x, flag, iter, workspace);
}
