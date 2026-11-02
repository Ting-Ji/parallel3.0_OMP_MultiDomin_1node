#ifndef DBEM_GMRES_CPU_H
#define DBEM_GMRES_CPU_H

#include <vector>

// Owned by one solve session. Never shared by concurrent solver invocations.
struct GmresCpuWorkspace {
    long n = 0, restart = 0;
    std::vector<double> V, W, r, u, ax, partial, h, q, f;
    std::vector<long> rowCuts;
    Wmatrix QT, R;
    bool uniqueOutput = false;
    int actualThreads = 0;
    double reduced = 0.0;
    double seconds[7] = {};
    void Prepare(long size, long basis);
    GmresCpuWorkspace() = default;
    GmresCpuWorkspace(const GmresCpuWorkspace&) = delete;
    GmresCpuWorkspace& operator=(const GmresCpuWorkspace&) = delete;
};

int gmres(DSquareElement*, CCSRMat&, Wvector&, PreConditioner&, long, double,
    long, Wvector&, Wvector&, int&, int*, GmresCpuWorkspace&, long step = -1);
#ifdef DBEM_GMRES_CPU_TESTS
void GmresCpuTestOperators(CCSRMat&, PreConditioner&, Wvector&, Wvector&, Wvector&, double&, double&);
#endif

#endif
