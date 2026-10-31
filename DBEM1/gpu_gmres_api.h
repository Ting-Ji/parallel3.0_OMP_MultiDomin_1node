#ifndef DBEM_GPU_GMRES_API_H
#define DBEM_GPU_GMRES_API_H
#include <stdint.h>

#if defined(_WIN32) && defined(DBEM_GPU_BUILD)
#define DBEM_GPU_API __declspec(dllexport)
#else
#define DBEM_GPU_API
#endif

// Stable C ABI. Both matrices are zero-based scalar CSR, in double precision.
struct DbemGpuStats {
    int flag, outer, inner, iterations;
    double relative_residual, seconds;
    uint64_t device_bytes;
};

extern "C" {
DBEM_GPU_API int dbem_gpu_abi_version();
DBEM_GPU_API int dbem_gpu_create(int n, int64_t nnz, const int* rows,
    const int* cols, const double* values, int64_t pnnz, const int* prows,
    const int* pcols, const double* pvalues, void** handle, char* error, int error_size);
DBEM_GPU_API int dbem_gpu_solve(void* handle, const double* rhs, const double* initial,
    int restart, int cycles, double tolerance, double* solution,
    DbemGpuStats* stats, char* error, int error_size);
DBEM_GPU_API void dbem_gpu_destroy(void* handle);
// Diagnostic operations: 0 = A*x, 1 = P*x, 2 = P*(A*x).
DBEM_GPU_API int dbem_gpu_apply(void* handle, int preconditioner, const double* input,
    double* output, char* error, int error_size);
}
#endif
