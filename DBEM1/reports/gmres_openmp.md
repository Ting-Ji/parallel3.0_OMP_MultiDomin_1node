# CCSR GMRES OpenMP backend

The original CCSR implementation remains available as the `legacy` reference.
`DBEM_GMRES_CPU_MODE=omp` selects the new backend; `auto` selects it for at least
65,536 scalar unknowns. Unset mode currently means `legacy`: enabling the default
requires both numerical gates and the 1.5x performance target on both large cases.
Other matrix formats and the CUDA branch are unchanged.

The backend retains FP64, left preconditioning, sequential modified Gram–Schmidt,
incremental W solution updates and an unpreconditioned true-residual test at each
inner iteration. The final restart-column combination uses exactly the active
columns, avoiding the old V/restart+1 versus q/restart mismatch. Nonfinite norms,
zero pivots and unsuccessful breakdown return failure, never convergence.

One OpenMP team is used per solve. CCSR row ranges are balanced by block-nonzero
count plus row count. Preconditioner application shares the team when output
indices are unique; overlapping outputs retain sequential accumulation. No
nested OpenMP region or additional runnable history-pool work is introduced.
Dot products and squared norms use fixed 1024-element blocks and an ordered
final sum. This is repeatable for identical inputs/builds, but changes rounding
relative to the original serial sum. Existing small Hessenberg/Givens work is
performed by one thread. No fast-math or reduced precision is used.

`GmresCpuWorkspace` is explicitly owned by the multidomain solve session. V is
column contiguous and W is row contiguous. Large arrays and small matrices are
resized only when dimensions/restart change. A legacy caller retains its original
function signature and gets a temporary workspace. Matrix values are not cached
by object address. The workspace is noncopyable and must not be used concurrently.

## Diagnostics

Set `DBEM_GMRES_CPU_PROFILE=1` to accumulate `gmres.cpu.spmv`,
`preconditioner`, `orthogonalization`, `norm`, `vector_update`, `small_system`,
and `workspace_setup` timings. Orthogonalization includes its dot/update passes;
norm timings do not include orthogonalization dots. Synchronization within an
operation is included; unrelated iteration control/logging is not. Substage sums
must not be added to enclosing `multidomain.gmres` times. Metadata records scalar
dimension, block nonzeros, actual team size, iterations, preconditioner fallback
and final relative residual. File logging occurs only on the serial coordinator.
Profile mode adds barriers and must be disabled for speed comparisons.

The isolated build additionally defines `DBEM_GMRES_CPU_VERIFY=1`: after **either**
backend returns, the original serial `SparseMul` and a serial norm independently
check the true residual. Every physical step is logged in
`gmres_cpu_validation.csv`; a failed/nonfinite verification cannot report success.
This audit cost is included equally in both benchmark modes. Production builds
without this define retain the in-iteration residual test, without duplicate audit.

## Reproduction

Use WSL for ordinary commands and installed Windows Python/MSVC for compilation.
From the repository root (paths supplied to Windows Python may be relative):

```sh
WIN_PY=/mnt/c/Users/33435/AppData/Local/Programs/Python/Python312/python.exe
"$WIN_PY" tools/build_gmres_openmp.py
(cd tmp/gmres_openmp/build && ./gmres_openmp_test.exe)
"$WIN_PY" tools/validate_gmres_openmp.py --suite small
"$WIN_PY" tools/validate_gmres_openmp.py --suite large-short
"$WIN_PY" tools/validate_gmres_openmp.py --suite large-full
"$WIN_PY" tools/validate_gmres_openmp.py --suite diagnostic
```

The scripts run all solver computations with the card and environment set to 40
threads. Cases/backends run sequentially. Validation output is isolated and no
production card, input, executable, or output is overwritten. Completed results
can be resumed only with matching executable hashes; incomplete/failed runs are
retained and require inspection rather than automatic deletion. The full suite
uses the same cube materials as the original cases; TSV global material alone is
changed to copper E=110, nu=0.34, rho=8.96, with domain materials retained.

The small gate covers 24-step one-/two-domain rods and 1000-step TSV 01/04/09.
The short large gate covers 1000-/2000-domain cubes for 100 steps, three alternating
legacy/OMP pairs each. The final gate covers the 1000-domain cube for 1000 steps.
State comparison is streamed in chunks; row identity, finite fields, per-group
relative L2 <=1e-6 and maximum difference <=1e-5 of the baseline group peak are
required for both displacement and traction (zero floor 1e-30). Each step must
also meet its original true-residual tolerance. Accuracy warnings in historical
cube results remain a separate issue from parallel/legacy consistency.

The serial reference includes its existing 40-thread preconditioner application.
Speedups are therefore relative to the existing partially parallel implementation,
not a single-thread baseline. The runner reports medians for GMRES, run.solve
and run.total plus sampled peak resident memory and iteration changes. No thread
sweep is performed. A 1.5x GMRES speedup on **both** large cases is required before
considering default automatic activation.
