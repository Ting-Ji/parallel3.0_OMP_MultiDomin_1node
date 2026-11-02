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
For at least 65,536 scalar unknowns, dot products and squared norms use fixed
1024-element blocks and an ordered final sum. Smaller forced-OMP systems retain
the original serial summation order: TSV 04's long-time traction failed the
strict field gate with changed reduction order despite passing every residual
test (the initial evidence is retained). Large-system fixed-block reduction is
repeatable for identical inputs/builds, but changes rounding relative to the
original serial sum. Existing small Hessenberg/Givens work is
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

## Initial validation and compatibility adjustment (2026-11-02)

The first build passed the 40-thread operator/solver tests, both 24-step rods,
and the full 1000-step TSV 01 comparison. TSV 04 passed every independent linear
residual check but failed the full-history traction gate: relative L2 7.0112e-5,
maximum absolute difference 5.5795e-7. Displacement relative L2 was 3.6150e-7.
The two cards were identical (SHA256
`1698256c7e901a2bdefa673039214291ecad8d11fef2b41e7f3e907ff5035d6f`),
with global copper E=110, nu=0.34, rho=8.96 and 40 threads. The printed c1/c2 in
DBEM.cpp uses materialContext.Get(0), the first silicon domain, not the global
card material. The copper global time-step path was retained.

The first differences were at roundoff scale in step 1, growing late in the
history. The small-system serial summation compatibility path was added; large
systems retain the planned parallel fixed-block reductions. The rebuilt tests
also exercise 65,541 unknowns, explicitly covering the fixed-block path. A fresh
100-step TSV 04 legacy/OMP pair then produced exactly equal displacement and
traction fields, with equal iteration totals and passing independent residuals.
This short check does not replace the full 1000-step gate.

Initial failed evidence is retained in `tmp/gmres_openmp/acceptance/`; the
compatibility check is in `tmp/gmres_openmp/quick_v2/`. The new complete acceptance
queue uses `tmp/gmres_openmp/acceptance_v2/` and the rebuilt solver in
`tmp/gmres_openmp/build_v2/`. Its live report is `acceptance_report.md` and its
machine-readable status is `acceptance_queue.json`. Run it from Windows Python:

```sh
"$WIN_PY" tools/run_gmres_openmp_acceptance.py \
  --exe tmp/gmres_openmp/build_v2/DBEM1_openmp.exe \
  --output tmp/gmres_openmp/acceptance_v2
```

The queue stops on any failed gate, retains all evidence and never changes the
default backend automatically. Default activation requires review after all
numerical gates and the large-case performance target have completed.
