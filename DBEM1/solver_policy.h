#ifndef DBEM_SOLVER_POLICY_H
#define DBEM_SOLVER_POLICY_H

#include <cstdio>
#include <cstdlib>

// Built-in policy: edit these defaults and rebuild. /D overrides are intended
// for isolated validation builds only; environment variables never enable them.
#ifndef DBEM_POLICY_VALIDATION_OUTPUT
#define DBEM_POLICY_VALIDATION_OUTPUT 0
#endif
#ifndef DBEM_POLICY_ASSEMBLY_CHECK
#define DBEM_POLICY_ASSEMBLY_CHECK 0
#endif
#ifndef DBEM_POLICY_SOLVE_DIAGNOSTICS
#define DBEM_POLICY_SOLVE_DIAGNOSTICS 0
#endif
#ifndef DBEM_POLICY_SERIAL_FALLBACK
#define DBEM_POLICY_SERIAL_FALLBACK 0
#endif
#ifndef DBEM_POLICY_REFERENCE_REPLACEMENT
#define DBEM_POLICY_REFERENCE_REPLACEMENT 0
#endif

namespace DBEMPolicy {
// Optional metrics/audits only. validation_state.csv is always written on success.
constexpr bool validationOutput = DBEM_POLICY_VALIDATION_OUTPUT != 0;
constexpr bool assemblyCrossCheck = DBEM_POLICY_ASSEMBLY_CHECK != 0;
constexpr bool solveDiagnostics = DBEM_POLICY_SOLVE_DIAGNOSTICS != 0;
constexpr bool allowSerialFallback = DBEM_POLICY_SERIAL_FALLBACK != 0;
constexpr bool allowReferenceReplacement = DBEM_POLICY_REFERENCE_REPLACEMENT != 0;

inline void Print()
{
    printf("DBEM policy: validationOutput=%d assemblyCrossCheck=%d solveDiagnostics=%d allowSerialFallback=%d allowReferenceReplacement=%d\n",
        validationOutput, assemblyCrossCheck, solveDiagnostics,
        allowSerialFallback, allowReferenceReplacement);
    const char* legacy[] = { "DBEM_VALIDATION_OUTPUT", "DBEM_RHS_DIAGNOSTIC",
        "DBEM_MULTIDOMAIN_GMRES_PTHREAD_CHECK", "DBEM_PARALLEL_KNOWN_RHS_CHECK" };
    for (const char* name : legacy) {
        char* value = nullptr;
        size_t size = 0;
        if (!_dupenv_s(&value, &size, name) && value)
            printf("DBEM policy: %s is ignored; configure solver_policy.h and rebuild.\n", name);
        free(value);
    }
}
}
#endif
