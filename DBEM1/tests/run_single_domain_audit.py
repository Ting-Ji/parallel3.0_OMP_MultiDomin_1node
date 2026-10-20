"""Build and audit legacy/modern one-domain GMRES in isolated directories.
Windows Python required. CSV comparisons cover the current even-step exports;
residual probes cover every solved step. Diagnostic timings are not benchmarks.
"""
from pathlib import Path
import argparse
import csv
import json
import math
import re
import time
import run_step0_policy_tests as h


def main():
    v = h.v
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path)
    parser.add_argument('--run-root', type=Path, default=h.ROOT / 'validation_runs/single_domain_audit_20261017')
    args = parser.parse_args()
    work = args.run_root.resolve()
    h.WORK = work
    exe = args.exe.resolve() if args.exe else h.build('diagnostic', diagnostics=True)
    summary = {'runs': [], 'comparisons': [], 'timings_are_diagnostic_not_benchmarks': True}
    run_cases = work / ('cases_' + time.strftime('%Y%m%d_%H%M%S'))
    def compare(a, b):
        def read(r):
            return {(int(x['step']), int(x['element'])): x for x in v.read_state(r.state_path)}
        left, right = read(a), read(b)
        assert left.keys() == right.keys()
        maximum = 0.
        passed = True
        worst = None
        for key in left:
            for field in v.STATE_FIELDS:
                x, y = left[key][field], right[key][field]
                assert math.isfinite(x) and math.isfinite(y)
                maximum = max(maximum, abs(x-y))
                if abs(x-y) > 1e-10 + 1e-7*max(abs(x), abs(y)):
                    passed = False
                    if worst is None or abs(x-y) > worst['difference']:
                        worst = dict(step=key[0], element=key[1], field=field, value_a=x, value_b=y, difference=abs(x-y))
        return dict(max_abs_difference=maximum, passed=passed, worst=worst)
    for label, mesh, steps in [('coarse', (1,1,1), 5), ('coarse_long', (1,1,1), 30), ('refined', (8,2,2), 30)]:
        v.generate_single_boundary_mesh.__defaults__ = mesh
        v.generate_multidomain_mesh.__defaults__ = mesh
        results = {}
        for modern in (False, True):
            for threads in (1,4):
                v.DEFAULT_THREADS = threads
                name = f'{label}_{"modern" if modern else "legacy"}_t{threads}'
                spec = v.CaseSpec(name, 1, steps, modern)
                case = v.prepare_case(run_cases, h.ROOT, exe, spec)
                started = time.perf_counter()
                r = v.run_case(case, spec, 180)
                seconds = time.perf_counter()-started
                assert r.returncode == 1, r.log_text[-4000:]
                elements = int(r.metrics['ElementCount'])
                assert {int(row['step']) for row in v.read_state(r.state_path)} == set(range(0,steps+1,2))
                assert len(v.read_state(r.state_path)) == (steps//2+1)*elements
                with (v.result_output_dir(r)/'residual_probe.csv').open(newline='') as stream:
                    probes = list(csv.DictReader(stream))
                assert [int(p['step']) for p in probes] == list(range(1,steps+1))
                residual = max(float(p['relativeResidual']) for p in probes)
                assert math.isfinite(residual) and residual <= 1e-8, residual
                if modern:
                    assert r.metrics['MultiDomainFinalFlag'] == '0'
                    assert int(r.metrics['MultiDomainSolvedSteps']) == steps
                    assert r.metrics['MultiDomainMatrixStructureValid'] == '1'
                checks = re.findall(r'maxDiff=([^;\s]+)', r.log_text)
                assert all(float(x) <= 1e-10 for x in checks), checks
                summary['runs'].append(dict(name=name, elements=elements, steps=steps, threads=threads,
                    max_relative_residual=residual, wall_seconds=seconds,
                    modern_inner_iterations=r.metrics.get('MultiDomainTotalInnerIterations'),
                    assembly_check_max_difference=max(map(float,checks),default=None)))
                results[modern,threads] = r
                print(f'PASS {name}: elements={elements}, residual={residual:.3e}, wall={seconds:.3f}s',flush=True)
        for title,a,b in [('legacy threads',(False,1),(False,4)),('modern threads',(True,1),(True,4)),
                          ('legacy/modern t1',(False,1),(True,1)),('legacy/modern t4',(False,4),(True,4))]:
            summary['comparisons'].append(dict(name=f'{label}: {title}', **compare(results[a],results[b])))
    summary['passed'] = all(c['passed'] for c in summary['comparisons'])
    summary['case_directory'] = str(run_cases)
    (work/'single_domain_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2),flush=True)

    return 0 if summary['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
