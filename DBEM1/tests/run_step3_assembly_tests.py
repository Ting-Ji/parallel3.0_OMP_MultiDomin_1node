"""Unified assembly queue: small two-domain five-step regression only."""
import json
import shutil
import run_step0_policy_tests as tests

tests.WORK = tests.ROOT / 'validation_runs/step3_assembly'
prior = tests.ROOT / 'validation_runs/step2_step0/case_parallel/case_md_2domain_n5/output/rod_validation_5'
baseline = tests.WORK / 'baseline_step2'
names = ['validation_state.csv'] + [f'TecValueFile_1node_{i}.dat' for i in range(6)]
assert all((prior / name).is_file() for name in names), 'Step 2 baseline required.'
baseline.mkdir(parents=True, exist_ok=True)
for name in names:
    if not (baseline / name).exists():
        shutil.copy2(prior / name, baseline / name)
tests.build('unit_strict', unit=True)
tests.build('unit_fallback', unit=True, fallback=True)
exe = tests.build('default')
outputs = [tests.run_small(exe, 'case_single', False, threads=1),
           tests.run_small(exe, 'case_parallel', False, threads=4),
           tests.run_small(tests.build('diagnostic', diagnostics=True), 'case_diagnostic', True)]
for name in names:
    expected = (baseline / name).read_bytes()
    for output in outputs:
        assert (output / name).read_bytes() == expected, (output, name)
summary = {'domains': 2, 'steps': 5, 'threads': [1, 4],
           'step2_baseline_byte_identical': True, 'cross_check_max_difference': 0,
           'default_optional_diagnostics_disabled': True, 'state_always_written': True,
           'cross_domain_queue_pool_reuse_and_failure_units': 'passed',
           'strict_and_explicit_resource_fallback': 'passed'}
(tests.WORK / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary), flush=True)
