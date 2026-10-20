"""Step 2: one two-domain five-step case, 1/4 workers, plus bounded unit checks."""
import json
import os
import subprocess
import shutil
import run_step0_policy_tests as tests

tests.WORK = tests.ROOT / 'validation_runs/step2_step0'
baseline = tests.ROOT / 'validation_runs/step1_cache/case_validation/case_md_2domain_n5/output/rod_validation_5'
assert (baseline / 'validation_state.csv').is_file(), 'Capture the step 1 baseline first.'
tests.build('unit', unit=True)
exe = tests.build('default')
single = tests.run_small(exe, 'case_single', False, threads=1)
parallel = tests.run_small(exe, 'case_parallel', False, threads=4)
diagnostic = tests.run_small(tests.build('diagnostic', diagnostics=True), 'case_diagnostic', True)
for name in ['validation_state.csv'] + [f'TecValueFile_1node_{i}.dat' for i in range(6)]:
    expected = (baseline / name).read_bytes()
    for output in (single,parallel,diagnostic):
        assert (output / name).read_bytes() == expected, (output,name)
log = (diagnostic.parents[1] / 'run.log').read_text(encoding='utf-8')
assert 'Step0 BCID immutability check passed.' in log

# Reuse only a copied case for a deliberate failed rerun: old success output
# must be invalidated before the solver reports an error, with no final/temp CSV.
source = parallel.parents[1]
failed_case = tests.WORK / 'case_failed'
shutil.copytree(source, failed_case, dirs_exist_ok=True)
(tests.WORK / 'output').mkdir(exist_ok=True)  # Legacy ../output/log.txt.
env = os.environ.copy()
env['DBEM_MULTIDOMAIN_TREE_ROOT'] = 'invalid-test-mode'
run = subprocess.run([str(failed_case / 'DBEM1.exe')], cwd=failed_case, env=env,
                     capture_output=True, timeout=120)
(failed_case / 'run.log').write_bytes(run.stdout + run.stderr)
assert run.returncode not in (0,1)
assert b'Invalid DBEM_MULTIDOMAIN_TREE_ROOT' in run.stdout + run.stderr
assert not (failed_case / 'output/rod_validation_5/validation_state.csv').exists()
assert not (failed_case / 'output/rod_validation_5/validation_state.csv.tmp').exists()
summary = {'domains':2, 'steps':5, 'threads':[1,4], 'step1_baseline_byte_identical':True,
           'default_state_always_written':True, 'cross_check_max_difference':0,
           'pool_reuse_slots_materials_and_failure_tests':'passed', 'failed_run_stale_state_removed':True}
(tests.WORK / 'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary),flush=True)
