"""One two-domain five-step case plus cache units, isolated from step 0 results."""
from pathlib import Path
import argparse
import json
import run_step0_policy_tests as tests

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--baseline', type=Path, help='Saved pre-change case output directory')
args = parser.parse_args()
saved = tests.ROOT / 'validation_runs/step0_policy/case_validation/case_md_2domain_n5/output/rod_validation_5'
baseline = args.baseline or (saved if saved.is_dir() else None)
if baseline is not None:
    assert (baseline / 'validation_state.csv').is_file(), 'Invalid baseline output directory'
tests.WORK = tests.ROOT / 'validation_runs/step1_cache'
tests.build('unit', unit=True)
default = tests.run_small(tests.build('default'), 'case_default', False)
validation = tests.run_small(tests.build('validation', diagnostics=True), 'case_validation', True)
files = list(default.glob('TecValueFile_1node_*.dat'))
assert len(files) == 6
for path in files:
    assert path.read_bytes() == (validation / path.name).read_bytes(), path.name
    if baseline is not None:
        assert path.read_bytes() == (baseline / path.name).read_bytes(), path.name
if baseline is not None:
    assert (validation / 'validation_state.csv').read_bytes() == (baseline / 'validation_state.csv').read_bytes()
summary = {'two_domain_steps': 5, 'cache_units': 'passed',
           'step0_baseline_byte_identical': True if baseline is not None else None,
           'default_vs_validation_byte_identical': True,
           'cross_check_max_difference': 0}
(tests.WORK / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary), flush=True)
