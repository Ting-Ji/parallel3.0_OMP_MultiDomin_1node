"""Small Windows/MSVC regression: policy units and one two-domain five-step case.

All generated projects, executables and cases stay in validation_runs/step0_policy.
Optional --baseline points to a saved pre-change case output directory.
"""
from pathlib import Path
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / 'validation_runs/step0_policy'
sys.path.insert(0, str(ROOT / 'tools'))
import run_multidomain_validation as v


def build(name, *, unit=False, diagnostics=False, fallback=False):
    ns = {'ms': 'http://schemas.microsoft.com/developer/msbuild/2003'}
    ET.register_namespace('', ns['ms'])
    tree = ET.parse(ROOT / 'DBEM1/DBEM1.vcxproj')
    project = tree.getroot()
    work = WORK / name
    work.mkdir(parents=True, exist_ok=True)
    for helper in ('gaussposition', 'gaussweight'):
        shutil.copy2(ROOT / 'DBEM1' / helper, work / helper)
    for item in project.findall('.//ms:ItemGroup/ms:ClCompile', ns):
        source = item.attrib['Include']
        item.set('Include', str(ROOT / 'DBEM1' / source))
        if unit and source == 'multidomain.cpp':
            item.set('Include', str(ROOT / 'DBEM1/tests/history_rows_test.cpp'))
        if unit and source == 'DBEM.cpp':
            ET.SubElement(item, '{' + ns['ms'] + '}PreprocessorDefinitions').text = (
                'main=dbem_application_main;%(PreprocessorDefinitions)')
    defines = [f'DBEM_POLICY_SERIAL_FALLBACK={int(fallback)}']
    if diagnostics:
        defines += ['DBEM_POLICY_VALIDATION_OUTPUT=1', 'DBEM_POLICY_ASSEMBLY_CHECK=1',
                    'DBEM_POLICY_SOLVE_DIAGNOSTICS=1']
    for group in project.findall('.//ms:ItemDefinitionGroup', ns):
        for item in group.findall('.//ms:AdditionalIncludeDirectories', ns):
            item.text = str(ROOT / 'include') + ';%(AdditionalIncludeDirectories)'
        for item in group.findall('.//ms:AdditionalLibraryDirectories', ns):
            item.text = str(ROOT / 'lib/x64') + ';%(AdditionalLibraryDirectories)'
        for compile_group in group.findall('ms:ClCompile', ns):
            node = compile_group.find('ms:PreprocessorDefinitions', ns)
            if node is None:
                node = ET.SubElement(compile_group, '{' + ns['ms'] + '}PreprocessorDefinitions')
            node.text = ';'.join(defines) + ';' + (node.text or '%(PreprocessorDefinitions)')
    path = work / (name + '.vcxproj')
    tree.write(path, encoding='utf-8', xml_declaration=True)
    msbuild = Path(os.environ.get('ProgramFiles', 'C:/Program Files')) / (
        'Microsoft Visual Studio/2022/Community/MSBuild/Current/Bin/MSBuild.exe')
    result = subprocess.run([str(msbuild), str(path), '/p:Configuration=Release',
                             '/p:Platform=x64', '/m:1', '/v:quiet', '/clp:ErrorsOnly'],
                            cwd=work, capture_output=True)
    (work / 'build.log').write_bytes(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError((result.stdout + result.stderr).decode(errors='replace'))
    exe = work / 'x64/Release' / (name + '.exe')
    if unit:
        result = subprocess.run([str(exe)], cwd=work, capture_output=True, timeout=60)
        (work / 'unit.log').write_bytes(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f'Unit exit code {result.returncode}: ' +
                               (result.stdout + result.stderr).decode(errors='replace'))
        print(f'PASS {name}', flush=True)
    return exe


def run_small(exe, name, diagnostics, threads=4):
    v.DEFAULT_THREADS = threads
    spec = next(s for s in v.case_specs(False, False) if s.name == 'md_2domain_n5')
    case = v.prepare_case(WORK / name, ROOT, exe, spec)
    result = v.run_case(case, spec, 120, require_validation=diagnostics)
    assert result.returncode == 1, result.log_text[-4000:]
    flags = re.findall(r'MultiDomain GMRES step (\d+): flag=(-?\d+)', result.log_text)
    assert flags == [(str(i), '0') for i in range(1, 6)], flags
    assert 'explicitly enabled' not in result.log_text
    output = v.result_output_dir(result)
    if diagnostics:
        assert v.check_md_run(result).passed
        assert v.check_interface(result, 1e-10, 1e-8).passed
        residual = v.check_residual_probe(result, abs_tol=1e-10, rel_tol=1e-8)
        assert residual.passed, residual
        checks = [line for line in result.log_text.splitlines() if 'pthread check ' in line]
        assert any('Step0' in line for line in checks)
        assert any('History' in line for line in checks)
        assert all(float(x) == 0 for line in checks for x in re.findall(r'maxDiff=([^;\s]+)', line))
    else:
        assert sorted(p.name for p in output.glob('*.csv')) == ['validation_state.csv']
        assert not (output / 'validation_metrics.txt').exists()
        assert 'pthread check ' not in result.log_text
        assert 'GMRES candidate' not in result.log_text
        assert 'self-check' not in result.log_text
    assert (output / 'validation_state.csv').is_file()
    print(f'PASS {name}: two domains, five steps', flush=True)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--units-only', action='store_true')
    args = parser.parse_args()
    if os.name != 'nt':
        raise SystemExit('Run with Windows Python: MSVC and Windows executables required.')
    build('unit_strict', unit=True)
    build('unit_fallback', unit=True, fallback=True)
    if args.units_only:
        return
    default = run_small(build('default'), 'case_default', False)
    validation = run_small(build('validation', diagnostics=True), 'case_validation', True)
    tec = list(default.glob('TecValueFile_1node_*.dat'))
    assert len(tec) == 6
    for file in tec:
        assert file.read_bytes() == (validation / file.name).read_bytes(), file.name
        if args.baseline:
            assert file.read_bytes() == (args.baseline / file.name).read_bytes(), file.name
    if args.baseline:
        assert (validation / 'validation_state.csv').read_bytes() == (args.baseline / 'validation_state.csv').read_bytes()
    summary = {'two_domain_steps': 5, 'default_vs_validation_byte_identical': True,
               'baseline_byte_identical': True if args.baseline else None,
               'cross_check_max_difference': 0, 'strict_and_opt_in_fallback_units': 'passed'}
    (WORK / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
