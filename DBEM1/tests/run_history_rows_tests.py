"""Windows/MSVC regression runner: tiny C++ tests plus three small cases once.

Reuses validation_runs/history_rhs/cases as the existing numerical baseline.
Never rebuilds or runs the baseline executable.
"""
from pathlib import Path
import json
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / 'validation_runs/history_rows'


def prepare_project():
    ns = {'ms': 'http://schemas.microsoft.com/developer/msbuild/2003'}
    ET.register_namespace('', ns['ms'])
    tree = ET.parse(ROOT / 'DBEM1/DBEM1.vcxproj')
    project = tree.getroot()
    for item in project.findall('.//ms:ItemGroup/ms:ClCompile', ns):
        name = item.attrib['Include']
        item.set('Include', str(ROOT / 'DBEM1' / name))
        if name == 'multidomain.cpp':
            item.set('Include', str(ROOT / 'DBEM1/tests/history_rows_test.cpp'))
        if name == 'DBEM.cpp':
            ET.SubElement(item, '{' + ns['ms'] + '}PreprocessorDefinitions').text = 'main=dbem_application_main;%(PreprocessorDefinitions)'
    for group in project.findall('.//ms:ItemDefinitionGroup', ns):
        for item in group.findall('.//ms:AdditionalIncludeDirectories', ns):
            item.text = str(ROOT / 'include') + ';%(AdditionalIncludeDirectories)'
        for item in group.findall('.//ms:AdditionalLibraryDirectories', ns):
            item.text = str(ROOT / 'lib/x64') + ';%(AdditionalLibraryDirectories)'
    work = WORK / 'unit'
    work.mkdir(parents=True, exist_ok=True)
    for name in ('gaussposition', 'gaussweight'):
        shutil.copy2(ROOT / 'DBEM1' / name, work / name)
    path = work / 'history_rows_test.vcxproj'
    tree.write(path, encoding='utf-8', xml_declaration=True)
    return path


def run_unit_tests():
    project = prepare_project()
    msbuild = Path(os.environ.get('ProgramFiles', 'C:/Program Files')) / 'Microsoft Visual Studio/2022/Community/MSBuild/Current/Bin/MSBuild.exe'
    build = subprocess.run([str(msbuild), str(project), '/p:Configuration=Release', '/p:Platform=x64', '/m:1', '/v:quiet', '/clp:ErrorsOnly'],
                           cwd=project.parent, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (WORK / 'unit_build.log').write_bytes(build.stdout)
    if build.returncode:
        print(build.stdout.decode(errors='replace'))
        raise RuntimeError('C++ test build failed')
    exe = project.parent / 'x64/Release/history_rows_test.exe'
    run = subprocess.run([str(exe)], cwd=project.parent, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
    (WORK / 'unit.log').write_bytes(run.stdout)
    print(run.stdout.decode(errors='replace'))
    if run.returncode:
        raise RuntimeError(f'C++ tests failed: {run.returncode}')


def run_small_cases():
    sys.path.insert(0, str(ROOT / 'tools'))
    import run_multidomain_validation as v
    baseline_root = ROOT / 'validation_runs/history_rhs/cases'
    specs = [spec for spec in v.case_specs(True, False)
             if spec.name in ('md_1domain_n5', 'md_2domain_n5', 'md_10domain_n5')]
    os.environ['DBEM_MULTIDOMAIN_GMRES_PTHREAD'] = '1'
    # Use an executable compiled with assemblyCrossCheck and validationOutput.
    report = []
    for spec in specs:
        v.DEFAULT_THREADS = 1 if spec.domain_count == 1 else 4
        case = v.prepare_case(WORK / 'cases', ROOT, ROOT / 'DBEM1/x64/Release/DBEM1.exe', spec)
        result = v.run_case(case, spec, 120)
        if 'DBEM policy:' in result.log_text and 'assemblyCrossCheck=1' not in result.log_text:
            raise RuntimeError('Rebuild with DBEM_POLICY_ASSEMBLY_CHECK=1; environment switches no longer enable checks.')
        checks = [v.check_md_run(result, require_self_check=spec.domain_count == 1)]
        if spec.domain_count > 1:
            checks.append(v.check_interface(result, 1e-10, 1e-8))
        if not all(check.passed for check in checks):
            raise AssertionError(str(checks))
        old = baseline_root / ('case_' + spec.name) / 'output' / f'rod_validation_{spec.nstep}' / 'validation_state.csv'
        identical = old.read_bytes() == result.state_path.read_bytes()
        if not identical:
            raise AssertionError(f'{spec.name}: state differs from saved baseline')
        history_checks = [line for line in result.log_text.splitlines() if 'pthread check History' in line]
        if spec.domain_count > 1:
            if not history_checks or not all(all(float(x) == 0 for x in re.findall(r'maxDiff=([^;\s]+)', line)) for line in history_checks):
                raise AssertionError(f'{spec.name}: history matrix mismatch')
            if 'falling back' in result.log_text or 'pthread_create failed' in result.log_text:
                raise AssertionError(f'{spec.name}: unexpected fallback')
            timing = (v.result_output_dir(result) / 'timing.txt').read_text()
            if 'multidomain.assembly.history.output.fill' not in timing:
                raise AssertionError('Parallel output path did not execute')
        report.append({'case': spec.name, 'threads': v.DEFAULT_THREADS, 'state_byte_identical': identical,
                       'history_checks': history_checks, 'checks': [check.message for check in checks]})
        print(f'PASS: {spec.name}, threads={v.DEFAULT_THREADS}, state identical, history checks={len(history_checks)}')
    (WORK / 'summary.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print('Known cross-domain-vs-single comparisons were not rerun; existing baseline reports are retained.')


if __name__ == '__main__':
    if os.name != 'nt':
        raise SystemExit('Run with Windows Python (MSVC and Windows executable required).')
    WORK.mkdir(parents=True, exist_ok=True)
    run_unit_tests()
    run_small_cases()
