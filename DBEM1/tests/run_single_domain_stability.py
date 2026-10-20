"""New-framework one-domain stability audit. Run using Windows Python.
All cases are isolated; production source/card/executable are not modified.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, os, re, shutil, sys, time
import run_step0_policy_tests as h
v = h.v

def triangle(tau):
    if tau <= 0:
        return 0.0
    phase = tau % 4.0
    return phase if phase <= 2 else 4-phase

def run_one(work, exe, config):
    name = config['name']
    mesh = tuple(config.get('mesh', (8,1,1)))
    nu = config.get('nu', 0.0)
    e = config.get('e', 20000.)
    rho = config.get('rho', 8.)
    beta = config.get('beta', .4)
    target = config.get('tau', 12.)
    threads = config.get('threads', 4)
    amplitude = config.get('amplitude', 1.)
    load = config.get('load', 'step')
    c = math.sqrt(e/rho)
    c1 = math.sqrt(e*(1-nu)/((1+nu)*(1-2*nu)*rho))
    half = config.get('half_width', .1)
    nx,ny,nz = mesh
    a=half
    builder=v.RodMeshBuilder()
    faces=[
        ('left',(0,-a,a),(0,a,a),(0,a,-a),(0,-a,-a),ny,nz),
        ('bottom',(0,-a,-a),(1,-a,-a),(1,-a,a),(0,-a,a),nx,nz),
        ('top',(0,a,-a),(0,a,a),(1,a,a),(1,a,-a),nz,nx),
        ('front',(0,-a,-a),(0,a,-a),(1,a,-a),(1,-a,-a),ny,nx),
        ('back',(0,-a,a),(1,-a,a),(1,a,a),(0,a,a),nx,ny),
        ('right',(1,-a,-a),(1,a,-a),(1,a,a),(1,-a,a),ny,nz)]
    for face in faces:
        builder.add_subdivided_face(*face)
    hmin=min(math.dist(builder.points[ele[j]-1],builder.points[ele[(j+1)%4]-1])
             for ele in builder.elements for j in range(4))
    dt_expected = beta*hmin/c1
    steps = int(math.ceil(target/(c*dt_expected)))
    v.DEFAULT_E, v.DEFAULT_V, v.DEFAULT_ROU = e, nu, rho
    v.DEFAULT_THREADS = threads
    v.DEFAULT_PRE_LEAF = config.get('leaf', 10)
    v.DEFAULT_ERROR = config.get('tolerance', 1e-8)
    spec = v.CaseSpec(name, 1, steps, True, beta_override=beta)
    case = work / ('case_'+name)
    (case/'input').mkdir(parents=True)
    (case/'output').mkdir()
    (work/'output').mkdir(exist_ok=True)
    shutil.copy2(exe, case/'DBEM1.exe')
    for helper in ('gaussposition','gaussweight'):
        shutil.copy2(h.ROOT/'DBEM1'/helper, case/helper)
    v.write_model(case/'input/rod_validation.model', builder)
    v.write_datacard(case/'BEM_DATACARD.DAT', 'rod_validation', spec, builder)
    card=case/'BEM_DATACARD.DAT'
    lines=card.read_text().splitlines()
    lines[17]=f'1 {e:.17g} {nu:.17g} {rho:.17g}'
    card.write_text('\n'.join(lines)+'\n',encoding='ascii')
    count = len(builder.elements)
    def value(step):
        tau = step*dt_expected*c
        if step == 0 or load == 'zero':
            return 0.
        if load == 'pulse' and tau > 1.:
            return 0.
        return amplitude
    with (case/'input/rod_validation.bd').open('w') as fp:
        for face in builder.face_names:
            fp.write('123\n' if face == 'left' else '456\n')
        for step in range(steps+1):
            for face in builder.face_names:
                # Local direction 3 is the outward normal; right end normal = +x.
                q = value(step) if face == 'right' else 0.
                fp.write(f'0 0 {q:.17g}\n')
    saved = os.environ.get('DBEM_MULTIDOMAIN_GMRES_PTHREAD')
    if config.get('serial_assembly'):
        os.environ['DBEM_MULTIDOMAIN_GMRES_PTHREAD'] = '0'
    else:
        os.environ.pop('DBEM_MULTIDOMAIN_GMRES_PTHREAD',None)
    start = time.perf_counter()
    try:
        result = v.run_case(case, spec, 300, require_validation=True)
    finally:
        if saved is None:
            os.environ.pop('DBEM_MULTIDOMAIN_GMRES_PTHREAD',None)
        else:
            os.environ['DBEM_MULTIDOMAIN_GMRES_PTHREAD'] = saved
    seconds = time.perf_counter()-start
    output = v.result_output_dir(result)
    record = dict(config, steps=steps, elements=count, returncode=result.returncode,
                  seconds=seconds, case=str(case), output=str(output), expected_dt=dt_expected)
    flags = re.findall(r'MultiDomain GMRES step (\d+): flag=(-?\d+) outer=(\d+) inner=(\d+)', result.log_text)
    record['flag_count'] = len(flags)
    record['nonzero_flags'] = [x for x in flags if x[1]!='0']
    record['solved_all_steps'] = (len(flags)==steps and all(x[1]=='0' for x in flags))
    record['max_inner_iterations'] = max((int(x[3]) for x in flags),default=None)
    record['assembly_max_difference'] = max((float(x) for x in re.findall(r'maxDiff=([^;\s]+)',result.log_text)),default=None)
    record['matrix_structure_valid'] = result.metrics.get('MultiDomainMatrixStructureValid')
    record['tree_mode'] = result.metrics.get('TreeRootMode')
    if not result.state_path.exists():
        record['passed_algebra'] = False
        return record, None
    dt = float(result.metrics['Dt'])
    record['dt'] = dt
    assert abs(dt-dt_expected) < 1e-12*dt_expected
    rows = v.read_state(result.state_path)
    record['state_rows'] = len(rows)
    record['finite_state'] = all(math.isfinite(row[f]) for row in rows for f in v.STATE_FIELDS)
    assert len(rows)==count*(steps//2+1), (name,len(rows),steps,count)
    right = {i for i,f in enumerate(builder.face_names) if f=='right'}
    left = {i for i,f in enumerate(builder.face_names) if f=='left'}
    by_step = {}
    bc_error = 0.
    for row in rows:
        step, ele = int(row['step']),int(row['element'])
        if ele in right:
            by_step.setdefault(step,[]).append(row['ux'])
        for f in ('ux','uy','uz') if ele in left else ('tx','ty','tz'):
            expected = value(step) if ele in right and f=='tx' else 0.
            bc_error = max(bc_error,abs(row[f]-expected))
    record['known_boundary_max_error'] = bc_error
    scale = 1/e
    record['max_u_over_unit_static'] = max(abs(row[f])/scale for row in rows for f in ('ux','uy','uz'))
    record['max_t_over_unit_load'] = max(abs(row[f]) for row in rows for f in ('tx','ty','tz'))
    curve = []
    for step,vals in sorted(by_step.items()):
        tau = c*step*dt
        expected = 0. if load=='zero' else amplitude*(triangle(tau)-(triangle(tau-1) if load=='pulse' else 0))
        numeric = sum(vals)/len(vals)/scale
        curve.append(dict(step=step,tau=tau,tip_over_unit_static=numeric,
                          analytic_nu0=expected if nu==0 else None))
    with (output/'stability_tip.csv').open('w',newline='') as fp:
        writer=csv.DictWriter(fp,fieldnames=curve[0].keys())
        writer.writeheader();writer.writerows(curve)
    record['early_tip_peak'] = max(abs(x['tip_over_unit_static']) for x in curve if x['tau']<=4.)
    late=[abs(x['tip_over_unit_static']) for x in curve if x['tau']>=max(4.,target-4)]
    record['late_tip_peak'] = max(late,default=None)
    record['growth_late_over_early'] = record['late_tip_peak']/max(record['early_tip_peak'],1e-300) if late else None
    if nu==0 and load!='zero':
        errors = [x['tip_over_unit_static']-x['analytic_nu0'] for x in curve]
        record['analytic_max_error_over_static'] = max(map(abs,errors))
        record['analytic_rmse_over_static'] = math.sqrt(sum(x*x for x in errors)/len(errors))
        record['bounded_tip_screen'] = max(abs(x['tip_over_unit_static']) for x in curve)<=10*abs(amplitude)
    with (output/'residual_probe.csv').open(newline='') as fp:
        probes = list(csv.DictReader(fp))
    record['max_relative_residual'] = max(float(x['relativeResidual']) for x in probes)
    record['max_abs_residual'] = max(float(x['maxAbsResidual']) for x in probes)
    record['residual_all_steps'] = [int(x['step']) for x in probes] == list(range(1,steps+1))
    record['passed_algebra'] = (result.returncode==1 and record['solved_all_steps'] and record['finite_state']
        and record['matrix_structure_valid']=='1' and record['residual_all_steps']
        and record['max_relative_residual']<=max(1.01*v.DEFAULT_ERROR,1e-12)
        and bc_error<1e-12)
    threshold=10*abs(amplitude)
    record['stability_screen_passed'] = (record['passed_algebra'] and
        record['max_u_over_unit_static'] <= threshold)
    record['stability_screen_definition'] = 'finite/converged and max displacement <= 10*abs(load)*L/E; zero load must stay exactly zero'
    record['state_sha256'] = hashlib.sha256(result.state_path.read_bytes()).hexdigest()
    return record, rows

def compare(a,b,scale=1.,traction_scale=None):
    assert len(a)==len(b)
    maximum=0.
    passed=True
    for p,q in zip(a,b):
        assert (p['step'],p['element'])==(q['step'],q['element'])
        for f in v.STATE_FIELDS:
            factor=traction_scale if f.startswith('t') and traction_scale is not None else scale
            x,y=p[f]*factor,q[f]
            diff=abs(x-y)
            maximum=max(maximum,diff)
            passed=passed and math.isfinite(diff) and diff<=1e-10+1e-7*max(abs(x),abs(y))
    return dict(passed=passed,max_abs_difference=maximum)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',type=Path)
    parser.add_argument('--run-root',type=Path,default=h.ROOT/'validation_runs/single_domain_stability_20261017')
    parser.add_argument('--phase',choices=('base','extended'),default='base')
    args=parser.parse_args()
    work=args.run_root.resolve()
    work.mkdir(parents=True,exist_ok=True)
    h.WORK=work
    exe=args.exe.resolve() if args.exe else h.build('diagnostic',diagnostics=True)
    run_dir=work/('runs_'+time.strftime('%Y%m%d_%H%M%S'))
    run_dir.mkdir()
    if args.phase=='base':
        configs=[
            dict(name='axial_t1',threads=1),
            dict(name='axial_t4',threads=4),
            dict(name='axial_t8',threads=8),
            dict(name='axial_t4_repeat',threads=4),
            dict(name='axial_serial_assembly',threads=4,serial_assembly=True),
            dict(name='zero_136',mesh=(16,2,2),tau=8.,load='zero'),
            dict(name='pulse',load='pulse'),
            dict(name='beta02',beta=.2),
            dict(name='beta06',beta=.6),
            dict(name='beta08',beta=.8),
            dict(name='refined136',mesh=(16,2,2)),
            dict(name='poisson025',nu=.25),
            dict(name='poisson045',nu=.45,tau=4.),
            dict(name='long_tau24',tau=24.),
            dict(name='load_negative',amplitude=-1.),
            dict(name='material_scale',e=80000.,rho=32.,amplitude=4.),
        ]
    else:
        configs=[dict(name='refined_beta02',mesh=(16,2,2),beta=.2),
                 dict(name='refined_beta06',mesh=(16,2,2),beta=.6),
                 dict(name='refined_long',mesh=(16,2,2),tau=24.),
                 dict(name='cube96',mesh=(4,4,4),half_width=.5,tau=12.),
                 dict(name='tolerance10',tolerance=1e-10),
                 dict(name='leaf4',leaf=4),dict(name='leaf20',leaf=20)]
    summary=dict(executable=str(exe),phase=args.phase,runs=[],comparisons=[])
    states={}
    for config in configs:
        try:
            record,rows=run_one(run_dir,exe,config)
            if rows is not None: states[config['name']]=rows
        except Exception as exc:
            record=dict(config,error=repr(exc),passed_algebra=False)
        summary['runs'].append(record)
        (run_dir/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
        print(json.dumps({k:record.get(k) for k in ('name','passed_algebra','steps','elements','seconds','max_relative_residual','max_u_over_unit_static','analytic_rmse_over_static','growth_late_over_early','error')}),flush=True)
    if args.phase=='base':
        for name in ('axial_t4','axial_t8','axial_t4_repeat','axial_serial_assembly','load_negative','material_scale'):
            if 'axial_t1' in states and name in states:
                factor=-1. if name=='load_negative' else 1.
                summary['comparisons'].append(dict(name='axial_t1 vs '+name,**compare(states['axial_t1'],states[name],factor,4. if name=='material_scale' else None)))
    summary['algebra_passed']=all(x['passed_algebra'] for x in summary['runs'])
    summary['stability_screen_passed']=all(x.get('stability_screen_passed',False) for x in summary['runs'])
    summary['comparison_passed']=all(x['passed'] for x in summary['comparisons'])
    (run_dir/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print('SUMMARY '+str(run_dir/'summary.json'),flush=True)
    return 0 if summary['algebra_passed'] and summary['comparison_passed'] and summary['stability_screen_passed'] else 1
if __name__=='__main__':
    raise SystemExit(main())
