"""Retained legacy entry: finite/residual and physical-growth checks, kept separate."""
from pathlib import Path
import json,time,csv,math,sys,argparse
sys.path.insert(0,str(Path.cwd()/"DBEM1/tests"))
import run_single_domain_stability as a
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--case",action="append",choices=("legacy_34_n480","legacy_136_n960","legacy_nu025_short","legacy_nu045_short"),help="Run only the named fixture; may be repeated")
parser.add_argument("--exe",type=Path,default=Path.cwd()/"validation_runs/e08_nan_diagnosis_20261019/after/base/x64/Release/base.exe")
args=parser.parse_args()
work=Path.cwd()/"validation_runs/legacy_completion_20261019"/("local_"+time.strftime("%Y%m%d_%H%M%S"))
work.mkdir(parents=True)
exe=args.exe.resolve()
assert exe.is_file(),"Build first: python.exe DBEM1/tests/singular_initialization_race.py --label after"
original=a.v.run_case
def legacy_run(case,*args,**kwargs):
    card=case/"BEM_DATACARD.DAT"
    lines=card.read_text().splitlines()
    assert len(lines)==19
    card.write_text("\n".join(lines[:16])+"\n")
    run_args=list(args)
    if len(run_args)>1:run_args[1]=1800
    else:kwargs["timeout"]=1800
    return original(case,*run_args,**kwargs)
a.v.run_case=legacy_run
records=[]
for config in [
 dict(name="legacy_34_n480",threads=4,tau=24.),
 dict(name="legacy_136_n960",threads=4,mesh=(16,2,2),tau=24.),
 dict(name="legacy_nu025_short",threads=4,nu=.25,tau=1.),
 dict(name="legacy_nu045_short",threads=4,nu=.45,tau=1.),
]:
    if args.case and config["name"] not in args.case:continue
    record,rows=a.run_one(work,exe,config)
    # The reused metrics parser targets multidomain; replace its route-dependent verdict.
    output=Path(record["output"])
    with (output/"residual_probe.csv").open() as fp:
        residual=list(csv.DictReader(fp))
    assert all(r["solver"]=="single" for r in residual)
    assert all(math.isfinite(float(r[k])) for r in residual for k in ("relativeResidual","maxAbsResidual","l2Residual","l2Rhs"))
    record["framework"]="retained legacy entry (no domain configuration)"
    for key in ("flag_count","nonzero_flags","solved_all_steps","matrix_structure_valid","tree_mode","assembly_max_difference","max_inner_iterations"):
        record.pop(key,None)
    record["passed_algebra"]=(record["returncode"]==1 and record["finite_state"]
        and record["residual_all_steps"] and record["max_relative_residual"]<=1.01e-8
        and record["known_boundary_max_error"]<1e-12)
    record["stability_screen_passed"]=(record["passed_algebra"] and record["max_u_over_unit_static"]<=10.)
    records.append(record)
    (work/"summary.json").write_text(json.dumps(records,indent=2))
    print(json.dumps(record),flush=True)
    assert record["passed_algebra"],record
print(work/"summary.json")
