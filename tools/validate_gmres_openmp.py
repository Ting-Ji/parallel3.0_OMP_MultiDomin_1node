"""Sequential, resumable 40-thread physical validation of legacy/OMP GMRES.

Run with Windows Python. All inputs/results are isolated under --output.
Suites are gated: small -> large-short (one pair/case) -> large-full.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import numpy as np
import pandas as pd
import psutil
import run_multidomain_validation as rod

ROOT = Path(__file__).resolve().parents[1]
CUBES = ROOT / "paper/E15/cube_xgrowth_1000_1500_2000_20261024/inputs"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def save(path, data):
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    temp.replace(path)


def prepare(directory, case, steps):
    run = directory / "DBEM1"
    run.mkdir(parents=True)
    (run / "input").mkdir()
    (directory / "output").mkdir()
    if case.startswith("rod"):
        count = int(case[3:])
        title = case
        spec = rod.CaseSpec(title, count, steps, True, beta_override=.8)
        mesh = rod.generate_multidomain_mesh(count)
        rod.write_model(run / "input" / (title + ".model"), mesh)
        rod.write_bd(run / "input" / (title + ".bd"), mesh, steps, count, True)
        rod.write_datacard(run / "BEM_DATACARD.DAT", title, spec, mesh)
        card = (run / "BEM_DATACARD.DAT").read_text().splitlines()
        origin = None
    else:
        source = CUBES / case if case.startswith("cube") else ROOT / "DBEM1/input"
        origin = source / (case + ".DATACARD")
        card = origin.read_text().splitlines()
        title = card[0]
        for extension in ["model", "bd"]:
            shutil.copy2(source / (case + "." + extension), run / "input" / (title + "." + extension))
        if case.startswith("tsv"):
            card[6:9] = ["110", "0.34", "8.96"]
    card[1] = "40"
    card[4] = str(steps)
    (run / "BEM_DATACARD.DAT").write_text("\n".join(card) + "\n", encoding="ascii")
    for name in ["gaussposition", "gaussweight"]:
        shutil.copy2(ROOT / "DBEM1" / name, run / name)
    manifest = {"case": case, "steps": steps, "threads": 40, "tolerance": float(card[12]),
                "material": card[6:9], "original_card_sha256": sha(origin) if origin else None,
                "input_sha256": {str(p.relative_to(run)): sha(p) for p in [run / "BEM_DATACARD.DAT", *sorted((run / "input").iterdir())]}}
    save(directory / "input_manifest.json", manifest)
    return run, manifest


def execute(args, case, steps, mode, repeat=0, profile=False):
    directory = args.output / case / (f"{mode}_{steps}_{repeat}" + ("_profile" if profile else ""))
    result_file = directory / "result.json"
    executable_hash = sha(args.exe)
    if result_file.exists():
        result = json.loads(result_file.read_text(encoding="utf-8"))
        if result["exe_sha256"] != executable_hash:
            raise RuntimeError(f"Binary changed; select another output directory: {directory}")
        if not result["success"]:
            raise RuntimeError(f"Existing failed result: {directory}")
        return result
    if directory.exists():
        raise RuntimeError(f"Inspect unfinished run before retry: {directory}")
    run, manifest = prepare(directory, case, steps)
    env = dict(os.environ)
    env.update(DBEM_GMRES_CPU_MODE=mode, DBEM_GMRES_CPU_PROFILE="1" if profile else "0", OMP_NUM_THREADS="40", OMP_DYNAMIC="FALSE")
    print(f"START {case} {steps} {mode} repeat={repeat}", flush=True)
    begin = time.perf_counter()
    peak = 0
    with (directory / "run.log").open("wb") as log:
        child = subprocess.Popen([str(args.exe)], cwd=run, env=env, stdin=subprocess.DEVNULL,
                                 stdout=log, stderr=subprocess.STDOUT, creationflags=subprocess.CREATE_NO_WINDOW)
        monitor = psutil.Process(child.pid)
        save(directory / "running.json", {"pid": child.pid, "start_unix": time.time(), "exe_sha256": executable_hash})
        while child.poll() is None:
            try:
                peak = max(peak, monitor.memory_info().rss)
            except psutil.NoSuchProcess:
                pass
            if time.perf_counter() - begin > args.timeout:
                child.kill()
                child.wait()
                raise TimeoutError(f"Timed out: {directory}")
            time.sleep(1)
    elapsed = time.perf_counter() - begin
    log = (directory / "run.log").read_bytes().decode("mbcs", errors="replace")
    records = [tuple(map(int, m)) for m in re.findall(r"MultiDomain GMRES step (\d+): flag=(-?\d+) outer=(\d+) inner=(\d+)", log)]
    timing_files = list((run / "output").rglob("timing.txt"))
    timings = {}
    if len(timing_files) == 1:
        with timing_files[0].open() as f:
            for row in csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t"):
                timings[row["stage"]] = float(row["cumulative_seconds"])
    state_files = list((run / "output").rglob("validation_state.csv"))
    checks = list((run / "output").rglob("gmres_cpu_validation.csv"))
    residuals = []
    if len(checks) == 1:
        with checks[0].open() as f:
            residuals = list(csv.DictReader(f))
    residual_ok = len(residuals) == steps and all(int(r["flag"]) == 0 and math.isfinite(float(r["relative_residual"])) and
        float(r["relative_residual"]) <= float(r["tolerance"]) for r in residuals)
    success = child.returncode in [0, 1] and len(records) == steps and all(r[1] == 0 for r in records) and len(state_files) == 1 and residual_ok
    result = {"case": case, "steps": steps, "mode": mode, "repeat": repeat, "success": success,
              "exe_sha256": executable_hash, "returncode": child.returncode, "elapsed_s": elapsed,
              "peak_sampled_rss_GiB": peak/2**30, "stage_timings": timings,
              "solved_steps": len(records), "inner_iterations": sum(r[3] for r in records),
              "max_relative_residual": max((float(r["relative_residual"]) for r in residuals if math.isfinite(float(r["relative_residual"]))), default=None),
              "state_file": str(state_files[0]) if state_files else None, "manifest": manifest}
    save(result_file, result)
    print(f"END {case} {mode} success={success} elapsed={elapsed:.1f}s gmres={timings.get('multidomain.gmres')}", flush=True)
    if not success:
        raise RuntimeError(f"Numerical/execution gate failed: {directory}")
    return result


def compare(a, b):
    if a["manifest"]["input_sha256"] != b["manifest"]["input_sha256"]:
        raise RuntimeError("Inputs differ between modes")
    fields = ["ux", "uy", "uz", "tx", "ty", "tz"]
    keys = ["step", "element", "localNode", "nodeId"]
    arrays = [pd.read_csv(r["state_file"], usecols=keys+fields, chunksize=200000) for r in [a,b]]
    sums = np.zeros((2,4)); rows = 0
    while True:
        ca, cb = next(arrays[0], None), next(arrays[1], None)
        if ca is None or cb is None:
            if ca is not None or cb is not None:
                raise RuntimeError("State row counts differ")
            break
        if not np.array_equal(ca[keys].to_numpy(), cb[keys].to_numpy()):
            raise RuntimeError("State keys differ")
        va, vb = ca[fields].to_numpy(), cb[fields].to_numpy()
        if not np.isfinite(va).all() or not np.isfinite(vb).all():
            raise RuntimeError("Nonfinite physical state")
        for g in range(2):
            ref = va[:,3*g:3*g+3]; diff = vb[:,3*g:3*g+3] - ref
            sums[g,0] += np.sum(diff*diff); sums[g,1] += np.sum(ref*ref)
            sums[g,2] = max(sums[g,2], np.max(np.abs(diff)))
            sums[g,3] = max(sums[g,3], np.max(np.abs(ref)))
        rows += len(ca)
    groups = {}
    for name, values in zip(["displacement", "traction"], sums):
        d2, r2, maximum, peak = map(float, values)
        relative = math.sqrt(d2)/max(math.sqrt(r2),1e-30)
        groups[name] = {"relative_L2": relative, "max_abs_difference": maximum, "reference_peak": peak,
                        "passed": relative <= 1e-6 and maximum <= 1e-5*max(peak,1e-30)}
    return {"passed": all(g["passed"] for g in groups.values()), "rows": rows, "groups": groups,
            "inner_iteration_difference": b["inner_iterations"]-a["inner_iterations"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=["small", "large-short", "large-full", "diagnostic"], required=True)
    parser.add_argument("--exe", type=Path, default=ROOT / "tmp/gmres_openmp/build/DBEM1_openmp.exe")
    parser.add_argument("--output", type=Path, default=ROOT / "tmp/gmres_openmp/validation")
    parser.add_argument("--timeout", type=float, default=21600)
    args = parser.parse_args()
    args.exe = args.exe.resolve(); args.output = args.output.resolve(); args.output.mkdir(parents=True, exist_ok=True)
    if args.suite == "small":
        cases = [("rod1",24),("rod2",24)] + [(f"tsv_c{i:02d}_md{2*i:02d}_orig_b060_n1000",1000) for i in [1,4,9]]
        repeats = 1
    else:
        prerequisite = "small" if args.suite == "large-short" else "large-short"
        if args.suite != "diagnostic":
            prior = json.loads((args.output / (prerequisite + ".json")).read_text(encoding="utf-8"))
            if not prior["passed"] or prior["exe_sha256"] != sha(args.exe):
                raise RuntimeError("Previous numerical gate not passed for this executable")
        cases = [(f"cube_md{n}_x{n//100}_face16", 1000 if args.suite == "large-full" else 100) for n in ([1000] if args.suite in ["large-full","diagnostic"] else [1000,2000])]
        repeats = 1
    summary = {"suite": args.suite, "exe_sha256": sha(args.exe), "sampling": "single_run_per_backend",
               "passed": False, "comparisons": [], "performance": {}}
    for case, steps in cases:
        pairs = []
        for repeat in range(repeats):
            if args.suite == "diagnostic":
                execute(args,case,steps,"omp",profile=True)
                continue
            order = ["legacy","omp"] if repeat % 2 == 0 else ["omp","legacy"]
            pair = {mode: execute(args,case,steps,mode,repeat) for mode in order}
            check = compare(pair["legacy"], pair["omp"])
            summary["comparisons"].append({"case": case,"repeat": repeat, **check})
            save(args.output / (args.suite + ".json"), summary)
            if not check["passed"]:
                raise RuntimeError(f"Field comparison failed: {case}; evidence retained")
            pairs.append(pair)
        if pairs:
            performance = {}
            for stage in ["multidomain.gmres","run.solve","run.total"]:
                serial = pairs[0]["legacy"]["stage_timings"][stage]
                parallel = pairs[0]["omp"]["stage_timings"][stage]
                performance[stage] = {"legacy_s": serial,"omp_s": parallel,"samples_per_backend": 1,"speedup": serial/parallel}
            summary["performance"][case] = performance
            save(args.output / (args.suite + ".json"), summary)
    summary["passed"] = True
    if args.suite == "large-short":
        summary["performance_gate_passed"] = all(p["multidomain.gmres"]["speedup"] >= 1.5 for p in summary["performance"].values())
    save(args.output / (args.suite + ".json"), summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
