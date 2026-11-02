"""Gated 40-thread RHS acceptance against read-only completed GMRES OMP runs.

Only the new RHS parallel solver is launched, once per case/step combination.
Incomplete/failed runs are retained and never automatically restarted.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

import psutil
import validate_gmres_openmp as validation

ROOT = validation.ROOT
CASES = {
    "small": [("rod1", 24), ("rod2", 24)] +
             [(f"tsv_c{i:02d}_md{2*i:02d}_orig_b060_n1000", 1000) for i in [1, 4, 9]],
    "large-short": [(f"cube_md{n}_x{n//100}_face16", 100) for n in [1000, 2000]],
    "large-full": [("cube_md1000_x10_face16", 1000)],
}
STAGES = ["multidomain.rhs", "multidomain.gmres", "run.solve", "run.total"]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def normalized_sources(manifest):
    return {name.replace("\\", "/"): digest for name, digest in manifest["source_sha256"].items()}


def audit(args):
    """Finish compatibility checks for every old result before launching any solver."""
    old_build = read(args.baseline_build / "build_manifest.json")
    new_build = read(args.exe.parent / "build_manifest.json")
    old_hash = old_build["executables"]["DBEM1_openmp.exe"]
    new_hash = validation.sha(args.exe)
    if new_build["executables"][args.exe.name] != new_hash:
        raise RuntimeError("New executable differs from its build manifest")
    if old_build["flags"] != new_build["flags"]:
        raise RuntimeError("Build flags differ; baseline not compatible; no baseline rerun")
    before, after = normalized_sources(old_build), normalized_sources(new_build)
    unexpected = [name for name, digest in before.items()
                  if name != "DBEM1/multidomain.cpp" and after.get(name) != digest]
    additions = [name for name in after if name not in before and not name.startswith("DBEM1/rhs_cpu_")]
    if unexpected or additions:
        raise RuntimeError(f"Changes outside RHS implementation: {unexpected + additions}")
    evidence = {"passed": False, "baseline_root": str(args.baseline),
                "baseline_exe_sha256": old_hash, "new_exe_sha256": new_hash,
                "build_flags_match": True, "non_rhs_sources_match": True,
                "baseline_gmres_mode": "omp", "baseline_rhs": "original pre-RHS implementation",
                "new_gmres_mode": "omp", "new_rhs_mode": "parallel", "threads": 40,
                "sampling": "one_new_run_per_case_reuse_old_repeat_0", "cases": []}
    baselines = {}
    for suite, cases in CASES.items():
        old_suite = read(args.baseline / (suite + ".json"))
        if not old_suite["passed"] or old_suite["exe_sha256"] != old_hash:
            raise RuntimeError(f"Old {suite} gate unavailable or incompatible; no rerun")
        for case, steps in cases:
            directory = args.baseline / case / f"omp_{steps}_0"
            result = read(directory / "result.json")
            if not (result["success"] and result["mode"] == "omp" and result["repeat"] == 0 and
                    result["steps"] == steps and result["exe_sha256"] == old_hash and
                    result["manifest"]["threads"] == 40 and result["solved_steps"] == steps):
                raise RuntimeError(f"Old run incompatible: {directory}; no rerun")
            proof = args.output / "input_audit" / f"{case}_{steps}"
            if not proof.exists():
                validation.prepare(proof, case, steps)
            expected_manifest = read(proof / "input_manifest.json")
            if expected_manifest != result["manifest"]:
                raise RuntimeError(f"Generated inputs differ from baseline: {case}/{steps}")
            for relative, digest in result["manifest"]["input_sha256"].items():
                if validation.sha(directory / "DBEM1" / relative) != digest:
                    raise RuntimeError(f"Old input evidence changed: {directory}/{relative}")
            state = Path(result["state_file"])
            if not state.is_file() or state.stat().st_size == 0:
                raise RuntimeError(f"Old full state missing: {state}; no rerun")
            checks = list((directory / "DBEM1/output").rglob("gmres_cpu_validation.csv"))
            if len(checks) != 1:
                raise RuntimeError(f"Old residual evidence missing: {directory}")
            with checks[0].open() as stream:
                residuals = list(csv.DictReader(stream))
            if len(residuals) != steps or not all(int(r["flag"]) == 0 and
                    validation.math.isfinite(float(r["relative_residual"])) and
                    float(r["relative_residual"]) <= float(r["tolerance"]) for r in residuals):
                raise RuntimeError(f"Old residual evidence invalid: {directory}")
            if any(result["stage_timings"].get(stage, 0) <= 0 for stage in STAGES):
                raise RuntimeError(f"Old timing evidence missing: {directory}")
            baselines[(case, steps)] = result
            evidence["cases"].append({"case": case, "steps": steps, "suite": suite,
                "baseline_result": str(directory / "result.json"), "input_sha256": expected_manifest["input_sha256"],
                "state_file": str(state), "state_bytes": state.stat().st_size,
                "residual_rows": len(residuals), "compatible": True})
    evidence["passed"] = True
    validation.save(args.output / "baseline_match.json", evidence)
    return baselines, evidence


def report(args, state):
    lines = ["# RHS parallel acceptance", "", f"Status: **{state['status']}**", "",
             f"Current suite/case: {state.get('current_suite', '-')} / {state.get('current_case', '-')}", "",
             "40 threads; GMRES=omp; RHS=parallel; one new run per case. Old GMRES OMP repeat=0 reused read-only.",
             "Cross-round single-run timings; machine-state variation is not excluded. No averages or medians.", "",
             f"New executable SHA256: `{state['exe_sha256']}`", "",
             f"Baseline executable SHA256: `{state.get('baseline_exe_sha256', 'pending audit')}`", "",
             "| Suite / case | Steps | Field gate | RHS speedup | GMRES ratio | Compute speedup | Total speedup | Iteration delta |",
             "|---|---:|---|---:|---:|---:|---:|---:|"]
    for suite in CASES:
        path = args.output / (suite + ".json")
        if not path.exists():
            continue
        data = read(path)
        for entry in data.get("comparisons", []):
            perf = data.get("performance", {}).get(entry["case"], {})
            ratios = [f"{perf[s]['speedup']:.3f}" if s in perf else "pending" for s in STAGES]
            lines.append(f"| {suite} / {entry['case']} | {entry['steps']} | {entry['passed']} | " +
                         " | ".join(ratios) + f" | {entry['inner_iteration_difference']} |")
    if "performance_gate_passed" in state:
        lines += ["", f"Both short large cases RHS >=1.5x target: **{state['performance_gate_passed']}**"]
    if state.get("error"):
        lines += ["", "Stopped: " + state["error"]]
    lines += ["", "Default RHS and GMRES backends remain legacy. All profile flags are disabled for acceptance.",
              "See baseline_match.json, suite JSON, per-run inputs/result/run.log/timing and full state CSV.",
              "workspace_setup time is reported separately and also remains inside run.solve/run.total.",
              "Historical physical accuracy warnings are separate from parallel consistency."]
    (args.output / "acceptance_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=ROOT / "tmp/rhs_parallel/build/DBEM1_openmp.exe")
    parser.add_argument("--output", type=Path, default=ROOT / "tmp/rhs_parallel/acceptance")
    parser.add_argument("--baseline", type=Path, default=ROOT / "tmp/gmres_openmp/acceptance_v2")
    parser.add_argument("--baseline-build", type=Path, default=ROOT / "tmp/gmres_openmp/build_v2")
    parser.add_argument("--timeout", type=float, default=21600)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--detach", action="store_true", help="start one hidden Windows queue and return")
    args = parser.parse_args()
    for name in ["exe", "output", "baseline", "baseline_build"]:
        setattr(args, name, getattr(args, name).resolve())
    if (args.output == args.baseline or args.output.is_relative_to(args.baseline) or
            args.output.is_relative_to(args.baseline_build)):
        raise RuntimeError("Output must be isolated from old evidence")
    args.output.mkdir(parents=True, exist_ok=True)
    if args.detach:
        if (args.output / "queue.lock").exists():
            raise RuntimeError("Queue lock already exists; inspect it instead of launching another queue")
        command = [sys.executable, str(Path(__file__).resolve()), *[x for x in sys.argv[1:] if x != "--detach"]]
        with (args.output / "acceptance_queue.log").open("ab") as log:
            child = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.DEVNULL, stdout=log,
                stderr=subprocess.STDOUT, creationflags=subprocess.CREATE_NO_WINDOW | subprocess.CREATE_NEW_PROCESS_GROUP)
        print(json.dumps({"queue_pid": child.pid, "output": str(args.output)}, ensure_ascii=False), flush=True)
        return
    state = {"status": "auditing", "started_utc": datetime.now(timezone.utc).isoformat(),
             "exe_sha256": validation.sha(args.exe), "queue_pid": os.getpid(),
             "queue_created_unix": psutil.Process().create_time(), "threads": 40,
             "completed_suites": [], "completed_cases": 0, "total_cases": 8}

    def update():
        state["updated_utc"] = datetime.now(timezone.utc).isoformat()
        validation.save(args.output / "acceptance_queue.json", state)
        report(args, state)

    lock = args.output / "queue.lock"
    # Deliberately exclusive: a crash leaves the lock/evidence for inspection.
    with lock.open("x") as stream:
        json.dump({"pid": os.getpid(), "created_unix": state["queue_created_unix"]}, stream)
    try:
        update()
        baselines, evidence = audit(args)
        state["baseline_exe_sha256"] = evidence["baseline_exe_sha256"]
        if args.audit_only:
            state["status"] = "audit_passed"; update(); return
        state["status"] = "running"; update()
        os.environ.update(DBEM_RHS_CPU_MODE="parallel", DBEM_RHS_CPU_PROFILE="0")
        for suite, cases in CASES.items():
            state["current_suite"] = suite
            path = args.output / (suite + ".json")
            summary = read(path) if path.exists() else {"suite": suite, "passed": False,
                "exe_sha256": state["exe_sha256"], "baseline_exe_sha256": evidence["baseline_exe_sha256"],
                "sampling": "one_new_run_reuse_old_baseline", "comparisons": [], "performance": {}}
            if summary["exe_sha256"] != state["exe_sha256"]:
                raise RuntimeError("New executable changed; choose another output directory")
            for case, steps in cases:
                state["current_case"] = case; state["current_steps"] = steps; update()
                print(f"CASE {suite} {case} steps={steps}; reuse old omp baseline", flush=True)
                old = baselines[(case, steps)]
                new = validation.execute(args, case, steps, "omp", repeat=0, profile=False)
                run = args.output / case / f"omp_{steps}_0"
                log = (run / "run.log").read_bytes().decode("mbcs", errors="replace")
                if "RHS CPU mode=parallel threads=40" not in log:
                    new["success"] = False; validation.save(run / "result.json", new)
                    raise RuntimeError(f"Parallel RHS 40-thread execution not confirmed: {run}")
                new.update(rhs_mode="parallel", rhs_profile=False, gmres_mode="omp")
                validation.save(run / "result.json", new)
                check = validation.compare(old, new)
                entry = {"case": case, "steps": steps, **check,
                         "baseline_result": str(args.baseline / case / f"omp_{steps}_0/result.json"),
                         "new_result": str(run / "result.json"),
                         "workspace_setup_s": new["stage_timings"].get("rhs.cpu.workspace_setup", 0),
                         "baseline_peak_sampled_rss_GiB": old["peak_sampled_rss_GiB"],
                         "new_peak_sampled_rss_GiB": new["peak_sampled_rss_GiB"]}
                summary["comparisons"] = [x for x in summary["comparisons"] if x["case"] != case] + [entry]
                summary["performance"][case] = {
                    stage: {"baseline_s": old["stage_timings"][stage], "parallel_s": new["stage_timings"][stage],
                            "speedup": old["stage_timings"][stage] / new["stage_timings"][stage], "samples_per_version": 1}
                    for stage in STAGES}
                validation.save(path, summary); update()
                if not check["passed"]:
                    raise RuntimeError(f"Field gate failed: {case}; evidence retained; no retry")
                state["completed_cases"] += 1; update()
            summary["passed"] = True
            if suite == "large-short":
                summary["performance_gate_passed"] = all(
                    values["multidomain.rhs"]["speedup"] >= 1.5 for values in summary["performance"].values())
                state["performance_gate_passed"] = summary["performance_gate_passed"]
            validation.save(path, summary)
            state["completed_suites"].append(suite); update()
        state["status"] = "completed"; update()
    except BaseException as error:
        state["status"] = "stopped_at_failed_gate"; state["error"] = str(error); update(); raise
    finally:
        lock.unlink()


if __name__ == "__main__":
    main()
