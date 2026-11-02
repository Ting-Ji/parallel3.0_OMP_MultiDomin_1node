"""One sequential 2000-domain long-time comparison, both using 40-thread OMP GMRES.

Run the validated GMRES-only executable, then the validated RHS-parallel one.
No solver rebuild, repeat sampling or default-backend change is performed.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import psutil
import validate_gmres_openmp as validation

ROOT = validation.ROOT
CASE = "cube_md2000_x20_face16"
VERSIONS = [
    ("gmres_only", "legacy", ROOT / "tmp/gmres_openmp/build_v2/DBEM1_openmp.exe",
     "41d545684ffe41d326848ffcdfec559d4467426bb07248acf7e784cdcf6d9611"),
    ("gmres_rhs", "parallel", ROOT / "tmp/rhs_parallel/build/DBEM1_openmp.exe",
     "cd72bccf91ddbe5ad43f1c8b4773afac36327b58237039c3c9e1d8f673e3843b"),
]
STAGES = ["multidomain.rhs", "multidomain.gmres", "run.solve", "run.total"]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def report(output, state):
    steps = state["steps"]
    lines = [f"# 2000-domain / {steps}-step GMRES and RHS comparison", "",
             f"Status: **{state['status']}**", "",
             f"Current version: {state.get('current_version', '-')}", "",
             "Geometry: 20 x 10 x 10; partition 20 x 10 x 10; 192000 elements; 576000 scalar unknowns.",
             f"{steps} steps; dt=0.03 ms; {steps * .03:g} ms total; 40 threads; sequential single runs.",
             "Both versions: GMRES=omp, profiles disabled, full state CSV and independent true residual audit.",
             "Original material E=200, nu=0, rho=8 and original mesh/load/tolerance retained.", "",
             "| Version | RHS seconds | GMRES seconds | Compute seconds | Total seconds | Inner iterations | Sampled RSS GiB |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for name, result in state.get("results", {}).items():
        values = [result["stage_timings"][stage] for stage in STAGES]
        lines.append("| " + name + " | " + " | ".join(f"{x:.3f}" for x in values) +
                     f" | {result['inner_iterations']} | {result['peak_sampled_rss_GiB']:.3f} |")
    if state.get("memory_available_GiB") is not None:
        lines += ["", f"Prelaunch memory available: {state['memory_available_GiB']:.2f} GiB; "
                  f"minimum: {state['minimum_available_GiB']:.2f} GiB."]
    if state.get("comparison"):
        check = state["comparison"]
        lines += ["", f"Full-field gate: **{check['passed']}**; compared rows: {check['rows']}; "
                  f"inner iteration delta: {check['inner_iteration_difference']}.", "",
                  "| Field | Relative L2 | Max absolute difference | Passed |", "|---|---:|---:|---|"]
        for field, values in check["groups"].items():
            lines.append(f"| {field} | {values['relative_L2']:.12g} | {values['max_abs_difference']:.12g} | {values['passed']} |")
    if state.get("speedups"):
        lines += ["", "Speedup = GMRES-only time / GMRES+RHS time:"]
        for stage, value in state["speedups"].items():
            lines += [f"- {stage}: {value:.4f}x"]
    if state.get("error"):
        lines += ["", "Stopped: " + state["error"]]
    lines += ["", f"One run/version; no averages/medians, no old baseline rerun beyond this explicitly requested new {steps}-step pair.",
              "Single-run timings do not exclude machine-state variation. Historical physical accuracy warnings remain separate.",
              "Default backends unchanged. Raw inputs, logs, result.json, timing, residual and full-field CSV are retained."]
    (output / "comparison_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def preflight(output, steps):
    source = validation.CUBES / CASE
    metadata = read(source / (CASE + ".json"))
    assert metadata["domains"] == 2000 and metadata["size"] == [20, 10, 10]
    assert metadata["elements"] == 192000 and metadata["nstep"] >= steps
    for suffix in ["model", "bd"]:
        if validation.sha(source / (CASE + "." + suffix)) != metadata[suffix + "_sha256"]:
            raise RuntimeError("Original input hash mismatch: " + suffix)
    card = (source / (CASE + ".DATACARD")).read_text().splitlines()
    if card[0] != CASE or int(card[4]) < steps or list(map(float, card[6:9])) != [200, 0, 8]:
        raise RuntimeError("Unexpected source card")
    builds = []
    for name, mode, executable, digest in VERSIONS:
        if validation.sha(executable) != digest:
            raise RuntimeError("Validated executable changed: " + name)
        manifest = read(executable.parent / "build_manifest.json")
        if manifest["executables"][executable.name] != digest:
            raise RuntimeError("Executable manifest mismatch: " + name)
        builds.append(manifest)
    if builds[0]["flags"] != builds[1]["flags"]:
        raise RuntimeError("Build flags differ")
    before, after = [{k.replace("\\", "/"): v for k, v in m["source_sha256"].items()} for m in builds]
    if any(after.get(k) != v for k, v in before.items() if k != "DBEM1/multidomain.cpp"):
        raise RuntimeError("Unexpected non-RHS source differences between validated builds")
    # Only a 100-step, 40-thread pair exists locally; old 64-thread historical
    # 1000-step outputs do not provide a matching performance baseline.
    evidence = {"case": CASE, "steps": steps, "threads": 40, "metadata": metadata,
                "source_card_sha256": validation.sha(source / (CASE + ".DATACARD")),
                "source_card_changes": {"threads": [int(card[1]), 40], "steps": [int(card[4]), steps]},
                "build_flags_match": True, "non_rhs_sources_match": True,
                "versions": {name: {"exe": str(exe), "sha256": digest, "gmres_mode": "omp", "rhs_mode": mode}
                             for name, mode, exe, digest in VERSIONS}}
    validation.save(output / "comparison_manifest.json", evidence)


def wait_for_memory(args, state, update):
    while True:
        available = psutil.virtual_memory().available / 2**30
        state["memory_available_GiB"] = available
        state["minimum_available_GiB"] = args.min_available_gib
        if available >= args.min_available_gib:
            state["status"] = "running"; update(); return
        state["status"] = "waiting_for_memory"; update()
        time.sleep(30)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, choices=[800, 1000], default=800)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--timeout", type=float, default=21600, help="maximum seconds per solver")
    parser.add_argument("--min-available-gib", type=float, help="default: 37 GiB for 800 steps, 41 for 1000")
    parser.add_argument("--detach", action="store_true")
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    args.output = (args.output or ROOT / f"tmp/rhs_parallel/cube2000_{args.steps}").resolve()
    if args.min_available_gib is None:
        args.min_available_gib = 37.0 if args.steps == 800 else 41.0
    if args.min_available_gib < 0:
        raise ValueError("Minimum memory cannot be negative")
    args.output.mkdir(parents=True, exist_ok=True)
    if args.detach:
        if (args.output / "queue.lock").exists():
            raise RuntimeError("Queue lock exists; inspect the running/unfinished queue")
        with (args.output / "comparison_queue.log").open("ab") as log:
            child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
                                      *[x for x in sys.argv[1:] if x != "--detach"]], cwd=ROOT,
                stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NO_WINDOW | subprocess.CREATE_NEW_PROCESS_GROUP)
        print(json.dumps({"queue_pid": child.pid, "output": str(args.output)}, ensure_ascii=True)); return
    state = {"status": "auditing", "case": CASE, "steps": args.steps, "threads": 40,
             "queue_pid": os.getpid(), "queue_created_unix": psutil.Process().create_time(),
             "started_utc": datetime.now(timezone.utc).isoformat(), "results": {}, "completed_versions": []}
    def update():
        state["updated_utc"] = datetime.now(timezone.utc).isoformat()
        validation.save(args.output / "comparison_queue.json", state); report(args.output, state)
    lock = args.output / "queue.lock"
    with lock.open("x") as stream:
        json.dump({"pid": os.getpid(), "created_unix": state["queue_created_unix"]}, stream)
    try:
        update(); preflight(args.output, args.steps)
        if args.audit_only:
            state["status"] = "audit_passed"; update(); return
        if psutil.disk_usage(str(args.output)).free < 100 * 2**30:
            raise RuntimeError("Less than 100 GiB free for both complete state exports")
        for name, rhs_mode, executable, digest in VERSIONS:
            state["current_version"] = name; update()
            path = args.output / name / CASE / f"omp_{args.steps}_0" / "result.json"
            if not path.exists():
                wait_for_memory(args, state, update)
            os.environ.update(DBEM_RHS_CPU_MODE=rhs_mode, DBEM_RHS_CPU_PROFILE="0")
            run_args = SimpleNamespace(output=args.output / name, exe=executable, timeout=args.timeout)
            result = validation.execute(run_args, CASE, args.steps, "omp", repeat=0, profile=False)
            if name == "gmres_rhs":
                log = (path.parent / "run.log").read_bytes().decode("mbcs", errors="replace")
                if "RHS CPU mode=parallel threads=40" not in log:
                    result["success"] = False; validation.save(path, result)
                    raise RuntimeError("40-thread RHS parallel path not confirmed")
            result.update(gmres_mode="omp", rhs_mode=rhs_mode, profiles_enabled=False)
            validation.save(path, result)
            state["results"][name] = result; state["completed_versions"].append(name); update()
        state["status"] = "comparing_full_fields"; update()
        old, new = state["results"]["gmres_only"], state["results"]["gmres_rhs"]
        check = validation.compare(old, new)
        if check["rows"] != (args.steps + 1) * 192000:
            raise RuntimeError("Unexpected complete state row count")
        state["comparison"] = check
        state["speedups"] = {stage: old["stage_timings"][stage] / new["stage_timings"][stage] for stage in STAGES}
        update()
        if not check["passed"]:
            raise RuntimeError("Full-field gate failed; evidence retained; no retry")
        state["status"] = "completed"; update()
    except BaseException as error:
        state["status"] = "stopped_at_failed_gate"; state["error"] = str(error); update(); raise
    finally:
        lock.unlink()


if __name__ == "__main__":
    main()
