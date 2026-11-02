"""Run the gated acceptance queue sequentially and keep an inspectable report.

--wait-small-pid attaches to an already-running small suite without launching
another solver. A failed gate stops the queue; this script never changes defaults.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time

import psutil
from validate_gmres_openmp import ROOT, save, sha


def report(root, state):
    lines = ["# GMRES OpenMP acceptance", "", f"Status: **{state['status']}**", "",
             "All solver runs use 40 threads and run sequentially. Default CPU mode remains legacy.", "",
             f"Executable SHA256: `{state['exe_sha256']}`", "",
             "| Suite / case | Numerical gate | GMRES speedup | Compute speedup | Total speedup |",
             "|---|---|---:|---:|---:|"]
    for suite in ["small", "large-short", "large-full"]:
        file = root / (suite + ".json")
        if not file.exists():
            continue
        data = json.loads(file.read_text(encoding="utf-8"))
        for case, perf in data.get("performance", {}).items():
            checks = [x for x in data["comparisons"] if x["case"] == case]
            passed = bool(checks) and all(x["passed"] for x in checks)
            ratios = [perf[name]["speedup"] for name in ["multidomain.gmres", "run.solve", "run.total"]]
            lines.append(f"| {suite} / {case} | {passed} | " + " | ".join(f"{x:.3f}" for x in ratios) + " |")
    if state.get("error"):
        lines += ["", "Queue stopped: " + state["error"]]
    lines += ["", "Detailed evidence: per-suite JSON, per-run result.json/input_manifest.json, run.log,",
              "timing.txt, gmres_cpu_validation.csv and full validation_state.csv in this directory.",
              "Large-case historical physical-accuracy warnings are not resolved by a backend comparison."]
    (root / "acceptance_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "tmp/gmres_openmp/acceptance")
    parser.add_argument("--exe", type=Path, default=ROOT / "tmp/gmres_openmp/build/DBEM1_openmp.exe")
    parser.add_argument("--wait-small-pid", type=int)
    args = parser.parse_args()
    args.output = args.output.resolve(); args.exe = args.exe.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    state = {"status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
             "exe_sha256": sha(args.exe), "threads": 40, "current_suite": "small", "completed_suites": []}

    def update():
        save(args.output / "acceptance_queue.json", state)
        report(args.output, state)

    def gate(suite):
        file = args.output / (suite + ".json")
        data = json.loads(file.read_text(encoding="utf-8"))
        if data["exe_sha256"] != state["exe_sha256"] or not data["passed"]:
            raise RuntimeError(f"{suite}: numerical gate did not pass for this executable")

    update()
    try:
        if args.wait_small_pid:
            process = psutil.Process(args.wait_small_pid)
            if "validate_gmres_openmp.py" not in " ".join(process.cmdline()):
                raise RuntimeError("Attach PID is not the small validation runner")
            while process.is_running():
                update()
                try:
                    process.wait(timeout=10)
                except psutil.TimeoutExpired:
                    pass
            gate("small")
            state["completed_suites"].append("small")
        for suite in (["large-short", "large-full", "diagnostic"] if args.wait_small_pid else ["small", "large-short", "large-full", "diagnostic"]):
            state["current_suite"] = suite
            update()
            command = [sys.executable, str(ROOT / "tools/validate_gmres_openmp.py"), "--suite", suite,
                       "--exe", str(args.exe), "--output", str(args.output)]
            with (args.output / (suite + "_queue.log")).open("ab") as log:
                child = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.DEVNULL,
                    stdout=log, stderr=subprocess.STDOUT, creationflags=subprocess.CREATE_NO_WINDOW)
                state["runner_pid"] = child.pid
                while child.poll() is None:
                    update()
                    time.sleep(10)
            if child.returncode:
                raise RuntimeError(f"{suite}: runner exited {child.returncode}; see {suite}_queue.log")
            if suite != "diagnostic":
                gate(suite)
            state["completed_suites"].append(suite)
            update()
        state["status"] = "completed"
        state["performance_gate_passed"] = json.loads((args.output / "large-short.json").read_text(encoding="utf-8"))["performance_gate_passed"]
    except Exception as error:
        state["status"] = "stopped_at_failed_gate"
        state["error"] = str(error)
        raise
    finally:
        state["updated_utc"] = datetime.now(timezone.utc).isoformat()
        update()


if __name__ == "__main__":
    main()
