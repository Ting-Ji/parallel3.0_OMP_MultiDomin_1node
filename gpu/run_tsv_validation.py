"""Full-length, isolated CPU/CUDA comparisons of existing TSV models.

Run with Windows Python (may be launched from WSL). Copied cards use the case's
copper global material and the requested thread count (default 40). Domain
materials, boundary history, 1000 steps and solver tolerance are preserved.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import psutil

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    temporary.replace(path)


def prepare(case, backend, output, threads):
    source = ROOT / "DBEM1/input"
    card = source / (case + ".DATACARD")
    config = json.loads((source / (case + ".config.json")).read_text())
    lines = card.read_text().splitlines()
    assert lines[0] == case and int(lines[4]) == 1000
    copper = config["copper"]
    modified = lines.copy()
    modified[1] = str(threads)
    modified[6:9] = [format(copper[key], ".17g") for key in ["E", "nu", "rho"]]
    assert modified[0] == lines[0] and modified[2:6] == lines[2:6] and modified[9:] == lines[9:]
    run = output / case / backend / "DBEM1"
    if run.exists():
        raise RuntimeError(f"Unfinished run already exists; inspect before retrying: {run}")
    (run / "input").mkdir(parents=True)
    (run.parent / "output").mkdir()
    (run / "BEM_DATACARD.DAT").write_text("\n".join(modified) + "\n")
    files = {}
    for extension in ["model", "bd"]:
        name = case + "." + extension
        shutil.copy2(source / name, run / "input" / name)
        files[name] = sha256(source / name)
        assert files[name] == sha256(run / "input" / name)
    for name in ["gaussposition", "gaussweight"]:
        shutil.copy2(ROOT / "DBEM1" / name, run / name)
    manifest = dict(case=case, backend=backend, original_card_sha256=sha256(card),
                    copied_card_sha256=sha256(run / "BEM_DATACARD.DAT"),
                    original_global_material=lines[6:9], global_material=copper,
                    unchanged_card_from_line_10=True, input_sha256=files,
                    threads=int(modified[1]), steps=int(modified[4]),
                    iterations=int(modified[11]), tolerance=float(modified[12]))
    save(run.parent / "input_manifest.json", manifest)
    return run, manifest


def parse_run(run, elapsed, returncode, peak_rss, samples, timeout):
    log = (run.parent / "run.log").read_bytes().decode("mbcs", errors="replace")
    steps = [tuple(map(int, row)) for row in re.findall(
        r"MultiDomain GMRES step (\d+): flag=(-?\d+) outer=(\d+) inner=(\d+)", log)]
    timings = {}
    timing_files = list((run / "output").rglob("timing.txt"))
    if len(timing_files) == 1:
        for row in timing_files[0].read_text().splitlines():
            values = row.split("\t")
            if len(values) == 5 and values[0] != "stage":
                timings[values[0]] = dict(calls=int(values[3]), seconds=float(values[4]))
    state_files = list((run / "output").rglob("validation_state.csv"))
    performance = re.search(r"GlobalUnknownBlocks=(\d+).*MTS0NonZeroBlocks=(\d+).*MatrixAssemblyTime=([\d.eE+-]+) GMRESTime=([\d.eE+-]+)", log)
    setup = re.findall(r"CUDA GMRES setup: n=(\d+) nnz=(\d+) pre_nnz=(\d+) seconds=([\d.eE+-]+)", log)
    gpu = re.findall(r"CUDA GMRES: flag=(-?\d+) outer=(\d+) inner=(\d+) residual_gpu=([\d.eE+-]+) residual_host=([\d.eE+-]+) solve_seconds=([\d.eE+-]+) device_MiB=([\d.eE+-]+)", log)
    success = returncode in [0, 1] and len(steps) == 1000 and all(s[1] == 0 for s in steps) and len(state_files) == 1
    result = dict(success=success, returncode=returncode, timeout=timeout,
                  process_seconds=elapsed, peak_sampled_process_rss_MiB=peak_rss / 2**20,
                  solved_steps=len(steps), failed_steps=[s for s in steps if s[1]],
                  inner_iteration_total=sum(s[3] for s in steps),
                  inner_iteration_max=max((s[3] for s in steps), default=0),
                  stage_timings=timings, setup_count=len(setup),
                  state_file=str(state_files[0]) if state_files else None,
                  gpu_system_memory_samples=samples)
    if performance:
        result.update(unknown_blocks=int(performance[1]), scalar_dimension=3 * int(performance[1]),
                      matrix_nonzero_blocks=int(performance[2]), assembly_seconds=float(performance[3]),
                      gmres_seconds=float(performance[4]))
    if setup:
        result["gpu_setup"] = dict(n=int(setup[0][0]), nnz=int(setup[0][1]),
                                   pre_nnz=int(setup[0][2]), seconds=float(setup[0][3]))
    if gpu:
        result.update(gpu_solve_count=len(gpu), gpu_host_residual_max=max(float(s[4]) for s in gpu),
                      gpu_internal_solve_seconds=sum(float(s[5]) for s in gpu),
                      gpu_explicit_device_MiB=max(float(s[6]) for s in gpu))
        result["success"] &= len(gpu) == 1000 and len(setup) == 1
    if not result["success"]:
        result["log_tail"] = log[-7000:]
    return result


def run_one(case, backend, output, timeout, threads):
    result_path = output / case / backend / "result.json"
    if result_path.exists():
        manifest = json.loads(result_path.with_name("input_manifest.json").read_text(encoding="utf-8"))
        if manifest["threads"] != threads:
            raise RuntimeError("Cached result uses a different thread count; select a new output directory")
        return json.loads(result_path.read_text(encoding="utf-8"))
    run, manifest = prepare(case, backend, output, threads)
    env = dict(os.environ)
    env["DBEM_GMRES_BACKEND"] = backend
    build = ROOT / "gpu/build"
    env["DBEM_CUDA_GMRES_DLL"] = str(build / ("dbem_cuda_gmres.dll" if backend == "cuda" else "intentionally-missing.dll"))
    cuda = Path(env.get("CUDA_PATH", r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6"))
    env["PATH"] = str(cuda / "bin") + os.pathsep + env.get("PATH", "")
    gpu_samples = []
    peak_rss = 0
    timed_out = False
    print(f"START {case} {backend}", flush=True)
    begin = time.perf_counter()
    with (run.parent / "run.log").open("wb") as log:
        process = subprocess.Popen([str(build / "DBEM1_cuda_optional.exe")], cwd=run, env=env,
                                   stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                   creationflags=subprocess.CREATE_NO_WINDOW)
        save(run.parent / "running.json", dict(pid=process.pid, started_unix=time.time(), manifest=manifest))
        monitored = psutil.Process(process.pid)
        last_gpu = -100
        while process.poll() is None:
            elapsed = time.perf_counter() - begin
            try:
                peak_rss = max(peak_rss, monitored.memory_info().rss)
            except psutil.NoSuchProcess:
                pass
            if elapsed - last_gpu >= 15:
                last_gpu = elapsed
                try:
                    query = subprocess.run(["nvidia-smi", "--query-gpu=memory.used,memory.free,utilization.gpu", "--format=csv,noheader,nounits"],
                                           capture_output=True, text=True, timeout=5,
                                           creationflags=subprocess.CREATE_NO_WINDOW)
                    gpu_samples.append(dict(seconds=elapsed, system_gpu_csv=query.stdout.strip()))
                    save(run.parent / "gpu_samples.json", gpu_samples)
                except (OSError, subprocess.TimeoutExpired):
                    pass
            if elapsed > timeout:
                timed_out = True
                process.terminate()
                process.wait(timeout=20)
                break
            try:
                process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                pass
        returncode = process.wait()
    elapsed = time.perf_counter() - begin
    result = parse_run(run, elapsed, returncode, peak_rss, gpu_samples, timed_out)
    save(result_path, result)
    (run.parent / "running.json").unlink()
    print(f"DONE {case} {backend}: success={result['success']} process={elapsed:.3f}s steps={result['solved_steps']} GMRES={result.get('gmres_seconds')}", flush=True)
    return result


def compare(cpu, gpu, elements):
    keys = ["step", "element", "localNode", "nodeId", "domain", "surfaceType", "bcid"]
    values = ["ux", "uy", "uz", "tx", "ty", "tz"]
    chunks = [pd.read_csv(run["state_file"], usecols=keys + values, chunksize=100000, float_precision="round_trip") for run in [cpu, gpu]]
    totals = {prefix: dict(reference_squared=0., difference_squared=0., max_absolute_difference=0., reference_max=0.) for prefix in ["u", "t"]}
    rows = 0
    step_rows = []
    for a, b in itertools.zip_longest(*chunks):
        assert a is not None and b is not None and a[keys].equals(b[keys]), "state row identity mismatch"
        rows += len(a)
        for prefix in ["u", "t"]:
            cols = [prefix + axis for axis in "xyz"]
            x, y = a[cols].to_numpy(), b[cols].to_numpy()
            assert np.isfinite(x).all() and np.isfinite(y).all(), "non-finite state output"
            difference = x - y
            total = totals[prefix]
            total["reference_squared"] += float(np.sum(x*x))
            total["difference_squared"] += float(np.sum(difference*difference))
            total["max_absolute_difference"] = max(total["max_absolute_difference"], float(np.max(np.abs(difference))))
            total["reference_max"] = max(total["reference_max"], float(np.max(np.abs(x))))
            partial = pd.DataFrame(dict(step=a["step"].to_numpy(), reference=np.sum(x*x, axis=1),
                                        difference=np.sum(difference*difference, axis=1)))
            grouped = partial.groupby("step").sum()
            for step, row in grouped.iterrows():
                step_rows.append(dict(step=int(step), field=prefix, **row.to_dict()))
    assert rows == elements * 1001, (rows, elements)
    per_step = pd.DataFrame(step_rows).groupby(["step", "field"]).sum().reset_index()
    per_step["relative_difference"] = np.sqrt(per_step["difference"]) / np.maximum(np.sqrt(per_step["reference"]), 1e-300)
    for prefix, total in totals.items():
        total["relative_difference"] = (total["difference_squared"] / max(total["reference_squared"], 1e-300))**.5
    return dict(rows=rows, finite=True, fields=totals,
                passed=all(t["relative_difference"] < 1e-5 for t in totals.values())), per_step


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", type=int, default=[4, 1, 9, 16])
    parser.add_argument("--timeout", type=float, default=14400, help="per solver process, seconds")
    parser.add_argument("--threads", type=int, default=40)
    parser.add_argument("--output", type=Path, default=ROOT / "gpu/validation/tsv_copper_t40_20261029")
    args = parser.parse_args()
    if os.name != "nt":
        parser.error("Use Windows Python to launch the existing Windows solver")
    args.output.mkdir(parents=True, exist_ok=True)
    source_files = sorted((ROOT / "DBEM1").glob("*.cpp")) + sorted((ROOT / "DBEM1").glob("*.h"))
    source_files += [ROOT / "gpu/gmres_cuda.cu", ROOT / "gpu/build/DBEM1_cuda_optional.exe", ROOT / "gpu/build/dbem_cuda_gmres.dll"]
    identity = {str(f.relative_to(ROOT)): sha256(f) for f in source_files}
    identity_path = args.output / "source_identity.json"
    if identity_path.exists():
        assert json.loads(identity_path.read_text()) == identity, "Source/binary changed since the previous run"
    else:
        save(identity_path, identity)
    report = dict(global_material=dict(E=110., nu=.34, rho=8.96), steps=1000,
                  cpu_threads=args.threads, logical_processors=os.cpu_count(), cases=[],
                  notes=["Single full-length run per backend; no statistical speedup confidence interval.",
                         "GPU memory samples are whole-device usage including desktop applications.",
                         "Process time includes assembly, initialization, 1000 time steps and output writing.",
                         "State comparison checks backend consistency, not physical convergence."])
    for index in args.cases:
        case = f"tsv_c{index:02d}_md{index*2:02d}_orig_b060_n1000"
        check = json.loads((ROOT / "DBEM1/input" / (case + ".selfcheck.json")).read_text())
        record = dict(case=case, domains=check["domain_count"], elements=check["element_count"], runs={})
        # Alternate first backend by size; all solver processes run sequentially.
        for backend in (["cuda", "cpu"] if index in [4, 9] else ["cpu", "cuda"]):
            record["runs"][backend] = run_one(case, backend, args.output, args.timeout, args.threads)
            save(args.output / "latest.json", dict(case=case, backend=backend, result=record["runs"][backend]))
        cpu, gpu = record["runs"]["cpu"], record["runs"]["cuda"]
        if cpu["success"] and gpu["success"]:
            record["comparison"], per_step = compare(cpu, gpu, record["elements"])
            per_step.to_csv(args.output / case / "per_step_comparison.csv", index=False)
            record["process_speedup"] = cpu["process_seconds"] / gpu["process_seconds"]
            record["gmres_speedup"] = cpu["gmres_seconds"] / gpu["gmres_seconds"]
        report["cases"].append(record)
        save(args.output / "summary.json", report)
        print(f"CASE {case}: process_speedup={record.get('process_speedup')} gmres_speedup={record.get('gmres_speedup')} comparison={record.get('comparison')}", flush=True)
    print(f"Report: {args.output / 'summary.json'}", flush=True)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
