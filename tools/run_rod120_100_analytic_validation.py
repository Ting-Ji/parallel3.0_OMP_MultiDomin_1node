"""Validate rod120_100 multi-domain material use against a 1D analytic baseline.

The test copies the rod120_100 case to an isolated run directory, intentionally
corrupts the global material line in BEM_DATACARD.DAT, keeps the multi-domain
material lines unchanged, runs DBEM1.exe, and compares the loaded-end response
with an analytic one-dimensional step-load estimate derived only from the
multi-domain material lines.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


CASE = "rod120_100"


@dataclass
class Material:
    e: float
    nu: float
    rho: float


@dataclass
class AnalyticResult:
    max_abs: float
    max_rel: float
    step: int
    time: float
    numeric: float
    analytic: float
    coord: Tuple[float, float, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--exe", type=Path, default=Path("x64/Release/DBEM1.exe"))
    parser.add_argument("--run-root", type=Path, default=Path("validation_runs/rod120_100_analytic"))
    parser.add_argument("--in-place", action="store_true",
                        help="temporarily run through repo/DBEM1/BEM_DATACARD.DAT and restore it afterwards")
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--abs-tol", type=float, default=2.5e-1)
    parser.add_argument("--rel-tol", type=float, default=3.0e-1)
    parser.add_argument("--pollute-e", type=float, default=200.0)
    parser.add_argument("--pollute-rho", type=float, default=20.0)
    return parser.parse_args()


def read_card(path: Path) -> List[str]:
    return path.read_text(encoding="ascii", errors="strict").splitlines()


def parse_materials(lines: Sequence[str]) -> Tuple[Material, List[Material], int]:
    global_mat = Material(float(lines[6]), float(lines[7]), float(lines[8]))
    domain_count = int(lines[16])
    mats: List[Material] = []
    for i in range(domain_count):
        parts = lines[17 + i].split()
        mats.append(Material(float(parts[1]), float(parts[2]), float(parts[3])))
    return global_mat, mats, domain_count


def write_polluted_card(src: Path, dst: Path, pollute_e: float, pollute_rho: float) -> None:
    lines = read_card(src)
    _, mats, domain_count = parse_materials(lines)
    lines[6] = f"{pollute_e:.17g}"
    lines[8] = f"{pollute_rho:.17g}"
    # Preserve multi-domain material rows exactly except for normalizing whitespace.
    for i, mat in enumerate(mats):
        lines[17 + i] = f"{i + 1} {mat.e:.17g} {mat.nu:.17g} {mat.rho:.17g}"
    if domain_count != 5:
        raise RuntimeError(f"expected rod120_100 DomainCount=5, got {domain_count}")
    dst.write_text("\n".join(lines) + "\n", encoding="ascii")


def copy_case(repo: Path, exe: Path, run_root: Path, pollute_e: float, pollute_rho: float) -> Path:
    stamp = time.strftime("%Y%m%d_%H%M%S")
    case_dir = run_root / stamp
    input_dir = case_dir / "input"
    output_dir = case_dir / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    # DBEM.cpp writes early log.txt/dt.txt to the parent output directory.
    (case_dir.parent / "output").mkdir(parents=True, exist_ok=True)
    shutil.copy2(repo / exe, case_dir / "DBEM1.exe")
    for suffix in (".model", ".bd"):
        shutil.copy2(repo / "DBEM1" / "input" / f"{CASE}{suffix}", input_dir / f"{CASE}{suffix}")
    write_polluted_card(
        repo / "DBEM1" / "input" / f"{CASE}.DATACARD",
        case_dir / "BEM_DATACARD.DAT",
        pollute_e,
        pollute_rho,
    )
    return case_dir


def case_output_dir(case_dir: Path) -> Path:
    output_root = case_dir / "output"
    candidates = sorted(path for path in output_root.glob(f"{CASE}_*") if path.is_dir())
    if candidates:
        return max(candidates, key=lambda path: path.stat().st_mtime)
    return output_root


def model_axis_and_length(path: Path) -> Tuple[int, float]:
    lines = path.read_text(encoding="ascii", errors="strict").splitlines()
    count = int(lines[0])
    coords = [[float(value) for value in lines[i].split()[:3]] for i in range(1, count + 1)]
    spans = [max(coord[i] for coord in coords) - min(coord[i] for coord in coords) for i in range(3)]
    axis = max(range(3), key=lambda i: spans[i])
    return axis, spans[axis]


def read_dt(case_dir: Path) -> float:
    candidates = [
        case_dir.parent / "output" / "dt.txt",
        case_output_dir(case_dir) / "dt.txt",
        case_dir / "output" / "dt.txt",
        case_dir / "dt.txt",
    ]
    for path in candidates:
        if path.exists():
            return float(path.read_text(encoding="ascii", errors="ignore").strip())
    raise RuntimeError("missing dt.txt")


def dominant_load_component(path: Path, element_count: int) -> int:
    sums = [0.0, 0.0, 0.0]
    with path.open("r", encoding="ascii", errors="strict") as fp:
        rows = [line.split() for line in fp.readlines()[element_count:]]
    for row in rows:
        if len(row) != 3:
            continue
        for i in range(3):
            sums[i] += abs(float(row[i]))
    component = max(range(3), key=lambda i: sums[i])
    if sums[component] <= 0.0:
        raise RuntimeError("rod120_100 boundary file has no non-zero load history")
    return component


def element_count_from_bd(path: Path) -> int:
    count = 0
    with path.open("r", encoding="ascii", errors="strict") as fp:
        for line in fp:
            parts = line.split()
            if len(parts) != 1:
                break
            count += 1
    return count


def analytic_displacement(material: Material, length: float, dt: float, step: int, load: float) -> float:
    wave_speed = math.sqrt(material.e / material.rho)
    t = step * dt
    static_tip = load * length / material.e
    correction = 0.0
    # Fixed left end and step traction at the right end:
    # u(L,t)=P*L/E - 2P/(E*L) * sum cos(c*k_n*t)/k_n^2,
    # k_n=(n+1/2)pi/L.
    for n in range(10000):
        k = (n + 0.5) * math.pi / length
        correction += math.cos(wave_speed * k * t) / (k * k)
    return static_tip - (2.0 * load / (material.e * length)) * correction


def read_domain_material_metrics(case_dir: Path) -> Dict[str, str]:
    metrics: Dict[str, str] = {}
    path = case_output_dir(case_dir) / "validation_metrics.txt"
    if not path.exists():
        return metrics
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            metrics[key] = value
    return metrics


def expected_wave_speeds(material: Material) -> Tuple[float, float]:
    shear = material.e / (2.0 * (1.0 + material.nu))
    c1 = math.sqrt(2.0 * shear * (1.0 - material.nu) /
                   (material.rho * (1.0 + material.nu) * (1.0 - 2.0 * material.nu)))
    c2 = math.sqrt(shear / material.rho)
    return c1, c2


def parse_named_float(text: str, name: str) -> float:
    match = re.search(rf"\b{name}\s*[= ]\s*([-+0-9.eE]+)", text)
    if not match:
        raise RuntimeError(f"missing {name} in material log line: {text}")
    return float(match.group(1))


def material_logs_match(case_dir: Path, material: Material, domain_count: int) -> bool:
    expected_c1, expected_c2 = expected_wave_speeds(material)
    rel_tol = 1.0e-6

    metrics = read_domain_material_metrics(case_dir)
    metric_lines = [value for key, value in metrics.items() if key.startswith("DomainMaterial[")]
    if len(metric_lines) != domain_count:
        return False
    for line in metric_lines:
        c1 = parse_named_float(line, "C1")
        c2 = parse_named_float(line, "C2")
        if abs(c1 - expected_c1) > rel_tol * max(1.0, abs(expected_c1)):
            return False
        if abs(c2 - expected_c2) > rel_tol * max(1.0, abs(expected_c2)):
            return False

    run_log = (case_dir / "run.log").read_text(encoding="utf-8", errors="replace")
    gmres_lines = [line for line in run_log.splitlines() if "MultiDomain GMRES material domain=" in line]
    if len(gmres_lines) != domain_count:
        return False
    for line in gmres_lines:
        c1 = parse_named_float(line, "C1")
        c2 = parse_named_float(line, "C2")
        if abs(c1 - expected_c1) > rel_tol * max(1.0, abs(expected_c1)):
            return False
        if abs(c2 - expected_c2) > rel_tol * max(1.0, abs(expected_c2)):
            return False
    return True


def compare_state(case_dir: Path, material: Material, length: float, axis_component: int) -> AnalyticResult:
    dt = read_dt(case_dir)
    state_path = case_output_dir(case_dir) / "validation_state.csv"
    if not state_path.exists():
        raise RuntimeError(f"missing {state_path}")
    component_field = ("ux", "uy", "uz")[axis_component]
    coord_fields = ("x", "y", "z")
    by_step: Dict[int, List[Tuple[float, Tuple[float, float, float]]]] = {}
    max_axis_coord = None
    with state_path.open("r", encoding="utf-8", errors="replace", newline="") as fp:
        rows = list(csv.DictReader(fp))
    for row in rows:
        value = float(row[coord_fields[axis_component]])
        max_axis_coord = value if max_axis_coord is None else max(max_axis_coord, value)
    if max_axis_coord is None:
        raise RuntimeError("validation_state.csv has no rows")
    for row in rows:
        if int(row["surfaceType"]) != 0 or int(row["bcid"]) != 456:
            continue
        axis_coord = float(row[coord_fields[axis_component]])
        if abs(axis_coord - max_axis_coord) > 1.0e-8:
            continue
        step = int(row["step"])
        coord = tuple(float(row[name]) for name in coord_fields)
        by_step.setdefault(step, []).append((float(row[component_field]), coord))
    if not by_step:
        raise RuntimeError("no loaded-end outer rows found in validation_state.csv")

    worst = AnalyticResult(0.0, 0.0, 0, 0.0, 0.0, 0.0, (0.0, 0.0, 0.0))
    for step, values in sorted(by_step.items()):
        numeric = sum(value for value, _ in values) / float(len(values))
        coord = values[0][1]
        expected = analytic_displacement(material, length, dt, step, 1.0)
        diff = abs(numeric - expected)
        rel = diff / max(abs(expected), 1.0e-12)
        if diff > worst.max_abs:
            worst = AnalyticResult(diff, rel, step, step * dt, numeric, expected, coord)
    return worst


def run_case(case_dir: Path, exe: Path, timeout: int) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["DBEM_VALIDATION_OUTPUT"] = "1"
    env["DBEM_MULTIDOMAIN_GMRES_PTHREAD_CHECK"] = "1"
    if os.name == "nt":
        command = [str(exe)]
        cwd = case_dir
    else:
        win_case_dir = subprocess.check_output(
            ["wslpath", "-w", str(case_dir)],
            text=True,
        ).strip()
        command = ["cmd.exe", "/c", f'cd /d "{win_case_dir}" && DBEM1.exe']
        cwd = None
    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
    )
    (case_dir / "run.log").write_text(result.stdout, encoding="utf-8", errors="replace")
    return result


def main() -> int:
    args = parse_args()
    repo = args.repo.resolve()
    backup_card: bytes | None = None
    card_path: Path | None = None
    if args.in_place:
        case_dir = repo / "DBEM1"
        exe = (repo / args.exe).resolve()
        card_path = case_dir / "BEM_DATACARD.DAT"
        backup_card = card_path.read_bytes() if card_path.exists() else None
        write_polluted_card(
            repo / "DBEM1" / "input" / f"{CASE}.DATACARD",
            card_path,
            args.pollute_e,
            args.pollute_rho,
        )
    else:
        case_dir = copy_case(repo, args.exe, repo / args.run_root, args.pollute_e, args.pollute_rho)
        exe = case_dir / "DBEM1.exe"
    print(f"case_dir={case_dir}")
    try:
        result = run_case(case_dir, exe, args.timeout)
        if result.returncode not in (0, 1):
            print(f"run failed: returncode={result.returncode}")
            return 2

        card_lines = read_card(case_dir / "BEM_DATACARD.DAT")
        global_mat, mats, domain_count = parse_materials(card_lines)
        if global_mat.e == mats[0].e and global_mat.rho == mats[0].rho:
            raise RuntimeError("global material was not polluted")
        if len({(m.e, m.nu, m.rho) for m in mats}) != 1:
            raise RuntimeError("rod120_100 analytic validation expects uniform multi-domain material")

        element_count = element_count_from_bd(case_dir / "input" / f"{CASE}.bd")
        load_component, length = model_axis_and_length(case_dir / "input" / f"{CASE}.model")
        worst = compare_state(case_dir, mats[0], length, load_component)

        material_ok = material_logs_match(case_dir, mats[0], domain_count)
        passed = worst.max_abs <= args.abs_tol or worst.max_rel <= args.rel_tol
        report = (
            f"max_abs={worst.max_abs:.6e}\n"
            f"max_rel={worst.max_rel:.6e}\n"
            f"step={worst.step}\n"
            f"time={worst.time:.17g}\n"
            f"numeric={worst.numeric:.17g}\n"
            f"analytic={worst.analytic:.17g}\n"
            f"coord={worst.coord}\n"
            f"rod_axis_component={load_component}\n"
            f"material_metrics_ok={int(material_ok)}\n"
            f"passed={int(passed and material_ok)}\n"
        )
        (case_dir / "analytic_validation.txt").write_text(report, encoding="ascii")
        print(report)
        return 0 if passed and material_ok else 1
    finally:
        if args.in_place and card_path:
            if backup_card is None:
                card_path.unlink(missing_ok=True)
            else:
                card_path.write_bytes(backup_card)


if __name__ == "__main__":
    raise SystemExit(main())
