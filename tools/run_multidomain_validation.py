#!/usr/bin/env python3
"""
Generate isolated rod cases and validate the multi-domain TDBEM path.

Default quick mode runs:
  - 1-domain old/new, NStep=2
  - 1-domain old/new, NStep=5
  - 2-domain same material, NStep=5

Use --ten to add 10-domain cases without 100-domain cases.
Use --full to add 10-domain, 100-domain, and 100-domain segmented-material cases.
Use --diagnose-two-domain to run only the 1-domain gate and 2-domain RHS diagnostics.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


DEFAULT_E = 20000.0
DEFAULT_V = 0.25
DEFAULT_ROU = 8.0
DEFAULT_BETA = 0.4
DEFAULT_GAP = 0.67
DEFAULT_THREADS = 1
DEFAULT_ITERATIONS = 300
DEFAULT_ERROR = 1.0e-8
DEFAULT_PRE_LEAF = 10
ROD_LENGTH_PER_DOMAIN = 1.0
ROD_HALF_WIDTH = 0.1
LOAD_VALUE = 1.0
DEFAULT_MESH_NX = 1
DEFAULT_MESH_NY = 1
DEFAULT_MESH_NZ = 1

STATE_FIELDS = ("ux", "uy", "uz", "tx", "ty", "tz")
RHS_FIELDS = ("knownRhs", "historyG", "historyT", "totalRhs", "Ax", "residual", "solution")
RHS_COMPARE_FIELDS = ("knownRhs", "historyG", "historyT", "totalRhs", "Ax", "residual")


@dataclass
class Material:
    e: float = DEFAULT_E
    v: float = DEFAULT_V
    rou: float = DEFAULT_ROU


@dataclass
class CaseSpec:
    name: str
    domain_count: int
    nstep: int
    multi_domain: bool
    single_baseline_segments: Optional[int] = None
    segmented_material: bool = False
    beta_override: Optional[float] = None
    rhs_diagnostic: bool = False


@dataclass
class RunResult:
    spec: CaseSpec
    case_dir: Path
    returncode: int
    metrics: Dict[str, str]
    log_text: str
    state_path: Path


@dataclass
class CompareResult:
    name: str
    passed: bool
    max_abs: float
    max_rel: float
    compared_values: int
    missing: int = 0
    message: str = ""


def is_success_returncode(returncode: int) -> bool:
    # This legacy program returns 1 from main() on a normal completed run.
    return returncode in (0, 1)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def find_exe(root: Path, explicit: Optional[str]) -> Path:
    if explicit:
        exe = Path(explicit)
        if not exe.is_absolute():
            exe = root / exe
        if not exe.exists():
            raise FileNotFoundError(f"DBEM1.exe not found: {exe}")
        return exe

    candidates = [
        root / "DBEM1" / "x64" / "Release" / "DBEM1.exe",
        root / "x64" / "Release" / "DBEM1.exe",
        root / "DBEM1" / "Release" / "DBEM1.exe",
    ]
    for exe in candidates:
        if exe.exists():
            return exe
    raise FileNotFoundError("DBEM1.exe not found. Build Release|x64 first.")


def clean_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def midpoint(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> Tuple[float, float, float]:
    return ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, (a[2] + b[2]) * 0.5)


def _lerp3(a: Tuple[float, float, float], b: Tuple[float, float, float], t: float) -> Tuple[float, float, float]:
    """Linear interpolate between two 3D points."""
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t)


class RodMeshBuilder:
    def __init__(self) -> None:
        self.points: List[Tuple[float, float, float]] = []
        self.elements: List[List[int]] = []
        self.face_names: List[str] = []

    def add_point(self, p: Tuple[float, float, float]) -> int:
        self.points.append(p)
        return len(self.points)

    def add_face(
        self,
        name: str,
        p0: Tuple[float, float, float],
        p1: Tuple[float, float, float],
        p2: Tuple[float, float, float],
        p3: Tuple[float, float, float],
    ) -> int:
        nodes = [
            p0,
            p1,
            p2,
            p3,
            midpoint(p0, p1),
            midpoint(p1, p2),
            midpoint(p2, p3),
            midpoint(p3, p0),
        ]
        self.elements.append([self.add_point(p) for p in nodes])
        self.face_names.append(name)
        return len(self.elements) - 1

    def add_subdivided_face(
        self,
        name: str,
        p0: Tuple[float, float, float],
        p1: Tuple[float, float, float],
        p2: Tuple[float, float, float],
        p3: Tuple[float, float, float],
        nu: int,
        nv: int,
    ) -> List[int]:
        """Add a subdivided quad face with nu * nv elements.

        The face is treated as a bilinear patch from p0..p3.  p0→p1 is the
        u-direction, p0→p3 is the v-direction.  Elements are looped u then v so
        that touching faces on adjacent domains pair correctly.
        """
        elem_ids: List[int] = []
        for iv in range(nv):
            for iu in range(nu):
                u0, u1 = iu / nu, (iu + 1) / nu
                v0, v1 = iv / nv, (iv + 1) / nv
                c00 = _lerp3(_lerp3(p0, p1, u0), _lerp3(p3, p2, u0), v0)
                c10 = _lerp3(_lerp3(p0, p1, u1), _lerp3(p3, p2, u1), v0)
                c11 = _lerp3(_lerp3(p0, p1, u1), _lerp3(p3, p2, u1), v1)
                c01 = _lerp3(_lerp3(p0, p1, u0), _lerp3(p3, p2, u0), v1)
                elem_ids.append(self.add_face(name, c00, c10, c11, c01))
        return elem_ids

    def face_element_indices(self, face_name: str) -> List[int]:
        """Return 0-based element indices for all elements with the given face name."""
        return [i for i, name in enumerate(self.face_names) if name == face_name]


def add_box_segment_faces(builder: RodMeshBuilder, x0: float, x1: float,
                          include_left: bool, include_right: bool,
                          nx: int = 1, ny: int = 1, nz: int = 1) -> None:
    """Add the six faces of a rectangular-prism segment, optionally subdivided.

    Face element order (consistent with generate_multidomain_rod.py):
      left (y-z plane, ny*nz), right (y-z plane, ny*nz),
      bottom / y-min (x-z plane, nx*nz), top / y-max (x-z plane, nx*nz),
      front / z-min (x-y plane, nx*ny), back / z-max (x-y plane, nx*ny).
    """
    a = ROD_HALF_WIDTH

    # Left face (x = x0, normal = -x, y-z plane): ny * nz elements
    if include_left:
        x = x0
        builder.add_subdivided_face("left", (x, -a, a), (x, a, a), (x, a, -a), (x, -a, -a), ny, nz)

    # Bottom / y-min face (x-z plane): nx * nz elements
    builder.add_subdivided_face("bottom", (x0, -a, -a), (x1, -a, -a), (x1, -a, a), (x0, -a, a), nx, nz)
    # Top / y-max face (x-z plane): nx * nz elements
    builder.add_subdivided_face("top", (x0, a, -a), (x0, a, a), (x1, a, a), (x1, a, -a), nx, nz)
    # Front / z-min face (x-y plane): nx * ny elements
    builder.add_subdivided_face("front", (x0, -a, -a), (x0, a, -a), (x1, a, -a), (x1, -a, -a), nx, ny)
    # Back / z-max face (x-y plane): nx * ny elements
    builder.add_subdivided_face("back", (x0, -a, a), (x1, -a, a), (x1, a, a), (x0, a, a), nx, ny)

    # Right face (x = x1, normal = +x, y-z plane): ny * nz elements
    if include_right:
        x = x1
        builder.add_subdivided_face("right", (x, -a, -a), (x, a, -a), (x, a, a), (x, -a, a), ny, nz)


def generate_single_boundary_mesh(segments: int,
                                   nx: int = DEFAULT_MESH_NX,
                                   ny: int = DEFAULT_MESH_NY,
                                   nz: int = DEFAULT_MESH_NZ) -> RodMeshBuilder:
    builder = RodMeshBuilder()
    for d in range(segments):
        x0 = d * ROD_LENGTH_PER_DOMAIN
        x1 = (d + 1) * ROD_LENGTH_PER_DOMAIN
        add_box_segment_faces(builder, x0, x1,
                              include_left=(d == 0), include_right=(d == segments - 1),
                              nx=nx, ny=ny, nz=nz)
    return builder


def generate_multidomain_mesh(domains: int,
                               nx: int = DEFAULT_MESH_NX,
                               ny: int = DEFAULT_MESH_NY,
                               nz: int = DEFAULT_MESH_NZ) -> RodMeshBuilder:
    builder = RodMeshBuilder()
    for d in range(domains):
        x0 = d * ROD_LENGTH_PER_DOMAIN
        x1 = (d + 1) * ROD_LENGTH_PER_DOMAIN
        add_box_segment_faces(builder, x0, x1, include_left=True, include_right=True,
                              nx=nx, ny=ny, nz=nz)
    return builder


def average_element_size(builder: RodMeshBuilder) -> float:
    total = 0.0
    for element in builder.elements:
        p0 = builder.points[element[0] - 1]
        p1 = builder.points[element[1] - 1]
        total += math.dist(p0, p1)
    return total / len(builder.elements)


def single_baseline_beta_for_matching_dt(segments: int) -> float:
    single_size = average_element_size(generate_single_boundary_mesh(segments))
    multi_size = average_element_size(generate_multidomain_mesh(segments))
    return DEFAULT_BETA * multi_size / single_size


def write_model(path: Path, builder: RodMeshBuilder) -> None:
    with path.open("w", newline="\n") as f:
        f.write(f"{len(builder.points)}\n")
        for x, y, z in builder.points:
            f.write(f"{x:.17g} {y:.17g} {z:.17g}\n")
        f.write(f"{len(builder.elements)}\n")
        for element in builder.elements:
            f.write(" ".join(str(node_id) for node_id in element))
            f.write("\n")
        f.write("0\n")


def is_left_end(face: str, element_index: int, domain_count: int, multi_domain: bool) -> bool:
    """Check if an element belongs to the left end cap (fixed BC, BCID=123).

    In single-domain mode, all "left" elements are the left end.
    In multi-domain mode, only the "left" face elements in the first domain
    are the true left end; interior "left" faces are interfaces.  We can't
    distinguish them with just the face name, so this function returns True
    for ALL "left" faces.  Callers must further filter by domain position.
    """
    return face == "left"


def is_right_end(face: str, element_index: int, domain_count: int, multi_domain: bool) -> bool:
    """Check if an element belongs to the right end cap (loaded face).

    Same caveat as is_left_end: returns True for ALL "right" faces.
    Callers must further filter by domain position.
    """
    return face == "right"


def _first_domain_left_element_count(builder: RodMeshBuilder, domain_count: int) -> int:
    """Return how many elements in the first domain have face name 'left'."""
    if domain_count <= 0:
        return 0
    total = len(builder.face_names)
    ele_per_domain = total // domain_count
    count = 0
    for i in range(ele_per_domain):
        if builder.face_names[i] == "left":
            count += 1
    return count


def _last_domain_right_element_range(builder: RodMeshBuilder, domain_count: int) -> Tuple[int, int]:
    """Return (first, last+1) 0-based element indices for the last domain's right face."""
    if domain_count <= 0:
        return (0, 0)
    total = len(builder.face_names)
    ele_per_domain = total // domain_count
    start = (domain_count - 1) * ele_per_domain
    first = -1
    last = -1
    for i in range(start, total):
        if builder.face_names[i] == "right":
            if first < 0:
                first = i
            last = i + 1
    return (first, last) if first >= 0 else (0, 0)


def write_bd(path: Path, builder: RodMeshBuilder, nstep: int, domain_count: int, multi_domain: bool) -> None:
    # Determine which elements are on the exterior ends.
    # In single-domain mode, all "left" / "right" faces are exterior.
    # In multi-domain mode, only the first domain's "left" and the last
    # domain's "right" are exterior.
    left_end_count = _first_domain_left_element_count(builder, domain_count)
    right_start, right_end = _last_domain_right_element_range(builder, domain_count)
    left_end_set: Set[int] = set()
    right_end_set: Set[int] = set()

    total = len(builder.face_names)
    ele_per_domain = total // domain_count if domain_count > 0 else total

    for i, face in enumerate(builder.face_names):
        if multi_domain:
            # First domain's left elements (indices 0..ele_per_domain-1)
            if face == "left" and i < ele_per_domain:
                left_end_set.add(i)
                continue
            # Last domain's right elements
            if face == "right" and right_start <= i < right_end:
                right_end_set.add(i)
                continue
        else:
            if face == "left":
                left_end_set.add(i)
            if face == "right":
                right_end_set.add(i)

    bcids: List[int] = []
    for i in range(total):
        if i in left_end_set:
            bcids.append(123)
        else:
            bcids.append(456)

    with path.open("w", newline="\n") as f:
        for bcid in bcids:
            f.write(f"{bcid}\n")
        for step in range(nstep + 1):
            for i in range(total):
                if i in left_end_set:
                    value = (0.0, 0.0, 0.0)
                elif i in right_end_set:
                    pulse = LOAD_VALUE if step > 0 else 0.0
                    value = (pulse, 0.0, 0.0)
                else:
                    value = (0.0, 0.0, 0.0)
                f.write(f"{value[0]:.17g} {value[1]:.17g} {value[2]:.17g}\n")


def material_for_domain(domain_index: int, segmented: bool) -> Material:
    if not segmented:
        return Material()
    group = domain_index // 20
    return Material(e=DEFAULT_E * (1.0 + 0.15 * group), v=DEFAULT_V, rou=DEFAULT_ROU)


def write_datacard(path: Path, title: str, spec: CaseSpec, builder: Optional[RodMeshBuilder] = None) -> None:
    beta = spec.beta_override if spec.beta_override is not None else DEFAULT_BETA
    lines = [
        title,
        str(DEFAULT_THREADS),
        "2",
        f"{beta:.17g}",
        str(spec.nstep),
        f"{DEFAULT_GAP:.17g}",
        f"{DEFAULT_E:.17g}",
        f"{DEFAULT_V:.17g}",
        f"{DEFAULT_ROU:.17g}",
        "1",
        "2",
        str(DEFAULT_ITERATIONS),
        f"{DEFAULT_ERROR:.17g}",
        str(DEFAULT_PRE_LEAF),
        "0",
        "0",
    ]

    if spec.multi_domain:
        lines.append(str(spec.domain_count))
        for d in range(spec.domain_count):
            mat = material_for_domain(d, spec.segmented_material)
            lines.append(f"{d + 1} {mat.e:.17g} {mat.v:.17g} {mat.rou:.17g}")
        if spec.domain_count == 1:
            lines.append("0")
        else:
            # Build interface pairs from the builder's face names when available
            interfaces = _build_interfaces_from_builder(builder, spec.domain_count)
            lines.append(str(len(interfaces)))
            for domain_a, ele_a, domain_b, ele_b, normal_sign in interfaces:
                lines.append(f"{domain_a} {ele_a} {domain_b} {ele_b} {normal_sign}")

    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def _build_interfaces_from_builder(builder: Optional[RodMeshBuilder],
                                    domain_count: int) -> List[Tuple[int, int, int, int, int]]:
    """Build interface tuples from the mesh builder's face names.

    Pairs each "right" face element in domain d with the corresponding
    "left" face element in domain d+1.  The builder must emit elements
    in a consistent face order (right face elements in the same relative
    positions as left face elements of the next domain).

    Falls back to the old 6-per-domain formula when no builder is available.
    """
    if builder is None or not builder.face_names:
        # Fallback for backwards compatibility (assumes 6 elements per domain)
        interfaces: List[Tuple[int, int, int, int, int]] = []
        for d in range(domain_count - 1):
            ele_a = d * 6 + 6
            ele_b = (d + 1) * 6 + 1
            interfaces.append((d + 1, ele_a, d + 2, ele_b, -1))
        return interfaces

    total = len(builder.face_names)
    ele_per_domain = total // domain_count
    interfaces: List[Tuple[int, int, int, int, int]] = []

    for d in range(domain_count - 1):
        domain_a = d + 1
        domain_b = d + 2
        # Collect "right" elements of domain d and "left" elements of domain d+1
        right_elems: List[int] = []
        left_elems: List[int] = []
        dom_a_start = d * ele_per_domain
        dom_a_end = dom_a_start + ele_per_domain
        dom_b_start = (d + 1) * ele_per_domain
        dom_b_end = dom_b_start + ele_per_domain

        for i in range(dom_a_start, dom_a_end):
            if builder.face_names[i] == "right":
                right_elems.append(i + 1)  # 1-based
        for i in range(dom_b_start, dom_b_end):
            if builder.face_names[i] == "left":
                left_elems.append(i + 1)  # 1-based

        if len(right_elems) != len(left_elems):
            raise RuntimeError(
                f"Interface element count mismatch at domain boundary {d}: "
                f"right={len(right_elems)} left={len(left_elems)}"
            )
        for ele_a, ele_b in zip(right_elems, left_elems):
            interfaces.append((domain_a, ele_a, domain_b, ele_b, -1))

    return interfaces


def prepare_case(run_root: Path, repo: Path, exe: Path, spec: CaseSpec) -> Path:
    case_dir = run_root / f"case_{clean_name(spec.name)}"
    input_dir = case_dir / "input"
    output_dir = case_dir / "output"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    (run_root / "output").mkdir(parents=True, exist_ok=True)

    shutil.copy2(exe, case_dir / "DBEM1.exe")
    for helper in ("gaussposition", "gaussweight"):
        shutil.copy2(repo / "DBEM1" / helper, case_dir / helper)

    title = "rod_validation"
    if spec.multi_domain:
        builder = generate_multidomain_mesh(spec.domain_count)
    else:
        segments = spec.single_baseline_segments or spec.domain_count
        builder = generate_single_boundary_mesh(segments)

    write_model(input_dir / f"{title}.model", builder)
    write_bd(input_dir / f"{title}.bd", builder, spec.nstep, spec.domain_count, spec.multi_domain)
    write_datacard(case_dir / "BEM_DATACARD.DAT", title, spec, builder)
    return case_dir


def is_wsl() -> bool:
    return platform.system().lower() == "linux" and "microsoft" in platform.release().lower()


def run_case(case_dir: Path, spec: CaseSpec, timeout: int) -> RunResult:
    env = os.environ.copy()
    env["DBEM_VALIDATION_OUTPUT"] = "1"
    if spec.rhs_diagnostic:
        env["DBEM_RHS_DIAGNOSTIC"] = "1"
    exe = case_dir / "DBEM1.exe"
    command = [str(exe)]
    if is_wsl():
        command = [f"./{exe.name}"]

    completed = subprocess.run(
        command,
        cwd=case_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        errors="replace",
        timeout=timeout,
    )
    log_text = completed.stdout + "\n--- STDERR ---\n" + completed.stderr
    (case_dir / "run.log").write_text(log_text, encoding="utf-8", errors="replace")
    metrics_path = case_dir / "output" / "validation_metrics.txt"
    metrics = parse_metrics(metrics_path)
    return RunResult(
        spec=spec,
        case_dir=case_dir,
        returncode=completed.returncode,
        metrics=metrics,
        log_text=log_text,
        state_path=case_dir / "output" / "validation_state.csv",
    )


def parse_metrics(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    if not path.exists():
        return result
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.strip()] = value.strip()
    return result


def read_state(path: Path) -> List[Dict[str, float]]:
    if not path.exists():
        raise FileNotFoundError(f"missing validation state: {path}")
    rows: List[Dict[str, float]] = []
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            row: Dict[str, float] = {}
            for key, value in raw.items():
                if key in ("step", "element", "localNode", "nodeId", "domain", "surfaceType", "bcid"):
                    row[key] = int(value)
                else:
                    row[key] = parse_float(value)
            rows.append(row)
    return rows


def read_rhs_breakdown(path: Path) -> List[Dict[str, object]]:
    if not path.exists():
        raise FileNotFoundError(f"missing RHS diagnostic: {path}")
    rows: List[Dict[str, object]] = []
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            row: Dict[str, object] = {}
            for key, value in raw.items():
                if key in ("step", "rowNode", "rowElement", "localNode", "domain", "surfaceType", "bcid"):
                    row[key] = int(value)
                elif key in ("component", "solver"):
                    row[key] = value
                else:
                    row[key] = parse_float(value)
            rows.append(row)
    return rows


def read_residual_probe(path: Path) -> List[Dict[str, object]]:
    if not path.exists():
        raise FileNotFoundError(f"missing residual diagnostic: {path}")
    rows: List[Dict[str, object]] = []
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            row: Dict[str, object] = {}
            for key, value in raw.items():
                if key == "step":
                    row[key] = int(value)
                elif key == "solver":
                    row[key] = value
                else:
                    row[key] = parse_float(value)
            rows.append(row)
    return rows


def parse_float(value: str) -> float:
    text = value.strip().lower()
    if "nan" in text:
        return float("nan")
    if text in ("inf", "+inf", "infinity", "+infinity"):
        return float("inf")
    if text in ("-inf", "-infinity"):
        return float("-inf")
    return float(value)


def rel_error(a: float, b: float) -> float:
    denom = max(abs(a), abs(b), 1.0e-300)
    return abs(a - b) / denom


def compare_rows(
    name: str,
    a_rows: Sequence[Dict[str, float]],
    b_rows: Sequence[Dict[str, float]],
    abs_tol: float,
    rel_tol: float,
    key_fields: Sequence[str],
    only_outer_b: bool = False,
) -> CompareResult:
    a_map: Dict[Tuple[object, ...], Dict[str, float]] = {}
    for row in a_rows:
        key = make_key(row, key_fields)
        a_map[key] = row

    max_abs = 0.0
    max_rel = 0.0
    count = 0
    missing = 0
    max_detail = ""

    for b in b_rows:
        if only_outer_b and int(b["surfaceType"]) != 0:
            continue
        key = make_key(b, key_fields)
        a = a_map.get(key)
        if a is None:
            missing += 1
            continue
        for field in STATE_FIELDS:
            diff = abs(float(a[field]) - float(b[field]))
            if math.isnan(diff) or math.isinf(diff):
                max_abs = math.inf
                max_rel = math.inf
                count += 1
                continue
            max_abs = max(max_abs, diff)
            rel = rel_error(float(a[field]), float(b[field]))
            max_rel = max(max_rel, rel)
            if diff >= max_abs:
                max_detail = (
                    f"field={field}, key={key}, "
                    f"a={float(a[field]):.17g}, b={float(b[field]):.17g}, "
                    f"abs={diff:.6e}, rel={rel:.6e}"
                )
            count += 1

    passed = missing == 0 and count > 0 and (max_abs <= abs_tol or max_rel <= rel_tol)
    msg = f"max_abs={max_abs:.6e}, max_rel={max_rel:.6e}, values={count}, missing={missing}"
    if max_detail:
        msg += f", top={max_detail}"
    return CompareResult(name, passed, max_abs, max_rel, count, missing, msg)


def make_key(row: Dict[str, float], key_fields: Sequence[str]) -> Tuple[object, ...]:
    key: List[object] = []
    for field in key_fields:
        if field in ("cx", "cy", "cz"):
            key.append(round(coord_value(row, field), 10))
        elif field in ("nx", "ny", "nz"):
            key.append(round(float(row.get(field, 0.0)), 10))
        else:
            key.append(int(row[field]))
    return tuple(key)


def coord_value(row: Dict[str, float], field: str) -> float:
    """Return a coordinate value from either center or nodal CSV columns.

    Older diagnostics wrote cx/cy/cz.  The current validation_state.csv emitted
    by DBEM1 writes x/y/z.  The comparison key only needs a stable spatial
    coordinate, so support both formats.
    """
    if field in row:
        return float(row[field])
    fallback = {"cx": "x", "cy": "y", "cz": "z"}.get(field)
    if fallback and fallback in row:
        return float(row[fallback])
    raise KeyError(field)


def check_md_run(result: RunResult, require_self_check: bool = False) -> CompareResult:
    if not is_success_returncode(result.returncode):
        return CompareResult(
            result.spec.name + " run",
            False,
            math.inf,
            math.inf,
            0,
            message=f"returncode={result.returncode}",
        )
    if require_self_check and "MultiDomain matrix self-check passed for DomainCount=1" not in result.log_text:
        return CompareResult(
            result.spec.name + " self-check",
            False,
            math.inf,
            math.inf,
            0,
            message="DomainCount=1 matrix self-check did not pass",
        )
    flag = result.metrics.get("MultiDomainFinalFlag")
    if result.spec.multi_domain and flag is not None and flag.strip() != "0":
        return CompareResult(
            result.spec.name + " gmres",
            False,
            math.inf,
            math.inf,
            0,
            message=f"MultiDomainFinalFlag={flag}",
        )
    if result.spec.multi_domain and result.metrics.get("MultiDomainMetricsValid") != "1":
        return CompareResult(
            result.spec.name + " metrics",
            False,
            math.inf,
            math.inf,
            0,
            message="MultiDomainMetricsValid is not 1",
        )
    return CompareResult(result.spec.name + " run", True, 0.0, 0.0, 0, message="ok")


def check_interface(result: RunResult, abs_tol: float, rel_tol: float) -> CompareResult:
    max_u = float(result.metrics.get("InterfaceOverallMaxAbsU", "inf"))
    max_t = float(result.metrics.get("InterfaceOverallMaxAbsTBalance", "inf"))
    scale = 0.0
    try:
        rows = read_state(result.state_path)
        for row in rows:
            for field in STATE_FIELDS:
                scale = max(scale, abs(float(row[field])))
    except Exception:
        scale = 0.0
    rel = max(max_u, max_t) / max(scale, 1.0e-300)
    passed = max(max_u, max_t) <= abs_tol or rel <= rel_tol
    msg = f"max_u={max_u:.6e}, max_t_balance={max_t:.6e}, rel={rel:.6e}"
    return CompareResult(result.spec.name + " interface", passed, max(max_u, max_t), rel, 0, message=msg)


def rhs_key(row: Dict[str, object]) -> Tuple[object, ...]:
    return (
        int(row["step"]),
        int(row.get("localNode", 0)),
        round(float(row["cx"]), 10),
        round(float(row["cy"]), 10),
        round(float(row["cz"]), 10),
        round(float(row.get("nx", 0.0)), 10),
        round(float(row.get("ny", 0.0)), 10),
        round(float(row.get("nz", 0.0)), 10),
        str(row["component"]),
    )


def check_residual_probe(result: RunResult, abs_tol: float, rel_tol: float) -> CompareResult:
    path = result.case_dir / "output" / "residual_probe.csv"
    try:
        rows = read_residual_probe(path)
    except Exception as exc:
        return CompareResult(result.spec.name + " residual probe", False, math.inf, math.inf, 0, message=str(exc))

    max_abs = 0.0
    max_rel = 0.0
    max_detail = ""
    for row in rows:
        abs_value = abs(float(row["maxAbsResidual"]))
        rel_value = abs(float(row["relativeResidual"]))
        if abs_value > max_abs or rel_value > max_rel:
            max_detail = (
                f"step={int(row['step'])}, abs={abs_value:.6e}, "
                f"rel={rel_value:.6e}, l2Rhs={float(row['l2Rhs']):.6e}"
            )
        max_abs = max(max_abs, abs_value)
        max_rel = max(max_rel, rel_value)
    passed = len(rows) > 0 and (max_abs <= abs_tol or max_rel <= rel_tol)
    return CompareResult(
        result.spec.name + " residual probe",
        passed,
        max_abs,
        max_rel,
        len(rows),
        message=f"max_abs={max_abs:.6e}, max_rel={max_rel:.6e}, top={max_detail}",
    )


def compare_rhs_breakdown(single: RunResult, md: RunResult) -> CompareResult:
    try:
        single_rows = read_rhs_breakdown(single.case_dir / "output" / "rhs_breakdown.csv")
        md_rows = read_rhs_breakdown(md.case_dir / "output" / "rhs_breakdown.csv")
    except Exception as exc:
        return CompareResult("2domain RHS outer diagnostic", False, math.inf, math.inf, 0, message=str(exc))

    single_map = {rhs_key(row): row for row in single_rows}
    diffs: List[Tuple[float, float, str, Dict[str, object], Dict[str, object]]] = []
    missing = 0
    compared = 0
    max_abs = 0.0
    max_rel = 0.0
    field_max: Dict[str, float] = {field: 0.0 for field in RHS_COMPARE_FIELDS}

    for md_row in md_rows:
        if int(md_row["surfaceType"]) != 0:
            continue
        key = rhs_key(md_row)
        single_row = single_map.get(key)
        if single_row is None:
            missing += 1
            continue
        for field in RHS_COMPARE_FIELDS:
            a = float(single_row[field])
            b = float(md_row[field])
            diff = abs(a - b)
            rel = rel_error(a, b)
            diffs.append((diff, rel, field, single_row, md_row))
            compared += 1
            max_abs = max(max_abs, diff)
            max_rel = max(max_rel, rel)
            field_max[field] = max(field_max[field], diff)

    diffs.sort(key=lambda item: item[0], reverse=True)
    diff_path = md.case_dir / "output" / "rhs_outer_diff.csv"
    with diff_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "rank", "field", "step", "localNode", "cx", "cy", "cz", "nx", "ny", "nz", "component",
            "single", "multidomain", "absDiff", "relDiff",
        ])
        for rank, (diff, rel, field, single_row, md_row) in enumerate(diffs[:50], start=1):
            writer.writerow([
                rank,
                field,
                int(md_row["step"]),
                int(md_row.get("localNode", -1)),
                f"{float(md_row['cx']):.17g}",
                f"{float(md_row['cy']):.17g}",
                f"{float(md_row['cz']):.17g}",
                f"{float(md_row.get('nx', 0.0)):.17g}",
                f"{float(md_row.get('ny', 0.0)):.17g}",
                f"{float(md_row.get('nz', 0.0)):.17g}",
                str(md_row["component"]),
                f"{float(single_row[field]):.17g}",
                f"{float(md_row[field]):.17g}",
                f"{diff:.17g}",
                f"{rel:.17g}",
            ])

    top_msg = "none"
    if diffs:
        diff, rel, field, single_row, md_row = diffs[0]
        top_msg = (
            f"field={field}, step={int(md_row['step'])}, "
            f"coord=({float(md_row['cx']):.10g},{float(md_row['cy']):.10g},{float(md_row['cz']):.10g}), "
            f"component={md_row['component']}, single={float(single_row[field]):.17g}, "
            f"md={float(md_row[field]):.17g}, abs={diff:.6e}, rel={rel:.6e}"
        )
    field_msg = ", ".join(f"{field}={field_max[field]:.3e}" for field in RHS_COMPARE_FIELDS)
    abs_tol = 1.0e-8
    rel_tol = 1.0e-5
    passed = missing == 0 and compared > 0 and (max_abs <= abs_tol or max_rel <= rel_tol)
    return CompareResult(
        "2domain RHS outer diagnostic",
        passed,
        max_abs,
        max_rel,
        compared,
        missing,
        f"max_abs={max_abs:.6e}, max_rel={max_rel:.6e}, missing={missing}, "
        f"top={top_msg}, by_field=[{field_msg}], report={diff_path}",
    )


def write_state_outer_diff_report(single: RunResult, md: RunResult) -> Path:
    single_rows = read_state(single.state_path)
    md_rows = read_state(md.state_path)
    outer_key = ("step", "localNode", "cx", "cy", "cz", "nx", "ny", "nz")
    single_map = {make_key(row, outer_key): row for row in single_rows}
    diffs: List[Tuple[float, float, str, Dict[str, float], Dict[str, float]]] = []
    for md_row in md_rows:
        if int(md_row["surfaceType"]) != 0:
            continue
        key = make_key(md_row, outer_key)
        single_row = single_map.get(key)
        if single_row is None:
            continue
        for field in STATE_FIELDS:
            a = float(single_row[field])
            b = float(md_row[field])
            diffs.append((abs(a - b), rel_error(a, b), field, single_row, md_row))

    diffs.sort(key=lambda item: item[0], reverse=True)
    diff_path = md.case_dir / "output" / "state_outer_diff.csv"
    with diff_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "rank", "field", "step", "localNode", "cx", "cy", "cz", "nx", "ny", "nz",
            "single", "multidomain", "absDiff", "relDiff",
        ])
        for rank, (diff, rel, field, single_row, md_row) in enumerate(diffs[:100], start=1):
            writer.writerow([
                rank,
                field,
                int(md_row["step"]),
                int(md_row.get("localNode", -1)),
                f"{coord_value(md_row, 'cx'):.17g}",
                f"{coord_value(md_row, 'cy'):.17g}",
                f"{coord_value(md_row, 'cz'):.17g}",
                f"{float(md_row.get('nx', 0.0)):.17g}",
                f"{float(md_row.get('ny', 0.0)):.17g}",
                f"{float(md_row.get('nz', 0.0)):.17g}",
                f"{float(single_row[field]):.17g}",
                f"{float(md_row[field]):.17g}",
                f"{diff:.17g}",
                f"{rel:.17g}",
            ])
    return diff_path


def case_specs(include_ten: bool, full: bool) -> List[CaseSpec]:
    specs = [
        CaseSpec("old_1domain_n2", 1, 2, False, single_baseline_segments=1),
        CaseSpec("md_1domain_n2", 1, 2, True),
        CaseSpec("old_1domain_n5", 1, 5, False, single_baseline_segments=1),
        CaseSpec("md_1domain_n5", 1, 5, True),
        CaseSpec("single_for_2domain_n5", 1, 5, False, single_baseline_segments=2,
                 beta_override=single_baseline_beta_for_matching_dt(2)),
        CaseSpec("md_2domain_n5", 2, 5, True),
    ]
    if include_ten or full:
        specs.extend(
            [
                CaseSpec("single_for_10domain_n5", 1, 5, False, single_baseline_segments=10,
                         beta_override=single_baseline_beta_for_matching_dt(10)),
                CaseSpec("md_10domain_n5", 10, 5, True),
            ]
        )
    if full:
        specs.extend(
            [
                CaseSpec("single_for_100domain_n5", 1, 5, False, single_baseline_segments=100,
                         beta_override=single_baseline_beta_for_matching_dt(100)),
                CaseSpec("md_100domain_n5", 100, 5, True),
                CaseSpec("md_100domain_segmented_material_n5", 100, 5, True, segmented_material=True),
            ]
        )
    return specs


def diagnostic_case_specs() -> List[CaseSpec]:
    return [
        CaseSpec("old_1domain_n5", 1, 5, False, single_baseline_segments=1),
        CaseSpec("md_1domain_n5", 1, 5, True),
        CaseSpec(
            "single_for_2domain_n5",
            1,
            5,
            False,
            single_baseline_segments=2,
            beta_override=single_baseline_beta_for_matching_dt(2),
            rhs_diagnostic=True,
        ),
        CaseSpec("md_2domain_n5", 2, 5, True, rhs_diagnostic=True),
    ]


def compare_suite(results: Dict[str, RunResult], include_ten: bool, full: bool) -> List[CompareResult]:
    checks: List[CompareResult] = []

    for name in ("md_1domain_n2", "md_1domain_n5", "md_2domain_n5"):
        checks.append(check_md_run(results[name], require_self_check=name.startswith("md_1domain")))

    for suffix in ("n2", "n5"):
        old_rows = read_state(results[f"old_1domain_{suffix}"].state_path)
        md_rows = read_state(results[f"md_1domain_{suffix}"].state_path)
        checks.append(
            compare_rows(
                f"1domain old/new {suffix}",
                old_rows,
                md_rows,
                abs_tol=1.0e-8,
                rel_tol=1.0e-6,
                key_fields=("step", "element", "localNode"),
            )
        )

    single2_rows = read_state(results["single_for_2domain_n5"].state_path)
    md2_rows = read_state(results["md_2domain_n5"].state_path)
    state_diff_path = write_state_outer_diff_report(results["single_for_2domain_n5"], results["md_2domain_n5"])
    checks.append(
        compare_rows(
            "2domain outer vs single",
            single2_rows,
            md2_rows,
            abs_tol=1.0e-7,
            rel_tol=1.0e-5,
            key_fields=("step", "localNode", "cx", "cy", "cz", "nx", "ny", "nz"),
            only_outer_b=True,
        )
    )
    checks[-1].message += f", report={state_diff_path}"
    checks.append(check_interface(results["md_2domain_n5"], abs_tol=1.0e-8, rel_tol=1.0e-6))

    if include_ten or full:
        checks.append(check_md_run(results["md_10domain_n5"]))

        single10_rows = read_state(results["single_for_10domain_n5"].state_path)
        md10_rows = read_state(results["md_10domain_n5"].state_path)
        checks.append(
            compare_rows(
                "10domain outer vs single",
                single10_rows,
                md10_rows,
                abs_tol=1.0e-7,
                rel_tol=1.0e-5,
                key_fields=("step", "localNode", "cx", "cy", "cz", "nx", "ny", "nz"),
                only_outer_b=True,
            )
        )
        checks.append(check_interface(results["md_10domain_n5"], abs_tol=1.0e-8, rel_tol=1.0e-6))

    if full:
        for name in ("md_100domain_n5", "md_100domain_segmented_material_n5"):
            checks.append(check_md_run(results[name]))

        single100_rows = read_state(results["single_for_100domain_n5"].state_path)
        md100_rows = read_state(results["md_100domain_n5"].state_path)
        checks.append(
            compare_rows(
                "100domain outer vs single",
                single100_rows,
                md100_rows,
                abs_tol=1.0e-7,
                rel_tol=1.0e-4,
                key_fields=("step", "localNode", "cx", "cy", "cz", "nx", "ny", "nz"),
                only_outer_b=True,
            )
        )
        checks.append(check_interface(results["md_100domain_n5"], abs_tol=1.0e-7, rel_tol=1.0e-5))
        checks.append(check_interface(results["md_100domain_segmented_material_n5"], abs_tol=1.0e-7, rel_tol=1.0e-5))

        material_lines = [
            value
            for key, value in results["md_100domain_segmented_material_n5"].metrics.items()
            if key.startswith("DomainMaterial[")
        ]
        material_ok = len(material_lines) == 100 and len(set(material_lines)) >= 5
        checks.append(
            CompareResult(
                "100domain segmented material output",
                material_ok,
                0.0,
                0.0,
                len(material_lines),
                message=f"domain_material_lines={len(material_lines)}, unique={len(set(material_lines))}",
            )
        )

    return checks


def compare_diagnostic_suite(results: Dict[str, RunResult]) -> List[CompareResult]:
    checks: List[CompareResult] = []

    checks.append(check_md_run(results["old_1domain_n5"]))
    checks.append(check_md_run(results["md_1domain_n5"], require_self_check=True))
    checks.append(check_md_run(results["single_for_2domain_n5"]))
    checks.append(check_md_run(results["md_2domain_n5"]))

    old_rows = read_state(results["old_1domain_n5"].state_path)
    md_rows = read_state(results["md_1domain_n5"].state_path)
    checks.append(
        compare_rows(
            "1domain old/new n5",
            old_rows,
            md_rows,
            abs_tol=1.0e-8,
            rel_tol=1.0e-6,
            key_fields=("step", "element", "localNode"),
        )
    )

    single2_rows = read_state(results["single_for_2domain_n5"].state_path)
    md2_rows = read_state(results["md_2domain_n5"].state_path)
    state_diff_path = write_state_outer_diff_report(results["single_for_2domain_n5"], results["md_2domain_n5"])
    checks.append(
        compare_rows(
            "2domain outer vs single",
            single2_rows,
            md2_rows,
            abs_tol=1.0e-7,
            rel_tol=1.0e-5,
                key_fields=("step", "localNode", "cx", "cy", "cz", "nx", "ny", "nz"),
            only_outer_b=True,
        )
    )
    checks[-1].message += f", report={state_diff_path}"
    checks.append(check_interface(results["md_2domain_n5"], abs_tol=1.0e-8, rel_tol=1.0e-6))
    checks.append(check_residual_probe(results["single_for_2domain_n5"], abs_tol=1.0e-10, rel_tol=1.0e-8))
    checks.append(check_residual_probe(results["md_2domain_n5"], abs_tol=1.0e-10, rel_tol=1.0e-8))
    checks.append(compare_rhs_breakdown(results["single_for_2domain_n5"], results["md_2domain_n5"]))
    return checks


def check_name_zh(name: str) -> str:
    mapping = {
        "old_1domain_n5 run": "旧单域路径 NStep=5 运行",
        "md_1domain_n2 run": "多域路径 1 域 NStep=2 运行",
        "md_1domain_n5 run": "多域路径 1 域 NStep=5 运行",
        "single_for_2domain_n5 run": "2 域对照用单域基准 NStep=5 运行",
        "md_2domain_n5 run": "多域路径 2 域 NStep=5 运行",
        "single_for_10domain_n5 run": "10 域对照用单域基准 NStep=5 运行",
        "md_10domain_n5 run": "多域路径 10 域 NStep=5 运行",
        "1domain old/new n2": "1 域新旧路径 NStep=2 对比",
        "1domain old/new n5": "1 域新旧路径 NStep=5 对比",
        "2domain outer vs single": "2 域外边界对单域基准",
        "10domain outer vs single": "10 域外边界对单域基准",
        "md_2domain_n5 interface": "2 域界面连续和平衡",
        "md_10domain_n5 interface": "10 域界面连续和平衡",
        "single_for_2domain_n5 residual probe": "2 域单域基准残差诊断",
        "md_2domain_n5 residual probe": "2 域多域路径残差诊断",
        "2domain RHS outer diagnostic": "2 域外边界 RHS 分解诊断",
        "all cases completed": "所有算例运行完成",
    }
    return mapping.get(name, name)


def write_summary(run_root: Path, results: Dict[str, RunResult], checks: Sequence[CompareResult]) -> None:
    summary = {
        "run_root": str(run_root),
        "cases": {
            name: {
                "case_dir": str(result.case_dir),
                "returncode": result.returncode,
                "metrics": result.metrics,
            }
            for name, result in results.items()
        },
        "checks": [
            {
                "name": check.name,
                "passed": check.passed,
                "max_abs": check.max_abs,
                "max_rel": check.max_rel,
                "compared_values": check.compared_values,
                "missing": check.missing,
                "message": check.message,
            }
            for check in checks
        ],
    }
    (run_root / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = ["# 多域验证汇总", ""]
    lines.append("## 算例")
    for name, result in results.items():
        lines.append(f"- `{name}`：返回码 `{result.returncode}`，目录 `{result.case_dir}`")
    lines.append("")
    lines.append("## 检查结果")
    for check in checks:
        status = "通过" if check.passed else "未通过"
        lines.append(f"- {status}: {check_name_zh(check.name)} ({check.message})")
    lines.append("")
    lines.append(f"运行目录：`{run_root}`")
    (run_root / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run DBEM multi-domain validation cases.")
    parser.add_argument("--ten", action="store_true", help="run 10-domain cases but skip 100-domain cases")
    parser.add_argument("--full", action="store_true", help="run 10-domain and 100-domain cases too")
    parser.add_argument(
        "--diagnose-two-domain",
        action="store_true",
        help="run only the 1-domain gate and 2-domain RHS/residual diagnostics",
    )
    parser.add_argument("--exe", help="path to DBEM1.exe; defaults to DBEM1/x64/Release/DBEM1.exe")
    parser.add_argument("--timeout", type=int, default=600, help="timeout per case in seconds")
    parser.add_argument("--run-root", help="custom validation run directory")
    args = parser.parse_args(argv)

    root = repo_root()
    exe = find_exe(root, args.exe)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = Path(args.run_root) if args.run_root else root / "validation_runs" / stamp
    run_root.mkdir(parents=True, exist_ok=True)

    include_ten = args.ten or args.full
    specs = diagnostic_case_specs() if args.diagnose_two_domain else case_specs(include_ten, args.full)
    results: Dict[str, RunResult] = {}
    for spec in specs:
        print(f"[run] {spec.name}")
        case_dir = prepare_case(run_root, root, exe, spec)
        result = run_case(case_dir, spec, args.timeout)
        results[spec.name] = result
        print(f"      returncode={result.returncode} dir={case_dir}")
        if not is_success_returncode(result.returncode) or not result.state_path.exists():
            print(f"      failed; see {case_dir / 'run.log'}")
            break

    expected = {spec.name for spec in specs}
    if set(results) != expected:
        checks = [
            CompareResult(
                "all cases completed",
                False,
                math.inf,
                math.inf,
                0,
                message=f"completed={sorted(results)}, expected={sorted(expected)}",
            )
        ]
        write_summary(run_root, results, checks)
        print(f"Summary: {run_root / 'summary.md'}")
        return 1

    checks = compare_diagnostic_suite(results) if args.diagnose_two_domain else compare_suite(results, include_ten, args.full)
    write_summary(run_root, results, checks)
    failed = [check for check in checks if not check.passed]
    for check in checks:
        status = "PASS" if check.passed else "FAIL"
        print(f"[{status}] {check.name}: {check.message}")
    print(f"Summary: {run_root / 'summary.md'}")
    return 0 if not failed else 2


if __name__ == "__main__":
    sys.exit(main())
