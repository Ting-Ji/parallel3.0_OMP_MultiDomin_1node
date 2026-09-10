#!/usr/bin/env python3
"""Audit existing DBEM1 output case directories without modifying them.

The structured validation files are treated as the source of truth.  A case
passes only when the solver reports a complete, converged run, the matrix and
interface invariants pass, every expected Tecplot step exists, and the full
validation state is finite and complete.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path


TEC_PATTERNS = {
    "1node": re.compile(r"^TecValueFile_1node_(\d+)\.dat$"),
    "realsingle": re.compile(r"^TecValueFile_realsingle_(\d+)\.dat$"),
    "legacy": re.compile(r"^TecValueFile_(\d+)\.dat$"),
}
REQUIRED_STATE_COLUMNS = (
    "step",
    "element",
    "localNode",
    "domain",
    "surfaceType",
    "bcid",
    "x",
    "y",
    "z",
    "ux",
    "uy",
    "uz",
    "tx",
    "ty",
    "tz",
)


@dataclass
class StateStats:
    rows: int = 0
    min_step: int | None = None
    max_step: int | None = None
    nonfinite: int = 0
    max_abs_u: float = 0.0
    max_abs_t: float = 0.0


@dataclass
class CaseResult:
    name: str
    status: str
    reasons: list[str]
    nstep: int | None
    solved_steps: int | None
    state_rows: int
    expected_state_rows: int | None
    tec_steps: int
    max_abs_u: float
    max_abs_t: float
    tec_family: str | None


def read_metrics(path: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            metrics[key.strip()] = value.strip()
    return metrics


def metric_int(metrics: dict[str, str], key: str) -> int | None:
    value = metrics.get(key)
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def scan_state(path: Path) -> tuple[StateStats, list[str]]:
    stats = StateStats()
    reasons: list[str] = []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as stream:
        reader = csv.reader(stream)
        try:
            header = tuple(next(reader))
        except StopIteration:
            return stats, ["validation_state.csv is empty"]
        missing_columns = [name for name in REQUIRED_STATE_COLUMNS if name not in header]
        if missing_columns:
            return stats, [f"validation_state.csv is missing columns {missing_columns}"]
        indices = {name: header.index(name) for name in REQUIRED_STATE_COLUMNS}

        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                reasons.append(f"state row {line_number} has {len(row)} columns (expected {len(header)})")
                continue
            stats.rows += 1
            try:
                step = int(row[indices["step"]])
            except ValueError:
                reasons.append(f"state row {line_number} has an invalid step")
                continue
            if stats.min_step is None or step < stats.min_step:
                stats.min_step = step
            if stats.max_step is None or step > stats.max_step:
                stats.max_step = step

            for field in ("ux", "uy", "uz", "tx", "ty", "tz"):
                index = indices[field]
                try:
                    value = float(row[index])
                except ValueError:
                    stats.nonfinite += 1
                    continue
                if not math.isfinite(value):
                    stats.nonfinite += 1
                    continue
                absolute = abs(value)
                if field.startswith("u"):
                    stats.max_abs_u = max(stats.max_abs_u, absolute)
                else:
                    stats.max_abs_t = max(stats.max_abs_t, absolute)
    return stats, reasons


def last_nonempty_line(path: Path) -> str:
    with path.open("rb") as stream:
        stream.seek(0, 2)
        position = stream.tell()
        data = b""
        while position > 0 and data.count(b"\n") < 2:
            chunk_size = min(8192, position)
            position -= chunk_size
            stream.seek(position)
            data = stream.read(chunk_size) + data
    lines = [line.strip() for line in data.decode("utf-8", errors="replace").splitlines() if line.strip()]
    return lines[-1] if lines else ""


def audit_case(case_dir: Path, max_abs_result: float) -> CaseResult:
    reasons: list[str] = []
    metrics_path = case_dir / "validation_metrics.txt"
    state_path = case_dir / "validation_state.csv"
    audit_path = case_dir / "interface_transfer_audit.csv"

    metrics: dict[str, str] = {}
    if not metrics_path.is_file():
        reasons.append("missing validation_metrics.txt")
    else:
        metrics = read_metrics(metrics_path)

    nstep = metric_int(metrics, "NStep")
    node_count = metric_int(metrics, "NodeCount")
    solved_steps = metric_int(metrics, "MultiDomainSolvedSteps")
    required_metrics = {
        "ValidationOutput": "1",
        "MultiDomainMetricsValid": "1",
        "MultiDomainFinalFlag": "0",
        "MultiDomainMatrixStructureValid": "1",
        "MultiDomainMatrixZeroRows": "0",
        "MultiDomainMatrixZeroCols": "0",
        "MultiDomainMatrixNonFiniteBlocks": "0",
        "MultiDomainInterfaceDofInvariantValid": "1",
        "MultiDomainInterfaceDofInvariantFailures": "0",
    }
    for key, expected in required_metrics.items():
        actual = metrics.get(key)
        if actual != expected:
            reasons.append(f"{key}={actual or 'MISSING'} (expected {expected})")
    if nstep is None:
        reasons.append("NStep is missing or invalid")
    if nstep is not None and solved_steps != nstep:
        reasons.append(f"MultiDomainSolvedSteps={solved_steps or 'MISSING'} (expected {nstep})")

    tec_families: dict[str, list[int]] = {name: [] for name in TEC_PATTERNS}
    for path in case_dir.iterdir():
        if not path.is_file():
            continue
        for family, pattern in TEC_PATTERNS.items():
            if match := pattern.match(path.name):
                tec_families[family].append(int(match.group(1)))
                break
    for indices in tec_families.values():
        indices.sort()
    tec_family = None
    tec_indices: list[int] = []
    if nstep is not None:
        expected_indices = list(range(nstep + 1))
        complete_families = [name for name, indices in tec_families.items() if indices == expected_indices]
        if complete_families:
            tec_family = complete_families[0]
            tec_indices = tec_families[tec_family]
        else:
            tec_family, tec_indices = max(tec_families.items(), key=lambda item: len(item[1]))
            missing = sorted(set(expected_indices) - set(tec_indices))
            extra = sorted(set(tec_indices) - set(expected_indices))
            reasons.append(
                f"Tecplot steps incomplete in all families; best={tec_family}, count={len(tec_indices)}, "
                f"missing={missing[:8]}, extra={extra[:8]}"
            )
    else:
        tec_family, tec_indices = max(tec_families.items(), key=lambda item: len(item[1]))
        if not tec_indices:
            reasons.append("no recognized TecValueFile_*.dat files")

    stats = StateStats()
    state_reasons: list[str] = []
    if not state_path.is_file():
        reasons.append("missing validation_state.csv")
    else:
        stats, state_reasons = scan_state(state_path)
        reasons.extend(state_reasons)
    expected_rows = None
    if nstep is not None and node_count is not None:
        expected_rows = (nstep + 1) * node_count
        if stats.rows != expected_rows:
            reasons.append(f"state rows={stats.rows} (expected {expected_rows})")
    if nstep is not None and (stats.min_step != 0 or stats.max_step != nstep):
        reasons.append(f"state step range={stats.min_step}..{stats.max_step} (expected 0..{nstep})")
    if stats.nonfinite:
        reasons.append(f"state contains {stats.nonfinite} non-finite/invalid result values")
    if max(stats.max_abs_u, stats.max_abs_t) > max_abs_result:
        reasons.append(
            f"catastrophic result magnitude exceeds {max_abs_result:.9g} "
            f"(max|u|={stats.max_abs_u:.9g}, max|t|={stats.max_abs_t:.9g})"
        )

    if not audit_path.is_file() or audit_path.stat().st_size == 0:
        reasons.append("missing or empty interface_transfer_audit.csv")
    elif nstep is not None:
        last_fields = last_nonempty_line(audit_path).split(",")
        try:
            audit_last_step = int(last_fields[1])
        except (IndexError, ValueError):
            reasons.append("interface audit final row is invalid")
        else:
            if audit_last_step != nstep:
                reasons.append(f"interface audit ends at step {audit_last_step} (expected {nstep})")

    return CaseResult(
        name=case_dir.name,
        status="PASS" if not reasons else "FAIL",
        reasons=reasons,
        nstep=nstep,
        solved_steps=solved_steps,
        state_rows=stats.rows,
        expected_state_rows=expected_rows,
        tec_steps=len(tec_indices),
        max_abs_u=stats.max_abs_u,
        max_abs_t=stats.max_abs_t,
        tec_family=tec_family,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="DBEM1/output directory")
    parser.add_argument("--csv", type=Path, help="optional CSV report path")
    parser.add_argument(
        "--max-abs-result",
        type=float,
        default=1.0e8,
        help="fail cases whose displacement or traction exceeds this divergence guard (default: 1e8)",
    )
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    if not output_dir.is_dir():
        parser.error(f"not a directory: {output_dir}")
    case_dirs = sorted((path for path in output_dir.iterdir() if path.is_dir()), key=lambda path: path.name)
    results: list[CaseResult] = []
    for case_dir in case_dirs:
        result = audit_case(case_dir, args.max_abs_result)
        results.append(result)
        reason = "; ".join(result.reasons) if result.reasons else "all checks passed"
        print(
            f"{result.status}\t{result.name}\tNStep={result.nstep}\t"
            f"state_rows={result.state_rows}\tTecSteps={result.tec_steps}({result.tec_family})\t"
            f"max|u|={result.max_abs_u:.9g}\tmax|t|={result.max_abs_t:.9g}\t{reason}",
            flush=True,
        )

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    "status",
                    "case",
                    "nstep",
                    "solved_steps",
                    "state_rows",
                    "expected_state_rows",
                    "tec_steps",
                    "tec_family",
                    "max_abs_u",
                    "max_abs_t",
                    "reasons",
                ]
            )
            for result in results:
                writer.writerow(
                    [
                        result.status,
                        result.name,
                        result.nstep,
                        result.solved_steps,
                        result.state_rows,
                        result.expected_state_rows,
                        result.tec_steps,
                        result.tec_family,
                        f"{result.max_abs_u:.17g}",
                        f"{result.max_abs_t:.17g}",
                        "; ".join(result.reasons),
                    ]
                )

    failures = sum(result.status == "FAIL" for result in results)
    print(f"SUMMARY\tcases={len(results)}\tpassed={len(results) - failures}\tfailed={failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
