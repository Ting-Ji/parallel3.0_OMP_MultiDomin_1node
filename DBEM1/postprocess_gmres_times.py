#!/usr/bin/env python
"""Collect GMRES timing files from batch cases into one table.

Default input:
    output/case_*/GMRESTime.txt

Default outputs:
    output/gmres_time_summary.csv
    output/gmres_time_summary.txt
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from statistics import mean


CASE_RE = re.compile(r"case_(\d+)$", re.IGNORECASE)
TIME_PATTERNS = {
    "assembly_time": re.compile(r"(?:组装时间|装配时间|assembly)\s*[:：]\s*([-+0-9.eE]+)", re.IGNORECASE),
    "solution_time": re.compile(r"(?:求解时间|solution)\s*[:：]\s*([-+0-9.eE]+)", re.IGNORECASE),
    "gmres_time": re.compile(r"(?:GMRES求解|GMRES\s*solve|gmres)\s*[:：]\s*([-+0-9.eE]+)", re.IGNORECASE),
}


def read_text(path: Path) -> str:
    for encoding in ("utf-8-sig", "gbk", "mbcs"):
        try:
            return path.read_text(encoding=encoding)
        except UnicodeDecodeError:
            continue
    return path.read_text(errors="replace")


def parse_time_file(path: Path) -> dict[str, float | int | str]:
    case_match = CASE_RE.search(path.parent.name)
    if not case_match:
        raise ValueError(f"Cannot parse case id from directory name: {path.parent.name}")

    text = read_text(path)
    row: dict[str, float | int | str] = {
        "case_id": int(case_match.group(1)),
        "case_name": path.parent.name,
        "time_file": str(path),
    }

    for key, pattern in TIME_PATTERNS.items():
        match = pattern.search(text)
        if not match:
            raise ValueError(f"Cannot parse {key} from {path}")
        row[key] = float(match.group(1))

    row["total_recorded_time"] = (
        float(row["assembly_time"])
        + float(row["solution_time"])
        + float(row["gmres_time"])
    )
    return row


def collect_rows(output_dir: Path) -> list[dict[str, float | int | str]]:
    files = sorted(
        output_dir.glob("case_*/GMRESTime.txt"),
        key=lambda p: int(CASE_RE.search(p.parent.name).group(1)) if CASE_RE.search(p.parent.name) else 10**12,
    )
    rows = [parse_time_file(path) for path in files]
    if not rows:
        raise FileNotFoundError(f"No GMRESTime.txt files found under {output_dir / 'case_*'}")

    first_assembly = float(rows[0]["assembly_time"])
    cumulative_time = 0.0
    for row in rows:
        assembly_time = float(row["assembly_time"])
        total_recorded_time = float(row["total_recorded_time"])
        cumulative_time += total_recorded_time
        row["assembly_speedup_vs_case_001"] = first_assembly / assembly_time if assembly_time > 0 else ""
        row["assembly_time_saved_vs_case_001"] = first_assembly - assembly_time
        row["cumulative_recorded_time"] = cumulative_time

    return rows


def write_csv(rows: list[dict[str, float | int | str]], path: Path) -> None:
    fieldnames = [
        "case_id",
        "case_name",
        "assembly_time",
        "solution_time",
        "gmres_time",
        "total_recorded_time",
        "assembly_speedup_vs_case_001",
        "assembly_time_saved_vs_case_001",
        "cumulative_recorded_time",
        "time_file",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(rows: list[dict[str, float | int | str]], path: Path) -> None:
    first = rows[0]
    rest = rows[1:]
    assembly_values = [float(row["assembly_time"]) for row in rows]
    solution_values = [float(row["solution_time"]) for row in rows]
    gmres_values = [float(row["gmres_time"]) for row in rows]
    total_values = [float(row["total_recorded_time"]) for row in rows]

    first_assembly = float(first["assembly_time"])
    rest_assembly_mean = mean(float(row["assembly_time"]) for row in rest) if rest else 0.0
    speedup = first_assembly / rest_assembly_mean if rest_assembly_mean > 0 else 0.0

    lines = [
        "GMRES batch timing summary",
        f"case_count: {len(rows)}",
        f"first_case: {first['case_name']}",
        f"first_case_assembly_time: {first_assembly:.6f}",
        f"subsequent_case_mean_assembly_time: {rest_assembly_mean:.6f}",
        f"estimated_assembly_speedup_first_vs_subsequent_mean: {speedup:.6f}",
        f"mean_assembly_time_all_cases: {mean(assembly_values):.6f}",
        f"mean_solution_time_all_cases: {mean(solution_values):.6f}",
        f"mean_gmres_time_all_cases: {mean(gmres_values):.6f}",
        f"mean_total_recorded_time_all_cases: {mean(total_values):.6f}",
        f"sum_total_recorded_time_all_cases: {sum(total_values):.6f}",
        f"min_assembly_time: {min(assembly_values):.6f}",
        f"max_assembly_time: {max(assembly_values):.6f}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Collect GMRES batch timing files into CSV and text summaries.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "output",
        help="Directory containing case_*/GMRESTime.txt. Default: DBEM1/output",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="CSV output path. Default: <output-dir>/gmres_time_summary.csv",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="Text summary output path. Default: <output-dir>/gmres_time_summary.txt",
    )
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    csv_path = args.csv.resolve() if args.csv else output_dir / "gmres_time_summary.csv"
    summary_path = args.summary.resolve() if args.summary else output_dir / "gmres_time_summary.txt"

    rows = collect_rows(output_dir)
    write_csv(rows, csv_path)
    write_summary(rows, summary_path)

    print(f"Collected {len(rows)} GMRES timing files.")
    print(f"CSV: {csv_path}")
    print(f"Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
