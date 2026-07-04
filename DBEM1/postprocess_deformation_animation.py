#!/usr/bin/env python3
"""Generate a deformation GIF from DBEM per-element Tecplot DAT files.

The C++ Tecplot output duplicates each element's 8 geometric nodes for display
and stores the corresponding discontinuous-node result values on those nodes.
This keeps the visual geometry closed while preserving validation_state.csv
semantics. This script reads the requested FEPOINT zone, normally ``deformed``,
and does not merge, average, smooth, or re-index coincident nodes.

The optional scale parameter applies an extra display-only displacement factor
on top of the DAT deformed geometry:

    Xplot = Xdeformed + (scale - 1) * ux
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


ZONE_RE = re.compile(r'ZONE\s+T="([^"]+)".*?\bN\s*=\s*(\d+).*?\bE\s*=\s*(\d+)', re.IGNORECASE)
STEP_RE = re.compile(r"_(\d+)\.dat$", re.IGNORECASE)


@dataclass
class FrameData:
    step: int
    path: Path
    coords: np.ndarray
    disp: np.ndarray
    disp_mag: np.ndarray
    conn: np.ndarray
    plot_coords: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a deformation animation from DBEM TecValueFile_*.dat output."
    )
    parser.add_argument("--input-dir", default="output", help="Directory containing TecValueFile_*.dat files.")
    parser.add_argument("--pattern", default="TecValueFile_*.dat", help="Input DAT filename pattern.")
    parser.add_argument("--output", default="output/deformation_animation.gif", help="Output GIF path.")
    parser.add_argument("--fps", type=float, default=4.0, help="Animation frames per second.")
    parser.add_argument("--scale", type=float, default=1.0, help="Extra displacement scale applied to deformed zone.")
    parser.add_argument("--zone", default="deformed", help="Tecplot zone title to animate.")
    parser.add_argument("--view", default="20,-60", help="3D view as elev,azim.")
    parser.add_argument("--dpi", type=int, default=140, help="Output image DPI.")
    parser.add_argument("--max-step", type=int, default=None, help="Maximum step to include. Defaults to NStep in validation_metrics.txt when available.")
    return parser.parse_args()


def step_from_path(path: Path) -> int:
    match = STEP_RE.search(path.name)
    if not match:
        return 0
    return int(match.group(1))


def resolve_input_dir(input_dir: str) -> Path:
    path = Path(input_dir)
    if path.exists():
        return path
    script_relative = Path(__file__).resolve().parent / input_dir
    if script_relative.exists():
        return script_relative
    return path


def parse_zone_header(line: str) -> tuple[str, int, int] | None:
    match = ZONE_RE.search(line)
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3))


def read_tecplot_zone(path: Path, zone_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = path.read_text(errors="ignore").splitlines()
    wanted = zone_name.lower()

    line_id = 0
    while line_id < len(lines):
        parsed = parse_zone_header(lines[line_id])
        if not parsed:
            line_id += 1
            continue

        title, node_count, element_count = parsed
        data_begin = line_id + 1
        data_end = data_begin + node_count
        conn_end = data_end + element_count
        if conn_end > len(lines):
            raise ValueError(f"{path}: zone {title!r} is truncated.")

        if title.lower() == wanted:
            node_rows = []
            for raw in lines[data_begin:data_end]:
                values = [float(item) for item in raw.split()]
                if len(values) < 6:
                    raise ValueError(f"{path}: node row has fewer than 6 columns: {raw!r}")
                node_rows.append(values)

            conn_rows = []
            for raw in lines[data_end:conn_end]:
                values = [int(item) - 1 for item in raw.split()]
                if len(values) != 4:
                    raise ValueError(f"{path}: only quadrilateral connectivity is supported: {raw!r}")
                conn_rows.append(values)

            data = np.asarray(node_rows, dtype=float)
            coords = data[:, 0:3]
            disp = data[:, 3:6]
            conn = np.asarray(conn_rows, dtype=np.int64)
            if not np.all(np.isfinite(data[:, 0:6])):
                raise ValueError(f"{path}: zone {title!r} contains non-finite coordinates or displacement values.")
            if conn.size and (int(conn.min()) < 0 or int(conn.max()) >= node_count):
                raise ValueError(f"{path}: zone {title!r} connectivity references nodes outside 1..{node_count}.")
            return coords, disp, conn

        line_id = conn_end

    raise ValueError(f"{path}: zone {zone_name!r} was not found.")


def parse_view(value: str) -> tuple[float, float]:
    parts = [item.strip() for item in value.split(",")]
    if len(parts) != 2:
        raise ValueError("--view must be formatted as elev,azim, for example 20,-60.")
    return float(parts[0]), float(parts[1])


def read_default_max_step(input_dir: Path) -> int | None:
    metrics = input_dir / "validation_metrics.txt"
    if not metrics.exists():
        return None
    for line in metrics.read_text(errors="ignore").splitlines():
        if line.startswith("NStep="):
            try:
                return int(line.split("=", 1)[1].strip())
            except ValueError:
                return None
    return None


def load_frames(input_dir: Path, pattern: str, zone: str, scale: float, max_step: int | None) -> list[FrameData]:
    files = sorted(input_dir.glob(pattern), key=step_from_path)
    if max_step is not None:
        files = [path for path in files if step_from_path(path) <= max_step]
    if not files:
        raise FileNotFoundError(f"No files matched {input_dir / pattern}.")

    frames: list[FrameData] = []
    expected_conn: np.ndarray | None = None
    for path in files:
        coords, disp, conn = read_tecplot_zone(path, zone)
        plot_coords = coords + (scale - 1.0) * disp
        disp_mag = np.linalg.norm(disp, axis=1)

        if expected_conn is None:
            expected_conn = conn
        elif conn.shape != expected_conn.shape or not np.array_equal(conn, expected_conn):
            raise ValueError(f"{path}: connectivity differs from the first frame.")

        frames.append(
            FrameData(
                step=step_from_path(path),
                path=path,
                coords=coords,
                disp=disp,
                disp_mag=disp_mag,
                conn=conn,
                plot_coords=plot_coords,
            )
        )
    return frames


def equal_axes_limits(all_coords: np.ndarray) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    mins = np.nanmin(all_coords, axis=0)
    maxs = np.nanmax(all_coords, axis=0)
    center = 0.5 * (mins + maxs)
    span = float(np.max(maxs - mins))
    if not math.isfinite(span) or span <= 0.0:
        span = 1.0
    half = 0.55 * span
    return (
        (center[0] - half, center[0] + half),
        (center[1] - half, center[1] + half),
        (center[2] - half, center[2] + half),
    )


def render_frame(
    frame: FrameData,
    axes_limits: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    color_limits: tuple[float, float],
    view: tuple[float, float],
    scale: float,
    dpi: int,
) -> np.ndarray:
    fig = plt.figure(figsize=(7.2, 5.8), dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")

    faces = frame.plot_coords[frame.conn]
    face_values = frame.disp_mag[frame.conn].mean(axis=1)
    vmin, vmax = color_limits
    if vmax <= vmin:
        vmax = vmin + 1.0
    colors = cm.viridis((face_values - vmin) / (vmax - vmin))

    collection = Poly3DCollection(faces, facecolors=colors, edgecolors=(0.2, 0.2, 0.2, 0.35), linewidths=0.25)
    ax.add_collection3d(collection)

    ax.set_xlim(*axes_limits[0])
    ax.set_ylim(*axes_limits[1])
    ax.set_zlim(*axes_limits[2])
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=view[0], azim=view[1])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(f"Step {frame.step}   scale={scale:g}")

    scalar_map = cm.ScalarMappable(cmap="viridis")
    scalar_map.set_clim(vmin, vmax)
    cbar = fig.colorbar(scalar_map, ax=ax, shrink=0.72, pad=0.08)
    cbar.set_label("|u|")

    fig.tight_layout()
    fig.canvas.draw()
    image = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
    plt.close(fig)
    return image


def save_gif(images: list[np.ndarray], output_path: Path, fps: float) -> None:
    from PIL import Image

    if fps <= 0.0:
        raise ValueError("--fps must be positive.")
    duration_ms = int(round(1000.0 / fps))
    pil_images = [Image.fromarray(image) for image in images]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pil_images[0].save(
        output_path,
        save_all=True,
        append_images=pil_images[1:],
        duration=duration_ms,
        loop=0,
    )


def main() -> int:
    args = parse_args()
    input_dir = resolve_input_dir(args.input_dir)
    output_path = Path(args.output)
    if not output_path.is_absolute() and not output_path.parent.exists():
        output_path = input_dir.parent / output_path

    view = parse_view(args.view)
    max_step = args.max_step
    if max_step is None:
        max_step = read_default_max_step(input_dir)
    frames = load_frames(input_dir, args.pattern, args.zone, args.scale, max_step)
    all_plot_coords = np.vstack([frame.plot_coords for frame in frames])
    axes_limits = equal_axes_limits(all_plot_coords)
    disp_values = np.concatenate([frame.disp_mag for frame in frames])
    color_limits = (float(np.nanmin(disp_values)), float(np.nanmax(disp_values)))

    images = [
        render_frame(frame, axes_limits, color_limits, view, args.scale, args.dpi)
        for frame in frames
    ]
    save_gif(images, output_path, args.fps)

    print(f"frames={len(frames)}")
    print(f"zone={args.zone}")
    print("topology=per-element-duplicated-geometry-node")
    print(f"scale={args.scale:g}")
    if max_step is not None:
        print(f"max_step={max_step}")
    print(f"output={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
