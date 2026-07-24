#!/usr/bin/env python3
"""Generate a conforming multi-domain rectangular rod case for DBEM1.

The current multi-domain reader divides elements into domains by element order:
domain d owns a contiguous block of EleNum / DomainCount elements.  Each domain
contains six faces, and each face is subdivided into nx*ny, nx*nz, or ny*nz
8-node quadrilateral boundary elements depending on the face orientation.

Element order within a domain (face-major):
  1. left  cap/interface face  (y-z plane): ny * nz elements, loop y then z
  2. right cap/interface face  (y-z plane): ny * nz elements, loop y then z
  3. y_min side face           (x-z plane): nx * nz elements, loop x then z
  4. y_max side face           (x-z plane): nx * nz elements, loop x then z
  5. z_min side face           (x-y plane): nx * ny elements, loop x then y
  6. z_max side face           (x-y plane): nx * ny elements, loop x then y

Total per domain: 2*ny*nz + 2*nx*nz + 2*nx*ny elements.

Adjacent domains share interface pairs on their touching y-z faces.  With
subdivision, each face has ny*nz elements, paired in the same y/z loop order.

The .bd format remains the legacy element-level format: one BCID per element,
then one triplet per element for every time step.  DBEM1 expands each triplet to
the eight discontinuous physical nodes internally.

The script also writes a matching single-domain outer-surface case.  That model
uses the same total rod geometry and boundary condition convention, but it does
not include any internal interface faces.  It is compatible with the legacy
constant-element solver path: the model still stores 8-node geometric
quadrilateral elements, while the solver creates one physical unknown node at
each element center internally (NodeNum = EleNum).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

Point = Tuple[float, float, float]
Element = Tuple[int, int, int, int, int, int, int, int]


def vec_sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vec_dot(a: Point, b: Point) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vec_cross(a: Point, b: Point) -> Point:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vec_norm(a: Point) -> float:
    return math.sqrt(vec_dot(a, a))


def vec_normalize(a: Point) -> Point:
    length = vec_norm(a)
    if length <= 0.0:
        raise ValueError("zero-length vector")
    return (a[0] / length, a[1] / length, a[2] / length)


def point_distance(a: Point, b: Point) -> float:
    return vec_norm(vec_sub(a, b))


def quad_basis(points: Sequence[Point], element: Element) -> Tuple[Point, Point, Point]:
    """Return the local t1, t2 and normal basis used by DSquareElement."""
    p0 = points[element[0] - 1]
    p1 = points[element[1] - 1]
    p3 = points[element[3] - 1]
    t1 = vec_normalize(vec_sub(p1, p0))
    raw_t2 = vec_normalize(vec_sub(p3, p0))
    normal = vec_normalize(vec_cross(t1, raw_t2))
    t2 = vec_cross(normal, t1)
    return t1, t2, normal


def global_to_local_components(points: Sequence[Point],
                               elements: Sequence[Element],
                               element_id: int,
                               global_value: Point) -> Point:
    t1, t2, normal = quad_basis(points, elements[element_id - 1])
    return (
        vec_dot(t1, global_value),
        vec_dot(t2, global_value),
        vec_dot(normal, global_value),
    )


def midpoint(a: Point, b: Point) -> Point:
    return ((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, (a[2] + b[2]) * 0.5)


def format_float(value: float) -> str:
    if abs(value) < 1.0e-14:
        value = 0.0
    return f"{value:.12g}"


class PointTable:
    def __init__(self, decimals: int = 12) -> None:
        self.decimals = decimals
        self.points: List[Point] = []
        self.ids: Dict[Tuple[float, float, float], int] = {}

    def add(self, p: Point) -> int:
        key = tuple(round(v, self.decimals) for v in p)
        point_id = self.ids.get(key)
        if point_id is None:
            point_id = len(self.points) + 1
            self.ids[key] = point_id
            self.points.append(key)  # type: ignore[arg-type]
        return point_id


def add_quad(table: PointTable, elements: List[Element], corners: Sequence[Point]) -> int:
    """Add one 8-node quadrilateral and return its 1-based element id.

    Local node order follows the common serendipity convention used by the
    solver: four corners followed by four edge midpoints.
    """
    if len(corners) != 4:
        raise ValueError("a quadrilateral requires four corner points")
    p0, p1, p2, p3 = corners
    nodes = (
        table.add(p0),
        table.add(p1),
        table.add(p2),
        table.add(p3),
        table.add(midpoint(p0, p1)),
        table.add(midpoint(p1, p2)),
        table.add(midpoint(p2, p3)),
        table.add(midpoint(p3, p0)),
    )
    elements.append(nodes)
    return len(elements)


def build_domain(table: PointTable,
                 elements: List[Element],
                 x0: float,
                 x1: float,
                 half_width: float,
                 half_height: float,
                 nx: int = 1,
                 ny: int = 1,
                 nz: int = 1) -> Dict[str, List[int]]:
    """Build one rectangular-prism domain with subdivided faces.

    Parameters
    ----------
    nx : int
        Subdivisions along the rod axis (x-direction) on side faces.
    ny : int
        Subdivisions along the y-direction on all faces.
    nz : int
        Subdivisions along the z-direction on all faces.

    Returns
    -------
    Dict[str, List[int]]
        Face name -> list of 1-based element IDs, in the order described in the
        module docstring.  This ordering is deterministic and matches the solver's
        assumption that adjacent domains' touch-faces are element-wise aligned.
    """
    y0, y1 = -half_width, half_width
    z0, z1 = -half_height, half_height

    dx = (x1 - x0) / max(nx, 1)
    dy = (y1 - y0) / max(ny, 1)
    dz = (z1 - z0) / max(nz, 1)

    ids: Dict[str, List[int]] = {
        "left": [],
        "right": [],
        "y_min": [],
        "y_max": [],
        "z_min": [],
        "z_max": [],
    }

    # ------------------------------------------------------------------
    # Left face  (x = x0, normal = -x, y-z plane): ny * nz elements
    # Winding (a, b, c, d) 鈥?counterclockwise viewed from -x.
    # ------------------------------------------------------------------
    for iy in range(ny):
        for iz in range(nz):
            cy0 = y0 + iy * dy
            cy1 = y0 + (iy + 1) * dy
            cz0 = z0 + iz * dz
            cz1 = z0 + (iz + 1) * dz
            a = (x0, cy0, cz0)
            b = (x0, cy1, cz0)
            c = (x0, cy1, cz1)
            d = (x0, cy0, cz1)
            ids["left"].append(add_quad(table, elements, (a, d, c, b)))

    # ------------------------------------------------------------------
    # Right face (x = x1, normal = +x, y-z plane): ny * nz elements
    # Winding (e, f, g, h) 鈥?counterclockwise viewed from +x.
    # Uses the same (iy, iz) loop order as the left face so that interface
    # pairing is element-wise correct.
    # ------------------------------------------------------------------
    for iy in range(ny):
        for iz in range(nz):
            cy0 = y0 + iy * dy
            cy1 = y0 + (iy + 1) * dy
            cz0 = z0 + iz * dz
            cz1 = z0 + (iz + 1) * dz
            e = (x1, cy0, cz0)
            f = (x1, cy1, cz0)
            g = (x1, cy1, cz1)
            h = (x1, cy0, cz1)
            ids["right"].append(add_quad(table, elements, (e, f, g, h)))

    # ------------------------------------------------------------------
    # y_min face (y = -half_width, normal = -y, x-z plane): nx * nz elements
    # Winding (a, e, h, d) 鈥?counterclockwise viewed from -y.
    # ------------------------------------------------------------------
    for ix in range(nx):
        for iz in range(nz):
            cx0 = x0 + ix * dx
            cx1 = x0 + (ix + 1) * dx
            cz0 = z0 + iz * dz
            cz1 = z0 + (iz + 1) * dz
            a = (cx0, y0, cz0)
            e = (cx1, y0, cz0)
            h = (cx1, y0, cz1)
            d = (cx0, y0, cz1)
            ids["y_min"].append(add_quad(table, elements, (a, e, h, d)))

    # ------------------------------------------------------------------
    # y_max face (y = +half_width, normal = +y, x-z plane): nx * nz elements
    # Winding (b, f, g, c) 鈥?counterclockwise viewed from +y.
    # ------------------------------------------------------------------
    for ix in range(nx):
        for iz in range(nz):
            cx0 = x0 + ix * dx
            cx1 = x0 + (ix + 1) * dx
            cz0 = z0 + iz * dz
            cz1 = z0 + (iz + 1) * dz
            b = (cx0, y1, cz0)
            f = (cx1, y1, cz0)
            g = (cx1, y1, cz1)
            c = (cx0, y1, cz1)
            ids["y_max"].append(add_quad(table, elements, (b, c, g, f)))

    # ------------------------------------------------------------------
    # z_min face (z = -half_height, normal = -z, x-y plane): nx * ny elements
    # Winding (a, e, f, b) 鈥?counterclockwise viewed from -z.
    # ------------------------------------------------------------------
    for ix in range(nx):
        for iy in range(ny):
            cx0 = x0 + ix * dx
            cx1 = x0 + (ix + 1) * dx
            cy0 = y0 + iy * dy
            cy1 = y0 + (iy + 1) * dy
            a = (cx0, cy0, z0)
            e = (cx1, cy0, z0)
            f = (cx1, cy1, z0)
            b = (cx0, cy1, z0)
            ids["z_min"].append(add_quad(table, elements, (a, b, f, e)))

    # ------------------------------------------------------------------
    # z_max face (z = +half_height, normal = +z, x-y plane): nx * ny elements
    # Winding (d, h, g, c) 鈥?counterclockwise viewed from +z.
    # ------------------------------------------------------------------
    for ix in range(nx):
        for iy in range(ny):
            cx0 = x0 + ix * dx
            cx1 = x0 + (ix + 1) * dx
            cy0 = y0 + iy * dy
            cy1 = y0 + (iy + 1) * dy
            d = (cx0, cy0, z1)
            h = (cx1, cy0, z1)
            g = (cx1, cy1, z1)
            c = (cx0, cy1, z1)
            ids["z_max"].append(add_quad(table, elements, (d, h, g, c)))

    return ids


def elems_per_domain(nx: int, ny: int, nz: int) -> int:
    """Return the number of boundary elements generated per domain."""
    return 2 * ny * nz + 2 * nx * nz + 2 * nx * ny


def _match_element_nodes(points: Sequence[Point],
                         elements: Sequence[Element],
                         ele_a: int,
                         ele_b: int,
                         tolerance: float = 1.0e-8) -> Tuple[List[int], List[int]]:
    node_points_a = [points[node_id - 1] for node_id in elements[ele_a - 1]]
    node_points_b = [points[node_id - 1] for node_id in elements[ele_b - 1]]
    scale = max(1.0, *(abs(c) for pt in node_points_a + node_points_b for c in pt))
    tol = tolerance * scale
    local_b_for_a = [-1] * 8
    local_a_for_b = [-1] * 8
    used_b: Set[int] = set()
    for local_a, point_a in enumerate(node_points_a):
        best_b = -1
        best_dist = float("inf")
        for local_b, point_b in enumerate(node_points_b):
            if local_b in used_b:
                continue
            dist = point_distance(point_a, point_b)
            if dist < best_dist:
                best_dist = dist
                best_b = local_b
        if best_b < 0 or best_dist > tol:
            raise RuntimeError(
                f"interface node mismatch: eleA={ele_a} localA={local_a} "
                f"eleB={ele_b} bestLocalB={best_b} distance={best_dist:.6e} tol={tol:.6e}"
            )
        local_b_for_a[local_a] = best_b
        local_a_for_b[best_b] = local_a
        used_b.add(best_b)
    if any(v < 0 for v in local_a_for_b):
        raise RuntimeError(f"interface node mapping is not one-to-one: eleA={ele_a} eleB={ele_b}")
    return local_b_for_a, local_a_for_b


def validate_rod_geometry(points: Sequence[Point],
                          elements: Sequence[Element],
                          faces_by_domain: Sequence[Dict[str, List[int]]],
                          interfaces: Sequence[Tuple[int, int, int, int, int]]) -> None:
    expected_normals: Dict[str, Point] = {
        "left": (-1.0, 0.0, 0.0),
        "right": (1.0, 0.0, 0.0),
        "y_min": (0.0, -1.0, 0.0),
        "y_max": (0.0, 1.0, 0.0),
        "z_min": (0.0, 0.0, -1.0),
        "z_max": (0.0, 0.0, 1.0),
    }
    for domain_id, faces in enumerate(faces_by_domain, start=1):
        for face_name, expected in expected_normals.items():
            for ele_id in faces[face_name]:
                normal = quad_basis(points, elements[ele_id - 1])[2]
                alignment = vec_dot(normal, expected)
                if alignment < 1.0 - 1.0e-10:
                    raise RuntimeError(
                        f"bad outward normal: domain={domain_id} face={face_name} ele={ele_id} "
                        f"normal={normal} expected={expected} dot={alignment:.12g}"
                    )
    for _, ele_a, _, ele_b, _ in interfaces:
        _match_element_nodes(points, elements, ele_a, ele_b)


def validate_single_domain_outer_surface(points: Sequence[Point],
                                         elements: Sequence[Element],
                                         length: float,
                                         width: float,
                                         height: float,
                                         tolerance: float = 1.0e-9) -> None:
    """Validate that a single-domain model contains only the rod exterior faces."""
    half_length = 0.5 * length
    half_width = 0.5 * width
    half_height = 0.5 * height
    scale = max(1.0, length, width, height)
    tol = tolerance * scale
    counts = {"left": 0, "right": 0, "y_min": 0, "y_max": 0, "z_min": 0, "z_max": 0}

    for ele_id, element in enumerate(elements, start=1):
        node_points = [points[node_id - 1] for node_id in element]
        center = (
            sum(p[0] for p in node_points[:4]) * 0.25,
            sum(p[1] for p in node_points[:4]) * 0.25,
            sum(p[2] for p in node_points[:4]) * 0.25,
        )
        normal = quad_basis(points, element)[2]
        axis = max(range(3), key=lambda idx: abs(normal[idx]))
        if abs(normal[axis]) < 1.0 - 1.0e-10:
            raise RuntimeError(f"single-domain element {ele_id} is not axis-aligned: normal={normal}")

        if axis == 0:
            if normal[0] < 0.0 and abs(center[0] + half_length) <= tol:
                counts["left"] += 1
            elif normal[0] > 0.0 and abs(center[0] - half_length) <= tol:
                counts["right"] += 1
            else:
                raise RuntimeError(
                    f"single-domain model contains an internal x-normal face: "
                    f"ele={ele_id} center={center} normal={normal}"
                )
        elif axis == 1:
            if normal[1] < 0.0 and abs(center[1] + half_width) <= tol:
                counts["y_min"] += 1
            elif normal[1] > 0.0 and abs(center[1] - half_width) <= tol:
                counts["y_max"] += 1
            else:
                raise RuntimeError(f"bad y-face location: ele={ele_id} center={center} normal={normal}")
        else:
            if normal[2] < 0.0 and abs(center[2] + half_height) <= tol:
                counts["z_min"] += 1
            elif normal[2] > 0.0 and abs(center[2] - half_height) <= tol:
                counts["z_max"] += 1
            else:
                raise RuntimeError(f"bad z-face location: ele={ele_id} center={center} normal={normal}")

    missing = [name for name, count in counts.items() if count <= 0]
    if missing:
        raise RuntimeError(f"single-domain model is missing exterior faces: {missing}")


def build_rod(domain_count: int,
              length: float,
              width: float,
              height: float,
              nx: int = 1,
              ny: int = 1,
              nz: int = 1) -> Tuple[List[Point], List[Element], List[Tuple[int, int, int, int, int]]]:
    """Build a rectangular rod from *domain_count* axially stacked domains.

    Returns
    -------
    points : List[Point]
        Point coordinates (0-indexed list; point IDs are 1-based).
    elements : List[Element]
        8-node quadrilateral elements (0-indexed list; element IDs are 1-based).
    interfaces : List[Tuple[int, int, int, int, int]]
        Each tuple is (domain_a, ele_a, domain_b, ele_b, normal_sign).
    """
    if domain_count < 1:
        raise ValueError("domain_count must be >= 1")
    if length <= 0.0 or width <= 0.0 or height <= 0.0:
        raise ValueError("length, width and height must be positive")
    if nx < 1 or ny < 1 or nz < 1:
        raise ValueError("nx, ny, nz must be >= 1")

    table = PointTable()
    elements: List[Element] = []
    faces_by_domain: List[Dict[str, List[int]]] = []
    dx = length / domain_count
    x_start = -0.5 * length

    for d in range(domain_count):
        x0 = x_start + d * dx
        x1 = x0 + dx
        faces_by_domain.append(
            build_domain(table, elements, x0, x1, width * 0.5, height * 0.5, nx, ny, nz)
        )

    # Interface pairing: right face of domain d <-> left face of domain d+1.
    # Both faces use the same (iy, iz) loop order, so element-wise alignment
    # is automatic.
    interfaces: List[Tuple[int, int, int, int, int]] = []
    for d in range(domain_count - 1):
        domain_a = d + 1
        domain_b = d + 2
        right_elems = faces_by_domain[d]["right"]
        left_elems = faces_by_domain[d + 1]["left"]
        if len(right_elems) != len(left_elems):
            raise RuntimeError(
                f"face element count mismatch at domain boundary {d}: "
                f"right={len(right_elems)} left={len(left_elems)}"
            )
        for ele_a, ele_b in zip(right_elems, left_elems):
            interfaces.append((domain_a, ele_a, domain_b, ele_b, -1))

    validate_rod_geometry(table.points, elements, faces_by_domain, interfaces)

    return table.points, elements, interfaces


def write_model(path: Path, points: Sequence[Point], elements: Sequence[Element]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as fp:
        fp.write(f"{len(points)}\n")
        for p in points:
            fp.write(f"{format_float(p[0])}\t{format_float(p[1])}\t{format_float(p[2])}\n")
        fp.write(f"{len(elements)}\n")
        for ele in elements:
            fp.write("\t".join(str(node_id) for node_id in ele) + "\n")
        fp.write("0\n")


def _right_end_element_set(domain_count: int, ele_per_domain: int, nx: int, ny: int, nz: int) -> Set[int]:
    """Return the set of element IDs belonging to the right-end face of the last domain."""
    face_count = ny * nz  # elements on the right face of one domain
    first_element_of_last_domain = (domain_count - 1) * ele_per_domain + 1
    # face order within domain: left(ny*nz), right(ny*nz), y_min, y_max, z_min, z_max
    right_face_start = first_element_of_last_domain + face_count
    right_face_end = right_face_start + face_count - 1
    return set(range(right_face_start, right_face_end + 1))


def _left_end_element_set(domain_count: int, ele_per_domain: int, ny: int, nz: int) -> Set[int]:
    """Return the set of element IDs belonging to the left-end face of the first domain."""
    face_count = ny * nz  # elements on the left face of one domain
    # face order within domain: left(ny*nz), right(ny*nz), y_min, y_max, z_min, z_max
    left_face_start = 1
    left_face_end = face_count
    return set(range(left_face_start, left_face_end + 1))


def boundary_value(step: int,
                   nstep: int,
                   element_id: int,
                   right_end_elements: Set[int],
                   points: Sequence[Point],
                   elements: Sequence[Element],
                   amplitude: float,
                   load_axis: str,
                   time_function: str) -> Tuple[float, float, float]:
    if amplitude == 0.0:
        scale = 0.0
    elif time_function == "constant":
        # Step 0 is the zero initial state used by the dynamic history terms.
        # Apply the constant load on physical time steps 1..nstep.
        scale = amplitude if step > 0 else 0.0
    elif time_function == "ramp":
        scale = amplitude * (step / nstep if nstep > 0 else 1.0)
    elif time_function == "sine":
        scale = amplitude * math.sin(math.pi * step / nstep if nstep > 0 else math.pi)
    else:
        raise ValueError(f"unsupported time function: {time_function}")

    # Only the right end face of the last domain receives the load.
    if element_id not in right_end_elements:
        return (0.0, 0.0, 0.0)

    if load_axis == "x":
        global_value = (scale, 0.0, 0.0)
    elif load_axis == "y":
        global_value = (0.0, scale, 0.0)
    elif load_axis == "z":
        global_value = (0.0, 0.0, scale)
    else:
        raise ValueError(f"unsupported load axis: {load_axis}")
    return global_to_local_components(points, elements, element_id, global_value)


def write_bd(path: Path,
             points: Sequence[Point],
             elements: Sequence[Element],
             domain_count: int,
             nstep: int,
             amplitude: float,
             load_axis: str,
             time_function: str,
             nx: int = 1,
             ny: int = 1,
             nz: int = 1) -> None:
    element_count = len(elements)
    ele_per_domain = element_count // domain_count

    left_end_elements = _left_end_element_set(domain_count, ele_per_domain, ny, nz)
    right_end_elements = _right_end_element_set(domain_count, ele_per_domain, nx, ny, nz)

    with path.open("w", encoding="ascii", newline="\n") as fp:
        # BCID per element: 123 = fixed (left end), 456 = traction (everything else)
        for ele in range(1, element_count + 1):
            if ele in left_end_elements:
                fp.write("123\n")
            else:
                fp.write("456\n")
        for step in range(nstep + 1):
            for ele in range(1, element_count + 1):
                # Left end: always zero (fixed displacement)
                if ele in left_end_elements:
                    value = (0.0, 0.0, 0.0)
                else:
                    value = boundary_value(
                        step, nstep, ele, right_end_elements,
                        points, elements, amplitude, load_axis, time_function)
                fp.write("\t".join(format_float(v) for v in value) + "\n")


def write_datacard(path: Path,
                   case_name: str,
                   domain_count: int,
                   interfaces: Sequence[Tuple[int, int, int, int, int]],
                   cfg: Config) -> None:
    with path.open("w", encoding="ascii", newline="\n") as fp:
        fp.write(f"{case_name}\n")
        fp.write(f"{cfg.threads}\n")
        fp.write(f"{cfg.flag_dyna}\n")
        fp.write(f"{format_float(cfg.beta)}\n")
        fp.write(f"{cfg.nstep}\n")
        fp.write(f"{format_float(cfg.gap)}\n")
        fp.write(f"{format_float(cfg.E)}\n")
        fp.write(f"{format_float(cfg.nu)}\n")
        fp.write(f"{format_float(cfg.rho)}\n")
        fp.write(f"{format_float(cfg.amplify)}\n")
        fp.write("2\n")
        fp.write(f"{cfg.iterations}\n")
        fp.write(f"{format_float(cfg.error)}\n")
        fp.write(f"{cfg.premaxleafpointnum}\n")
        fp.write("0\n")
        fp.write("0\n")
        fp.write(f"{domain_count}\n")
        for d in range(1, domain_count + 1):
            fp.write(f"{d} {format_float(cfg.E)} {format_float(cfg.nu)} {format_float(cfg.rho)}\n")
        fp.write(f"{len(interfaces)}\n")
        for domain_a, ele_a, domain_b, ele_b, normal_sign in interfaces:
            fp.write(f"{domain_a} {ele_a} {domain_b} {ele_b} {normal_sign}\n")


def write_legacy_datacard(path: Path, case_name: str, cfg: Config) -> None:
    """Write a datacard without the optional multi-domain tail."""
    with path.open("w", encoding="ascii", newline="\n") as fp:
        fp.write(f"{case_name}\n")
        fp.write(f"{cfg.threads}\n")
        fp.write(f"{cfg.flag_dyna}\n")
        fp.write(f"{format_float(cfg.beta)}\n")
        fp.write(f"{cfg.nstep}\n")
        fp.write(f"{format_float(cfg.gap)}\n")
        fp.write(f"{format_float(cfg.E)}\n")
        fp.write(f"{format_float(cfg.nu)}\n")
        fp.write(f"{format_float(cfg.rho)}\n")
        fp.write(f"{format_float(cfg.amplify)}\n")
        fp.write("2\n")
        fp.write(f"{cfg.iterations}\n")
        fp.write(f"{format_float(cfg.error)}\n")
        fp.write(f"{cfg.premaxleafpointnum}\n")
        fp.write("0\n")
        fp.write("0\n")


# ============================================================
# Editable default parameters. Run this script to regenerate the current rod120 case.
# ============================================================
class Config:
    domains: int = 10              # axial domain count
    case: str = "10domain_odd_5"           # case name
    out_dir: str = "DBEM1/input"   # output directory
    card_out: str = ""             # datacard output path; empty = input/<case>.DATACARD
    write_single_domain_case: bool = False  # write a single-domain outer-surface case
    single_case_suffix: str = ""
    single_card_out: str = ""      # empty = input/<case>_single.DATACARD

    length: float = 5           # rod length
    width: float = 0.5             # cross-section width
    height: float = 0.5            # cross-section height

    # Mesh subdivisions per domain face.
    nx: int = 5                    # x-axis subdivisions on side faces
    ny: int = 5                    # y subdivisions
    nz: int = 5                    # z subdivisions

    nstep: int = 1000                # time steps
    load_axis: str = "x"           # global load direction: "x", "y", "z"
    load_amplitude: float = 0.2   # right-end traction amplitude
    time_function: str = "constant"  # "constant", "ramp", "sine"

    threads: int = 20
    flag_dyna: int = 2
    beta: float = 0.6
    gap: float = 0.4
    E: float = 200
    nu: float = 0
    rho: float = 8
    amplify: float = 1.0
    iterations: int = 200
    error: float = 1.0e-5
    premaxleafpointnum: int = 10

def main() -> None:
    cfg = Config()

    nx, ny, nz = cfg.nx, cfg.ny, cfg.nz
    if nx < 1 or ny < 1 or nz < 1:
        raise SystemExit("nx, ny, nz must be >= 1")

    e_per_domain = elems_per_domain(nx, ny, nz)
    total_elements = cfg.domains * e_per_domain
    case_name = cfg.case or f"rod{total_elements}"
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    points, elements, interfaces = build_rod(
        cfg.domains, cfg.length, cfg.width, cfg.height, nx, ny, nz
    )
    if len(elements) % cfg.domains != 0:
        raise RuntimeError("element count must be divisible by domain count")
    if len(elements) != cfg.domains * e_per_domain:
        raise RuntimeError(
            f"unexpected element count: {len(elements)} vs {cfg.domains * e_per_domain}"
        )

    model_path = out_dir / f"{case_name}.model"
    bd_path = out_dir / f"{case_name}.bd"
    card_path = Path(cfg.card_out) if cfg.card_out else out_dir / f"{case_name}.DATACARD"

    write_model(model_path, points, elements)
    write_bd(bd_path, points, elements, cfg.domains, cfg.nstep,
             cfg.load_amplitude, cfg.load_axis, cfg.time_function,
             nx, ny, nz)
    write_datacard(card_path, case_name, cfg.domains, interfaces, cfg)

    single_case_name = f"{case_name}{cfg.single_case_suffix}"
    single_model_path = single_bd_path = single_card_path = None
    single_points: List[Point] = []
    single_elements: List[Element] = []
    single_constant_node_count = 0
    if cfg.write_single_domain_case:
        single_nx = cfg.domains * nx
        single_points, single_elements, single_interfaces = build_rod(
            1, cfg.length, cfg.width, cfg.height, single_nx, ny, nz
        )
        if single_interfaces:
            raise RuntimeError("single-domain outer-surface case must not contain interfaces")
        validate_single_domain_outer_surface(single_points, single_elements, cfg.length, cfg.width, cfg.height)
        single_constant_node_count = len(single_elements)
        single_model_path = out_dir / f"{single_case_name}.model"
        single_bd_path = out_dir / f"{single_case_name}.bd"
        single_card_path = (
            Path(cfg.single_card_out)
            if cfg.single_card_out
            else out_dir / f"{single_case_name}.DATACARD"
        )
        write_model(single_model_path, single_points, single_elements)
        write_bd(
            single_bd_path, single_points, single_elements, 1, cfg.nstep,
            cfg.load_amplitude, cfg.load_axis, cfg.time_function,
            single_nx, ny, nz
        )
        write_legacy_datacard(single_card_path, single_case_name, cfg)

    face_nyz = ny * nz
    face_nxz = nx * nz
    face_nxy = nx * ny

    print(f"case: {case_name}")
    print(f"mesh divisions: nx={nx} ny={ny} nz={nz}")
    print(f"points: {len(points)}")
    print(f"elements: {len(elements)} ({e_per_domain} per domain)")
    print(f"  left/right faces: {face_nyz} elements each")
    print(f"  y-min/y-max faces: {face_nxz} elements each")
    print(f"  z-min/z-max faces: {face_nxy} elements each")
    print(f"interfaces: {len(interfaces)} ({face_nyz} per domain boundary)")
    print(f"model: {model_path}")
    print(f"bd: {bd_path}")
    print(f"datacard: {card_path}")
    if cfg.write_single_domain_case:
        print("single-domain outer-surface case:")
        print(f"  case: {single_case_name}")
        print(f"  mesh divisions: nx={cfg.domains * nx} ny={ny} nz={nz}")
        print(f"  geometry points: {len(single_points)}")
        print(f"  elements: {len(single_elements)} (no internal interface faces)")
        print(f"  constant-element physical nodes: {single_constant_node_count}")
        print(f"  model: {single_model_path}")
        print(f"  bd: {single_bd_path}")
        print(f"  datacard: {single_card_path}")

    # Show first few interfaces as a sanity check
    max_show = min(5, len(interfaces))
    for itf in interfaces[:max_show]:
        print(f"interface: domain {itf[0]} ele {itf[1]} <-> domain {itf[2]} ele {itf[3]}, normalSign {itf[4]}")
    if len(interfaces) > max_show:
        print(f"  ... and {len(interfaces) - max_show} more interface pairs")


if __name__ == "__main__":
    main()


