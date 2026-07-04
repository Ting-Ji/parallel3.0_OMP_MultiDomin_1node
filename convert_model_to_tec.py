#!/usr/bin/env python
"""将 DBEM .model 文件转换为 Tecplot FEPOINT .dat 格式
   8节点四边形 → 拆分为4个4节点四边形 + 单元中心点

用法: python convert_model_to_tec.py DBEM1/input/rod5.model
"""

import sys
from pathlib import Path


def convert_model_to_tecplot(model_path: str, output_path: str = None):
    model_path = Path(model_path)
    if output_path is None:
        output_path = model_path.with_suffix(".dat")
    else:
        output_path = Path(output_path)

    with open(model_path, "r") as f:
        lines = f.readlines()

    lines = [l.strip() for l in lines if l.strip()]

    idx = 0
    num_nodes = int(lines[idx]); idx += 1

    # 读取节点坐标
    nodes = []
    for _ in range(num_nodes):
        parts = lines[idx].split()
        nodes.append(tuple(float(p) for p in parts[:3]))
        idx += 1

    num_elems = int(lines[idx]); idx += 1

    # 读取单元连接 (1-based → 0-based)
    elems = []
    for _ in range(num_elems):
        parts = lines[idx].split()
        elems.append([int(p) - 1 for p in parts[:8]])
        idx += 1

    # 计算每个单元的中心点坐标
    # 参照 write_tec.cpp 中 cal_center 的算法
    centers = []
    for e in elems:
        cx = cy = cz = 0.0
        for k in range(4):       # 角节点权重 -0.25
            n = e[k]
            cx += -0.25 * nodes[n][0]
            cy += -0.25 * nodes[n][1]
            cz += -0.25 * nodes[n][2]
        for k in range(4, 8):    # 中节点权重 0.5
            n = e[k]
            cx += 0.5 * nodes[n][0]
            cy += 0.5 * nodes[n][1]
            cz += 0.5 * nodes[n][2]
        centers.append((cx, cy, cz))

    # 拆分 8 节点四边形为 4 个标准四边形
    # 单元局部节点: 角点1-4, 中节点5-8
    # 4个子四边形:
    #   Q1: 角1, 中5, 中心, 中8
    #   Q2: 角2, 中6, 中心, 中5
    #   Q3: 角3, 中7, 中心, 中6
    #   Q4: 角4, 中8, 中心, 中7
    sub_quads = [
        (0, 4, 8, 7),   # 角1, 中5, 中心, 中8
        (1, 5, 8, 4),   # 角2, 中6, 中心, 中5
        (2, 6, 8, 5),   # 角3, 中7, 中心, 中6
        (3, 7, 8, 6),   # 角4, 中8, 中心, 中7
    ]

    total_nodes = num_nodes + num_elems  # 原始节点 + 中心点
    total_elems = num_elems * 4          # 每个单元拆成4个

    with open(output_path, "w") as f:
        # === Tecplot 文件头 ===
        f.write('TITLE = "BEM Model"\n')
        f.write('VARIABLES = "X", "Y", "Z"\n')
        f.write(f'ZONE T="mesh", N={total_nodes}, E={total_elems}, F=FEPOINT, ET=QUADRILATERAL\n')

        # 输出所有原始节点坐标
        for x, y, z in nodes:
            f.write(f"{x}\t{y}\t{z}\n")

        # 输出单元中心点坐标
        for cx, cy, cz in centers:
            f.write(f"{cx}\t{cy}\t{cz}\n")

        # 输出拆分后的单元连接 (1-based)
        center_start = num_nodes + 1  # 中心点起始编号
        for ei, e in enumerate(elems):
            center_id = center_start + ei
            # 构建局部→全局节点映射
            local = [
                e[0] + 1, e[1] + 1, e[2] + 1, e[3] + 1,   # 角节点 (1-based)
                e[4] + 1, e[5] + 1, e[6] + 1, e[7] + 1,   # 中节点
                center_id                                     # 中心点
            ]
            for sq in sub_quads:
                f.write(f"{local[sq[0]]}\t{local[sq[1]]}\t{local[sq[2]]}\t{local[sq[3]]}\n")

    print(f"转换完成: {output_path}")
    print(f"  原始节点数: {num_nodes}")
    print(f"  中心点数量: {num_elems}")
    print(f"  总节点数:   {total_nodes}")
    print(f"  总单元数:   {num_elems} × 4 = {total_elems}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python convert_model_to_tec.py <model文件路径> [输出路径]")
        print("示例: python convert_model_to_tec.py DBEM1/input/rod5.model")
        sys.exit(1)
    convert_model_to_tecplot(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
