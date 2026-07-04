#!/usr/bin/env python
"""Convert DBEM text .model/.bd inputs to binary files.

Examples:
  python convert_inputs_to_binary.py --model DBEM1/input/earth1000000.model
  python convert_inputs_to_binary.py --bd DBEM1/input/earth1000000.bd --ele-num 1000 --nstep 5
  python convert_inputs_to_binary.py --bd-dir DBEM1/input/bd_400 --pattern "*.bd" --ele-num 16384 --nstep 63
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path


MAGIC = b"DBEMBIN\0"
VERSION = 1
KIND_MODEL = 1
KIND_BD = 2
ENDIAN_MARK = 0x01020304


def tokens(path: Path):
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            for item in line.split():
                yield item


def write_header(f, kind: int) -> None:
    f.write(struct.pack("<8sIIII", MAGIC, VERSION, kind, ENDIAN_MARK, 0))


def write_chunked(f, values, fmt: str, chunk_size: int = 65536) -> None:
    chunk = []
    for value in values:
        chunk.append(value)
        if len(chunk) >= chunk_size:
            f.write(struct.pack("<" + fmt * len(chunk), *chunk))
            chunk.clear()
    if chunk:
        f.write(struct.pack("<" + fmt * len(chunk), *chunk))


def convert_model(path: Path) -> Path:
    it = tokens(path)
    out = path.with_suffix(path.suffix + ".bin")
    with out.open("wb") as f:
        write_header(f, KIND_MODEL)

        point_num = int(next(it))
        f.write(struct.pack("<q", point_num))
        write_chunked(f, (float(next(it)) for _ in range(point_num * 3)), "d")

        ele_num = int(next(it))
        f.write(struct.pack("<q", ele_num))
        write_chunked(f, (int(next(it)) - 1 for _ in range(ele_num * 8)), "q")

        try:
            inf_ele_num = int(next(it))
        except StopIteration:
            inf_ele_num = 0
        f.write(struct.pack("<q", inf_ele_num))
        write_chunked(f, (int(next(it)) - 1 for _ in range(inf_ele_num * 4)), "q")

    return out


def convert_bd(path: Path, ele_num: int, nstep: int) -> Path:
    it = tokens(path)
    out = path.with_suffix(path.suffix + ".bin")
    with out.open("wb") as f:
        write_header(f, KIND_BD)
        f.write(struct.pack("<qq", ele_num, nstep))

        write_chunked(f, (int(next(it)) for _ in range(ele_num)), "i")

        value_count = (nstep + 1) * ele_num * 3
        write_chunked(f, (float(next(it)) for _ in range(value_count)), "d")

    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert DBEM .model/.bd text files to binary inputs.")
    parser.add_argument("--model", type=Path, action="append", default=[], help="Text .model file to convert.")
    parser.add_argument("--bd", type=Path, action="append", default=[], help="Text .bd file to convert.")
    parser.add_argument("--bd-dir", type=Path, default=None, help="Directory containing .bd files to convert.")
    parser.add_argument("--pattern", default="*.bd", help="Pattern for --bd-dir. Default: *.bd")
    parser.add_argument("--ele-num", type=int, default=None, help="Element count required for .bd conversion.")
    parser.add_argument("--nstep", type=int, default=None, help="NStep required for .bd conversion.")
    args = parser.parse_args()

    outputs: list[Path] = []

    for model in args.model:
        outputs.append(convert_model(model))

    bd_files = list(args.bd)
    if args.bd_dir:
        bd_files.extend(sorted(args.bd_dir.glob(args.pattern)))

    if bd_files and (args.ele_num is None or args.nstep is None):
        parser.error("--ele-num and --nstep are required when converting .bd files")

    for bd in bd_files:
        outputs.append(convert_bd(bd, args.ele_num, args.nstep))

    for out in outputs:
        print(out)
    print(f"Converted {len(outputs)} file(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
