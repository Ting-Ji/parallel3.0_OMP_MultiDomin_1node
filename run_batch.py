"""Run DBEM1 cases listed in a text file, one process at a time."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
SUCCESS_CODES = {0, 1}  # DBEM1 currently returns 1 after a successful run.


@dataclass
class Case:
    number: int
    requested: str
    model: Path | None = None
    relative_name: str = ""
    card_bytes: bytes = b""
    error: str = ""


def read_list(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig") as stream:
        return [line.strip() for line in stream if line.strip() and not line.lstrip().startswith("#")]


def resolve_cases(input_dir: Path, names: list[str]) -> list[Case]:
    # Index only models. A missing companion is reported for that particular case.
    models: dict[str, list[Path]] = {}
    for model in input_dir.rglob("*.model"):
        models.setdefault(model.stem, []).append(model)
    resolved: list[Case] = []
    for number, requested in enumerate(names, 1):
        case = Case(number, requested)
        resolved.append(case)
        normalized = requested.replace("\\", "/")
        parts = normalized.split("/")
        if (not normalized or any(part in ("", ".", "..") for part in parts)
                or ":" in normalized or requested.endswith((".model", ".bd", ".DATACARD"))):
            case.error = "模型名必须是 input 下的相对路径，且不带扩展名"
            continue
        try:
            encoded_name = normalized.encode("ascii")
        except UnicodeEncodeError:
            case.error = "模型名须为 ASCII；当前 DBEM1 的参数卡读取方式不支持非 ASCII 模型名"
            continue
        if len(encoded_name) > 48 or any(chr(byte).isspace() for byte in encoded_name):
            case.error = "模型名不能含空格且不能超过 48 字节（DBEM1 读取限制）"
            continue
        if len(parts) > 1:
            candidates = [Path(str(input_dir.joinpath(*parts)) + ".model")]
        else:
            candidates = models.get(requested, [])
        if len(candidates) > 1:
            paths = ", ".join(sorted(path.relative_to(input_dir).as_posix() for path in candidates))
            case.error = f"存在同名模型，请在清单中写相对路径：{paths}"
            continue
        if not candidates or not candidates[0].is_file():
            case.error = "找不到 .model"
            continue
        model = candidates[0]
        if not model.resolve().is_relative_to(input_dir.resolve()):
            case.error = "模型路径超出 input 目录"
            continue
        missing = [suffix for suffix in (".bd", ".DATACARD")
                   if not model.with_suffix(suffix).is_file()]
        if missing:
            case.error = "缺少 " + "、".join(missing)
            continue
        card = model.with_suffix(".DATACARD")
        try:
            raw = card.read_bytes()
            token = re.search(rb"\S+", raw)
            if token is None or raw.startswith(b"\xef\xbb\xbf"):
                raise ValueError("参数卡为空或包含 UTF-8 BOM")
            case.card_bytes = raw[:token.start()] + encoded_name + raw[token.end():]
        except (OSError, ValueError) as exc:
            case.error = f"无法读取参数卡：{exc}"
            continue
        case.model = model
        case.relative_name = model.relative_to(input_dir).with_suffix("").as_posix()
    return resolved


def atomic_card_write(path: Path, content: bytes) -> None:
    temporary = path.with_name(path.name + ".batch.tmp")
    try:
        temporary.write_bytes(content)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def stop_child(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()


def run_batch(project_root: Path, list_path: Path, exe: Path, dry_run: bool = False) -> int:
    dbem_dir = project_root / "DBEM1"
    input_dir = dbem_dir / "input"
    card_path = dbem_dir / "BEM_DATACARD.DAT"
    public_log = project_root / "output" / "log.txt"
    try:
        names = read_list(list_path)
    except OSError as exc:
        print(f"无法读取清单：{exc}", file=sys.stderr)
        return 2
    if not names:
        print("清单中没有模型名", file=sys.stderr)
        return 2
    cases = resolve_cases(input_dir, names)
    missing_globals = [str(path) for path in (exe, dbem_dir / "gaussposition",
                       dbem_dir / "gaussweight", card_path) if not path.is_file()]
    if missing_globals:
        print("缺少运行所需文件：" + "，".join(missing_globals), file=sys.stderr)
        return 2
    for case in cases:
        label = case.relative_name or case.requested
        print(f"[{case.number}/{len(cases)}] {label}: " + (case.error or "文件齐全"), flush=True)
    if dry_run:
        failures = sum(bool(case.error) for case in cases)
        print(f"检查完毕：可运行 {len(cases) - failures}，有问题 {failures}；未修改文件")
        return 1 if failures else 0

    output_root = project_root / "output"
    output_root.mkdir(exist_ok=True)
    lock_path = output_root / ".dbem_batch.lock"
    try:
        lock_fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        print(f"已有批处理锁：{lock_path}。请确认旧进程已停止，并按锁内备份路径恢复参数卡。", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"无法创建批处理锁：{exc}", file=sys.stderr)
        return 2

    batch_dir = output_root / ("batch_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
    backup = batch_dir / "BEM_DATACARD.original.DAT"
    results: list[dict[str, object]] = []
    interrupted = False
    backup_complete = False
    original_stat: os.stat_result | None = None
    try:
        with os.fdopen(lock_fd, "w", encoding="utf-8") as stream:
            json.dump({"pid": os.getpid(), "backup": str(backup), "list": str(list_path)}, stream,
                      ensure_ascii=False)
            stream.flush()
            os.fsync(stream.fileno())
        batch_dir.mkdir()
        original_stat = card_path.stat()
        backup.write_bytes(card_path.read_bytes())
        backup_complete = True
        for case in cases:
            start = time.monotonic()
            code: int | None = None
            message = case.error
            status = "FAILED" if message else "SUCCESS"
            if not message:
                stdout_log = batch_dir / f"{case.number:04d}_console.log"
                log_before = public_log.stat().st_mtime_ns if public_log.is_file() else None
                started = False
                process: subprocess.Popen[bytes] | None = None
                try:
                    atomic_card_write(card_path, case.card_bytes)
                    with stdout_log.open("wb") as stream:
                        process = subprocess.Popen([str(exe)], cwd=dbem_dir,
                                                   stdout=stream, stderr=subprocess.STDOUT)
                        started = True
                        try:
                            code = process.wait()
                        except KeyboardInterrupt:
                            stop_child(process)
                            status = "INTERRUPTED"
                            message = "用户中断"
                            interrupted = True
                    if not interrupted and code not in SUCCESS_CODES:
                        status = "FAILED"
                        message = f"进程返回码 {code}"
                except KeyboardInterrupt:
                    if process is not None:
                        stop_child(process)
                    status = "INTERRUPTED"
                    message = "用户中断"
                    interrupted = True
                except (OSError, subprocess.SubprocessError) as exc:
                    status = "FAILED"
                    message = f"启动或运行失败：{exc}"
                if started and public_log.is_file() and public_log.stat().st_mtime_ns != log_before:
                    try:
                        shutil.copy2(public_log, batch_dir / f"{case.number:04d}_dbem.log")
                    except OSError as exc:
                        message = (message + "；" if message else "") + f"复制公共日志失败：{exc}"
            results.append({"序号": case.number, "清单模型名": case.requested,
                            "匹配模型名": case.relative_name, "状态": status,
                            "返回码": "" if code is None else code,
                            "耗时秒": f"{time.monotonic() - start:.2f}", "说明": message})
            print(f"[{case.number}/{len(cases)}] {case.requested}: {status}"
                  + (f"（{message}）" if message else ""), flush=True)
            if interrupted:
                break
        with (batch_dir / "summary.csv").open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["序号", "清单模型名", "匹配模型名",
                                                       "状态", "返回码", "耗时秒", "说明"])
            writer.writeheader()
            writer.writerows(results)
        failures = sum(row["状态"] != "SUCCESS" for row in results)
        print(f"汇总：成功 {len(results) - failures}，失败/中断 {failures}；日志：{batch_dir}")
        return 130 if interrupted else (1 if failures else 0)
    finally:
        try:
            if backup_complete:
                atomic_card_write(card_path, backup.read_bytes())
                if original_stat is not None:
                    os.chmod(card_path, stat.S_IMODE(original_stat.st_mode))
                    os.utime(card_path, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
        except OSError as exc:
            print(f"参数卡恢复失败！请手动从 {backup} 恢复：{exc}", file=sys.stderr)
            raise
        else:
            lock_path.unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="按清单顺序调用 DBEM1.exe")
    parser.add_argument("cases", type=Path, help="模型清单；相对路径以脚本所在目录为基准")
    parser.add_argument("--exe", type=Path, default=Path("DBEM1/x64/Release/DBEM1.exe"),
                        help="DBEM1 可执行文件；相对路径以脚本所在目录为基准")
    parser.add_argument("--dry-run", action="store_true", help="只检查文件，不修改参数卡或运行 EXE")
    args = parser.parse_args(argv)
    list_path = (PROJECT_ROOT / args.cases).resolve()
    exe = (PROJECT_ROOT / args.exe).resolve()
    return run_batch(PROJECT_ROOT, list_path, exe, args.dry_run)


if __name__ == "__main__":
    sys.exit(main())
