"""Build RHS parallel solver/tests with the same MSVC flags as the saved baseline."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if __name__ == "__main__":
    extra = sys.argv[1:]
    if "--output" not in extra:
        extra = ["--output", str(ROOT / "tmp/rhs_parallel/build")] + extra
    subprocess.run([sys.executable, str(ROOT / "tools/build_gmres_openmp.py"), "--rhs-tests", *extra], check=True)
