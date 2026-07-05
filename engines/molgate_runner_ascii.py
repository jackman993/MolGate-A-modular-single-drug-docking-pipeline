"""
ASCII-safe wrapper for MolGate runner.
Runs the pretty runner and strips non-ASCII chars from output.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
PYTHON = sys.executable


def _to_ascii_line(line: str) -> str:
    return line.encode("ascii", errors="replace").decode("ascii", errors="replace")


def main(argv: list[str] | None = None) -> int:
    args = argv if argv is not None else sys.argv[1:]
    cmd = [PYTHON, str(BASE_DIR / "molgate_runner_pretty.py"), *args]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(_to_ascii_line(line))
    return proc.wait()


if __name__ == "__main__":
    raise SystemExit(main())
