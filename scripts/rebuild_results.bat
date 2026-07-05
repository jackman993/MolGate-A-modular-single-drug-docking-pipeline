#!/usr/bin/env python3
"""Re-export catalog-300 release into ./results (wrapper; source of truth in moltgate_UI)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UI_SCRIPT = ROOT.parent / "moltgate_UI" / "scripts" / "export_catalog_300_release.py"
OUT = ROOT / "results"

if not UI_SCRIPT.is_file():
    print(f"Missing {UI_SCRIPT}")
    print("Clone moltgate_UI beside MolGate-CMD-Release, or run export_catalog_300_release.py manually.")
    sys.exit(1)

cmd = [sys.executable, str(UI_SCRIPT), "--out-dir", str(OUT)]
raise SystemExit(subprocess.call(cmd))
