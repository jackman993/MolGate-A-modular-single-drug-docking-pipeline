#!/usr/bin/env python3
"""Sync master_index + drug_catalog_300 from moltgate_UI into CMD engines (Phase C)."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UI = ROOT.parent / "moltgate_UI"
UI_CONFIG = UI / "config"
ENGINES = ROOT / "engines"
OUT_CATALOG = ENGINES / "config" / "drug_catalog_300.csv"
OUT_INDEX = ENGINES / "master_index.json"
RESULTS_CATALOG = ROOT / "results" / "catalog_300_registry.csv"


def merge_registry_csv() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for p in sorted(UI_CONFIG.glob("drug_catalog_batch*.csv")):
        rows.extend(csv.DictReader(p.open(encoding="utf-8-sig")))
    by_id: dict[int, dict[str, str]] = {}
    for r in rows:
        try:
            by_id[int(r["drug_id"])] = r
        except (KeyError, ValueError):
            continue
    return [by_id[k] for k in sorted(by_id)]


def write_catalog(rows: list[dict[str, str]], path: Path) -> None:
    if not rows:
        raise SystemExit("No catalog rows to write")
    fields = list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def count_index_keys(path: Path) -> int:
    raw = json.loads(path.read_text(encoding="utf-8"))
    return sum(1 for k in raw if not str(k).startswith("_") and isinstance(raw[k], dict))


def main() -> int:
    ap = argparse.ArgumentParser(description="Sync catalog-300 + master_index from moltgate_UI")
    ap.add_argument("--ui-home", default=str(UI), help="Path to moltgate_UI")
    args = ap.parse_args()
    ui = Path(args.ui_home).resolve()
    ui_config = ui / "config"
    ui_index = ui_config / "master_index.json"

    if not ui_index.is_file():
        print(f"Missing {ui_index}")
        return 1

    rows = merge_registry_csv() if ui_config.is_dir() else []
    if not rows and RESULTS_CATALOG.is_file():
        rows = list(csv.DictReader(RESULTS_CATALOG.open(encoding="utf-8-sig")))
        print(f"[sync] catalog from existing {RESULTS_CATALOG}")
    elif rows:
        print(f"[sync] catalog from {ui_config}/drug_catalog_batch*.csv ({len(rows)} rows)")
    else:
        print("No catalog source found")
        return 1

    shutil.copy2(ui_index, OUT_INDEX)
    write_catalog(rows, OUT_CATALOG)
    write_catalog(rows, RESULTS_CATALOG)

    print(f"[sync] master_index -> {OUT_INDEX} ({count_index_keys(OUT_INDEX)} keys)")
    print(f"[sync] catalog -> {OUT_CATALOG} ({len(rows)} rows)")
    print(f"[sync] catalog -> {RESULTS_CATALOG} (mirror)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
