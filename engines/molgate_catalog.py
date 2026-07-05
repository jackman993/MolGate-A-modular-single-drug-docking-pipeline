"""Load drug catalog-300 CSV for --catalog-drug-id runner (Phase C)."""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

_PKG_ROOT = Path(__file__).resolve().parent
_CATALOG_CANDIDATES = (
    _PKG_ROOT / "config" / "drug_catalog_300.csv",
    _PKG_ROOT.parent / "results" / "catalog_300_registry.csv",
)


def catalog_path() -> Path:
    for p in _CATALOG_CANDIDATES:
        if p.is_file():
            return p
    raise FileNotFoundError(
        "找不到 catalog CSV。請執行 scripts/sync_catalog_from_ui.bat "
        "或確認 engines/config/drug_catalog_300.csv 存在。"
    )


@lru_cache(maxsize=1)
def load_catalog() -> dict[int, dict[str, str]]:
    path = catalog_path()
    rows: dict[int, dict[str, str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            try:
                did = int(row["drug_id"])
            except (KeyError, TypeError, ValueError):
                continue
            rows[did] = row
    return rows


def catalog_row_to_drug(row: dict[str, str]) -> dict:
    """Convert catalog CSV row to runner drug dict."""
    notes = (row.get("notes") or "").strip()
    anchor = (row.get("anchor_mode") or "").strip()
    tier = (row.get("tier") or "").strip()
    extra = f" | tier={tier}" if tier else ""
    if anchor:
        extra += f" | anchor={anchor}"
    return {
        "drug": (row.get("canonical_name") or "").strip(),
        "index_key": (row.get("index_key") or "").strip(),
        "smiles": (row.get("smiles") or "").strip(),
        "target": (row.get("target_label") or row.get("index_key") or "").strip(),
        "pdb_id": (row.get("pdb_id") or "").strip().upper(),
        "het_id": (row.get("het_id") or "").strip().upper(),
        "catalog_drug_id": int(row["drug_id"]),
        "note": (notes + extra).strip(" |") or f"catalog drug_id={row['drug_id']}",
        "_from_catalog": True,
    }


def drug_by_catalog_id(catalog_drug_id: int) -> dict:
    rows = load_catalog()
    try:
        row = rows[int(catalog_drug_id)]
    except KeyError as exc:
        raise KeyError(catalog_drug_id) from exc
    return catalog_row_to_drug(row)
