"""Wire $f_route into CMD single-drug runner."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

from molgate_anchor_routing import (
    RouteResolution,
    apply_resolution_to_session_index,
    load_catalog_rows,
    resolve_anchor_route,
    set_allow_network_fetch,
    write_route_resolution,
)
from molgate_pdb_cache import PDB_DIR, ensure_pdb_cached, ensure_runtime_dirs


def build_route_entry(drug: dict[str, Any], load_master_index: Callable[[], dict[str, dict]]) -> dict[str, Any]:
    key = str(drug.get("index_key") or drug.get("target") or "").strip().upper()
    idx = load_master_index().get(key, {})
    entry: dict[str, Any] = {**idx}
    entry["index_key"] = key
    entry["pdb_id"] = (drug.get("pdb_id") or idx.get("pdb_id") or "").strip().upper()
    entry["het_id"] = (drug.get("het_id") or idx.get("het_id") or idx.get("cocrystal_het") or "").strip().upper()
    entry["target_label"] = drug.get("target") or idx.get("target_label") or idx.get("target") or ""
    entry["target_drug"] = drug.get("drug") or idx.get("target_drug") or idx.get("common_name") or ""
    if drug.get("_from_catalog") and drug.get("catalog_drug_id") is not None:
        from molgate_catalog import load_catalog

        row = load_catalog().get(int(drug["catalog_drug_id"]), {})
        for field in ("anchor_mode", "pocket_rule", "tier", "target_label"):
            if row.get(field):
                entry[field] = row[field]
    return entry


def apply_anchor_routing(
    drug: dict[str, Any],
    load_master_index: Callable[[], dict[str, dict]],
    *,
    allow_network: bool,
) -> tuple[dict[str, Any], RouteResolution]:
    ensure_runtime_dirs()
    set_allow_network_fetch(allow_network)
    entry = build_route_entry(drug, load_master_index)
    pdb_id = str(entry.get("pdb_id") or "").strip().upper()
    if pdb_id:
        ensure_pdb_cached(pdb_id, allow_network=allow_network)

    catalog_rows = load_catalog_rows()
    resolution = resolve_anchor_route(entry, PDB_DIR, catalog_rows=catalog_rows)

    if resolution.fallback_applied:
        eff = resolution.effective_entry
        drug = {
            **drug,
            "pdb_id": str(eff.get("pdb_id") or drug.get("pdb_id") or "").strip().upper(),
            "het_id": str(
                eff.get("het_id") or eff.get("cocrystal_het") or drug.get("het_id") or ""
            ).strip().upper(),
        }
        eff_pdb = str(drug.get("pdb_id") or "").strip().upper()
        if eff_pdb:
            ensure_pdb_cached(eff_pdb, allow_network=allow_network)
        print(f"  [route] {resolution.status}: {resolution.reason}")
    elif resolution.status in ("direct_ok", "proxy_ok"):
        print(f"  [route] {resolution.status}: {resolution.reason}")
    elif resolution.status == "no_route":
        print(f"  [route] WARNING {resolution.status}: {resolution.reason}")

    return drug, resolution


def attach_route_to_session(session_dir: Path, resolution: RouteResolution) -> None:
    write_route_resolution(session_dir, resolution)
    apply_resolution_to_session_index(session_dir, resolution)
