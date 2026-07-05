"""
Anchor route resolution ($f_route) — verify catalog HET in PDB; fallback Direct → Proxy.

Ported from molgate_ui/anchor_routing.py for CMD release (DrugOps paper).
"""
from __future__ import annotations

import csv
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from molgate_pdb_cache import PDB_DIR, ensure_pdb_cached, find_cached_pdb

_CONFIG_DIR = Path(__file__).resolve().parent / "config"
_CATALOG_PATH = _CONFIG_DIR / "drug_catalog_300.csv"

_BUFFER_IONS = {
    "HOH", "WAT", "GOL", "EDO", "PEG", "MPD", "DMS", "FMT", "ACT", "TRS", "MES",
    "EPE", "BME", "DTT", "TCE", "IMD", "SO4", "PO4", "CIT", "ACY", "EOH", "IPA",
    "MOH", "BOG", "CL", "NA", "MG", "ZN", "CA", "K", "FE", "MN", "CU", "NI", "CO",
    "CD", "HG", "PT", "AU", "AG", "SR", "BA", "CS", "RB", "IOD", "BR", "F", "LI",
    "AL", "SE",
}
_GLYCANS = {
    "NAG", "NGA", "BMA", "MAN", "GAL", "GLC", "GCS", "MAL", "MAB", "MBG", "BGC",
    "FUC", "FUL", "RAM", "LFR", "GLO", "XYP", "XYS", "XYF", "XUL", "SIA", "SLB",
    "A2G", "M6D", "MMA", "MGL", "NGK", "IDR",
}
_MW_CUTOFF = 150.0

SESSION_ROUTE_KEYS = (
    "pdb_id",
    "het_id",
    "cocrystal_het",
    "pocket_rule",
    "anchor_mode",
    "state",
    "target_drug",
    "common_name",
    "tier",
    "_route_fallback",
    "_original_pdb_id",
    "_original_het_id",
    "_proxy_source_key",
    "_route_resolution_status",
)

_SURROGATE_RE = re.compile(r"^(?:surrogate_[^:]+|cocrystal_ligand):([A-Z0-9]{3})$", re.I)

_allow_network_fetch = False


def set_allow_network_fetch(flag: bool) -> None:
    global _allow_network_fetch
    _allow_network_fetch = bool(flag)


@dataclass
class HetResidue:
    resname: str
    chain: str
    resseq: int
    heavy_atom_count: int
    mw_estimate: float

    def as_candidate(self, *, role: str = "candidate") -> dict[str, Any]:
        return {
            "resname": self.resname,
            "chain": self.chain,
            "resseq": self.resseq,
            "heavy_atom_count": self.heavy_atom_count,
            "mw_estimate": round(self.mw_estimate, 1),
            "ligand_role": role,
        }


@dataclass
class RouteResolution:
    status: str
    het_found: bool
    fallback_applied: bool
    reason: str
    original_entry: dict[str, Any]
    effective_entry: dict[str, Any]
    pdb_ligand_resnames: list[str] = field(default_factory=list)
    drug_like_resnames: list[str] = field(default_factory=list)
    proxy_source_key: str | None = None
    suggested_ccd: str | None = None

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


def load_catalog_rows(config_dir: Path | None = None) -> list[dict[str, str]]:
    root = config_dir or _CONFIG_DIR
    catalog = root / "drug_catalog_300.csv"
    if catalog.is_file():
        with catalog.open(encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.DictReader(fh))
        rows.sort(key=lambda r: int(r["drug_id"]))
        return rows
    rows: list[dict[str, str]] = []
    for path in sorted(root.glob("drug_catalog_batch*.csv")):
        with path.open(encoding="utf-8-sig", newline="") as fh:
            rows.extend(csv.DictReader(fh))
    rows.sort(key=lambda r: int(r["drug_id"]))
    return rows


def resolve_pdb_path(pdb_id: str, pdb_dir: Path | None = None) -> Path | None:
    root = pdb_dir or PDB_DIR
    pid = pdb_id.strip().upper()
    for name in (f"{pid}.pdb", f"{pid.lower()}.pdb"):
        p = root / name
        if p.is_file() and p.stat().st_size > 0:
            return p
    cached = find_cached_pdb(pid)
    if cached:
        return cached
    if _allow_network_fetch:
        return ensure_pdb_cached(pid, allow_network=True)
    return None


def parse_pocket_rule_het(pocket_rule: str) -> str:
    m = _SURROGATE_RE.match((pocket_rule or "").strip())
    return m.group(1).upper() if m else ""


def expected_het_from_entry(entry: dict[str, Any]) -> str:
    for key in ("het_id", "cocrystal_het"):
        v = entry.get(key)
        if isinstance(v, str) and v.strip():
            return v.strip().upper()
    return parse_pocket_rule_het(str(entry.get("pocket_rule") or ""))


def scan_pdb_het_residues(pdb_path: Path) -> list[HetResidue]:
    if not pdb_path.is_file():
        return []
    seen: set[tuple[str, str, int]] = set()
    heavy: dict[tuple[str, str, int], int] = {}

    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM") or len(line) < 22:
            continue
        resname = line[17:20].strip().upper()
        chain = line[21].strip() or " "
        try:
            resseq = int(line[22:26].strip())
        except ValueError:
            continue
        key = (resname, chain, resseq)
        if key not in seen:
            seen.add(key)
            heavy[key] = 0
        element = line[76:78].strip() if len(line) >= 78 else ""
        if element and element.upper() != "H":
            heavy[key] = heavy.get(key, 0) + 1
        elif not element:
            heavy[key] = heavy.get(key, 0) + 1

    out: list[HetResidue] = []
    for key in sorted(seen):
        resname, chain, resseq = key
        ha = heavy.get(key, 0)
        mw = float(ha * 13.0 if ha else 10.0)
        out.append(HetResidue(resname, chain, resseq, ha, mw))
    return out


def is_drug_like_resname(resname: str, *, mw: float | None = None) -> bool:
    rsn = resname.upper()
    if rsn in _BUFFER_IONS or rsn in _GLYCANS:
        return False
    if mw is not None and mw < _MW_CUTOFF:
        return False
    return True


def drug_like_residues(residues: list[HetResidue]) -> list[HetResidue]:
    return [r for r in residues if is_drug_like_resname(r.resname, mw=r.mw_estimate)]


def find_proxy_template(
    catalog_rows: list[dict[str, str]],
    *,
    target_label: str,
    exclude_key: str,
    pdb_dir: Path | None = None,
) -> tuple[dict[str, str], str] | None:
    label = target_label.strip().upper()
    excl = exclude_key.strip().upper()
    candidates: list[tuple[int, dict[str, str], str]] = []

    for row in catalog_rows:
        if (row.get("target_label") or "").strip().upper() != label:
            continue
        key = (row.get("index_key") or "").strip().upper()
        if key == excl:
            continue
        het = (row.get("het_id") or "").strip().upper()
        pdb_id = (row.get("pdb_id") or "").strip().upper()
        if not pdb_id:
            continue
        pdb_path = resolve_pdb_path(pdb_id, pdb_dir)
        if not pdb_path:
            continue
        residues = scan_pdb_het_residues(pdb_path)
        resnames = {r.resname for r in residues}
        effective_het = ""
        if het and het in resnames:
            effective_het = het
        else:
            drug_like = drug_like_residues(residues)
            if drug_like:
                effective_het = max(drug_like, key=lambda r: (r.mw_estimate, r.heavy_atom_count)).resname
        if not effective_het:
            continue
        tier = int(row.get("tier") or "9")
        score = tier * 10 + (0 if (row.get("anchor_mode") or "") == "direct" else 1)
        candidates.append((score, row, effective_het))

    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0])
    row, het = candidates[0][1], candidates[0][2]
    return row, het


def _apply_proxy_template(
    entry: dict[str, Any],
    template: dict[str, str],
    effective_het: str,
    reason: str,
) -> RouteResolution:
    orig = dict(entry)
    het = effective_het.strip().upper()
    pdb_id = (template.get("pdb_id") or "").strip().upper()
    effective = dict(entry)
    effective["pdb_id"] = pdb_id
    effective["het_id"] = het
    effective["cocrystal_het"] = het
    effective["anchor_mode"] = "proxy"
    effective["state"] = "bound_proxy"
    base_rule = str(template.get("pocket_rule") or "")
    if base_rule.startswith("surrogate_"):
        prefix = base_rule.split(":", 1)[0]
        effective["pocket_rule"] = f"{prefix}:{het}"
    else:
        effective["pocket_rule"] = f"surrogate_template_ligand:{het}"
    effective["_route_fallback"] = "DIRECT_TO_CROSS_PDB_PROXY"
    effective["_original_pdb_id"] = orig.get("pdb_id")
    effective["_original_het_id"] = expected_het_from_entry(orig)
    effective["_proxy_source_key"] = (template.get("index_key") or "").strip().upper()
    effective["_route_resolution_status"] = "cross_pdb_proxy"
    if template.get("tier"):
        effective["tier"] = template.get("tier")
    return RouteResolution(
        status="cross_pdb_proxy",
        het_found=False,
        fallback_applied=True,
        reason=reason,
        original_entry=orig,
        effective_entry=effective,
        proxy_source_key=effective["_proxy_source_key"],
        suggested_ccd=het,
    )


def _apply_same_pdb_surrogate(
    entry: dict[str, Any],
    surrogate: HetResidue,
    reason: str,
) -> RouteResolution:
    orig = dict(entry)
    effective = dict(entry)
    effective["het_id"] = surrogate.resname
    effective["cocrystal_het"] = surrogate.resname
    effective["anchor_mode"] = "proxy"
    effective["state"] = "bound_proxy"
    effective["pocket_rule"] = f"surrogate_auto_ligand:{surrogate.resname}"
    effective["_route_fallback"] = "DIRECT_TO_SAME_PDB_SURROGATE"
    effective["_original_het_id"] = expected_het_from_entry(orig)
    effective["_route_resolution_status"] = "same_pdb_surrogate"
    effective["suggested_ccd"] = surrogate.resname
    return RouteResolution(
        status="same_pdb_surrogate",
        het_found=False,
        fallback_applied=True,
        reason=reason,
        original_entry=orig,
        effective_entry=effective,
        suggested_ccd=surrogate.resname,
        drug_like_resnames=[surrogate.resname],
    )


def resolve_anchor_route(
    entry: dict[str, Any],
    pdb_dir: Path | None = None,
    *,
    catalog_rows: list[dict[str, str]] | None = None,
) -> RouteResolution:
    orig = dict(entry)
    pdb_id = str(entry.get("pdb_id") or "").strip().upper()
    expected = expected_het_from_entry(entry)
    anchor_mode = str(entry.get("anchor_mode") or "direct").strip().lower()
    pocket_rule = str(entry.get("pocket_rule") or "")
    index_key = str(entry.get("index_key") or entry.get("_index_key") or "").strip().upper()
    target_label = str(entry.get("target_label") or entry.get("target") or "").strip()

    pdb_path = resolve_pdb_path(pdb_id, pdb_dir) if pdb_id else None
    residues = scan_pdb_het_residues(pdb_path) if pdb_path else []
    resnames = sorted({r.resname for r in residues})
    drug_like = drug_like_residues(residues)
    drug_names = sorted({r.resname for r in drug_like})

    if expected and expected in resnames:
        eff = dict(entry)
        eff.setdefault("_route_resolution_status", "direct_ok")
        return RouteResolution(
            status="direct_ok",
            het_found=True,
            fallback_applied=False,
            reason=f"{expected} present in {pdb_id}",
            original_entry=orig,
            effective_entry=eff,
            pdb_ligand_resnames=resnames,
            drug_like_resnames=drug_names,
        )

    if anchor_mode != "direct" or pocket_rule.startswith("surrogate_"):
        rule_het = parse_pocket_rule_het(pocket_rule)
        check_het = rule_het or expected
        if check_het and check_het in resnames:
            eff = dict(entry)
            eff.setdefault("_route_resolution_status", "proxy_ok")
            return RouteResolution(
                status="proxy_ok",
                het_found=True,
                fallback_applied=False,
                reason=f"proxy het {check_het} present in {pdb_id}",
                original_entry=orig,
                effective_entry=eff,
                pdb_ligand_resnames=resnames,
                drug_like_resnames=drug_names,
            )

    rows = catalog_rows if catalog_rows is not None else load_catalog_rows()

    if drug_like:
        best = max(drug_like, key=lambda r: (r.mw_estimate, r.heavy_atom_count))
        return _apply_same_pdb_surrogate(
            entry,
            best,
            reason=(
                f"catalog het {expected or '?'} missing in {pdb_id}; "
                f"auto surrogate {best.resname} (~{best.mw_estimate:.0f} Da)"
            ),
        )

    template_hit = find_proxy_template(
        rows,
        target_label=target_label,
        exclude_key=index_key,
        pdb_dir=pdb_dir,
    )
    if template_hit:
        template, eff_het = template_hit
        return _apply_proxy_template(
            entry,
            template,
            eff_het,
            reason=(
                f"catalog het {expected or '?'} missing in {pdb_id} and no drug-like HET; "
                f"cross-PDB proxy {template.get('index_key')} "
                f"({template.get('pdb_id')}/{eff_het})"
            ),
        )

    suggested = drug_names[0] if drug_names else (resnames[0] if resnames else None)
    return RouteResolution(
        status="no_route",
        het_found=False,
        fallback_applied=False,
        reason=f"het {expected or '?'} missing in {pdb_id}; no fallback template found",
        original_entry=orig,
        effective_entry=dict(entry),
        pdb_ligand_resnames=resnames,
        drug_like_resnames=drug_names,
        suggested_ccd=suggested,
    )


def write_route_resolution(session_dir: Path, resolution: RouteResolution) -> Path:
    path = session_dir / "route_resolution.json"
    path.write_text(
        json.dumps(resolution.to_json(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


def apply_resolution_to_session_index(session_dir: Path, resolution: RouteResolution) -> None:
    idx_path = session_dir / "session_index.json"
    base: dict[str, Any] = {}
    if idx_path.is_file():
        try:
            base = json.loads(idx_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            base = {}
    for key in SESSION_ROUTE_KEYS:
        if key in resolution.effective_entry:
            base[key] = resolution.effective_entry[key]
    base["_route_resolution"] = {
        "status": resolution.status,
        "fallback_applied": resolution.fallback_applied,
        "reason": resolution.reason,
        "proxy_source_key": resolution.proxy_source_key,
    }
    idx_path.write_text(json.dumps(base, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def merge_session_route_overrides(entry: dict[str, Any], session_dir: Path) -> dict[str, Any]:
    """Overlay session_index / route_resolution onto a master_index entry."""
    out = dict(entry)
    session_path = session_dir / "session_index.json"
    if session_path.is_file():
        try:
            session_entry = json.loads(session_path.read_text(encoding="utf-8"))
            for key in SESSION_ROUTE_KEYS:
                val = session_entry.get(key)
                if val not in (None, ""):
                    out[key] = val
        except (json.JSONDecodeError, OSError):
            pass
    route_path = session_dir / "route_resolution.json"
    if route_path.is_file():
        try:
            route = json.loads(route_path.read_text(encoding="utf-8"))
            eff = route.get("effective_entry") or {}
            for key in SESSION_ROUTE_KEYS:
                val = eff.get(key)
                if val not in (None, ""):
                    out[key] = val
        except (json.JSONDecodeError, OSError):
            pass
    return out


def m4_recover_empty_candidates(
    session_dir: Path,
    entry: dict[str, Any],
    *,
    pdb_dir: Path | None = None,
) -> dict[str, Any] | None:
    """
    Stage 22 recovery when Stage 21 has zero candidates.
    Re-scan cleaned/raw PDB and pick surrogate or cross-pdb proxy metadata.
    """
    pdb_id = str(entry.get("pdb_id") or "").strip().upper()
    cleaned = session_dir / "protein" / "cleaned.pdb"
    scan_path = cleaned if cleaned.is_file() else None
    if scan_path is None and pdb_id:
        scan_path = resolve_pdb_path(pdb_id, pdb_dir)
    if not scan_path or not scan_path.is_file():
        return None

    residues = scan_pdb_het_residues(scan_path)
    drug_like = drug_like_residues(residues)
    if drug_like:
        best = max(drug_like, key=lambda r: (r.mw_estimate, r.heavy_atom_count))
        resolution = _apply_same_pdb_surrogate(
            entry,
            best,
            reason=f"M4 recovery: no Stage-21 candidates; use {best.resname} from {scan_path.name}",
        )
        write_route_resolution(session_dir, resolution)
        apply_resolution_to_session_index(session_dir, resolution)
        return {
            "effective_entry": resolution.effective_entry,
            "candidates": [best.as_candidate(role="target")],
            "resolution": resolution,
        }

    resolution = resolve_anchor_route(entry, pdb_dir)
    if resolution.status == "cross_pdb_proxy":
        write_route_resolution(session_dir, resolution)
        apply_resolution_to_session_index(session_dir, resolution)
        return {
            "effective_entry": resolution.effective_entry,
            "candidates": [],
            "resolution": resolution,
            "needs_rerun_from_m1": True,
        }
    return None
