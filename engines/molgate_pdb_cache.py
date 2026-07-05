"""Shared PDB cache for anchor routing and Module 1 (CMD)."""
from __future__ import annotations

from pathlib import Path

_PKG_ROOT = Path(__file__).resolve().parent
PDB_DIR = _PKG_ROOT / "data" / "pdb"
SESSIONS = _PKG_ROOT / "molgate_sessions"


def ensure_runtime_dirs() -> None:
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    SESSIONS.mkdir(parents=True, exist_ok=True)


def find_cached_pdb(pdb_id: str) -> Path | None:
    """Return engines/data/pdb/{PDBID}.pdb if present."""
    pid = pdb_id.strip().upper()
    for name in (f"{pid}.pdb", f"{pid.lower()}.pdb"):
        p = PDB_DIR / name
        if p.is_file() and p.stat().st_size > 0:
            return p
    return None


def find_pdb_for_session(pdb_id: str, session_dir: Path) -> Path | None:
    pid = pdb_id.strip().upper()
    session_path = session_dir / "pdb" / f"{pid}.pdb"
    if session_path.is_file() and session_path.stat().st_size > 0:
        return session_path
    cached = find_cached_pdb(pid)
    if cached:
        return cached
    if SESSIONS.is_dir():
        for p in sorted(
            SESSIONS.glob(f"*/pdb/{pid}.pdb"),
            key=lambda x: x.stat().st_mtime,
            reverse=True,
        ):
            if p.is_file() and p.stat().st_size > 0:
                return p
    return None


def ensure_pdb_cached(pdb_id: str, *, allow_network: bool = False) -> Path | None:
    """Ensure PDB is in shared cache; optionally download from RCSB."""
    ensure_runtime_dirs()
    pid = pdb_id.strip().upper()
    if not pid:
        return None
    hit = find_cached_pdb(pid)
    if hit:
        return hit
    if not allow_network:
        return None
    from molgate_module1 import _download_pdb_text

    body, _url, _diag = _download_pdb_text(pid)
    if not body:
        return None
    dest = PDB_DIR / f"{pid}.pdb"
    dest.write_text(body, encoding="utf-8")
    return dest if dest.stat().st_size > 0 else None
