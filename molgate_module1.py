"""
MolGate — Module 1: Input validation
Stage 01: SMILES standardization (RDKit)
Stage 02: PDB download & validation (RCSB HTTPS + optional Biopython)
Stage 03: Master Index lookup (JSON registry)

Session artifacts: ./molgate_sessions/<session_id>/stage_*.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import socket
import sys
import time
from datetime import datetime
from pathlib import Path

# 同目錄模組（無論從何處執行）
_PKG_ROOT = Path(__file__).resolve().parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import ProxyHandler, Request, build_opener, urlopen

from rdkit import Chem
from rdkit.Chem import Descriptors, rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize

from molgate_manifest import append_registry, finalize_module1_manifest

# ── Paths ─────────────────────────────────────────────────────────
BASE_DIR = _PKG_ROOT / "molgate_sessions"
DEFAULT_INDEX_PATH = _PKG_ROOT / "master_index.json"

SCHEMA_VERSION = "0.3"

# ── Master Index (loaded at runtime) ─────────────────────────────
_MASTER_INDEX: dict[str, Any] = {}


def load_master_index(path: Path | None = None) -> dict[str, Any]:
    """Load Master Index JSON; merge with empty dict. Keys are uppercased for lookup."""
    p = path or DEFAULT_INDEX_PATH
    if not p.is_file():
        raise FileNotFoundError(
            f"找不到 Master Index：{p}\n"
            "請在 MolGate 目錄放置 master_index.json，或使用 --index 指定路徑。"
        )
    raw = json.loads(p.read_text(encoding="utf-8"))
    ver = raw.pop("_schema_version", None)
    raw.pop("_comment", None)
    out: dict[str, Any] = {}
    for k, v in raw.items():
        if not isinstance(v, dict):
            continue
        entry = dict(v)
        entry["_index_schema"] = ver or SCHEMA_VERSION
        out[k.upper()] = entry
    return out


def get_master_index(path: Path | None = None) -> dict[str, Any]:
    global _MASTER_INDEX
    if not _MASTER_INDEX or path is not None:
        _MASTER_INDEX = load_master_index(path)
    return _MASTER_INDEX


# ── Session ───────────────────────────────────────────────────────
def make_session_id(smiles: str, target: str, pdb_id: str) -> str:
    raw = f"{smiles}|{target}|{pdb_id}|{SCHEMA_VERSION}"
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]


def save_stage(session_dir: Path, stage: int, name: str, data: dict) -> dict:
    data = dict(data)
    data["_stage"] = stage
    data["_name"] = name
    data["_timestamp"] = datetime.now().isoformat()
    data["_schema_version"] = SCHEMA_VERSION
    session_dir.mkdir(parents=True, exist_ok=True)
    path = session_dir / f"stage_{stage:02d}_{name}.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return data


def load_stage(session_dir: Path, stage: int, name: str) -> dict | None:
    path = session_dir / f"stage_{stage:02d}_{name}.json"
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


# ── Stage 01 輔助：立體未定義掃描（進 Module 2 前預警）────────────────
def scan_undefined_stereo(mol: Chem.Mol) -> dict[str, Any]:
    """
    偵測「可能具立體但 SMILES 未標記」之手性中心與雙鍵 E/Z。
    在標準化後分子上執行，供下游決策是否先補齊立體再進 Module 2。
    """
    rdmolops.AssignStereochemistry(mol, cleanIt=True, force=True)

    n_chiral = 0
    n_ez = 0
    chiral_atoms: list[int] = []
    ez_bonds: list[list[int]] = []
    method = ""

    try:
        from rdkit.Chem import rdMolDescriptors

        if hasattr(rdMolDescriptors, "CalcNumUnspecifiedAtomStereoCenters"):
            n_chiral = int(rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol))
        if hasattr(rdMolDescriptors, "CalcNumUnspecifiedBondStereoCenters"):
            n_ez = int(rdMolDescriptors.CalcNumUnspecifiedBondStereoCenters(mol))
        if hasattr(rdMolDescriptors, "CalcNumUnspecifiedAtomStereoCenters") or hasattr(
            rdMolDescriptors, "CalcNumUnspecifiedBondStereoCenters"
        ):
            method = "rdMolDescriptors"
    except Exception:
        method = ""

    if not method:
        try:
            centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            for c in centers:
                if len(c) > 1 and c[1] == "?":
                    n_chiral += 1
                    chiral_atoms.append(int(c[0]))
        except Exception:
            pass
        method = "FindMolChiralCenters"

    # 細節與後備計數：FindPotentialStereo（亦補齊舊版缺少的 rd 函式時）
    try:
        from rdkit.Chem import FindPotentialStereo
        from rdkit.Chem.rdchem import StereoSpecified, StereoType

        for entry in FindPotentialStereo(mol):
            if entry.specified != StereoSpecified.Unspecified:
                continue
            if entry.type == StereoType.Atom_Tetrahedral:
                a = int(entry.centeredOn)
                if a not in chiral_atoms:
                    chiral_atoms.append(a)
            elif entry.type == StereoType.Bond_Double:
                try:
                    bidx = int(entry.centeredOn)
                    bnd = mol.GetBondWithIdx(bidx)
                    a1, a2 = bnd.GetBeginAtomIdx(), bnd.GetEndAtomIdx()
                    pair = [a1, a2]
                    if pair not in ez_bonds and [a2, a1] not in ez_bonds:
                        ez_bonds.append(pair)
                except Exception:
                    continue
    except Exception:
        pass

    chiral_atoms = sorted(set(chiral_atoms))
    if n_chiral == 0 and chiral_atoms:
        n_chiral = len(chiral_atoms)
    if n_ez == 0 and ez_bonds:
        n_ez = len(ez_bonds)

    return {
        "n_undefined_chiral_centers": n_chiral,
        "n_undefined_ez_bonds": n_ez,
        "undefined_chiral_atom_indices": chiral_atoms,
        "undefined_ez_bond_atom_pairs": ez_bonds,
        "stereo_scan_method": method or "mixed",
    }


# ── Stage 01: SMILES ─────────────────────────────────────────────
def stage01_smiles(smiles: str, session_dir: Path) -> dict:
    print("\n[Stage 01] SMILES 標準化...")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result = {"status": "FAIL", "reason": "無法解析 SMILES", "input_smiles": smiles}
        save_stage(session_dir, 1, "smiles", result)
        raise ValueError(f"Stage 01 FAIL: {result['reason']}")

    try:
        canonical = rdMolStandardize.StandardizeSmiles(smiles)
    except Exception as e:
        result = {
            "status": "FAIL",
            "reason": f"StandardizeSmiles 失敗: {e}",
            "input_smiles": smiles,
        }
        save_stage(session_dir, 1, "smiles", result)
        raise ValueError(f"Stage 01 FAIL: {result['reason']}")

    mol_std = Chem.MolFromSmiles(canonical)
    if mol_std is None:
        result = {"status": "FAIL", "reason": "標準化後 SMILES 無效", "input_smiles": smiles}
        save_stage(session_dir, 1, "smiles", result)
        raise ValueError(f"Stage 01 FAIL: {result['reason']}")

    mw = Descriptors.MolWt(mol_std)
    n_atoms = mol_std.GetNumAtoms()

    stereo = scan_undefined_stereo(mol_std)

    result: dict[str, Any] = {
        "status": "PASS",
        "input_smiles": smiles,
        "canonical_smiles": canonical,
        "mol_weight": round(mw, 2),
        "num_atoms": n_atoms,
        "is_small_molecule": mw < 900,
        "n_undefined_chiral_centers": stereo["n_undefined_chiral_centers"],
        "n_undefined_ez_bonds": stereo["n_undefined_ez_bonds"],
        "undefined_chiral_atom_indices": stereo["undefined_chiral_atom_indices"],
        "undefined_ez_bond_atom_pairs": stereo["undefined_ez_bond_atom_pairs"],
        "stereo_scan_method": stereo["stereo_scan_method"],
        "recommend_resolve_stereo_before_module2": False,
    }

    warn_msgs: list[str] = []
    if not result["is_small_molecule"]:
        result["status"] = "WARN"
        warn_msgs.append("分子量 > 900 Da，可能不適用標準小分子對接流程")

    if stereo["n_undefined_chiral_centers"] > 0:
        result["status"] = "WARN"
        warn_msgs.append(
            f"未定義手性中心：{stereo['n_undefined_chiral_centers']} 處"
            f"（atom_idx: {stereo['undefined_chiral_atom_indices']}）"
            " — 建議補齊 @ 立體再進 Module 2"
        )

    if stereo["n_undefined_ez_bonds"] > 0:
        result["status"] = "WARN"
        warn_msgs.append(
            f"未定義 E/Z 雙鍵：{stereo['n_undefined_ez_bonds']} 處"
            f"（atom pairs: {stereo['undefined_ez_bond_atom_pairs']}）"
            " — 建議補齊 / 與 \\ 再進 Module 2"
        )

    if stereo["n_undefined_chiral_centers"] > 0 or stereo["n_undefined_ez_bonds"] > 0:
        result["recommend_resolve_stereo_before_module2"] = True

    if warn_msgs:
        result["warnings"] = warn_msgs

    save_stage(session_dir, 1, "smiles", result)
    print(f"  ✓ Canonical SMILES: {canonical}")
    print(f"  ✓ MW: {mw:.1f} Da | Atoms: {n_atoms}")
    print(f"  ✓ 立體掃描: 未定義手性 {stereo['n_undefined_chiral_centers']} | 未定義 E/Z {stereo['n_undefined_ez_bonds']}")
    if stereo["n_undefined_chiral_centers"] > 0:
        print(
            f"  ⚠ WARN: 有未定義手性中心（indices: {stereo['undefined_chiral_atom_indices']}）"
        )
    if stereo["n_undefined_ez_bonds"] > 0:
        print(f"  ⚠ WARN: 有未定義 E/Z（pairs: {stereo['undefined_ez_bond_atom_pairs']}）")
    if result.get("recommend_resolve_stereo_before_module2"):
        print("  → 建議：先處理立體再進 Module 2（仍為 WARN，非強制中止）")
    print(f"  → Status: {result['status']}")
    return result


# ── Stage 02: PDB ─────────────────────────────────────────────────
_PDB_SOURCE_URLS = [
    ("rcsb_download", "https://files.rcsb.org/download/{id}.pdb"),
    ("rcsb_view", "https://files.rcsb.org/view/{id}.pdb"),
    ("pdbe_ent", "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{id}.ent"),
]


def _validate_pdb_content(text: str, pdb_id: str) -> dict:
    lines = text.splitlines()
    has_atom = any(l.startswith(("ATOM  ", "HETATM")) for l in lines)
    has_cryst = any(l.startswith("CRYST1") for l in lines)
    atom_n = sum(1 for l in lines if l.startswith("ATOM"))
    het_n = sum(1 for l in lines if l.startswith("HETATM"))
    chains = sorted({l[21] for l in lines if len(l) > 21 and l.startswith("ATOM")})
    if not has_atom and not any(l.startswith("HETATM") for l in lines):
        return {"ok": False, "reason": "檔案無 ATOM/HETATM 記錄"}
    if len(text.strip()) < 200:
        return {"ok": False, "reason": "檔案過短，可能非有效 PDB"}
    return {
        "ok": True,
        "pdb_id": pdb_id.upper(),
        "atom_count": atom_n,
        "hetatm_count": het_n,
        "chains": chains,
        "has_hetatm": het_n > 0,
        "has_cryst1": has_cryst,
    }


def _format_network_reason(err: Exception) -> str:
    """Normalize network exceptions to readable diagnostics."""
    if isinstance(err, HTTPError):
        return f"HTTP {err.code}: {err.reason}"
    if isinstance(err, URLError):
        r = err.reason
        if isinstance(r, OSError):
            if getattr(r, "winerror", None) == 10061:
                return "連線被拒絕 (WinError 10061，可能是防火牆/代理/公司網路限制)"
            return f"網路錯誤: {r}"
        return f"網路錯誤: {r}"
    return str(err)


def _download_pdb_text(pid: str, *, max_retries: int = 2) -> tuple[str | None, str | None, list[str]]:
    """
    Try multiple mirrors with retries.
    Return (body, source_url, diagnostics).
    """
    diagnostics: list[str] = []
    for source_name, pattern in _PDB_SOURCE_URLS:
        url = pattern.format(id=pid)
        for attempt in range(1, max_retries + 1):
            req = Request(url, headers={"User-Agent": "MolGate/0.2 (Module1; research)"})
            for mode in ("proxy_env", "direct"):
                try:
                    if mode == "proxy_env":
                        with urlopen(req, timeout=45) as resp:
                            body = resp.read().decode("utf-8", errors="replace")
                    else:
                        opener = build_opener(ProxyHandler({}))
                        with opener.open(req, timeout=45) as resp:
                            body = resp.read().decode("utf-8", errors="replace")
                    if "404" in body[:200] or len(body) < 200:
                        diagnostics.append(f"{source_name} try#{attempt} {mode}: 內容異常或為空")
                        continue
                    return body, url, diagnostics
                except Exception as e:  # noqa: BLE001
                    diagnostics.append(f"{source_name} try#{attempt} {mode}: {_format_network_reason(e)}")
            time.sleep(0.6 * attempt)
    return None, None, diagnostics


def _find_local_pdb(pid: str, session_dir: Path) -> Path | None:
    """Find local PDB by priority: current session, shared repo, previous sessions."""
    local_candidates = [
        session_dir / "pdb" / f"{pid}.pdb",
        _PKG_ROOT / "pdb" / f"{pid}.pdb",
        _PKG_ROOT / "tool" / "pdb" / f"{pid}.pdb",
    ]
    for c in local_candidates:
        if c.is_file():
            return c

    # Reuse from previous sessions if available.
    if BASE_DIR.is_dir():
        for p in sorted(BASE_DIR.glob(f"*/pdb/{pid}.pdb"), key=lambda x: x.stat().st_mtime, reverse=True):
            if p.is_file():
                return p
    return None


def stage02_pdb(pdb_id: str, session_dir: Path, *, allow_network: bool = False) -> dict:
    print(f"\n[Stage 02] PDB 載入與驗證：{pdb_id}...")

    pid = pdb_id.strip().upper()
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}", pid):
        result = {"status": "FAIL", "pdb_id": pid, "reason": "PDB ID 格式應為 4 字元（例：7DFP）"}
        save_stage(session_dir, 2, "pdb", result)
        raise ValueError(f"Stage 02 FAIL: {result['reason']}")

    pdb_dir = session_dir / "pdb"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    out_path = pdb_dir / f"{pid}.pdb"
    url = None

    # Offline-first: prefer local prepared PDB sources.
    src = _find_local_pdb(pid, session_dir)
    if src is not None:
        body = src.read_text(encoding="utf-8", errors="replace")
        v = _validate_pdb_content(body, pid)
        if not v.get("ok"):
            result = {
                "status": "FAIL",
                "pdb_id": pid,
                "pdb_path": str(src.resolve()),
                "reason": f"本地 PDB 驗證失敗: {v.get('reason')}",
            }
            save_stage(session_dir, 2, "pdb", result)
            raise ValueError(f"Stage 02 FAIL: {result['reason']}")

        out_path.write_text(body, encoding="utf-8")
        result = {
            "status": "PASS",
            "pdb_id": pid,
            "pdb_path": str(out_path.resolve()),
            "source": "local",
            "source_path": str(src.resolve()),
            "download_url": None,
            "atom_count": v["atom_count"],
            "hetatm_count": v["hetatm_count"],
            "chains": v["chains"],
            "has_hetatm": v["has_hetatm"],
            "has_cryst1": v["has_cryst1"],
        }
        save_stage(session_dir, 2, "pdb", result)
        print(f"  ✓ 已由本地載入：{src}")
        print(f"  ✓ 複製到 session：{out_path}")
        print(
            f"  ✓ ATOM: {v['atom_count']} | HETATM: {v['hetatm_count']} | Chains: {v['chains']}"
        )
        print(f"  → Status: {result['status']}")
        return result

    if not allow_network:
        result = {
            "status": "FAIL",
            "pdb_id": pid,
            "reason": (
                "離線模式未找到本地 PDB。"
                "請先放置於 ./pdb 或 ./tool/pdb，或先跑過一個含此 PDB 的 session，"
                "或改用 --allow-network。"
            ),
        }
        save_stage(session_dir, 2, "pdb", result)
        raise ValueError(f"Stage 02 FAIL: {result['reason']}")

    body, url, diagnostics = _download_pdb_text(pid)
    if body is None:
        host_diag = ""
        try:
            socket.getaddrinfo("files.rcsb.org", 443)
        except Exception as e:  # noqa: BLE001
            host_diag = f" | DNS 檢查失敗(files.rcsb.org): {e}"
        proxy_hint = (
            os.environ.get("HTTPS_PROXY")
            or os.environ.get("https_proxy")
            or os.environ.get("HTTP_PROXY")
            or os.environ.get("http_proxy")
        )
        proxy_diag = f" | 代理設定: {proxy_hint}" if proxy_hint else " | 代理設定: <none>"
        result = {
            "status": "FAIL",
            "pdb_id": pid,
            "reason": "；".join(diagnostics) + host_diag + proxy_diag,
        }
        save_stage(session_dir, 2, "pdb", result)
        raise ValueError(f"Stage 02 FAIL: {result['reason']}")

    out_path.write_text(body, encoding="utf-8")
    v = _validate_pdb_content(body, pid)
    if not v.get("ok"):
        result = {"status": "FAIL", "pdb_id": pid, "pdb_path": str(out_path), "reason": v.get("reason")}
        save_stage(session_dir, 2, "pdb", result)
        raise ValueError(f"Stage 02 FAIL: {result['reason']}")

    result = {
        "status": "PASS",
        "pdb_id": pid,
        "pdb_path": str(out_path.resolve()),
        "source": "network",
        "source_path": None,
        "download_url": url,
        "network_attempts": diagnostics,
        "atom_count": v["atom_count"],
        "hetatm_count": v["hetatm_count"],
        "chains": v["chains"],
        "has_hetatm": v["has_hetatm"],
        "has_cryst1": v["has_cryst1"],
    }
    save_stage(session_dir, 2, "pdb", result)
    print(f"  ✓ 已寫入：{out_path}")
    print(
        f"  ✓ ATOM: {v['atom_count']} | HETATM: {v['hetatm_count']} | Chains: {v['chains']}"
    )
    print(f"  → Status: {result['status']}")
    return result


# ── Stage 03: Master Index ────────────────────────────────────────
def stage03_index(
    target: str,
    pdb_id: str,
    session_dir: Path,
    index_path: Path | None = None,
) -> dict:
    print(f"\n[Stage 03] Master Index 查詢：target={target}...")

    idx = get_master_index(index_path)
    key = target.strip().upper()
    entry = idx.get(key)

    if entry is None:
        result = {
            "status": "FAIL",
            "target": key,
            "reason": f"Target '{key}' 不在 Master Index，請先於 master_index.json 策編",
        }
        save_stage(session_dir, 3, "index", result)
        raise ValueError(f"Stage 03 FAIL: {result['reason']}")

    if not entry.get("allowed", True):
        result = {
            "status": "BLOCKED",
            "target": key,
            "reason": f"此 Target 已禁用：{entry.get('notes', '')}",
        }
        save_stage(session_dir, 3, "index", result)
        raise ValueError(f"Stage 03 BLOCKED: {result['reason']}")

    rec_pdb = str(entry.get("pdb_id", "")).upper()
    in_pdb = pdb_id.strip().upper()

    if rec_pdb and rec_pdb != in_pdb:
        status = "WARN"
        warning = f"建議 PDB 為 Master Index 指定之 {rec_pdb}，目前輸入為 {in_pdb}"
    else:
        status = "PASS"
        warning = None

    result = {
        "status": status,
        "target": key,
        "pdb_id": rec_pdb or in_pdb,
        "input_pdb": in_pdb,
        "state": entry.get("state"),
        "family": entry.get("family"),
        "pocket_rule": entry.get("pocket_rule"),
        "cocrystal_het": entry.get("cocrystal_het"),
        "override": entry.get("override") or {},
        "notes": entry.get("notes"),
    }
    if warning:
        result["warning"] = warning

    save_stage(session_dir, 3, "index", result)
    print(f"  ✓ Target: {key} | Family: {entry.get('family')} | State: {entry.get('state')}")
    print(f"  ✓ Pocket rule: {entry.get('pocket_rule')} | HET: {entry.get('cocrystal_het')}")
    if entry.get("override"):
        print(f"  ✓ Override: {entry.get('override')}")
    print(f"  → Status: {result['status']}")
    if warning:
        print(f"  ⚠ {warning}")
    return result


# ── Module 1 runner ───────────────────────────────────────────────
def run_module1(
    smiles: str,
    target: str,
    pdb_id: str,
    index_path: Path | None = None,
    session_id: str | None = None,
    drug_id: str | None = None,
    drug_label: str | None = None,
    allow_network: bool = False,
) -> dict | None:
    print("=" * 58)
    print("  MolGate — Module 1: Input validation")
    print("=" * 58)

    BASE_DIR.mkdir(parents=True, exist_ok=True)
    sid = session_id or make_session_id(smiles, target, pdb_id)
    session_dir = BASE_DIR / sid
    session_dir.mkdir(parents=True, exist_ok=True)

    meta = {
        "session_id": sid,
        "module": 1,
        "created": datetime.now().isoformat(),
        "smiles_input": smiles,
        "target": target,
        "pdb_id": pdb_id,
        "index_path": str((index_path or DEFAULT_INDEX_PATH).resolve()),
        "drug_id": drug_id,
        "drug_label": drug_label,
    }
    (session_dir / "session_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print(f"\n  Session ID: {sid}")
    print(f"  Session dir: {session_dir}")

    try:
        r1 = stage01_smiles(smiles, session_dir)
        r2 = stage02_pdb(pdb_id, session_dir, allow_network=allow_network)
        r3 = stage03_index(target, pdb_id, session_dir, index_path=index_path)
    except ValueError as e:
        print(f"\n❌ Module 1 中止: {e}")
        return None

    finalize_module1_manifest(
        session_dir,
        meta,
        r1,
        r2,
        r3,
        drug_id=drug_id,
        drug_label=drug_label,
    )
    append_registry(
        BASE_DIR,
        {
            "event": "module1_complete",
            "session_id": sid,
            "session_dir": str(session_dir.resolve()),
            "drug_id": drug_id,
            "drug_label": drug_label,
            "target": target,
            "pdb_id": r2.get("pdb_id") or pdb_id,
            "timestamp": datetime.now().isoformat(),
        },
    )

    print("\n" + "=" * 58)
    print("  Module 1 完成")
    print(f"  日誌目錄: {session_dir}")
    print(f"  Manifest: {session_dir / 'molgate_run_manifest.json'}")
    print("=" * 58)

    return {"session_id": sid, "session_dir": str(session_dir), "stage01": r1, "stage02": r2, "stage03": r3}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="MolGate Module 1 — input validation")
    p.add_argument("--smiles", required=True, help="輸入 SMILES")
    p.add_argument("--target", required=True, help="Master Index 鍵名，例：D2R_7DFP、COX2")
    p.add_argument("--pdb-id", required=True, dest="pdb_id", help="RCSB PDB ID，例：7DFP")
    p.add_argument("--index", type=Path, default=None, help="master_index.json 路徑")
    p.add_argument("--drug-id", default=None, help="藥物管理用 ID（例：M01、IBU-001）")
    p.add_argument("--drug-label", default=None, help="藥物顯示名稱（選填）")
    p.add_argument(
        "--allow-network",
        action="store_true",
        help="允許 Stage 02 連線下載 PDB（預設離線，只用本地資料）",
    )
    args = p.parse_args(argv)

    r = run_module1(
        args.smiles,
        args.target,
        args.pdb_id,
        index_path=args.index,
        drug_id=args.drug_id,
        drug_label=args.drug_label,
        allow_network=args.allow_network,
    )
    return 0 if r is not None else 1


if __name__ == "__main__":
    sys.exit(main())
