"""
MolGate — Module 4: Pocket Definition (Engine A)
Stage 21: HETATM scan & filter

輸入：protein/cleaned.pdb（Module 3 Engine A）、prep_decisions.json（須 Module 3 Engine B 完成並
      ready_for_engine_e=true，與 Module 3 Engine C 前置條件一致）
輸出：stage_21_hetatm.json（候選配體列表；含 ligand_role=target/candidate）
規則：Master Index 之 het_id／cocrystal_het 白名單優先；糖基化／輔因子黑名單（輔因子可於 Index
      設 keep_cofactors 保留）；其餘套用緩衝離子黑名單與 MW 門檻。
依賴：biopython, rdkit；讀 Index 時與 molgate_module3_engineA 共用載入邏輯。
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from Bio.PDB import PDBParser

_PKG_ROOT = Path(__file__).resolve().parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

SCHEMA_VERSION = "0.2"

# ── 黑名單：緩衝劑、離子、溶劑 ───────────────────────
BLACKLIST = {
    "HOH","WAT","GOL","EDO","PEG","MPD","DMS","FMT",
    "ACT","TRS","MES","EPE","BME","DTT","TCE","IMD",
    "SO4","PO4","CIT","ACY","EOH","IPA","MOH","BOG",
    "CL","NA","MG","ZN","CA","K","FE","MN","CU","NI",
    "CO","CD","HG","PT","AU","AG","SR","BA","CS","RB",
    "IOD","BR","F","LI","AL","SE",
}

# ── 糖基化／寡糖常見殘基（非口袋藥物配體）────────────────
GLYCAN_BLACKLIST = {
    "NAG", "NGA", "BMA", "MAN", "GAL", "GLC", "GCS", "MAL", "MAB", "MBG", "BGC",
    "FUC", "FUL", "RAM", "LFR", "GLO", "XYP", "XYS", "XYF", "XUL",
    "SIA", "SLB", "A2G", "M6D", "MMA", "MGL", "NGK", "IDR",
}

# ── 常見輔因子／核苷酸／金屬簇（非口袋配體；除非 Index 標記 keep_cofactors）──
COFACTOR_BLACKLIST = {
    "HEM", "HEA", "HEB", "HEC", "HAS", "CHL", "BCL", "CLA", "FAD", "FMN",
    "NAD", "NAP", "NAI", "NAO", "NAM", "NDP", "NDG", "NT2", "NAJ", "NAC",
    "ATP", "ADP", "AMP", "ANP", "GTP", "GDP", "GMP", "GNP",
    "UTP", "UDP", "UMP", "CTP", "CDP", "CMP", "TTP", "TDP", "TMP",
    "COA", "SAM", "SAH", "ACO", "SCA",
    "SF4", "FES", "F3S", "CLF", "CFM", "FSO", "FS4",
    "PLP", "PLH", "PMP", "PSP", "MGD", "GTN", "Q10", "U10", "3CO", "CUA", "CU1",
}

MW_CUTOFF = 150.0  # Da 以下排除（Master Index 標記之目標配體不受此限）

# ── 工具函式 ──────────────────────────────────────────
def save_stage(session_dir: Path, stage: int, name: str, data: dict[str, Any]) -> dict[str, Any]:
    data = dict(data)
    data["_stage"] = stage
    data["_name"] = name
    data["_timestamp"] = datetime.now().isoformat()
    data["_schema_version"] = SCHEMA_VERSION
    path = session_dir / f"stage_{stage:02d}_{name}.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return data

def load_prep_decisions(session_dir: Path) -> dict[str, Any]:
    path = session_dir / "prep_decisions.json"
    if not path.is_file():
        raise FileNotFoundError(
            "找不到 prep_decisions.json，請先完成 Module 3 Engine B（人機決策 → ready_for_engine_e）"
        )
    return json.loads(path.read_text(encoding="utf-8"))


def validate_module3_prerequisites(
    session_dir: Path,
    target: str | None,
    *,
    skip_ready_check: bool = False,
    skip_manifest_target_check: bool = False,
) -> dict[str, Any]:
    """
    與 molgate_module3_engineC 對齊：須 prep_decisions、ready_for_engine_e；
    若提供 --target 則與 manifest / prep_decisions.target 交叉驗證。
    """
    if target and not skip_manifest_target_check:
        try:
            from molgate_manifest import validate_target_matches_manifest

            validate_target_matches_manifest(session_dir, target)
        except FileNotFoundError:
            pass
        except ValueError as e:
            print(f"❌ {e}")
            sys.exit(1)

    try:
        prep = load_prep_decisions(session_dir)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    if not skip_ready_check and not prep.get("ready_for_engine_e"):
        print(
            "❌ prep_decisions.json 尚未標記 ready_for_engine_e=true，請先完成 Module 3 Engine B"
        )
        sys.exit(1)

    if target:
        pt = (prep.get("target") or "").strip().upper()
        if pt and pt != target.strip().upper():
            print(
                f"❌ prep_decisions.target={prep.get('target')!r} 與 --target={target!r} 不一致"
            )
            sys.exit(1)

    pdb_id = (prep.get("pdb_id") or "").strip().upper()
    if not pdb_id:
        print("❌ prep_decisions 缺少 pdb_id（請先完成 Module 3 Engine B）")
        sys.exit(1)

    return prep


def load_cleaned_pdb(session_dir: Path) -> Path:
    p = session_dir / "protein" / "cleaned.pdb"
    if not p.is_file():
        raise FileNotFoundError(
            "找不到 protein/cleaned.pdb，請先完成 Module 3 Engine A（產出 cleaned.pdb）"
        )
    return p


def resolve_index_path(session_dir: Path, cli_index: Path | None) -> Path:
    if cli_index is not None:
        return cli_index
    meta = session_dir / "session_meta.json"
    if meta.is_file():
        try:
            data = json.loads(meta.read_text(encoding="utf-8"))
            p = data.get("index_path")
            if isinstance(p, str) and p.strip():
                return Path(p)
        except (json.JSONDecodeError, OSError):
            pass
    mp = session_dir / "molgate_run_manifest.json"
    if mp.is_file():
        try:
            data = json.loads(mp.read_text(encoding="utf-8"))
            p = data.get("index_path")
            if isinstance(p, str) and p.strip():
                return Path(p)
        except (json.JSONDecodeError, OSError):
            pass
    return _PKG_ROOT / "master_index.json"


def load_target_index_entry(
    session_dir: Path,
    prep: dict[str, Any],
    index_path: Path | None,
) -> dict[str, Any]:
    """Master Index 條目（與 Module 3 對齊）；失敗時以 prep 的 cocrystal_het 等為後備。"""
    target = (prep.get("target") or "").strip()
    if not target:
        return {}
    resolved = resolve_index_path(session_dir, index_path)
    try:
        from molgate_module3_engineA import build_engine_entry, load_master_index_file

        master = load_master_index_file(resolved)
        return build_engine_entry(master, target, session_dir)
    except (FileNotFoundError, KeyError, OSError):
        out: dict[str, Any] = {}
        for k in ("cocrystal_het", "het_id"):
            v = prep.get(k)
            if isinstance(v, str) and v.strip():
                out[k] = v.strip()
        ov = prep.get("override")
        if isinstance(ov, dict):
            out["override"] = ov
        return out


def get_primary_het_ids(entry: dict[str, Any], prep: dict[str, Any]) -> set[str]:
    """Master Index / prep 標記之共晶配體 ID（het_id 或 cocrystal_het）。"""
    s: set[str] = set()
    for src in (entry, prep):
        if not isinstance(src, dict):
            continue
        for key in ("het_id", "cocrystal_het"):
            v = src.get(key)
            if isinstance(v, str) and v.strip():
                s.add(v.strip().upper())
            elif isinstance(v, list):
                for x in v:
                    if isinstance(x, str) and x.strip():
                        s.add(x.strip().upper())
    return s


def get_keep_cofactors(entry: dict[str, Any]) -> set[str]:
    """除非列於此，否則套用 COFACTOR_BLACKLIST。可寫在條目或 override.keep_cofactors。"""
    s: set[str] = set()
    for key in ("keep_cofactors", "retain_cofactors", "cofactor_whitelist"):
        v = entry.get(key)
        if isinstance(v, list):
            s |= {x.strip().upper() for x in v if isinstance(x, str) and x.strip()}
    ov = entry.get("override")
    if isinstance(ov, dict):
        for key in ("keep_cofactors", "retain_cofactors", "cofactor_whitelist"):
            v = ov.get(key)
            if isinstance(v, list):
                s |= {x.strip().upper() for x in v if isinstance(x, str) and x.strip()}
    return s


def sort_stage21_candidates(cands: list[dict[str, Any]]) -> list[dict[str, Any]]:
    targets = [c for c in cands if c.get("ligand_role") == "target"]
    rest = [c for c in cands if c.get("ligand_role") != "target"]
    targets.sort(key=lambda c: (c["chain"], c["resseq"]))
    rest.sort(key=lambda c: (-float(c["mw_estimate"]), c["chain"], c["resseq"]))
    return targets + rest


def get_het_atoms(structure, resname: str, chain_id: str, resseq: int) -> list:
    atoms = []
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for res in chain:
                if (res.get_resname().strip() == resname and
                        res.get_id()[1] == resseq):
                    atoms.extend(res.get_atoms())
    return atoms

def estimate_mw_from_pdb(resname: str) -> float:
    """
    從 PDB HETATM 估算 MW。
    優先用 RDKit CCD 查詢，失敗則用原子數粗估。
    """
    try:
        from rdkit.Chem import AllChem, Descriptors
        from rdkit.Chem.rdmolfiles import MolFromPDBBlock
        # 嘗試用 CCD SMILES 查詢（需要網路或本地 CCD）
        # 這裡用原子計數粗估作為 fallback
    except ImportError:
        pass

    # 粗估：平均有機分子每原子約 13 Da
    # 實際上我們用 PDB 結構裡的原子數來估
    return None  # 由呼叫端用原子數估算

# ── Stage 21：HETATM 掃描過濾 ─────────────────────────
def stage21_hetatm(
    session_dir: Path,
    prep: dict[str, Any],
    index_path: Path | None,
) -> dict[str, Any]:
    print("\n[Stage 21] HETATM 掃描過濾...")

    idx_entry = load_target_index_entry(session_dir, prep, index_path)
    primary_ids = get_primary_het_ids(idx_entry, prep)
    keep_cof = get_keep_cofactors(idx_entry)

    print(f"  → Master Index 目標 het_id / cocrystal_het：{sorted(primary_ids) or '（無）'}")
    if keep_cof:
        print(f"  → 輔因子保留（keep_cofactors）：{sorted(keep_cof)}")

    cleaned_path = load_cleaned_pdb(session_dir)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(cleaned_path))

    # 掃描所有 HETATM 殘基
    raw_het = []
    seen = set()

    for model in structure:
        for chain in model:
            for res in chain:
                hetflag, resseq, _ = res.get_id()
                resname = res.get_resname().strip().upper()

                # 只處理 HETATM（hetflag 非空白）
                if not hetflag.strip():
                    continue

                key = (resname, chain.id, resseq)
                if key in seen:
                    continue
                seen.add(key)

                atoms = list(res.get_atoms())
                atom_count = len(atoms)

                # 用原子數估算 MW（粗估，平均有機分子每非氫原子約 13 Da）
                # 更準確的做法是查 CCD，這裡用保守估算
                heavy_atoms = [a for a in atoms if a.element and a.element != 'H']
                mw_estimate = len(heavy_atoms) * 13.0 if heavy_atoms else atom_count * 10.0

                raw_het.append({
                    "resname": resname,
                    "chain": chain.id,
                    "resseq": resseq,
                    "atom_count": atom_count,
                    "heavy_atom_count": len(heavy_atoms),
                    "mw_estimate": round(mw_estimate, 1),
                })

    print(f"  ✓ 掃描到 HETATM 殘基：{len(raw_het)} 個")

    # 過濾：白名單（Index het）優先 → 糖 → 輔因子 → 緩衝/離子 → MW
    excluded_blacklist: list[str] = []
    excluded_glycan: list[str] = []
    excluded_cofactor: list[str] = []
    excluded_mw: list[str] = []
    candidates: list[dict[str, Any]] = []

    for het in raw_het:
        resname = het["resname"]
        mw = float(het["mw_estimate"])
        rsn = resname.upper()

        if rsn in primary_ids:
            row = dict(het)
            row["ligand_role"] = "target"
            candidates.append(row)
            continue

        if rsn in GLYCAN_BLACKLIST:
            excluded_glycan.append(rsn)
            continue

        if rsn in COFACTOR_BLACKLIST and rsn not in keep_cof:
            excluded_cofactor.append(rsn)
            continue

        if rsn in BLACKLIST:
            excluded_blacklist.append(rsn)
            continue

        if mw < MW_CUTOFF:
            excluded_mw.append(f"{resname}({mw:.0f}Da)")
            continue

        row = dict(het)
        row["ligand_role"] = "candidate"
        candidates.append(row)

    candidates = sort_stage21_candidates(candidates)

    print(f"  ✓ 目標配體（Index 白名單）：{sum(1 for c in candidates if c.get('ligand_role') == 'target')} 個")
    print(f"  ✓ 糖基化殘基排除：{len(excluded_glycan)} 個")
    print(f"  ✓ 輔因子排除：{len(excluded_cofactor)} 個")
    print(f"  ✓ 緩衝/離子/溶劑黑名單排除：{len(excluded_blacklist)} 個")
    print(f"  ✓ MW < {MW_CUTOFF} Da 排除：{len(excluded_mw)} 個")
    print(f"  ✓ 候選配體（目標 + 其餘）：{len(candidates)} 個")

    for c in candidates:
        tag = "【目標】" if c.get("ligand_role") == "target" else ""
        print(
            f"    {tag}{c['resname']} chain {c['chain']} res {c['resseq']} "
            f"(~{c['mw_estimate']:.0f} Da, {c['heavy_atom_count']} heavy atoms)"
        )

    if not candidates:
        status = "WARN"
        print(f"  ⚠ WARN: 無候選配體，可能是 apo 結構或 HET 已被移除")
    else:
        status = "PASS"

    print(f"  → Status: {status}")

    result = {
        "status": status,
        "raw_hetatm_count": len(raw_het),
        "excluded_glycan_count": len(excluded_glycan),
        "excluded_cofactor_count": len(excluded_cofactor),
        "excluded_blacklist_count": len(excluded_blacklist),
        "excluded_mw_count": len(excluded_mw),
        "candidate_count": len(candidates),
        "candidates": candidates,
        "mw_cutoff": MW_CUTOFF,
        "filter_rules": {
            "primary_het_ids": sorted(primary_ids),
            "keep_cofactors": sorted(keep_cof),
            "glycan_blacklist_size": len(GLYCAN_BLACKLIST),
            "cofactor_blacklist_size": len(COFACTOR_BLACKLIST),
        },
        "module3_prep": {
            "pdb_id": (prep.get("pdb_id") or "").strip().upper(),
            "target": prep.get("target"),
            "ready_for_engine_e": bool(prep.get("ready_for_engine_e")),
        },
    }
    return save_stage(session_dir, 21, "hetatm", result)

# ── 主流程 ────────────────────────────────────────────
def run_engine_a(
    session_dir: Path,
    *,
    target: str | None = None,
    index_path: Path | None = None,
    skip_ready_check: bool = False,
    skip_manifest_target_check: bool = False,
):
    print("=" * 55)
    print("  MolGate — Module 4 Engine A")
    print("  Stage 21: HETATM 掃描過濾")
    print("=" * 55)

    prep = validate_module3_prerequisites(
        session_dir,
        target,
        skip_ready_check=skip_ready_check,
        skip_manifest_target_check=skip_manifest_target_check,
    )
    pdb_id = (prep.get("pdb_id") or "").strip().upper()
    tgt = (prep.get("target") or target or "").strip().upper()
    print(f"\n  Module 3 銜接：PDB={pdb_id}" + (f" | Target={tgt}" if tgt else ""))

    try:
        result = stage21_hetatm(session_dir, prep, index_path)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Pipeline stopped: {e}")
        sys.exit(1)

    print("\n" + "=" * 55)
    if result["status"] == "WARN":
        print("  Engine A 完成（WARN — 無候選配體）")
        print("  → 若為 apo 結構，請在 Master Index 標記 pocket_rule: manual_coords")
    else:
        print(f"  Engine A 完成 ✓ — {result['candidate_count']} 個候選配體")
    print(f"  輸出：{session_dir / 'stage_21_hetatm.json'}")
    print("=" * 55)

    return result

# ── CLI ──────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="MolGate Module 4 Engine A")
    ap.add_argument("--session-dir", required=True, help="Session 目錄路徑")
    ap.add_argument(
        "--target",
        default=None,
        help="標靶 ID（建議與 Module 3 一致；若提供則與 prep_decisions / manifest 交叉驗證）",
    )
    ap.add_argument(
        "--index",
        default=None,
        help="master_index.json 路徑（預設：session_meta / manifest 或 MolGate 目錄下）",
    )
    ap.add_argument(
        "--skip-ready-check",
        action="store_true",
        help="不檢查 ready_for_engine_e（僅除錯用）",
    )
    ap.add_argument(
        "--skip-manifest-target-check",
        action="store_true",
        help="不與 molgate_run_manifest.json 的 target 交叉驗證",
    )
    args = ap.parse_args()

    sd = Path(args.session_dir)
    if not sd.is_dir():
        print(f"❌ Session 目錄不存在：{sd}")
        sys.exit(1)

    idx = Path(args.index) if args.index else None
    run_engine_a(
        sd,
        target=args.target,
        index_path=idx,
        skip_ready_check=args.skip_ready_check,
        skip_manifest_target_check=args.skip_manifest_target_check,
    )
