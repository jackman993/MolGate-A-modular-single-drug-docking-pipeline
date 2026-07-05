"""
MolGate — Module 3: Protein Preparation (Engine A)
Stage 10: Family classification
Stage 11: State confirmation
Stage 12: Contaminant chain removal
Stage 13: Missing loop assessment
Stage 14: Water removal

全自動層，輸出 cleaned.pdb
Master Index 與 Module 1 共用 master_index.json（單一真相來源）。
可選合併 session 內 stage_03_index.json（與該次 run 對齊）。
依賴：biopython
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa

_PKG_ROOT = Path(__file__).resolve().parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

SCHEMA_VERSION = "0.2"
DEFAULT_INDEX_PATH = _PKG_ROOT / "master_index.json"
INDEX_SCHEMA_FALLBACK = "0.2"


def load_master_index_file(path: Path | None = None) -> dict[str, Any]:
    """與 molgate_module1.load_master_index 對齊的純 JSON 載入（不依賴 RDKit）。"""
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
        entry["_index_schema"] = ver or INDEX_SCHEMA_FALLBACK
        out[k.upper()] = entry
    return out


def load_stage_03(session_dir: Path) -> dict[str, Any] | None:
    p = session_dir / "stage_03_index.json"
    if not p.is_file():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def build_engine_entry(
    master: dict[str, Any],
    target: str,
    session_dir: Path,
) -> dict[str, Any]:
    """
    以 master_index 為底，若有 stage_03 則以該次 run 的 pdb_id / state / family 等覆寫，
    避免與 Module 1 實際寫入的 stage_03 漂移。
    """
    key = target.strip().upper()
    base = master.get(key)
    if base is None:
        raise KeyError(key)
    entry: dict[str, Any] = dict(base)
    s3 = load_stage_03(session_dir)
    if s3:
        for k in (
            "pdb_id",
            "state",
            "family",
            "notes",
            "cocrystal_het",
            "override",
            "pocket_rule",
        ):
            if k in s3 and s3[k] is not None:
                entry[k] = s3[k]
    # 可選欄位（Index 可擴充）；預設與舊內嵌 Index 行為一致
    entry.setdefault("contaminant_chains", [])
    entry.setdefault("keep_water", False)
    if not entry["contaminant_chains"] and isinstance(entry.get("override"), dict):
        oc = entry["override"].get("contaminant_chains")
        if isinstance(oc, list):
            entry["contaminant_chains"] = oc
    if "keep_water" not in base and isinstance(entry.get("override"), dict):
        kw = entry["override"].get("keep_water")
        if isinstance(kw, bool):
            entry["keep_water"] = kw
    return entry

FAMILY_PREP_PROFILE = {
    "cyclooxygenase": "標準準備，無特殊需求",
    "GPCR":           "注意 Asp3.32 override，移除奈米抗體 chain",
    "kinase":         "確認 DFG-in/out 構型，注意 Mg2+ 輔因子",
    "protease":       "確認催化三聯體 His/Asp/Ser 質子化態",
    "nuclear_rec":    "確認配體結合口袋 His 狀態",
    "default":        "標準準備流程",
}

# ── 工具函式 ──────────────────────────────────────────
def save_stage(session_dir: Path, stage: int, name: str, data: dict) -> dict:
    data["_stage"] = stage
    data["_name"] = name
    data["_timestamp"] = datetime.now().isoformat()
    data["_schema_version"] = SCHEMA_VERSION
    path = session_dir / f"stage_{stage:02d}_{name}.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return data

def find_pdb_file(session_dir: Path, pdb_id: str) -> Path:
    pdb_dir = session_dir / "pdb"
    candidates = (list(pdb_dir.glob(f"*{pdb_id.lower()}*")) +
                  list(pdb_dir.glob(f"*{pdb_id.upper()}*")))
    if not candidates:
        raise FileNotFoundError(f"找不到 PDB 檔案：{pdb_id}，請先跑 Module 1")
    return candidates[0]

# ── Stage 10：家族分類 ────────────────────────────────
def stage10_family(entry: dict, session_dir: Path) -> dict:
    print("\n[Stage 10] 蛋白家族分類...")
    family = entry.get("family", "default")
    profile = FAMILY_PREP_PROFILE.get(family, FAMILY_PREP_PROFILE["default"])
    print(f"  ✓ Family: {family}")
    print(f"  ✓ Prep profile: {profile}")
    print(f"  → Status: AUTO")
    return save_stage(session_dir, 10, "family", {
        "status": "AUTO",
        "family": family,
        "prep_profile": profile,
    })

# ── Stage 11：State 確認 ──────────────────────────────
def stage11_state(entry: dict, session_dir: Path) -> dict:
    print("\n[Stage 11] 結構態確認...")
    state = entry.get("state", "unknown")
    notes = entry.get("notes", "")
    status = "WARN" if state == "unknown" else "AUTO"
    if status == "WARN":
        print(f"  ⚠ WARN: state = unknown，結果可信度降低")
    print(f"  ✓ State: {state}")
    if notes:
        print(f"  ✓ Notes: {notes}")
    print(f"  → Status: {status}")
    return save_stage(session_dir, 11, "state", {
        "status": status,
        "state": state,
        "notes": notes,
    })

# ── Stage 12：污染鏈移除 ──────────────────────────────
def stage12_chains(entry: dict, structure, session_dir: Path):
    print("\n[Stage 12] 污染鏈移除...")
    contaminant = entry.get("contaminant_chains", [])
    all_chains = [c.id for m in structure for c in m]
    kept = [c for c in all_chains if c not in contaminant]
    if contaminant:
        print(f"  ⚠ 移除污染鏈：{contaminant}")
    else:
        print(f"  ✓ 無污染鏈")
    print(f"  ✓ 保留 chains：{kept}")
    print(f"  → Status: AUTO")
    r = save_stage(session_dir, 12, "chains", {
        "status": "AUTO",
        "all_chains": all_chains,
        "contaminant_chains": contaminant,
        "kept_chains": kept,
    })
    return r, kept

# ── Stage 13：缺失環評估 ──────────────────────────────
def stage13_loops(structure, kept_chains: list, session_dir: Path) -> dict:
    print("\n[Stage 13] 缺失環評估...")
    missing = []
    for model in structure:
        for chain in model:
            if chain.id not in kept_chains:
                continue
            residues = [r for r in chain if is_aa(r, standard=True)]
            if not residues:
                continue
            seq_ids = sorted([r.get_id()[1] for r in residues])
            for i in range(len(seq_ids) - 1):
                gap = seq_ids[i+1] - seq_ids[i]
                if gap > 1:
                    missing.append({
                        "chain": chain.id,
                        "after": seq_ids[i],
                        "before": seq_ids[i+1],
                        "gap_size": gap - 1,
                    })
    status = "WARN" if missing else "AUTO"
    if missing:
        print(f"  ⚠ 偵測到 {len(missing)} 個缺失片段")
        for seg in missing[:3]:
            print(f"    Chain {seg['chain']}: 殘基 {seg['after']}～{seg['before']} 缺失 {seg['gap_size']} 個")
        if len(missing) > 3:
            print(f"    ...（共 {len(missing)} 個）")
    else:
        print(f"  ✓ 無缺失環，結構完整")
    print(f"  → Status: {status}")
    return save_stage(session_dir, 13, "loops", {
        "status": status,
        "missing_count": len(missing),
        "missing_segments": missing[:20],
    })

# ── Stage 14：去結晶水 ────────────────────────────────
def stage14_water(entry: dict, structure, kept_chains: list, session_dir: Path) -> dict:
    print("\n[Stage 14] 去結晶水...")
    keep_water = entry.get("keep_water", False)
    water_count = sum(
        1 for m in structure for c in m
        if c.id in kept_chains
        for r in c if r.get_resname().strip() in ["HOH", "WAT"]
    )
    if keep_water:
        status = "CONFIRM"
        action = "keep_all"
        print(f"  ⚠ Master Index 標記保留水（共 {water_count} 個）")
    else:
        status = "AUTO"
        action = "remove_all"
        print(f"  ✓ 移除全部結晶水（共 {water_count} 個）")
    print(f"  → Status: {status}")
    return save_stage(session_dir, 14, "water", {
        "status": status,
        "action": action,
        "water_count": water_count,
        "keep_water": keep_water,
    })

# ── PDB 輸出 Select ───────────────────────────────────
class CleanSelect(Select):
    def __init__(self, kept_chains, remove_water):
        self.kept_chains = kept_chains
        self.remove_water = remove_water

    def accept_chain(self, chain):
        return chain.id in self.kept_chains

    def accept_residue(self, residue):
        if self.remove_water and residue.get_resname().strip() in ["HOH", "WAT"]:
            return False
        return True

def _update_manifest_engine_a(
    session_dir: Path,
    cleaned_path: Path,
    flags: list[str],
) -> None:
    """若存在 molgate_run_manifest.json，寫入 cleaned.pdb 與 module3_engine_a 狀態。"""
    try:
        from molgate_manifest import read_manifest, rel_path, write_manifest
    except ImportError:
        return
    mp = session_dir / "molgate_run_manifest.json"
    if not mp.is_file():
        return
    try:
        m = read_manifest(session_dir)
        paths = dict(m.get("paths") or {})
        rmap = dict(paths.get("relative") or {})
        rmap["cleaned_pdb"] = rel_path(session_dir, cleaned_path)
        paths["relative"] = rmap
        paths["cleaned_pdb_absolute"] = str(cleaned_path.resolve())
        ns = dict(paths.get("next_step") or {})
        ns["after_module3_engine_a"] = {
            "read_cleaned_pdb_from": rmap["cleaned_pdb"],
            "hint": "下一步：Engine B（人機決策 15/17 + Flag 16/18 → prep_decisions）",
        }
        paths["next_step"] = ns
        m["paths"] = paths
        pl = dict(m.get("pipeline") or {})
        pl["module3_engine_a"] = {
            "status": "complete",
            "completed_at": datetime.now().isoformat(),
            "stages": ["10", "11", "12", "13", "14"],
        }
        m["pipeline"] = pl
        summ = dict(m.get("summary") or {})
        summ["module3_engine_a_flags"] = flags
        summ["cleaned_pdb_relative"] = rmap["cleaned_pdb"]
        m["summary"] = summ
        write_manifest(session_dir, m)
    except Exception:
        pass


# ── 主流程 ────────────────────────────────────────────
def run_engine_a(
    session_dir: Path,
    target: str,
    *,
    index_path: Path | None = None,
    skip_manifest_target_check: bool = False,
):
    print("=" * 55)
    print("  MolGate — Module 3 Engine A（全自動層）")
    print("  Stage 10/11/12/13/14 → cleaned.pdb")
    print("=" * 55)

    try:
        master = load_master_index_file(index_path)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    try:
        entry = build_engine_entry(master, target, session_dir)
    except KeyError:
        print(f"❌ Target '{target}' 不在 Master Index")
        sys.exit(1)

    if entry.get("allowed") is False:
        print(f"❌ Target '{target}' 在 Master Index 標記為 allowed=false，已中止")
        sys.exit(1)

    if not skip_manifest_target_check:
        try:
            from molgate_manifest import validate_target_matches_manifest

            validate_target_matches_manifest(session_dir, target)
        except FileNotFoundError:
            pass
        except ValueError as e:
            print(f"❌ {e}")
            sys.exit(1)

    pdb_id = (entry.get("pdb_id") or "").upper()
    if not pdb_id:
        print("❌ 條目缺少 pdb_id")
        sys.exit(1)
    print(f"\n  Target：{target.upper()} | PDB：{pdb_id}")

    try:
        pdb_path = find_pdb_file(session_dir, pdb_id)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    print(f"  PDB：{pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(pdb_path))

    flags = []

    try:
        r10 = stage10_family(entry, session_dir)
        r11 = stage11_state(entry, session_dir)
        if r11["status"] == "WARN":
            flags.append("WARN_state_unknown")

        r12, kept_chains = stage12_chains(entry, structure, session_dir)
        r13 = stage13_loops(structure, kept_chains, session_dir)
        if r13["status"] == "WARN":
            flags.append(f"WARN_missing_loops:{r13['missing_count']}")

        r14 = stage14_water(entry, structure, kept_chains, session_dir)

    except Exception as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)

    # 輸出 cleaned.pdb
    protein_dir = session_dir / "protein"
    protein_dir.mkdir(exist_ok=True)
    cleaned_path = protein_dir / "cleaned.pdb"

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(cleaned_path), CleanSelect(kept_chains, r14["action"] == "remove_all"))

    with open(cleaned_path) as f:
        lines = f.readlines()
    atom_count = len([l for l in lines if l.startswith("ATOM")])
    hetatm_count = len([l for l in lines if l.startswith("HETATM")])

    print(f"\n  ✓ cleaned.pdb 完成")
    print(f"  ✓ ATOM: {atom_count} | HETATM: {hetatm_count} | Chains: {kept_chains}")

    _update_manifest_engine_a(session_dir, cleaned_path, flags)

    print("\n" + "=" * 55)
    print("  Engine A 完成 ✓")
    if flags:
        print(f"  Flags：{' | '.join(flags)}")
    print(f"  輸出：{cleaned_path}")
    print(f"  Session log：{session_dir}")
    print("=" * 55)

    return {
        "stage10": r10, "stage11": r11,
        "stage12": r12, "stage13": r13, "stage14": r14,
        "cleaned_pdb": str(cleaned_path),
        "kept_chains": kept_chains,
        "flags": flags,
    }

# ── CLI ──────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MolGate Module 3 Engine A")
    parser.add_argument("--session-dir", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument(
        "--index",
        type=Path,
        default=None,
        help="Master Index JSON（預設：MolGate 目錄下 master_index.json）",
    )
    parser.add_argument(
        "--skip-manifest-target-check",
        action="store_true",
        help="不與 molgate_run_manifest.json 的 target 交叉驗證",
    )
    args = parser.parse_args()

    session_dir = Path(args.session_dir)
    if not session_dir.exists():
        print(f"❌ Session 目錄不存在：{session_dir}")
        sys.exit(1)

    run_engine_a(
        session_dir,
        args.target,
        index_path=args.index,
        skip_manifest_target_check=args.skip_manifest_target_check,
    )
