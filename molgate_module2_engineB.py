"""
MolGate — Module 2: Ligand Enumeration (Engine B)
Stage 04: Protonation state (Dimorphite-DL @ pH 7.4)
Stage 05: Tautomer (RDKit canonical)

依賴：
  pip install dimorphite-dl rdkit
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

_PKG_ROOT = Path(__file__).resolve().parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

from molgate_manifest import (
    append_registry,
    finalize_module2_engine_b_manifest,
    try_load_manifest_or_rebuild,
    validate_target_matches_manifest,
)

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# ── 嘗試載入 Dimorphite-DL ───────────────────────────
try:
    from dimorphite_dl import protonate_smiles
    DIMORPHITE_OK = True
except ImportError:
    DIMORPHITE_OK = False
    print("[WARN] dimorphite-dl 未安裝，Stage 04 將降級為原始 SMILES")

# Stage 04/05 產物版本（與 Module 1 的 stage_01 schema 分開）
SCHEMA_VERSION = "0.2"
# Module 1 Stage 01 目前可能為 0.2、0.3…
MODULE1_STAGE01_SCHEMAS = frozenset({"0.2", "0.3"})

# ── Master Index risk tier 查詢 ───────────────────────
# 實際使用時從 JSON 或 DB 讀取，這裡示範內嵌
LIGAND_RISK = {
    "COX2": {
        "stage04_tier": "A",
        "stage04_rule": "neutral_only",
        "stage04_note": "羧酸藥，pH 7.4 主要去質子化，單一形態足夠",
        "stage05_tier": "A",
        "stage05_note": "無高風險互變核心",
    },
    "D2R": {
        "stage04_tier": "C",
        "stage04_rule": "manual_select",
        "stage04_note": "CNS 胺基藥，pKa 接近 7.4，需手選質子化態",
        "stage05_tier": "A",
        "stage05_note": "無高風險互變核心",
    },
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

def load_stage(session_dir: Path, stage: int, name: str) -> dict:
    """讀取本引擎產生的 Stage（04/05），schema 必須與 SCHEMA_VERSION 一致。"""
    path = session_dir / f"stage_{stage:02d}_{name}.json"
    if not path.exists():
        raise FileNotFoundError(f"找不到 Stage {stage} 暫存：{path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("_schema_version") != SCHEMA_VERSION:
        raise ValueError(f"Schema version mismatch at stage {stage}，請重跑此關卡")
    return data


def load_stage01_from_module1(session_dir: Path) -> dict:
    """
    讀取 molgate_module1 產生的 stage_01_smiles.json。
    不與 Engine B 的 SCHEMA_VERSION 綁定（Module 1 獨立演進 0.2、0.3…）。
    """
    path = session_dir / "stage_01_smiles.json"
    if not path.exists():
        raise FileNotFoundError(f"找不到 Module 1 輸出：{path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    ver = data.get("_schema_version")
    if ver is not None and ver not in MODULE1_STAGE01_SCHEMAS:
        print(
            f"[WARN] Stage 01 schema 為 {ver!r}，未在已知清單 {sorted(MODULE1_STAGE01_SCHEMAS)}，仍嘗試繼續"
        )
    if "canonical_smiles" not in data:
        raise ValueError("stage_01_smiles.json 缺少 canonical_smiles")
    return data

def prompt_select(candidates: list[str], label: str) -> str:
    """人工選擇介面"""
    print(f"\n  ── 需要人工選擇：{label} ──")
    for i, c in enumerate(candidates):
        print(f"  [{i}] {c}")
    while True:
        try:
            idx = int(input("  請輸入編號選擇："))
            if 0 <= idx < len(candidates):
                return candidates[idx]
        except (ValueError, KeyboardInterrupt):
            pass
        print("  請輸入有效編號")

# ── Stage 04：質子化狀態 ──────────────────────────────
def stage04_protonation(smiles: str, target: str, session_dir: Path) -> dict:
    print("\n[Stage 04] 質子化狀態（pH 7.4）...")

    risk = LIGAND_RISK.get(target.upper(), {})
    tier = risk.get("stage04_tier", "B")
    rule = risk.get("stage04_rule", "dimorphite")
    note = risk.get("stage04_note", "")

    print(f"  → Tier: {tier} | Rule: {rule}")
    if note:
        print(f"  → Note: {note}")

    # Dimorphite-DL 產生候選
    candidates = []

    if DIMORPHITE_OK:
        try:
            from dimorphite_dl import protonate_smiles
            results = protonate_smiles(smiles, ph_min=7.35, ph_max=7.45)
            candidates = list(results) if results else [smiles]
        except Exception as e:
            print(f"  ⚠ Dimorphite-DL 執行錯誤：{e}，降級為原始 SMILES")
            candidates = [smiles]
    else:
        # 降級：直接用輸入 SMILES
        candidates = [smiles]

    print(f"  ✓ 候選數：{len(candidates)}")
    for c in candidates:
        print(f"    {c}")

    # 選擇邏輯
    if tier == "C":
        # 強制人工選擇
        selected = prompt_select(candidates, "質子化形態")
        status = "MANUAL"
    elif len(candidates) > 2:
        # 候選太多，自動升級到 C
        print(f"  ⚠ 候選數 > 2，自動升級為手動選擇")
        selected = prompt_select(candidates, "質子化形態")
        status = "MANUAL"
        tier = "C_auto"
    else:
        # 自動取第一個（主要形態）
        selected = candidates[0]
        status = "AUTO" if tier == "A" else "CONFIRM"

    print(f"  ✓ 選定：{selected}")
    print(f"  → Status: {status}")

    result = {
        "status": status,
        "tier": tier,
        "input_smiles": smiles,
        "candidates": candidates,
        "selected_smiles": selected,
        "ph": 7.4,
        "note": note,
    }
    return save_stage(session_dir, 4, "protonation", result)

# ── Stage 05：互變異構體 ──────────────────────────────
def stage05_tautomer(smiles: str, target: str, session_dir: Path) -> dict:
    print("\n[Stage 05] 互變異構體（canonical tautomer）...")

    risk = LIGAND_RISK.get(target.upper(), {})
    tier = risk.get("stage05_tier", "A")
    note = risk.get("stage05_note", "")

    print(f"  → Tier: {tier}")
    if note:
        print(f"  → Note: {note}")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result = {"status": "FAIL", "reason": "無法解析 Stage 04 輸出的 SMILES"}
        save_stage(session_dir, 5, "tautomer", result)
        raise ValueError(f"Stage 05 FAIL: {result['reason']}")

    # RDKit canonical tautomer
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(mol)
    canon_smiles = Chem.MolToSmiles(canon)

    # 列出所有候選
    tauts = enumerator.Enumerate(mol)
    all_tauts = [Chem.MolToSmiles(t) for t in tauts]

    print(f"  ✓ Canonical tautomer：{canon_smiles}")
    print(f"  ✓ 所有候選數：{len(all_tauts)}")

    # tier C 才停下來手選
    if tier == "C":
        print(f"  ⚠ 高風險互變標記，需手動確認")
        for t in all_tauts:
            print(f"    {t}")
        selected = prompt_select(all_tauts, "互變異構體")
        status = "MANUAL"
    else:
        selected = canon_smiles
        status = "AUTO" if tier == "A" else "CONFIRM"

    print(f"  ✓ 選定：{selected}")
    print(f"  → Status: {status}")

    result = {
        "status": status,
        "tier": tier,
        "input_smiles": smiles,
        "canonical_tautomer": canon_smiles,
        "all_tautomers": all_tauts,
        "selected_smiles": selected,
        "tautomer_count": len(all_tauts),
        "note": note,
    }
    return save_stage(session_dir, 5, "tautomer", result)

# ── 主流程 ────────────────────────────────────────────
def run_engine_b(session_dir: Path, target: str):
    print("=" * 55)
    print("  MolGate — Module 2 Engine B")
    print("  Stage 04: Protonation | Stage 05: Tautomer")
    print("=" * 55)

    try:
        try_load_manifest_or_rebuild(session_dir)
        validate_target_matches_manifest(session_dir, target)
    except ValueError as e:
        print(f"❌ {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    # 讀取 Stage 01 輸出（Module 1，schema 與本引擎分離）
    try:
        s1 = load_stage01_from_module1(session_dir)
    except FileNotFoundError:
        print("❌ 找不到 Module 1 輸出，請先跑 molgate_module1.py")
        sys.exit(1)

    smiles = s1["canonical_smiles"]
    print(f"\n  輸入 SMILES（來自 Stage 01）：{smiles}")
    print(f"  Target：{target}")

    if s1.get("recommend_resolve_stereo_before_module2"):
        print(
            "\n  ⚠ WARN：Stage 01 標記「建議先處理立體再進 Module 2」"
            f"（未定義手性 {s1.get('n_undefined_chiral_centers', '?')} | "
            f"未定義 E/Z {s1.get('n_undefined_ez_bonds', '?')}）"
        )
        print("     仍可繼續 Engine B；結果之解讀風險自負。\n")

    try:
        r4 = stage04_protonation(smiles, target, session_dir)
        r5 = stage05_tautomer(r4["selected_smiles"], target, session_dir)
    except ValueError as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)

    finalize_module2_engine_b_manifest(session_dir, r4, r5)
    append_registry(
        session_dir.parent,
        {
            "event": "module2_engine_b_complete",
            "session_id": session_dir.name,
            "session_dir": str(session_dir.resolve()),
            "target": target,
            "final_smiles": r5["selected_smiles"],
            "timestamp": datetime.now().isoformat(),
        },
    )

    print("\n" + "=" * 55)
    print("  Engine B 完成 ✓")
    print(f"  最終 SMILES → Stage 09 使用：{r5['selected_smiles']}")
    print(f"  Session log：{session_dir}")
    print(f"  Manifest：{session_dir / 'molgate_run_manifest.json'}")
    print("=" * 55)

    return {
        "stage04": r4,
        "stage05": r5,
        "final_smiles": r5["selected_smiles"],
    }

# ── CLI ──────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MolGate Module 2 Engine B")
    parser.add_argument("--session-dir", required=True, help="Module 1 的 session 目錄路徑")
    parser.add_argument("--target", required=True, help="Target 名稱（如 COX2）")
    args = parser.parse_args()

    session_dir = Path(args.session_dir)
    if not session_dir.exists():
        print(f"❌ Session 目錄不存在：{session_dir}")
        sys.exit(1)

    run_engine_b(session_dir, args.target)
