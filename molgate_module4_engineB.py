"""
MolGate — Module 4: Pocket Definition (Engine B)
Stage 22: Ligand identity confirmation (人機)

輸入：stage_21_hetatm.json + Master Index（het_id / cocrystal_het、pocket_rule、可選 common_name、target_drug）
輸出：stage_22_ligand.json（選定配體）
依賴：與 molgate_module4_engineA 共用 Index 解析（molgate_module3_engineA）
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

_PKG_ROOT = Path(__file__).resolve().parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

from molgate_module4_engineA import load_target_index_entry

SCHEMA_VERSION = "0.2"


def load_prep_overlay(session_dir: Path, target: str) -> dict[str, Any]:
    """合併 prep_decisions.json 與 CLI target（CLI 優先）。"""
    prep: dict[str, Any] = {"target": target.strip()}
    p = session_dir / "prep_decisions.json"
    if p.is_file():
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
            if isinstance(data, dict):
                for k, v in data.items():
                    if k.startswith("_"):
                        continue
                    prep[k] = v
        except (json.JSONDecodeError, OSError):
            pass
    prep["target"] = target.strip()
    return prep


def expected_het_id(entry: dict[str, Any]) -> str:
    """比對用單一代碼：優先 het_id，否則 cocrystal_het。"""
    for key in ("het_id", "cocrystal_het"):
        v = entry.get(key)
        if isinstance(v, str) and v.strip():
            return v.strip().upper()
    return ""


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


def load_stage21(session_dir: Path) -> dict[str, Any]:
    path = session_dir / "stage_21_hetatm.json"
    if not path.exists():
        raise FileNotFoundError("找不到 stage_21_hetatm.json，請先跑 Engine A")
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("_schema_version") != SCHEMA_VERSION:
        raise ValueError("Stage 21 schema version mismatch，請重跑 Engine A")
    return data


def prompt_select(candidates: list[dict], label: str) -> dict:
    print(f"\n  ── {label} ──")
    for i, c in enumerate(candidates):
        tag = "【目標】 " if c.get("ligand_role") == "target" else ""
        print(
            f"  [{i}] {tag}{c['resname']} chain {c['chain']} "
            f"res {c['resseq']} (~{c['mw_estimate']:.0f} Da, "
            f"{c['heavy_atom_count']} heavy atoms)"
        )
    while True:
        try:
            raw = input("  請輸入編號選擇：").strip()
            idx = int(raw)
            if 0 <= idx < len(candidates):
                return candidates[idx]
        except (ValueError, KeyboardInterrupt):
            pass
        print("  請輸入有效編號")


# ── Stage 22：配體身份確認 ────────────────────────────
def stage22_ligand(
    session_dir: Path,
    target: str,
    index_path: Path | None,
) -> dict[str, Any]:
    print("\n[Stage 22] 配體身份確認...")

    s21 = load_stage21(session_dir)
    candidates = s21.get("candidates", [])

    if not candidates:
        result = {
            "status": "FAIL",
            "reason": "Stage 21 無候選配體，無法確認配體身份",
        }
        save_stage(session_dir, 22, "ligand", result)
        raise ValueError(f"Stage 22 FAIL: {result['reason']}")

    prep = load_prep_overlay(session_dir, target)
    entry = load_target_index_entry(session_dir, prep, index_path)
    if not entry:
        result = {
            "status": "FAIL",
            "reason": "無法載入 Master Index 條目（請確認 --target 與 master_index.json）",
        }
        save_stage(session_dir, 22, "ligand", result)
        raise ValueError(f"Stage 22 FAIL: {result['reason']}")

    het_id = expected_het_id(entry)
    pocket_rule = (entry.get("pocket_rule") or "").strip()
    common_name = (entry.get("common_name") or "").strip()
    target_drug = (entry.get("target_drug") or "").strip()
    pdb_id = (entry.get("pdb_id") or prep.get("pdb_id") or "").strip().upper()

    print(f"  → Master Index expected het（het_id／cocrystal_het）：{het_id or '（未設定）'}")
    if common_name:
        print(f"  → 共晶配體俗名（common_name）：{common_name}")
    if target_drug:
        print(f"  → Docking 目標藥物（target_drug）：{target_drug}")
    print(f"  → Pocket rule: {pocket_rule}")
    print(f"  → 候選配體數：{len(candidates)}")

    if not het_id:
        print(f"  ⚠ 無法自動比對（Index 缺 het_id / cocrystal_het），請手動選擇")
        selected = prompt_select(candidates, "請選擇用於口袋定義的配體")
        status = "MANUAL"
        print(f"  ✓ 選定：{selected['resname']} chain {selected['chain']}")
    else:
        matched = [c for c in candidates if c["resname"].upper() == het_id]

        if len(matched) == 1:
            selected = matched[0]
            status = "AUTO"
            print(f"  ✓ 自動比對成功：{selected['resname']} chain {selected['chain']}")

        elif len(matched) > 1:
            print(f"  ⚠ 多個 {het_id} 存在於不同 chain，需要確認")
            selected = prompt_select(matched, f"選擇 {het_id} 所在 chain")
            status = "MANUAL"
            print(f"  ✓ 選定：{selected['resname']} chain {selected['chain']}")

        elif len(candidates) == 1:
            print(
                f"  ⚠ WARN: 唯一候選 {candidates[0]['resname']} "
                f"與 Index 期望 {het_id} 不符"
            )
            print(f"  → 可能是 PDB／CCD 代碼與 Index 不一致")
            confirm = input(
                f"  是否使用 {candidates[0]['resname']} 作為口袋定義配體？[y/N]："
            ).strip().lower()
            if confirm == "y":
                selected = candidates[0]
                status = "MANUAL"
                print(f"  ✓ 手動確認：{selected['resname']}")
            else:
                result = {
                    "status": "FAIL",
                    "reason": f"用戶拒絕使用 {candidates[0]['resname']}，無法繼續",
                }
                save_stage(session_dir, 22, "ligand", result)
                raise ValueError(f"Stage 22 FAIL: {result['reason']}")

        else:
            print(f"  ⚠ 無法自動比對，候選列表如下：")
            selected = prompt_select(candidates, "請選擇用於口袋定義的配體")
            status = "MANUAL"
            print(f"  ✓ 選定：{selected['resname']} chain {selected['chain']}")

    print(f"  → Status: {status}")

    result: dict[str, Any] = {
        "status": status,
        "selected_het_id": selected["resname"],
        "chain": selected["chain"],
        "resseq": selected["resseq"],
        "mw_estimate": selected["mw_estimate"],
        "heavy_atom_count": selected["heavy_atom_count"],
        "master_index_expected_het": het_id,
        "master_index_het_id": het_id,
        "match": selected["resname"].upper() == het_id if het_id else None,
        "pocket_rule": pocket_rule,
        "ready_for_engine_c": True,
    }
    if pdb_id:
        result["pdb_id"] = pdb_id
    if common_name:
        result["common_name"] = common_name
    if target_drug:
        result["target_drug"] = target_drug
    return save_stage(session_dir, 22, "ligand", result)


# ── 主流程 ────────────────────────────────────────────
def run_engine_b(session_dir: Path, target: str, index_path: Path | None = None):
    print("=" * 55)
    print("  MolGate — Module 4 Engine B")
    print("  Stage 22: 配體身份確認（人機）")
    print("=" * 55)

    prep = load_prep_overlay(session_dir, target)
    entry = load_target_index_entry(session_dir, prep, index_path)
    if not entry:
        print(f"❌ Target '{target}' 無法從 Master Index 載入條目（檔案或鍵名錯誤）")
        sys.exit(1)

    pdb_id = (entry.get("pdb_id") or prep.get("pdb_id") or "?").strip()
    het_hint = expected_het_id(entry) or "?"
    print(f"\n  Target：{target.upper()} | PDB：{pdb_id}")
    print(f"  Index 期望共晶代碼：{het_hint}")

    try:
        result = stage22_ligand(session_dir, target, index_path)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)

    print("\n" + "=" * 55)
    print("  Engine B 完成 ✓")
    print(
        f"  選定配體：{result['selected_het_id']} "
        f"chain {result['chain']} res {result['resseq']}"
    )
    print(f"  ready_for_engine_c = True")
    print(f"  輸出：{session_dir / 'stage_22_ligand.json'}")
    print("=" * 55)

    return result


# ── CLI ──────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="MolGate Module 4 Engine B")
    ap.add_argument("--session-dir", required=True, help="Session 目錄路徑")
    ap.add_argument("--target", required=True, help="Target 名稱（如 COX2，須對應 master_index.json）")
    ap.add_argument(
        "--index",
        default=None,
        help="master_index.json 路徑（預設：session_meta / manifest 或 MolGate 目錄下）",
    )
    args = ap.parse_args()

    sd = Path(args.session_dir)
    if not sd.is_dir():
        print(f"❌ Session 目錄不存在：{sd}")
        sys.exit(1)

    idx = Path(args.index) if args.index else None
    run_engine_b(sd, args.target, index_path=idx)
