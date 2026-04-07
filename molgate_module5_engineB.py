"""
MolGate — Module 5: Docking & Output (Engine B)
Stage 27: 最終人工審核

輸入：docking_result.json（Engine A 契約）
輸出：stage_27_review.json、final_report.json

決策：
  ACCEPT  → 寫入 final_report.json，pipeline 完成
  RETRY   → 標記需重跑的 Module，不寫最終報告
  FLAG    → 標記問題，結果存檔但不視為完成

依賴：無額外依賴
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

SCHEMA_VERSION = "0.2"

RETRY_TARGETS = {
    "1": ("Module 1", "SMILES / PDB 輸入問題"),
    "2": ("Module 2", "配體準備問題（質子化、手性）"),
    "3": ("Module 3", "蛋白準備問題（去水、殘基 override）"),
    "4": ("Module 4", "口袋定義問題（配體選擇、座標）"),
    "5": ("Module 5", "Docking 參數問題（exhaustiveness、box size）"),
}

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

def load_docking_result(session_dir: Path) -> dict[str, Any]:
    path = session_dir / "docking_result.json"
    if not path.exists():
        raise FileNotFoundError("找不到 docking_result.json，請先跑 Engine A")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not data.get("ready_for_review"):
        raise ValueError("docking_result.json 未標記 ready_for_review=True，請重跑 Engine A")
    sv = str(data.get("schema_version") or "").strip()
    if sv and sv != SCHEMA_VERSION:
        print(f"  ⚠ docking_result.schema_version={sv!r}，建議重跑 Module 5 Engine A（{SCHEMA_VERSION}）")
    return data

def display_result(dr: dict[str, Any]):
    """顯示 docking 結果摘要。"""
    print("\n" + "─" * 55)
    print("  DOCKING 結果摘要")
    print("─" * 55)

    target = dr.get("target", "?")
    pdb_id = dr.get("pdb_id", "?")
    het_id = dr.get("het_id", "?")
    chain  = dr.get("chain", "?")
    td     = dr.get("target_drug", "")
    cn     = dr.get("common_name", "")

    print(f"  Target：{target} | PDB：{pdb_id}")
    if td:
        print(f"  Docking 藥物：{td}")
    if cn:
        print(f"  共晶配體：{cn} ({het_id} chain {chain})")
    else:
        print(f"  共晶配體：{het_id} chain {chain}")

    print(f"\n  Best pose：Mode {dr['best_pose_mode']} | "
          f"{dr['best_affinity']:.3f} kcal/mol")
    print(f"  Confidence tier：{dr['confidence_tier']}")
    print(f"  QC pass：{dr['qc_pass']}")

    cutoff = float(dr.get("out_of_pocket_rmsd_lb_cutoff", 8.0))
    oop_modes = dr.get("out_of_pocket_modes") or []
    poses = dr.get("all_poses", [])
    if poses:
        print(f"\n  所有 pose：")
        for p in poses:
            rmsd_warn = " ⚠ 口袋外" if p["rmsd_lb"] > cutoff else ""
            print(f"    Mode {p['mode']:2d}: {p['affinity']:7.3f} kcal/mol | "
                  f"RMSD lb={p['rmsd_lb']:.2f}{rmsd_warn}")
    if oop_modes:
        print(f"  ⚠ OUT_OF_POCKET cutoff：{cutoff:.1f} Å | modes={oop_modes}")

    flags = dr.get("flags", [])
    if flags:
        print(f"\n  Flags：")
        for f in flags:
            print(f"    ⚠ {f}")

    print("─" * 55)

def prompt_decision() -> tuple[str, str]:
    """人工審核決策介面。"""
    print("\n  ── Stage 27：最終審核 ──")
    print("  [1] ACCEPT  → 接受結果，寫入最終報告")
    print("  [2] RETRY   → 標記問題，指定重跑 Module")
    print("  [3] FLAG    → 標記問題，結果存檔但不視為完成")
    print()

    while True:
        try:
            raw = input("  請選擇 [1/2/3]：").strip()
            if raw in ("1", "2", "3"):
                break
        except KeyboardInterrupt:
            sys.exit(0)
        print("  請輸入 1、2 或 3")

    note = ""
    if raw == "1":
        try:
            note = input("  備註（選填，直接 Enter 跳過）：").strip()
        except KeyboardInterrupt:
            pass
        return "ACCEPT", note

    if raw == "2":
        print("\n  指定需重跑的 Module：")
        for k, (mod, reason) in RETRY_TARGETS.items():
            print(f"  [{k}] {mod}：{reason}")
        retry_mod = ""
        retry_reason = ""
        while True:
            try:
                m = input("  請選擇 Module [1-5]：").strip()
                if m in RETRY_TARGETS:
                    retry_mod, retry_reason = RETRY_TARGETS[m]
                    break
            except KeyboardInterrupt:
                sys.exit(0)
            print("  請輸入 1～5")
        try:
            note = input(f"  問題說明（選填）：").strip()
        except KeyboardInterrupt:
            pass
        return "RETRY", f"重跑 {retry_mod}：{retry_reason}。{note}".strip()

    if raw == "3":
        try:
            note = input("  FLAG 原因（選填）：").strip()
        except KeyboardInterrupt:
            pass
        return "FLAG", note

    return "FLAG", ""


def build_non_interactive_note(
    decision: str,
    note: str | None,
    retry_module: str | None = None,
) -> str:
    """API/自動化模式下組裝 note。"""
    base = (note or "").strip()
    if decision == "RETRY":
        if retry_module and retry_module in RETRY_TARGETS:
            mod, reason = RETRY_TARGETS[retry_module]
            prefix = f"重跑 {mod}：{reason}"
            return f"{prefix}。{base}".strip()
        return base or "重跑（未指定模組）"
    return base

def write_final_report(
    session_dir: Path,
    dr: dict[str, Any],
    note: str,
) -> Path:
    """寫入最終報告。"""
    report = {
        "schema_version": SCHEMA_VERSION,
        "status": "COMPLETE",
        "target": dr.get("target", ""),
        "pdb_id": dr.get("pdb_id", ""),
        "het_id": dr.get("het_id", ""),
        "chain":  dr.get("chain", ""),
        "pocket_rule": dr.get("pocket_rule", ""),
        "best_affinity": dr["best_affinity"],
        "best_pose_mode": dr["best_pose_mode"],
        "confidence_tier": dr["confidence_tier"],
        "qc_pass": dr["qc_pass"],
        "flags": dr.get("flags", []),
        "all_poses": dr.get("all_poses", []),
        "output_pdbqt": dr.get("output_pdbqt", ""),
        "pocket_center": dr.get("pocket_center", {}),
        "out_of_pocket_rmsd_lb_cutoff": dr.get("out_of_pocket_rmsd_lb_cutoff"),
        "out_of_pocket_modes": dr.get("out_of_pocket_modes", []),
        "reviewer_note": note,
        "reviewed_at": datetime.now().isoformat(),
    }
    for key in ("common_name", "target_drug"):
        if dr.get(key):
            report[key] = dr[key]

    path = session_dir / "final_report.json"
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return path

def _update_manifest_module5_engine_b(
    session_dir: Path,
    decision: str,
    report_path: Path | None,
) -> None:
    try:
        from molgate_manifest import read_manifest, rel_path, write_manifest
    except ImportError:
        return
    mp = session_dir / "molgate_run_manifest.json"
    if not mp.is_file():
        return
    try:
        m = read_manifest(session_dir)
        pl = dict(m.get("pipeline") or {})
        pl["module5_engine_b"] = {
            "status": decision.lower(),
            "completed_at": datetime.now().isoformat(),
            "stages": ["27"],
        }
        m["pipeline"] = pl
        if report_path:
            paths = dict(m.get("paths") or {})
            rmap = dict(paths.get("relative") or {})
            rmap["final_report"] = rel_path(session_dir, report_path)
            paths["relative"] = rmap
            m["paths"] = paths
            summ = dict(m.get("summary") or {})
            summ["pipeline_complete"] = True
            summ["final_report"] = rmap["final_report"]
            m["summary"] = summ
        write_manifest(session_dir, m)
    except Exception:
        pass

# ── Stage 27：最終人工審核 ────────────────────────────
def stage27_review(
    session_dir: Path,
    dr: dict[str, Any],
    *,
    decision_override: str | None = None,
    note_override: str | None = None,
    retry_module: str | None = None,
) -> dict[str, Any]:
    print("\n[Stage 27] 最終人工審核...")

    display_result(dr)
    if decision_override:
        decision = decision_override.strip().upper()
        if decision not in ("ACCEPT", "RETRY", "FLAG"):
            raise ValueError(f"不支援 decision={decision_override!r}，需為 ACCEPT/RETRY/FLAG")
        note = build_non_interactive_note(decision, note_override, retry_module)
    else:
        decision, note = prompt_decision()

    report_path = None
    status = decision

    if decision == "ACCEPT":
        report_path = write_final_report(session_dir, dr, note)
        print(f"\n  ✓ 結果接受")
        print(f"  ✓ 最終報告：{report_path}")
    elif decision == "RETRY":
        print(f"\n  ↩ 標記重跑：{note}")
    elif decision == "FLAG":
        print(f"\n  ⚑ 標記問題：{note or '（無說明）'}")

    print(f"  → Status: {status}")

    result = {
        "status": status,
        "decision": decision,
        "note": note,
        "final_report": str(report_path.as_posix()) if report_path else None,
        "pipeline_complete": decision == "ACCEPT",
    }
    r = save_stage(session_dir, 27, "review", result)
    _update_manifest_module5_engine_b(session_dir, decision, report_path)
    return r

# ── 主流程 ────────────────────────────────────────────
def run_engine_b(
    session_dir: Path,
    *,
    decision: str | None = None,
    note: str | None = None,
    retry_module: str | None = None,
):
    print("=" * 55)
    print("  MolGate — Module 5 Engine B")
    print("  Stage 27：最終人工審核")
    print("=" * 55)

    try:
        dr = load_docking_result(session_dir)
    except (FileNotFoundError, ValueError) as e:
        print(f"❌ {e}")
        sys.exit(1)

    try:
        result = stage27_review(
            session_dir,
            dr,
            decision_override=decision,
            note_override=note,
            retry_module=retry_module,
        )
    except Exception as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)

    print("\n" + "=" * 55)
    decision = result["decision"]
    if decision == "ACCEPT":
        print("  Pipeline 完成 ✓")
        print(f"  最終報告：{session_dir / 'final_report.json'}")
    elif decision == "RETRY":
        print(f"  ↩ 需重跑：{result['note']}")
    else:
        print(f"  ⚑ 已標記：{result['note'] or '（無說明）'}")
    print("=" * 55)

    return result

# ── CLI ──────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="MolGate Module 5 Engine B")
    ap.add_argument("--session-dir", required=True)
    ap.add_argument(
        "--decision",
        choices=("ACCEPT", "RETRY", "FLAG"),
        default=None,
        help="非互動模式決策（供 API/UI 使用）",
    )
    ap.add_argument(
        "--note",
        default=None,
        help="非互動模式備註（可選）",
    )
    ap.add_argument(
        "--retry-module",
        choices=("1", "2", "3", "4", "5"),
        default=None,
        help="decision=RETRY 時指定重跑模組（1-5）",
    )
    args = ap.parse_args()

    sd = Path(args.session_dir)
    if not sd.is_dir():
        print(f"❌ Session 目錄不存在：{sd}")
        sys.exit(1)

    run_engine_b(
        sd,
        decision=args.decision,
        note=args.note,
        retry_module=args.retry_module,
    )
