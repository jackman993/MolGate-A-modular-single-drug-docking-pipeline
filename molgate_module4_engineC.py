"""
MolGate — Module 4: Pocket Definition (Engine C)
Stage 23: Pocket center calculation

輸入：stage_22_ligand.json + protein/cleaned.pdb + Master Index（可選 box_padding）
輸出：stage_23_pocket_center.json、pocket_config.json（Module 5 / Vina 契約）

依賴：biopython、numpy；Index 載入與 Engine A/B 對齊。
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from Bio.PDB import PDBParser

_PKG_ROOT = Path(__file__).resolve().parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

from molgate_module4_engineA import load_target_index_entry
from molgate_module4_engineB import load_prep_overlay

SCHEMA_VERSION = "0.2"
STAGE22_SCHEMA_VERSION = "0.2"
DEFAULT_PADDING = 10.0  # Å
GPCR_DEFAULT_PADDING = 12.0


def resolve_box_padding(entry: dict[str, Any]) -> tuple[float, str]:
    """Index 條目或 override.box_padding；GPCR family 預設較大 padding。"""
    ov = entry.get("override") if isinstance(entry.get("override"), dict) else {}
    for label, src in (("Master Index", entry), ("override", ov)):
        v = src.get("box_padding")
        if v is not None:
            return float(v), f"{label} box_padding"
    fam = str(entry.get("family") or "").lower()
    if "gpcr" in fam:
        return GPCR_DEFAULT_PADDING, "family 含 GPCR → 預設 12 Å"
    return DEFAULT_PADDING, "預設 10 Å"


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


def load_stage22(session_dir: Path) -> dict[str, Any]:
    path = session_dir / "stage_22_ligand.json"
    if not path.exists():
        raise FileNotFoundError("找不到 stage_22_ligand.json，請先跑 Engine B")
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("_schema_version") != STAGE22_SCHEMA_VERSION:
        raise ValueError(
            f"stage_22_ligand.json schema 須為 {STAGE22_SCHEMA_VERSION}，請重跑 Module 4 Engine B"
        )
    if not data.get("ready_for_engine_c"):
        raise ValueError("stage_22_ligand.json 未標記 ready_for_engine_c=True，請重跑 Engine B")
    return data


def load_cleaned_pdb(session_dir: Path) -> Path:
    p = session_dir / "protein" / "cleaned.pdb"
    if not p.is_file():
        raise FileNotFoundError("找不到 protein/cleaned.pdb，請先跑 Module 3 Engine A")
    return p


def get_ligand_coords(
    structure,
    het_id: str,
    chain_id: str,
    resseq: int,
) -> np.ndarray:
    """從 PDB 結構抓取指定配體所有重原子座標。"""
    coords = []
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for res in chain:
                rn = res.get_resname().strip().upper()
                rs = res.get_id()[1]
                if rn == het_id.upper() and rs == resseq:
                    for atom in res.get_atoms():
                        if atom.element and atom.element.upper() != "H":
                            coords.append(atom.get_vector().get_array())
    if not coords:
        raise ValueError(
            f"在 cleaned.pdb 找不到配體：{het_id} chain {chain_id} res {resseq}\n"
            f"  → 請確認 cleaned.pdb 是否保留了 HETATM，"
            f"或 Stage 22 / Index het_id 與 resseq 是否正確"
        )
    return np.array(coords)


def prompt_confirm(label: str, value: Any) -> bool:
    print(f"\n  ── CONFIRM：{label} ──")
    if isinstance(value, dict):
        for k, v in value.items():
            if not k.startswith("_"):
                print(f"    {k}: {v}")
    else:
        print(f"    {value}")
    while True:
        try:
            raw = input("  接受此結果？[Y/n]：").strip().lower()
            if raw in ("", "y", "yes"):
                return True
            if raw in ("n", "no"):
                return False
        except KeyboardInterrupt:
            sys.exit(0)
        print("  請輸入 y 或 n")


def _update_manifest_module4_engine_c(session_dir: Path, pocket_config_path: Path) -> None:
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
        rmap["pocket_config"] = rel_path(session_dir, pocket_config_path)
        paths["relative"] = rmap
        paths["pocket_config_absolute"] = str(pocket_config_path.resolve())
        ns = dict(paths.get("next_step") or {})
        ns["after_module4_engine_c"] = {
            "read_pocket_config_from": rmap["pocket_config"],
            "hint": "Module 5：AutoDock Vina 使用 center / box_size",
        }
        paths["next_step"] = ns
        m["paths"] = paths
        pl = dict(m.get("pipeline") or {})
        pl["module4_engine_c"] = {
            "status": "complete",
            "completed_at": datetime.now().isoformat(),
            "stages": ["23"],
        }
        m["pipeline"] = pl
        write_manifest(session_dir, m)
    except Exception:
        pass


# ── Stage 23：口袋中心計算 ────────────────────────────
def stage23_pocket_center(session_dir: Path, entry: dict[str, Any]) -> dict[str, Any]:
    print("\n[Stage 23] 口袋中心計算...")

    s22 = load_stage22(session_dir)
    het_id = s22["selected_het_id"]
    chain = s22["chain"]
    resseq = s22["resseq"]
    pocket_rule = (s22.get("pocket_rule") or "").strip()

    padding, padding_src = resolve_box_padding(entry)

    print(f"  → 配體：{het_id} chain {chain} res {resseq}")
    print(f"  → Pocket rule：{pocket_rule}")
    print(f"  → Padding：{padding} Å（{padding_src}）")

    cleaned_path = load_cleaned_pdb(session_dir)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(cleaned_path))

    try:
        coords = get_ligand_coords(structure, het_id, chain, resseq)
    except ValueError as e:
        result = {
            "status": "FAIL",
            "reason": str(e),
        }
        save_stage(session_dir, 23, "pocket_center", result)
        raise

    atom_count = len(coords)
    print(f"  ✓ 配體重原子數：{atom_count}")

    center = coords.mean(axis=0)
    cx, cy, cz = float(center[0]), float(center[1]), float(center[2])

    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    span = max_coords - min_coords

    box_x = round(float(span[0]) + padding * 2, 1)
    box_y = round(float(span[1]) + padding * 2, 1)
    box_z = round(float(span[2]) + padding * 2, 1)

    box_x = max(box_x, 20.0)
    box_y = max(box_y, 20.0)
    box_z = max(box_z, 20.0)

    print(f"  ✓ 重心：({cx:.3f}, {cy:.3f}, {cz:.3f})")
    print(f"  ✓ Box size：{box_x} × {box_y} × {box_z} Å")

    confirm_data = {
        "center_x": round(cx, 3),
        "center_y": round(cy, 3),
        "center_z": round(cz, 3),
        "box_size_x": box_x,
        "box_size_y": box_y,
        "box_size_z": box_z,
        "padding": padding,
    }

    accepted = prompt_confirm("口袋座標與 Box Size", confirm_data)

    if not accepted:
        print("\n  手動輸入座標（直接 Enter 保留計算值）：")

        def ask_float(label: str, default: float) -> float:
            try:
                raw = input(f"    {label} [{default}]：").strip()
                return float(raw) if raw else default
            except (ValueError, KeyboardInterrupt):
                return default

        cx = ask_float("center_x", round(cx, 3))
        cy = ask_float("center_y", round(cy, 3))
        cz = ask_float("center_z", round(cz, 3))
        box_x = ask_float("box_size_x", box_x)
        box_y = ask_float("box_size_y", box_y)
        box_z = ask_float("box_size_z", box_z)
        padding = ask_float("padding", padding)
        status = "MANUAL"
        print(f"  ✓ 手動座標接受")
    else:
        status = "CONFIRM"
        print(f"  ✓ 計算結果接受")

    print(f"  → Status: {status}")

    pdb_id = (s22.get("pdb_id") or "").strip().upper()
    result: dict[str, Any] = {
        "status": status,
        "het_id": het_id,
        "chain": chain,
        "resseq": resseq,
        "atom_count": atom_count,
        "center_x": round(cx, 3),
        "center_y": round(cy, 3),
        "center_z": round(cz, 3),
        "box_size_x": box_x,
        "box_size_y": box_y,
        "box_size_z": box_z,
        "padding": padding,
        "padding_source": padding_src,
        "pocket_rule": pocket_rule,
    }
    if pdb_id:
        result["pdb_id"] = pdb_id
    cn = (s22.get("common_name") or "").strip()
    td = (s22.get("target_drug") or "").strip()
    if cn:
        result["common_name"] = cn
    if td:
        result["target_drug"] = td
    return save_stage(session_dir, 23, "pocket_center", result)


def write_pocket_config(
    session_dir: Path,
    r23: dict[str, Any],
    *,
    target: str,
) -> Path:
    config: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "target": target.strip().upper(),
        "het_id": r23["het_id"],
        "chain": r23["chain"],
        "resseq": r23["resseq"],
        "center_x": r23["center_x"],
        "center_y": r23["center_y"],
        "center_z": r23["center_z"],
        "box_size_x": r23["box_size_x"],
        "box_size_y": r23["box_size_y"],
        "box_size_z": r23["box_size_z"],
        "padding": r23["padding"],
        "atom_count": r23["atom_count"],
        "pocket_rule": r23.get("pocket_rule", ""),
        "ready_for_docking": r23["status"] in ("CONFIRM", "MANUAL"),
        "_timestamp": datetime.now().isoformat(),
    }
    if r23.get("pdb_id"):
        config["pdb_id"] = r23["pdb_id"]
    if r23.get("common_name"):
        config["common_name"] = r23["common_name"]
    if r23.get("target_drug"):
        config["target_drug"] = r23["target_drug"]
    path = session_dir / "pocket_config.json"
    path.write_text(json.dumps(config, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def run_engine_c(
    session_dir: Path,
    target: str,
    *,
    index_path: Path | None = None,
    skip_manifest_target_check: bool = False,
    skip_prep_target_check: bool = False,
):
    print("=" * 55)
    print("  MolGate — Module 4 Engine C")
    print("  Stage 23: 口袋中心計算 → pocket_config.json")
    print("=" * 55)

    if not skip_manifest_target_check:
        try:
            from molgate_manifest import validate_target_matches_manifest

            validate_target_matches_manifest(session_dir, target)
        except FileNotFoundError:
            pass
        except ValueError as e:
            print(f"❌ {e}")
            sys.exit(1)

    prep = load_prep_overlay(session_dir, target)
    if not skip_prep_target_check:
        pt = (prep.get("target") or "").strip().upper()
        if pt and pt != target.strip().upper():
            print(
                f"❌ prep_decisions.target={prep.get('target')!r} 與 --target={target!r} 不一致"
            )
            sys.exit(1)

    entry = load_target_index_entry(session_dir, prep, index_path)
    if not entry:
        print(f"❌ 無法載入 Master Index 條目（--target 或 master_index.json）")
        sys.exit(1)

    print(f"\n  Target：{target.upper()}")

    try:
        r23 = stage23_pocket_center(session_dir, entry)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)

    if r23.get("status") == "FAIL":
        print(f"\n❌ Stage 23 失敗：{r23.get('reason', '')}")
        sys.exit(1)

    config_path = write_pocket_config(session_dir, r23, target=target)
    _update_manifest_module4_engine_c(session_dir, config_path)

    print("\n" + "=" * 55)
    print("  Engine C 完成 ✓")
    print(f"  中心：({r23['center_x']}, {r23['center_y']}, {r23['center_z']})")
    print(f"  Box：{r23['box_size_x']} × {r23['box_size_y']} × {r23['box_size_z']} Å")
    print(f"  契約：{config_path}")
    print(f"  ready_for_docking = {r23['status'] in ('CONFIRM', 'MANUAL')}")
    print("=" * 55)

    return r23


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="MolGate Module 4 Engine C")
    ap.add_argument("--session-dir", required=True, help="Session 目錄路徑")
    ap.add_argument("--target", required=True, help="Target 名稱（如 COX2）")
    ap.add_argument(
        "--index",
        default=None,
        help="master_index.json 路徑（預設：session_meta / manifest 或 MolGate 目錄）",
    )
    ap.add_argument(
        "--skip-manifest-target-check",
        action="store_true",
        help="不與 molgate_run_manifest.json 的 target 交叉驗證",
    )
    ap.add_argument(
        "--skip-prep-target-check",
        action="store_true",
        help="不檢查 prep_decisions.target 與 --target 一致",
    )
    args = ap.parse_args()

    sd = Path(args.session_dir)
    if not sd.is_dir():
        print(f"❌ Session 目錄不存在：{sd}")
        sys.exit(1)

    idx = Path(args.index) if args.index else None
    run_engine_c(
        sd,
        args.target,
        index_path=idx,
        skip_manifest_target_check=args.skip_manifest_target_check,
        skip_prep_target_check=args.skip_prep_target_check,
    )
