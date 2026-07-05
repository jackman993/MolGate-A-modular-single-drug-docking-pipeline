"""
MolGate — Module 3: Protein Preparation (Engine B)
人機決策層（單一腳本跑完）：
  Stage 15：有偵測到輔因子／金屬才停，否則一路向下
  Stage 17：僅 GPCR／金屬酶相關標靶才停，否則一路向下
  Flag 16／18：PROPKA／鹽橋占位（SKIPPED）
最後寫入 prep_decisions.json，並設 ready_for_engine_e = true（供 Engine C 銜接）

依賴：biopython；Index 載入與 Engine A 共用（import molgate_module3_engineA）
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

from molgate_module3_engineA import build_engine_entry, load_master_index_file

SCHEMA_VERSION = "0.3"
PREP_DECISIONS_SCHEMA = "0.3"

CRYSTALLIZATION_HET: set[str] = {
    "HOH", "WAT", "GOL", "EDO", "PEG", "MPD",
    "SO4", "PO4", "ACT", "DMS", "FMT", "TRS", "MES", "EPE",
    "CL", "NA", "K", "BR", "IOD", "NO3", "NH4",
}

METAL_RESNAMES: set[str] = {
    "ZN", "MG", "MN", "FE", "FE2", "FES", "CU", "NI", "CO", "MO", "W", "V",
    "CD", "HG", "PT", "AU", "AG", "PB", "CA", "SR", "BA", "RB", "CS", "LI",
    "CR", "EU", "GD", "SM", "IR", "RU", "RH", "OS", "RE", "TA",
}


def family_needs_stage17_hitl(entry: dict[str, Any]) -> bool:
    """GPCR 或金屬酶／含金屬輔因子語意之 family 才需 Stage 17 人機確認。"""
    if entry.get("metalloenzyme") is True:
        return True
    fam_raw = (entry.get("family") or "").strip()
    if not fam_raw:
        return False
    fup = fam_raw.upper()
    if fup == "GPCR" or "GPCR" in fup:
        return True
    flo = fam_raw.lower()
    for key in ("metall", "metallo", "zinc", "metal"):
        if key in flo:
            return True
    return False


def save_stage(session_dir: Path, stage: int, name: str, data: dict) -> dict:
    data = dict(data)
    data["_stage"] = stage
    data["_name"] = name
    data["_timestamp"] = datetime.now().isoformat()
    data["_schema_version"] = SCHEMA_VERSION
    path = session_dir / f"stage_{stage:02d}_{name}.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return data


def load_cleaned_pdb(session_dir: Path) -> Path:
    path = session_dir / "protein" / "cleaned.pdb"
    if not path.is_file():
        raise FileNotFoundError("找不到 cleaned.pdb，請先跑 Engine A")
    return path


def infer_kept_chains(session_dir: Path, structure) -> list[str]:
    s12 = session_dir / "stage_12_chains.json"
    if s12.is_file():
        try:
            data = json.loads(s12.read_text(encoding="utf-8"))
            kc = data.get("kept_chains")
            if isinstance(kc, list) and kc:
                return kc
        except (json.JSONDecodeError, OSError, KeyError):
            pass
    return sorted({c.id for m in structure for c in m})


def load_prep_decisions(session_dir: Path) -> dict[str, Any]:
    path = session_dir / "prep_decisions.json"
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    return {
        "schema_version": PREP_DECISIONS_SCHEMA,
        "decisions": {},
        "flags": [],
    }


def save_prep_decisions(session_dir: Path, decisions: dict[str, Any]) -> None:
    path = session_dir / "prep_decisions.json"
    decisions = dict(decisions)
    decisions["_timestamp"] = datetime.now().isoformat()
    path.write_text(json.dumps(decisions, indent=2, ensure_ascii=False), encoding="utf-8")


def prompt_choice(options: list[str], label: str, default: int = 0) -> int:
    print(f"\n  ── {label} ──")
    for i, opt in enumerate(options):
        marker = "（預設）" if i == default else ""
        print(f"  [{i}] {opt}{marker}")
    while True:
        try:
            raw = input(f"  請輸入編號（直接 Enter = {default}）：").strip()
            if raw == "":
                return default
            idx = int(raw)
            if 0 <= idx < len(options):
                return idx
        except (ValueError, KeyboardInterrupt):
            pass
        print("  請輸入有效編號")


def prompt_yes_no(question: str, default: bool = True) -> bool:
    default_str = "Y/n" if default else "y/N"
    while True:
        try:
            raw = input(f"  {question} [{default_str}]：").strip().lower()
            if raw == "":
                return default
            if raw in ("y", "yes"):
                return True
            if raw in ("n", "no"):
                return False
        except KeyboardInterrupt:
            sys.exit(0)
        print("  請輸入 y 或 n")


def _cocrystal_resnames(entry: dict[str, Any]) -> set[str]:
    names: set[str] = set()
    for key in ("cocrystal_het", "het_id"):
        v = entry.get(key)
        if isinstance(v, str) and v.strip():
            names.add(v.strip().upper())
    return names


def _normalize_override(raw: Any) -> dict[str, str]:
    if not isinstance(raw, dict):
        return {}
    out: dict[str, str] = {}
    for k, v in raw.items():
        if isinstance(k, str) and v is not None:
            out[k.strip()] = str(v).strip()
    return out


def stage15_cofactor(
    entry: dict[str, Any],
    structure,
    kept_chains: list[str],
    session_dir: Path,
    *,
    non_interactive: bool,
    auto_action: str,
) -> dict[str, Any]:
    print("\n[Stage 15] 輔因子／金屬處理...")

    cocrystal = _cocrystal_resnames(entry)
    cofactor_index = entry.get("cofactor")
    metal_index = entry.get("metal")

    detected: list[dict[str, Any]] = []
    seen: set[tuple[str, str, Any]] = set()

    for model in structure:
        for chain in model:
            if chain.id not in kept_chains:
                continue
            for res in chain:
                resname = res.get_resname().strip().upper()
                hetflag = res.get_id()[0].strip()
                if not (hetflag == "H_" or hetflag.startswith("H")):
                    continue
                if resname in cocrystal:
                    continue
                if resname in METAL_RESNAMES:
                    key = (resname, chain.id, res.get_id()[1])
                    if key not in seen:
                        seen.add(key)
                        detected.append({
                            "resname": resname,
                            "chain": chain.id,
                            "resseq": res.get_id()[1],
                            "kind": "metal",
                        })
                    continue
                if resname in CRYSTALLIZATION_HET:
                    continue
                key = (resname, chain.id, res.get_id()[1])
                if key not in seen:
                    seen.add(key)
                    detected.append({
                        "resname": resname,
                        "chain": chain.id,
                        "resseq": res.get_id()[1],
                        "kind": "het",
                    })

    index_flagged: list[str] = []
    if cofactor_index:
        index_flagged.append(f"cofactor:{cofactor_index}")
    if metal_index:
        index_flagged.append(f"metal:{metal_index}")

    needs_confirm = len(detected) > 0 or len(index_flagged) > 0

    if not needs_confirm:
        print("  ✓ 無需確認的輔因子／金屬候選，一路向下（AUTO）")
        print("  → Status: AUTO")
        return save_stage(session_dir, 15, "cofactor", {
            "status": "AUTO",
            "detected": [],
            "index_flagged": [],
            "action": "none",
            "kept": [],
            "removed": [],
        })

    print("  ⚠ 偵測到輔因子／金屬候選（人機確認）：")
    for d in detected:
        print(f"    [{d.get('kind', 'het')}] {d['resname']} chain {d['chain']} res {d['resseq']}")
    for f in index_flagged:
        print(f"    Master Index 標記：{f}")

    if non_interactive:
        action = auto_action if auto_action in ("remove_all", "keep_all") else "remove_all"
        if action == "remove_all":
            kept: list[str] = []
            removed = [d["resname"] for d in detected] + index_flagged
            print(f"  ✓ --non-interactive：{action}")
        else:
            kept = [d["resname"] for d in detected] + index_flagged
            removed = []
            print(f"  ✓ --non-interactive：{action}")
        status = "AUTO"
        return save_stage(session_dir, 15, "cofactor", {
            "status": status,
            "detected": detected,
            "index_flagged": index_flagged,
            "action": action,
            "kept": kept,
            "removed": removed,
        })

    idx = prompt_choice(
        ["移除全部（預設安全）", "全部保留", "逐一確認"],
        label="輔因子／金屬處理方式",
        default=0,
    )

    if idx == 0:
        action = "remove_all"
        kept = []
        removed = [d["resname"] for d in detected] + index_flagged
        print("  ✓ 選擇：移除全部")
        status = "MANUAL"
    elif idx == 1:
        action = "keep_all"
        kept = [d["resname"] for d in detected] + index_flagged
        removed = []
        print("  ✓ 選擇：保留全部")
        status = "MANUAL"
    else:
        kept = []
        removed = []
        for d in detected:
            keep = prompt_yes_no(
                f"保留 {d['resname']} chain {d['chain']}？",
                default=False,
            )
            if keep:
                kept.append(d["resname"])
            else:
                removed.append(d["resname"])
        action = "selective"
        status = "MANUAL"

    print(f"  → Status: {status}")
    return save_stage(session_dir, 15, "cofactor", {
        "status": status,
        "detected": detected,
        "index_flagged": index_flagged,
        "action": action,
        "kept": kept,
        "removed": removed,
    })


def flag16_propka(session_dir: Path) -> dict[str, Any]:
    print("\n[Flag 16] PROPKA 滴定態（占位，併入 Engine B）...")
    print("  → Flag: PROPKA_SKIPPED")
    print("  → Status: SKIPPED")
    return save_stage(session_dir, 16, "propka", {
        "status": "SKIPPED",
        "flag": "PROPKA_SKIPPED",
        "note": "快速版占位。正式版需 PROPKA/PDB2PQR 並重跑。",
    })


def stage17_override(
    entry: dict[str, Any],
    session_dir: Path,
    *,
    non_interactive: bool,
) -> dict[str, Any]:
    print("\n[Stage 17] 關鍵殘基 Override...")

    family = entry.get("family", "default")
    override = _normalize_override(entry.get("override"))
    needs_hitl = family_needs_stage17_hitl(entry) or len(override) > 0

    if not needs_hitl:
        print(f"  ✓ family={family}，非 GPCR／金屬酶優先類別且無 Index override，一路向下（AUTO）")
        print("  → Status: AUTO")
        return save_stage(session_dir, 17, "override", {
            "status": "AUTO",
            "family": family,
            "override_rules": {},
        })

    print(f"  ⚠ 需確認：family={family}（GPCR／金屬酶相關或具 Index override）")

    if override:
        print(f"  Master Index override 建議（{len(override)} 項）：")
        for residue, rule in override.items():
            print(f"    {residue} → {rule}")
    else:
        print("  Master Index 無預設 override，請自行確認質子化／關鍵殘基敘述")

    if non_interactive:
        final_rules = dict(override)
        print("  ✓ --non-interactive：採用 Index 之 override（無則空）")
        return save_stage(session_dir, 17, "override", {
            "status": "AUTO",
            "family": family,
            "override_rules": final_rules,
        })

    final_rules: dict[str, str] = {}
    status = "MANUAL"

    if override:
        accepted = prompt_yes_no("接受上述 override 建議？", default=True)
        if accepted:
            final_rules = dict(override)
            print("  ✓ 接受建議 override")
        else:
            print("  手動輸入 override（格式：殘基鍵=狀態，例 HIS114=protonated）")
            print("  輸入 done 或空白行結束")
            while True:
                try:
                    raw = input("  輸入：").strip()
                    if raw.lower() == "done" or raw == "":
                        break
                    if "=" in raw:
                        k, v = raw.split("=", 1)
                        final_rules[k.strip()] = v.strip()
                        print(f"  ✓ 加入：{k.strip()} → {v.strip()}")
                except KeyboardInterrupt:
                    break
    else:
        print("  手動輸入 override（殘基鍵=狀態）；直接 Enter 結束＝無 override")
        while True:
            try:
                raw = input("  輸入：").strip()
                if raw.lower() == "done" or raw == "":
                    break
                if "=" in raw:
                    k, v = raw.split("=", 1)
                    final_rules[k.strip()] = v.strip()
                    print(f"  ✓ 加入：{k.strip()} → {v.strip()}")
            except KeyboardInterrupt:
                break

    print(f"  → Status: {status}")
    return save_stage(session_dir, 17, "override", {
        "status": status,
        "family": family,
        "override_rules": final_rules,
    })


def flag18_saltbridge(session_dir: Path) -> dict[str, Any]:
    print("\n[Flag 18] 鹽橋分析（占位，併入 Engine B）...")
    print("  → Flag: SALT_BRIDGE_SKIPPED")
    print("  → Status: SKIPPED")
    return save_stage(session_dir, 18, "saltbridge", {
        "status": "SKIPPED",
        "flag": "SALT_BRIDGE_SKIPPED",
        "note": "快速版占位。正式版需 PROPKA 完成後再分析。",
    })


def _merge_flags(existing: list[Any], new_flags: list[str]) -> list[str]:
    out: list[str] = []
    for x in existing or []:
        if isinstance(x, str) and x not in out:
            out.append(x)
    for f in new_flags:
        if f not in out:
            out.append(f)
    return out


def _update_manifest_engine_b(session_dir: Path, prep_path: Path) -> None:
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
        rmap["prep_decisions"] = rel_path(session_dir, prep_path)
        paths["relative"] = rmap
        paths["prep_decisions_absolute"] = str(prep_path.resolve())
        ns = dict(paths.get("next_step") or {})
        ns["after_module3_engine_b"] = {
            "read_prep_decisions_from": rmap["prep_decisions"],
            "hint": "下一步：Engine C（receptor.pdbqt，Stage 19/20）；需 ready_for_engine_e=true",
        }
        paths["next_step"] = ns
        m["paths"] = paths
        pl = dict(m.get("pipeline") or {})
        pl["module3_engine_b"] = {
            "status": "complete",
            "completed_at": datetime.now().isoformat(),
            "stages": ["15", "16", "17", "18"],
        }
        m["pipeline"] = pl
        summ = dict(m.get("summary") or {})
        summ["prep_decisions_relative"] = rmap["prep_decisions"]
        summ["ready_for_engine_e"] = True
        m["summary"] = summ
        write_manifest(session_dir, m)
    except Exception:
        pass


def run_engine_b(
    session_dir: Path,
    target: str,
    *,
    index_path: Path | None = None,
    skip_manifest_target_check: bool = False,
    non_interactive: bool = False,
    auto_hets: str = "remove_all",
) -> dict[str, Any]:
    print("=" * 55)
    print("  MolGate — Module 3 Engine B（人機決策層）")
    print("  Stage 15 → Flag 16 → Stage 17 → Flag 18 → prep_decisions.json")
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

    try:
        cleaned_path = load_cleaned_pdb(session_dir)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(cleaned_path))
    kept_chains = infer_kept_chains(session_dir, structure)

    print(f"\n  Target：{target.upper()} | PDB：{pdb_id}")
    print(f"  cleaned.pdb：{cleaned_path}")
    print(f"  kept_chains：{kept_chains}")

    try:
        r15 = stage15_cofactor(
            entry,
            structure,
            kept_chains,
            session_dir,
            non_interactive=non_interactive,
            auto_action=auto_hets,
        )
        r16 = flag16_propka(session_dir)
        r17 = stage17_override(
            entry,
            session_dir,
            non_interactive=non_interactive,
        )
        r18 = flag18_saltbridge(session_dir)
    except Exception as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)

    prep = load_prep_decisions(session_dir)
    prep["schema_version"] = PREP_DECISIONS_SCHEMA
    prep["target"] = target.upper()
    prep["pdb_id"] = pdb_id
    ch = entry.get("cocrystal_het") or entry.get("het_id")
    prep["cocrystal_het"] = ch.strip().upper() if isinstance(ch, str) and ch.strip() else None

    dec = dict(prep.get("decisions") or {})
    dec["stage15_cofactor"] = {
        "action": r15["action"],
        "kept": r15["kept"],
        "removed": r15["removed"],
        "detected": r15["detected"],
        "index_flagged": r15["index_flagged"],
        "status": r15["status"],
    }
    dec["stage16_propka"] = {
        "status": r16["status"],
        "flag": r16.get("flag"),
        "note": r16.get("note"),
    }
    dec["stage17_override"] = {
        "rules": r17["override_rules"],
        "family": r17["family"],
        "status": r17["status"],
    }
    dec["stage18_saltbridge"] = {
        "status": r18["status"],
        "flag": r18.get("flag"),
        "note": r18.get("note"),
    }
    prep["decisions"] = dec

    prep["engines"] = dict(prep.get("engines") or {})
    prep["engines"]["engine_b"] = {
        "status": "complete",
        "completed_at": datetime.now().isoformat(),
        "stages": ["15", "16", "17", "18"],
    }

    flags: list[str] = []
    if r16.get("flag"):
        flags.append(str(r16["flag"]))
    if r18.get("flag"):
        flags.append(str(r18["flag"]))
    if r15["status"] == "MANUAL":
        flags.append("STAGE15_MANUAL_COFACTOR")
    if r17["status"] == "MANUAL":
        flags.append("STAGE17_MANUAL_OVERRIDE")
    prep["flags"] = _merge_flags(prep.get("flags"), flags)

    prep["ready_for_engine_e"] = True

    prep_path = session_dir / "prep_decisions.json"
    save_prep_decisions(session_dir, prep)
    _update_manifest_engine_b(session_dir, prep_path)

    print("\n" + "=" * 55)
    print("  Engine B 完成 ✓")
    print("  ready_for_engine_e = true")
    print(f"  契約：{prep_path}")
    print(f"  Session：{session_dir}")
    print("=" * 55)

    return {
        "stage15": r15,
        "stage16": r16,
        "stage17": r17,
        "stage18": r18,
        "prep_decisions": str(prep_path),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="MolGate Module 3 Engine B")
    ap.add_argument("--session-dir", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument(
        "--index",
        type=Path,
        default=None,
        help="Master Index JSON（預設：MolGate 目錄下 master_index.json）",
    )
    ap.add_argument(
        "--skip-manifest-target-check",
        action="store_true",
        help="不與 molgate_run_manifest.json 的 target 交叉驗證",
    )
    ap.add_argument(
        "--non-interactive",
        action="store_true",
        help="Stage 15／17 不詢問（15 依 --auto-hets；17 採 Index override）",
    )
    ap.add_argument(
        "--auto-hets",
        choices=("remove_all", "keep_all"),
        default="remove_all",
        help="搭配 --non-interactive 的 Stage 15（預設 remove_all）",
    )
    args = ap.parse_args()

    sd = Path(args.session_dir)
    if not sd.is_dir():
        print(f"❌ Session 目錄不存在：{sd}")
        sys.exit(1)

    run_engine_b(
        sd,
        args.target,
        index_path=args.index,
        skip_manifest_target_check=args.skip_manifest_target_check,
        non_interactive=args.non_interactive,
        auto_hets=args.auto_hets,
    )
