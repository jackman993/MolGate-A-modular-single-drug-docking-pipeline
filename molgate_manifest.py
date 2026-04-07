"""
MolGate — Run manifest & registry（對齊各 Module 輸入/輸出路徑）

每個 session 目錄一個 molgate_run_manifest.json，記錄：
  - 藥物/標靶識別、各 stage 檔案相對路徑、pipeline 完成狀態
  - 下一步應讀哪個檔案（給 Module 2/3/… 用）

另可選寫入 molgate_sessions/registry.jsonl（一行一筆 run，方便外部管理）
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

MANIFEST_FILENAME = "molgate_run_manifest.json"
REGISTRY_FILENAME = "registry.jsonl"
MANIFEST_SCHEMA = "1.0"


def _now() -> str:
    return datetime.now().isoformat()


def manifest_path(session_dir: Path) -> Path:
    return session_dir / MANIFEST_FILENAME


def read_manifest(session_dir: Path) -> dict[str, Any]:
    p = manifest_path(session_dir)
    if not p.is_file():
        raise FileNotFoundError(f"找不到 manifest：{p}")
    return json.loads(p.read_text(encoding="utf-8"))


def write_manifest(session_dir: Path, data: dict[str, Any]) -> Path:
    data = dict(data)
    data.setdefault("_schema_version", MANIFEST_SCHEMA)
    data["updated_at"] = _now()
    p = manifest_path(session_dir)
    session_dir.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return p


def rel_path(session_dir: Path, file_path: Path | str) -> str:
    fp = Path(file_path)
    try:
        return fp.resolve().relative_to(session_dir.resolve()).as_posix()
    except ValueError:
        return fp.as_posix()


def build_paths_block(
    session_dir: Path,
    pdb_id: str,
    *,
    stage01: dict[str, Any] | None = None,
    stage02: dict[str, Any] | None = None,
    stage03: dict[str, Any] | None = None,
    stage04: dict[str, Any] | None = None,
    stage05: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """統一 paths / next_step 結構。"""
    sd = session_dir.resolve()
    rel: dict[str, str] = {
        "stage_01_smiles": "stage_01_smiles.json",
        "stage_02_pdb": "stage_02_pdb.json",
        "stage_03_index": "stage_03_index.json",
        "session_meta": "session_meta.json",
        "manifest": MANIFEST_FILENAME,
    }
    pdb_fn = f"pdb/{pdb_id.upper()}.pdb"
    rel["receptor_pdb"] = pdb_fn

    if stage04 is not None:
        rel["stage_04_protonation"] = "stage_04_protonation.json"
    if stage05 is not None:
        rel["stage_05_tautomer"] = "stage_05_tautomer.json"

    next_m2: dict[str, Any] = {
        "read_canonical_smiles_from": "stage_01_smiles.json",
        "field": "canonical_smiles",
        "read_receptor_pdb_from": "stage_02_pdb.json",
        "field_pdb_path": "pdb_path",
    }
    next_m2b_done: dict[str, Any] | None = None
    if stage05 is not None:
        next_m2b_done = {
            "ligand_smiles_for_stage09": {
                "file": "stage_05_tautomer.json",
                "field": "selected_smiles",
            }
        }

    out: dict[str, Any] = {
        "session_dir": str(sd),
        "relative": rel,
        "absolute_pdb_path": str(sd / pdb_fn) if pdb_id else None,
        "next_step": {
            "after_module1": {
                "module2_engine_b": next_m2,
            },
        },
    }
    if next_m2b_done:
        out["next_step"]["after_module2_engine_b"] = next_m2b_done
    return out


def finalize_module1_manifest(
    session_dir: Path,
    meta: dict[str, Any],
    r1: dict[str, Any],
    r2: dict[str, Any],
    r3: dict[str, Any],
    *,
    drug_id: str | None = None,
    drug_label: str | None = None,
) -> dict[str, Any]:
    sid = meta.get("session_id", session_dir.name)
    pdb_id = (r2.get("pdb_id") or r3.get("pdb_id") or meta.get("pdb_id") or "").upper()
    did = drug_id if drug_id is not None else meta.get("drug_id")
    dlab = drug_label if drug_label is not None else meta.get("drug_label")

    manifest: dict[str, Any] = {
        "_schema_version": MANIFEST_SCHEMA,
        "run_id": sid,
        "drug_id": did,
        "drug_label": dlab,
        "target": r3.get("target") or meta.get("target"),
        "pdb_id": pdb_id,
        "index_path": meta.get("index_path"),
        "created_at": meta.get("created", _now()),
        "updated_at": _now(),
        "pipeline": {
            "module1": {
                "status": "complete",
                "completed_at": _now(),
                "stages": ["01", "02", "03"],
            },
            "module2_engine_b": {"status": "pending", "stages": []},
            "module2_engine_c": {"status": "pending", "stages": []},
        },
        "paths": build_paths_block(session_dir, pdb_id, stage01=r1, stage02=r2, stage03=r3),
        "summary": {
            "stage01_status": r1.get("status"),
            "recommend_resolve_stereo_before_module2": r1.get(
                "recommend_resolve_stereo_before_module2", False
            ),
            "stage02_status": r2.get("status"),
            "stage03_status": r3.get("status"),
            "cocrystal_het": r3.get("cocrystal_het"),
        },
    }
    write_manifest(session_dir, manifest)
    return manifest


def finalize_module2_engine_b_manifest(
    session_dir: Path,
    r4: dict[str, Any],
    r5: dict[str, Any],
) -> dict[str, Any]:
    m = read_manifest(session_dir)
    pdb_id = m.get("pdb_id") or ""
    m["paths"] = build_paths_block(
        session_dir,
        pdb_id,
        stage04=r4,
        stage05=r5,
    )
    m["pipeline"] = dict(m.get("pipeline") or {})
    m["pipeline"]["module2_engine_b"] = {
        "status": "complete",
        "completed_at": _now(),
        "stages": ["04", "05"],
    }
    if "module2_engine_c" not in m["pipeline"]:
        m["pipeline"]["module2_engine_c"] = {"status": "pending", "stages": []}
    m["summary"] = dict(m.get("summary") or {})
    m["summary"]["stage04_status"] = r4.get("status")
    m["summary"]["stage05_status"] = r5.get("status")
    m["summary"]["final_smiles_post_m2"] = r5.get("selected_smiles")
    write_manifest(session_dir, m)
    return m


def finalize_module2_engine_c_manifest(
    session_dir: Path,
    stage09: dict[str, Any],
    ligand_pdbqt_file: Path,
) -> dict[str, Any]:
    """Stage 09：配體 PDBQT 產物寫入 manifest，供 Vina / Module 3 讀取。"""
    m = read_manifest(session_dir)
    ligand_pdbqt_file = ligand_pdbqt_file.resolve()
    rel = rel_path(session_dir, ligand_pdbqt_file)

    paths = dict(m.get("paths") or {})
    rmap = dict(paths.get("relative") or {})
    rmap["ligand_pdbqt"] = rel
    rmap["stage_09_ligand_pdbqt"] = "stage_09_ligand_pdbqt.json"
    paths["relative"] = rmap
    paths["ligand_pdbqt_absolute"] = str(ligand_pdbqt_file)

    ns = dict(paths.get("next_step") or {})
    ns["after_module2_engine_c"] = {
        "ligand_pdbqt": {
            "relative": rel,
            "absolute": str(ligand_pdbqt_file),
        },
        "stage_record": "stage_09_ligand_pdbqt.json",
        "hint": "下一步：受體 PDBQT（Module 3）或 AutoDock Vina --ligand",
    }
    paths["next_step"] = ns
    m["paths"] = paths

    m["pipeline"] = dict(m.get("pipeline") or {})
    m["pipeline"]["module2_engine_c"] = {
        "status": "complete",
        "completed_at": _now(),
        "stages": ["09"],
    }
    m["summary"] = dict(m.get("summary") or {})
    m["summary"]["stage09_status"] = stage09.get("status")
    m["summary"]["stage09_engine"] = stage09.get("engine")
    m["summary"]["ligand_pdbqt_relative"] = rel
    write_manifest(session_dir, m)
    return m


def append_registry(
    base_sessions_dir: Path,
    entry: dict[str, Any],
) -> None:
    """追加一行 JSON 到 molgate_sessions/registry.jsonl"""
    base_sessions_dir.mkdir(parents=True, exist_ok=True)
    path = base_sessions_dir / REGISTRY_FILENAME
    line = json.dumps(entry, ensure_ascii=False) + "\n"
    with open(path, "a", encoding="utf-8") as f:
        f.write(line)


def try_load_manifest_or_rebuild(session_dir: Path) -> dict[str, Any]:
    """若無 manifest，嘗試用 session_meta + 既有 stage 檔重建（舊 session 相容）。"""
    p = manifest_path(session_dir)
    if p.is_file():
        return read_manifest(session_dir)

    meta_p = session_dir / "session_meta.json"
    if not meta_p.is_file():
        raise FileNotFoundError("無 manifest 且無 session_meta.json，無法對齊路徑")
    meta = json.loads(meta_p.read_text(encoding="utf-8"))
    r1 = json.loads((session_dir / "stage_01_smiles.json").read_text(encoding="utf-8"))
    r2 = json.loads((session_dir / "stage_02_pdb.json").read_text(encoding="utf-8"))
    r3 = json.loads((session_dir / "stage_03_index.json").read_text(encoding="utf-8"))
    return finalize_module1_manifest(session_dir, meta, r1, r2, r3)


def validate_target_matches_manifest(session_dir: Path, target: str) -> None:
    m = read_manifest(session_dir)
    mt = (m.get("target") or "").strip().upper()
    t = target.strip().upper()
    if mt and t and mt != t:
        raise ValueError(
            f"CLI --target={target!r} 與 manifest 的 target={m.get('target')!r} 不一致，"
            "請使用同一標靶或換新 session。"
        )
