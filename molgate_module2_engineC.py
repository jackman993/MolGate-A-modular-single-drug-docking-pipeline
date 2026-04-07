"""
MolGate — Module 2 Engine C：配體 PDBQT（Stage 09）

預設（--engine meeko，快速）：僅 Meeko。
備援：--engine auto 為 Meeko → obabel；--engine obabel 僅 Open Babel。

輸入 SMILES：
  - 優先 stage_05_tautomer.json → selected_smiles（若已跑 Engine B）
  - 否則 stage_01_smiles.json → canonical_smiles

輸出：
  - ligand/ligand.pdbqt
  - stage_09_ligand_pdbqt.json
  - 更新 molgate_run_manifest.json
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

_PKG_ROOT = Path(__file__).resolve().parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

from molgate_manifest import (
    append_registry,
    finalize_module2_engine_c_manifest,
    try_load_manifest_or_rebuild,
    validate_target_matches_manifest,
)

from rdkit import Chem
from rdkit.Chem import AllChem

SCHEMA_VERSION = "0.1"


def save_stage(session_dir: Path, stage: int, name: str, data: dict[str, Any]) -> dict[str, Any]:
    data = dict(data)
    data["_stage"] = stage
    data["_name"] = name
    data["_timestamp"] = datetime.now().isoformat()
    data["_schema_version"] = SCHEMA_VERSION
    session_dir.mkdir(parents=True, exist_ok=True)
    path = session_dir / f"stage_{stage:02d}_{name}.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return data


def load_json(p: Path) -> dict[str, Any]:
    return json.loads(p.read_text(encoding="utf-8"))


def resolve_input_smiles(session_dir: Path) -> tuple[str, str]:
    p5 = session_dir / "stage_05_tautomer.json"
    if p5.is_file():
        d = load_json(p5)
        smi = d.get("selected_smiles")
        if smi:
            return smi, "stage_05_tautomer.json:selected_smiles"
    p1 = session_dir / "stage_01_smiles.json"
    if not p1.is_file():
        raise FileNotFoundError("找不到 stage_05_tautomer.json 或 stage_01_smiles.json")
    d = load_json(p1)
    smi = d.get("canonical_smiles")
    if not smi:
        raise ValueError("stage_01_smiles.json 缺少 canonical_smiles")
    return smi, "stage_01_smiles.json:canonical_smiles (未跑 Engine B)"


def build_3d_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("無法自 SMILES 建分子")
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    try:
        params = AllChem.ETKDGv3()
    except AttributeError:
        params = AllChem.ETKDG()
    params.randomSeed = 0xC0FFEE
    err = AllChem.EmbedMolecule(mol, params)
    if err < 0:
        err = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if err < 0:
        raise RuntimeError("RDKit EmbedMolecule 失敗，無法產生 3D 構象")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass
    return mol


def _coerce_pdbqt_string(raw: Any) -> str:
    """
    Meeko 0.7+ 的 write_string 可能回傳 str，或 (pdbqt_str, …) 元組。
    """
    if raw is None:
        return ""
    if isinstance(raw, str):
        return raw
    if isinstance(raw, bytes):
        return raw.decode("utf-8", errors="replace")
    if isinstance(raw, (tuple, list)) and len(raw) > 0:
        first = raw[0]
        if isinstance(first, str):
            return first
        if isinstance(first, bytes):
            return first.decode("utf-8", errors="replace")
    raise TypeError(f"PDBQT 輸出型別無法轉成字串: {type(raw)!r}")


def engine_meeko(mol: Chem.Mol, out_pdbqt: Path) -> dict[str, Any]:
    """
    Meeko 0.6+：MoleculePreparation(mol) 回傳 setup 列表，PDBQT 用 PDBQTWriterLegacy.write_string。
    較舊版：prepare(mol) + setup 上的 write_* 方法。
    """
    try:
        from meeko import MoleculePreparation
    except ImportError as e:
        msg = str(e)
        hint = ""
        if "gemmi" in msg.lower():
            hint = " 請執行：pip install gemmi（Meeko 依賴 gemmi）。"
        raise RuntimeError(f"Meeko 未安裝或無法 import: {e}.{hint}") from e

    preparator = MoleculePreparation()
    setups = None
    api_route = ""

    # 0.6+ 慣用 __call__(mol)；舊版用 prepare(mol)
    try:
        setups = preparator(mol)
        api_route = "__call__"
    except Exception:
        pass
    if not setups and hasattr(preparator, "prepare"):
        try:
            setups = preparator.prepare(mol)
            api_route = "prepare"
        except Exception as e:
            raise RuntimeError(f"Meeko prepare 失敗: {e}") from e

    if not setups:
        raise RuntimeError("Meeko 未回傳任何 molecule setup（檢查分子是否含氫與 3D）")

    setup = setups[0]
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    pdbqt_text: str | None = None
    detail = ""

    try:
        from meeko import PDBQTWriterLegacy

        pdbqt_text = PDBQTWriterLegacy.write_string(setup)
        detail = "PDBQTWriterLegacy.write_string"
    except Exception as e_legacy:
        err_legacy = str(e_legacy)
        try:
            from meeko import PDBQTWriter

            pdbqt_text = PDBQTWriter.write_string(setup)
            detail = "PDBQTWriter.write_string"
        except Exception:
            wfile = getattr(setup, "write_pdbqt_file", None)
            if callable(wfile):
                wfile(str(out_pdbqt))
                return {
                    "engine": "meeko",
                    "detail": "write_pdbqt_file",
                    "api": api_route,
                }
            writer = getattr(setup, "write_pdbqt_string", None)
            if writer is None:
                writer = getattr(setup, "write_string", None)
            if callable(writer):
                pdbqt_text = writer()
                detail = "setup.write_pdbqt_string/write_string"
            else:
                raise RuntimeError(
                    "Meeko 無法輸出 PDBQT：請升級 meeko（建議 pip install -U meeko）"
                    f" 並確認可用 PDBQTWriterLegacy。Legacy 錯誤：{err_legacy}"
                ) from e_legacy

    try:
        pdbqt_text = _coerce_pdbqt_string(pdbqt_text)
    except TypeError as e:
        raise RuntimeError(f"Meeko PDBQT 輸出格式異常: {e}") from e

    if not pdbqt_text or not pdbqt_text.strip():
        raise RuntimeError("Meeko 產生的 PDBQT 字串為空")
    out_pdbqt.write_text(pdbqt_text, encoding="utf-8")
    return {"engine": "meeko", "detail": detail, "api": api_route}


def _obabel_run(obabel: str, argv: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [obabel, *argv],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def _write_temp_smi(smiles: str) -> Path:
    """單行 SMILES 寫入暫存 .smi，避免 stdin 與參數順序在 OB 3.x 上觸發錯誤。"""
    tf = tempfile.NamedTemporaryFile(
        mode="w", suffix=".smi", delete=False, encoding="utf-8", newline="\n"
    )
    try:
        tf.write(smiles.strip() + "\n")
        path = Path(tf.name)
    finally:
        tf.close()
    return path


def _obabel_try_convert(
    obabel: str,
    argv: list[str],
    out_pdbqt: Path,
    *,
    strategy: str,
) -> dict[str, Any] | None:
    proc = _obabel_run(obabel, argv)
    if proc.returncode == 0 and out_pdbqt.is_file() and out_pdbqt.stat().st_size >= 50:
        return {
            "engine": "obabel",
            "obabel_path": obabel,
            "strategy": strategy,
            "argv": argv,
        }
    return None


def engine_obabel(smiles: str, out_pdbqt: Path) -> dict[str, Any]:
    """
    Open Babel → PDBQT。

    不同 OB 版本對「--gen3d / -O」順序要求不一，故依序嘗試多組 argv。
    """
    obabel = shutil.which("obabel")
    if not obabel:
        raise RuntimeError("找不到 obabel（請安裝 Open Babel 並加入 PATH）")

    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    err_parts: list[str] = []
    sp = str(out_pdbqt.resolve())

    smi_path: Path | None = None
    try:
        smi_path = _write_temp_smi(smiles)
        si = str(smi_path.resolve())

        # 策略 A：多種 OB 3.x 可接受的參數順序
        smi_to_pdbqt_variants = [
            [si, "--gen3d", "-O", sp],
            [si, "-O", sp, "--gen3d"],
            ["-ismi", si, "-O", sp, "--gen3d"],
            ["-ismi", si, "--gen3d", "-O", sp],
        ]
        for argv in smi_to_pdbqt_variants:
            ok = _obabel_try_convert(obabel, argv, out_pdbqt, strategy="smi_to_pdbqt")
            if ok:
                return ok
            proc = _obabel_run(obabel, argv)
            err_parts.append((proc.stderr or proc.stdout or "").strip()[:400])

        # 策略 B：smi → 暫存 pdb → pdbqt
        tmp_pdb = out_pdbqt.with_name(out_pdbqt.stem + "_obabel_tmp.pdb")
        tp = str(tmp_pdb.resolve())
        try:
            smi_to_pdb_variants = [
                [si, "--gen3d", "-O", tp],
                [si, "-O", tp, "--gen3d"],
                ["-ismi", si, "-O", tp, "--gen3d"],
            ]
            proc1_ok = False
            for argv1 in smi_to_pdb_variants:
                proc1 = _obabel_run(obabel, argv1)
                if proc1.returncode == 0 and tmp_pdb.is_file():
                    proc1_ok = True
                    break
                err_parts.append((proc1.stderr or proc1.stdout or "").strip()[:400])
            if proc1_ok:
                pdb_to_pdbqt_variants = [
                    [tp, "-O", sp],
                    ["-ipdb", tp, "-O", sp],
                ]
                for argv2 in pdb_to_pdbqt_variants:
                    ok = _obabel_try_convert(
                        obabel, argv2, out_pdbqt, strategy="smi_to_pdb_then_pdbqt"
                    )
                    if ok:
                        return ok
                    proc2 = _obabel_run(obabel, argv2)
                    err_parts.append((proc2.stderr or proc2.stdout or "").strip()[:400])
            else:
                err_parts.append("smi→pdb 所有 argv 變體皆失敗")
        finally:
            if tmp_pdb.is_file():
                try:
                    tmp_pdb.unlink()
                except OSError:
                    pass
    finally:
        if smi_path is not None and smi_path.is_file():
            try:
                smi_path.unlink()
            except OSError:
                pass

    err_tail = " | ".join(e for e in err_parts if e)
    raise RuntimeError(
        "obabel 無法產生 PDBQT（已嘗試暫存 .smi 與 pdb 中繼）。"
        " 可能原因：Open Babel 未含 pdbqt 外掛，或 Python 3.13 與 conda-forge openbabel 版本不相容。"
        " 建議：① pip install gemmi meeko 後使用 --engine meeko；② 另建 Python 3.11/3.12 conda 環境再裝 openbabel；③ Windows 安裝官方 Open Babel 並確認 obabel 在 PATH。"
        f" 末次訊息：{err_tail}"
    )


def run_engine_c(
    session_dir: Path,
    target: str,
    engine: str = "meeko",
) -> dict[str, Any] | None:
    print("=" * 58)
    print("  MolGate — Module 2 Engine C（Stage 09：配體 PDBQT）")
    print("=" * 58)

    try:
        try_load_manifest_or_rebuild(session_dir)
        validate_target_matches_manifest(session_dir, target)
    except (ValueError, FileNotFoundError) as e:
        print(f"❌ {e}")
        return None

    smiles, smi_src = resolve_input_smiles(session_dir)
    print(f"\n  SMILES 來源：{smi_src}")
    print(f"  SMILES：{smiles[:120]}{'...' if len(smiles) > 120 else ''}")

    s1_path = session_dir / "stage_01_smiles.json"
    if s1_path.is_file():
        s1 = load_json(s1_path)
        if s1.get("recommend_resolve_stereo_before_module2"):
            print(
                "\n  ⚠ WARN：Stage 01 曾標記立體未補齊；PDBQT 仍會產生，解讀風險自負。\n"
            )

    lig_dir = session_dir / "ligand"
    lig_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = lig_dir / "ligand.pdbqt"

    mol_3d = build_3d_mol(smiles)

    print(f"\n  引擎模式：{engine}" + ("（預設快速：僅 Meeko）" if engine == "meeko" else ""))

    engines_order: list[str]
    if engine == "auto":
        engines_order = ["meeko", "obabel"]
    elif engine == "meeko":
        engines_order = ["meeko"]
    elif engine == "obabel":
        engines_order = ["obabel"]
    else:
        print(f"❌ 未知 --engine：{engine}")
        return None

    last_err: str | None = None
    meta: dict[str, Any] = {}
    for eng in engines_order:
        try:
            if eng == "meeko":
                meta = engine_meeko(mol_3d, out_pdbqt)
                break
            if eng == "obabel":
                meta = engine_obabel(smiles, out_pdbqt)
                break
        except Exception as e:
            last_err = f"{eng}: {e}"
            print(f"  ⚠ {last_err}")
            continue
    else:
        result = {
            "status": "FAIL",
            "reason": "所有引擎皆失敗",
            "last_error": last_err,
            "tried": engines_order,
        }
        save_stage(session_dir, 9, "ligand_pdbqt", result)
        print(f"\n❌ Engine C 失敗：{last_err}")
        return None

    result = {
        "status": "PASS",
        "input_smiles": smiles,
        "smiles_source": smi_src,
        "engine": meta.get("engine"),
        "engine_detail": meta,
        "ligand_pdbqt": str(out_pdbqt.resolve()),
    }
    save_stage(session_dir, 9, "ligand_pdbqt", result)

    finalize_module2_engine_c_manifest(session_dir, result, out_pdbqt)
    append_registry(
        session_dir.parent,
        {
            "event": "module2_engine_c_complete",
            "session_id": session_dir.name,
            "session_dir": str(session_dir.resolve()),
            "target": target,
            "ligand_pdbqt": str(out_pdbqt.resolve()),
            "engine": result.get("engine"),
            "timestamp": datetime.now().isoformat(),
        },
    )

    print("\n" + "=" * 58)
    print("  Engine C 完成 ✓")
    print(f"  配體 PDBQT：{out_pdbqt}")
    print(f"  引擎：{result.get('engine')}")
    print(f"  Manifest：{session_dir / 'molgate_run_manifest.json'}")
    print("=" * 58)
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description="MolGate Module 2 Engine C — ligand PDBQT")
    ap.add_argument("--session-dir", required=True, type=Path, help="session 目錄（與 Module 1 相同）")
    ap.add_argument("--target", required=True, help="須與 manifest 一致，如 COX2")
    ap.add_argument(
        "--engine",
        choices=("auto", "meeko", "obabel"),
        default="meeko",
        help="預設 meeko（快速）；auto＝Meeko 失敗再 obabel；obabel＝僅 Open Babel",
    )
    args = ap.parse_args()
    sd = args.session_dir
    if not sd.is_dir():
        print(f"❌ session 目錄不存在：{sd}")
        return 1
    r = run_engine_c(sd, args.target, engine=args.engine)
    return 0 if r is not None else 1


if __name__ == "__main__":
    raise SystemExit(main())
