"""
MolGate — Module 3: Protein Preparation（Engine C）
Stage 19：自 cleaned.pdb 產出 apo 受體 receptor.pdbqt（auto：預設先 Meeko，再 Open Babel，最後 MGLTools）
Stage 20：幾何／結構 QC（輕量占位）

輸入：protein/cleaned.pdb、prep_decisions.json（須 ready_for_engine_e = true）
依賴：biopython、RDKit（Meeko 路徑）；可選 PATH 有 obabel（部分 Windows 安裝無 pdbqt 寫出）；可選 MGLTools prepare_receptor4.py
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa

_PKG_ROOT = Path(__file__).resolve().parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(errors="replace")

SCHEMA_VERSION = "0.3"


class ApoProteinSelect(Select):
    """保留胺基酸主體／常見修飾殘基，排除水與非蛋白 HET。"""

    def accept_residue(self, residue):
        rn = residue.get_resname().strip().upper()
        if rn in ("HOH", "WAT"):
            return False
        return is_aa(residue, standard=False)


def save_stage(session_dir: Path, stage: int, name: str, data: dict[str, Any]) -> dict[str, Any]:
    data = dict(data)
    data["_stage"] = stage
    data["_name"] = name
    data["_timestamp"] = datetime.now().isoformat()
    data["_schema_version"] = SCHEMA_VERSION
    path = session_dir / f"stage_{stage:02d}_{name}.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return data


def load_cleaned_pdb(session_dir: Path) -> Path:
    p = session_dir / "protein" / "cleaned.pdb"
    if not p.is_file():
        raise FileNotFoundError("找不到 protein/cleaned.pdb，請先跑 Engine A")
    return p


def load_prep_decisions(session_dir: Path) -> dict[str, Any]:
    path = session_dir / "prep_decisions.json"
    if not path.is_file():
        raise FileNotFoundError("找不到 prep_decisions.json，請先跑 Engine B")
    return json.loads(path.read_text(encoding="utf-8"))


def _write_apo_pdb(structure, out_pdb: Path) -> None:
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb), ApoProteinSelect())


def _run_cmd(argv: list[str], *, cwd: Path | None = None, timeout: int = 300) -> tuple[int, str]:
    try:
        r = subprocess.run(
            argv,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
            cwd=str(cwd) if cwd else None,
        )
    except (OSError, subprocess.TimeoutExpired) as e:
        return -1, str(e)
    out = ((r.stderr or "") + "\n" + (r.stdout or "")).strip()
    return r.returncode, out


def _obabel_format_probe(exe: str) -> str:
    """回傳 obabel -L formats 摘要，用於診斷是否支援 pdbqt 寫出。"""
    code, msg = _run_cmd([exe, "-L", "formats"], timeout=30)
    if code != 0 or not msg:
        return f"無法列出格式（exit {code}）：{msg[:400]}"
    low = msg.lower()
    if "pdbqt" in low:
        return "此 obabel 的格式清單內含 pdbqt（理論上可寫 pdbqt）；若仍失敗請檢查參數或改用 --engine meeko / mgltools。"
    return (
        "此 obabel 的格式清單內未見 pdbqt — 常見於精簡編譯，無法用 obabel 寫 receptor.pdbqt。"
        " 請改用：--engine meeko，或安裝含 pdbqt 的 Open Babel，或使用 MGLTools prepare_receptor4.py。"
    )


def _obabel_executables() -> list[str]:
    out: list[str] = []
    for name in ("obabel", "babel"):
        p = shutil.which(name)
        if p and p not in out:
            out.append(p)
    return out


def _obabel_pdb_to_pdbqt(pdb_in: Path, pdbqt_out: Path) -> tuple[bool, str]:
    """多種 obabel/babel 寫法；Windows 若編譯未含 pdbqt 會回 False。"""
    exes = _obabel_executables()
    if not exes:
        return False, "PATH 中找不到 obabel 或 babel"

    pdb_in = pdb_in.resolve()
    pdbqt_out = pdbqt_out.resolve()
    pdbqt_out.parent.mkdir(parents=True, exist_ok=True)
    if pdbqt_out.is_file():
        pdbqt_out.unlink()

    # Open Babel 3.x：輸出格式建議用 -of pdbqt，且 -O 指定檔名；順序依官方說明調整
    variants: list[list[str]] = []
    for obabel in exes:
        variants.extend(
            [
                [obabel, "-ipdb", str(pdb_in), "-O", str(pdbqt_out), "-of", "pdbqt"],
                [obabel, str(pdb_in), "-O", str(pdbqt_out), "-of", "pdbqt"],
                # 舊式 -opdbqt（部分環境仍吃）
                [obabel, "-ipdb", str(pdb_in), "-opdbqt", "-O", str(pdbqt_out)],
                [obabel, str(pdb_in), "-O", str(pdbqt_out), "-opdbqt"],
                [obabel, "-ipdb", str(pdb_in), "-O", str(pdbqt_out), "-opdbqt"],
                [obabel, str(pdb_in), "-O", str(pdbqt_out)],
                [obabel, "-ipdb", str(pdb_in), "-O", str(pdbqt_out)],
            ]
        )

    last_err = ""
    primary_exe = exes[0]
    for argv in variants:
        code, msg = _run_cmd(argv)
        if pdbqt_out.is_file() and pdbqt_out.stat().st_size >= 20:
            return True, ""
        last_err = msg or f"exit {code}"

    probe = _obabel_format_probe(primary_exe)
    return False, f"{last_err}\n[obabel 診斷] {probe}"


def _meeko_try_write(
    pdb_in: Path,
    pdbqt_out: Path,
    *,
    allow_bad_res: bool,
) -> tuple[bool, str]:
    argv = [
        sys.executable,
        "-m",
        "meeko.cli.mk_prepare_receptor",
        "--read_pdb",
        str(pdb_in),
        "--write_pdbqt",
        str(pdbqt_out),
    ]
    if allow_bad_res:
        argv.append("--allow_bad_res")
    code, msg = _run_cmd(argv, cwd=pdbqt_out.parent)
    if pdbqt_out.is_file() and pdbqt_out.stat().st_size >= 20:
        return True, ""
    stem = pdbqt_out.with_suffix("")
    alt = Path(str(stem) + "_rigid.pdbqt")
    if alt.is_file() and alt.stat().st_size >= 20:
        alt.replace(pdbqt_out)
        return True, ""
    return False, msg or f"meeko exit {code}"


def _meeko_pdb_to_pdbqt(pdb_in: Path, pdbqt_out: Path) -> tuple[bool, str]:
    """Meeko mk_prepare_receptor（不依賴 obabel pdbqt）。失敗時以 --allow_bad_res 重試（交聯殘基等）。"""
    pdb_in = pdb_in.resolve()
    pdbqt_out = pdbqt_out.resolve()
    pdbqt_out.parent.mkdir(parents=True, exist_ok=True)

    def _cleanup_out() -> None:
        if pdbqt_out.is_file():
            pdbqt_out.unlink()
        stem = pdbqt_out.with_suffix("")
        alt = Path(str(stem) + "_rigid.pdbqt")
        if alt.is_file():
            alt.unlink()

    _cleanup_out()
    ok, err = _meeko_try_write(pdb_in, pdbqt_out, allow_bad_res=False)
    if ok:
        return True, ""

    _cleanup_out()
    ok2, err2 = _meeko_try_write(pdb_in, pdbqt_out, allow_bad_res=True)
    if ok2:
        print(
            "  ⚠ Meeko 已使用 --allow_bad_res（模板無法匹配的殘基已略過， docking 前請確認）",
        )
        return True, ""

    # 2HU4 這類多聚體偶發跨殘基錯誤時，嘗試排除問題鏈再重試
    bond_pairs = re.findall(r"\(([A-Za-z0-9]):\d+,\s*([A-Za-z0-9]):\d+\)", f"{err}\n{err2}")
    bad_chains = sorted({c for a, b in bond_pairs for c in (a, b)})
    if bad_chains:
        sanitized_pdb = pdbqt_out.parent / "apo_for_docking_sanitized.pdb"
        try:
            kept = 0
            with pdb_in.open("r", encoding="utf-8", errors="replace") as fin, sanitized_pdb.open(
                "w",
                encoding="utf-8",
            ) as fout:
                for line in fin:
                    rec = line[:6]
                    if rec.startswith(("ATOM  ", "HETATM", "ANISOU")) and len(line) > 21:
                        if line[21] in bad_chains:
                            continue
                        kept += 1
                    fout.write(line)
            if kept > 0:
                _cleanup_out()
                ok3, err3 = _meeko_try_write(sanitized_pdb, pdbqt_out, allow_bad_res=False)
                if ok3:
                    print(
                        f"  ⚠ Meeko fallback：排除問題鏈 {bad_chains} 後成功產生 receptor.pdbqt",
                    )
                    return True, ""
                _cleanup_out()
                ok4, err4 = _meeko_try_write(sanitized_pdb, pdbqt_out, allow_bad_res=True)
                if ok4:
                    print(
                        f"  ⚠ Meeko fallback：排除問題鏈 {bad_chains} + --allow_bad_res 後成功",
                    )
                    return True, ""
                err2 = f"{err2}\nmeeko（sanitize chains {bad_chains}）：{err3}\nmeeko（sanitize+allow_bad_res）：{err4}"
        except OSError as e:
            err2 = f"{err2}\nmeeko sanitize fallback 失敗: {e}"

    return False, f"meeko：{err}\nmeeko（--allow_bad_res）：{err2}"


def _find_mgltools_prepare_script() -> tuple[Path, Path] | None:
    """
    回傳 (python.exe, prepare_receptor4.py)；找不到則 None。
    可用環境變數 MGLTOOLS_ROOT 覆寫安裝根目錄。
    """
    roots: list[Path] = []
    env_root = os.environ.get("MGLTOOLS_ROOT", "").strip()
    if env_root:
        roots.append(Path(env_root))
    roots.extend(
        [
            Path(r"C:\Program Files (x86)\MGLTools-1.5.7"),
            Path(r"C:\Program Files\MGLTools-1.5.7"),
        ]
    )
    script_names = (
        Path("Lib") / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
        Path("lib") / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
        Path("MGLToolsPckgs") / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py",
    )
    seen: set[Path] = set()
    for root in roots:
        try:
            root = root.resolve()
        except OSError:
            continue
        if root in seen:
            continue
        seen.add(root)
        py = root / "python.exe"
        if not py.is_file():
            continue
        for rel in script_names:
            script = root / rel
            if script.is_file():
                return py, script
    return None


def _mgltools_pdb_to_pdbqt(pdb_in: Path, pdbqt_out: Path) -> tuple[bool, str]:
    found = _find_mgltools_prepare_script()
    if not found:
        return False, "找不到 MGLTools（prepare_receptor4.py）。可設定 MGLTOOLS_ROOT 或安裝 MGLTools 1.5.x。"
    py_exe, script = found
    pdb_in = pdb_in.resolve()
    pdbqt_out = pdbqt_out.resolve()
    pdbqt_out.parent.mkdir(parents=True, exist_ok=True)
    if pdbqt_out.is_file():
        pdbqt_out.unlink()
    argv = [
        str(py_exe),
        str(script),
        "-r",
        str(pdb_in),
        "-o",
        str(pdbqt_out),
        "-A",
        "hydrogens",
    ]
    code, msg = _run_cmd(argv, timeout=600)
    if pdbqt_out.is_file() and pdbqt_out.stat().st_size >= 20:
        return True, ""
    return False, msg or f"MGLTools prepare_receptor4 exit {code}"


def _pdb_to_pdbqt(
    pdb_in: Path,
    pdbqt_out: Path,
    *,
    engine: str,
) -> tuple[bool, str, str]:
    """
    engine: auto | obabel | meeko | mgltools
    回傳 (ok, error_message, used_engine_label)
    """
    if engine == "obabel":
        ok, err = _obabel_pdb_to_pdbqt(pdb_in, pdbqt_out)
        return ok, err, "obabel"
    if engine == "meeko":
        ok, err = _meeko_pdb_to_pdbqt(pdb_in, pdbqt_out)
        return ok, err, "meeko"
    if engine == "mgltools":
        ok, err = _mgltools_pdb_to_pdbqt(pdb_in, pdbqt_out)
        return ok, err, "mgltools"

    # auto：meeko → obabel → MGLTools（與多數 Windows 環境一致：先求 Vina 相容 PDBQT，再備援）
    ok, err = _meeko_pdb_to_pdbqt(pdb_in, pdbqt_out)
    if ok:
        return True, "", "meeko"
    ok2, err2 = _obabel_pdb_to_pdbqt(pdb_in, pdbqt_out)
    if ok2:
        return True, "", "obabel"
    ok3, err3 = _mgltools_pdb_to_pdbqt(pdb_in, pdbqt_out)
    if ok3:
        return True, "", "mgltools"
    return False, f"meeko: {err}\nobabel: {err2}\nmgltools: {err3}", ""


def stage19_receptor(
    session_dir: Path,
    pdb_id: str,
    cleaned_path: Path,
    *,
    engine: str = "auto",
) -> dict[str, Any]:
    print("\n[Stage 19] Receptor 準備（apo → PDBQT）...")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(cleaned_path))

    apo_pdb = session_dir / "protein" / "apo_for_docking.pdb"
    _write_apo_pdb(structure, apo_pdb)
    print(f"  ✓ 寫入 apo PDB：{apo_pdb}")

    out_pdbqt = session_dir / "protein" / "receptor.pdbqt"
    ok, err, used = _pdb_to_pdbqt(apo_pdb, out_pdbqt, engine=engine)
    if ok:
        print(f"  ✓ receptor.pdbqt：{out_pdbqt}（引擎：{used}）")
        return save_stage(session_dir, 19, "receptor_pdbqt", {
            "status": "PASS",
            "engine": used,
            "apo_pdb": str(apo_pdb.as_posix()),
            "receptor_pdbqt": str(out_pdbqt.as_posix()),
        })

    print(f"  ❌ PDBQT 轉換失敗：{err}")
    return save_stage(session_dir, 19, "receptor_pdbqt", {
        "status": "FAIL",
        "engine": used or engine,
        "apo_pdb": str(apo_pdb.as_posix()),
        "receptor_pdbqt": None,
        "error": err,
    })


def stage20_geometry_qc(session_dir: Path, pdb_id: str, cleaned_path: Path) -> dict[str, Any]:
    print("\n[Stage 20] 受體幾何 QC（輕量）...")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(cleaned_path))

    n_chain = 0
    n_aa = 0
    for model in structure:
        for chain in model:
            n_chain += 1
            for res in chain:
                if is_aa(res, standard=False) and res.get_resname().strip().upper() not in (
                    "HOH",
                    "WAT",
                ):
                    n_aa += 1

    status = "PASS" if n_aa > 0 else "WARN"
    if status == "WARN":
        print("  ⚠ 可辨識胺基酸殘基過少，請檢查 cleaned.pdb")
    else:
        print(f"  ✓ Chains（模型內列舉）：{n_chain} | 胺基酸殘基（約）：{n_aa}")

    print("  → Status:", status)
    print("  （正式版可接 RDKit／clash 檢查／fpocket 等）")

    return save_stage(session_dir, 20, "receptor_geometry_qc", {
        "status": status,
        "chain_count": n_chain,
        "approx_aa_residues": n_aa,
        "note": "占位：僅計數。深度幾何 QC 待擴充。",
    })


def _update_manifest_engine_c(session_dir: Path, pdbqt_path: Path | None) -> None:
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
        if pdbqt_path and pdbqt_path.is_file():
            rmap["receptor_pdbqt"] = rel_path(session_dir, pdbqt_path)
            paths["receptor_pdbqt_absolute"] = str(pdbqt_path.resolve())
        paths["relative"] = rmap
        ns = dict(paths.get("next_step") or {})
        ns["after_module3_engine_c"] = {
            "receptor_pdbqt": rmap.get("receptor_pdbqt"),
            "hint": "可搭配 Module 2 之 ligand_pdbqt 執行 AutoDock Vina",
        }
        paths["next_step"] = ns
        m["paths"] = paths
        pl = dict(m.get("pipeline") or {})
        pl["module3_engine_c"] = {
            "status": "complete",
            "completed_at": datetime.now().isoformat(),
            "stages": ["19", "20"],
        }
        m["pipeline"] = pl
        summ = dict(m.get("summary") or {})
        if pdbqt_path and pdbqt_path.is_file():
            summ["receptor_pdbqt_relative"] = rmap.get("receptor_pdbqt")
        m["summary"] = summ
        write_manifest(session_dir, m)
    except Exception:
        pass


def run_engine_c(
    session_dir: Path,
    target: str,
    *,
    skip_ready_check: bool = False,
    skip_manifest_target_check: bool = False,
    stage19_engine: str = "auto",
) -> dict[str, Any]:
    print("=" * 55)
    print("  MolGate — Module 3 Engine C（Receptor 準備）")
    print("  Stage 19/20 → protein/receptor.pdbqt")
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

    prep = load_prep_decisions(session_dir)
    if not skip_ready_check and not prep.get("ready_for_engine_e"):
        print("❌ prep_decisions.json 尚未標記 ready_for_engine_e=true，請先完成 Engine B")
        sys.exit(1)

    pt = (prep.get("target") or "").strip().upper()
    if pt and pt != target.strip().upper():
        print(f"❌ prep_decisions.target={prep.get('target')!r} 與 --target={target!r} 不一致")
        sys.exit(1)

    pdb_id = (prep.get("pdb_id") or "").strip().upper()
    if not pdb_id:
        print("❌ prep_decisions 缺少 pdb_id")
        sys.exit(1)

    cleaned_path = load_cleaned_pdb(session_dir)
    print(f"\n  Target：{target.upper()} | PDB：{pdb_id}")
    print(f"  cleaned.pdb：{cleaned_path}")

    try:
        r19 = stage19_receptor(
            session_dir,
            pdb_id,
            cleaned_path,
            engine=stage19_engine,
        )
        r20 = stage20_geometry_qc(session_dir, pdb_id, cleaned_path)
    except Exception as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)

    pdbqt_path = session_dir / "protein" / "receptor.pdbqt"
    if r19.get("status") == "PASS":
        _update_manifest_engine_c(session_dir, pdbqt_path)

    print("\n" + "=" * 55)
    print("  Engine C 完成")
    if r19.get("status") == "PASS":
        print(f"  輸出：{pdbqt_path}")
    else:
        print(
            "  Stage 19 未產出 receptor.pdbqt（可試：obabel 含 pdbqt 編譯、"
            "--engine meeko、--engine mgltools，或手動轉 apo_for_docking.pdb）"
        )
    print("=" * 55)

    if r19.get("status") != "PASS":
        sys.exit(1)

    return {"stage19": r19, "stage20": r20}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="MolGate Module 3 Engine C")
    ap.add_argument("--session-dir", required=True)
    ap.add_argument("--target", required=True)
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
    ap.add_argument(
        "--engine",
        choices=("auto", "obabel", "meeko", "mgltools"),
        default="auto",
        help=(
            "PDB→PDBQT：auto＝meeko → obabel → MGLTools；"
            "obabel 須為含 pdbqt 寫出之編譯；mgltools 需安裝 MGLTools（可設 MGLTOOLS_ROOT）"
        ),
    )
    args = ap.parse_args()

    sd = Path(args.session_dir)
    if not sd.is_dir():
        print(f"❌ Session 目錄不存在：{sd}")
        sys.exit(1)

    run_engine_c(
        sd,
        args.target,
        skip_ready_check=args.skip_ready_check,
        skip_manifest_target_check=args.skip_manifest_target_check,
        stage19_engine=args.engine,
    )
