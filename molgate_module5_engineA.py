"""
MolGate — Module 5: Docking & Output (Engine A)
Stage 24: AutoDock Vina 執行
Stage 25: Pose QC
Stage 26: 結果寫入 docking_result.json

輸入：
  pocket_config.json          ← Module 4 Engine C
  ligand/ligand.pdbqt         ← Module 2 Engine C
  protein/receptor.pdbqt      ← Module 3 Engine C

輸出：
  docking/output.pdbqt        ← Vina 所有 pose
  stage_24_vina.json
  stage_25_qc.json
  stage_26_result.json
  docking_result.json         ← Engine B 契約

依賴：AutoDock Vina 1.2.x（PATH 或 --vina-bin 指定）
契約對齊：pocket_config.json（Module 4 Engine C，含 target / pdb_id / common_name / target_drug / resseq）
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

_PKG_ROOT = Path(__file__).resolve().parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

SCHEMA_VERSION = "0.2"
POCKET_CONFIG_SCHEMA_MIN = "0.2"

# QC 閾值
QC_SCORE_THRESHOLD = -6.0  # kcal/mol，低於此值才算合理（較佳結合）
RMSD_OUT_OF_POCKET_CUTOFF = 8.0  # Å，超過視為偏離主口袋

_FLOAT_RE = re.compile(r"^[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?$")


def _is_float_token(s: str) -> bool:
    return bool(_FLOAT_RE.match(s.strip()))


def _normalize_resname_3(raw: str) -> str:
    """MGLTools 偶發把欄位併成 'A1LYS' 類；取末三碼為殘基名。"""
    r = raw.strip().upper()
    if len(r) >= 3 and r[-3:].isalpha():
        return r[-3:]
    return (r + "XXX")[:3]


def _split_merged_atom_resname(merged: str) -> tuple[str, str] | None:
    """
    MGLTools 常把原子名與三字母殘基併成單一 token（如 HH1A1ARG、NH1AARG）。
    自尾端取三個字母為殘基名，其餘為原子名（供 _format 再截成 4 字元）。
    """
    m = merged.strip().upper()
    if len(m) < 4:
        return None
    r3 = m[-3:]
    if not r3.isalpha():
        return None
    atom = m[:-3].strip()
    if not atom:
        return None
    return atom, r3


def _pdbqt_coord_slices_valid(line: str) -> bool:
    """Vina 依 PDB 固定欄讀 x/y/z（30–54）；欄內不得含空白，否則會報座標解析錯誤。"""
    if len(line) < 54:
        return False
    for a, b in ((30, 38), (38, 46), (46, 54)):
        t = line[a:b].strip()
        if not _is_float_token(t):
            return False
    return True


def _format_pdbqt_atom_line(
    rec: str,
    serial: int,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    occ: float,
    temp_f: float,
    charge: float,
    ad_type: str,
) -> str:
    """輸出 Vina 可讀的 PDBQT 行（欄寬對齊）。rec 為 'ATOM' 或 'HETATM'。"""
    an = (atom_name.strip()[:4]).ljust(4)
    rn = (_normalize_resname_3(resname))[:3].rjust(3)
    ch = (chain_id.strip() or "A")[:1]
    typ = ad_type.strip()
    prefix = (rec + "  ")[:6]
    # PDB：13–16 atom、17 altLoc（空白）、18–20 resName、21 chain、22–26 resSeq
    return (
        f"{prefix}{serial:5d} {an} {rn} {ch}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{temp_f:6.2f}    {charge:8.3f} {typ}\n"
    )


def _parse_pdbqt_atom_line_tokens(line: str) -> dict[str, Any] | None:
    """
    以空白切分還原 ATOM 欄位（處理 MGLTools 把 resname 寫成 A1LYS 等情形）。
    期待尾端：x y z occ temp charge ad_type
    """
    raw = line.rstrip("\n\r")
    if not (raw.startswith("ATOM") or raw.startswith("HETATM")):
        return None
    parts = raw.split()
    if len(parts) < 12:
        return None

    # 尾端錨定：最後 7 token 為 x,y,z,occ,temp,charge,ad_type（支援 12-token 合併欄位）
    tail = parts[-7:]
    if len(tail) == 7 and all(_is_float_token(t) for t in tail[:-1]) and not _is_float_token(tail[-1]):
        x, y, z = float(tail[0]), float(tail[1]), float(tail[2])
        occ, temp_f, charge = float(tail[3]), float(tail[4]), float(tail[5])
        ad_type = tail[6]
        head = parts[:-7]
        if not head or head[0] not in ("ATOM", "HETATM"):
            return None
        if len(head) == 6:
            try:
                serial = int(head[1])
            except ValueError:
                return None
            atom_name = head[2]
            resname = head[3]
            chain_id = head[4]
            try:
                resseq = int(head[5])
            except ValueError:
                return None
            return {
                "serial": serial,
                "atom_name": atom_name,
                "resname": resname,
                "chain_id": chain_id,
                "resseq": resseq,
                "x": x,
                "y": y,
                "z": z,
                "occ": occ,
                "temp_f": temp_f,
                "charge": charge,
                "ad_type": ad_type,
            }
        if len(head) == 5:
            try:
                serial = int(head[1])
            except ValueError:
                return None
            merged = head[2]
            chain_id = head[3]
            try:
                resseq = int(head[4])
            except ValueError:
                return None
            split = _split_merged_atom_resname(merged)
            if split is None:
                return None
            atom_name, resname = split
            return {
                "serial": serial,
                "atom_name": atom_name,
                "resname": resname,
                "chain_id": chain_id,
                "resseq": resseq,
                "x": x,
                "y": y,
                "z": z,
                "occ": occ,
                "temp_f": temp_f,
                "charge": charge,
                "ad_type": ad_type,
            }

    if len(parts) < 13:
        return None
    # 自尾端找連續三個座標
    for i in range(len(parts) - 6, 5, -1):
        if i + 2 >= len(parts):
            continue
        if not (
            _is_float_token(parts[i])
            and _is_float_token(parts[i + 1])
            and _is_float_token(parts[i + 2])
        ):
            continue
        if not _is_float_token(parts[i + 3]) or not _is_float_token(parts[i + 4]):
            continue
        if not _is_float_token(parts[i + 5]):
            continue
        ad_type = parts[-1]
        charge = float(parts[i + 5])
        temp_f = float(parts[i + 4])
        occ = float(parts[i + 3])
        x, y, z = float(parts[i]), float(parts[i + 1]), float(parts[i + 2])
        head = parts[:i]
        if len(head) < 6:
            return None
        try:
            serial = int(head[1])
        except ValueError:
            return None
        atom_name = head[2]
        resname = head[3]
        chain_id = head[4]
        try:
            resseq = int(head[5])
        except ValueError:
            return None
        return {
            "serial": serial,
            "atom_name": atom_name,
            "resname": resname,
            "chain_id": chain_id,
            "resseq": resseq,
            "x": x,
            "y": y,
            "z": z,
            "occ": occ,
            "temp_f": temp_f,
            "charge": charge,
            "ad_type": ad_type,
        }
    return None


def _try_read_pdb_fixed_columns(line: str) -> dict[str, Any] | None:
    """標準 PDB 欄位讀取（座標在 31–54 欄）；成功表示此行已對齊。"""
    if len(line) < 68:
        return None
    try:
        rec = line[0:6].strip()
        if rec not in ("ATOM", "HETATM"):
            return None
        serial = int(line[6:11])
        atom_name = line[12:16]
        resname = line[17:20]
        chain_id = line[21]
        resseq = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        occ = float(line[54:60])
        temp_f = float(line[60:66])
        tail = line[66:].strip().split()
        if len(tail) < 2:
            return None
        charge = float(tail[0])
        ad_type = tail[1]
        return {
            "rec": rec,
            "serial": serial,
            "atom_name": atom_name,
            "resname": resname,
            "chain_id": chain_id,
            "resseq": resseq,
            "x": x,
            "y": y,
            "z": z,
            "occ": occ,
            "temp_f": temp_f,
            "charge": charge,
            "ad_type": ad_type,
        }
    except (ValueError, IndexError):
        return None


def _ad4_type_is_hydrogen(ad_type: str) -> bool:
    t = ad_type.strip().upper()
    if not t:
        return False
    # AD4 氫類型多以 H 開頭；受體中極少見 Hg（若遇到可再擴充白名單）
    if t == "H" or (t[0] == "H" and len(t) == 2):
        return True
    return False


def _sanitize_receptor_pdbqt_for_vina(
    receptor_in: Path,
    receptor_out: Path,
    *,
    drop_hydrogens: bool,
) -> tuple[int, int]:
    """
    重寫 receptor PDBQT，修正欄位錯位。
    回傳 (輸出的 ATOM/HETATM 行數, 重新格式化行數)。
    """
    text = receptor_in.read_text(encoding="utf-8", errors="replace")
    out_lines: list[str] = []
    n_atom = 0
    n_reformatted = 0
    for line in text.splitlines(keepends=True):
        bare = line.rstrip("\n\r")
        if not (bare.startswith("ATOM") or bare.startswith("HETATM")):
            out_lines.append(line if line.endswith("\n") else line + "\n")
            continue
        padded = bare + " " * 120
        col = _try_read_pdb_fixed_columns(padded)
        if col is not None and _pdbqt_coord_slices_valid(padded):
            if drop_hydrogens and _ad4_type_is_hydrogen(col["ad_type"]):
                continue
            out_lines.append(line if line.endswith("\n") else line + "\n")
            n_atom += 1
            continue
        parsed = _parse_pdbqt_atom_line_tokens(bare)
        if parsed is None:
            out_lines.append(line if line.endswith("\n") else line + "\n")
            n_atom += 1
            continue
        if drop_hydrogens and _ad4_type_is_hydrogen(parsed["ad_type"]):
            continue
        rec = "HETATM" if bare.startswith("HETATM") else "ATOM"
        out_lines.append(
            _format_pdbqt_atom_line(
                rec,
                parsed["serial"],
                parsed["atom_name"],
                parsed["resname"],
                parsed["chain_id"],
                parsed["resseq"],
                parsed["x"],
                parsed["y"],
                parsed["z"],
                parsed["occ"],
                parsed["temp_f"],
                parsed["charge"],
                parsed["ad_type"],
            )
        )
        n_atom += 1
        n_reformatted += 1
    receptor_out.parent.mkdir(parents=True, exist_ok=True)
    receptor_out.write_text("".join(out_lines), encoding="utf-8")
    return n_atom, n_reformatted


def _prepare_receptor_for_vina(session_dir: Path, receptor: Path) -> Path:
    """
    產生供 Vina 使用的 receptor：修正 MGLTools 等造成的欄位錯位。
    輸出 protein/receptor_for_vina.pdbqt（至少會複製並掃描；有修正時會提示）。
    """
    fixed = session_dir / "protein" / "receptor_for_vina.pdbqt"
    n_atom, n_ref = _sanitize_receptor_pdbqt_for_vina(receptor, fixed, drop_hydrogens=False)
    if n_atom <= 0:
        return receptor
    if n_ref > 0:
        print(f"  → receptor PDBQT 已正規化供 Vina：{fixed}（修正 {n_ref} 行，共 {n_atom} 行原子）")
    else:
        print(f"  → 使用正規化副本供 Vina：{fixed}（{n_atom} 行，欄位已對齊無需改寫）")
    return fixed


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

def load_pocket_config(session_dir: Path) -> dict[str, Any]:
    path = session_dir / "pocket_config.json"
    if not path.exists():
        raise FileNotFoundError("找不到 pocket_config.json，請先跑 Module 4 Engine C")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not data.get("ready_for_docking"):
        raise ValueError("pocket_config.json 未標記 ready_for_docking=True，請重跑 Module 4")
    sv = (data.get("schema_version") or "").strip()
    if sv and sv != POCKET_CONFIG_SCHEMA_MIN:
        print(
            f"  ⚠ pocket_config.schema_version={sv!r}，建議與 Module 4 Engine C（{POCKET_CONFIG_SCHEMA_MIN}）一致後再跑"
        )
    return data


def load_prep_decisions(session_dir: Path) -> dict[str, Any] | None:
    p = session_dir / "prep_decisions.json"
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None


def resolve_and_validate_target(
    session_dir: Path,
    pocket: dict[str, Any],
    cli_target: str | None,
    *,
    skip_manifest_target_check: bool = False,
    skip_prep_target_check: bool = False,
) -> str:
    """
    標靶：CLI --target 優先，否則 pocket_config.target；
    與 prep_decisions.target、molgate_run_manifest.json 交叉驗證（可 skip）。
    """
    pt_pocket = (pocket.get("target") or "").strip().upper()
    cli = (cli_target or "").strip().upper()
    target = cli or pt_pocket
    if not target:
        raise ValueError(
            "缺少 target：請使用 --target，或確認 pocket_config.json 已寫入 target（Module 4 Engine C）"
        )
    if pt_pocket and cli and pt_pocket != cli:
        raise ValueError(
            f"pocket_config.target={pocket.get('target')!r} 與 --target={cli_target!r} 不一致"
        )

    if not skip_prep_target_check:
        prep = load_prep_decisions(session_dir)
        if prep:
            pt = (prep.get("target") or "").strip().upper()
            if pt and pt != target:
                raise ValueError(
                    f"prep_decisions.target={prep.get('target')!r} 與標靶 {target!r} 不一致"
                )

    if not skip_manifest_target_check:
        try:
            from molgate_manifest import validate_target_matches_manifest

            validate_target_matches_manifest(session_dir, target)
        except FileNotFoundError:
            pass
        except ValueError:
            raise

    return target


def _update_manifest_module5_engine_a(
    session_dir: Path,
    docking_result_path: Path,
    output_pdbqt: Path,
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
        paths = dict(m.get("paths") or {})
        rmap = dict(paths.get("relative") or {})
        rmap["docking_result"] = rel_path(session_dir, docking_result_path)
        rmap["docking_output_pdbqt"] = rel_path(session_dir, output_pdbqt)
        paths["relative"] = rmap
        paths["docking_result_absolute"] = str(docking_result_path.resolve())
        paths["docking_output_pdbqt_absolute"] = str(output_pdbqt.resolve())
        ns = dict(paths.get("next_step") or {})
        ns["after_module5_engine_a"] = {
            "read_docking_result_from": rmap["docking_result"],
            "hint": "Engine B：審核 docking_result.json",
        }
        paths["next_step"] = ns
        m["paths"] = paths
        pl = dict(m.get("pipeline") or {})
        pl["module5_engine_a"] = {
            "status": "complete",
            "completed_at": datetime.now().isoformat(),
            "stages": ["24", "25", "26"],
        }
        m["pipeline"] = pl
        write_manifest(session_dir, m)
    except Exception:
        pass

def find_vina(vina_bin: str | None) -> str:
    """找 AutoDock Vina 執行檔。"""
    if vina_bin:
        p = Path(vina_bin)
        if p.is_file():
            return str(p)
        raise FileNotFoundError(f"找不到指定的 Vina：{vina_bin}")
    for name in ("vina", "vina.exe", "autodock_vina", "autodock_vina.exe"):
        found = shutil.which(name)
        if found:
            return found
    raise FileNotFoundError(
        "找不到 AutoDock Vina，請安裝後加入 PATH，或用 --vina-bin 指定路徑"
    )

def parse_vina_output(output_text: str) -> list[dict[str, Any]]:
    """解析 Vina stdout 的 pose 分數列表。"""
    poses = []
    for line in output_text.splitlines():
        # Vina 1.2.x 輸出格式：
        # mode |   affinity | dist from best mode
        #      | (kcal/mol) | rmsd l.b.| rmsd u.b.
        #   1   -8.612       0.000      0.000
        m = re.match(r'\s*(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)', line)
        if m:
            poses.append({
                "mode":     int(m.group(1)),
                "affinity": float(m.group(2)),
                "rmsd_lb":  float(m.group(3)),
                "rmsd_ub":  float(m.group(4)),
            })
    return poses

def run_cmd(argv: list[str], timeout: int = 600) -> tuple[int, str, str]:
    try:
        r = subprocess.run(
            argv,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        return r.returncode, r.stdout or "", r.stderr or ""
    except subprocess.TimeoutExpired:
        return -1, "", "Timeout"
    except OSError as e:
        return -1, "", str(e)

# ── Stage 24：AutoDock Vina 執行 ──────────────────────
def stage24_vina(
    session_dir: Path,
    pocket: dict[str, Any],
    vina_bin: str,
    exhaustiveness: int = 8,
    num_modes: int = 9,
) -> dict[str, Any]:
    print("\n[Stage 24] AutoDock Vina 執行...")

    receptor = session_dir / "protein" / "receptor.pdbqt"
    ligand   = session_dir / "ligand"  / "ligand.pdbqt"

    if not receptor.is_file():
        raise FileNotFoundError(f"找不到 receptor.pdbqt：{receptor}")
    if not ligand.is_file():
        raise FileNotFoundError(f"找不到 ligand.pdbqt：{ligand}")

    out_dir = session_dir / "docking"
    out_dir.mkdir(exist_ok=True)
    out_pdbqt = out_dir / "output.pdbqt"
    log_file  = out_dir / "vina.log"

    receptor_use = _prepare_receptor_for_vina(session_dir, receptor)

    print(f"  → Receptor：{receptor_use}")
    print(f"  → Ligand：{ligand}")
    print(f"  → Center：({pocket['center_x']}, {pocket['center_y']}, {pocket['center_z']})")
    print(f"  → Box：{pocket['box_size_x']} × {pocket['box_size_y']} × {pocket['box_size_z']} Å")
    print(f"  → Exhaustiveness：{exhaustiveness} | Modes：{num_modes}")
    print(f"  → 執行中（請稍候）...")

    def _run_vina(rec_path: Path) -> tuple[int, str, str, str]:
        argv = [
            vina_bin,
            "--receptor",    str(rec_path),
            "--ligand",      str(ligand),
            "--center_x",    str(pocket["center_x"]),
            "--center_y",    str(pocket["center_y"]),
            "--center_z",    str(pocket["center_z"]),
            "--size_x",      str(pocket["box_size_x"]),
            "--size_y",      str(pocket["box_size_y"]),
            "--size_z",      str(pocket["box_size_z"]),
            "--out",         str(out_pdbqt),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes",   str(num_modes),
        ]
        code, stdout, stderr = run_cmd(argv)
        combined = stdout + "\n" + stderr
        return code, stdout, stderr, combined

    code, stdout, stderr, combined = _run_vina(receptor_use)
    poses = parse_vina_output(combined)

    pdbqt_parse_fail = (
        code != 0
        and not poses
        and (
            "PDBQT parsing" in combined
            or "parsing error" in combined.lower()
            or "not valid" in combined.lower()
        )
    )
    if pdbqt_parse_fail:
        noh_path = session_dir / "protein" / "receptor_for_vina_noh.pdbqt"
        n_kept, n_ref = _sanitize_receptor_pdbqt_for_vina(
            receptor_use, noh_path, drop_hydrogens=True
        )
        print(
            f"  ⚠ Vina 無法解析 receptor PDBQT（常見於 MGLTools 欄位錯位）；"
            f"已剔除 AD4 氫原子列後重試（保留 {n_kept} 行）"
        )
        receptor_use = noh_path
        code, stdout, stderr, combined = _run_vina(receptor_use)
        poses = parse_vina_output(combined)

    # 寫 log
    log_file.write_text(combined, encoding="utf-8")

    # 若 stdout 沒有，嘗試從 log 再解析
    if not poses:
        poses = parse_vina_output(combined)

    if code != 0 and not poses:
        result = {
            "status": "FAIL",
            "reason": f"Vina 執行失敗（exit {code}）",
            "stderr": stderr[:500],
        }
        save_stage(session_dir, 24, "vina", result)
        raise ValueError(f"Stage 24 FAIL: {result['reason']}\n{stderr[:200]}")

    print(f"  ✓ Vina 完成，找到 {len(poses)} 個 pose")
    for p in poses:
        print(f"    Mode {p['mode']:2d}: {p['affinity']:7.3f} kcal/mol "
              f"| RMSD lb={p['rmsd_lb']:.2f} ub={p['rmsd_ub']:.2f}")

    result = {
        "status": "PASS",
        "vina_bin": vina_bin,
        "exhaustiveness": exhaustiveness,
        "num_modes": num_modes,
        "pose_count": len(poses),
        "poses": poses,
        "output_pdbqt": str(out_pdbqt.as_posix()),
        "log_file": str(log_file.as_posix()),
        "pocket_het_id": pocket.get("het_id", ""),
        "pocket_chain":  pocket.get("chain", ""),
    }
    return save_stage(session_dir, 24, "vina", result)

# ── Stage 25：Pose QC ─────────────────────────────────
def stage25_qc(
    session_dir: Path,
    r24: dict[str, Any],
    pocket: dict[str, Any],
) -> dict[str, Any]:
    print("\n[Stage 25] Pose QC...")

    poses = r24.get("poses", [])
    if not poses:
        result = {
            "status": "FAIL",
            "reason": "Stage 24 無 pose 輸出",
        }
        save_stage(session_dir, 25, "qc", result)
        raise ValueError("Stage 25 FAIL: 無 pose 可評估")

    best = poses[0]
    score = best["affinity"]

    # QC 判斷
    flags = []
    qc_pass = True

    # 分數檢查
    if score > QC_SCORE_THRESHOLD:
        flags.append(f"SCORE_WEAK: best pose {score:.3f} kcal/mol > 閾值 {QC_SCORE_THRESHOLD}")
        qc_pass = False
        print(f"  ⚠ 分數偏弱：{score:.3f} kcal/mol（閾值 {QC_SCORE_THRESHOLD}）")
    else:
        print(f"  ✓ 分數合理：{score:.3f} kcal/mol")

    # RMSD 檢查（若多個 pose 都一樣，可能有問題）
    if len(poses) > 1:
        rmsds = [p["rmsd_lb"] for p in poses[1:]]
        if all(r < 0.5 for r in rmsds):
            flags.append("RMSD_COLLAPSED: 所有 pose RMSD 相近，搜索空間可能太小")
            print(f"  ⚠ 所有 pose RMSD 相近，建議增大 box size 或 exhaustiveness")

    # RMSD lb 偏大：僅當口袋外 pose 佔多數才降級
    out_of_pocket_modes = [p["mode"] for p in poses if float(p["rmsd_lb"]) > RMSD_OUT_OF_POCKET_CUTOFF]
    out_of_pocket_count = len(out_of_pocket_modes)
    total_poses = len(poses)
    if out_of_pocket_count > 0 and total_poses > 0:
        if (out_of_pocket_count / total_poses) > 0.5:
            flags.append("OUT_OF_POCKET_MAJORITY")
            qc_pass = False
            print(
                f"  ⚠ OUT_OF_POCKET_MAJORITY：{out_of_pocket_count}/{total_poses} "
                f"poses RMSD lb > {RMSD_OUT_OF_POCKET_CUTOFF:.1f} Å"
            )
        else:
            flags.append(f"OUT_OF_POCKET_MINOR: {out_of_pocket_count}/{total_poses} poses outside")
            print(
                f"  ℹ OUT_OF_POCKET_MINOR：{out_of_pocket_count}/{total_poses} "
                f"poses RMSD lb > {RMSD_OUT_OF_POCKET_CUTOFF:.1f} Å（不降級）"
            )

    status = "PASS" if qc_pass else "WARN"
    if qc_pass:
        print(f"  ✓ QC 通過")
    else:
        print(f"  ⚠ QC 有警告，建議確認結果")

    print(f"  → Status: {status}")

    result = {
        "status": status,
        "best_pose_mode": best["mode"],
        "best_affinity": score,
        "best_rmsd_lb": best["rmsd_lb"],
        "best_rmsd_ub": best["rmsd_ub"],
        "score_threshold": QC_SCORE_THRESHOLD,
        "out_of_pocket_rmsd_lb_cutoff": RMSD_OUT_OF_POCKET_CUTOFF,
        "out_of_pocket_modes": out_of_pocket_modes,
        "qc_pass": qc_pass,
        "flags": flags,
        "all_poses": poses,
    }
    return save_stage(session_dir, 25, "qc", result)

# ── Stage 26：結果寫入 ────────────────────────────────
def stage26_result(
    session_dir: Path,
    r24: dict[str, Any],
    r25: dict[str, Any],
    pocket: dict[str, Any],
    target: str,
) -> dict[str, Any]:
    print("\n[Stage 26] 結果寫入...")

    flags_all = r25.get("flags", [])

    # confidence tier
    if r25["qc_pass"] and not flags_all:
        tier = "A"
    elif r25["qc_pass"]:
        tier = "B"
    else:
        tier = "C"

    result: dict[str, Any] = {
        "status": "PASS",
        "target": target.upper(),
        "pdb_id": pocket.get("pdb_id", ""),
        "het_id": pocket.get("het_id", ""),
        "chain": pocket.get("chain", ""),
        "resseq": pocket.get("resseq"),
        "pocket_rule": pocket.get("pocket_rule", ""),
        "best_affinity": r25["best_affinity"],
        "best_pose_mode": r25["best_pose_mode"],
        "confidence_tier": tier,
        "qc_pass": r25["qc_pass"],
        "flags": flags_all,
        "all_poses": r25["all_poses"],
        "output_pdbqt": r24["output_pdbqt"],
        "pocket_center": {
            "x": pocket["center_x"],
            "y": pocket["center_y"],
            "z": pocket["center_z"],
        },
        "ready_for_review": True,
    }
    cn = (pocket.get("common_name") or "").strip()
    td = (pocket.get("target_drug") or "").strip()
    if cn:
        result["common_name"] = cn
    if td:
        result["target_drug"] = td

    r = save_stage(session_dir, 26, "result", result)

    # 寫 docking_result.json（Engine B 契約）
    docking_result_path = session_dir / "docking_result.json"
    docking_result: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "target": result["target"],
        "pdb_id": result["pdb_id"],
        "het_id": result["het_id"],
        "chain": result.get("chain", ""),
        "resseq": result.get("resseq"),
        "pocket_rule": result.get("pocket_rule", ""),
        "best_affinity": result["best_affinity"],
        "best_pose_mode": result["best_pose_mode"],
        "confidence_tier": tier,
        "qc_pass": result["qc_pass"],
        "flags": flags_all,
        "all_poses": result["all_poses"],
        "output_pdbqt": result["output_pdbqt"],
        "pocket_center": result["pocket_center"],
        "ready_for_review": True,
        "_timestamp": datetime.now().isoformat(),
    }
    if "common_name" in result:
        docking_result["common_name"] = result["common_name"]
    if "target_drug" in result:
        docking_result["target_drug"] = result["target_drug"]
    docking_result_path.write_text(
        json.dumps(docking_result, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(f"  ✓ Confidence tier：{tier}")
    print(f"  ✓ docking_result.json 寫入完成")
    print(f"  → Status: PASS")

    return r

# ── 主流程 ────────────────────────────────────────────
def run_engine_a(
    session_dir: Path,
    target: str | None,
    vina_bin: str | None = None,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    *,
    skip_manifest_target_check: bool = False,
    skip_prep_target_check: bool = False,
):
    print("=" * 55)
    print("  MolGate — Module 5 Engine A")
    print("  Stage 24/25/26：Vina + QC + 結果寫入")
    print("=" * 55)

    try:
        pocket = load_pocket_config(session_dir)
    except (FileNotFoundError, ValueError) as e:
        print(f"❌ {e}")
        sys.exit(1)

    try:
        resolved_target = resolve_and_validate_target(
            session_dir,
            pocket,
            target,
            skip_manifest_target_check=skip_manifest_target_check,
            skip_prep_target_check=skip_prep_target_check,
        )
    except ValueError as e:
        print(f"❌ {e}")
        sys.exit(1)

    try:
        vina = find_vina(vina_bin)
        print(f"\n  Vina：{vina}")
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    td = (pocket.get("target_drug") or "").strip()
    drug_note = f" | target_drug：{td}" if td else ""
    print(
        f"  Target：{resolved_target} | PDB：{pocket.get('pdb_id', '?')}"
        f"{drug_note}"
    )

    try:
        r24 = stage24_vina(session_dir, pocket, vina, exhaustiveness, num_modes)
        r25 = stage25_qc(session_dir, r24, pocket)
        r26 = stage26_result(session_dir, r24, r25, pocket, resolved_target)
    except ValueError as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Pipeline stopped: {e}")
        sys.exit(1)

    dr_path = session_dir / "docking_result.json"
    out_pdbqt = Path(r24["output_pdbqt"])
    _update_manifest_module5_engine_a(session_dir, dr_path, out_pdbqt)

    print("\n" + "=" * 55)
    print("  Engine A 完成 ✓")
    print(f"  Best pose：Mode {r25['best_pose_mode']} | "
          f"{r25['best_affinity']:.3f} kcal/mol")
    print(f"  Confidence tier：{r26['confidence_tier']}")
    if r25.get("flags"):
        print(f"  Flags：{' | '.join(r25['flags'])}")
    print(f"  契約：{session_dir / 'docking_result.json'}")
    print(f"  ready_for_review = True")
    print("=" * 55)

    return {"stage24": r24, "stage25": r25, "stage26": r26}

# ── CLI ──────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="MolGate Module 5 Engine A")
    ap.add_argument("--session-dir", required=True)
    ap.add_argument(
        "--target",
        default=None,
        help="標靶 ID（建議指定；若省略則使用 pocket_config.json 的 target）",
    )
    ap.add_argument(
        "--vina-bin",
        default=None,
        help="AutoDock Vina 執行檔路徑（預設從 PATH 找）",
    )
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument(
        "--skip-manifest-target-check",
        action="store_true",
        help="不與 molgate_run_manifest.json 的 target 交叉驗證",
    )
    ap.add_argument(
        "--skip-prep-target-check",
        action="store_true",
        help="不檢查 prep_decisions.target 與標靶一致",
    )
    args = ap.parse_args()

    sd = Path(args.session_dir)
    if not sd.is_dir():
        print(f"❌ Session 目錄不存在：{sd}")
        sys.exit(1)

    run_engine_a(
        sd,
        args.target,
        vina_bin=args.vina_bin,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        skip_manifest_target_check=args.skip_manifest_target_check,
        skip_prep_target_check=args.skip_prep_target_check,
    )
