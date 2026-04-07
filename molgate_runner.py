"""
MolGate — Single Drug Runner
一次跑一個藥的完整 pipeline（Module 1 → Module 5）
人機協作正常進行，遇到決策點自然停下來

用法：
  python molgate_runner.py
  python molgate_runner.py --drug 0        # 直接選第 0 個藥（目前 0–11）
  python molgate_runner.py --from-stage 3  # 從 Module 3 繼續（已有 session）
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

PYTHON = sys.executable
BASE_DIR = Path(__file__).resolve().parent

if hasattr(sys.stdout, "reconfigure"):
    # Avoid Windows console codec crashes on non-encodable CJK chars.
    sys.stdout.reconfigure(errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(errors="replace")

# ── 藥物清單 ──────────────────────────────────────────
DRUGS = [
    {
        "drug":    "Ibuprofen",
        "index_key": "COX2",
        "smiles":  "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "target":  "COX2",
        "pdb_id":  "5IKT",
        "het_id":  "TLF",
        "note":    "COX-2 抑制劑，NSAID 代表藥物",
    },
    {
        "drug":    "Celecoxib",
        "index_key": "COX2_CLX",
        "smiles":  "CC1=CC=C(C=C1)C1=CC(=NN1C1=CC=C(C=C1)S(N)(=O)=O)C(F)(F)F",
        "target":  "COX2",
        "pdb_id":  "3LN1",
        "het_id":  "CLX",
        "note":    "COX-2 選擇性抑制劑",
    },
    {
        "drug":    "Oseltamivir",
        "index_key": "NEURAMINIDASE",
        "smiles":  "CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@@H](NC(C)=O)[C@@H](N)C1",
        "target":  "Neuraminidase",
        "pdb_id":  "2HU4",
        "het_id":  "G39",
        "note":    "流感神經氨酸酶抑制劑（克流感）",
    },
    {
        "drug":    "Imatinib",
        "index_key": "ABLABL",
        "smiles":  "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        "target":  "ABLABL",
        "pdb_id":  "1IEP",
        "het_id":  "STI",
        "note":    "BCR-ABL 激酶抑制劑（基利克）",
    },
    {
        "drug":    "Methotrexate",
        "index_key": "DHFR",
        "smiles":  "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1",
        "target":  "DHFR",
        "pdb_id":  "1U72",
        "het_id":  "MTX",
        "note":    "二氫葉酸還原酶抑制劑，抗癌/抗類風溼",
    },
    {
        "drug":    "Fluconazole",
        "index_key": "CYP51",
        "smiles":  "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F",
        "target":  "CYP51",
        "pdb_id":  "4WMZ",
        "het_id":  "FLC",
        "note":    "CYP51 抑制劑，抗黴菌藥",
    },
    {
        "drug":    "Sildenafil",
        "index_key": "PDE5",
        "smiles":  "CCCC1=NN(C)C(=O)c2[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc21",
        "target":  "PDE5",
        "pdb_id":  "2H42",
        "het_id":  "SLD",
        "note":    "PDE5 抑制劑（威而鋼）",
    },
    {
        "drug":    "Acetazolamide",
        "index_key": "CA2",
        "smiles":  "CC(=O)Nc1nnc(S(N)(=O)=O)s1",
        "target":  "CA2",
        "pdb_id":  "3HS4",
        "het_id":  "AZM",
        "note":    "碳酸酐酶 II 抑制劑，利尿劑",
    },
    {
        "drug":    "Gefitinib",
        "index_key": "EGFR",
        "smiles":  "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
        "target":  "EGFR",
        "pdb_id":  "4WKQ",
        "het_id":  "IRE",
        "note":    "EGFR 激酶抑制劑（艾瑞莎）",
    },
    {
        "drug":    "Warfarin",
        "index_key": "VKOR",
        "smiles":  "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
        "target":  "VKOR",
        "pdb_id":  "6WV3",
        "het_id":  "SWF",
        "note":    "人源 VKOR + warfarin（6WV3，共晶配體 CCD：SWF）；勿用 4KP6（該為 PDE4B / 1S1）",
    },
    {
        "drug":    "Tacrolimus",
        "index_key": "FKBP12_TACROLIMUS",
        "smiles":  r"C[C@@H]1C[C@@H]([C@@H]2[C@H](C[C@H]([C@@](O2)(C(=O)C(=O)N3CCCC[C@H]3C(=O)O[C@@H]([C@@H]([C@H](CC(=O)[C@@H](/C=C(/C1)\C)CC=C)O)C)/C(=C/[C@@H]4CC[C@H]([C@@H](C4)OC)O)/C)O)C)OC)OC",
        "target":  "FKBP12",
        "pdb_id":  "1FKJ",
        "het_id":  "FK5",
        "note":    "FKBP12–FK506（1FKJ，CCD FK5）；direct tacrolimus anchor，適合大環／立體佔據基準",
    },
    {
        "drug":    "Quetiapine",
        "index_key": "DRD2_QUETIAPINE_PROXY",
        "smiles":  "C1CN(CCN1CCOCCO)C2=NC3=CC=CC=C3SC4=CC=CC=C42",
        "target":  "DRD2",
        "pdb_id":  "6CM4",
        "het_id":  "8NU",
        "note":    "Proxy：無 quetiapine 直接受體共晶；以 DRD2/6CM4 口袋（共晶 8NU）當暫代，勿當 direct benchmark",
    },
]

_SESSION_ID_RE = re.compile(r"^[a-fA-F0-9]{16}$")

# ── 工具函式 ──────────────────────────────────────────
def run(cmd: list[str]) -> int:
    """執行指令，繼承 stdin/stdout（保留人機介入）。"""
    print(f"\n  $ {' '.join(cmd)}")
    print("  " + "─" * 50)
    env = dict(os.environ)
    env["PYTHONUTF8"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    result = subprocess.run(cmd, env=env)
    return result.returncode

def banner(msg: str):
    print("\n" + "═" * 55)
    print(f"  {msg}")
    print("═" * 55)


def print_startup_banner() -> None:
    """啟動儀式橫幅（ASCII；SESSION START 為實際執行時間）。"""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(
        f"""
================================================================================
   __  __       _  _____       _       
  |  \/  | ___ | |/ ____| __ _| |_ ___ 
  | |\/| |/ _ \\| | |  __ / _` | __/ _ \\
  | |  | | (_) | | | |_ | (_| | ||  __/
  |_|  |_|\\___/|_|\\____|\\__,_|\\__\\___/

               S I N G L E   D R U G   R U N N E R
--------------------------------------------------------------------------------
  [BOOT] Initializing molecular workflow...
  [MODE] Structured docking pipeline
  [CORE] Input -> Ligand -> Receptor -> Pocket -> Docking -> Review
================================================================================

  SESSION START : {ts}
  ENGINE        : MolGate
  PIPELINE      : Single Drug Runner

================================================================================
"""
    )


def target_key(drug: dict) -> str:
    """Master Index 查詢鍵（與顯示 target 分離）。"""
    return str(drug.get("index_key") or drug["target"]).strip().upper()


def load_master_index() -> dict[str, dict]:
    p = BASE_DIR / "master_index.json"
    if not p.is_file():
        return {}
    try:
        raw = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}
    out: dict[str, dict] = {}
    for k, v in raw.items():
        if k.startswith("_"):
            continue
        if isinstance(v, dict):
            out[k.strip().upper()] = v
    return out


def validate_drug_index_entry(drug: dict) -> tuple[bool, str]:
    """跑前檢查：index key 是否存在，且 pdb/het 與 runner 設定是否一致。"""
    key = target_key(drug)
    idx = load_master_index()
    ent = idx.get(key)
    if ent is None:
        return False, f"Master Index 缺少 key={key}（請先補 master_index.json）"

    mism = []
    p_runner = str(drug.get("pdb_id", "")).strip().upper()
    p_idx = str(ent.get("pdb_id", "")).strip().upper()
    if p_runner and p_idx and p_runner != p_idx:
        mism.append(f"pdb_id runner={p_runner} vs index={p_idx}")

    h_runner = str(drug.get("het_id", "")).strip().upper()
    h_idx = str(ent.get("het_id") or ent.get("cocrystal_het") or "").strip().upper()
    if h_runner and h_idx and h_runner != h_idx:
        mism.append(f"het_id runner={h_runner} vs index={h_idx}")

    if mism:
        return False, "；".join(mism)
    return True, "OK"

def select_drug(default: int | None) -> tuple[int, dict]:
    if default is not None:
        d = DRUGS[default]
        print(f"\n  選定藥物：[{default}] {d['drug']} / {d['target']} (index_key={target_key(d)})")
        return default, d

    print("\n  ── 選擇藥物 ──")
    for i, d in enumerate(DRUGS):
        print(f"  [{i:2d}] {d['drug']:<18} → {d['target']:<12} [{target_key(d):<14}] {d['note']}")

    while True:
        try:
            raw = input("\n  請輸入編號：").strip()
            idx = int(raw)
            if 0 <= idx < len(DRUGS):
                return idx, DRUGS[idx]
        except (ValueError, KeyboardInterrupt):
            pass
        print("  請輸入有效編號")

def get_session_id(session_dir: Path) -> str | None:
    """從 session 目錄找最新的 session ID。"""
    if not session_dir.exists():
        return None
    sessions = sorted(session_dir.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
    if sessions:
        return sessions[0].name
    return None


def list_valid_sessions(session_base: Path) -> list[Path]:
    """只列出合法 session 目錄（16 位 hex 且含 session_meta.json）。"""
    if not session_base.exists():
        return []
    out: list[Path] = []
    for p in session_base.iterdir():
        if not p.is_dir():
            continue
        if not _SESSION_ID_RE.match(p.name):
            continue
        if not (p / "session_meta.json").is_file():
            continue
        out.append(p)
    out.sort(key=lambda x: x.stat().st_mtime, reverse=True)
    return out

def find_latest_session(drug: dict) -> Path | None:
    """找這個藥最近的 session。"""
    session_base = BASE_DIR / "molgate_sessions"
    if not session_base.exists():
        return None
    # 找有 stage_03_index.json（Module 1 完成）且 index key 符合的 session
    key = target_key(drug)
    for s in list_valid_sessions(session_base):
        idx = s / "stage_03_index.json"
        if idx.exists():
            try:
                data = json.loads(idx.read_text(encoding="utf-8"))
                if data.get("target", "").upper() == key:
                    return s
            except Exception:
                pass
    return None

# ── Pipeline 各 Module ────────────────────────────────
def run_module1(drug: dict) -> Path | None:
    banner("MODULE 1 — Input Validation")
    code = run([
        PYTHON, str(BASE_DIR / "molgate_module1.py"),
        "--smiles",  drug["smiles"],
        "--target",  target_key(drug),
        "--pdb-id",  drug["pdb_id"],
        "--allow-network",
    ])
    if code != 0:
        print("❌ Module 1 失敗")
        return None

    # 找剛建立且與本次 target/pdb 對應的 session（避免誤撿其他藥物最新 session）
    session_base = BASE_DIR / "molgate_sessions"
    key = target_key(drug)
    expect_pdb = str(drug.get("pdb_id", "")).strip().upper()
    for s in list_valid_sessions(session_base):
        idx = s / "stage_03_index.json"
        if not idx.is_file():
            continue
        try:
            data = json.loads(idx.read_text(encoding="utf-8"))
        except Exception:
            continue
        t = str(data.get("target", "")).strip().upper()
        p = str(data.get("pdb_id", "")).strip().upper()
        if t == key and (not expect_pdb or p == expect_pdb):
            return s
    return None

def run_module2(sd: Path, drug: dict) -> bool:
    banner("MODULE 2 — Ligand Enumeration")
    for script, label in [
        ("molgate_module2_engineB.py", "Engine B：質子化 + 互變異構"),
        ("molgate_module2_engineC.py", "Engine C：PDBQT 轉換"),
    ]:
        print(f"\n  [{label}]")
        code = run([
            PYTHON, str(BASE_DIR / script),
            "--session-dir", str(sd),
            "--target", target_key(drug),
        ])
        if code != 0:
            print(f"❌ Module 2 {script} 失敗")
            return False
    return True

def run_module3(sd: Path, drug: dict) -> bool:
    banner("MODULE 3 — Protein Preparation")
    for script, label in [
        ("molgate_module3_engineA.py", "Engine A：PDB 清理"),
        ("molgate_module3_engineB.py", "Engine B：人機決策"),
        ("molgate_module3_engineC.py", "Engine C：Receptor PDBQT"),
    ]:
        print(f"\n  [{label}]")
        code = run([
            PYTHON, str(BASE_DIR / script),
            "--session-dir", str(sd),
            "--target", target_key(drug),
        ])
        if code != 0:
            print(f"❌ Module 3 {script} 失敗")
            return False
    return True

def run_module4(sd: Path, drug: dict) -> bool:
    banner("MODULE 4 — Pocket Definition")
    for script, label in [
        ("molgate_module4_engineA.py", "Engine A：HETATM 掃描"),
        ("molgate_module4_engineB.py", "Engine B：配體確認"),
        ("molgate_module4_engineC.py", "Engine C：口袋中心計算"),
    ]:
        print(f"\n  [{label}]")
        code = run([
            PYTHON, str(BASE_DIR / script),
            "--session-dir", str(sd),
            "--target", target_key(drug),
        ])
        if code != 0:
            print(f"❌ Module 4 {script} 失敗")
            return False
    return True

def run_module5(sd: Path, drug: dict, vina_bin: str | None) -> bool:
    banner("MODULE 5 — Docking & Review")

    resolved_vina_bin = vina_bin
    if not resolved_vina_bin:
        for cand in (
            BASE_DIR / "tool" / "vina" / "vina.exe",
            BASE_DIR / "tool" / "vina" / "vina.exe.exe",
        ):
            if cand.is_file():
                resolved_vina_bin = str(cand)
                break

    cmd_a = [
        PYTHON, str(BASE_DIR / "molgate_module5_engineA.py"),
        "--session-dir", str(sd),
        "--target", target_key(drug),
    ]
    if resolved_vina_bin:
        cmd_a += ["--vina-bin", resolved_vina_bin]
        print(f"  [Vina binary] {resolved_vina_bin}")

    print("\n  [Engine A：Vina + QC]")
    code = run(cmd_a)
    if code != 0:
        print("❌ Module 5 Engine A 失敗")
        return False

    print("\n  [Engine B：人工審核]")
    code = run([
        PYTHON, str(BASE_DIR / "molgate_module5_engineB.py"),
        "--session-dir", str(sd),
    ])
    if code != 0:
        print("❌ Module 5 Engine B 失敗")
        return False

    return True

# ── 主流程 ────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="MolGate Single Drug Runner")
    ap.add_argument("--drug",       type=int, default=None, help="藥物編號（0-11）")
    ap.add_argument("--from-stage", type=int, default=1,    help="從哪個 Module 開始（1-5）")
    ap.add_argument("--session-dir",default=None,           help="指定 session 目錄（from-stage > 1 時使用）")
    ap.add_argument("--vina-bin",   default=None,           help="AutoDock Vina 執行檔路徑")
    args = ap.parse_args()

    print_startup_banner()

    # 選藥
    drug_idx, drug = select_drug(args.drug)
    print(f"\n  藥物：{drug['drug']}")
    print(
        f"  Target：{drug['target']} | index_key：{target_key(drug)} "
        f"| PDB：{drug['pdb_id']} | HET：{drug['het_id']}"
    )
    print(f"  說明：{drug['note']}")

    ok, msg = validate_drug_index_entry(drug)
    if not ok:
        print(f"❌ Index 前置檢查失敗：{msg}")
        print("   請先補齊/修正 master_index.json 後再跑")
        sys.exit(1)

    from_stage = args.from_stage

    # 確認 session
    if from_stage > 1:
        if args.session_dir:
            sd = Path(args.session_dir)
        else:
            sd = find_latest_session(drug)
            if sd is None:
                print("❌ 找不到此藥的既有 session，請從 --from-stage 1 開始")
                sys.exit(1)
        print(f"\n  繼續 session：{sd}")
    else:
        sd = None

    print(f"\n  從 Module {from_stage} 開始\n")
    input("  按 Enter 開始...")

    # ── Module 1 ──
    if from_stage <= 1:
        sd = run_module1(drug)
        if sd is None:
            sys.exit(1)
        print(f"\n  Session 建立：{sd}")

    # ── Module 2 ──
    if from_stage <= 2:
        if not run_module2(sd, drug):
            sys.exit(1)

    # ── Module 3 ──
    if from_stage <= 3:
        if not run_module3(sd, drug):
            sys.exit(1)

    # ── Module 4 ──
    if from_stage <= 4:
        if not run_module4(sd, drug):
            sys.exit(1)

    # ── Module 5 ──
    if from_stage <= 5:
        if not run_module5(sd, drug, args.vina_bin):
            sys.exit(1)

    # ── 完成 ──
    banner("Pipeline 完成 ✓")
    print(f"  藥物：{drug['drug']} / {drug['target']}")
    print(f"  Session：{sd}")
    report = sd / "final_report.json"
    if report.exists():
        r = json.loads(report.read_text(encoding="utf-8"))
        print(f"  Best affinity：{r.get('best_affinity','?')} kcal/mol")
        print(f"  Confidence tier：{r.get('confidence_tier','?')}")
        print(f"  Decision：{r.get('status','?')}")
    print(f"\n  結束時間：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    main()
