"""
MolGate — Module 4 本機面板後端（方向 2）
- 提供靜態 molgate_module4_panel.html
- GET  /api/module4/status?session=…  讀 stage_21/22/23、pocket_config
- POST /api/module4/run-engine-a     subprocess 執行 Engine A（B/C 仍用終端）

使用方式（在 MolGate 目錄）：
  py -3 molgate_module4_server.py
瀏覽器：http://127.0.0.1:8766/?session=你的16位hex&target=COX2
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import urllib.parse
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from typing import Any

_PKG_ROOT = Path(__file__).resolve().parent
_SESSIONS = _PKG_ROOT / "molgate_sessions"
_PANEL = _PKG_ROOT / "molgate_module4_panel.html"

_SESSION_ID_RE = re.compile(r"^[a-fA-F0-9]{16}$")


def _read_json(p: Path) -> dict[str, Any] | None:
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None


def _safe_session_dir(session_id: str) -> Path | None:
    sid = (session_id or "").strip()
    if not sid or not _SESSION_ID_RE.match(sid):
        return None
    p = (_SESSIONS / sid).resolve()
    try:
        p.relative_to(_SESSIONS.resolve())
    except ValueError:
        return None
    return p if p.is_dir() else None


def _run_subprocess(argv: list[str], cwd: Path) -> tuple[int, str]:
    try:
        r = subprocess.run(
            argv,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=3600,
        )
    except subprocess.TimeoutExpired:
        return -1, "subprocess timeout (3600s)"
    except OSError as e:
        return -1, str(e)
    out = (r.stdout or "") + ("\n" if r.stdout and r.stderr else "") + (r.stderr or "")
    return r.returncode, out.strip() or f"(exit {r.returncode})"


def build_m4_status(session_id: str) -> dict[str, Any]:
    err: str | None = None
    sd = _safe_session_dir(session_id)
    if not session_id.strip():
        err = "請提供 session（16 位 hex）"
    elif not sd:
        err = f"找不到 session：{session_id!r}"

    base: dict[str, Any] = {
        "ok": sd is not None,
        "error": err,
        "session_id": session_id.strip() if session_id else "",
        "session_dir_abs": str(sd.resolve()) if sd else "",
        "target": "",
        "pdb_id": "",
        "prep_ready": False,
        "engine_a": {"complete": False, "summary": "", "status": "", "raw": None},
        "engine_b": {"complete": False, "summary": "", "status": "", "raw": None},
        "engine_c": {"complete": False, "summary": "", "raw": None, "pocket_config": None},
        "cli": {"engine_b": "", "engine_c": "", "chdir": str(_PKG_ROOT.resolve())},
    }
    if not sd:
        return base

    meta = _read_json(sd / "session_meta.json")
    prep = _read_json(sd / "prep_decisions.json")
    target = (
        (prep or {}).get("target")
        or (meta or {}).get("target")
        or ""
    )
    if isinstance(target, str):
        target = target.strip()
    pdb_id = (
        (prep or {}).get("pdb_id")
        or (meta or {}).get("pdb_id")
        or ""
    )
    if isinstance(pdb_id, str):
        pdb_id = pdb_id.strip().upper()

    base["target"] = target or ""
    base["pdb_id"] = pdb_id or ""
    base["prep_ready"] = bool(prep and prep.get("ready_for_engine_e"))

    s21 = _read_json(sd / "stage_21_hetatm.json")
    base["engine_a"]["raw"] = s21
    if s21:
        nc = s21.get("candidate_count", 0)
        st = (s21.get("status") or "").upper()
        base["engine_a"]["summary"] = f"{nc} candidates · {st}" if st else str(nc)
        base["engine_a"]["status"] = st
        base["engine_a"]["complete"] = True
    else:
        base["engine_a"]["summary"] = "—"

    s22 = _read_json(sd / "stage_22_ligand.json")
    base["engine_b"]["raw"] = s22
    if s22 and s22.get("ready_for_engine_c"):
        het = s22.get("selected_het_id", "—")
        ch = s22.get("chain", "")
        rs = s22.get("resseq", "")
        st = (s22.get("status") or "").upper()
        base["engine_b"]["summary"] = f"{het} chain {ch} res {rs}"
        base["engine_b"]["status"] = st
        base["engine_b"]["complete"] = True
    else:
        base["engine_b"]["summary"] = "—"

    s23 = _read_json(sd / "stage_23_pocket_center.json")
    pc = _read_json(sd / "pocket_config.json")
    base["engine_c"]["raw"] = s23
    base["engine_c"]["pocket_config"] = pc
    if pc and pc.get("ready_for_docking"):
        cx, cy, cz = pc.get("center_x"), pc.get("center_y"), pc.get("center_z")
        base["engine_c"]["summary"] = f"({cx}, {cy}, {cz})"
        base["engine_c"]["complete"] = True
    elif s23 and (s23.get("status") in ("CONFIRM", "MANUAL")):
        base["engine_c"]["summary"] = (
            f"({s23.get('center_x')}, {s23.get('center_y')}, {s23.get('center_z')})"
        )
        base["engine_c"]["complete"] = bool(pc and pc.get("ready_for_docking"))
    else:
        base["engine_c"]["summary"] = "—"

    py = sys.executable
    tq = target or "YOUR_TARGET"
    sd_abs = str(sd.resolve())
    base["cli"]["engine_b"] = (
        f'"{py}" "{_PKG_ROOT / "molgate_module4_engineB.py"}" '
        f'--session-dir "{sd_abs}" --target {tq}'
    )
    base["cli"]["engine_c"] = (
        f'"{py}" "{_PKG_ROOT / "molgate_module4_engineC.py"}" '
        f'--session-dir "{sd_abs}" --target {tq}'
    )

    return base


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt: str, *args: Any) -> None:
        sys.stderr.write(f"[molgate-m4] {args[0]}\n")

    def _send(self, code: int, body: bytes, content_type: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_json(self, code: int, obj: Any) -> None:
        data = json.dumps(obj, ensure_ascii=False).encode("utf-8")
        self._send(code, data, "application/json; charset=utf-8")

    def do_GET(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path in ("/", "/molgate_module4_panel.html"):
            if not _PANEL.is_file():
                self._send(404, b"panel html missing", "text/plain; charset=utf-8")
                return
            self._send(200, _PANEL.read_bytes(), "text/html; charset=utf-8")
            return

        if path == "/api/module4/status":
            qs = urllib.parse.parse_qs(parsed.query)
            sid = (qs.get("session") or [""])[0]
            self._send_json(200, build_m4_status(sid))
            return

        self._send(404, b"not found", "text/plain; charset=utf-8")

    def do_POST(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path != "/api/module4/run-engine-a":
            self._send(404, b"not found", "text/plain; charset=utf-8")
            return

        length = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(length) if length else b"{}"
        try:
            payload = json.loads(raw.decode("utf-8"))
        except json.JSONDecodeError:
            self._send_json(400, {"ok": False, "error": "invalid JSON"})
            return

        session_id = (payload.get("session") or "").strip()
        sd = _safe_session_dir(session_id)
        if not sd:
            self._send_json(400, {"ok": False, "error": "invalid or missing session id"})
            return

        target = (payload.get("target") or "").strip()
        index = (payload.get("index") or "").strip()

        argv = [
            sys.executable,
            str(_PKG_ROOT / "molgate_module4_engineA.py"),
            "--session-dir",
            str(sd),
        ]
        if target:
            argv.extend(["--target", target])
        if index:
            argv.extend(["--index", index])

        code, output = _run_subprocess(argv, _PKG_ROOT)
        status = build_m4_status(session_id)
        self._send_json(
            200,
            {
                "ok": code == 0,
                "returncode": code,
                "log": output,
                "status": status,
            },
        )


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="MolGate Module 4 panel server")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8766)
    args = ap.parse_args()
    httpd = HTTPServer((args.host, args.port), Handler)
    print(f"MolGate Module 4 panel → http://{args.host}:{args.port}/")
    print("範例：http://127.0.0.1:8766/?session=你的session&target=COX2")
    print("按 Ctrl+C 結束")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n已停止")


if __name__ == "__main__":
    main()
