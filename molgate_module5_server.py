"""
MolGate — Module 5 本機面板後端
- 提供靜態 molgate_module5_panel.html
- GET  /api/module5/status?session=…
- POST /api/module5/run-engine-a
- POST /api/module5/run-engine-b

使用方式（在 MolGate 目錄）:
  py -3 molgate_module5_server.py
瀏覽器:
  http://127.0.0.1:8767/?session=你的16位hex&target=COX2
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
_PANEL = _PKG_ROOT / "molgate_module5_panel.html"
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


def build_m5_status(session_id: str) -> dict[str, Any]:
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
        "drug": "",
        "engine_a": {"complete": False, "summary": "", "status": "", "raw24": None, "raw25": None, "raw26": None},
        "engine_b": {"complete": False, "summary": "", "status": "", "raw27": None},
        "docking_result": None,
        "final_report": None,
    }
    if not sd:
        return base

    meta = _read_json(sd / "session_meta.json")
    dr = _read_json(sd / "docking_result.json")
    fr = _read_json(sd / "final_report.json")
    s24 = _read_json(sd / "stage_24_vina.json")
    s25 = _read_json(sd / "stage_25_qc.json")
    s26 = _read_json(sd / "stage_26_result.json")
    s27 = _read_json(sd / "stage_27_review.json")

    target = (dr or {}).get("target") or (meta or {}).get("target") or ""
    pdb_id = (dr or {}).get("pdb_id") or (meta or {}).get("pdb_id") or ""
    drug = (dr or {}).get("target_drug") or ""

    base["target"] = str(target).strip()
    base["pdb_id"] = str(pdb_id).strip().upper()
    base["drug"] = str(drug).strip()

    base["docking_result"] = dr
    base["final_report"] = fr
    base["engine_a"]["raw24"] = s24
    base["engine_a"]["raw25"] = s25
    base["engine_a"]["raw26"] = s26
    base["engine_b"]["raw27"] = s27

    if dr and dr.get("ready_for_review"):
        aff = dr.get("best_affinity")
        tier = dr.get("confidence_tier", "")
        qc = "PASS" if dr.get("qc_pass") else "WARN"
        base["engine_a"]["summary"] = f"{aff} kcal/mol · tier {tier} · {qc}"
        base["engine_a"]["status"] = "PASS" if dr.get("qc_pass") else "WARN"
        base["engine_a"]["complete"] = True

    if s27:
        dec = str(s27.get("decision") or s27.get("status") or "").upper()
        base["engine_b"]["summary"] = dec or "—"
        base["engine_b"]["status"] = dec or ""
        base["engine_b"]["complete"] = dec == "ACCEPT"

    return base


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt: str, *args: Any) -> None:
        sys.stderr.write(f"[molgate-m5] {args[0]}\n")

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

        if path in ("/", "/molgate_module5_panel.html"):
            if not _PANEL.is_file():
                self._send(404, b"panel html missing", "text/plain; charset=utf-8")
                return
            self._send(200, _PANEL.read_bytes(), "text/html; charset=utf-8")
            return

        if path == "/api/module5/status":
            qs = urllib.parse.parse_qs(parsed.query)
            sid = (qs.get("session") or [""])[0]
            self._send_json(200, build_m5_status(sid))
            return

        self._send(404, b"not found", "text/plain; charset=utf-8")

    def do_POST(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
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

        if parsed.path == "/api/module5/run-engine-a":
            target = (payload.get("target") or "").strip()
            vina_bin = (payload.get("vina_bin") or "").strip()
            argv = [
                sys.executable,
                str(_PKG_ROOT / "molgate_module5_engineA.py"),
                "--session-dir",
                str(sd),
            ]
            if target:
                argv.extend(["--target", target])
            if vina_bin:
                argv.extend(["--vina-bin", vina_bin])
            code, output = _run_subprocess(argv, _PKG_ROOT)
            status = build_m5_status(session_id)
            self._send_json(200, {"ok": code == 0, "returncode": code, "log": output, "status": status})
            return

        if parsed.path == "/api/module5/run-engine-b":
            decision = str(payload.get("decision") or "").strip().upper()
            note = str(payload.get("note") or "").strip()
            retry_module = str(payload.get("retry_module") or "").strip()
            if decision not in ("ACCEPT", "RETRY", "FLAG"):
                self._send_json(400, {"ok": False, "error": "decision must be ACCEPT/RETRY/FLAG"})
                return
            argv = [
                sys.executable,
                str(_PKG_ROOT / "molgate_module5_engineB.py"),
                "--session-dir",
                str(sd),
                "--decision",
                decision,
            ]
            if note:
                argv.extend(["--note", note])
            if decision == "RETRY" and retry_module in ("1", "2", "3", "4", "5"):
                argv.extend(["--retry-module", retry_module])
            code, output = _run_subprocess(argv, _PKG_ROOT)
            status = build_m5_status(session_id)
            self._send_json(200, {"ok": code == 0, "returncode": code, "log": output, "status": status})
            return

        self._send(404, b"not found", "text/plain; charset=utf-8")


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="MolGate Module 5 panel server")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8767)
    args = ap.parse_args()
    httpd = HTTPServer((args.host, args.port), Handler)
    print(f"MolGate Module 5 panel → http://{args.host}:{args.port}/")
    print("範例：http://127.0.0.1:8767/?session=你的session&target=COX2")
    print("按 Ctrl+C 結束")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n已停止")


if __name__ == "__main__":
    main()

