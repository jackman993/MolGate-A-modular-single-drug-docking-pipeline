"""
MolGate — Module 3 本機面板後端
提供靜態 panel HTML、session 狀態 JSON、subprocess 執行 Engine A/B/C，
以及 Ollama HTTP API 代理（/api/ollama/chat，避免瀏覽器 CORS）。

環境變數：OLLAMA_HOST（預設 http://127.0.0.1:11434）

使用方式（在 MolGate 目錄）：
  py -3.11 molgate_module3_server.py
瀏覽器開啟 http://127.0.0.1:8765/
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from typing import Any

_PKG_ROOT = Path(__file__).resolve().parent
_OLLAMA_BASE = (os.environ.get("OLLAMA_HOST") or "http://127.0.0.1:11434").rstrip("/")
_SESSIONS = _PKG_ROOT / "molgate_sessions"
_PANEL = _PKG_ROOT / "molgate_module3_panel.html"

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


def _stage_status_class(rec: dict[str, Any] | None) -> tuple[str, str, str]:
    """回傳 (顯示文字, status CSS class, 狀態標籤)"""
    if not rec:
        return "—", "", ""
    st = (rec.get("status") or "").upper()
    if st == "AUTO":
        return st, "s-auto", "AUTO"
    if st == "WARN":
        return st, "s-warn", "WARN"
    if st == "MANUAL":
        return st, "s-manual", "MANUAL"
    if st == "SKIPPED":
        return "SKIPPED", "s-skip", "SKIPPED"
    if st == "PASS":
        return "PASS", "s-auto", "PASS"
    if st == "FAIL":
        return "FAIL", "s-fail", "FAIL"
    return st or "—", "s-auto", st or "—"


def _summarize_stage_a_row(n: int, rec: dict[str, Any] | None) -> dict[str, Any]:
    if not rec:
        return {"num": n, "val": "—", "statusClass": "", "statusText": "", "raw": None}
    val = "—"
    if n == 10:
        val = rec.get("family", "—")
    elif n == 11:
        val = rec.get("state", "—")
    elif n == 12:
        kc = rec.get("kept_chains") or []
        val = f"保留 {','.join(kc)}" if kc else "—"
    elif n == 13:
        mc = rec.get("missing_count", 0)
        val = f"{mc} 處缺失" if mc else "無缺失環"
    elif n == 14:
        wc = rec.get("water_count", 0)
        act = rec.get("action", "")
        val = f"水 {wc} 個 · {act}"
    cls, sc, stxt = _stage_status_class(rec)
    return {"num": n, "val": val, "statusClass": sc, "statusText": stxt, "raw": rec}


def build_status(session_id: str) -> dict[str, Any]:
    err: str | None = None
    sd = _safe_session_dir(session_id)
    if not session_id.strip():
        err = "請輸入 session id（16 位 hex）"
    elif not sd:
        err = f"找不到 session：{session_id!r}"

    out: dict[str, Any] = {
        "ok": sd is not None,
        "error": err,
        "session_id": session_id.strip() if session_id else "",
        "target": "",
        "pdb_id": "",
        "paths": {},
        "engine_a": {"complete": False, "stages": [], "cleaned_exists": False},
        "engine_b": {"complete": False, "stages": [], "prep": None, "ready_for_engine_e": False},
        "engine_c": {"complete": False, "stages": [], "receptor_exists": False},
    }
    if not sd:
        return out

    s3 = _read_json(sd / "stage_03_index.json")
    man = _read_json(sd / "molgate_run_manifest.json")
    out["target"] = (s3 or {}).get("target") or (man or {}).get("target") or ""
    out["pdb_id"] = (s3 or {}).get("pdb_id") or (man or {}).get("pdb_id") or ""
    out["paths"] = {
        "cleaned_pdb": str((sd / "protein" / "cleaned.pdb").as_posix()),
        "prep_decisions": str((sd / "prep_decisions.json").as_posix()),
        "receptor_pdbqt": str((sd / "protein" / "receptor.pdbqt").as_posix()),
    }

    stages_a = []
    for i, name in enumerate(
        ["family", "state", "chains", "loops", "water"],
        start=10,
    ):
        rec = _read_json(sd / f"stage_{i:02d}_{name}.json")
        stages_a.append(_summarize_stage_a_row(i, rec))
    out["engine_a"]["stages"] = stages_a
    cleaned = sd / "protein" / "cleaned.pdb"
    out["engine_a"]["cleaned_exists"] = cleaned.is_file()
    out["engine_a"]["complete"] = cleaned.is_file() and _read_json(sd / "stage_14_water.json") is not None

    b_stages = []
    for i, name in [(15, "cofactor"), (16, "propka"), (17, "override"), (18, "saltbridge")]:
        rec = _read_json(sd / f"stage_{i:02d}_{name}.json")
        val = "—"
        if rec:
            if i == 15:
                val = f"{rec.get('action', '—')} · {rec.get('status', '')}"
            elif i == 17:
                rules = rec.get("override_rules") or {}
                val = f"{len(rules)} 條規則" if rules else rec.get("status", "—")
            else:
                val = rec.get("flag") or rec.get("status", "—")
        cls, sc, stxt = _stage_status_class(rec)
        b_stages.append({"num": i, "val": val, "statusClass": sc, "statusText": stxt, "raw": rec})
    out["engine_b"]["stages"] = b_stages
    prep = _read_json(sd / "prep_decisions.json")
    out["engine_b"]["prep"] = prep
    rdy = bool(prep and prep.get("ready_for_engine_e"))
    out["engine_b"]["ready_for_engine_e"] = rdy
    out["engine_b"]["complete"] = rdy

    c_stages = []
    for i, name in [(19, "receptor_pdbqt"), (20, "receptor_geometry_qc")]:
        rec = _read_json(sd / f"stage_{i:02d}_{name}.json")
        val = "—"
        if rec:
            if i == 19:
                val = rec.get("engine") or rec.get("status", "—")
                if rec.get("receptor_pdbqt"):
                    val = f"{val} → receptor.pdbqt"
            elif i == 20:
                val = f"AA≈{rec.get('approx_aa_residues', '—')} · {rec.get('status', '')}"
        cls, sc, stxt = _stage_status_class(rec)
        c_stages.append({"num": i, "val": val, "statusClass": sc, "statusText": stxt, "raw": rec})
    out["engine_c"]["stages"] = c_stages
    rp = sd / "protein" / "receptor.pdbqt"
    out["engine_c"]["receptor_exists"] = rp.is_file()
    s19 = _read_json(sd / "stage_19_receptor_pdbqt.json")
    out["engine_c"]["complete"] = bool(
        s19 and s19.get("status") == "PASS" and rp.is_file()
    )

    # Summary 用
    s14 = _read_json(sd / "stage_14_water.json")
    s12 = _read_json(sd / "stage_12_chains.json")
    out["summary"] = {
        "cleaned_rel": "protein/cleaned.pdb" if cleaned.is_file() else "—",
        "chains": ",".join(s12.get("kept_chains") or []) if s12 else "—",
        "water_removed": (
            f"{s14.get('water_count', '—')} 個" if s14 else "—"
        ),
        "receptor_rel": "protein/receptor.pdbqt" if rp.is_file() else "—",
        "engine_c": (s19 or {}).get("engine") or "—",
        "flags": (prep or {}).get("flags") if prep else [],
    }
    return out


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


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt: str, *args: Any) -> None:
        sys.stderr.write(f"[molgate-m3] {args[0]}\n")

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

        if path in ("/", "/molgate_module3_panel.html"):
            if not _PANEL.is_file():
                self._send(404, b"panel html missing", "text/plain; charset=utf-8")
                return
            body = _PANEL.read_bytes()
            self._send(200, body, "text/html; charset=utf-8")
            return

        if path == "/api/status":
            qs = urllib.parse.parse_qs(parsed.query)
            sid = (qs.get("session") or [""])[0]
            self._send_json(200, build_status(sid))
            return

        if path == "/api/sessions":
            if not _SESSIONS.is_dir():
                self._send_json(200, {"ok": True, "sessions": []})
                return
            dirs = sorted(
                p.name for p in _SESSIONS.iterdir() if p.is_dir() and _SESSION_ID_RE.match(p.name)
            )
            self._send_json(200, {"ok": True, "sessions": dirs})
            return

        self._send(404, b"not found", "text/plain; charset=utf-8")

    def _handle_ollama_chat(self, payload: dict[str, Any]) -> None:
        """轉發至 Ollama /api/chat；stream=true 時以 SSE 回傳 NDJSON 行。"""
        messages = payload.get("messages")
        if not isinstance(messages, list) or not messages:
            self._send_json(400, {"ok": False, "error": "messages 必須為非空陣列"})
            return
        model = (payload.get("model") or "gemma3").strip()
        stream = bool(payload.get("stream", True))
        body_obj: dict[str, Any] = {
            "model": model,
            "messages": messages,
            "stream": stream,
        }
        if "options" in payload and isinstance(payload["options"], dict):
            body_obj["options"] = payload["options"]
        body = json.dumps(body_obj, ensure_ascii=False).encode("utf-8")
        url = f"{_OLLAMA_BASE}/api/chat"
        req = urllib.request.Request(
            url,
            data=body,
            method="POST",
            headers={"Content-Type": "application/json; charset=utf-8"},
        )
        if not stream:
            try:
                with urllib.request.urlopen(req, timeout=600) as resp:
                    raw_out = resp.read().decode("utf-8")
                self._send_json(200, json.loads(raw_out))
            except urllib.error.HTTPError as e:
                err_body = e.read().decode("utf-8", errors="replace")
                self._send_json(502, {"ok": False, "error": f"ollama HTTP {e.code}: {err_body}"})
            except Exception as e:
                self._send_json(502, {"ok": False, "error": str(e)})
            return

        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream; charset=utf-8")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("X-Accel-Buffering", "no")
        self.end_headers()
        try:
            with urllib.request.urlopen(req, timeout=600) as resp:
                while True:
                    line = resp.readline()
                    if not line:
                        break
                    line = line.strip()
                    if not line:
                        continue
                    self.wfile.write(b"data: " + line + b"\n\n")
                    self.wfile.flush()
        except Exception as e:
            err = json.dumps({"error": str(e)}, ensure_ascii=False).encode("utf-8")
            self.wfile.write(b"data: " + err + b"\n\n")
            self.wfile.flush()

    def do_POST(self) -> None:
        parsed = urllib.parse.urlparse(self.path)
        length = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(length) if length else b"{}"
        try:
            payload = json.loads(raw.decode("utf-8"))
        except json.JSONDecodeError:
            self._send_json(400, {"ok": False, "error": "invalid JSON"})
            return

        if parsed.path == "/api/ollama/chat":
            self._handle_ollama_chat(payload)
            return

        if parsed.path != "/api/run":
            self._send(404, b"not found", "text/plain; charset=utf-8")
            return

        engine = (payload.get("engine") or "").lower().strip()
        session_id = (payload.get("session") or "").strip()
        target = (payload.get("target") or "").strip()
        sd = _safe_session_dir(session_id)
        if not sd:
            self._send_json(400, {"ok": False, "error": "invalid or missing session id"})
            return
        if not target:
            self._send_json(400, {"ok": False, "error": "missing target"})
            return

        py = sys.executable
        cwd = _PKG_ROOT

        if engine == "a":
            argv = [
                py,
                str(_PKG_ROOT / "molgate_module3_engineA.py"),
                "--session-dir",
                str(sd),
                "--target",
                target,
            ]
        elif engine == "b":
            auto_hets = (payload.get("auto_hets") or "remove_all").strip()
            if auto_hets not in ("remove_all", "keep_all"):
                auto_hets = "remove_all"
            argv = [
                py,
                str(_PKG_ROOT / "molgate_module3_engineB.py"),
                "--session-dir",
                str(sd),
                "--target",
                target,
                "--non-interactive",
                "--auto-hets",
                auto_hets,
            ]
        elif engine == "c":
            eng = (payload.get("engine_c") or "auto").strip()
            if eng not in ("auto", "obabel", "meeko"):
                eng = "auto"
            argv = [
                py,
                str(_PKG_ROOT / "molgate_module3_engineC.py"),
                "--session-dir",
                str(sd),
                "--target",
                target,
                "--engine",
                eng,
            ]
        else:
            self._send_json(400, {"ok": False, "error": "engine must be a, b, or c"})
            return

        code, output = _run_subprocess(argv, cwd)
        status = build_status(session_id)
        self._send_json(
            200,
            {
                "ok": code == 0,
                "returncode": code,
                "output": output,
                "status": status,
            },
        )


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="MolGate Module 3 panel server")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    args = ap.parse_args()
    httpd = HTTPServer((args.host, args.port), Handler)
    print(f"MolGate Module 3 panel → http://{args.host}:{args.port}/")
    print(f"Ollama 代理 → {_OLLAMA_BASE}（可設 OLLAMA_HOST）")
    print("按 Ctrl+C 結束")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n已停止")


if __name__ == "__main__":
    main()
