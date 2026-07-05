@echo off
setlocal

set "ROOT=%~dp0.."
cd /d "%ROOT%"

echo [MolGate] Starting Module5 server on 8767
echo ROOT=%CD%
echo.

python "molgate_module5_server.py" --host 127.0.0.1 --port 8767

endlocal

