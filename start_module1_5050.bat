@echo off
setlocal

set "ROOT=%~dp0.."
cd /d "%ROOT%"

echo [MolGate] Starting Module1 server on 5050
echo ROOT=%CD%
echo.

set "MOLGATE_PORT=5050"
python "module1\molgate_server.py"

endlocal

