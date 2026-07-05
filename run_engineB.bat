@echo off
setlocal

set "ROOT=%~dp0"
cd /d "%ROOT%"

if "%~2"=="" (
  echo Usage: run_engineB.bat "SESSION_DIR" "TARGET"
  echo Example: run_engineB.bat "molgate_sessions\abc123" "COX2"
  exit /b 1
)

python "engines\molgate_module2_engineB.py" --session-dir "%~1" --target "%~2"

endlocal
