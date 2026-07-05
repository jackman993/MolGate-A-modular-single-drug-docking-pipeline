@echo off
setlocal

set "ROOT=%~dp0"
cd /d "%ROOT%"

if "%~2"=="" (
  echo Usage: run_engineC.bat "SESSION_DIR" "TARGET"
  echo Example: run_engineC.bat "molgate_sessions\abc123" "COX2"
  exit /b 1
)

python "engines\molgate_module2_engineC.py" --session-dir "%~1" --target "%~2" --engine meeko

endlocal
