@echo off
setlocal
cd /d "%~dp0"
chcp 65001 >nul
set PYTHONUTF8=1
python "engines\molgate_runner_pretty.py" %*
endlocal
