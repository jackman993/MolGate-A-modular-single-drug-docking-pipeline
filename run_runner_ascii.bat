@echo off
setlocal
cd /d "%~dp0"
python "engines\molgate_runner_ascii.py" %*
endlocal
