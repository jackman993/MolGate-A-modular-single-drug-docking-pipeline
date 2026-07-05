@echo off
setlocal
cd /d "%~dp0.."
python "%~dp0sync_catalog_from_ui.py" %*
endlocal
