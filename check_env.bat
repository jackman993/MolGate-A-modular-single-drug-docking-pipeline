@echo off
setlocal
cd /d "%~dp0"

echo [MolGate ENV CHECK]
echo ROOT=%CD%
echo.

where python
python -c "import sys; print('python_exe=', sys.executable); print('python_ver=', sys.version)"
echo.

python -c "mods=['flask','rdkit','meeko','scipy','gemmi','Bio']; import importlib.util as u; print({m: bool(u.find_spec(m)) for m in mods})"
if exist "tool\vina\vina.exe" (echo vina=OK: %CD%\tool\vina\vina.exe) else (echo vina=MISSING)
echo.

echo [Done]
endlocal
