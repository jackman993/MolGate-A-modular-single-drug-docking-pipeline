@echo off
setlocal
cd /d "%~dp0"

echo [MolGate] Full pipeline from Module 1 (drug index 0)
echo ROOT=%CD%
echo.

call run_runner_pretty.bat --drug 0 --from-stage 1 --allow-network --non-interactive --vina-bin "%CD%\tool\vina\vina.exe"
set "RC=%ERRORLEVEL%"

echo.
echo [Done] exit code=%RC%
endlocal
exit /b %RC%
