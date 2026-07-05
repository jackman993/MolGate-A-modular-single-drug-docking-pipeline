@echo off
setlocal
cd /d "%~dp0"
if "%~1"=="" (
  echo Usage: run_catalog_drug.bat DRUG_ID
  echo Example: run_catalog_drug.bat 26
  exit /b 1
)
call run_runner_pretty.bat --catalog-drug-id %1 --from-stage 1 --allow-network --non-interactive --vina-bin "%CD%\tool\vina\vina.exe"
endlocal
