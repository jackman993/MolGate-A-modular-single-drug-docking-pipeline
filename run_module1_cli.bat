@echo off
setlocal

set "ROOT=%~dp0"
cd /d "%ROOT%"

if "%~3"=="" (
  echo Usage: run_module1_cli.bat "SMILES" "TARGET" "PDB_ID"
  echo Example: run_module1_cli.bat "CC(=O)OC1=CC=CC=C1C(=O)O" "COX2" "7DFP"
  exit /b 1
)

set "SMILES=%~1"
set "TARGET=%~2"
set "PDBID=%~3"

python "engines\molgate_module1.py" --smiles "%SMILES%" --target "%TARGET%" --pdb-id "%PDBID%" --index "engines\master_index.json"

endlocal
