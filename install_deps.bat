@echo off
setlocal

echo [MolGate] Installing Python dependencies...
python -m pip install --upgrade pip
python -m pip install flask meeko scipy gemmi dimorphite-dl biopython
conda install -y -c conda-forge rdkit

echo [Done]
endlocal
