# MolGate-A-modular-single-drug-docking-pipeline
A modular single-drug docking pipeline for structure-based screening and reproducible workflow experiments. > Current release focuses on stability and reproducibility. > Runtime console logs may still contain Traditional Chinese messages in some modules.
---
## Highlights
- Modular 5-stage pipeline (Input -> Ligand -> Receptor -> Pocket -> Docking/Review)
- Session-based outputs for reproducibility
- Meeko-first strategy with fallback options (Open Babel / MGLTools, depending on module)
- Supports direct anchor and proxy anchor index design (`anchor_mode`)
---
## Repository Structure
MolGate/
├── molgate_runner.py
├── molgate_module1.py
├── molgate_module2_engineB.py
├── molgate_module2_engineC.py
├── molgate_module3_engineA.py
├── molgate_module3_engineB.py
├── molgate_module3_engineC.py
├── molgate_module4_engineA.py
├── molgate_module4_engineB.py
├── molgate_module4_engineC.py
├── molgate_module5_engineA.py
├── molgate_module5_engineB.py
├── config/
├── docs/
├── data/
├── molgate_sessions/
└── tools/

## Pipeline Overview
Module 1: Input validation, PDB loading, index lookup
Module 2: Ligand protonation/tautomer handling and ligand PDBQT generation
Module 3: Receptor preparation and receptor PDBQT generation
Module 4: Pocket ligand identification and docking box generation
Module 5: Vina docking, QC, and review decision
See docs/PIPELINE_OVERVIEW.md for details.
Prerequisites
Python 3.11+
AutoDock Vina (recommended local binary path under tools/vina/)
Optional tools depending on chosen engine:
meeko
Open Babel (obabel)
MGLTools (prepare_receptor4.py) for legacy fallback

## Installation
Quick Start
## 1) Prepare index files
Copy sample config files and edit for your targets:
config/master_index.sample.json -> your master_index.json
config/drug_library.sample.json (optional)
## 2) Run a single drug with runner
## 3) Resume from a later module (example)

## Engine Policy (Current)
Ligand path (Module 2 Engine C): default meeko
Receptor path (Module 3 Engine C): auto = Meeko -> Open Babel -> MGLTools
Recommended practice: use Meeko first; fallback only on failure or quality concerns

## Index Design Notes
master_index.json is currently the critical control plane.
Recommended required fields per entry:
pdb_id
het_id / cocrystal_het
pocket_rule
target_drug
common_name
notes
anchor_mode (direct or proxy, recommended)
See config/index_schema.md and docs/USER_INDEX_GUIDE.md.

## Known Limitations
Index generation is not yet fully automated from remote APIs
Some modules still output Traditional Chinese console messages
PDB/CCD mapping may require manual validation for difficult targets
Legacy toolchain fallback may produce format edge cases requiring cleanup

## Troubleshooting
See docs/TROUBLESHOOTING.md.
Common issues include:
Vina PDBQT parsing errors
Windows console encoding problems
Meeko --allow_bad_res behavior and receptor quality checks

## Reproducibility
Each run creates a session folder with staged JSON artifacts
Keep molgate_sessions/ out of git for clean repository history

## Roadmap
Full English runtime logs
Better automated index assistant (API + validation layer)
Stronger receptor/ligand quality gates before docking
Batch mode for high-throughput runs
## License
[Choose and add your license here, e.g. MIT]
