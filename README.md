# MolGate: A Modular Single-Drug Docking Pipeline (DrugOps)

Standalone **command-line** release for the DrugOps / MolGate virtual-screening pipeline (Windows CMD, no web UI). This release includes the catalog-300 registry, runtime anchor routing (`f_route`), automated QC taxonomy, and per-entry batch outcomes used by the DrugOps paper.

## Repository layout

| Path | Description |
|------|-------------|
| `engines/` | Module 1–5 Python pipeline scripts |
| `engines/config/drug_catalog_300.csv` | **300-entry** drug–target catalog (`--catalog-drug-id`) |
| `engines/master_index.json` | Structural anchor registry (**304** index keys) |
| `engines/molgate_anchor_routing.py` | Runtime anchor routing (`f_route`): Direct -> same-PDB surrogate -> cross-PDB proxy |
| `engines/molgate_route_runner.py` | Runner integration for anchor routing |
| `examples/repaglinide_route_showcase.json` | Paper showcase case for runtime routing (`drug_id=62`) |
| `results/` | **300-case** registry + batch run outcomes (paper Data Availability) |
| `tool/vina/vina.exe` | AutoDock Vina (Windows) |
| `*.bat` | One-click launch scripts |
| `run_from_scratch.bat` | Full 27-step pipeline demo (DRUGS index 0, Ibuprofen) |

---

<img width="1077" height="652" alt="574605006-13eadd18-29af-4672-aa2c-48a2b1d9e252" src="https://github.com/user-attachments/assets/c36ae9c6-5d50-4330-bade-c697c3552533" />


## Quick start (Windows CMD)

```bat
cd /d C:\path\to\MolGate-A-modular-single-drug-docking-pipeline
check_env.bat
install_deps.bat
run_from_scratch.bat
```

---

## Running a single drug

### CreaDrug24 demo (`--drug 0` … `--drug 23`)

Built-in array of 24 seed drugs (maps to catalog `drug_id` 1–24):

```bat
run_runner_pretty.bat --drug 0 --from-stage 1 --allow-network --non-interactive --vina-bin "%CD%\tool\vina\vina.exe"
```

ASCII-safe output (any CMD code page):

```bat
run_runner_ascii.bat --drug 0 --from-stage 1 --allow-network --non-interactive --vina-bin "%CD%\tool\vina\vina.exe"
```

### Catalog-300 (`--catalog-drug-id 1` … `300`)

All **300** curated entries are in `engines/config/drug_catalog_300.csv`. Run any entry by catalog ID:

```bat
run_catalog_drug.bat 26
```

Equivalent:

```bat
run_runner_pretty.bat --catalog-drug-id 26 --from-stage 1 --allow-network --non-interactive --vina-bin "%CD%\tool\vina\vina.exe"
```

| Catalog range | How to run |
|---------------|------------|
| `drug_id` 1–24 | `--drug 0` … `--drug 23` **or** `--catalog-drug-id 1` … `24` |
| `drug_id` 25–300 | `--catalog-drug-id 25` … `300` |

**Sync catalog from UI workspace** (after editing `moltgate_UI/config/`):

```bat
scripts\sync_catalog_from_ui.bat
```

### Runtime anchor routing (`f_route`)

Before Module 1, the runner verifies whether the catalog HET is present in the declared PDB. If not, it applies the same routing policy used in the catalog-300 study:

```text
Direct -> same-PDB surrogate -> cross-PDB proxy
```

The routing decision is written to each session as `route_resolution.json` and merged into `session_index.json` before downstream modules run.

Paper showcase: **Repaglinide** (`drug_id=62`, `REPAGLINIDE_KATP`):

```bat
run_catalog_drug.bat 62
```

Expected route behavior:

| Field | Value |
|-------|-------|
| Catalog declaration | `3VE0` / `REP` |
| Effective route | `cross_pdb_proxy` |
| Effective PDB | `3SYQ` |
| Proxy source | `GLYBURIDE_KATP` |

Use `--skip-routing` only for debugging legacy behavior without `f_route`.

---

## Catalog-300 results (paper release)

The `results/` folder ships the **300-entry structural anchor registry** and **276 automated batch pipeline outcomes** (`drug_id` 25–300).

### Files

| File | Description |
|------|-------------|
| `results/catalog_300_registry.csv` | All **300** curated entries: SMILES, `index_key`, PDB, HET, tier, `anchor_mode` |
| `results/catalog_300_run_summary.csv` | Per-entry outcomes: ΔG, `qc_pass`, `paper_outcome`, route, `session_id` |
| `results/catalog_300_run_summary.json` | Same data, machine-readable |
| `results/stats_300.json` | Aggregate statistics |

### Coverage

| Range | Count | Source |
|-------|------:|--------|
| `drug_id` 1–24 | 24 | **CreaDrug24 seed set** — in registry; batch runner skipped these. Reproduce via `--drug 0` … `--drug 23`. |
| `drug_id` 25–300 | 276 | **Batch pipeline runs** — rows in `catalog_300_run_summary.*` |
| **Registry total** | **300** | Full catalog in `catalog_300_registry.csv` |

### Batch statistics (`stats_300.json`, 276 runs)

| Metric | Value |
|--------|------:|
| Completed (`status=ok`) | 232 |
| Pipeline-fail | 44 |
| **Pass** (`qc_pass=true`) | **181** |
| **QC-flagged** (`qc_pass=false`) | **51** |
| Tier 1 / Tier 2 | 151 / 125 |
| Best affinity range (kcal/mol) | −11.66 … −2.73 (mean −7.96) |

**Route resolution** (implemented in this CMD release via `engines/molgate_anchor_routing.py`):

| `route_status` | Count |
|----------------|------:|
| `same_pdb_surrogate` | 137 |
| `direct_ok` | 71 |
| `no_route` | 34 |
| `cross_pdb_proxy` | 17 |
| `unknown` | 17 |

Regenerate results bundle from the UI workspace:

```bat
python ..\moltgate_UI\scripts\export_catalog_300_release.py --out-dir results
```

---

## Outcome taxonomy (DrugOps paper)

The **paper** classifies each run using **automated QC** (`qc_pass` in `docking_result.json` / `final_report.json`), **not** the Stage-27 review label.

| Paper class | Condition |
|-------------|-----------|
| **Pass** | Docking complete **and** `qc_pass=true` |
| **QC-flagged** | Docking complete **and** `qc_pass=false` |
| **Pipeline-fail** | No complete docking outcome |

With `--non-interactive`, Stage 27 is **forced to ACCEPT** so every completed run writes `final_report.json`. That **ACCEPT** means *artifacts archived* — it does **not** mean QC passed.

Validated QC-flagged example: **Aspirin** (`drug_id=28`, `ASPIRIN_COX2`).

```bat
run_catalog_drug.bat 28
```

Observed result in this CMD release:

| Field | Value |
|-------|-------|
| Route | `direct_ok` (`TLF` present in `5IKT`) |
| Best affinity | `-6.688 kcal/mol` |
| Confidence tier | `C` |
| `qc_pass` | `false` |
| QC flag | `OUT_OF_POCKET_MAJORITY` |
| Stage-27 decision | `ACCEPT` |
| Paper outcome | **QC-flagged** |

This illustrates the intended separation: Stage-27 `ACCEPT` archives artifacts in batch mode, while `qc_pass=false` remains **QC-flagged** for the paper.

Check classification:

```bat
findstr qc_pass engines\molgate_sessions\YOUR_SESSION\docking_result.json
```

---

## Smoke tests before release

Run these from the repository root:

```bat
:: direct_ok + QC-flagged case
run_catalog_drug.bat 28

:: cross-PDB proxy routing showcase
run_catalog_drug.bat 62

:: built-in CreaDrug24 seed demo
run_runner_pretty.bat --drug 0 --from-stage 1 --allow-network --non-interactive --vina-bin "%CD%\tool\vina\vina.exe"
```

Expected artifacts:

| Test | Key artifact |
|------|--------------|
| Aspirin `drug_id=28` | `docking_result.json` with `qc_pass=false` and `OUT_OF_POCKET_MAJORITY` |
| Repaglinide `drug_id=62` | `route_resolution.json` with `status=cross_pdb_proxy` |
| Ibuprofen `--drug 0` | Complete 27-stage session and final report |

---

## Data Availability (paper citation text)

> The 300-entry structural anchor registry, runtime anchor-routing implementation (`f_route`), and per-entry docking summary statistics are available in the MolGate command-line release (`engines/molgate_anchor_routing.py`, `results/catalog_300_registry.csv`, `results/catalog_300_run_summary.csv`, `results/stats_300.json`) at GitHub tag **v1.0.0-catalog300**. The command-line pipeline is in the same repository (`engines/`). Interactive browsing is provided separately (Hugging Face Space). Full session artifacts are available on request.

Suggested release tag: **`v1.0.0-catalog300`**

---

## Dependencies

- Python 3.10+ (conda recommended)
- rdkit, meeko, scipy, gemmi, biopython
- Run `install_deps.bat` to install most packages

---

## Publishing to GitHub

1. Create a new GitHub repository.
2. Push the **entire** release folder (`.bat` files included).
3. Do **not** push run outputs under `engines/molgate_sessions/` (listed in `.gitignore`).

```bat
cd /d C:\path\to\MolGate-A-modular-single-drug-docking-pipeline
git init
git add .
git commit -m "Add MolGate CMD release with catalog-300 results"
git branch -M main
git remote add origin https://github.com/YOUR_USER/YOUR_REPO.git
git push -u origin main
git tag v1.0.0-catalog300
git push --tags
```
