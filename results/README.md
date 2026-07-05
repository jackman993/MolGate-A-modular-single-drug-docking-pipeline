# DrugOps catalog-300 results

Per-entry docking outcomes from the **300-entry structural anchor registry** (DrugOps paper). This folder is the **Data Availability** bundle for GitHub release tag `v1.0.0-catalog300`.

---

## Files

| File | Description |
|------|-------------|
| `catalog_300_registry.csv` | All **300** curated entries: `drug_id`, SMILES, `index_key`, PDB, HET, tier, `anchor_mode` |
| `catalog_300_run_summary.csv` | Pipeline outcomes: best affinity, `qc_pass`, `paper_outcome`, route, `session_id` |
| `catalog_300_run_summary.json` | Same as CSV (machine-readable) |
| `stats_300.json` | Aggregate statistics |

---

## Coverage (300 registry vs 276 run rows)

| Range | Count | Source |
|-------|------:|--------|
| `drug_id` 1–24 | 24 | **CreaDrug24 seed set** — listed in `catalog_300_registry.csv`; batch runner skipped these (`--skip-creadrug`). Reproduce via CMD: `--drug 0` … `--drug 23`. |
| `drug_id` 25–300 | 276 | **Batch pipeline runs** — rows in `catalog_300_run_summary.*` |
| **Registry total** | **300** | Full catalog in `catalog_300_registry.csv` |

---

## Batch statistics (`stats_300.json`)

Summary of **276** automated batch runs (`drug_id` 25–300):

| Metric | Value |
|--------|------:|
| Completed (`status=ok`) | 232 |
| Pipeline-fail | 44 |
| **Pass** (`qc_pass=true`) | **181** |
| **QC-flagged** (`qc_pass=false`) | **51** |
| Tier 1 / Tier 2 | 151 / 125 |
| Best affinity (kcal/mol) | −11.66 … −2.73 (mean −7.96) |

**Route resolution** (from UI batch runs):

| `route_status` | Count |
|----------------|------:|
| `same_pdb_surrogate` | 137 |
| `direct_ok` | 71 |
| `no_route` | 34 |
| `cross_pdb_proxy` | 17 |
| `unknown` | 17 |

---

## Paper outcome taxonomy

| Class | Rule |
|-------|------|
| **Pass** | `status=ok` and `qc_pass=true` |
| **QC-flagged** | `status=ok` and `qc_pass=false` |
| **Pipeline-fail** | `status` not `ok` |

Stage-27 **ACCEPT** in batch CMD mode archives artifacts only; it does **not** override `qc_pass`.

---

## Reproduce one entry

From the repository root:

```bat
run_catalog_drug.bat 26
```

Or:

```bat
run_runner_pretty.bat --catalog-drug-id 26 --from-stage 1 --allow-network --non-interactive --vina-bin "%CD%\tool\vina\vina.exe"
```

Requires `engines/master_index.json` and `engines/config/drug_catalog_300.csv` (sync from UI if needed):

```bat
scripts\sync_catalog_from_ui.bat
```

---

## Regenerate this bundle

From a sibling `moltgate_UI` workspace:

```bat
python ..\moltgate_UI\scripts\export_catalog_300_release.py --out-dir results
```

Or from the repo root:

```bat
scripts\rebuild_results.bat
```

---

## Citation

Cite GitHub release tag **`v1.0.0-catalog300`** and the DrugOps paper Data Availability statement. See the root `README.md` for the full citation paragraph.
