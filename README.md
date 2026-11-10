# Soybean yield responses to ENSO in South America

Public research-code and archived-results repository accompanying the manuscript by **Rogério de Souza Nóia Júnior, Bruno Fardim Christo, Felipe Schwerz and Henrique Boriolo Dias**, revised 10 November 2026 for *Agricultural and Forest Meteorology*.

The study covers 2,146 administrative units in Argentina and Brazil during 1974–2024. It separates mean ENSO responses, temporal recurrence, conditional predictive skill and spatial aggregation. This is a manuscript-revision archive, not a claim of an accepted publication or an operational seasonal forecast service.

## What this release reproduces

| Component | Available here | Reproduction boundary |
|---|---|---|
| ENSO indices and classification | Frozen monthly/seasonal series, annual windows and standard-library classifier | Reproduces all 51 harvest assignments: 17 El Niño, 19 La Niña and 15 neutral |
| Agricultural/ENSO analysis | Quality-controlled yields, reconstructed local panel, archived targets, predictions, folds, tuning, original R scripts and Python replay/sensitivity code | ENSO replay and additional linear-trend validation can run offline |
| Revision analyses | Antecedent-ENSO, continuous-water and matched-climate code; returned predictions, diagnostics, repair logs and audits | Returned results can be audited offline; the full water/climate inputs are supplied in the release assets listed below |
| Figures | Final main and supplementary figures, original R figure scripts, scientific diagnostic code and reporting-boundary code/data | Figure versions and limitations are described below |

The full original scientific project was subsequently supplied by the corresponding author on 10 November 2026 and is included as **three release assets**: 2,150 files (1.61 GB before compression). This includes the daily NASA POWER caches and the full grid/year climate panel that were absent from the initial revision archive. Download and restore these separately using `restore_original_data.py` below.

The original R scripts are preserved, but their presence does not imply that every original seven-index hierarchical climate/PCA stage was rerun for this public release. Some archived exploratory figure scripts use earlier predictor sets; do not substitute those outputs for the final seven-index figures. Original inputs/figures and final revision outputs are kept separately to preserve provenance.

## Quick start

Use Python 3.12 or 3.13. Install the numerical dependencies in a separate environment:

```bash
git clone https://github.com/rogeriosouzanoia/soybean-enso-south-america.git
cd soybean-enso-south-america
python -m venv .venv
# Activate .venv using the command for your operating system.
python -m pip install -r requirements.txt
python unpack_data.py
python enso/reproduce_classification.py
python revision/verify_revision.py
python validation/audit_recovered_data.py
python revision/audit_returned_results.py
```

The ZIP files are ordinary versioned research data, not opaque executables. `unpack_data.py` verifies their SHA-256 digests, rejects unsafe member paths and restores the expected directory structure. Download all repository files (clone or GitHub's Download ZIP), rather than downloading a single script. Audits write their reports into the corresponding results folders.

### Additional model reruns

```bash
python validation/run_validation.py
python revision/run_lifecycle.py
```

These commands refit models and may take substantially longer. Numerical agreement can depend on Python/NumPy/BLAS versions. Historical environment records are included; publication checks for this release are recorded separately in `PUBLICATION_CHECKS.json`.

To stage recovered inputs for the original R workflow:

```bash
python validation/restore_enso_project.py
```

The original R scripts are in `validation/archived_R_scripts/`. In R, set `Sys.setenv(ENSO_SOYBEAN_PROJECT_DIR = "/absolute/path/to/enso_project")`, or pass the project path as the first argument for scripts that support command-line arguments. The initial download script now also honors this environment variable. Package requirements are listed in `R_DEPENDENCIES.txt`; an exact original R lockfile was not available. Re-downloading live data may produce different source values from the frozen analysis inputs.

### Restore the original data and rerun water/climate analyses

The [versioned release](https://github.com/rogeriosouzanoia/soybean-enso-south-america/releases/tag/v2026.11.10) includes:

- `original_project_core_2026-11-10.zip`: original agricultural inputs, detrending outputs, ENSO records, climate-index panels, model outputs, metadata and historical figures.
- `original_power_daily_grid_2026-11-10.zip`: 882 daily meteorological grid files.
- `original_power_et0_inputs_2026-11-10.zip`: 882 active ET0 input files and three pre-correction backups.

```bash
python restore_original_data.py --download
python revision/run_water_memory.py --project original_project --output reruns/water
python revision/run_climate_revision.py --project original_project --water reruns/water/water_initialization_local.csv.gz --output reruns/climate
```

The restore command uses only the Python standard library, checks archive and member SHA-256 values and restores the original directory structure under `original_project/`. It refuses to overwrite modified files. For a manual download, use `python restore_original_data.py --archive-dir /path/to/downloaded/assets`. Allow space for the compressed assets and the 1.61 GB extracted project. The complete original weather cache covers 1981 through July 2024; the revision requires 1984-01-01 through 2024-07-31, audited here for all eight used weather variables.

`ORIGINAL_DATA_MANIFEST.json` records all 2,150 files and their source/public hashes. `ORIGINAL_DATA_CHECKS.json` records complete expected weather coverage without missing/invalid required values and agreement of all 82,259 original/recovered agricultural rows after normalizing administrative-code padding. Only three metadata files had personal computer paths replaced; numerical values were not changed. The three ET0 backups preserve the pre-correction values. The four original off-season radiation estimates recorded in the revision remain estimates, not NASA observations.

The optional `revision/EXECUTAR_REVISAO.R` wrapper runs the preflight and numerical workflow through reticulate. Set `ENSO_SOYBEAN_PROJECT_DIR` to the restored `original_project` and `ENSO_REVISION_DIR` to this repository's `revision` directory. It may install dependencies and repair eligible isolated off-season radiation gaps, preserving backups and logs; read its behavior before running it on a different cache. The released active cache already passes the required-variable audit. These checks establish input availability and consistency; they are not a new end-to-end R or climate-model rerun.

## ENSO definition

Phase requires ONI ≥ +0.50 °C or ≤ −0.50 °C for at least five consecutive overlapping three-month seasons. Each harvest inherits its DJF phase. Intensity uses the maximum absolute ONI over the **whole qualifying episode**, with cutoffs 1.00, 1.50 and 2.00 °C; strong and very strong are pooled for display. Thresholds use the frozen two-decimal values without rounding beforehand.

The episode may peak after the focal harvest. These are retrospective episode labels, not concurrent DJF intensity or information available before sowing. Continuous regressions use actual seasonal ONI, including neutral years. The frozen ONI and monthly Niño 3.4 files correspond to ERSSTv6; monthly Niño 3.4 uses the 1981–2010 reference. See `enso/README.txt` for sources and provenance.

## Repository contents

- `validation/`: original R sources and independent ENSO replay/sensitivities; input/result data unpack from `data_validation_*.zip`.
- `revision/`: completed revision-analysis Python/R sources; inputs, archived outputs and diagnostics unpack from `data_revision_*.zip`.
- `enso/`: frozen series, annual metadata and classification code.
- `figures/`: reporting-boundary reconstruction script; final figures unpack from `figures_geometry_*.zip`.
- `geometry/`: frozen generalized Natural Earth 1:50 million boundaries, restored by the unpacker.
- `SOURCE_PROVENANCE.json`, `PUBLICATION_CHANGES.csv`: source/public hashes and explicitly documented portability edits.
- `DATA_MANIFEST.json`: inventory and checksums for the archived revision data and figures.
- `ORIGINAL_DATA_MANIFEST.json`, `ORIGINAL_DATA_CHECKS.json`, `ORIGINAL_DATA_CHANGES.json`, `ORIGINAL_WEATHER_AUDIT.csv`: original-data inventory, verification and portability changes.
- `restore_original_data.py`: verified restoration of the three large release assets.

The reporting map in Figures 1 and 3 outlines Argentina and unions of complete Brazilian states assigned to each reporting group; darker polygons identify sampled administrative units. It is not a historical soybean land-cover map. Rebuild it with `python figures/rebuild_reporting_boundaries.py` after unpacking the data. Final raster figures retain their original native resolution; TIFF metadata does not add detail. The locator also has a vector PDF. `revision/response_illustrations/` contains scientific diagnostic plots used during revision; confidential peer-review correspondence and editorial documents are not part of this research repository.

## Sources and citation

Agricultural records originate from [IBGE SIDRA table 5457](https://sidra.ibge.gov.br/tabela/5457) and Argentina's [Agricultural Estimates](https://datosestimaciones.magyp.gob.ar/). Weather originates from [NASA POWER](https://power.larc.nasa.gov/). ENSO series originate from [NOAA CPC](https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt) and [NOAA PSL](https://psl.noaa.gov/enso/mei/). Reporting boundaries originate from [Natural Earth](https://www.naturalearthdata.com/). Source-provider attribution and applicable terms remain relevant; this repository does not assign a new blanket license to third-party data.

Use `CITATION.cff` to cite this repository and record the exact Git commit used. No article DOI or repository DOI has been assigned here. Public availability on GitHub is not a claim of a permanent DOI-backed archive. Authorship and analysis responsibility remain with the manuscript authors. AI assistance in code development, computational reproduction and editorial preparation is disclosed in the manuscript.
