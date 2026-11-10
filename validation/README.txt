HISTORICAL RECORD — 14 October 2026. See ../README.md (repository root README) for the current public release and scope. Some archived paths were made portable; see PUBLICATION_CHANGES.csv. Statements below describe the October recovery, not the current deposit.

ENSO SOYBEAN REPRODUCIBILITY AND ADDITIONAL VALIDATION
Manuscript: Soybean yield responses to ENSO in South America
Prepared 14 October 2026. Corresponds to Supplement S10 and Table S13.

COMPLETED
Recovered quality-controlled annual yields: 82,259 rows, 263 Argentine departments and 1,883 Brazilian municipalities. Reconstructed the exact prepared ENSO modelling targets from archived local anomalies plus original yields and frozen NOAA indices. Reproduced regional targets, the primary field ridge-model fitting/tuning/predictions, and all 4,292 local phase means/medians from original yields. Added three linear-trend sensitivities, including preprocessing repeated within validation splits and past-only annual tests for 1995-2024. All 410 validation checks and 8 recovered-data checks passed.

RUN THE INDEPENDENT ANALYSES
Use Python 3.12 with the pinned NumPy and pandas versions in requirements.txt.
1. python -m pip install -r requirements.txt
2. python run_validation.py
3. python audit_recovered_data.py
Results are written to results/. The scripts require no network connection. The frozen protocol is protocol.json. Raw agricultural downloads and original R computations are not repeated by these Python scripts.

PRIMARY MODEL VERSUS SENSITIVITIES
published_replay refits the original field model to its archived prepared regional targets, reproducing all held-out predictions to about 1e-14 percentage points. Its bootstrap intervals are independently recomputed with a different seed and do not replace Table S8.
linear_loo_control uses full-record leave-one-year-out linear trends and the original episode folds.
linear_nested re-estimates a prespecified linear trend inside every outer and inner episode split.
linear_forward trains only through t-2 for each test harvest t in 1995-2024, with three past-only inner blocks, and training-derived regional reference areas.
The last three designs use linear trends, not the original selection among eight trend candidates. Their target transformations, periods and weighting can differ; score differences are not isolated estimates of leakage bias. All candidate and baseline losses use common finite support. Observed contemporaneous DJF ONI, realized within-region harvested-area weighting and the fixed historical cohort mean these are retrospective tests, not forecasts issued before the season.

CONTENTS
inputs/: frozen quality-controlled yield CSV; archived field predictions; original local prepared anomalies extracted without duplicate model rows; frozen NOAA series; a reconstructed merged local ENSO panel (gzip); original phase summaries and spatial/calendars inventory.
results/: all held-out predictions, annual error components, scores/intervals, every inner-tuning candidate, outer/inner year memberships, computational checks, numerical replay comparisons and measured environment versions.
archived_R_scripts/: recovered original analysis/figure scripts, including source-download, harmonization and detrending stages. These retain their original defaults and dependencies. Their presence does not imply that all original R stages were rerun.
climate_recovery/: recovered AWC100 water-balance records, country climate summaries, crop calendars and original NASA request/file manifests. Country summaries are not local climate inputs. Original local Windows paths in manifests are provenance, not accessible files in this archive.
restore_enso_project.py: stages the recovered local ENSO panel, yield input, NOAA series, calendars and R sources into a conventional enso_project directory. It does not recreate the missing local climate panel or download meteorology.
file_inventory.csv: sizes and SHA-256 hashes of the packaged files, excluding this inventory itself.

REMAINING CLIMATE RELEASE
The local soybean_climate_indices_grid_year.csv.gz and exact seven-index climate-model/PCA outputs were not found in the available files. The older figure_data.zip and several saved climate tables use eleven predictors or five mechanism groups. They must not be relabelled as outputs from the seven-index model. NASA API retrieval was attempted but timed out in this environment. Complete the exact local weather/climate panel and matching seven-index model outputs before claiming complete article reproducibility or finalizing a public data/code accession. No public deposit was performed.

PROVENANCE
Yield input: soybean_Brazil_Argentina_analysis_ready_yield(1).csv, generated from the study IBGE and Argentine agricultural records. Original anomalies were extracted from ENSO_spatially_resolved_oof_predictions.csv, SHA-256 recorded in recovery_provenance.json. The merged panel is reconstructed and labelled accordingly; it is not presented as a recovered untouched original panel file. Source redistribution terms remain those of the original providers; no new data license or public accession has been assigned here.
