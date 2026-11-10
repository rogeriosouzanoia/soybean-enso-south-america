Supplementary Data S1
Soybean yield responses to ENSO in South America
AGRFORMET-D-26-03447

Frozen ENSO exposures and retrospective classification

The eight input CSVs preserve the values supplied from the original project. All 40 annual ENSO fields, including phase and intensity, match the archived yield-model input for all 51 harvest years, 1974–2024. Their data values have not been replaced by a current online download.

ONI source: https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt
Series: ERSSTv6, confirmed against https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/oni/v6/
Niño 3.4 monthly source: https://psl.noaa.gov/data/correlation/nina34.anom.data
Series: ERSSTv6; anomalies relative to 1981–2010, as identified in the source-file footer.
MEI.v2 source: https://psl.noaa.gov/enso/mei/data/meiv2.data
Source descriptions checked 10 November 2026. The original exact download timestamp is not encoded in these CSVs. Source URLs are live and may subsequently change; use these frozen files to reproduce the published exposure assignments.

Run python3 reproduce_classification.py from this folder. Only the Python standard library is required. It verifies phase and intensity assignments, exports Annual_ENSO_episode_classification.csv and reports the 17 El Niño, 19 La Niña and 15 neutral harvest years. Month-centered dates label the middle month of an ONI three-month season. For example, NDJ 2015 means November 2015 to January 2016.

Phase rule: ONI at least +0.50 °C or at most −0.50 °C for five consecutive overlapping seasons; otherwise neutral. Each harvest takes its DJF phase. Intensity is the maximum absolute ONI over the full qualifying run, using the unrounded two-decimal archived values: weak <1.00, moderate <1.50, strong <2.00, very strong ≥2.00 °C. Strong and very strong are pooled for manuscript displays. The full-episode peak can occur after the focal harvest. It is not the concurrent DJF intensity or information assumed available before sowing. Repeated harvests within a run are not independent episodes. Continuous regressions retain the actual seasonal ONI values, including neutral years.

The annual metadata table records the qualifying run, its first and last seasonal center dates, all seasons reaching its absolute peak, the original strength and the pooled display class. Neutral rows have no qualifying episode or intensity peak. Raw candidate-run lengths remain available in ENSO_ONI_seasonal.csv.

ENSO_indices_harmonized_monthly.csv aligns the monthly index records. ENSO_indices_by_soybean_harvest_year_1974_2024.csv contains the annual windows. The three QC tables summarize these data. This archive documents ENSO exposure preparation; it is not a release of the complete crop, weather or fitted-model environment.
