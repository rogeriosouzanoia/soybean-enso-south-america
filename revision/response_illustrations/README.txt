RESPONSE FIGURES

The letter contains 13 complete figures: seven new response-only diagnostics
(R1–R7), complete manuscript Figures 1, 2 and 5, and complete supplementary
Figures S1, S7 and S8. Response-only figures are not additions to the article.

R1  Lifecycle counts and national antecedent-ONI increments
R2  Original seasonal-reset CWS versus continuous storage and local changes
R3  All 32 ENSO increments by preprocessing, initialization and scoring scale
R4  Climate-index importance before and after initialization changes
R5  Original primary ENSO field-skill estimates and confidence intervals
R6  All retained years versus the ≥80%-coverage subset in Table S6
R7  Every MATOPIBA harvest contributing to the intensity-class means

The numerical analyses, main article, tracked article and supplement are
unchanged in this update. No new model fitting or uncertainty estimation was
performed. Inferential intervals are taken from the completed or archived
outputs; descriptive displays are identified explicitly in their captions.
Counts of local seasons or harvests are not treated as independent ENSO events.

REPRODUCE THE NEW FIGURES
From the package root, with its existing graphical dependencies installed:
  python response_illustrations/make_response_illustrations.py

The script reads the completed results folders and the two archived input
CSVs, plus the reported Table S6 values copied to inputs/. It writes PNG and
vector PDF exports to figures/. The source-data hashes and all captions are
in figure_manifest.json. The coverage figure in source_figures/ is the exact
unmodified image embedded in the Supplementary Material.

response_analysis_coverage.csv maps each empirical reviewer request to its
visual evidence. illustration_verification.json records the numerical and
document checks. The letter was rendered and all 23 pages reviewed.
