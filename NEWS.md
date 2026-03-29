# CoMorTrax 0.99.0 (2026-03-23)

## Initial Bioconductor Submission

- Initial Bioconductor submission
- Full comorbidity-aware pipeline implemented across 8 modules
- Multi-label disease encoding with automatic comparison enumeration
- Comorbidity-aware differential expression (Wilcoxon, edgeR, MAST)
- Interaction scoring module: permutation-based synergy/antagonism detection
- Pathway-level comorbidity enrichment (GO, KEGG, Reactome)
- Cell-type vulnerability composite scoring
- Network rewiring analysis between disease states
- ML classification layer (Random Forest, LASSO, Elastic Net)
- Spatial transcriptomics mapping of comorbidity scores
- Simulation engine with injected synergistic/antagonistic signals
- Multi-layer validation (biological, robustness, cross-dataset)
- Unified visualization interface with 10 plot types
- Structured result export (CSV + RDS)

# CoMorTrax 0.99.1 (2026-03-29)

## New Features

- Added `CoMorTrax_AD_T2D_example` bundled dataset: compact 400-cell x 150-gene
  example representative of Mathys et al. 2019 (AD, GSE125050) and Segerstolpe
  et al. 2016 (T2D, E-MTAB-5061).
- Added `inst/extdata/real_dataset_analysis.R`: complete code to download and
  process 5 published GEO datasets (GSE125050, E-MTAB-5061, GSE157783,
  GSE174188, GSE213380).
- Added `inst/extdata/real_dataset_benchmark.txt`: expected results on all 5
  real datasets with GEO accessions and literature references.
- Updated DESCRIPTION: LazyData=true, corrected email, expanded Description.
- Updated README.md: Bioconductor-style with real dataset benchmark table.

## Bug Fixes

- Fixed `synThreshold`/`antThreshold` params inheritance from object store
  (was using old 0.5 default instead of new 0.25, causing recall=0.275).
- Fixed `.is_binary_col()` false positive on string disease columns.
- Fixed 3-disease simulation to create all pairwise comorbid groups.
- Fixed adaptive correlation threshold in `buildComorbidityNetworks()`.
