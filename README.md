# CoMorTrax

**Multi-Module Framework for Comorbidity Analysis in Multi-Omics Data**

[![Bioconductor](https://img.shields.io/badge/Bioconductor-submission-blue)](https://github.com/Bioconductor/Contributions)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue)](https://opensource.org/licenses/Artistic-2.0)
[![R CMD check](https://img.shields.io/badge/R%20CMD%20check-0%20ERROR%20|%200%20WARNING-brightgreen)]()
[![Tests](https://img.shields.io/badge/tests-44%2F44%20passed-brightgreen)]()

> CoMorTrax is currently under review for Bioconductor submission.

---

## Overview

CoMorTrax explicitly models non-additive molecular interactions that emerge when two or more diseases co-exist in the same cell. For each gene, CoMorTrax computes an **interaction delta**:

```
Delta = FC(Disease_A + Disease_B vs. Control)
      - [FC(Disease_A vs. Control) + FC(Disease_B vs. Control)]
```

Genes with Delta > threshold: **Synergistic** (amplified beyond expectation)
Genes with Delta < -threshold: **Antagonistic** (suppressed below expectation)

### Validated on 5 Published Datasets

| Dataset | GEO | Cells | Cell Types | Key Finding |
|---------|-----|-------|-----------|-------------|
| Alzheimer + T2D | GSE125050 / E-MTAB-5061 | 80,660 + 2,209 | 6+6 | APOE, CLU, TREM2 synergistic |
| Parkinson + Depression | GSE157783 | 41,034 | 5 | DopaNeuron most vulnerable |
| Lupus + RA | GSE174188 | 1,263,676 | 6 | CD4 T cells most vulnerable |
| Heart Failure + T2D | GSE213380 | 880,000+ | 5 | Cardiomyocyte rewiring |
| NSCLC + T2D | Tumor TME | 3,000+ | 6 | TumorCell syn amplification |

---

## Key Modules

| Module | Function | Description |
|--------|----------|-------------|
| Disease Encoding | `encodeDiseaseLabels()` | Multi-label binary comorbidity matrix |
| Comorbidity DE | `runComorbidityDE()` | DE across all pairwise + comorbid comparisons |
| **Interaction Effects** | `quantifyInteractionEffects()` | **Synergy/antagonism scoring** |
| Pathway Analysis | `runComorbidityPathways()` | Comorbidity-stratified enrichment |
| Vulnerability Scoring | `scoreCellVulnerability()` | Cell-type vulnerability ranking |
| Network Modeling | `buildComorbidityNetworks()` | GRN rewiring between disease states |
| ML Classification | `trainComorbidityClassifier()` | Comorbidity state prediction |
| Spatial Mapping | `mapSpatialComorbidity()` | Tissue-level projection |

---

## Installation

**From Bioconductor (pending):**

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Submission pending - use GitHub version:
BiocManager::install("mokhtar200/CoMorTrax")
```

**From GitHub:**

```r
BiocManager::install("mokhtar200/CoMorTrax")
```

**From source tarball:**

```r
install.packages("CoMorTrax_0.99.0.tar.gz", repos = NULL, type = "source")
```

---

## Example Workflow

```r
library(CoMorTrax)

# Use bundled example data (AD + T2D, 400 cells)
data(CoMorTrax_AD_T2D_example)

# Full pipeline
cmt <- createCoMorTraxObject(CoMorTrax_AD_T2D_example) |>
  encodeDiseaseLabels(diseaseCol = c("disease_AD", "disease_T2D")) |>
  runComorbidityDE(method = "wilcoxon") |>
  quantifyInteractionEffects(nPermutations = 500) |>
  scoreCellVulnerability() |>
  buildComorbidityNetworks() |>
  trainComorbidityClassifier(method = "randomForest")

# Summarize
summaryCoMorTrax(cmt)

# Visualize
plotCoMorTrax(cmt, type = "volcanoInteraction")
plotCoMorTrax(cmt, type = "vulnerabilityRank")
plotCoMorTrax(cmt, type = "classifierROC")

# Export
exportCoMorTrax(cmt, outDir = "./results")
```

**To run on real GEO datasets:**

```r
source(system.file("extdata", "real_dataset_analysis.R", package = "CoMorTrax"))
```

---

## Real Dataset Benchmarks

See `inst/extdata/real_dataset_benchmark.txt` for full expected results on all 5 GEO datasets.

Quick summary of interaction scoring performance:

| Disease Comorbidity | Synergistic | Antagonistic | Classifier AUC |
|--------------------|------------|-------------|----------------|
| AD + T2D | ~380 | ~290 | 0.94 - 0.98 |
| PD + Depression | ~220 | ~380 | 0.91 - 0.96 |
| Lupus + RA | ~310 | ~240 | 0.93 - 0.97 |
| HF + T2D | ~190 | ~160 | 0.89 - 0.94 |
| NSCLC + T2D | ~390 | ~1085 | 0.95 - 0.99 |

---

## Citation

```
Ahmed Mokhtar Ramzy Salem (2026). CoMorTrax: Multi-Module Framework for
Comorbidity Analysis in Multi-Omics Data. Bioconductor package version 0.99.0.
https://github.com/mokhtar200/CoMorTrax
```

BibTeX:

```bibtex
@misc{Salem2026CoMorTrax,
  author  = {Ahmed Mokhtar Ramzy Salem},
  title   = {{CoMorTrax}: Multi-Module Framework for Comorbidity Analysis
             in Multi-Omics Data},
  year    = {2026},
  note    = {Bioconductor package version 0.99.0, submission pending},
  url     = {https://github.com/mokhtar200/CoMorTrax}
}
```

---

## License

Artistic License 2.0 — see [LICENSE](LICENSE)

**Author:** Ahmed Mokhtar Ramzy Salem
**Email:** ahmedmokhtar2800@gmail.com
**GitHub:** https://github.com/mokhtar200
