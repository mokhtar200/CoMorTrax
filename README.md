# CoMorTrax

**Comorbidity-Aware Multi-Scale Transcriptomic Analysis Framework**

[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue)](https://opensource.org/licenses/Artistic-2.0)
[![R CMD check](https://img.shields.io/badge/R%20CMD%20check-0%20ERROR%20%7C%200%20WARNING-brightgreen)]()

> CoMorTrax is currently under review at Bioconductor.

---

## What CoMorTrax Does

Comorbidities affect the majority of clinical patients, yet existing single-cell
tools analyze each disease in isolation. CoMorTrax is the first Bioconductor
framework dedicated to modeling **non-additive molecular interactions** that
emerge when two or more diseases co-exist in the same cell.

For each gene, CoMorTrax computes an **interaction delta**:

```
Delta = FC(Disease_A + Disease_B vs. Control)
      - [FC(Disease_A vs. Control) + FC(Disease_B vs. Control)]
```

Genes with a large positive Delta are **synergistic** (amplified beyond
expectation). Genes with a large negative Delta are **antagonistic**
(suppressed below expectation).

---

## Features

| Module | Function | Description |
|--------|----------|-------------|
| Disease Encoding | `encodeDiseaseLabels()` | Multi-label binary disease matrix |
| Comorbidity DE | `runComorbidityDE()` | DE across all pairwise + comorbid comparisons |
| **Interaction Effects** | `quantifyInteractionEffects()` | **Synergy/antagonism scoring** |
| Pathway Analysis | `runComorbidityPathways()` | Comorbidity-stratified enrichment |
| Vulnerability Scoring | `scoreCellVulnerability()` | Cell-type vulnerability ranking |
| Network Modeling | `buildComorbidityNetworks()` | GRN rewiring between disease states |
| ML Classification | `trainComorbidityClassifier()` | Comorbidity state prediction |
| Spatial Mapping | `mapSpatialComorbidity()` | Tissue-level projection |

---

## Installation

**From Bioconductor (pending review):**

```r
# Submission pending to Bioconductor.
# Once accepted:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CoMorTrax")
```

**From GitHub (development version):**

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mokhtar200/CoMorTrax")
```

---

## Example Workflow

```r
library(CoMorTrax)
set.seed(42)

# 1. Simulate a comorbidity dataset with known interaction signals
sce <- simulateCoMorTrax(
  n_cells  = 400,
  n_genes  = 300,
  n_syn    = 20,
  n_ant    = 15,
  diseases = c("AD", "T2D")
)

# 2. Create object and encode disease labels
cmt <- createCoMorTraxObject(sce) |>
  encodeDiseaseLabels(diseaseCol = c("disease_AD", "disease_T2D"))

# 3. Comorbidity-aware differential expression
cmt <- runComorbidityDE(cmt, method = "wilcoxon")

# 4. Quantify interaction effects (core module)
cmt <- quantifyInteractionEffects(cmt, nPermutations = 1000)
cat("Synergistic genes:", length(getSynergisticGenes(cmt)), "\n")
cat("Antagonistic genes:", length(getAntagonisticGenes(cmt)), "\n")

# 5. Downstream analyses
cmt <- scoreCellVulnerability(cmt)
cmt <- buildComorbidityNetworks(cmt)
cmt <- trainComorbidityClassifier(cmt, method = "randomForest")

# 6. Summarise and export
summaryCoMorTrax(cmt)
exportCoMorTrax(cmt, outDir = "./results")
```

---

## Citation

If you use CoMorTrax in your research, please cite:

```
Ahmed Mokhtar (2026). CoMorTrax: A Comorbidity-Aware Multi-Scale
Transcriptomic Analysis Framework for Single-Cell Disease Modeling.
Bioconductor package version 0.99.0.
https://github.com/mokhtar200/CoMorTrax
```

BibTeX:

```bibtex
@misc{Mokhtar2026CoMorTrax,
  author  = {Ahmed Mokhtar},
  title   = {{CoMorTrax}: A Comorbidity-Aware Multi-Scale Transcriptomic
             Analysis Framework},
  year    = {2026},
  note    = {Bioconductor package version 0.99.0, submission pending},
  url     = {https://github.com/mokhtar200/CoMorTrax}
}
```

---

## License

This package is licensed under the [Artistic License 2.0](LICENSE).

---

## Author

**Ahmed Mokhtar**
Bioinformatics Research Group
Email: ahmed.salem@bioinformatics.org
GitHub: https://github.com/mokhtar200

---

## Links

- GitHub repository: https://github.com/mokhtar200/CoMorTrax
- Bug reports: https://github.com/mokhtar200/CoMorTrax/issues
- Bioconductor submission: pending review
