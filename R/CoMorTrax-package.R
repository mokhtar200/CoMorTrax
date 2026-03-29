#' @title CoMorTrax: Comorbidity-Aware Multi-Omics Disease Modeling Engine
#'
#' @description
#' CoMorTrax is a Bioconductor framework that explicitly models comorbid
#' disease interactions at single-cell resolution. The package provides a
#' complete analytical pipeline from raw multi-omics data to clinically
#' interpretable comorbidity signatures.
#'
#' @section Core Modules:
#' \describe{
#'   \item{\strong{Multi-Label Disease Encoding}}{
#'     Encodes multi-disease metadata into binary disease vectors per cell,
#'     supporting any number of concurrent diseases.
#'     See \code{\link{encodeDiseaseLabels}}.
#'   }
#'   \item{\strong{Comorbidity-Aware Differential Expression}}{
#'     Extends standard DE analysis to model all pairwise and multi-way
#'     disease comparisons simultaneously.
#'     See \code{\link{runComorbidityDE}}.
#'   }
#'   \item{\strong{Interaction Effect Quantification}}{
#'     Detects non-additive gene expression effects under comorbidity,
#'     classifying genes as synergistic, antagonistic, or additive.
#'     See \code{\link{quantifyInteractionEffects}}.
#'   }
#'   \item{\strong{Pathway-Level Comorbidity Analysis}}{
#'     Performs enrichment analysis separately for synergistic and
#'     antagonistic gene sets to reveal comorbidity-specific pathways.
#'     See \code{\link{runComorbidityPathways}}.
#'   }
#'   \item{\strong{Cell-Type Vulnerability Scoring}}{
#'     Ranks cell types by their susceptibility to comorbid pathology
#'     using interaction gene counts and pathway disruption scores.
#'     See \code{\link{scoreCellVulnerability}}.
#'   }
#'   \item{\strong{Gene Regulatory Network Modeling}}{
#'     Builds and compares gene regulatory networks across disease states
#'     to identify network rewiring under comorbidity.
#'     See \code{\link{buildComorbidityNetworks}}.
#'   }
#'   \item{\strong{Machine Learning Prediction}}{
#'     Trains classifiers to predict comorbidity state from single-cell
#'     gene expression profiles.
#'     See \code{\link{trainComorbidityClassifier}}.
#'   }
#'   \item{\strong{Spatial Mapping}}{
#'     Maps comorbidity signals onto spatial transcriptomics data to
#'     provide anatomical context.
#'     See \code{\link{mapSpatialComorbidity}}.
#'   }
#' }
#'
#' @section Workflow:
#' The recommended analysis pipeline:
#' \enumerate{
#'   \item \code{\link{createCoMorTraxObject}} - construct the core object
#'   \item \code{\link{preprocessCoMorTrax}} - QC, normalization, batch correction
#'   \item \code{\link{encodeDiseaseLabels}} - multi-label disease encoding
#'   \item \code{\link{runComorbidityDE}} - comorbidity-aware DE
#'   \item \code{\link{quantifyInteractionEffects}} - interaction scoring
#'   \item \code{\link{runComorbidityPathways}} - pathway enrichment
#'   \item \code{\link{scoreCellVulnerability}} - vulnerability ranking
#'   \item \code{\link{buildComorbidityNetworks}} - network analysis
#'   \item \code{\link{trainComorbidityClassifier}} - ML prediction
#'   \item \code{\link{plotCoMorTrax}} - comprehensive visualization
#' }
#'
#' @section Author:
#' \strong{Ahmed Mokhtar Ramzy Salem}
#'
#' @references
#' Salem AMR (2024). CoMorTrax: A Comorbidity-Aware Multi-Omics Framework
#' for Modeling Disease Interactions at Single-Cell Resolution.
#' \emph{Bioinformatics}.
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/ahmedramzy/CoMorTrax}
#'   \item Report bugs: \url{https://github.com/ahmedramzy/CoMorTrax/issues}
#' }
#'
#' @docType package
#' @name CoMorTrax-package
#' @aliases CoMorTrax
"_PACKAGE"

## usethis namespace: start
#' @import methods
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @importFrom S4Vectors DataFrame SimpleList metadata metadata<-
#' @importFrom BiocGenerics counts colnames rownames
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam
#' @importFrom Matrix Matrix sparseMatrix rowSums colSums
## usethis namespace: end
NULL

#' @keywords internal
.CoMorTraxEnv <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
    .CoMorTraxEnv$version <- utils::packageVersion("CoMorTrax")
    invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
    ver <- utils::packageVersion("CoMorTrax")
    msg <- paste0(
        "\n",
        "CoMorTrax v", ver, "\n",
        "Comorbidity-Aware Multi-Omics Disease Modeling\n",
        "Author: Ahmed Mokhtar Ramzy Salem | Bioconductor Package\n"
    )
    packageStartupMessage(msg)
    invisible(NULL)
}
