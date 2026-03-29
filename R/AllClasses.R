## =============================================================================
## AllClasses.R - S4 Class Definitions for CoMorTrax
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @title CoMorTraxObject: Core S4 class for comorbidity analysis
#'
#' @description
#' The \code{CoMorTraxObject} is the central data container for all CoMorTrax
#' analyses. It extends \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' to store disease labels, interaction scores, network models, and classifier
#' results alongside standard single-cell data.
#'
#' @slot diseaseMatrix A \code{DataFrame} of binary disease indicators per cell.
#'   Rows = cells, columns = diseases. Values are 0 (absent) or 1 (present).
#' @slot diseaseGroups A named \code{list} mapping group labels (e.g.,
#'   "AD_only", "AD_DM", "control") to cell barcodes.
#' @slot deResults A named \code{list} of differential expression result
#'   \code{data.frame}s, one per comparison. Created by
#'   \code{\link{runComorbidityDE}}.
#' @slot interactionScores A \code{DataFrame} of per-gene interaction scores
#'   (observed, expected, delta, classification). Created by
#'   \code{\link{quantifyInteractionEffects}}.
#' @slot pathwayResults A named \code{list} of pathway enrichment results
#'   for synergistic, antagonistic, and additive gene sets. Created by
#'   \code{\link{runComorbidityPathways}}.
#' @slot cellVulnerability A \code{DataFrame} with vulnerability scores per
#'   cell type. Created by \code{\link{scoreCellVulnerability}}.
#' @slot networks A named \code{list} of \code{igraph} gene regulatory network
#'   objects, one per disease state. Created by
#'   \code{\link{buildComorbidityNetworks}}.
#' @slot networkDiff A \code{list} containing network rewiring statistics
#'   between disease states.
#' @slot classifier A named \code{list} holding trained model objects and
#'   performance metrics from \code{\link{trainComorbidityClassifier}}.
#' @slot spatialResults A \code{list} of spatial mapping results (optional).
#'   Created by \code{\link{mapSpatialComorbidity}}.
#' @slot params A named \code{list} of analysis parameters used throughout
#'   the workflow for reproducibility.
#' @slot coMorTraxVersion A \code{character} storing the package version at
#'   object creation.
#'
#' @author Ahmed Mokhtar Ramzy Salem
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @exportClass CoMorTraxObject
#'
#' @examples
#' # Create a minimal CoMorTraxObject
#' library(SingleCellExperiment)
#' counts_mat <- matrix(rpois(500, 5), nrow = 50, ncol = 10)
#' rownames(counts_mat) <- paste0("Gene", seq_len(50))
#' colnames(counts_mat) <- paste0("Cell", seq_len(10))
#' sce <- SingleCellExperiment(assays = list(counts = counts_mat))
#' obj <- createCoMorTraxObject(sce)
#' obj
setClass(
    "CoMorTraxObject",
    contains = "SingleCellExperiment",
    slots = list(
        diseaseMatrix    = "DataFrame",
        diseaseGroups    = "list",
        deResults        = "list",
        interactionScores = "DataFrame",
        pathwayResults   = "list",
        cellVulnerability = "DataFrame",
        networks         = "list",
        networkDiff      = "list",
        classifier       = "list",
        spatialResults   = "list",
        params           = "list",
        coMorTraxVersion = "character"
    ),
    prototype = list(
        diseaseMatrix    = S4Vectors::DataFrame(),
        diseaseGroups    = list(),
        deResults        = list(),
        interactionScores = S4Vectors::DataFrame(),
        pathwayResults   = list(),
        cellVulnerability = S4Vectors::DataFrame(),
        networks         = list(),
        networkDiff      = list(),
        classifier       = list(),
        spatialResults   = list(),
        params           = list(),
        coMorTraxVersion = "0.99.0"
    )
)

## -----------------------------------------------------------------------------
## Validity method
## -----------------------------------------------------------------------------

setValidity("CoMorTraxObject", function(object) {
    errors <- character()

    # Check diseaseMatrix dimensions match cell count
    dm <- object@diseaseMatrix
    if (nrow(dm) > 0 && nrow(dm) != ncol(object)) {
        errors <- c(errors,
            paste0("diseaseMatrix has ", nrow(dm), " rows but object has ",
                   ncol(object), " cells. They must match."))
    }

    # Check diseaseMatrix values are 0/1
    if (nrow(dm) > 0) {
        vals <- unlist(as.list(dm))
        vals <- vals[!is.na(vals)]
        if (!all(vals %in% c(0L, 1L, 0, 1))) {
            errors <- c(errors,
                "diseaseMatrix must contain only binary values (0 or 1).")
        }
    }

    # interactionScores: if non-empty, must have required columns
    is_df <- object@interactionScores
    if (nrow(is_df) > 0) {
        req_cols <- c("gene", "observed", "expected", "interactionDelta",
                      "classification")
        missing_cols <- setdiff(req_cols, colnames(is_df))
        if (length(missing_cols) > 0) {
            errors <- c(errors,
                paste0("interactionScores is missing required columns: ",
                       paste(missing_cols, collapse = ", ")))
        }
    }

    if (length(errors) == 0) TRUE else errors
})
