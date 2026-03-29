## =============================================================================
## spatial_mapping.R - Spatial Comorbidity Mapping (Module 8)
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @describeIn mapSpatialComorbidity Map comorbidity signatures to spatial data.
#'
#' @description
#' \strong{Module 8: Spatial Comorbidity Mapping}
#'
#' Projects comorbidity interaction signatures (synergistic, antagonistic, and
#' additive gene sets) onto spatial transcriptomics data, providing anatomical
#' context for where comorbid pathology manifests in tissue.
#'
#' The approach:
#' \enumerate{
#'   \item Compute per-spot module scores for each interaction gene set.
#'   \item Map the comorbidity vulnerability score across tissue regions.
#'   \item Overlay disease-specific signals to identify anatomically
#'     convergent zones.
#' }
#'
#' @param object A \code{\link{CoMorTraxObject}} with interaction scores.
#' @param spatialData A \code{SingleCellExperiment} or compatible object
#'   containing spatial transcriptomics data (e.g., Visium). Must include
#'   spatial coordinates in \code{colData} columns specified by
#'   \code{coordCols}.
#' @param coordCols Character vector (length 2); column names in
#'   \code{colData(spatialData)} holding x/y spatial coordinates.
#'   Default \code{c("x_coord", "y_coord")}.
#' @param scoreMethod Character; spot scoring method. One of
#'   \code{"mean"} (mean log-expression of gene set),
#'   \code{"ssGSEA"} (requires \pkg{GSVA}).
#' @param geneSetTypes Character vector; which interaction classes to project.
#'   Default: all three.
#' @param verbose Logical. Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return An updated \code{\link{CoMorTraxObject}} with \code{spatialResults}
#'   populated. Each element contains per-spot scores for one interaction class.
#'
#' @author Ahmed Mokhtar Ramzy Salem
#' @export
#' @examples
#' \dontrun{
#' obj <- mapSpatialComorbidity(obj,
#'                              spatialData = visium_sce,
#'                              coordCols = c("x_coord", "y_coord"))
#' }
setMethod("mapSpatialComorbidity", "CoMorTraxObject",
    function(object,
             spatialData   = NULL,
             coordCols     = c("x_coord", "y_coord"),
             scoreMethod   = "mean",
             geneSetTypes  = c("synergistic", "antagonistic", "additive"),
             verbose       = TRUE,
             ...) {

        if (nrow(interactionScores(object)) == 0)
            stop("Run quantifyInteractionEffects() first.")
        if (is.null(spatialData))
            stop("spatialData must be provided.")

        if (verbose)
            .msg("Module 8: Spatial Comorbidity Mapping ...", "step")

        cd_sp <- SummarizedExperiment::colData(spatialData)
        if (!all(coordCols %in% colnames(cd_sp)))
            stop("coordCols '", paste(coordCols, collapse = "', '"),
                 "' not found in spatialData colData.")

        coords <- as.data.frame(cd_sp[, coordCols])
        colnames(coords) <- c("x", "y")

        # Get expression matrix for spatial data
        sp_expr <- .get_norm_matrix_from_sce(spatialData)

        spatial_scores <- list()

        for (gs_type in geneSetTypes) {
            genes <- switch(gs_type,
                synergistic  = getSynergisticGenes(object),
                antagonistic = getAntagonisticGenes(object),
                additive     = getAdditiveGenes(object)
            )

            # Find overlapping genes
            genes_in_sp <- intersect(genes, rownames(sp_expr))
            if (length(genes_in_sp) < 3) {
                if (verbose)
                    message(sprintf("  Too few %s genes in spatial data; skipping.",
                                    gs_type))
                next
            }

            if (verbose)
                .msg(sprintf("  Scoring %d %s genes across %d spots ...",
                             length(genes_in_sp), gs_type,
                             ncol(sp_expr)), "info")

            scores <- if (scoreMethod == "mean") {
                colMeans(sp_expr[genes_in_sp, , drop = FALSE])
            } else if (scoreMethod == "ssGSEA") {
                .require_pkg("GSVA")
                gs_list <- list(geneSet = genes_in_sp)
                res <- GSVA::gsva(sp_expr, gs_list, method = "ssgsea",
                                  verbose = FALSE)
                as.numeric(res[1, ])
            } else {
                stop("scoreMethod must be 'mean' or 'ssGSEA'.")
            }

            spatial_scores[[gs_type]] <- data.frame(
                spot     = colnames(spatialData),
                x        = coords$x,
                y        = coords$y,
                score    = as.numeric(scores),
                geneSet  = gs_type,
                nGenes   = length(genes_in_sp),
                stringsAsFactors = FALSE
            )
        }

        object@spatialResults <- spatial_scores

        if (verbose)
            .msg(sprintf("Spatial mapping complete: %d gene sets projected.",
                         length(spatial_scores)), "success")
        object
    }
)

#' @keywords internal
.get_norm_matrix_from_sce <- function(sce) {
    if ("logcounts" %in% SummarizedExperiment::assayNames(sce))
        as.matrix(SingleCellExperiment::logcounts(sce))
    else
        as.matrix(SummarizedExperiment::assay(sce, 1))
}

#' @title Plot Spatial Comorbidity Map
#' @description Visualize comorbidity scores mapped onto tissue coordinates.
#' @param object A \code{\link{CoMorTraxObject}} with spatial results.
#' @param geneSet Character; which gene set to plot. One of
#'   \code{"synergistic"}, \code{"antagonistic"}, \code{"additive"}.
#' @param pointSize Numeric; spot size. Default \code{1.5}.
#' @param colors Character vector of length 2; low/high colors for the score.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#' plotSpatialComorbidity(obj, geneSet = "synergistic")
#' }
plotSpatialComorbidity <- function(object, geneSet = "synergistic",
                                    pointSize = 1.5,
                                    colors = c("#F7F7F7", "#B2182B")) {
    sp_res <- object@spatialResults
    if (length(sp_res) == 0) stop("Run mapSpatialComorbidity() first.")
    if (!geneSet %in% names(sp_res))
        stop("'", geneSet, "' not in spatialResults. Available: ",
             paste(names(sp_res), collapse = ", "))

    df <- sp_res[[geneSet]]

    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = score)) +
    ggplot2::geom_point(size = pointSize) +
    ggplot2::scale_color_gradient(low = colors[1], high = colors[2],
                                   name = "Score") +
    ggplot2::labs(
        title = paste("Spatial Comorbidity Map:", geneSet),
        subtitle = sprintf("%d genes projected", unique(df$nGenes)),
        x = "Spatial X", y = "Spatial Y") +
    ggplot2::theme_classic() +
    ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10, color = "grey40"),
        axis.text     = ggplot2::element_blank(),
        axis.ticks    = ggplot2::element_blank())
}
