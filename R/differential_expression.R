## =============================================================================
## differential_expression.R - Comorbidity-Aware DE Engine (Module 2)
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @describeIn runComorbidityDE Run comorbidity-aware differential expression.
#'
#' @description
#' \strong{Module 2: Comorbidity-Aware Differential Expression Engine}
#'
#' Extends standard single-disease DE analysis to model ALL pairwise and
#' multi-way disease comparisons simultaneously. For diseases AD and DM,
#' this produces comparisons:
#' \enumerate{
#'   \item AD_only vs control
#'   \item DM_only vs control
#'   \item AD_DM vs control
#'   \item AD_DM vs AD_only
#'   \item AD_DM vs DM_only
#' }
#' Results are stratified by cell type (via \code{cellTypeCol}) when provided.
#'
#' @param object A \code{\link{CoMorTraxObject}} with encoded disease labels.
#' @param cellTypeCol Character or \code{NULL}; column in \code{colData}
#'   holding cell-type annotations. If \code{NULL}, all cells are treated
#'   as one population.
#' @param controlLabel Character; group label for control cells.
#' @param method Character; DE method. One of \code{"wilcoxon"} (default,
#'   via \pkg{limma}/\pkg{scran}), \code{"MAST"}, \code{"edgeR"},
#'   \code{"DESeq2"}.
#' @param fdrThreshold Numeric; FDR threshold for significance calls.
#' @param logFCThreshold Numeric; minimum absolute log2 fold-change.
#' @param minCellsDE Integer; minimum cells required in each group to
#'   attempt DE. Default \code{10L}.
#' @param customComparisons Named list; optionally supply explicit comparisons
#'   as \code{list(name = c(numerator, denominator))} to override the
#'   automatic enumeration.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#' @param verbose Logical. Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return An updated \code{\link{CoMorTraxObject}} with
#'   \code{deResults(object)} populated. Each element is a \code{data.frame}
#'   with columns: \code{gene}, \code{logFC}, \code{pvalue}, \code{padj},
#'   \code{meanExpr_num}, \code{meanExpr_denom}, \code{significant},
#'   \code{direction}, \code{cellType}, \code{comparison}.
#'
#' @author Ahmed Mokhtar Ramzy Salem
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom SingleCellExperiment logcounts
#' @importFrom stats wilcox.test p.adjust
#' @examples
#' \dontrun{
#' obj <- runComorbidityDE(obj, cellTypeCol = "cellType",
#'                         method = "wilcoxon", fdrThreshold = 0.05)
#' }
setMethod("runComorbidityDE", "CoMorTraxObject",
    function(object,
             cellTypeCol      = NULL,
             controlLabel     = "control",
             method           = "wilcoxon",
             fdrThreshold     = NULL,
             logFCThreshold   = NULL,
             minCellsDE       = 5L,
             customComparisons = NULL,
             BPPARAM          = BiocParallel::SerialParam(),
             verbose          = TRUE,
             ...) {

        .check_disease_encoded(object)
        p <- object@params
        fdrThreshold   <- .param(fdrThreshold,   .param(p$fdrThreshold,   0.05))
        logFCThreshold <- .param(logFCThreshold, .param(p$logFCThreshold, 0.5))

        if (verbose) .msg("Module 2: Comorbidity-Aware DE Analysis ...", "step")

        cd    <- SummarizedExperiment::colData(object)
        grps  <- diseaseGroups(object)

        # Build comparison list
        if (!is.null(customComparisons)) {
            comps <- customComparisons
        } else {
            comp_df <- listComorbidityComparisons(object, controlLabel)
            comps   <- setNames(
                lapply(seq_len(nrow(comp_df)), function(i)
                    c(comp_df$numerator[i], comp_df$denominator[i])),
                comp_df$comparison
            )
        }

        # Cell-type stratification
        if (!is.null(cellTypeCol)) {
            if (!cellTypeCol %in% colnames(cd))
                stop("cellTypeCol '", cellTypeCol, "' not in colData.")
            cell_types <- levels(factor(cd[[cellTypeCol]]))
        } else {
            cell_types <- "all_cells"
        }

        # Expression matrix
        expr_mat <- .get_norm_matrix(object)

        all_results <- list()

        for (ct in cell_types) {
            for (comp_name in names(comps)) {
                num_grp  <- comps[[comp_name]][1]
                den_grp  <- comps[[comp_name]][2]

                # Cell indices
                if (!is.null(cellTypeCol)) {
                    ct_cells  <- which(cd[[cellTypeCol]] == ct)
                    num_cells <- intersect(grps[[num_grp]], ct_cells)
                    den_cells <- intersect(grps[[den_grp]], ct_cells)
                } else {
                    num_cells <- grps[[num_grp]]
                    den_cells <- grps[[den_grp]]
                }

                # Check group sizes
                if (length(num_cells) < minCellsDE ||
                    length(den_cells) < minCellsDE) {
                    if (verbose)
                        message(sprintf(
                            "  Skipping %s | %s - insufficient cells (%d / %d)",
                            ct, comp_name,
                            length(num_cells), length(den_cells)))
                    next
                }

                res_key <- paste0(ct, "__", comp_name)

                if (verbose)
                    .msg(sprintf("  %s | %s ...", ct, comp_name), "info")

                de_res <- switch(method,
                    wilcoxon = .de_wilcoxon(
                        expr_mat, num_cells, den_cells,
                        fdrThreshold, logFCThreshold),
                    edgeR    = .de_edgeR(
                        object, num_cells, den_cells,
                        fdrThreshold, logFCThreshold),
                    MAST     = .de_MAST(
                        object, num_cells, den_cells,
                        fdrThreshold, logFCThreshold),
                    stop("method must be 'wilcoxon', 'edgeR', or 'MAST'.")
                )

                de_res$cellType   <- ct
                de_res$comparison <- comp_name
                all_results[[res_key]] <- de_res
            }
        }

        deResults(object) <- all_results
        n_total <- sum(vapply(all_results, nrow, integer(1)))
        n_sig   <- sum(vapply(all_results, function(r)
            sum(r$significant, na.rm = TRUE), integer(1)))

        if (verbose)
            .msg(sprintf("DE complete: %d comparisons, %d significant genes.",
                         length(all_results), n_sig), "success")
        object
    }
)

## --- Internal DE engines -----------------------------------------------------

#' @keywords internal
.de_wilcoxon <- function(expr_mat, num_idx, den_idx,
                         fdrThreshold, logFCThreshold) {
    genes    <- rownames(expr_mat)
    n_genes  <- length(genes)
    num_mat  <- expr_mat[, num_idx, drop = FALSE]
    den_mat  <- expr_mat[, den_idx, drop = FALSE]

    logFC    <- numeric(n_genes)
    pvals    <- numeric(n_genes)
    mean_num <- numeric(n_genes)
    mean_den <- numeric(n_genes)

    for (i in seq_len(n_genes)) {
        x <- num_mat[i, ]
        y <- den_mat[i, ]
        mean_num[i] <- mean(x)
        mean_den[i] <- mean(y)
        logFC[i]    <- log2((mean(x) + 1e-6) / (mean(y) + 1e-6))

        # Wilcoxon rank-sum
        wt           <- suppressWarnings(stats::wilcox.test(x, y,
                            exact = FALSE, alternative = "two.sided"))
        pvals[i]     <- wt$p.value
    }

    padj <- stats::p.adjust(pvals, method = "BH")

    data.frame(
        gene        = genes,
        logFC       = logFC,
        pvalue      = pvals,
        padj        = padj,
        meanExpr_num  = mean_num,
        meanExpr_denom = mean_den,
        significant = (padj < fdrThreshold &
                       abs(logFC) >= logFCThreshold),
        direction   = ifelse(logFC > 0, "up", "down"),
        stringsAsFactors = FALSE
    )
}

#' @keywords internal
.de_edgeR <- function(object, num_idx, den_idx,
                      fdrThreshold, logFCThreshold) {
    .require_pkg("edgeR")
    counts_mat <- as.matrix(
        SummarizedExperiment::assay(object, "counts")[,
            c(num_idx, den_idx)])
    group <- factor(c(rep("num", length(num_idx)),
                      rep("den", length(den_idx))))
    dge   <- edgeR::DGEList(counts = counts_mat, group = group)
    dge   <- edgeR::filterByExpr(dge)
    dge   <- edgeR::normLibSizes(dge)
    design <- stats::model.matrix(~group)
    dge   <- edgeR::estimateDisp(dge, design)
    fit   <- edgeR::glmQLFit(dge, design)
    qlf   <- edgeR::glmQLFTest(fit)
    tbl   <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table

    data.frame(
        gene        = rownames(tbl),
        logFC       = tbl$logFC,
        pvalue      = tbl$PValue,
        padj        = tbl$FDR,
        meanExpr_num  = NA_real_,
        meanExpr_denom = NA_real_,
        significant = (tbl$FDR < fdrThreshold &
                       abs(tbl$logFC) >= logFCThreshold),
        direction   = ifelse(tbl$logFC > 0, "up", "down"),
        stringsAsFactors = FALSE
    )
}

#' @keywords internal
.de_MAST <- function(object, num_idx, den_idx,
                     fdrThreshold, logFCThreshold) {
    .require_pkg("MAST")
    expr_mat <- .get_norm_matrix(object)
    sub_mat  <- expr_mat[, c(num_idx, den_idx)]
    cond     <- factor(c(rep("num", length(num_idx)),
                         rep("den", length(den_idx))))
    fdat     <- data.frame(primerid = rownames(sub_mat),
                           row.names = rownames(sub_mat))
    cdat     <- data.frame(wellKey = colnames(sub_mat),
                           cond    = cond,
                           row.names = colnames(sub_mat))
    sca      <- MAST::FromMatrix(sub_mat, cdat, fdat)
    zlmfit   <- MAST::zlm(~cond, sca)
    summ     <- MAST::summary(zlmfit, doLRT = "condnum")
    summ_dt  <- summ$datatable
    fcHurdle <- merge(
        summ_dt[summ_dt$contrast == "condnum" & summ_dt$component == "H", ,
                drop = FALSE][, c("primerid", "Pr(>Chisq)")],
        summ_dt[summ_dt$contrast == "condnum" & summ_dt$component == "logFC", ,
                drop = FALSE][, c("primerid", "coef")],
        by = "primerid"
    )
    colnames(fcHurdle) <- c("gene", "pvalue", "logFC")
    fcHurdle$padj      <- stats::p.adjust(fcHurdle$pvalue, method = "BH")
    fcHurdle$significant <- (fcHurdle$padj < fdrThreshold &
                              abs(fcHurdle$logFC) >= logFCThreshold)
    fcHurdle$direction <- ifelse(fcHurdle$logFC > 0, "up", "down")
    fcHurdle$meanExpr_num  <- NA_real_
    fcHurdle$meanExpr_denom <- NA_real_
    fcHurdle
}

## --- Convenience: extract a merged DE table ----------------------------------

#' @title Get Merged DE Results Table
#' @description Combine all DE results into a single \code{data.frame}.
#' @param object A \code{\link{CoMorTraxObject}} with DE results.
#' @param significantOnly Logical; if \code{TRUE} return only significant hits.
#' @return A \code{data.frame}.
#' @export
#' @examples
#' \dontrun{
#' de_all <- getDeTable(obj)
#' de_sig <- getDeTable(obj, significantOnly = TRUE)
#' }
getDeTable <- function(object, significantOnly = FALSE) {
    if (length(deResults(object)) == 0)
        stop("No DE results found. Run runComorbidityDE() first.")
    tbl <- do.call(rbind, lapply(names(deResults(object)), function(nm) {
        d <- deResults(object)[[nm]]
        d$resultKey <- nm
        d
    }))
    rownames(tbl) <- NULL
    if (significantOnly) tbl <- tbl[tbl$significant, ]
    tbl
}
