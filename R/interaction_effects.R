## =============================================================================
## interaction_effects.R - Interaction Effect Quantification (Module 3)
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @describeIn quantifyInteractionEffects Quantify comorbidity interaction effects.
#'
#' @description
#' \strong{Module 3: Interaction Effect Quantification - Core Novelty}
#'
#' Detects non-additive gene expression effects under comorbidity using the
#' fundamental interaction formula:
#'
#' \deqn{
#'   \Delta_{\text{interaction}} =
#'   \text{Expr}(D_1 + D_2) - [\text{Expr}(D_1) + \text{Expr}(D_2) -
#'   \text{Expr}(\text{ctrl})]
#' }
#'
#' Genes are classified as:
#' \describe{
#'   \item{\strong{Synergistic}}{
#'     \eqn{\Delta > \theta_{\text{syn}}} - expression is \emph{amplified}
#'     beyond the sum of individual effects. These genes drive comorbidity
#'     pathology beyond what either disease alone causes.
#'   }
#'   \item{\strong{Antagonistic}}{
#'     \eqn{\Delta < \theta_{\text{ant}}} - expression is \emph{suppressed}
#'     relative to the additive expectation. These genes reflect disease
#'     cross-inhibition.
#'   }
#'   \item{\strong{Additive}}{
#'     \eqn{\theta_{\text{ant}} \le \Delta \le \theta_{\text{syn}}} -
#'     comorbid expression is the simple sum of individual disease effects.
#'   }
#' }
#'
#' Statistical significance is assessed via permutation testing.
#'
#' @param object A \code{\link{CoMorTraxObject}} with DE results.
#' @param disease1 Character; name of the first disease group
#'   (e.g., \code{"AD"}).
#' @param disease2 Character; name of the second disease group
#'   (e.g., \code{"DM"}).
#' @param comorbidGroup Character; name of the comorbid group
#'   (e.g., \code{"AD_DM"}). If \code{NULL}, inferred automatically.
#' @param controlLabel Character; control group label.
#' @param cellTypeCol Character or \code{NULL}; if provided, compute
#'   interaction scores per cell type.
#' @param synThreshold Numeric; \eqn{\Delta} threshold above which a gene is
#'   classified as synergistic (on log2FC scale). Default \code{0.5}.
#' @param antThreshold Numeric; \eqn{\Delta} threshold below which a gene is
#'   classified as antagonistic. Default \code{-0.5}.
#' @param nPermutations Integer; permutations for null distribution estimation.
#'   Default \code{1000L}.
#' @param seed Integer; random seed for permutation reproducibility.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#' @param verbose Logical. Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return An updated \code{\link{CoMorTraxObject}} with \code{interactionScores}
#'   populated. The \code{DataFrame} contains:
#'   \describe{
#'     \item{\code{gene}}{Gene identifier.}
#'     \item{\code{observed}}{Mean log2FC in the comorbid group vs control.}
#'     \item{\code{expected}}{Additive expectation: FC(D1) + FC(D2) - FC(ctrl).}
#'     \item{\code{interactionDelta}}{Observed minus expected.}
#'     \item{\code{pvalue}}{Permutation p-value.}
#'     \item{\code{padj}}{BH-adjusted p-value.}
#'     \item{\code{classification}}{One of \code{"synergistic"},
#'       \code{"antagonistic"}, \code{"additive"}.}
#'     \item{\code{cellType}}{Cell type (if \code{cellTypeCol} provided).}
#'   }
#'
#' @author Ahmed Mokhtar Ramzy Salem
#' @export
#' @importFrom stats p.adjust
#' @importFrom BiocParallel bplapply
#' @examples
#' \dontrun{
#' obj <- quantifyInteractionEffects(obj,
#'                                   disease1 = "AD",
#'                                   disease2 = "DM",
#'                                   comorbidGroup = "AD_DM",
#'                                   nPermutations = 1000)
#' }
setMethod("quantifyInteractionEffects", "CoMorTraxObject",
    function(object,
             disease1      = NULL,
             disease2      = NULL,
             comorbidGroup = NULL,
             controlLabel  = "control",
             cellTypeCol   = NULL,
             synThreshold  = NULL,
             antThreshold  = NULL,
             nPermutations = NULL,
             seed          = NULL,
             BPPARAM       = BiocParallel::SerialParam(),
             verbose       = TRUE,
             ...) {

        .check_disease_encoded(object)
        if (length(deResults(object)) == 0)
            stop("Run runComorbidityDE() first.")

        p <- object@params
        synThreshold  <- .param(synThreshold,  .param(p$synThreshold,   0.25))
        antThreshold  <- .param(antThreshold,  .param(p$antThreshold,  -0.25))
        nPermutations <- .param(nPermutations, .param(p$nPermutations, 100L))
        seed          <- .param(seed,          .param(p$seed,           42L))
        alpha         <- .param(p$interactionAlpha, 0.05)

        if (verbose)
            .msg("Module 3: Interaction Effect Quantification ...", "step")

        grp_names <- names(diseaseGroups(object))
        dm        <- as.data.frame(diseaseMatrix(object))
        cd        <- SummarizedExperiment::colData(object)

        # Infer disease and comorbid groups from available groups
        if (is.null(disease1) || is.null(disease2)) {
            # First try exact 2-disease group names
            multi_grps <- grp_names[vapply(grp_names, function(g) {
                g != controlLabel &&
                length(strsplit(g, "_")[[1]]) == 2
            }, logical(1))]

            # Fallback: any multi-disease group (3-way, etc.)
            if (length(multi_grps) == 0) {
                multi_grps <- grp_names[vapply(grp_names, function(g) {
                    g != controlLabel &&
                    length(strsplit(g, "_")[[1]]) >= 2
                }, logical(1))]
            }

            if (length(multi_grps) == 0)
                stop("No comorbid groups found. ",
                     "Please specify disease1, disease2, comorbidGroup.")

            if (is.null(comorbidGroup)) comorbidGroup <- multi_grps[1]
            parts    <- strsplit(comorbidGroup, "_")[[1]]
            disease1 <- .param(disease1, parts[1])
            disease2 <- .param(disease2, parts[2])
        }

        if (is.null(comorbidGroup))
            comorbidGroup <- paste(sort(c(disease1, disease2)), collapse = "_")

        # Find single-disease group names
        d1_grp <- .find_single_disease_group(grp_names, disease1, controlLabel)
        d2_grp <- .find_single_disease_group(grp_names, disease2, controlLabel)

        # For the comorbid group: use exact match if available, else derive
        # from disease matrix (cells positive for BOTH disease1 and disease2)
        if (comorbidGroup %in% grp_names) {
            cm_grp <- comorbidGroup
        } else if (ncol(dm) >= 2 && disease1 %in% colnames(dm) &&
                   disease2 %in% colnames(dm)) {
            # Build a virtual comorbid group from diseaseMatrix
            cm_cells <- which(dm[[disease1]] == 1 & dm[[disease2]] == 1)
            ctrl_cells <- diseaseGroups(object)[[controlLabel]]
            if (length(cm_cells) < 3)
                stop("Too few cells (", length(cm_cells),
                     ") with both ", disease1, " and ", disease2, ".")
            # Register the virtual group
            grps_new <- diseaseGroups(object)
            cm_grp   <- paste(sort(c(disease1, disease2)), collapse = "_")
            grps_new[[cm_grp]] <- cm_cells
            # Temporarily update the object's groups
            diseaseGroups(object) <- grps_new
            grp_names <- names(diseaseGroups(object))
            if (verbose)
                .msg(sprintf("  Virtual comorbid group '%s' derived from disease matrix (%d cells).",
                             cm_grp, length(cm_cells)), "info")
        } else {
            stop("Comorbid group '", comorbidGroup,
                 "' not found in diseaseGroups and cannot be derived. ",
                 "Available groups: ", paste(grp_names, collapse = ", "))
        }

        if (verbose)
            .msg(sprintf("  Computing: %s vs [%s + %s] vs %s",
                         cm_grp, d1_grp, d2_grp, controlLabel), "info")

        # Cell-type stratification
        cd <- SummarizedExperiment::colData(object)
        # Auto-detect cell type column if not provided
        if (is.null(cellTypeCol)) {
            cd_cols <- colnames(cd)
            for (guess in c("cellType", "cell_type", "CellType", "cluster",
                            "leiden_cluster", "seurat_clusters")) {
                if (guess %in% cd_cols) { cellTypeCol <- guess; break }
            }
        }

        if (!is.null(cellTypeCol)) {
            cell_types <- levels(factor(cd[[cellTypeCol]]))
        } else {
            cell_types <- "all_cells"
        }

        expr_mat <- .get_norm_matrix(object)
        all_scores <- list()

        for (ct in cell_types) {
            if (!is.null(cellTypeCol)) {
                ct_idx   <- which(cd[[cellTypeCol]] == ct)
                d1_cells <- intersect(diseaseGroups(object)[[d1_grp]], ct_idx)
                d2_cells <- intersect(diseaseGroups(object)[[d2_grp]], ct_idx)
                cm_cells <- intersect(diseaseGroups(object)[[cm_grp]], ct_idx)
                ct_cells <- intersect(diseaseGroups(object)[[controlLabel]], ct_idx)
            } else {
                d1_cells <- diseaseGroups(object)[[d1_grp]]
                d2_cells <- diseaseGroups(object)[[d2_grp]]
                cm_cells <- diseaseGroups(object)[[cm_grp]]
                ct_cells <- diseaseGroups(object)[[controlLabel]]
            }

            min_n <- min(length(d1_cells), length(d2_cells),
                         length(cm_cells), length(ct_cells))
            if (min_n < 3) {
                if (verbose)
                    message("  Skipping ", ct, " - too few cells.")
                next
            }

            # Compute mean log2FC per gene
            fc_d1 <- .mean_logFC(expr_mat, d1_cells, ct_cells)
            fc_d2 <- .mean_logFC(expr_mat, d2_cells, ct_cells)
            fc_cm <- .mean_logFC(expr_mat, cm_cells, ct_cells)

            # Interaction delta
            expected <- fc_d1 + fc_d2
            delta    <- fc_cm - expected

            # Permutation p-values
            if (verbose)
                .msg(sprintf("  Permutation testing (%d perm) ...",
                             nPermutations), "info")
            set.seed(seed)
            pvals <- .permutation_pvals(
                expr_mat, d1_cells, d2_cells, cm_cells, ct_cells,
                delta, nPermutations, BPPARAM)

            padj  <- stats::p.adjust(pvals, method = "BH")

            classification <- ifelse(
                delta > synThreshold  & padj < alpha,
                "synergistic",
                ifelse(
                    delta < antThreshold & padj < alpha,
                    "antagonistic",
                    "additive"
                )
            )

            score_df <- S4Vectors::DataFrame(
                gene              = rownames(expr_mat),
                observed          = fc_cm,
                expected          = expected,
                fc_d1             = fc_d1,
                fc_d2             = fc_d2,
                interactionDelta  = delta,
                pvalue            = pvals,
                padj              = padj,
                classification    = classification,
                cellType          = ct,
                disease1          = disease1,
                disease2          = disease2,
                comorbidGroup     = cm_grp
            )
            all_scores[[ct]] <- score_df
        }

        combined <- do.call(rbind, all_scores)
        interactionScores(object) <- combined

        tbl <- table(combined$classification)
        if (verbose) {
            .msg("  Classification summary:", "info")
            for (cls in names(tbl))
                .msg(sprintf("    %-15s: %d genes", cls, tbl[[cls]]), "info")
            .msg("Interaction effect quantification complete.", "success")
        }
        object
    }
)

## --- Helpers -----------------------------------------------------------------

#' @keywords internal
.mean_logFC <- function(expr_mat, num_idx, den_idx) {
    mean_num <- rowMeans(expr_mat[, num_idx, drop = FALSE])
    mean_den <- rowMeans(expr_mat[, den_idx, drop = FALSE])
    log2((mean_num + 1e-6) / (mean_den + 1e-6))
}

#' @keywords internal
.permutation_pvals <- function(expr_mat, d1_idx, d2_idx, cm_idx, ctrl_idx,
                                obs_delta, n_perm, BPPARAM) {

    all_idx <- c(d1_idx, d2_idx, cm_idx, ctrl_idx)
    sizes   <- c(length(d1_idx), length(d2_idx),
                 length(cm_idx), length(ctrl_idx))

    null_deltas <- BiocParallel::bplapply(seq_len(n_perm), function(i) {
        perm     <- sample(all_idx)
        p_d1     <- perm[seq_len(sizes[1])]
        p_d2     <- perm[sizes[1] + seq_len(sizes[2])]
        p_cm     <- perm[sizes[1] + sizes[2] + seq_len(sizes[3])]
        p_ctrl   <- perm[sizes[1] + sizes[2] + sizes[3] + seq_len(sizes[4])]
        p_fc_d1  <- .mean_logFC(expr_mat, p_d1, p_ctrl)
        p_fc_d2  <- .mean_logFC(expr_mat, p_d2, p_ctrl)
        p_fc_cm  <- .mean_logFC(expr_mat, p_cm, p_ctrl)
        p_fc_cm - (p_fc_d1 + p_fc_d2)
    }, BPPARAM = BPPARAM)

    null_mat <- do.call(cbind, null_deltas)   # genes x perms
    # Two-sided p-value
    vapply(seq_len(nrow(null_mat)), function(i) {
        (sum(abs(null_mat[i, ]) >= abs(obs_delta[i])) + 1) / (n_perm + 1)
    }, numeric(1))
}

#' @keywords internal
.find_single_disease_group <- function(grp_names, disease, controlLabel) {
    candidates <- grp_names[vapply(grp_names, function(g) {
        parts <- strsplit(g, "_")[[1]]
        g != controlLabel && length(parts) == 1 && parts == disease
    }, logical(1))]
    if (length(candidates) == 0)
        stop("Single-disease group for '", disease, "' not found. ",
             "Available groups: ", paste(grp_names, collapse = ", "))
    candidates[1]
}

## --- Convenience accessors ---------------------------------------------------

#' @title Get Synergistic Genes
#' @description Extract the list of synergistically expressed genes under
#'   comorbidity.
#' @param object A \code{\link{CoMorTraxObject}} with interaction scores.
#' @param cellType Character or \code{NULL}; restrict to one cell type.
#' @return Character vector of gene names.
#' @export
getSynergisticGenes <- function(object, cellType = NULL) {
    .get_classified_genes(object, "synergistic", cellType)
}

#' @title Get Antagonistic Genes
#' @description Extract the list of antagonistically expressed genes.
#' @param object A \code{\link{CoMorTraxObject}} with interaction scores.
#' @param cellType Character or \code{NULL}; restrict to one cell type.
#' @return Character vector of gene names.
#' @export
getAntagonisticGenes <- function(object, cellType = NULL) {
    .get_classified_genes(object, "antagonistic", cellType)
}

#' @title Get Additive Genes
#' @description Extract the list of additively expressed genes.
#' @param object A \code{\link{CoMorTraxObject}} with interaction scores.
#' @param cellType Character or \code{NULL}; restrict to one cell type.
#' @return Character vector of gene names.
#' @export
getAdditiveGenes <- function(object, cellType = NULL) {
    .get_classified_genes(object, "additive", cellType)
}

#' @keywords internal
.get_classified_genes <- function(object, cls, cellType) {
    df <- interactionScores(object)
    if (nrow(df) == 0)
        stop("Run quantifyInteractionEffects() first.")
    if (!is.null(cellType))
        df <- df[df$cellType == cellType, ]
    as.character(df$gene[df$classification == cls])
}
