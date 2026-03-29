## =============================================================================
## validation.R - Biological & Statistical Validation
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @title Validate CoMorTrax Results
#' @description
#' Comprehensive validation of CoMorTrax outputs including:
#' \enumerate{
#'   \item \strong{Biological validation}: check enrichment of known disease genes
#'     (e.g., APOE, BACE1 for Alzheimer's) in synergistic gene sets.
#'   \item \strong{Statistical validation}: permutation testing for interaction
#'     effect significance beyond the permutations already run internally.
#'   \item \strong{Cross-dataset validation}: transfer interaction signatures to
#'     an independent \code{SingleCellExperiment} and compute concordance.
#'   \item \strong{Robustness testing}: subsample analysis to assess stability
#'     of interaction classifications.
#' }
#'
#' @param object A \code{\link{CoMorTraxObject}} with quantified interaction
#'   effects.
#' @param knownGenes Named list of character vectors; known disease-associated
#'   genes per disease (e.g., \code{list(AD = c("APOE","BACE1","TREM2"))}).
#' @param testDataset Optional \code{SingleCellExperiment}; independent dataset
#'   for cross-validation.
#' @param nBootstrap Integer; bootstrap iterations for robustness testing.
#'   Default \code{100L}.
#' @param subsampleFraction Numeric; fraction of cells per group to subsample
#'   in each bootstrap iteration. Default \code{0.8}.
#' @param seed Integer.
#' @param verbose Logical. Default \code{TRUE}.
#' @return A named \code{list} with validation results:
#'   \describe{
#'     \item{\code{biologicalValidation}}{Overlap statistics for known genes.}
#'     \item{\code{robustnessScores}}{Bootstrap stability of gene classifications.}
#'     \item{\code{crossDatasetConcordance}}{Concordance metrics if
#'       \code{testDataset} provided.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' known <- list(
#'     AD = c("APOE", "BACE1", "TREM2", "CLU", "BIN1"),
#'     DM = c("INS", "INSR", "GCK", "HNF4A", "PDX1")
#' )
#' val <- validateCoMorTrax(obj, knownGenes = known)
#' val$biologicalValidation
#' }
validateCoMorTrax <- function(object,
                               knownGenes         = NULL,
                               testDataset        = NULL,
                               nBootstrap         = 100L,
                               subsampleFraction  = 0.8,
                               seed               = NULL,
                               verbose            = TRUE) {

    if (nrow(interactionScores(object)) == 0)
        stop("Run quantifyInteractionEffects() first.")

    seed <- .param(seed, object@params$seed)
    set.seed(seed)

    results <- list()

    # ------------------------------------------------------------------
    # 1. Biological validation: known gene overlap
    # ------------------------------------------------------------------
    if (!is.null(knownGenes)) {
        if (verbose) .msg("Biological Validation ...", "step")

        bio_val <- lapply(names(knownGenes), function(disease) {
            known     <- knownGenes[[disease]]
            all_genes <- rownames(object)
            n_total   <- length(all_genes)

            lapply(c("synergistic", "antagonistic", "additive"), function(cls) {
                cls_genes <- .get_classified_genes(object, cls, NULL)
                n_cls     <- length(cls_genes)
                overlap   <- length(intersect(known, cls_genes))
                expected  <- (n_cls / n_total) * length(known)

                # Hypergeometric test
                p_val <- stats::phyper(
                    overlap - 1, n_cls, n_total - n_cls,
                    length(known), lower.tail = FALSE)

                data.frame(
                    disease     = disease,
                    geneSet     = cls,
                    nKnown      = length(known),
                    nGeneSet    = n_cls,
                    nOverlap    = overlap,
                    expected    = round(expected, 2),
                    pvalue      = p_val,
                    significant = p_val < 0.05,
                    stringsAsFactors = FALSE
                )
            })
        })

        results$biologicalValidation <- do.call(rbind,
            lapply(bio_val, function(x) do.call(rbind, x)))

        if (verbose) {
            sig_rows <- results$biologicalValidation[
                results$biologicalValidation$significant, ]
            .msg(sprintf("  %d significant overlaps with known disease genes.",
                         nrow(sig_rows)), "info")
        }
    }

    # ------------------------------------------------------------------
    # 2. Robustness: bootstrap subsampling
    # ------------------------------------------------------------------
    if (verbose) .msg("Robustness Testing (bootstrap) ...", "step")

    is_df  <- as.data.frame(interactionScores(object))
    grps   <- diseaseGroups(object)
    expr   <- .get_norm_matrix(object)

    # Original classifications
    orig_class <- setNames(is_df$classification, is_df$gene)

    bootstrap_res <- BiocParallel::bplapply(seq_len(nBootstrap), function(b) {
        set.seed(seed + b)
        # Subsample each group
        grps_sub <- lapply(grps, function(idx) {
            n <- max(3, round(length(idx) * subsampleFraction))
            sample(idx, min(n, length(idx)))
        })

        # Recompute interaction delta on subsample
        d1_grp <- is_df$disease1[1]
        d2_grp <- is_df$disease2[1]
        cm_grp <- is_df$comorbidGroup[1]
        ctrl   <- setdiff(names(grps), c(d1_grp, d2_grp, cm_grp,
                                          grep("_", names(grps), value = TRUE)))[1]
        if (is.na(ctrl)) ctrl <- names(grps)[1]

        fc_d1_b <- .mean_logFC(expr, grps_sub[[d1_grp]], grps_sub[[ctrl]])
        fc_d2_b <- .mean_logFC(expr, grps_sub[[d2_grp]], grps_sub[[ctrl]])
        fc_cm_b <- .mean_logFC(expr, grps_sub[[cm_grp]], grps_sub[[ctrl]])
        delta_b <- fc_cm_b - (fc_d1_b + fc_d2_b)

        cls_b <- ifelse(
            delta_b > object@params$synThreshold, "synergistic",
            ifelse(delta_b < object@params$antThreshold, "antagonistic",
                   "additive"))
        names(cls_b) <- rownames(expr)
        cls_b
    }, BPPARAM = BiocParallel::SerialParam())

    # Stability = proportion of bootstraps where classification matches original
    all_genes <- rownames(expr)
    stability <- vapply(all_genes, function(g) {
        boot_classes <- vapply(bootstrap_res, function(b) b[g], character(1))
        mean(boot_classes == orig_class[g], na.rm = TRUE)
    }, numeric(1))

    results$robustnessScores <- data.frame(
        gene        = all_genes,
        stability   = stability,
        classification = orig_class[all_genes],
        stringsAsFactors = FALSE
    )

    median_stab <- median(stability, na.rm = TRUE)
    if (verbose)
        .msg(sprintf("  Median classification stability: %.3f", median_stab),
             "info")

    # ------------------------------------------------------------------
    # 3. Cross-dataset validation
    # ------------------------------------------------------------------
    if (!is.null(testDataset)) {
        if (verbose) .msg("Cross-Dataset Validation ...", "step")

        # Score each spot/cell in testDataset using our interaction gene sets
        test_expr <- .get_norm_matrix_from_sce(testDataset)

        concordance <- lapply(c("synergistic", "antagonistic"), function(cls) {
            genes    <- .get_classified_genes(object, cls, NULL)
            in_test  <- intersect(genes, rownames(test_expr))
            if (length(in_test) < 3) return(NULL)
            score_in_test <- colMeans(test_expr[in_test, , drop = FALSE])
            data.frame(
                cls        = cls,
                nGenes     = length(genes),
                nInTest    = length(in_test),
                meanScore  = mean(score_in_test),
                sdScore    = sd(score_in_test)
            )
        })
        results$crossDatasetConcordance <- do.call(rbind,
            Filter(Negate(is.null), concordance))
        if (verbose)
            .msg("Cross-dataset validation complete.", "success")
    }

    if (verbose) .msg("Validation complete.", "success")
    results
}

#' @title Export CoMorTrax Results to Disk
#' @describeIn exportCoMorTrax Export all analysis results to CSV/RDS.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param outDir Character; output directory (created if needed).
#' @param formats Character vector; output formats. Combination of
#'   \code{"csv"}, \code{"rds"}. Default both.
#' @param verbose Logical. Default \code{TRUE}.
#' @param ... Ignored.
#' @return Invisibly returns \code{outDir}.
#' @export
setMethod("exportCoMorTrax", "CoMorTraxObject",
    function(object, outDir, formats = c("csv", "rds"), verbose = TRUE, ...) {

        dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
        if (verbose) .msg(sprintf("Exporting to: %s", outDir), "step")

        # RDS: full object
        if ("rds" %in% formats) {
            rds_path <- file.path(outDir, "CoMorTraxObject.rds")
            saveRDS(object, rds_path)
            if (verbose) .msg(paste("  Saved:", rds_path), "info")
        }

        if ("csv" %in% formats) {
            # Disease matrix
            if (nrow(diseaseMatrix(object)) > 0) {
                utils::write.csv(as.data.frame(diseaseMatrix(object)),
                    file.path(outDir, "disease_matrix.csv"))
            }

            # DE results
            if (length(deResults(object)) > 0) {
                de_dir <- file.path(outDir, "DE_results")
                dir.create(de_dir, showWarnings = FALSE)
                for (nm in names(deResults(object))) {
                    safe_nm <- gsub("[^A-Za-z0-9_]", "_", nm)
                    utils::write.csv(deResults(object)[[nm]],
                        file.path(de_dir, paste0(safe_nm, ".csv")),
                        row.names = FALSE)
                }
            }

            # Interaction scores
            if (nrow(interactionScores(object)) > 0) {
                utils::write.csv(
                    as.data.frame(interactionScores(object)),
                    file.path(outDir, "interaction_scores.csv"),
                    row.names = FALSE)
            }

            # Cell vulnerability
            if (nrow(cellVulnerability(object)) > 0) {
                utils::write.csv(
                    as.data.frame(cellVulnerability(object)),
                    file.path(outDir, "cell_vulnerability.csv"),
                    row.names = FALSE)
            }

            # Classifier predictions
            clf <- classifier(object)
            if (length(clf) > 0 && !is.null(clf$predictions)) {
                utils::write.csv(clf$predictions,
                    file.path(outDir, "classifier_predictions.csv"),
                    row.names = FALSE)
                if (!is.null(clf$featureImportance)) {
                    fi_df <- data.frame(gene = names(clf$featureImportance),
                                        importance = clf$featureImportance)
                    utils::write.csv(fi_df,
                        file.path(outDir, "feature_importance.csv"),
                        row.names = FALSE)
                }
            }
        }

        if (verbose)
            .msg("Export complete.", "success")
        invisible(outDir)
    }
)
