## =============================================================================
## preprocessing.R - Object Construction & Preprocessing
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

## --- createCoMorTraxObject ---------------------------------------------------

#' @describeIn createCoMorTraxObject Create from a SingleCellExperiment.
#' @param x A \code{SingleCellExperiment} object.
#' @param projectName Character; a label for this project/analysis.
#' @param params Named list of global parameters overriding defaults.
#' @param ... Ignored.
#' @return A \code{\link{CoMorTraxObject}}.
#' @export
#' @importFrom methods new validObject
#' @importFrom SummarizedExperiment colData rowData
#' @examples
#' library(SingleCellExperiment)
#' counts_mat <- matrix(rpois(200, 5), nrow = 20, ncol = 10)
#' rownames(counts_mat) <- paste0("Gene", seq_len(20))
#' colnames(counts_mat) <- paste0("Cell", seq_len(10))
#' sce <- SingleCellExperiment(assays = list(counts = counts_mat))
#' obj <- createCoMorTraxObject(sce, projectName = "TestProject")
setMethod("createCoMorTraxObject", "SingleCellExperiment",
    function(x, projectName = "CoMorTraxProject", params = list(), ...) {

        .msg("Creating CoMorTraxObject from SingleCellExperiment ...", "info")

        default_params <- list(
            projectName      = projectName,
            minCells         = 3L,
            minGenes         = 200L,
            maxMitoPercent   = 25,
            nHVGs            = 2000L,
            nPCs             = 30L,
            resolution       = 0.5,
            normMethod       = "scran",
            batchMethod      = "harmony",
            deMethod         = "wilcoxon",
            fdrThreshold     = 0.05,
            logFCThreshold   = 0.25,
            interactionAlpha = 0.05,
            synThreshold     = 0.25,
            antThreshold     = -0.25,
            nPermutations    = 1000L,
            seed             = 42L
        )

        # merge user params
        for (nm in names(params)) default_params[[nm]] <- params[[nm]]

        obj <- new("CoMorTraxObject",
            x,
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
            params           = default_params,
            coMorTraxVersion = as.character(
                utils::packageVersion("CoMorTrax"))
        )

        validObject(obj)
        .msg(sprintf("Object created: %d genes x %d cells.",
                     nrow(obj), ncol(obj)), "success")
        obj
    }
)

#' @describeIn createCoMorTraxObject Create from a raw count matrix.
#' @param x A numeric \code{matrix} or \code{dgCMatrix} of raw counts.
#'   Rows = genes, columns = cells.
#' @param ... Passed to the \code{SingleCellExperiment} method.
#' @export
setMethod("createCoMorTraxObject", "matrix",
    function(x, ...) {
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = x))
        createCoMorTraxObject(sce, ...)
    }
)

#' @describeIn createCoMorTraxObject Create from a sparse dgCMatrix.
#' @export
setMethod("createCoMorTraxObject", "dgCMatrix",
    function(x, ...) {
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = x))
        createCoMorTraxObject(sce, ...)
    }
)

## --- preprocessCoMorTrax -----------------------------------------------------

#' @describeIn preprocessCoMorTrax Preprocess a CoMorTraxObject.
#'
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param minCells Integer; minimum cells a gene must be expressed in to
#'   be retained. Default: taken from \code{object@params$minCells}.
#' @param minGenes Integer; minimum genes per cell for QC filter. Default:
#'   taken from \code{object@params$minGenes}.
#' @param maxMitoPercent Numeric [0-100]; maximum mitochondrial read percentage
#'   per cell. Default: \code{object@params$maxMitoPercent}.
#' @param normMethod Character; normalization strategy. One of
#'   \code{"scran"} (default), \code{"logNorm"}, \code{"scTransform"}.
#' @param nHVGs Integer; number of highly variable genes to select.
#' @param nPCs Integer; number of principal components.
#' @param batchVar Character or \code{NULL}; column in \code{colData} to use
#'   for batch correction. If \code{NULL}, batch correction is skipped.
#' @param batchMethod Character; batch correction method. One of
#'   \code{"harmony"} (requires the \pkg{harmony} package),
#'   \code{"fastMNN"}, \code{"none"}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#' @param verbose Logical; print progress messages. Default \code{TRUE}.
#' @param ... Ignored.
#' @return An updated \code{\link{CoMorTraxObject}} with QC metrics, normalized
#'   expression, reduced dimensions (\code{"PCA"}, \code{"UMAP"}), and
#'   cluster assignments stored in \code{colData}.
#' @export
#' @importFrom scuttle perCellQCMetrics perFeatureQCMetrics
#' @importFrom scran computeSumFactors modelGeneVar getTopHVGs
#' @importFrom scater runPCA runUMAP
#' @importFrom SingleCellExperiment reducedDim reducedDim<- logcounts logcounts<-
#' @importFrom SummarizedExperiment colData colData<-
#' @examples
#' \dontrun{
#' obj <- preprocessCoMorTrax(obj, minCells = 3, minGenes = 200,
#'                             maxMitoPercent = 25, normMethod = "scran",
#'                             nHVGs = 2000, nPCs = 30, batchVar = "batch")
#' }
setMethod("preprocessCoMorTrax", "CoMorTraxObject",
    function(object,
             minCells       = NULL,
             minGenes       = NULL,
             maxMitoPercent = NULL,
             normMethod     = NULL,
             nHVGs          = NULL,
             nPCs           = NULL,
             batchVar       = NULL,
             batchMethod    = NULL,
             BPPARAM        = BiocParallel::SerialParam(),
             verbose        = TRUE,
             ...) {

        p <- object@params
        minCells       <- .param(minCells,       p$minCells)
        minGenes       <- .param(minGenes,       p$minGenes)
        maxMitoPercent <- .param(maxMitoPercent, p$maxMitoPercent)
        normMethod     <- .param(normMethod,     p$normMethod)
        nHVGs          <- .param(nHVGs,          p$nHVGs)
        nPCs           <- .param(nPCs,           p$nPCs)
        batchMethod    <- .param(batchMethod,    p$batchMethod)

        if (verbose) .msg("Step 1: Quality Control ...", "step")

        # --- QC metrics -------------------------------------------------------
        # Detect mitochondrial genes (human: MT-, mouse: mt-)
        is_mito <- grepl("^MT-|^mt-", rownames(object))
        qc_metrics <- scuttle::perCellQCMetrics(object,
                                                subsets = list(Mito = is_mito),
                                                BPPARAM = BPPARAM)

        SummarizedExperiment::colData(object) <- cbind(
            SummarizedExperiment::colData(object), qc_metrics)

        # Filter cells
        keep_cells <- (
            qc_metrics$sum >= minGenes &
            qc_metrics$detected >= minGenes &
            (if (any(is_mito))
                qc_metrics$subsets_Mito_percent <= maxMitoPercent
             else TRUE)
        )

        n_removed_cells <- sum(!keep_cells)
        object <- object[, keep_cells]
        if (verbose)
            .msg(sprintf("  Removed %d low-quality cells. Retained: %d",
                         n_removed_cells, ncol(object)), "info")

        # Filter genes
        gene_counts <- Matrix::rowSums(
            SummarizedExperiment::assay(object, "counts") > 0)
        keep_genes <- gene_counts >= minCells
        n_removed_genes <- sum(!keep_genes)
        object <- object[keep_genes, ]
        if (verbose)
            .msg(sprintf("  Removed %d low-expressed genes. Retained: %d",
                         n_removed_genes, nrow(object)), "info")

        # --- Normalization ----------------------------------------------------
        if (verbose) .msg("Step 2: Normalization ...", "step")

        if (normMethod == "scran") {
            # Quick pre-clustering for pool-based size factors
            set.seed(p$seed)
            clusters <- scran::quickCluster(object, BPPARAM = BPPARAM)
            object   <- scran::computeSumFactors(object, clusters = clusters,
                                                 BPPARAM = BPPARAM)
            object   <- scater::logNormCounts(object)
        } else if (normMethod == "logNorm") {
            object <- scater::logNormCounts(object)
        } else if (normMethod == "scTransform") {
            # Pearson residuals stored as "pearson" assay
            if (requireNamespace("sctransform", quietly = TRUE)) {
                expr_mat <- as.matrix(
                    SummarizedExperiment::assay(object, "counts"))
                vst_res  <- sctransform::vst(expr_mat, verbosity = 0)
                SummarizedExperiment::assay(object, "pearson",
                                            withDimnames = TRUE) <-
                    vst_res$y
            } else {
                warning("sctransform not installed; falling back to logNorm.")
                object <- scater::logNormCounts(object)
            }
        } else {
            stop("normMethod must be 'scran', 'logNorm', or 'scTransform'.")
        }

        # --- HVG selection ---------------------------------------------------
        if (verbose) .msg("Step 3: Highly Variable Gene Selection ...", "step")
        norm_assay <- if (normMethod == "scTransform") "pearson" else "logcounts"
        gene_var   <- scran::modelGeneVar(object, assay.type = norm_assay,
                                          BPPARAM = BPPARAM)
        hvgs       <- scran::getTopHVGs(gene_var, n = nHVGs)
        if (verbose) .msg(sprintf("  Selected %d HVGs.", length(hvgs)), "info")

        # Store HVGs
        SummarizedExperiment::rowData(object)$isHVG <- rownames(object) %in% hvgs

        # --- PCA -------------------------------------------------------------
        if (verbose) .msg("Step 4: PCA ...", "step")
        object <- scater::runPCA(object, ncomponents = nPCs,
                                 subset_row = hvgs,
                                 exprs_values = norm_assay,
                                 BPPARAM = BPPARAM)

        # --- Batch correction ------------------------------------------------
        if (!is.null(batchVar) && batchVar != "none") {
            if (verbose)
                .msg(sprintf("Step 5: Batch Correction (%s) ...",
                             batchMethod), "step")

            cd <- SummarizedExperiment::colData(object)
            if (!batchVar %in% colnames(cd))
                stop("batchVar '", batchVar, "' not found in colData.")

            batch_labels <- cd[[batchVar]]

            if (batchMethod == "harmony") {
                .require_pkg("harmony")
                pca_mat <- SingleCellExperiment::reducedDim(object, "PCA")
                harmony_res <- harmony::HarmonyMatrix(
                    pca_mat, batch_labels,
                    do_pca = FALSE, verbose = FALSE)
                SingleCellExperiment::reducedDim(object, "HARMONY") <-
                    harmony_res

            } else if (batchMethod == "fastMNN") {
                .require_pkg("batchelor")
                lc_mat <- SingleCellExperiment::logcounts(object)[hvgs, ]
                mnn_out <- batchelor::fastMNN(lc_mat,
                                              batch = batch_labels,
                                              BPPARAM = BPPARAM)
                SingleCellExperiment::reducedDim(object, "MNN") <-
                    SingleCellExperiment::reducedDim(mnn_out, "corrected")
            }
        }

        # --- UMAP ------------------------------------------------------------
        if (verbose) .msg("Step 6: UMAP ...", "step")
        use_dimred <- if (!is.null(batchVar) && batchVar != "none") {
            if (batchMethod == "harmony") "HARMONY"
            else if (batchMethod == "fastMNN") "MNN"
            else "PCA"
        } else "PCA"

        set.seed(p$seed)
        object <- scater::runUMAP(object, dimred = use_dimred,
                                  BPPARAM = BPPARAM)

        # --- Clustering ------------------------------------------------------
        if (verbose) .msg("Step 7: Graph-based Clustering ...", "step")
        snn_graph <- scran::buildSNNGraph(object, use.dimred = use_dimred,
                                          BPPARAM = BPPARAM)
        clusters  <- igraph::cluster_leiden(snn_graph,
                                            resolution_parameter = p$resolution)
        SummarizedExperiment::colData(object)$cluster <-
            factor(igraph::membership(clusters))

        if (verbose)
            .msg(sprintf("  %d clusters identified.",
                         nlevels(
                             SummarizedExperiment::colData(object)$cluster)),
                 "success")

        # --- Update params ---------------------------------------------------
        object@params$minCells       <- minCells
        object@params$minGenes       <- minGenes
        object@params$maxMitoPercent <- maxMitoPercent
        object@params$normMethod     <- normMethod
        object@params$nHVGs          <- nHVGs
        object@params$nPCs           <- nPCs
        object@params$batchVar       <- batchVar
        object@params$batchMethod    <- batchMethod

        if (verbose) .msg("Preprocessing complete.", "success")
        object
    }
)
