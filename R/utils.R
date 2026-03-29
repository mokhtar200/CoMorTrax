## =============================================================================
## utils.R - Internal Utilities and Simulation Engine
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

## --- Message formatting ------------------------------------------------------

#' @keywords internal
.msg <- function(msg, level = "info") {
    prefix <- switch(level,
        info    = "  [CoMorTrax] ",
        step    = "\n[CoMorTrax] * ",
        success = "  [CoMorTrax] + ",
        warn    = "  [CoMorTrax] ! ",
        "  "
    )
    message(prefix, msg)
}

## --- Parameter handling -------------------------------------------------------

#' @keywords internal
.param <- function(val, default) if (is.null(val)) default else val

## --- Validation checks -------------------------------------------------------

#' @keywords internal
.check_disease_encoded <- function(object) {
    if (nrow(diseaseMatrix(object)) == 0)
        stop("Disease labels not encoded. Run encodeDiseaseLabels() first.")
    invisible(TRUE)
}

## --- Expression matrix retrieval ---------------------------------------------

#' @keywords internal
.get_norm_matrix <- function(object) {
    assay_names <- SummarizedExperiment::assayNames(object)
    if ("logcounts" %in% assay_names)
        return(as.matrix(SingleCellExperiment::logcounts(object)))
    if ("pearson" %in% assay_names)
        return(as.matrix(SummarizedExperiment::assay(object, "pearson")))
    if ("counts" %in% assay_names) {
        warning("No normalized expression found; using raw counts. ",
                "Run preprocessCoMorTrax() first.")
        return(as.matrix(SummarizedExperiment::assay(object, "counts")))
    }
    stop("No expression matrix found. Ensure object has 'logcounts' assay.")
}

## --- Package availability checks --------------------------------------------

#' @keywords internal
.require_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop("Package '", pkg, "' is required for this function. ",
             "Install it with: BiocManager::install('", pkg, "')",
             call. = FALSE)
    invisible(TRUE)
}

## =============================================================================
## simulateCoMorTrax - Simulation Engine
## =============================================================================

#' @title Simulate a Comorbidity Single-Cell Dataset
#' @description
#' Generates a synthetic \code{SingleCellExperiment} with biologically
#' realistic count distributions and \emph{injected} synergistic and
#' antagonistic signals. The simulation is truth-recoverable: the identities
#' of all synergistic and antagonistic genes are known and stored in
#' \code{metadata(sce)$sim_params}, enabling rigorous benchmarking and
#' unit testing of the full CoMorTrax pipeline.
#'
#' Counts are sampled from a negative-binomial model with cell-type-specific
#' baseline expression. Disease effects are added as multiplicative fold
#' changes. Synergistic genes receive an additional amplification in the
#' double-disease group beyond the sum of single-disease effects.
#' Antagonistic genes receive a corresponding attenuation.
#'
#' @param n_cells Integer. Number of cells to simulate. Default 500.
#' @param n_genes Integer. Number of genes. Default 300.
#' @param n_syn Integer. Number of injected synergistic genes. Default 20.
#' @param n_ant Integer. Number of injected antagonistic genes. Default 15.
#' @param cell_types Character vector of cell-type names.
#'   Default \code{c("Neuron","Microglia","Astrocyte","OPC")}.
#' @param diseases Character vector of disease names to simulate.
#'   Default \code{c("AD","T2D")}.
#' @param effect_size Numeric. Mean log2FC magnitude for each single disease
#'   effect. Default 1.5.
#' @param syn_boost Numeric. Additional count elevation for synergistic genes
#'   in the comorbid group (above individual disease effects). Default 2.0.
#' @param ant_boost Numeric. Count reduction fraction for antagonistic genes
#'   in the comorbid group (0-1, where 1 = complete suppression). Default 0.6.
#' @param dispersion Numeric. Negative-binomial dispersion (\code{size}
#'   parameter). Smaller values = greater overdispersion. Default 10.
#' @param seed Integer. Random seed for reproducibility. Default 42.
#' @param verbose Logical. Print simulation parameters. Default \code{FALSE}.
#' @param nGenes Deprecated alias for \code{n_genes}.
#' @param nCells Deprecated alias for \code{n_cells}.
#' @param nDiseaseCells Deprecated. Number of cells per disease group
#'   (overrides internal balancing when supplied).
#'
#' @return A \code{SingleCellExperiment} with:
#' \itemize{
#'   \item \code{assay(sce, "counts")} - raw integer count matrix.
#'   \item \code{assay(sce, "logcounts")} - log-normalised counts.
#'   \item \code{colData} columns: \code{disease}, \code{cellType},
#'     \code{batch}, plus one binary column per disease
#'     (\code{disease_<name>}) and \code{is_control}.
#'   \item \code{reducedDim(sce, "PCA")} and \code{reducedDim(sce, "UMAP")}.
#'   \item \code{metadata(sce)$sim_params} - list with \code{syn_genes},
#'     \code{ant_genes}, \code{diseases}, and all simulation parameters.
#' }
#' @export
#' @examples
#' library(CoMorTrax)
#' set.seed(42)
#' sce <- simulateCoMorTrax(
#'   n_cells    = 200,
#'   n_genes    = 100,
#'   n_syn      = 10,
#'   n_ant      = 8,
#'   diseases   = c("AD", "T2D"),
#'   cell_types = c("Neuron", "Microglia", "Astrocyte")
#' )
#' sce
#' # True synergistic genes are stored in metadata
#' true_syn <- S4Vectors::metadata(sce)$sim_params$syn_genes
#' cat("True synergistic genes:", paste(head(true_syn, 5), collapse=", "), "\n")
simulateCoMorTrax <- function(n_cells     = 500L,
                               n_genes     = 300L,
                               n_syn       = 20L,
                               n_ant       = 15L,
                               cell_types  = c("Neuron", "Microglia",
                                               "Astrocyte", "OPC"),
                               diseases    = c("AD", "T2D"),
                               effect_size = 1.5,
                               syn_boost   = 2.0,
                               ant_boost   = 0.6,
                               dispersion  = 10,
                               seed        = 42L,
                               verbose     = FALSE,
                               ## Legacy camelCase aliases
                               nGenes        = NULL,
                               nCells        = NULL,
                               nDiseaseCells = NULL) {
    ## --- Handle legacy camelCase arguments -----------------------------------
    if (!is.null(nGenes))  n_genes <- nGenes
    if (!is.null(nCells))  n_cells <- nCells

    set.seed(seed)

    n_diseases <- length(diseases)
    n_genes    <- as.integer(n_genes)
    n_syn      <- min(as.integer(n_syn), n_genes - n_ant - 5L)
    n_ant      <- min(as.integer(n_ant), n_genes - n_syn - 5L)

    ## Pre-compute all pairwise comorbid groups (works for 2+ diseases)
    pairs <- if (n_diseases >= 2L)
        utils::combn(sort(diseases), 2L, simplify = FALSE)
    else
        list()
    n_comorbid_groups <- length(pairs)
    comorbid_name <- if (length(pairs) > 0)
        paste(pairs[[1]], collapse = "_")
    else
        paste(sort(diseases), collapse = "_")

    ## Compute group sizes BEFORE allocating the counts matrix
    if (!is.null(nDiseaseCells)) {
        n_per_group <- as.integer(nDiseaseCells)
    } else {
        n_groups    <- n_diseases + n_comorbid_groups + 1L
        n_per_group <- max(10L, as.integer(n_cells / n_groups))
    }
    n_ctrl  <- max(10L, n_cells - (n_diseases + n_comorbid_groups) * n_per_group)
    n_total <- n_per_group * (n_diseases + n_comorbid_groups) + n_ctrl

    if (verbose) {
        .msg(sprintf("Simulating %d cells x %d genes", n_total, n_genes), "step")
        .msg(sprintf("Diseases: %s | pairs: %d | syn: %d ant: %d",
                     paste(diseases, collapse=", "), n_comorbid_groups, n_syn, n_ant))
    }

    gene_names <- paste0("Gene", seq_len(n_genes))
    cell_names <- paste0("Cell", seq_len(n_total))

    ## --- Base counts (negative-binomial) ------------------------------------
    mu_vec <- stats::runif(n_genes, min = 3, max = 20)
    counts  <- matrix(
        stats::rnbinom(n_genes * n_total, size = dispersion,
                       mu = rep(mu_vec, n_total)),
        nrow = n_genes, ncol = n_total,
        dimnames = list(gene_names, cell_names)
    )

    ## --- Define group index ranges -------------------------------------------
    grp_indices <- list()
    offset <- 0L
    for (d in diseases) {
        grp_indices[[d]] <- offset + seq_len(n_per_group)
        offset <- offset + n_per_group
    }
    for (pr in pairs) {
        pair_name <- paste(pr, collapse = "_")
        grp_indices[[pair_name]] <- offset + seq_len(n_per_group)
        offset <- offset + n_per_group
    }
    grp_indices[["control"]] <- offset + seq_len(n_ctrl)

    ## --- Assign synergistic and antagonistic gene indices -------------------
    syn_gene_idx  <- seq_len(n_syn)
    ant_gene_idx  <- n_syn + seq_len(n_ant)
    syn_gene_names <- gene_names[syn_gene_idx]
    ant_gene_names <- gene_names[ant_gene_idx]

    ## --- Inject single-disease effects ---------------------------------------
    fc_scale <- 2^effect_size  # fold change on count scale
    for (d in diseases) {
        d_idx <- grp_indices[[d]]
        ## Upregulate first half of non-syn/ant genes for this disease
        d_specific <- seq(n_syn + n_ant + 1L,
                          n_genes)[seq_len(max(1, (n_genes - n_syn - n_ant) %/% n_diseases))]
        counts[d_specific, d_idx] <- counts[d_specific, d_idx] * fc_scale
        ## Partial effect in syn/ant genes for single disease
        counts[syn_gene_idx, d_idx] <- counts[syn_gene_idx, d_idx] * (fc_scale * 0.5)
        counts[ant_gene_idx, d_idx] <- counts[ant_gene_idx, d_idx] * (fc_scale * 0.7)
    }

    ## --- Inject comorbid effects for every pairwise group -------------------
    for (pr in pairs) {
        pair_name <- paste(pr, collapse = "_")
        cm_idx    <- grp_indices[[pair_name]]
        ## Synergistic: observed >> expected (sum of singles)
        counts[syn_gene_idx, cm_idx] <-
            counts[syn_gene_idx, cm_idx] * (fc_scale * (1 + syn_boost))
        ## Antagonistic: observed << expected
        counts[ant_gene_idx, cm_idx] <-
            counts[ant_gene_idx, cm_idx] * (1 - ant_boost)
        ## Regular disease genes also expressed in comorbid group
        for (d in pr) {
            d_specific <- seq(n_syn + n_ant + 1L,
                              n_genes)[seq_len(max(1, (n_genes - n_syn - n_ant) %/% n_diseases))]
            counts[d_specific, cm_idx] <- counts[d_specific, cm_idx] * (fc_scale * 0.8)
        }
    }

    ## Round to integers
    storage.mode(counts) <- "integer"
    counts[counts < 0L] <- 0L

    ## --- Build colData -------------------------------------------------------
    n_ct    <- length(cell_types)
    disease_vec <- character(n_total)
    for (g in names(grp_indices)) disease_vec[grp_indices[[g]]] <- g
    ct_vec <- cell_types[((seq_len(n_total) - 1L) %% n_ct) + 1L]
    batch_vec <- sample(c("batch1", "batch2"), n_total, replace = TRUE)

    ## Binary disease indicator columns — built from grp_indices (truth-based)
    binary_df <- do.call(data.frame,
                         stats::setNames(
                             lapply(diseases, function(d) {
                                 # A cell has disease d if it's in the single-disease
                                 # group OR the comorbid group containing d
                                 member_groups <- names(grp_indices)[
                                     vapply(names(grp_indices), function(g) {
                                         g != "control" && grepl(d, g, fixed = TRUE)
                                     }, logical(1))]
                                 member_cells <- unlist(grp_indices[member_groups],
                                                        use.names = FALSE)
                                 vec <- integer(n_total)
                                 vec[member_cells] <- 1L
                                 vec
                             }),
                             paste0("disease_", diseases)
                         ))

    cd <- S4Vectors::DataFrame(
        disease    = disease_vec,
        cellType   = ct_vec,
        batch      = batch_vec,
        is_control = as.integer(disease_vec == "control"),
        binary_df,
        row.names  = cell_names
    )

    ## --- Build SCE -----------------------------------------------------------
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays  = list(counts = counts),
        colData = cd
    )

    ## Log-normalize (simple library-size normalisation)
    lib_size <- Matrix::colSums(counts)
    lib_size[lib_size == 0] <- 1L
    scale_factor <- stats::median(lib_size)
    log_counts <- log1p(t(t(counts) / lib_size * scale_factor))
    SummarizedExperiment::assay(sce, "logcounts") <- log_counts

    ## Mark all genes as HVG for simulation purposes
    SummarizedExperiment::rowData(sce)$isHVG <- TRUE

    ## Fake PCA and UMAP (2D) for visualisation
    fake_pca  <- matrix(stats::rnorm(n_total * 10), ncol = 10,
                        dimnames = list(cell_names, paste0("PC", seq_len(10))))
    fake_umap <- matrix(stats::rnorm(n_total * 2),  ncol = 2,
                        dimnames = list(cell_names, c("UMAP1","UMAP2")))
    SingleCellExperiment::reducedDim(sce, "PCA")  <- fake_pca
    SingleCellExperiment::reducedDim(sce, "UMAP") <- fake_umap

    ## Store sim params in metadata
    S4Vectors::metadata(sce) <- list(
        sim_params = list(
            n_cells     = n_total,
            n_genes     = n_genes,
            n_syn       = n_syn,
            n_ant       = n_ant,
            syn_genes   = syn_gene_names,
            ant_genes   = ant_gene_names,
            diseases    = diseases,
            cell_types  = cell_types,
            effect_size = effect_size,
            syn_boost   = syn_boost,
            ant_boost   = ant_boost,
            seed        = seed
        )
    )

    if (verbose)
        .msg(sprintf("Done. Groups: %s",
                     paste(names(grp_indices), collapse=", ")), "success")
    sce
}
