## =============================================================================
## tests/testthat/test-CoMorTrax.R
## Comprehensive Unit Tests for the CoMorTrax Package
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

library(testthat)
library(CoMorTrax)
library(SingleCellExperiment)
library(S4Vectors)

## ---------------------------------------------------------------------------
## Helper: build minimal test object
## ---------------------------------------------------------------------------
.make_test_obj <- function(n_genes = 100, n_cells = 60,
                            n_disease = 15, seed = 42) {
    set.seed(seed)
    counts <- matrix(
        rnbinom(n_genes * n_cells, size = 10, mu = 5),
        nrow = n_genes, ncol = n_cells,
        dimnames = list(paste0("Gene", seq_len(n_genes)),
                        paste0("Cell", seq_len(n_cells)))
    )
    # Inject DE signal
    counts[1:10, 1:n_disease]                          <- counts[1:10, 1:n_disease] + 20
    counts[11:20, (n_disease+1):(2*n_disease)]         <- counts[11:20, (n_disease+1):(2*n_disease)] + 20
    counts[21:25, (2*n_disease+1):(3*n_disease)]       <- counts[21:25, (2*n_disease+1):(3*n_disease)] + 40
    counts[26:30, (2*n_disease+1):(3*n_disease)]       <- pmax(0, counts[26:30, (2*n_disease+1):(3*n_disease)] - 15)

    disease_vec <- c(
        rep("AD",      n_disease),
        rep("DM",      n_disease),
        rep("AD,DM",   n_disease),
        rep("control", n_cells - 3 * n_disease)
    )
    ct_vec    <- rep(c("TypeA","TypeB"), length.out = n_cells)
    batch_vec <- rep(c("batch1","batch2"), length.out = n_cells)

    cd  <- DataFrame(disease  = disease_vec,
                     cellType = ct_vec,
                     batch    = batch_vec)
    sce <- SingleCellExperiment(assays  = list(counts = counts),
                                colData = cd)
    sce <- scater::logNormCounts(sce)
    SummarizedExperiment::rowData(sce)$isHVG <- TRUE

    set.seed(seed)
    SingleCellExperiment::reducedDim(sce, "PCA")  <-
        matrix(rnorm(n_cells * 10), ncol = 10)
    SingleCellExperiment::reducedDim(sce, "UMAP") <-
        matrix(rnorm(n_cells * 2),  ncol = 2)

    createCoMorTraxObject(sce, projectName = "TestProject")
}

## ===========================================================================
## Test: createCoMorTraxObject
## ===========================================================================

test_that("createCoMorTraxObject works from SingleCellExperiment", {
    obj <- .make_test_obj()
    expect_s4_class(obj, "CoMorTraxObject")
    expect_equal(nrow(obj), 100)
    expect_equal(ncol(obj), 60)
    expect_true(length(obj@params) > 0)
    expect_equal(obj@params$projectName, "TestProject")
    expect_true(nchar(obj@coMorTraxVersion) > 0)
})

test_that("createCoMorTraxObject works from raw matrix", {
    set.seed(1)
    mat <- matrix(rpois(500, 5), nrow = 50, ncol = 10,
                  dimnames = list(paste0("G", seq_len(50)),
                                  paste0("C", seq_len(10))))
    obj <- createCoMorTraxObject(mat)
    expect_s4_class(obj, "CoMorTraxObject")
})

test_that("createCoMorTraxObject works from dgCMatrix", {
    set.seed(2)
    mat <- Matrix::sparseMatrix(
        i = sample(50, 100, replace = TRUE),
        j = sample(20, 100, replace = TRUE),
        x = rpois(100, 3),
        dims = c(50, 20),
        dimnames = list(paste0("G", seq_len(50)),
                        paste0("C", seq_len(20)))
    )
    obj <- createCoMorTraxObject(mat)
    expect_s4_class(obj, "CoMorTraxObject")
})

test_that("params can be overridden", {
    obj <- .make_test_obj()
    obj2 <- createCoMorTraxObject(
        as(obj, "SingleCellExperiment"),
        params = list(nHVGs = 500L))
    expect_equal(obj2@params$nHVGs, 500L)
})

## ===========================================================================
## Test: S4 validity
## ===========================================================================

test_that("CoMorTraxObject validity rejects bad diseaseMatrix dims", {
    obj <- .make_test_obj()
    bad_dm <- DataFrame(AD = rep(0, nrow(obj)))  # nrow != ncol
    expect_error(
        {obj@diseaseMatrix <- bad_dm; validObject(obj)},
        regexp = "rows"
    )
})

test_that("CoMorTraxObject validity rejects non-binary diseaseMatrix", {
    obj <- .make_test_obj()
    dm  <- DataFrame(AD = rep(2L, ncol(obj)))  # invalid: 2
    expect_error(
        {obj@diseaseMatrix <- dm; validObject(obj)},
        regexp = "binary"
    )
})

## ===========================================================================
## Test: encodeDiseaseLabels
## ===========================================================================

test_that("encodeDiseaseLabels populates diseaseMatrix correctly", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control",
                               multiLabelSep = ",",
                               verbose = FALSE)
    dm <- diseaseMatrix(obj)
    expect_true(nrow(dm) == ncol(obj))
    expect_true(all(colnames(dm) %in% c("AD", "DM")))
    # All values binary
    expect_true(all(unlist(as.list(dm)) %in% c(0L, 1L)))
})

test_that("encodeDiseaseLabels creates correct disease groups", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control",
                               verbose = FALSE)
    grps <- diseaseGroups(obj)
    expect_true("control" %in% names(grps))
    expect_true(any(grepl("AD", names(grps))))
    expect_true(any(grepl("DM", names(grps))))
    # Comorbid group exists
    expect_true(any(vapply(names(grps), function(g)
        all(c("AD","DM") %in% strsplit(g,"_")[[1]]), logical(1))))
})

test_that("colData diseaseGroup column is populated", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control",
                               verbose = FALSE)
    cd <- SummarizedExperiment::colData(obj)
    expect_true("diseaseGroup" %in% colnames(cd))
    expect_equal(length(cd$diseaseGroup), ncol(obj))
})

test_that("listComorbidityComparisons returns a data.frame", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control",
                               verbose = FALSE)
    comps <- listComorbidityComparisons(obj)
    expect_s3_class(comps, "data.frame")
    expect_true(nrow(comps) > 0)
    expect_true(all(c("comparison","numerator","denominator") %in%
                        colnames(comps)))
})

## ===========================================================================
## Test: runComorbidityDE
## ===========================================================================

test_that("runComorbidityDE produces results for all comparisons", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon",
                            fdrThreshold = 0.05, verbose = FALSE)
    de <- deResults(obj)
    expect_true(length(de) > 0)
    # Each result is a data.frame
    for (res in de) {
        expect_s3_class(res, "data.frame")
        expect_true(all(c("gene","logFC","pvalue","padj",
                          "significant","direction") %in% colnames(res)))
    }
})

test_that("runComorbidityDE stratifies by cell type", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, cellTypeCol = "cellType",
                            method = "wilcoxon", verbose = FALSE)
    de_names <- names(deResults(obj))
    expect_true(any(grepl("TypeA", de_names)))
    expect_true(any(grepl("TypeB", de_names)))
})

test_that("deResults accessor with comparison argument works", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    first_comp <- names(deResults(obj))[1]
    single_res <- deResults(obj, comparison = first_comp)
    expect_s3_class(single_res, "data.frame")
})

test_that("deResults accessor errors on unknown comparison", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    expect_error(deResults(obj, comparison = "NONEXISTENT"),
                 regexp = "not found")
})

test_that("getDeTable returns combined data.frame", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    all_de <- getDeTable(obj)
    expect_s3_class(all_de, "data.frame")
    expect_true(nrow(all_de) > 0)
})

## ===========================================================================
## Test: quantifyInteractionEffects
## ===========================================================================

test_that("quantifyInteractionEffects returns valid DataFrame", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj,
                                      nPermutations = 10L,
                                      verbose = FALSE)
    is_df <- interactionScores(obj)
    expect_true(nrow(is_df) > 0)
    expect_true(all(c("gene","observed","expected","interactionDelta",
                      "classification") %in% colnames(is_df)))
})

test_that("interaction classifications are restricted to valid values", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)
    cls_vals <- unique(interactionScores(obj)$classification)
    expect_true(all(cls_vals %in% c("synergistic","antagonistic","additive")))
})

test_that("getSynergisticGenes / getAntagonisticGenes return character vectors", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)
    syn <- getSynergisticGenes(obj)
    ant <- getAntagonisticGenes(obj)
    expect_type(syn, "character")
    expect_type(ant, "character")
    # No gene should be both synergistic and antagonistic
    expect_equal(length(intersect(syn, ant)), 0L)
})

test_that("interaction scores have p-values in [0,1]", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)
    pvals <- interactionScores(obj)$pvalue
    expect_true(all(pvals >= 0 & pvals <= 1, na.rm = TRUE))
})

## ===========================================================================
## Test: scoreCellVulnerability
## ===========================================================================

test_that("scoreCellVulnerability returns ranked DataFrame", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)
    obj <- scoreCellVulnerability(obj, cellTypeCol = "cellType",
                                  includePathwayScore = FALSE,
                                  verbose = FALSE)
    cv <- cellVulnerability(obj)
    expect_true(nrow(cv) > 0)
    req_cols <- c("cellType","nSynergistic","nAntagonistic",
                  "vulnerabilityScore","vulnerabilityRank")
    expect_true(all(req_cols %in% colnames(cv)))
    # Ranks are sequential
    expect_equal(sort(cv$vulnerabilityRank), seq_len(nrow(cv)))
})

test_that("vulnerability scores are in [0,1]", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)
    obj <- scoreCellVulnerability(obj, includePathwayScore = FALSE,
                                  verbose = FALSE)
    scores <- cellVulnerability(obj)$vulnerabilityScore
    expect_true(all(scores >= 0 & scores <= 1, na.rm = TRUE))
})

## ===========================================================================
## Test: buildComorbidityNetworks
## ===========================================================================

test_that("buildComorbidityNetworks creates igraph objects", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- buildComorbidityNetworks(obj, method = "correlation",
                                    maxGenes = 50L,
                                    computeRewiring = TRUE,
                                    verbose = FALSE)
    nets <- networks(obj)
    expect_true(length(nets) > 0)
    for (g in nets) expect_s3_class(g, "igraph")
})

test_that("network vertex count equals gene subset size", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- buildComorbidityNetworks(obj, maxGenes = 30L,
                                    computeRewiring = FALSE,
                                    verbose = FALSE)
    nets <- networks(obj)
    for (g in nets) {
        expect_lte(igraph::vcount(g), 30L)
    }
})

test_that("getNetworkHubs returns named numeric vector", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- buildComorbidityNetworks(obj, maxGenes = 30L,
                                    computeRewiring = FALSE,
                                    verbose = FALSE)
    state <- names(networks(obj))[1]
    hubs  <- getNetworkHubs(obj, state = state, topN = 5L)
    expect_type(hubs, "double")
    expect_lte(length(hubs), 5L)
    expect_true(!is.null(names(hubs)))
})

## ===========================================================================
## Test: trainComorbidityClassifier
## ===========================================================================

test_that("trainComorbidityClassifier runs without error (RF)", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)
    obj <- trainComorbidityClassifier(obj,
                                      method      = "randomForest",
                                      featureType = "interaction",
                                      cvFolds     = 3L,
                                      verbose     = FALSE)
    clf <- classifier(obj)
    expect_true(length(clf) > 0)
    expect_true(!is.null(clf$model))
    expect_true(!is.null(clf$performance))
    expect_true(!is.null(clf$predictions))
})

test_that("classifier performance auc is in [0,1] or NA", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)
    obj <- trainComorbidityClassifier(obj, method = "randomForest",
                                      cvFolds = 3L, verbose = FALSE)
    auc <- classifier(obj)$performance$auc
    expect_true(is.na(auc) || (auc >= 0 && auc <= 1))
})

test_that("classifier predictions have correct number of rows", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)
    obj <- trainComorbidityClassifier(obj, method = "randomForest",
                                      cvFolds = 3L, verbose = FALSE)
    preds <- classifier(obj)$predictions
    expect_equal(nrow(preds), ncol(obj))
})

## ===========================================================================
## Test: show and summaryCoMorTrax
## ===========================================================================

test_that("show method runs without error", {
    obj <- .make_test_obj()
    expect_output(show(obj), "CoMorTraxObject")
})

test_that("summaryCoMorTrax runs without error on empty object", {
    obj <- .make_test_obj()
    expect_output(summaryCoMorTrax(obj), "CoMorTrax Analysis Summary")
})

test_that("summaryCoMorTrax works after full pipeline", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)
    expect_output(summaryCoMorTrax(obj), "Interaction effects")
})

## ===========================================================================
## Test: simulateCoMorTrax
## ===========================================================================

test_that("simulateCoMorTrax returns a valid SingleCellExperiment", {
    # new snake_case API
    sce <- simulateCoMorTrax(n_cells = 80, n_genes = 100,
                              n_syn = 5, n_ant = 4, seed = 1)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_equal(nrow(sce), 100)
    # colData has disease column
    cd <- SummarizedExperiment::colData(sce)
    expect_true("disease" %in% colnames(cd))
    # logcounts present
    expect_true("logcounts" %in% SummarizedExperiment::assayNames(sce))
    # sim_params in metadata
    sp <- S4Vectors::metadata(sce)$sim_params
    expect_equal(length(sp$syn_genes), 5)
    expect_equal(length(sp$ant_genes), 4)

    # legacy camelCase aliases still work
    sce2 <- simulateCoMorTrax(nGenes = 80, nCells = 60,
                               nDiseaseCells = 15, seed = 2)
    expect_s4_class(sce2, "SingleCellExperiment")
    expect_equal(nrow(sce2), 80)
})

## ===========================================================================
## Test: exportCoMorTrax
## ===========================================================================

test_that("exportCoMorTrax writes files to disk", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    obj <- runComorbidityDE(obj, method = "wilcoxon", verbose = FALSE)
    obj <- quantifyInteractionEffects(obj, nPermutations = 10L, verbose = FALSE)

    out_dir <- tempfile("CoMorTrax_test_export_")
    exportCoMorTrax(obj, outDir = out_dir, formats = c("csv","rds"),
                    verbose = FALSE)

    expect_true(file.exists(file.path(out_dir, "CoMorTraxObject.rds")))
    expect_true(file.exists(file.path(out_dir, "disease_matrix.csv")))
    expect_true(file.exists(file.path(out_dir, "interaction_scores.csv")))
    unlink(out_dir, recursive = TRUE)
})

## ===========================================================================
## Test: Accessor replacement methods
## ===========================================================================

test_that("diseaseMatrix<- validates and stores correctly", {
    obj <- .make_test_obj()
    dm  <- DataFrame(AD = rep(0L, ncol(obj)),
                     DM = rep(1L, ncol(obj)))
    diseaseMatrix(obj) <- dm
    expect_equal(nrow(diseaseMatrix(obj)), ncol(obj))
})

test_that("interactionScores<- validates required columns", {
    obj  <- .make_test_obj()
    bad  <- DataFrame(gene = "G1", observed = 1.0)  # missing columns
    expect_error(
        {interactionScores(obj) <- bad; validObject(obj)},
        regexp = "missing"
    )
})

## ===========================================================================
## Test: edge cases
## ===========================================================================

test_that("runComorbidityDE skips groups with too few cells", {
    obj <- .make_test_obj(n_cells = 40, n_disease = 5)
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    # minCellsDE = 10 > 5, so many comparisons should be skipped
    expect_warning(
        runComorbidityDE(obj, method = "wilcoxon",
                         minCellsDE = 10L, verbose = FALSE),
        NA  # no error — just messages
    )
})

test_that("quantifyInteractionEffects errors without prior DE", {
    obj <- .make_test_obj()
    obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
                               controlLabel = "control", verbose = FALSE)
    expect_error(quantifyInteractionEffects(obj, verbose = FALSE),
                 regexp = "runComorbidityDE")
})

test_that("scoreCellVulnerability errors without interaction scores", {
    obj <- .make_test_obj()
    expect_error(scoreCellVulnerability(obj, verbose = FALSE),
                 regexp = "quantifyInteractionEffects")
})
