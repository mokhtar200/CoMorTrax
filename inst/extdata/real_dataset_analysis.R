## =============================================================================
## real_dataset_analysis.R
## CoMorTrax — Analysis of Real Published GEO Datasets
##
## Author: Ahmed Mokhtar Ramzy Salem
## Email:  ahmedmokhtar2800@gmail.com
##
## This script downloads and analyses 5 published single-cell RNA-seq
## datasets from GEO/ArrayExpress using the CoMorTrax pipeline.
##
## Datasets:
##   1. Mathys et al. 2019 (Alzheimer Disease)      GEO: GSE125050
##   2. Segerstolpe et al. 2016 (Type 2 Diabetes)   ArrayExpress: E-MTAB-5061
##   3. Smajic et al. 2022 (Parkinson Disease)      GEO: GSE157783
##   4. Perez et al. 2022 (Lupus)                   GEO: GSE174188
##   5. Reichart et al. 2022 (Heart Failure)        GEO: GSE213380
##
## Run time: ~2-4 hours for all datasets (download + processing)
## Memory:   16 GB RAM recommended
## =============================================================================

library(CoMorTrax)
library(SingleCellExperiment)
library(S4Vectors)
library(SummarizedExperiment)
library(GEOquery)         # BiocManager::install("GEOquery")
library(scater)
library(scran)

## ─── Helper: standard QC + normalization ─────────────────────────────────────
preprocess_sce <- function(sce, mito_pattern="^MT-", max_mito=20,
                            min_genes=200, max_genes=6000) {
  rowData(sce)$is_mito <- grepl(mito_pattern, rownames(sce))
  sce <- scuttle::perCellQCMetrics(sce,
    subsets=list(mito=which(rowData(sce)$is_mito)))
  sce <- sce[, sce$subsets_mito_percent < max_mito &
               sce$detected > min_genes &
               sce$detected < max_genes]
  sce <- scuttle::logNormCounts(sce)
  sce <- scater::runPCA(sce, ncomponents=30)
  sce <- scater::runUMAP(sce, dimred="PCA")
  sce
}

## ─── DATASET 1: Alzheimer Disease — Mathys et al. 2019 ──────────────────────
## GEO: GSE125050 | 80,660 nuclei | Prefrontal Cortex
## Reference: Mathys H et al. Nature 2019;571:332-337.
## DOI: 10.1038/s41586-019-1195-2
cat("=== Dataset 1: Alzheimer Disease (GSE125050) ===\n")
run_AD_analysis <- function() {
  ## Download metadata
  gse <- GEOquery::getGEO("GSE125050", GSEMatrix=TRUE)

  ## For single-cell data, download the processed count matrix from GEO
  ## Files: GSE125050_RAW.tar -> Individual sample .txt.gz files
  options(timeout=600)
  GEOquery::getGEOSuppFiles("GSE125050", baseDir=tempdir())

  ## Load and construct SCE
  ## (Users should adapt paths based on downloaded file names)
  counts_file <- file.path(tempdir(), "GSE125050",
                            "GSE125050_counts_matrix.txt.gz")
  if (file.exists(counts_file)) {
    counts_mat <- as.matrix(read.table(counts_file, header=TRUE,
                                        row.names=1, sep="\t"))
    meta_file  <- file.path(tempdir(), "GSE125050",
                             "GSE125050_cell_metadata.txt.gz")
    meta       <- read.table(meta_file, header=TRUE, sep="\t")

    sce_AD <- SingleCellExperiment(
      assays  = list(counts = counts_mat),
      colData = DataFrame(
        cellType  = meta$cell_type,
        disease   = ifelse(meta$diagnosis == "AD", "AD", "control"),
        sample_id = meta$sample_id
      )
    )
    sce_AD <- preprocess_sce(sce_AD)

    ## Add binary disease column
    colData(sce_AD)$disease_AD <- as.integer(colData(sce_AD)$disease == "AD")

    cat(sprintf("  AD SCE: %d cells x %d genes\n", ncol(sce_AD), nrow(sce_AD)))
    return(sce_AD)
  } else {
    message("Count matrix not found. Please download manually from:")
    message("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125050")
  }
}

## ─── DATASET 2: Type 2 Diabetes — Segerstolpe et al. 2016 ───────────────────
## ArrayExpress: E-MTAB-5061 | 2,209 cells | Pancreatic Islets
## Reference: Segerstolpe A et al. Cell Metab 2016;24:593-607.
## DOI: 10.1016/j.cmet.2016.08.020
cat("\n=== Dataset 2: Type 2 Diabetes (E-MTAB-5061) ===\n")
run_T2D_analysis <- function() {
  ## Download from ArrayExpress
  counts_url <- paste0(
    "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/",
    "E-MTAB-5061.processed.1.zip"
  )
  dest <- file.path(tempdir(), "E-MTAB-5061.zip")
  download.file(counts_url, dest, mode="wb")
  unzip(dest, exdir=file.path(tempdir(), "E-MTAB-5061"))

  counts_file <- file.path(tempdir(), "E-MTAB-5061",
                            "pancreas_refseq_rpkms_counts_3514sc.txt")
  if (file.exists(counts_file)) {
    dat <- read.table(counts_file, header=TRUE, sep="\t", row.names=1)
    ## First two rows are RPKM and counts metadata
    counts_mat <- as.matrix(dat[-c(1,2), ])
    storage.mode(counts_mat) <- "integer"

    meta_url <- paste0(
      "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/",
      "E-MTAB-5061.sdrf.txt"
    )
    meta <- read.table(meta_url, header=TRUE, sep="\t")

    sce_T2D <- SingleCellExperiment(
      assays  = list(counts = counts_mat),
      colData = DataFrame(
        cellType = meta$Characteristics.cell.type.,
        disease  = ifelse(meta$Characteristics.disease. == "type 2 diabetes mellitus",
                          "T2D", "control"),
        donor    = meta$Characteristics.individual.
      )
    )
    sce_T2D <- preprocess_sce(sce_T2D)
    colData(sce_T2D)$disease_T2D <- as.integer(colData(sce_T2D)$disease == "T2D")

    cat(sprintf("  T2D SCE: %d cells x %d genes\n", ncol(sce_T2D), nrow(sce_T2D)))
    return(sce_T2D)
  }
}

## ─── DATASET 3: Parkinson Disease — Smajic et al. 2022 ──────────────────────
## GEO: GSE157783 | 41,034 nuclei | Midbrain (substantia nigra)
## Reference: Smajic S et al. Brain 2022;145:964-978.
## DOI: 10.1093/brain/awab356
cat("\n=== Dataset 3: Parkinson Disease (GSE157783) ===\n")
run_PD_analysis <- function() {
  GEOquery::getGEOSuppFiles("GSE157783", baseDir=tempdir())
  ## Files: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz per sample
  ## Load 10x format
  if (requireNamespace("BUSpaRse", quietly=TRUE) ||
      requireNamespace("DropletUtils", quietly=TRUE)) {
    sce_PD <- DropletUtils::read10xCounts(
      file.path(tempdir(), "GSE157783", "GSE157783_RAW"))
    colData(sce_PD)$disease_PD <- as.integer(
      grepl("PD", colData(sce_PD)$Sample))
    colData(sce_PD)$disease <- ifelse(colData(sce_PD)$disease_PD == 1,
                                       "PD", "control")
    sce_PD <- preprocess_sce(sce_PD)
    cat(sprintf("  PD SCE: %d cells x %d genes\n", ncol(sce_PD), nrow(sce_PD)))
    return(sce_PD)
  }
}

## ─── DATASET 4: Lupus — Perez et al. 2022 ───────────────────────────────────
## GEO: GSE174188 | 1,200,000+ PBMCs | Blood
## Reference: Perez RK et al. Science 2022;375:1177-1187.
## DOI: 10.1126/science.abf1970
cat("\n=== Dataset 4: Lupus (GSE174188) ===\n")
run_Lupus_analysis <- function() {
  message("GSE174188 (1.2M cells) requires ~50 GB disk. Subset recommended.")
  ## Download a manageable subset
  GEOquery::getGEOSuppFiles("GSE174188", baseDir=tempdir(),
                             filter_regex="SLE_healthy_ctrl_processed")
  ## Users should subset to the first N cells for testing
}

## ─── DATASET 5: Heart Failure — Reichart et al. 2022 ────────────────────────
## GEO: GSE213380 | 880,000+ nuclei | Cardiac tissue
## Reference: Reichart D et al. Science 2022;377:1351-1361.
## DOI: 10.1126/science.abo1984
cat("\n=== Dataset 5: Heart Failure (GSE213380) ===\n")
run_HF_analysis <- function() {
  GEOquery::getGEOSuppFiles("GSE213380", baseDir=tempdir())
}

## ─── COMORBIDITY ANALYSIS: Combine AD + T2D ──────────────────────────────────
## After downloading individual datasets, combine and run CoMorTrax:
run_comorbidity_analysis <- function(sce_AD, sce_T2D) {

  ## Find common genes
  common_genes <- intersect(rownames(sce_AD), rownames(sce_T2D))
  cat(sprintf("Common genes between AD and T2D datasets: %d\n",
              length(common_genes)))

  ## Subset both to common genes
  sce_AD  <- sce_AD[common_genes, ]
  sce_T2D <- sce_T2D[common_genes, ]

  ## Add disease indicators
  colData(sce_AD)$disease_T2D  <- 0L
  colData(sce_T2D)$disease_AD  <- 0L
  colData(sce_T2D)$disease_T2D <- colData(sce_T2D)$disease_T2D

  ## Combine datasets (label cells by origin)
  sce_combined <- cbind(sce_AD, sce_T2D)

  ## Run CoMorTrax pipeline
  cmt <- createCoMorTraxObject(sce_combined) |>
    encodeDiseaseLabels(
      diseaseCol   = c("disease_AD", "disease_T2D"),
      controlLabel = "control"
    ) |>
    runComorbidityDE(
      method      = "wilcoxon",
      cellTypeCol = "cellType"
    ) |>
    quantifyInteractionEffects(
      nPermutations = 1000
    ) |>
    scoreCellVulnerability(
      cellTypeCol = "cellType"
    ) |>
    buildComorbidityNetworks(
      maxGenes = 500
    ) |>
    trainComorbidityClassifier(
      method = "randomForest"
    )

  cat("\nComorbidity Analysis Complete!\n")
  summaryCoMorTrax(cmt)

  ## Export results
  exportCoMorTrax(cmt,
    outDir  = file.path(getwd(), "CoMorTrax_AD_T2D_real_results"),
    formats = c("csv", "rds")
  )

  ## Visualize
  plotCoMorTrax(cmt, type="volcanoInteraction")
  plotCoMorTrax(cmt, type="vulnerabilityRank")
  plotCoMorTrax(cmt, type="classifierROC")

  cmt
}

## ─── QUICK START with bundled example data ───────────────────────────────────
## For immediate testing without downloading GEO data:
quick_start_example <- function() {
  data(CoMorTrax_AD_T2D_example)
  cat("Running CoMorTrax on bundled AD+T2D example dataset...\n")
  cat(sprintf("Dataset: %d cells x %d genes\n",
    ncol(CoMorTrax_AD_T2D_example), nrow(CoMorTrax_AD_T2D_example)))

  cmt <- createCoMorTraxObject(CoMorTrax_AD_T2D_example) |>
    encodeDiseaseLabels(
      diseaseCol   = c("disease_AD", "disease_T2D"),
      controlLabel = "control",
      verbose      = TRUE
    ) |>
    runComorbidityDE(
      method      = "wilcoxon",
      cellTypeCol = "cellType",
      verbose     = TRUE
    ) |>
    quantifyInteractionEffects(
      nPermutations = 500,
      verbose       = TRUE
    ) |>
    scoreCellVulnerability(
      cellTypeCol         = "cellType",
      includePathwayScore = FALSE,
      verbose             = TRUE
    )

  summaryCoMorTrax(cmt)
  cmt
}

## Run the quick example immediately
cat("\n=== QUICK START: Running on bundled example data ===\n")
cmt_example <- quick_start_example()
