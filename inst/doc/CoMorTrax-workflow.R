## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  warning = FALSE,
  message = FALSE,
  eval = FALSE
)

## ----install-bioc, eval=FALSE-------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("CoMorTrax")

## ----install-github, eval=FALSE-----------------------------------------------
#  BiocManager::install("AhmedSalem/CoMorTrax")

## ----load---------------------------------------------------------------------
#  library(CoMorTrax)
#  library(SingleCellExperiment)
#  library(ggplot2)

## ----simulate-----------------------------------------------------------------
#  set.seed(42)
#  sce <- simulateCoMorTrax(
#    n_cells    = 1000,
#    n_genes    = 500,
#    n_syn      = 30,     # synergistic genes (non-additive amplification)
#    n_ant      = 20,     # antagonistic genes (non-additive suppression)
#    cell_types = c("Neuron", "Microglia", "Astrocyte", "OPC"),
#    diseases   = c("AD", "T2D", "Dementia"),
#    seed       = 42
#  )
#  
#  sce

## ----show-sim, echo=FALSE-----------------------------------------------------
#  # Expected output:
#  # class: SingleCellExperiment
#  # dim: 500 1000
#  # metadata(1): params
#  # ...
#  # colData names(8): cell_type disease_AD disease_T2D disease_Dementia ...

## ----create-from-sce----------------------------------------------------------
#  cmt <- createCoMorTraxObject(
#    data   = sce,
#    params = list(
#      minCells   = 3,
#      minGenes   = 200,
#      maxMitoPCT = 20
#    )
#  )
#  cmt

## ----create-from-matrix, eval=FALSE-------------------------------------------
#  # Assuming 'counts_mat' is a genes × cells matrix
#  # and 'cell_meta' is a DataFrame with disease columns
#  
#  cmt <- createCoMorTraxObject(
#    data      = counts_mat,
#    colData   = cell_meta,
#    rowData   = gene_meta
#  )

## ----preprocess---------------------------------------------------------------
#  cmt <- preprocessCoMorTrax(
#    object         = cmt,
#    normMethod     = "logNorm",    # "logNorm", "scran", or "scTransform"
#    nHVG           = 2000,
#    nPCs           = 30,
#    batchCorrect   = FALSE,        # Set TRUE if 'batch' column present in colData
#    batchMethod    = "harmony",    # "harmony" or "fastMNN"
#    runUMAP        = TRUE,
#    verbose        = TRUE
#  )

## ----encode-binary------------------------------------------------------------
#  cmt <- encodeDiseaseLabels(
#    object       = cmt,
#    diseaseCol   = c("disease_AD", "disease_T2D", "disease_Dementia"),
#    controlLabel = "control",
#    verbose      = TRUE
#  )

## ----encode-categorical, eval=FALSE-------------------------------------------
#  # If colData has a 'diagnosis' column with values like "AD;T2D", "Control", etc.
#  cmt <- encodeDiseaseLabels(
#    object       = cmt,
#    diseaseCol   = "diagnosis",
#    separator    = ";",
#    controlLabel = "Control"
#  )

## ----inspect-disease----------------------------------------------------------
#  # View the binary disease matrix (cells × diseases)
#  head(diseaseMatrix(cmt))
#  
#  # View all comparison groups defined
#  str(diseaseGroups(cmt))
#  
#  # Enumerate all pairwise + comorbid comparisons
#  listComorbidityComparisons(cmt)

## ----run-de-------------------------------------------------------------------
#  cmt <- runComorbidityDE(
#    object      = cmt,
#    cellTypeCol = "cell_type",
#    comparisons = "auto",      # Run all enumerated comparisons
#    method      = "wilcoxon",  # "wilcoxon", "edgeR", or "MAST"
#    minCells    = 10,
#    logFCcutoff = 0.5,
#    adjPcutoff  = 0.05,
#    verbose     = TRUE
#  )

## ----access-de----------------------------------------------------------------
#  # Get DE table for a specific comparison and cell type
#  de_table <- getDeTable(
#    cmt,
#    comparison = "AD_T2D_vs_Control",
#    cellType   = "Neuron"
#  )
#  
#  head(de_table)
#  #     gene     logFC    pvalue      padj significant direction
#  # 1  APOE4  2.341    1.2e-08  3.4e-06        TRUE          up
#  # 2  CLU    1.892    4.5e-07  8.9e-05        TRUE          up
#  # 3  TREM2 -1.203    2.3e-06  3.1e-04        TRUE        down
#  # ...

## ----interaction-effects------------------------------------------------------
#  cmt <- quantifyInteractionEffects(
#    object         = cmt,
#    disease1       = "AD",
#    disease2       = "T2D",
#    cellTypeCol    = "cell_type",
#    synThreshold   = 0.5,      # |Δ| > 0.5 classified synergistic
#    antThreshold   = -0.5,     # Δ < -0.5 classified antagonistic
#    nPermutations  = 1000,     # Permutation-based p-values
#    adjPcutoff     = 0.05,
#    verbose        = TRUE
#  )

## ----interaction-genes--------------------------------------------------------
#  # Synergistic genes (amplified in comorbidity beyond expectation)
#  syn_genes <- getSynergisticGenes(cmt, cellType = "Neuron")
#  cat("Synergistic genes in Neurons:", length(syn_genes), "\n")
#  
#  # Antagonistic genes
#  ant_genes <- getAntagonisticGenes(cmt, cellType = "Microglia")
#  
#  # Additive genes
#  add_genes <- getAdditiveGenes(cmt)
#  
#  # Full interaction table
#  int_table <- interactionScores(cmt)
#  head(int_table[order(abs(int_table$interactionDelta), decreasing=TRUE), ])

## ----pathway-analysis---------------------------------------------------------
#  cmt <- runComorbidityPathways(
#    object     = cmt,
#    genesets   = c("GO_BP", "KEGG", "Reactome"),
#    organism   = "hsa",        # "hsa" for human, "mmu" for mouse
#    minGSSize  = 10,
#    maxGSSize  = 500,
#    pCutoff    = 0.05,
#    verbose    = TRUE
#  )

## ----pathway-comparison-------------------------------------------------------
#  # Side-by-side comparison: synergistic vs. antagonistic pathways
#  path_comparison <- comparePathways(
#    cmt,
#    cellType = "Neuron",
#    topN     = 20
#  )
#  
#  head(path_comparison)

## ----cell-vulnerability-------------------------------------------------------
#  cmt <- scoreCellVulnerability(
#    object      = cmt,
#    cellTypeCol = "cell_type",
#    weights     = list(
#      synergistic  = 0.35,   # weight for number of synergistic genes
#      antagonistic = 0.25,   # weight for number of antagonistic genes
#      pathway      = 0.25,   # weight for pathway disruption score
#      magnitude    = 0.15    # weight for mean |interaction delta|
#    ),
#    normalize   = TRUE
#  )
#  
#  # View vulnerability rankings
#  vul_scores <- cellVulnerability(cmt)
#  vul_scores[order(vul_scores$vulnerabilityScore, decreasing=TRUE), ]

## ----networks-----------------------------------------------------------------
#  cmt <- buildComorbidityNetworks(
#    object      = cmt,
#    method      = "correlation",   # "correlation", "WGCNA", or "GENIE3"
#    states      = c("AD", "T2D", "AD_T2D", "Control"),
#    corMethod   = "pearson",
#    corCutoff   = 0.6,
#    pCutoff     = 0.01,
#    verbose     = TRUE
#  )

## ----rewiring-----------------------------------------------------------------
#  # Compare network topology between two states
#  rewiring <- computeRewiring(
#    cmt,
#    state1   = "Control",
#    state2   = "AD_T2D",
#    topHubs  = 20
#  )
#  
#  # Get hub genes
#  hubs <- getNetworkHubs(cmt, state = "AD_T2D", top = 20)
#  print(hubs)

## ----ml-classifier------------------------------------------------------------
#  cmt <- trainComorbidityClassifier(
#    object          = cmt,
#    method          = "randomForest",   # "randomForest", "lasso", or "elasticNet"
#    featureSource   = "interaction",    # use interaction genes as features
#    nTrees          = 500,
#    kFold           = 5,
#    verbose         = TRUE
#  )
#  
#  # View classifier performance
#  clf_results <- classifier(cmt)
#  cat("Cross-validated AUC:", round(clf_results$cv_auc, 3), "\n")
#  cat("Balanced Accuracy:", round(clf_results$balanced_accuracy, 3), "\n")

## ----spatial, eval=FALSE------------------------------------------------------
#  cmt <- mapSpatialComorbidity(
#    object       = cmt,
#    spatialData  = spatial_sce,   # SpatialExperiment object
#    method       = "mean_expr",   # "mean_expr" or "ssgsea"
#    geneSet      = "synergistic"
#  )
#  
#  # Visualize on tissue
#  plotSpatialComorbidity(cmt, title = "AD+T2D Synergistic Score")

## ----plot-umap----------------------------------------------------------------
#  plotCoMorTrax(cmt, type = "umap",
#    colorBy = "disease_state",
#    title   = "Cells Colored by Comorbidity State"
#  )

## ----plot-heatmap-------------------------------------------------------------
#  plotCoMorTrax(cmt, type = "interactionHeatmap",
#    cellType = "Neuron",
#    topN     = 50,
#    title    = "Top 50 Interaction Genes — Neurons"
#  )

## ----plot-volcano-------------------------------------------------------------
#  plotCoMorTrax(cmt, type = "volcanoInteraction",
#    cellType     = "Neuron",
#    labelTopN    = 15,
#    synColor     = "#E63946",
#    antColor     = "#457B9D",
#    addColor     = "#6c757d"
#  )

## ----plot-pathways------------------------------------------------------------
#  plotCoMorTrax(cmt, type = "pathwayDotplot",
#    cellType  = "Neuron",
#    geneSet   = "synergistic",
#    topN      = 15
#  )

## ----plot-vulnerability-------------------------------------------------------
#  plotCoMorTrax(cmt, type = "vulnerabilityRank",
#    fillColor = "#E63946",
#    title     = "Cell-Type Vulnerability to AD+T2D Comorbidity"
#  )

## ----plot-roc-----------------------------------------------------------------
#  plotCoMorTrax(cmt, type = "classifierROC",
#    title = "Comorbidity Classifier Performance"
#  )

## ----plot-overview------------------------------------------------------------
#  plotCoMorTrax(cmt, type = "overview")

## ----validate-----------------------------------------------------------------
#  val_results <- validateCoMorTrax(
#    object       = cmt,
#    knownGenes   = c("APOE", "CLU", "TREM2", "BIN1",  # AD genes
#                     "TCF7L2", "PPARG", "KCNJ11"),      # T2D genes
#    nBootstrap   = 100,
#    verbose      = TRUE
#  )
#  
#  # Print validation summary
#  print(val_results)

## ----export-------------------------------------------------------------------
#  exportCoMorTrax(
#    object    = cmt,
#    outDir    = "./CoMorTrax_results",
#    format    = c("csv", "rds"),
#    compress  = TRUE
#  )

## ----full-example-------------------------------------------------------------
#  library(CoMorTrax)
#  set.seed(2024)
#  
#  # 1. Simulate data with known comorbidity signal
#  sce <- simulateCoMorTrax(
#    n_cells = 800, n_genes = 400,
#    n_syn = 25, n_ant = 15,
#    diseases = c("AD", "T2D")
#  )
#  
#  # 2–3. Create object + preprocess
#  cmt <- createCoMorTraxObject(sce) |>
#    preprocessCoMorTrax(normMethod = "logNorm", nHVG = 1500)
#  
#  # 4. Encode disease labels
#  cmt <- encodeDiseaseLabels(cmt, diseaseCol = c("disease_AD","disease_T2D"))
#  
#  # 5. Comorbidity-aware DE
#  cmt <- runComorbidityDE(cmt, method = "wilcoxon", minCells = 5)
#  
#  # 6. Interaction effects
#  cmt <- quantifyInteractionEffects(
#    cmt, disease1 = "AD", disease2 = "T2D",
#    nPermutations = 500
#  )
#  
#  # 7. Pathway enrichment
#  cmt <- runComorbidityPathways(cmt, genesets = "GO_BP")
#  
#  # 8. Vulnerability scoring
#  cmt <- scoreCellVulnerability(cmt)
#  
#  # 9. Network modeling
#  cmt <- buildComorbidityNetworks(cmt, method = "correlation")
#  
#  # 10. ML classifier
#  cmt <- trainComorbidityClassifier(cmt, method = "randomForest")
#  
#  # Full summary
#  summaryCoMorTrax(cmt)
#  
#  # Export
#  exportCoMorTrax(cmt, outDir = "./my_results")

## ----session-info-------------------------------------------------------------
#  sessionInfo()

