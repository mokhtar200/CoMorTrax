## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.align = "center",
  warning  = FALSE,
  message  = FALSE,
  eval     = FALSE
)

## ----install, eval=FALSE------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("CoMorTrax")

## ----install-dev, eval=FALSE--------------------------------------------------
#  BiocManager::install("AhmedSalem/CoMorTrax")

## ----load---------------------------------------------------------------------
#  library(CoMorTrax)
#  library(SingleCellExperiment)
#  library(S4Vectors)
#  library(SummarizedExperiment)

## ----simulate-----------------------------------------------------------------
#  set.seed(42)
#  sce <- simulateCoMorTrax(
#    n_cells    = 400,
#    n_genes    = 300,
#    n_syn      = 20,    # genes synergistically amplified under comorbidity
#    n_ant      = 15,    # genes antagonistically suppressed under comorbidity
#    diseases   = c("AD", "T2D"),
#    cell_types = c("Neuron", "Microglia", "Astrocyte")
#  )
#  
#  # True interaction gene identities are stored in metadata
#  true_syn <- metadata(sce)$sim_params$syn_genes
#  true_ant <- metadata(sce)$sim_params$ant_genes
#  cat("True synergistic genes:", paste(head(true_syn, 5), collapse = ", "), "\n")
#  cat("True antagonistic genes:", paste(head(true_ant, 5), collapse = ", "), "\n")
#  
#  sce

## ----create-------------------------------------------------------------------
#  cmt <- createCoMorTraxObject(sce)
#  cmt

## ----preprocess---------------------------------------------------------------
#  cmt <- preprocessCoMorTrax(
#    cmt,
#    normMethod   = "logNorm",
#    nHVG         = 2000,
#    nPCs         = 30,
#    runUMAP      = TRUE,
#    batchCorrect = FALSE,
#    verbose      = TRUE
#  )

## ----encode-------------------------------------------------------------------
#  cmt <- encodeDiseaseLabels(
#    cmt,
#    diseaseCol   = c("disease_AD", "disease_T2D"),
#    controlLabel = "control",
#    verbose      = TRUE
#  )
#  
#  # Inspect the binary disease matrix
#  head(diseaseMatrix(cmt))
#  
#  # All automatically enumerated comparisons
#  listComorbidityComparisons(cmt)

## ----de-----------------------------------------------------------------------
#  cmt <- runComorbidityDE(
#    cmt,
#    method      = "wilcoxon",
#    cellTypeCol = "cellType",
#    verbose     = TRUE
#  )
#  
#  # Access results for a specific comparison
#  de <- getDeTable(cmt)
#  head(de[de$significant, ])

## ----interaction--------------------------------------------------------------
#  cmt <- quantifyInteractionEffects(
#    cmt,
#    nPermutations = 1000,
#    synThreshold  = 0.25,
#    antThreshold  = -0.25,
#    verbose       = TRUE
#  )
#  
#  # Access results
#  int_df <- interactionScores(cmt)
#  table(int_df$classification)
#  
#  # Retrieve gene lists
#  syn_genes <- getSynergisticGenes(cmt)
#  ant_genes <- getAntagonisticGenes(cmt)
#  cat("Synergistic genes detected:", length(syn_genes), "\n")
#  
#  # Benchmark: compare detected vs true injected
#  true_syn <- metadata(sce)$sim_params$syn_genes
#  cat("Precision:", round(mean(syn_genes %in% true_syn), 3), "\n")
#  cat("Recall:",    round(mean(true_syn %in% syn_genes), 3), "\n")

## ----pathways-----------------------------------------------------------------
#  cmt <- runComorbidityPathways(
#    cmt,
#    genesets = c("GO_BP", "KEGG"),
#    organism = "hsa",
#    pCutoff  = 0.05,
#    verbose  = TRUE
#  )
#  
#  # Compare enriched pathways between synergistic and antagonistic sets
#  comparePathways(cmt, topN = 10)

## ----vulnerability------------------------------------------------------------
#  cmt <- scoreCellVulnerability(
#    cmt,
#    cellTypeCol         = "cellType",
#    includePathwayScore = TRUE,
#    weights             = list(synergistic  = 0.35,
#                                antagonistic = 0.25,
#                                pathway      = 0.25,
#                                magnitude    = 0.15),
#    verbose             = TRUE
#  )
#  
#  cellVulnerability(cmt)

## ----networks-----------------------------------------------------------------
#  cmt <- buildComorbidityNetworks(
#    cmt,
#    method          = "correlation",
#    corrThreshold   = 0.6,
#    computeRewiring = TRUE,
#    verbose         = TRUE
#  )
#  
#  # Retrieve hub genes for the comorbid state
#  hubs <- getNetworkHubs(cmt, state = "AD_T2D", topN = 10)
#  print(hubs)
#  
#  # Network rewiring between control and comorbid state
#  rew <- computeRewiring(cmt, state1 = "control", state2 = "AD_T2D")
#  cat("Edge Jaccard (control vs AD+T2D):", round(rew$jaccardEdges, 3), "\n")
#  cat("Gained edges:", rew$gainedEdges, "| Lost edges:", rew$lostEdges, "\n")

## ----ml-----------------------------------------------------------------------
#  cmt <- trainComorbidityClassifier(
#    cmt,
#    method      = "randomForest",
#    featureType = "interaction",
#    cvFolds     = 5L,
#    verbose     = TRUE
#  )
#  
#  clf <- classifier(cmt)
#  cat("CV AUC:", round(clf$performance$auc, 3), "\n")
#  head(clf$feature_importance)

## ----plots, eval=FALSE--------------------------------------------------------
#  # UMAP colored by disease state
#  plotCoMorTrax(cmt, type = "umap", colorBy = "disease_state")
#  
#  # Interaction heatmap of top genes
#  plotCoMorTrax(cmt, type = "interactionHeatmap", topN = 50)
#  
#  # Volcano: interaction delta vs -log10(p)
#  plotCoMorTrax(cmt, type = "volcanoInteraction", labelTopN = 15)
#  
#  # Cell vulnerability ranking
#  plotCoMorTrax(cmt, type = "vulnerabilityRank")
#  
#  # Classifier ROC curve
#  plotCoMorTrax(cmt, type = "classifierROC")
#  
#  # Full overview (all key plots)
#  plotCoMorTrax(cmt, type = "overview")

## ----validate-----------------------------------------------------------------
#  # Known disease genes from literature
#  known <- c("APOE", "CLU", "TREM2", "BIN1",   # AD GWAS hits
#             "TCF7L2", "PPARG", "KCNJ11")        # T2D GWAS hits
#  
#  val <- validateCoMorTrax(
#    cmt,
#    knownGenes   = known,
#    nBootstrap   = 100,
#    verbose      = TRUE
#  )
#  
#  cat("Hypergeometric p-value:", val$biological$pvalue, "\n")
#  cat("Bootstrap stability:   ", round(val$robustness, 3), "\n")

## ----export-------------------------------------------------------------------
#  exportCoMorTrax(
#    cmt,
#    outDir  = "./CoMorTrax_results",
#    formats = c("csv", "rds"),
#    verbose = TRUE
#  )

## ----full-example-------------------------------------------------------------
#  library(CoMorTrax)
#  set.seed(2026)
#  
#  # Simulate
#  sce <- simulateCoMorTrax(n_cells = 200, n_genes = 150,
#                            n_syn = 10, n_ant = 8,
#                            diseases = c("AD", "T2D"))
#  
#  # Full pipeline in one chain
#  cmt <- createCoMorTraxObject(sce)                           |>
#    encodeDiseaseLabels(diseaseCol = c("disease_AD",
#                                        "disease_T2D"))       |>
#    runComorbidityDE(method = "wilcoxon")                     |>
#    quantifyInteractionEffects(nPermutations = 200)           |>
#    scoreCellVulnerability(includePathwayScore = FALSE)       |>
#    buildComorbidityNetworks(maxGenes = 100)                  |>
#    trainComorbidityClassifier(method = "randomForest",
#                                cvFolds = 3L)
#  
#  # Summary
#  summaryCoMorTrax(cmt)

## ----session------------------------------------------------------------------
#  sessionInfo()

