## =============================================================================
## AllGenerics.R - S4 Generic Definitions for CoMorTrax
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

## --- Accessors ---------------------------------------------------------------

#' @title Accessor: Disease Matrix
#' @description Get or set the binary disease label matrix.
#' @param x A \code{\link{CoMorTraxObject}}.
#' @param value A \code{DataFrame} of binary disease indicators.
#' @return A \code{DataFrame} (getter) or a modified object (setter).
#' @export
setGeneric("diseaseMatrix",
    function(x) standardGeneric("diseaseMatrix"))

#' @export
#' @rdname diseaseMatrix
setGeneric("diseaseMatrix<-",
    function(x, value) standardGeneric("diseaseMatrix<-"))

#' @title Accessor: Disease Groups
#' @description Get or set disease group membership lists.
#' @param x A \code{\link{CoMorTraxObject}}.
#' @param value A named \code{list} mapping group IDs to cell indices.
#' @return A named \code{list} (getter) or a modified object (setter).
#' @export
setGeneric("diseaseGroups",
    function(x) standardGeneric("diseaseGroups"))

#' @export
#' @rdname diseaseGroups
setGeneric("diseaseGroups<-",
    function(x, value) standardGeneric("diseaseGroups<-"))

#' @title Accessor: DE Results
#' @description Get or set differential expression results.
#' @param x A \code{\link{CoMorTraxObject}}.
#' @param comparison Optional character; name of comparison to retrieve.
#' @param value A named \code{list} of DE result data.frames.
#' @return A named \code{list} or single \code{data.frame}.
#' @export
setGeneric("deResults",
    function(x, comparison = NULL) standardGeneric("deResults"))

#' @export
#' @rdname deResults
setGeneric("deResults<-",
    function(x, value) standardGeneric("deResults<-"))

#' @title Accessor: Interaction Scores
#' @description Get or set gene-level interaction effect scores.
#' @param x A \code{\link{CoMorTraxObject}}.
#' @param value A \code{DataFrame} of interaction scores.
#' @return A \code{DataFrame} (getter) or modified object (setter).
#' @export
setGeneric("interactionScores",
    function(x) standardGeneric("interactionScores"))

#' @export
#' @rdname interactionScores
setGeneric("interactionScores<-",
    function(x, value) standardGeneric("interactionScores<-"))

#' @title Accessor: Pathway Results
#' @description Get or set pathway enrichment results.
#' @param x A \code{\link{CoMorTraxObject}}.
#' @param geneSet Optional character; one of \code{"synergistic"},
#'   \code{"antagonistic"}, \code{"additive"}.
#' @param value A named \code{list} of enrichment results.
#' @return A \code{list} (getter) or modified object (setter).
#' @export
setGeneric("pathwayResults",
    function(x, geneSet = NULL) standardGeneric("pathwayResults"))

#' @export
#' @rdname pathwayResults
setGeneric("pathwayResults<-",
    function(x, value) standardGeneric("pathwayResults<-"))

#' @title Accessor: Cell Vulnerability Scores
#' @description Get or set cell-type vulnerability scores.
#' @param x A \code{\link{CoMorTraxObject}}.
#' @param value A \code{DataFrame} of vulnerability scores.
#' @return A \code{DataFrame} (getter) or modified object (setter).
#' @export
setGeneric("cellVulnerability",
    function(x) standardGeneric("cellVulnerability"))

#' @export
#' @rdname cellVulnerability
setGeneric("cellVulnerability<-",
    function(x, value) standardGeneric("cellVulnerability<-"))

#' @title Accessor: Networks
#' @description Get or set disease-state gene regulatory networks.
#' @param x A \code{\link{CoMorTraxObject}}.
#' @param state Optional character; disease state name.
#' @param value A named \code{list} of \code{igraph} objects.
#' @return A \code{list} or single \code{igraph} object.
#' @export
setGeneric("networks",
    function(x, state = NULL) standardGeneric("networks"))

#' @export
#' @rdname networks
setGeneric("networks<-",
    function(x, value) standardGeneric("networks<-"))

#' @title Accessor: Classifier
#' @description Get or set the trained comorbidity classifier.
#' @param x A \code{\link{CoMorTraxObject}}.
#' @param value A named \code{list} of model and performance objects.
#' @return A named \code{list} or modified object.
#' @export
setGeneric("classifier",
    function(x) standardGeneric("classifier"))

#' @export
#' @rdname classifier
setGeneric("classifier<-",
    function(x, value) standardGeneric("classifier<-"))

## --- Analysis generics -------------------------------------------------------

#' @title Create a CoMorTraxObject
#' @description Construct the core \code{\link{CoMorTraxObject}} from a
#'   \code{SingleCellExperiment} or raw count matrix.
#' @param x Input data: a \code{SingleCellExperiment}, \code{matrix}, or
#'   \code{dgCMatrix}.
#' @param ... Additional arguments passed to methods.
#' @return A \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("createCoMorTraxObject",
    function(x, ...) standardGeneric("createCoMorTraxObject"))

#' @title Preprocess CoMorTrax Data
#' @description Quality control, normalization, dimensionality reduction,
#'   and optional batch correction.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return An updated \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("preprocessCoMorTrax",
    function(object, ...) standardGeneric("preprocessCoMorTrax"))

#' @title Encode Disease Labels
#' @description Build the multi-label binary disease matrix from metadata.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return An updated \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("encodeDiseaseLabels",
    function(object, ...) standardGeneric("encodeDiseaseLabels"))

#' @title Run Comorbidity-Aware Differential Expression
#' @description Perform DE analysis for all disease combination comparisons.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return An updated \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("runComorbidityDE",
    function(object, ...) standardGeneric("runComorbidityDE"))

#' @title Quantify Interaction Effects
#' @description Compute per-gene interaction scores and classify genes.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return An updated \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("quantifyInteractionEffects",
    function(object, ...) standardGeneric("quantifyInteractionEffects"))

#' @title Run Comorbidity Pathway Analysis
#' @description Enrichment analysis for synergistic/antagonistic gene sets.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return An updated \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("runComorbidityPathways",
    function(object, ...) standardGeneric("runComorbidityPathways"))

#' @title Score Cell-Type Vulnerability
#' @description Rank cell types by comorbidity-driven molecular disruption.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return An updated \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("scoreCellVulnerability",
    function(object, ...) standardGeneric("scoreCellVulnerability"))

#' @title Build Comorbidity Networks
#' @description Construct gene regulatory networks per disease state.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return An updated \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("buildComorbidityNetworks",
    function(object, ...) standardGeneric("buildComorbidityNetworks"))

#' @title Train Comorbidity Classifier
#' @description Train and evaluate a machine learning comorbidity predictor.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return An updated \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("trainComorbidityClassifier",
    function(object, ...) standardGeneric("trainComorbidityClassifier"))

#' @title Map Spatial Comorbidity
#' @description Project comorbidity signatures onto spatial transcriptomics.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return An updated \code{\link{CoMorTraxObject}}.
#' @export
setGeneric("mapSpatialComorbidity",
    function(object, ...) standardGeneric("mapSpatialComorbidity"))

#' @title Plot CoMorTrax Results
#' @description Unified plotting interface for all CoMorTrax output types.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param type Character; the plot type to produce.
#' @param ... Additional arguments.
#' @return A \code{ggplot} object or list thereof.
#' @export
setGeneric("plotCoMorTrax",
    function(object, type, ...) standardGeneric("plotCoMorTrax"))

#' @title Summary of CoMorTraxObject
#' @description Print a structured overview of completed analysis modules.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param ... Additional arguments.
#' @return Invisibly returns \code{object}.
#' @export
setGeneric("summaryCoMorTrax",
    function(object, ...) standardGeneric("summaryCoMorTrax"))

#' @title Export CoMorTrax Results
#' @description Write analysis outputs to disk in standard formats.
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param outDir Character; output directory path.
#' @param ... Additional arguments.
#' @return Invisibly returns \code{outDir}.
#' @export
setGeneric("exportCoMorTrax",
    function(object, outDir, ...) standardGeneric("exportCoMorTrax"))
