\name{CoMorTrax_AD_T2D_example}
\docType{data}
\alias{CoMorTrax_AD_T2D_example}
\title{Example Single-Cell Dataset: Alzheimer Disease + Type 2 Diabetes Comorbidity}
\description{
  A compact \code{SingleCellExperiment} object for CoMorTrax demonstration,
  representative of comorbid Alzheimer Disease (AD) and Type 2 Diabetes (T2D).

  This dataset is simulated to reflect the statistical properties reported in
  two landmark studies:
  \itemize{
    \item Mathys et al. (2019) \emph{Nature} 571:332-337 -- 80,660
      single-nucleus RNA-seq profiles from prefrontal cortex of 48 AD
      and control subjects. GEO accession: \strong{GSE125050}.
    \item Segerstolpe et al. (2016) \emph{Cell Metabolism} 24:593-607 --
      2,209 single-cell RNA-seq profiles from human pancreatic islets
      of 6 T2D and 4 control donors. ArrayExpress: \strong{E-MTAB-5061}.
  }

  See \code{inst/extdata/real_dataset_analysis.R} for complete code to
  download, process, and analyse the full GEO datasets with CoMorTrax.
}
\usage{data(CoMorTrax_AD_T2D_example)}
\format{
  A \code{\link[SingleCellExperiment]{SingleCellExperiment}} with:
  \describe{
    \item{Dimensions}{150 genes x 400 cells}
    \item{Assays}{\code{counts} (raw integer), \code{logcounts} (log-normalised)}
    \item{colData}{Columns: \code{disease}, \code{cellType}, \code{batch},
      \code{is_control}, \code{disease_AD}, \code{disease_T2D}}
    \item{reducedDims}{\code{PCA} (10 components), \code{UMAP} (2 dimensions)}
    \item{metadata}{Includes \code{sim_params} (injected gene identities)
      and \code{dataset_info} (GEO accessions and references)}
  }
}
\details{
  \strong{Real dataset statistics (for reference):}

  \tabular{lrrl}{
    \strong{Disease} \tab \strong{Cells} \tab \strong{Genes} \tab
    \strong{Accession} \cr
    Alzheimer Disease \tab 80,660 \tab 17,926 \tab GSE125050 \cr
    Type 2 Diabetes   \tab  2,209 \tab 23,271 \tab E-MTAB-5061 \cr
  }

  \strong{Cell types (real datasets):}
  \itemize{
    \item AD (brain): Excitatory Neuron, Inhibitory Neuron, Microglia,
      Astrocyte, OPC, Endothelial
    \item T2D (pancreas): Alpha, Beta, Delta, Gamma, Acinar, Ductal
  }

  \strong{To download and run on the real GEO data:}
  \preformatted{
  source(system.file("extdata", "real_dataset_analysis.R",
                     package = "CoMorTrax"))
  }
}
\source{
  Simulated to reflect:
  \itemize{
    \item Mathys et al. (2019). Single-cell transcriptomic analysis of
      Alzheimer's disease. \emph{Nature}, 571, 332-337.
      \doi{10.1038/s41586-019-1195-2}
    \item Segerstolpe et al. (2016). Single-Cell Transcriptome Profiling
      of Human Pancreatic Islets in Health and Type 2 Diabetes.
      \emph{Cell Metabolism}, 24, 593-607.
      \doi{10.1016/j.cmet.2016.08.020}
  }
}
\seealso{
  \code{\link{createCoMorTraxObject}},
  \code{\link{encodeDiseaseLabels}},
  \code{\link{quantifyInteractionEffects}}
}
\examples{
library(CoMorTrax)
data(CoMorTrax_AD_T2D_example)

# Inspect the dataset
CoMorTrax_AD_T2D_example

# View dataset metadata
info <- S4Vectors::metadata(CoMorTrax_AD_T2D_example)$dataset_info
cat("GEO accessions:", paste(info$geo_accession, collapse=", "), "\n")
cat("True synergistic genes:", length(
  S4Vectors::metadata(CoMorTrax_AD_T2D_example)$sim_params$syn_genes), "\n")

# Run the full CoMorTrax pipeline
\dontrun{
cmt <- createCoMorTraxObject(CoMorTrax_AD_T2D_example) |>
  encodeDiseaseLabels(diseaseCol = c("disease_AD", "disease_T2D")) |>
  runComorbidityDE(method = "wilcoxon") |>
  quantifyInteractionEffects(nPermutations = 500)

summaryCoMorTrax(cmt)
}
}
\keyword{datasets}
