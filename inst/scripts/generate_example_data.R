## CoMorTrax — Example Data Generation Script
## Author: Ahmed Mokhtar Ramzy Salem
## Description: Generates example datasets for CoMorTrax vignettes and testing.
##
## Run this script ONCE to regenerate the data stored in inst/extdata/.
## Usage: Rscript inst/scripts/generate_example_data.R

message("CoMorTrax Example Data Generator")
message("Author: Ahmed Mokhtar Ramzy Salem")
message("=========================================")

library(CoMorTrax)

out_dir <- system.file("extdata", package = "CoMorTrax")
if (out_dir == "") out_dir <- file.path("inst", "extdata")

# ── Dataset 1: AD + T2D Comorbidity (small, for testing) ─────────────────────
message("\n[1/3] Simulating AD+T2D comorbidity dataset (n=500 cells)...")
set.seed(2024)

sce_small <- simulateCoMorTrax(
    n_cells    = 500,
    n_genes    = 300,
    n_syn      = 20,
    n_ant      = 15,
    cell_types = c("Neuron", "Microglia", "Astrocyte"),
    diseases   = c("AD", "T2D"),
    seed       = 2024
)

saveRDS(sce_small,
        file = file.path(out_dir, "example_sce_small.rds"),
        compress = "xz")
message("  Saved: example_sce_small.rds (", ncol(sce_small), " cells x ",
        nrow(sce_small), " genes)")

# ── Dataset 2: Three-way Comorbidity (AD + T2D + Dementia) ───────────────────
message("\n[2/3] Simulating three-way comorbidity dataset (n=800 cells)...")
set.seed(2024)

sce_3way <- simulateCoMorTrax(
    n_cells    = 800,
    n_genes    = 400,
    n_syn      = 30,
    n_ant      = 20,
    cell_types = c("Neuron", "Microglia", "Astrocyte", "OPC"),
    diseases   = c("AD", "T2D", "Dementia"),
    seed       = 2024
)

saveRDS(sce_3way,
        file = file.path(out_dir, "example_sce_3way.rds"),
        compress = "xz")
message("  Saved: example_sce_3way.rds (", ncol(sce_3way), " cells x ",
        nrow(sce_3way), " genes)")

# ── Dataset 3: Known synergistic gene list (for validation) ──────────────────
message("\n[3/3] Generating known gene sets for validation...")

# AD and T2D genes curated from literature (GWAS/DisGeNET)
known_ad_genes <- c(
    "APOE", "CLU", "PICALM", "CR1", "BIN1", "MS4A6A", "CD33",
    "ABCA7", "EPHA1", "CD2AP", "SORL1", "TREM2", "FERMT2",
    "SLC24A4", "HLA-DRB5", "PTK2B", "ZCWPW1", "CELF1", "INPP5D",
    "MEF2C", "NME8", "IGHV1-67", "CASS4", "FERMT2", "DSG2"
)

known_t2d_genes <- c(
    "TCF7L2", "PPARG", "KCNJ11", "NOTCH2", "WFS1", "CDKAL1",
    "IGF2BP2", "CDKN2A", "CDKN2B", "HHEX", "LOC387761",
    "SLC30A8", "JAZF1", "CDC123", "TSPAN8", "THADA", "ADAMTS9",
    "CAMK1D", "FANCL", "DGKB", "MNTR1B", "ZBED3", "KLF14",
    "DEXI", "CENTD2", "RASGRP1", "C2CD4A", "VPS26A", "HNF1B"
)

# Save as a named list
known_gene_sets <- list(
    AD_genes  = known_ad_genes,
    T2D_genes = known_t2d_genes
)

saveRDS(known_gene_sets,
        file = file.path(out_dir, "known_gene_sets.rds"),
        compress = "xz")
message("  Saved: known_gene_sets.rds (",
        length(known_ad_genes), " AD genes, ",
        length(known_t2d_genes), " T2D genes)")

# ── Summary ──────────────────────────────────────────────────────────────────
message("\n=========================================")
message("Data generation complete.")
message("Files written to: ", out_dir)
message("\nFiles:")
for (f in list.files(out_dir, full.names=FALSE)) {
    message("  - ", f)
}
