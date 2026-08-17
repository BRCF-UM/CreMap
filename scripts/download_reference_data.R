#!/usr/bin/env Rscript
## Build CreMAP default reference Seurat objects and save them under reference/*.rds
## (one-time; requires network + Bioconductor packages TabulaMurisSenisData, scRNAseq, etc.).
##
## Run from the project root:
##   Rscript scripts/download_reference_data.R
##
## Output:
##   reference/mouse_spleen_reference.rds
##   reference/human_spleen_reference.rds
##
## The Shiny app loads these files first (no ExperimentHub at runtime).

options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = "FALSE")

ca <- commandArgs(trailingOnly = FALSE)
fn <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
if (length(fn) != 1L) {
  stop("Run from the repo root: Rscript scripts/download_reference_data.R", call. = FALSE)
}
root <- normalizePath(file.path(dirname(fn), ".."))
setwd(root)

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

en <- new.env(parent = globalenv())
sys.source(file.path(root, "R", "real_spleen_data.R"), envir = en)

en$hub_allow_downloads_no_prompt()
refd <- en$cremap_reference_dir()
dir.create(refd, recursive = TRUE, showWarnings = FALSE)
message("Reference directory: ", refd)
message("Building mouse object from Tabula Muris Senis (may download via ExperimentHub)…")
mouse <- en$load_mouse_spleen_from_hub()
if (is.null(mouse)) {
  stop("Mouse hub load failed (TabulaMurisSenisData / network).", call. = FALSE)
}
message("Building human object from He organ atlas (may download)…")
human <- en$load_human_spleen_from_hub()
if (is.null(human)) {
  stop("Human hub load failed (scRNAseq / network).", call. = FALSE)
}

paths <- en$cremap_reference_paths()
message("Saving: ", paths$mouse)
saveRDS(mouse, paths$mouse, compress = "xz")
message("Saving: ", paths$human)
saveRDS(human, paths$human, compress = "xz")
message("Done. Restart the CreMAP app; it will read these files without calling hubs.")
