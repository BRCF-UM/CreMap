#!/usr/bin/env Rscript
## Install packages into the project renv library.
##
## requireNamespace() is NOT reliable here: renv can hide the system library while
## install_deps still "sees" global packages and skips installing. We only trust
## packages present under the project library path.

options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = "FALSE")

ca <- commandArgs(trailingOnly = FALSE)
fn <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
if (length(fn) != 1L) {
  stop("Run from the repo root: Rscript scripts/install_deps.R", call. = FALSE)
}
root <- normalizePath(file.path(dirname(fn), ".."))
setwd(root)

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

proj_lib <- .libPaths()[[1L]]
pkg_in_proj <- function(pkg) {
  file.exists(file.path(proj_lib, pkg))
}

pkgs <- c(
  "shiny", "bslib", "dplyr", "plotly", "httr2", "Matrix", "DT",
  "Seurat", "SeuratObject"
)

missing <- pkgs[!vapply(pkgs, pkg_in_proj, logical(1))]

if (length(missing)) {
  message("Project library: ", proj_lib)
  message("Installing (first run may take several minutes, especially Seurat): ",
    paste(missing, collapse = ", "))
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv", lib = proj_lib, dependencies = TRUE)
  }
  tryCatch(
    renv::install(missing, prompt = FALSE),
    error = function(e) {
      message("renv::install failed: ", conditionMessage(e))
      message("Falling back to install.packages(lib = project library) …")
      install.packages(missing, lib = proj_lib, dependencies = TRUE)
    }
  )
} else {
  message("All required packages are already in the project library.")
}

still <- pkgs[!vapply(pkgs, pkg_in_proj, logical(1))]
if (length(still)) {
  stop(
    "These packages are still missing from the project library:\n  ",
    paste(still, collapse = ", "),
    "\nTry: Rscript scripts/install_deps.R",
    call. = FALSE
  )
}

invisible(lapply(pkgs, function(p) {
  suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE))
}))

## Optional Bioconductor packages for reference mouse/human spleen defaults
bioc_pkgs <- c(
  "SingleCellExperiment", "SummarizedExperiment",
  "ExperimentHub", "TabulaMurisSenisData", "scRNAseq",
  "cellxgenedp", "zellkonverter"
)
if (!pkg_in_proj("BiocManager")) {
  message("Installing BiocManager (needed for reference spleen data packages)…")
  install.packages("BiocManager", lib = proj_lib, dependencies = NA)
}
bioc_missing <- bioc_pkgs[!vapply(bioc_pkgs, pkg_in_proj, logical(1))]
if (length(bioc_missing)) {
  message(
    "Installing Bioconductor packages for reference spleen data (first run can take a long time): ",
    paste(bioc_missing, collapse = ", ")
  )
  tryCatch(
    {
      loadNamespace("BiocManager", lib.loc = proj_lib)
      bioc_install <- getExportedValue("BiocManager", "install")
      bioc_install(bioc_missing, lib = proj_lib, update = FALSE, ask = FALSE)
    },
    error = function(e) {
      message("BiocManager::install failed: ", conditionMessage(e))
      message("Install later with: BiocManager::install(c('TabulaMurisSenisData','scRNAseq','cellxgenedp','zellkonverter'))")
    }
  )
}

if (requireNamespace("renv", quietly = TRUE)) {
  tryCatch(
    renv::snapshot(prompt = FALSE),
    error = function(e) message("(Optional) renv::snapshot skipped: ", conditionMessage(e))
  )
}

message("Done.")
