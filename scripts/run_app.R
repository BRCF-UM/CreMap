#!/usr/bin/env Rscript
## Launch the CreMAP Shiny app (uses renv if present).

Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = "FALSE")

ca <- commandArgs(trailingOnly = FALSE)
fn <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
if (length(fn) != 1L) {
  stop("Run from the repo root: Rscript scripts/run_app.R", call. = FALSE)
}
root <- normalizePath(file.path(dirname(fn), ".."))
setwd(root)

## Bioconductor hubs ask in the R console to create a cache; Shiny has no prompt for that.
## Set before any Hub-using package loads (see ?ExperimentHub::getExperimentHubOption).
options(
  EXPERIMENT_HUB_ASK = FALSE,
  ANNOTATION_HUB_ASK = FALSE
)

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

proj_lib <- .libPaths()[[1L]]
pkg_in_proj <- function(pkg) {
  file.exists(file.path(proj_lib, pkg))
}

need <- c("shiny", "Seurat")
miss <- need[!vapply(need, pkg_in_proj, logical(1))]
if (length(miss)) {
  stop(
    "Missing from project library (", proj_lib, "): ", paste(miss, collapse = ", "),
    "\nYour global R may have these packages, but renv uses only the project library.",
    "\nRun: Rscript scripts/install_deps.R",
    call. = FALSE
  )
}

host <- Sys.getenv("SHINY_HOST", unset = "127.0.0.1")
port <- as.integer(Sys.getenv("SHINY_PORT", unset = "8787"))
launch <- tolower(Sys.getenv("SHINY_LAUNCH_BROWSER", unset = "true")) %in% c("1", "true", "yes")

message("CreMAP: http://", host, ":", port)
shiny::runApp(".", host = host, port = port, launch.browser = launch)
