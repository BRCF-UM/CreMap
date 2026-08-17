#!/usr/bin/env Rscript
## Deploy CreMAP to shinyapps.io without placing credentials in this repository.
## Run with a vanilla R session so the deployment client is independent of renv:
##   R --vanilla -f scripts/deploy_shinyapps.R

ca <- commandArgs(trailingOnly = FALSE)
fn <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
if (length(fn) != 1L) {
  stop("Run from the repo root: R --vanilla -f scripts/deploy_shinyapps.R", call. = FALSE)
}
root <- normalizePath(file.path(dirname(fn), ".."))
setwd(root)

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop(
    "The rsconnect deployment client is missing. Install it with:\n",
    "  R --vanilla -e 'install.packages(\"rsconnect\")'",
    call. = FALSE
  )
}

## Load the app's project library after rsconnect is loaded from the user's
## deployment-tool library. Strict manifest generation needs the installed
## package metadata as well as renv.lock.
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

app_name <- Sys.getenv("SHINYAPPS_APP_NAME", unset = "cremap")
account <- Sys.getenv("SHINYAPPS_ACCOUNT", unset = "")

files <- rsconnect::listDeploymentFiles(root)
file_sizes <- file.info(file.path(root, files))$size
bundle_mb <- sum(file_sizes, na.rm = TRUE) / 1024^2
message("Deployment bundle: ", length(files), " files, ", round(bundle_mb, 1), " MB")

reference_files <- c(
  "reference/mouse_spleen_reference.rds",
  "reference/human_spleen_reference.rds"
)
if (!all(reference_files %in% files)) {
  message(
    "Reference RDS files are not both present in the bundle. ",
    "The hosted app will use synthetic demo data."
  )
}

deploy_args <- list(
  appDir = root,
  appName = app_name,
  appTitle = "CreMAP",
  server = "shinyapps.io",
  appFiles = files,
  launch.browser = TRUE,
  forceUpdate = TRUE,
  dependencyResolution = "strict"
)
if (nzchar(account)) {
  deploy_args$account <- account
}

do.call(rsconnect::deployApp, deploy_args)
