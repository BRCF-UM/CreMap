#' Load real reference spleen scRNA-seq for mouse and human defaults.
#'
#' **Preferred:** pre-built Seurat objects under **`reference/mouse_spleen_reference.rds`**
#' and **`reference/human_spleen_reference.rds`** (see `scripts/download_reference_data.R`).
#' Those load instantly with no ExperimentHub traffic.
#'
#' Otherwise: mouse from Tabula Muris Senis (`TabulaMurisSenisData`); human from He organ atlas
#' (`scRNAseq::HeOrganAtlasData`), with optional CELLxGENE if `CREMAP_ENABLE_CELLXGENE_HUMAN=1`.
#'
#' If packages or network fail, these return `NULL` and the app falls back to synthetic demos.

MAX_DEFAULT_CELLS <- 4500L

#' Whether live reference-atlas downloads are allowed
#'
#' Live ExperimentHub downloads are useful locally, but are a poor default on
#' ephemeral shinyapps.io workers. Set `CREMAP_ENABLE_REMOTE_REFERENCES=1` or
#' `0` to override. When unset, they are enabled locally and disabled when
#' `R_CONFIG_ACTIVE=shinyapps`.
cremap_remote_references_enabled <- function() {
  value <- trimws(tolower(Sys.getenv("CREMAP_ENABLE_REMOTE_REFERENCES", unset = "")))
  if (nzchar(value)) {
    return(value %in% c("1", "true", "yes", "on"))
  }
  !identical(tolower(Sys.getenv("R_CONFIG_ACTIVE", unset = "")), "shinyapps")
}

#' Allow hub cache creation and downloads without an R-console prompt (needed for Shiny).
#'
#' ExperimentHub / AnnotationHub normally ask interactively once to confirm the cache
#' directory. That question appears in the **terminal** running the app, not in the browser.
#' This function disables that prompt so reference data can load from the UI.
hub_allow_downloads_no_prompt <- function() {
  options(EXPERIMENT_HUB_ASK = FALSE, ANNOTATION_HUB_ASK = FALSE)
  if (requireNamespace("ExperimentHub", quietly = TRUE)) {
    suppressWarnings(
      try(ExperimentHub::setExperimentHubOption("ASK", FALSE), silent = TRUE)
    )
  }
  if (requireNamespace("AnnotationHub", quietly = TRUE)) {
    suppressWarnings(
      try(AnnotationHub::setAnnotationHubOption("ASK", FALSE), silent = TRUE)
    )
  }
  invisible(NULL)
}

#' Directory for shipped / pre-downloaded reference Seurat RDS files.
#'
#' Set **`CREMAP_REFERENCE_DIR`** to an absolute path to override the default
#' **`<project>/reference/`** (i.e. next to `app.R` when `getwd()` is the project root).
cremap_reference_dir <- function() {
  d <- Sys.getenv("CREMAP_REFERENCE_DIR", unset = "")
  if (nzchar(d)) {
    if (!dir.exists(d)) {
      warning("CREMAP_REFERENCE_DIR does not exist: ", d, call. = FALSE)
      return(file.path(normalizePath(getwd(), winslash = "/"), "reference"))
    }
    return(normalizePath(d, winslash = "/", mustWork = TRUE))
  }
  file.path(normalizePath(getwd(), winslash = "/"), "reference")
}

cremap_reference_paths <- function() {
  d <- cremap_reference_dir()
  list(
    mouse = file.path(d, "mouse_spleen_reference.rds"),
    human = file.path(d, "human_spleen_reference.rds")
  )
}

cremap_local_references_ready <- function() {
  p <- cremap_reference_paths()
  isTRUE(file.exists(p$mouse)) && isTRUE(file.exists(p$human))
}

#' Read a reference Seurat saved by `scripts/download_reference_data.R` (RNA + umap required).
read_reference_seurat_file <- function(path, species_label) {
  if (!is.character(path) || !nzchar(path) || !isTRUE(file.exists(path))) {
    return(NULL)
  }
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    return(NULL)
  }
  if (!"RNA" %in% names(obj@assays)) {
    warning("Local reference missing RNA assay: ", path, call. = FALSE)
    return(NULL)
  }
  if (!"umap" %in% names(obj@reductions)) {
    warning("Local reference missing umap reduction: ", path, call. = FALSE)
    return(NULL)
  }
  attr(obj, "cremap_default_source") <- "reference"
  attr(obj, "cremap_reference_origin") <- "local_rds"
  if (identical(species_label, "human")) {
    attr(obj, "cremap_human_reference") <- "local_rds"
  }
  obj
}

#' First non-empty metadata column from a preference list
detect_celltype_col <- function(cd) {
  pref <- c(
    "cell_ontology_class",
    "free_annotation",
    "reclustered.fine",
    "reclustered.broad",
    "cell_type",
    "celltype",
    "cell",
    "broad_cell",
    "cluster"
  )
  nm <- colnames(cd)
  hit <- pref[pref %in% nm]
  if (length(hit)) {
    return(hit[[1L]])
  }
  if ("leiden" %in% nm) {
    return("leiden")
  }
  NA_character_
}

#' Subset SingleCellExperiment columns
sce_sample_cells <- function(sce, max_cells, seed) {
  if (ncol(sce) <= max_cells) {
    return(sce)
  }
  set.seed(seed)
  sce[, sample.int(ncol(sce), max_cells)]
}

#' TabulaMurisSenisDroplet returns a named list of SCEs (or, rarely, one SCE).
tms_result_to_sce <- function(lst) {
  if (inherits(lst, "SingleCellExperiment")) {
    return(lst)
  }
  if (!is.list(lst) || !length(lst)) {
    return(NULL)
  }
  if (!is.null(lst$Spleen)) {
    return(lst$Spleen)
  }
  nm <- names(lst)
  if (length(nm)) {
    ii <- which(tolower(nm) == "spleen")
    if (length(ii)) {
      return(lst[[ii[[1L]]]])
    }
  }
  lst[[1L]]
}

#' Duplicate gene symbols break Seurat import; make.unique preserves order.
sce_dedup_rownames <- function(sce) {
  rn <- rownames(sce)
  if (!length(rn) || !anyDuplicated(rn)) {
    return(sce)
  }
  rownames(sce) <- make.unique(rn, sep = "_")
  sce
}

#' Fallback when [Seurat::as.Seurat] fails on some SCE objects.
sce_to_seurat_manual <- function(sce, species_tag, seed) {
  an <- SummarizedExperiment::assayNames(sce)
  pick <- function(cands) {
    hit <- cands[cands %in% an]
    if (length(hit)) {
      return(hit[[1L]])
    }
    an[[1L]]
  }
  count_name <- pick(c("counts", "count"))
  cm <- SummarizedExperiment::assay(sce, count_name)
  if (inherits(cm, "DelayedArray") || inherits(cm, "DelayedMatrix")) {
    cm <- as.matrix(cm)
  }
  if (anyDuplicated(rownames(cm))) {
    rownames(cm) <- make.unique(rownames(cm), sep = "_")
  }
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  ct_col <- detect_celltype_col(cd)
  if (is.na(ct_col)) {
    ct_vec <- rep("unknown", ncol(sce))
  } else {
    ct_vec <- as.character(cd[[ct_col]])
    ct_vec[!nzchar(ct_vec) | is.na(ct_vec)] <- "unknown"
  }
  obj <- Seurat::CreateSeuratObject(counts = cm, meta.data = cd)
  if ("logcounts" %in% an) {
    lc <- SummarizedExperiment::assay(sce, "logcounts")
    if (inherits(lc, "DelayedArray") || inherits(lc, "DelayedMatrix")) {
      lc <- as.matrix(lc)
    }
    if (anyDuplicated(rownames(lc))) {
      rownames(lc) <- make.unique(rownames(lc), sep = "_")
    }
    obj <- Seurat::SetAssayData(obj, layer = "data", new.data = lc)
  } else {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  }
  rd_names <- tryCatch(
    SingleCellExperiment::reducedDimNames(sce),
    error = function(e) character()
  )
  umap_name <- rd_names[tolower(rd_names) == "umap"]
  if (length(umap_name)) {
    emb <- as.matrix(SingleCellExperiment::reducedDim(sce, umap_name[[1L]]))
    colnames(emb) <- paste0("UMAP_", seq_len(ncol(emb)))
    rownames(emb) <- colnames(obj)
    obj[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = emb,
      key = "umap_",
      assay = Seurat::DefaultAssay(obj)
    )
  } else {
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = min(2000L, nrow(obj)), verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    npc <- min(30L, ncol(obj) - 1L, nrow(obj) - 1L)
    npc <- max(npc, 2L)
    obj <- Seurat::RunPCA(obj, npcs = npc, verbose = FALSE)
    obj <- Seurat::RunUMAP(obj, dims = 1:min(15L, npc), verbose = FALSE)
  }
  obj <- ensure_umap_named(obj)
  obj$cell_type <- ct_vec
  obj$species <- species_tag
  obj$tissue <- "spleen"
  Seurat::Idents(obj) <- obj$cell_type
  obj
}

#' Ensure Seurat has a reduction named `umap` (app requirement)
ensure_umap_named <- function(obj) {
  rns <- names(obj@reductions)
  if ("umap" %in% rns) {
    return(obj)
  }
  if ("UMAP" %in% rns) {
    obj[["umap"]] <- obj[["UMAP"]]
    return(obj)
  }
  for (nm in rns) {
    if (tolower(nm) == "umap") {
      obj[["umap"]] <- obj[[nm]]
      return(obj)
    }
  }
  obj
}

#' Hub-backed SCE objects often store counts in DelayedArray form; Seurat needs dense/sparse matrices.
sce_materialize_delayed_assays <- function(sce) {
  an <- SummarizedExperiment::assayNames(sce)
  if (!length(an)) {
    return(sce)
  }
  for (nm in an) {
    m <- SummarizedExperiment::assay(sce, nm, withDimnames = TRUE)
    if (inherits(m, "DelayedArray") || inherits(m, "DelayedMatrix")) {
      SummarizedExperiment::assay(sce, nm, withDimnames = TRUE) <- as.matrix(m)
    }
  }
  sce
}

#' Convert SCE subset to normalized Seurat with `cell_type` and `umap`
sce_to_seurat_default <- function(sce, species_tag, seed) {
  sce <- sce_dedup_rownames(sce)
  sce <- sce_materialize_delayed_assays(sce)
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  ct_col <- detect_celltype_col(cd)
  if (is.na(ct_col)) {
    ct_vec <- rep("unknown", ncol(sce))
  } else {
    ct_vec <- as.character(cd[[ct_col]])
    ct_vec[!nzchar(ct_vec) | is.na(ct_vec)] <- "unknown"
  }

  obj <- tryCatch(
    {
      o <- Seurat::as.Seurat(sce)
      ensure_umap_named(o)
    },
    error = function(e) {
      warning(
        "as.Seurat(SCE) failed (trying manual import): ",
        conditionMessage(e),
        call. = FALSE
      )
      tryCatch(
        sce_to_seurat_manual(sce, species_tag, seed),
        error = function(e2) {
          stop(
            "Could not convert SCE to Seurat (as.Seurat: ",
            conditionMessage(e), "; manual: ", conditionMessage(e2), ")",
            call. = FALSE
          )
        }
      )
    }
  )

  if (!"umap" %in% names(obj@reductions)) {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = min(2000L, nrow(obj)), verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    npc <- min(30L, ncol(obj) - 1L, nrow(obj) - 1L)
    npc <- max(npc, 2L)
    obj <- Seurat::RunPCA(obj, npcs = npc, verbose = FALSE)
    obj <- Seurat::RunUMAP(obj, dims = 1:min(15L, npc), verbose = FALSE)
    obj <- ensure_umap_named(obj)
  } else {
    nd <- suppressWarnings(try(Seurat::GetAssayData(obj, layer = "data"), silent = TRUE))
    if (inherits(nd, "try-error") || length(nd) == 0L) {
      obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    }
  }

  obj$cell_type <- ct_vec
  obj$species <- species_tag
  obj$tissue <- "spleen"
  Seurat::Idents(obj) <- obj$cell_type
  obj
}

#' Mouse spleen from Tabula Muris Senis (hub only); `NULL` on failure
load_mouse_spleen_from_hub <- function(max_cells = MAX_DEFAULT_CELLS, seed = 42L) {
  hub_allow_downloads_no_prompt()
  if (!requireNamespace("TabulaMurisSenisData", quietly = TRUE)) {
    return(NULL)
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    return(NULL)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    return(NULL)
  }
  lst <- tryCatch(
    TabulaMurisSenisData::TabulaMurisSenisDroplet(
      tissues = "Spleen",
      processedCounts = TRUE,
      reducedDims = TRUE,
      infoOnly = FALSE
    ),
    error = function(e) {
      warning("TabulaMurisSenisDroplet: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(lst) || !length(lst)) {
    return(NULL)
  }
  sce <- tms_result_to_sce(lst)
  if (is.null(sce)) {
    warning("TabulaMurisSenisDroplet: no SingleCellExperiment in result", call. = FALSE)
    return(NULL)
  }
  sce <- sce_sample_cells(sce, max_cells, seed)
  obj <- tryCatch(
    sce_to_seurat_default(sce, "mouse", seed),
    error = function(e) {
      warning("mouse SCE → Seurat: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (!is.null(obj)) {
    attr(obj, "cremap_default_source") <- "reference"
  }
  obj
}

#' Mouse spleen reference: local RDS first, then Tabula Muris Senis hub.
load_mouse_spleen_reference <- function(max_cells = MAX_DEFAULT_CELLS, seed = 42L) {
  p <- cremap_reference_paths()$mouse
  loc <- read_reference_seurat_file(p, "mouse")
  if (!is.null(loc)) {
    return(loc)
  }
  if (!cremap_remote_references_enabled()) {
    return(NULL)
  }
  load_mouse_spleen_from_hub(max_cells, seed)
}

#' Human spleen SCE from He organ atlas (several argument combinations).
load_human_he_organ_spleen_sce <- function() {
  if (!requireNamespace("scRNAseq", quietly = TRUE)) {
    return(NULL)
  }
  arg_lists <- list(
    list(tissue = "Spleen", ensembl = FALSE, location = FALSE, legacy = FALSE),
    list(tissue = "Spleen", ensembl = FALSE, location = FALSE, legacy = TRUE)
  )
  for (args in arg_lists) {
    sce <- tryCatch(
      do.call(scRNAseq::HeOrganAtlasData, args),
      error = function(e) NULL
    )
    if (inherits(sce, "SingleCellExperiment") && ncol(sce) > 0L) {
      return(sce)
    }
  }
  warning(
    "HeOrganAtlasData(spleen): all attempts failed; try scRNAseq update or CELLxGENE fallback.",
    call. = FALSE
  )
  NULL
}

#' Human spleen from CELLxGENE (Tabula Sapiens and other studies) via cellxgenedp.
#'
#' Picks a small human spleen `.h5ad` when possible to limit download size.
load_human_spleen_cellxgene_sce <- function() {
  if (!requireNamespace("cellxgenedp", quietly = TRUE)) {
    return(NULL)
  }
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    return(NULL)
  }
  hub_allow_downloads_no_prompt()
  tryCatch(
    {
      res <- NULL
      cgdb <- cellxgenedp::db()
      ds <- cellxgenedp::datasets(cgdb)
      has_org <- "organism" %in% colnames(ds) || "organism_ontology_term_id" %in% colnames(ds)
      if (NROW(ds) && has_org) {
        idx <- rep(TRUE, NROW(ds))
        if ("organism" %in% colnames(ds)) {
          oc <- ds[["organism"]]
          ft <- try(
            cellxgenedp::facets_filter(oc, "label", "Homo sapiens", exact = TRUE),
            silent = TRUE
          )
          if (!inherits(ft, "try-error") && length(ft) == NROW(ds)) {
            idx <- idx & ft
          } else {
            idx <- idx & grepl("Homo sapiens", as.character(oc), fixed = TRUE)
          }
        } else {
          ot <- ds[["organism_ontology_term_id"]]
          idx <- idx & grepl("9606", as.character(ot), fixed = TRUE)
        }
        tissue_col <- intersect(c("tissue", "tissue_type"), colnames(ds))
        if (length(tissue_col)) {
          tc <- ds[[tissue_col[[1L]]]]
          ft <- try(
            cellxgenedp::facets_filter(tc, "label", "spleen", exact = FALSE),
            silent = TRUE
          )
          if (!inherits(ft, "try-error") && length(ft) == NROW(ds)) {
            idx <- idx & ft
          } else {
            idx <- idx & grepl("spleen", as.character(tc), ignore.case = TRUE)
          }
        } else if ("title" %in% colnames(ds)) {
          idx <- idx & grepl("spleen", ds[["title"]], ignore.case = TRUE)
        } else if ("name" %in% colnames(ds)) {
          idx <- idx & grepl("spleen", ds[["name"]], ignore.case = TRUE)
        } else {
          idx <- rep(FALSE, NROW(ds))
        }
        ds_h <- ds[idx, , drop = FALSE]
        if (NROW(ds_h)) {
          if ("cell_count" %in% colnames(ds_h)) {
            cc <- suppressWarnings(as.numeric(ds_h[["cell_count"]]))
            o <- order(cc, na.last = TRUE)
            k <- min(8L, length(o))
            ds_h <- ds_h[o[seq_len(k)], , drop = FALSE]
          } else {
            ds_h <- ds_h[seq_len(min(5L, NROW(ds_h))), , drop = FALSE]
          }
          fl <- cellxgenedp::files(cgdb)
          fl <- fl[fl$dataset_id %in% ds_h$dataset_id & fl$filetype == "h5ad", , drop = FALSE]
          if (NROW(fl)) {
            if ("filesize" %in% colnames(fl)) {
              fs <- suppressWarnings(as.numeric(fl[["filesize"]]))
              fl <- fl[order(fs, na.last = TRUE), , drop = FALSE]
            }
            fp <- cellxgenedp::files_download(fl[1L, , drop = FALSE], dry.run = FALSE)
            if (length(fp) && nzchar(fp[[1L]]) && file.exists(fp[[1L]])) {
              res <- zellkonverter::readH5AD(fp[[1L]], use_h5 = TRUE)
            }
          }
        }
      }
      res
    },
    error = function(e) {
      warning("CELLxGENE human spleen: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}

#' Human spleen from hub only (He organ atlas, optional CELLxGENE); `NULL` on failure
#'
#' CELLxGENE is **off by default**. Set `CREMAP_ENABLE_CELLXGENE_HUMAN=1` to allow fallback.
load_human_spleen_from_hub <- function(max_cells = MAX_DEFAULT_CELLS, seed = 43L) {
  hub_allow_downloads_no_prompt()
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    return(NULL)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    return(NULL)
  }
  sce <- load_human_he_organ_spleen_sce()
  ref_tag <- "he_organ_atlas"
  if (is.null(sce) && identical(Sys.getenv("CREMAP_ENABLE_CELLXGENE_HUMAN", "0"), "1")) {
    sce <- load_human_spleen_cellxgene_sce()
    ref_tag <- "cellxgene"
  }
  if (is.null(sce)) {
    return(NULL)
  }
  sce <- sce_sample_cells(sce, max_cells, seed)
  obj <- tryCatch(
    sce_to_seurat_default(sce, "human", seed),
    error = function(e) {
      warning("human SCE → Seurat: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (!is.null(obj)) {
    attr(obj, "cremap_default_source") <- "reference"
    attr(obj, "cremap_human_reference") <- ref_tag
  }
  obj
}

#' Human spleen reference: local RDS first, then hub loaders.
load_human_spleen_reference <- function(max_cells = MAX_DEFAULT_CELLS, seed = 43L) {
  p <- cremap_reference_paths()$human
  loc <- read_reference_seurat_file(p, "human")
  if (!is.null(loc)) {
    return(loc)
  }
  if (!cremap_remote_references_enabled()) {
    return(NULL)
  }
  load_human_spleen_from_hub(max_cells, seed)
}
