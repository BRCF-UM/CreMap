#' Locate a Cell Ranger–style matrix directory under `root`.
#'
#' Looks for `matrix.mtx` or `matrix.mtx.gz` (recursively if needed), e.g.
#' `.../outs/filtered_feature_bc_matrix/` or the parent `outs/`.
#'
#' @param root Directory path (absolute recommended).
#' @return Normalized path to the folder containing the matrix, or `NA_character_`.
find_cellranger_matrix_dir <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  hits <- list.files(
    root,
    pattern = "^matrix\\.mtx(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(hits)) {
    return(normalizePath(dirname(hits[[1L]]), mustWork = FALSE))
  }
  for (nm in c("matrix.mtx", "matrix.mtx.gz")) {
    f <- file.path(root, nm)
    if (file.exists(f)) {
      return(normalizePath(root, mustWork = FALSE))
    }
  }
  NA_character_
}

#' Build a Seurat object from Cell Ranger matrix files and run a default QC + UMAP pipeline.
#'
#' Expects the usual trio under one directory (see [Seurat::Read10X()]):
#' `matrix.mtx` (or `.gz`), `barcodes.tsv` (or `.gz`), `features.tsv` or `genes.tsv` (or `.gz`).
#'
#' @param mtx_dir Directory that directly contains the matrix and barcode / feature files.
#' @param max_cells If set and counts exceed this, randomly subsample columns (for speed).
#' @param resolution Passed to [Seurat::FindClusters()].
#' @param npcs Max PCA dimensions (capped inside against matrix size).
#' @param verbose Passed to Seurat steps.
#' @return A Seurat object with `RNA`, `umap`, `seurat_clusters`, and `cell_type` (= clusters for convenience).
seurat_from_cellranger_mtx <- function(mtx_dir,
                                       max_cells = NULL,
                                       seed = 42L,
                                       resolution = 0.5,
                                       npcs = 50L,
                                       verbose = FALSE) {
  mtx_dir <- normalizePath(mtx_dir, mustWork = TRUE)
  counts <- Seurat::Read10X(mtx_dir)
  if (is.list(counts)) {
    counts <- counts[[1L]]
    warning("10X folder contained multiple matrices; using the first.", call. = FALSE)
  }
  if (!inherits(counts, "dgCMatrix")) {
    counts <- SeuratObject::as.sparse(counts)
  }

  if (is.numeric(max_cells) && max_cells > 0L && ncol(counts) > max_cells) {
    set.seed(seed)
    keep <- sample.int(ncol(counts), max_cells)
    counts <- counts[, keep, drop = FALSE]
  }

  obj <- Seurat::CreateSeuratObject(counts = counts, project = "CellRanger")
  obj <- Seurat::NormalizeData(obj, verbose = verbose)
  obj <- Seurat::FindVariableFeatures(
    obj,
    nfeatures = min(3000L, nrow(obj)),
    verbose = verbose
  )
  obj <- Seurat::ScaleData(obj, verbose = verbose)
  nc <- ncol(obj)
  ng <- nrow(obj)
  ndim <- min(
    max(2L, as.integer(npcs)),
    max(2L, nc - 1L),
    max(2L, ng - 1L)
  )
  obj <- Seurat::RunPCA(obj, npcs = ndim, verbose = verbose)
  n_use <- min(30L, ndim, ncol(obj[["pca"]]))
  n_use <- max(n_use, 2L)
  obj <- Seurat::RunUMAP(obj, dims = seq_len(n_use), verbose = verbose)
  obj <- Seurat::FindNeighbors(obj, dims = seq_len(n_use), verbose = verbose)
  obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = verbose)
  obj$cell_type <- obj$seurat_clusters
  Seurat::Idents(obj) <- obj$seurat_clusters
  obj
}

#' Unzip archive, find matrix directory, build Seurat object (then removes temp dir).
seurat_from_cellranger_zip <- function(zip_path, ...) {
  zip_info <- utils::unzip(zip_path, list = TRUE)
  if (!NROW(zip_info)) {
    stop("ZIP archive is empty.", call. = FALSE)
  }

  archive_names <- gsub("\\\\", "/", as.character(zip_info$Name))
  unsafe_name <- grepl("^/|^[A-Za-z]:/", archive_names) |
    vapply(strsplit(archive_names, "/", fixed = TRUE), function(parts) ".." %in% parts, logical(1))
  if (any(unsafe_name)) {
    stop("ZIP contains an unsafe absolute or parent-directory path.", call. = FALSE)
  }

  max_uncompressed_mb <- suppressWarnings(as.numeric(
    Sys.getenv("CREMAP_MAX_UNCOMPRESSED_MB", unset = "500")
  ))
  if (!is.finite(max_uncompressed_mb) || max_uncompressed_mb <= 0) {
    max_uncompressed_mb <- 500
  }
  expanded_bytes <- sum(as.numeric(zip_info$Length), na.rm = TRUE)
  if (expanded_bytes > max_uncompressed_mb * 1024^2) {
    stop(
      "ZIP expands to ", round(expanded_bytes / 1024^2),
      " MB; server limit is ", format(max_uncompressed_mb), " MB.",
      call. = FALSE
    )
  }

  td <- tempfile("cr_mtx_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::unzip(zip_path, exdir = td)
  mtx_dir <- find_cellranger_matrix_dir(td)
  if (is.na(mtx_dir)) {
    stop(
      "ZIP did not contain matrix.mtx / matrix.mtx.gz with barcodes and features.",
      call. = FALSE
    )
  }
  seurat_from_cellranger_mtx(mtx_dir, ...)
}
