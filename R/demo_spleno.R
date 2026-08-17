#' Build a small mouse splenocyte-like Seurat object for local demos.
#'
#' Uses realistic MGI gene symbols and **granular** splenic populations (B subsets,
#' T subsets, NK, monocyte, neutrophil, DC subsets, macrophages, stromal). Replace
#' with your own `Seurat` object (RDS) for publication-grade analysis.

SPLENO_MARKERS <- data.frame(
  cell_type = c(
    rep("B_follicular", 3L),
    rep("B_marginal_zone", 3L),
    rep("Plasma_cell", 3L),
    rep("T_CD4_naive", 3L),
    rep("T_CD4_memory", 3L),
    rep("T_CD8_naive", 3L),
    rep("T_CD8_effector", 3L),
    rep("Treg", 3L),
    rep("NK_cell", 3L),
    rep("Monocyte", 3L),
    rep("Neutrophil", 3L),
    rep("cDC1", 3L),
    rep("cDC2", 3L),
    rep("pDC", 3L),
    rep("Macrophage_RP", 3L),
    rep("Stromal_FRC", 3L)
  ),
  gene = c(
    "Ms4a1", "Cd19", "Bank1",
    "Ms4a1", "Cr2", "Cd1d1",
    "Jchain", "Mzb1", "Irf4",
    "Cd3d", "Cd4", "Sell",
    "Cd3d", "Cd4", "Cd44",
    "Cd3d", "Cd8a", "Sell",
    "Cd3d", "Cd8a", "Gzmb",
    "Cd3d", "Cd4", "Foxp3",
    "Ncr1", "Klrb1c", "Klrk1",
    "Ly6c2", "Ccr2", "Csf1r",
    "S100a8", "S100a9", "Ly6g",
    "Clec9a", "Xcr1", "Itgax",
    "Itgax", "Sirpa", "Fcgr1",
    "Siglech", "Irf7", "Tcf4",
    "Adgre1", "Marco", "Mrc1",
    "Col1a1", "Dcn", "Pdgfra"
  ),
  stringsAsFactors = FALSE
)

#' @param seed Integer RNG seed for reproducible demo data.
#' @return A normalized Seurat object with RNA assay, UMAP, and `cell_type`.
build_demo_spleno_seurat <- function(seed = 42L) {
  set.seed(seed)

  n_by <- c(
    B_follicular = 140L,
    B_marginal_zone = 100L,
    Plasma_cell = 80L,
    T_CD4_naive = 160L,
    T_CD4_memory = 180L,
    T_CD8_naive = 120L,
    T_CD8_effector = 140L,
    Treg = 90L,
    NK_cell = 130L,
    Monocyte = 120L,
    Neutrophil = 150L,
    cDC1 = 80L,
    cDC2 = 160L,
    pDC = 70L,
    Macrophage_RP = 140L,
    Stromal_FRC = 60L
  )
  types <- rep(names(n_by), n_by)
  n <- length(types)

  extra_n <- 380L
  extra_genes <- sprintf("Gm%05d", seq_len(extra_n))
  genes <- unique(c(SPLENO_MARKERS$gene, "Malat1", "Actb", "Gapdh", "Rps29", extra_genes))
  G <- length(genes)
  gi <- stats::setNames(seq_len(G), genes)

  counts <- Matrix::Matrix(0, nrow = G, ncol = n, sparse = TRUE)
  colnames(counts) <- paste0("cell_", seq_len(n))
  rownames(counts) <- genes

  lambda_bg <- stats::rgamma(G, shape = 2, rate = 2)
  lambda_bg[gi[c("Malat1", "Actb", "Gapdh", "Rps29")]] <-
    lambda_bg[gi[c("Malat1", "Actb", "Gapdh", "Rps29")]] * 4

  for (j in seq_len(n)) {
    ct <- types[j]
    lam <- lambda_bg
    mk <- SPLENO_MARKERS[SPLENO_MARKERS$cell_type == ct, "gene", drop = TRUE]
    for (g in mk) {
      lam[gi[[g]]] <- lam[gi[[g]]] + stats::rlnorm(1, meanlog = 2.2, sdlog = 0.25)
    }
    counts[, j] <- Matrix::sparseMatrix(
      i = seq_len(G),
      j = rep(1L, G),
      x = stats::rpois(G, pmax(lam, 1e-6)),
      dims = c(G, 1L)
    )
  }

  meta <- data.frame(
    cell_type = types,
    tissue = "spleen",
    species = "mouse",
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta, project = "MouseSpleenDemo")
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = min(2000L, G), verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = min(30L, G - 1L), verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, dims = 1:min(20L, ncol(obj[["pca"]])),
                         verbose = FALSE)
  Seurat::Idents(obj) <- obj$cell_type
  attr(obj, "cremap_default_source") <- "synthetic"
  obj
}

#' Score cells by simple marker panel (mean z-score of panel genes per cell type)
suggest_cell_types_from_markers <- function(obj, panel = SPLENO_MARKERS, assay = "RNA") {
  genes_use <- intersect(unique(panel$gene), rownames(obj))
  panel <- panel[panel$gene %in% genes_use, , drop = FALSE]
  if (!nrow(panel)) {
    return(rep(NA_character_, ncol(obj)))
  }

  mat <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = "data")
  )
  mat <- as.matrix(mat[genes_use, , drop = FALSE])
  z <- t(scale(t(mat)))

  types <- unique(panel$cell_type)
  scores <- sapply(types, function(ct) {
    g <- panel$gene[panel$cell_type == ct]
    if (!length(g)) {
      return(rep(0, ncol(obj)))
    }
    colMeans(z[g, , drop = FALSE], na.rm = TRUE)
  })
  types[max.col(scores, ties.method = "first")]
}
