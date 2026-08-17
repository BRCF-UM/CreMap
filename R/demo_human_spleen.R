#' Human spleen–like demo (synthetic scRNA-seq) for side-by-side comparison with mouse.
#'
#' Uses HGNC-style symbols and the same **granular** splenic cell-type scheme as the
#' mouse demo (B / T / NK / myeloid / DC / macrophage / stromal subsets).

HUMAN_SPLEEN_MARKERS <- data.frame(
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
    "MS4A1", "CD19", "BANK1",
    "MS4A1", "CR2", "CD1D",
    "JCHAIN", "MZB1", "IRF4",
    "CD3D", "CD4", "SELL",
    "CD3D", "CD4", "CD44",
    "CD3D", "CD8A", "SELL",
    "CD3D", "CD8A", "GZMB",
    "CD3D", "CD4", "FOXP3",
    "NCR1", "KLRD1", "KLRK1",
    "LY6C2", "CCR2", "CSF1R",
    "S100A8", "S100A9", "LY6G",
    "CLEC9A", "XCR1", "ITGAX",
    "ITGAX", "SIRPA", "FCGR1A",
    "IL3RA", "IRF7", "TCF4",
    "ADGRE1", "MARCO", "MRC1",
    "COL1A1", "DCN", "PDGFRA"
  ),
  stringsAsFactors = FALSE
)

#' @param seed RNG seed (use a different seed than mouse demo for distinct layout).
#' @return Normalized Seurat object with RNA, UMAP, `cell_type`.
build_demo_human_spleen_seurat <- function(seed = 43L) {
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
  extra_genes <- sprintf("LINC%05d", seq_len(extra_n))
  genes <- unique(c(
    HUMAN_SPLEEN_MARKERS$gene,
    "MALAT1", "ACTB", "GAPDH", "RPS29",
    extra_genes
  ))
  G <- length(genes)
  gi <- stats::setNames(seq_len(G), genes)

  counts <- Matrix::Matrix(0, nrow = G, ncol = n, sparse = TRUE)
  colnames(counts) <- paste0("HUMAN_", seq_len(n))
  rownames(counts) <- genes

  lambda_bg <- stats::rgamma(G, shape = 2, rate = 2)
  lambda_bg[gi[c("MALAT1", "ACTB", "GAPDH", "RPS29")]] <-
    lambda_bg[gi[c("MALAT1", "ACTB", "GAPDH", "RPS29")]] * 4

  for (j in seq_len(n)) {
    ct <- types[j]
    lam <- lambda_bg
    mk <- HUMAN_SPLEEN_MARKERS[HUMAN_SPLEEN_MARKERS$cell_type == ct, "gene", drop = TRUE]
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
    species = "human",
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta, project = "HumanSpleenDemo")
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
