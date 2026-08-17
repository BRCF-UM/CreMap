## CreMAP — mouse splenocyte scRNA-seq explorer with MouseMine (MGI) Cre lookup

## File uploads are 5 MB by default in Shiny, which is too small for Seurat/10x
## data. Keep the public default conservative; operators can override it.
cremap_max_upload_mb <- suppressWarnings(as.numeric(
  Sys.getenv("CREMAP_MAX_UPLOAD_MB", unset = "100")
))
if (!is.finite(cremap_max_upload_mb) || cremap_max_upload_mb <= 0) {
  cremap_max_upload_mb <- 100
}
options(shiny.maxRequestSize = cremap_max_upload_mb * 1024^2)

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(dplyr)
  library(plotly)
  library(httr2)
  library(Seurat)
  library(Matrix)
  library(DT)
})

## Hub packages (TabulaMurisSenisData, scRNAseq atlases) read this when loaded: no
## interactive cache prompt in the R console while Shiny runs in the browser.
options(EXPERIMENT_HUB_ASK = FALSE, ANNOTATION_HUB_ASK = FALSE)

source("R/mousemine.R", local = TRUE)
source("R/demo_spleno.R", local = TRUE)
source("R/demo_human_spleen.R", local = TRUE)
source("R/read10x.R", local = TRUE)
source("R/real_spleen_data.R", local = TRUE)

default_meta_col <- "cell_type"

`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
}

default_data_kind <- function(obj) {
  a <- attr(obj, "cremap_default_source", exact = TRUE)
  if (identical(a, "reference")) {
    if (identical(attr(obj, "cremap_reference_origin", exact = TRUE), "local_rds")) {
      return("reference atlas (local RDS)")
    }
    hr <- attr(obj, "cremap_human_reference", exact = TRUE)
    if (identical(hr, "cellxgene")) {
      return("reference atlas (human: CELLxGENE)")
    }
    if (identical(hr, "he_organ_atlas")) {
      return("reference atlas (human: He et al. 2020)")
    }
    return("reference atlas")
  }
  if (identical(a, "synthetic")) {
    "synthetic demo"
  } else {
    "unknown"
  }
}

umap_df <- function(obj, color_col = NULL, expr_gene = NULL) {
  emb <- Seurat::Embeddings(obj, reduction = "umap")
  vars <- unique(c(
    if (!is.null(color_col) && nzchar(color_col)) color_col,
    if (!is.null(expr_gene) && nzchar(expr_gene)) expr_gene
  ))
  vars <- vars[vars %in% c(colnames(obj@meta.data), rownames(obj))]
  if (length(vars)) {
    df <- Seurat::FetchData(obj, vars = vars)
  } else {
    df <- data.frame(row.names = colnames(obj))
  }
  if (!is.null(expr_gene) && nzchar(expr_gene) && !expr_gene %in% colnames(df)) {
    if (expr_gene %in% rownames(obj)) {
      df[[expr_gene]] <- as.numeric(Seurat::FetchData(obj, vars = expr_gene)[[1]])
    } else {
      df[[expr_gene]] <- NA_real_
    }
  }
  df$cell <- colnames(obj)
  df$UMAP_1 <- emb[, 1]
  df$UMAP_2 <- if (ncol(emb) >= 2) emb[, 2] else rep(0, nrow(emb))
  df
}

validate_user_seurat <- function(obj) {
  if (!inherits(obj, "Seurat")) {
    return("File must contain a Seurat object.")
  }
  if (!"umap" %in% names(obj@reductions)) {
    return("Seurat object must include a UMAP reduction named 'umap'.")
  }
  if (!"RNA" %in% names(obj@assays)) {
    return("Seurat object must include an RNA assay.")
  }
  NA_character_
}

plotly_no_data <- function(title) {
  plotly::plot_ly() |>
    plotly::layout(
      title = list(text = title, font = list(size = 13)),
      xaxis = list(visible = FALSE, range = c(0, 1)),
      yaxis = list(visible = FALSE, range = c(0, 1)),
      margin = list(t = 40)
    )
}

## Height (px) for horizontal Cre driver bar chart: room per cell type + title/axes
cre_bar_plot_height_px <- function(n_types) {
  n <- max(1L, suppressWarnings(as.integer(n_types)))
  max(240L, min(960L, 72L + 22L * n))
}

## Stable y-axis order for Cre bar charts (alphabetical by cell type, not by mean)
cre_bar_df_fixed_levels <- function(df) {
  if (!nrow(df)) {
    return(df)
  }
  lev <- sort(unique(as.character(df$ct)))
  df$ct <- factor(df$ct, levels = lev)
  df
}

normalize_gb_query_id <- function(x) {
  x <- trimws(as.character(x))
  if (!length(x) || is.na(x) || !nzchar(x)) {
    return(NA_character_)
  }
  if (grepl("^MGI:\\d+", x, ignore.case = TRUE)) {
    return(sub("^mgi:", "MGI:", x, ignore.case = TRUE))
  }
  resolve_allele_mgi_id(x)
}

plot_gb_locus <- function(loci, overlap, symbol = "") {
  if (is.null(loci) || !nrow(loci)) {
    return(plotly_no_data("No integration coordinates in MGI for this allele"))
  }
  loc <- loci[1, , drop = FALSE]
  chr <- loc$chromosome[[1]]
  a0 <- loc$start[[1]]
  b0 <- loc$end[[1]]
  if (is.na(a0) || is.na(b0)) {
    return(plotly_no_data("Invalid genomic coordinates"))
  }
  pad <- max(5000L, as.integer(0.15 * (b0 - a0 + 1L)))
  xmin <- max(1L, a0 - pad)
  xmax <- b0 + pad
  title <- if (nzchar(symbol)) paste0(symbol, " — chr", chr) else paste0("chr", chr)
  p <- plotly::plot_ly() |>
    plotly::layout(
      title = list(text = title, font = list(size = 14)),
      xaxis = list(
        title = paste0("Position (", loc$assembly[[1]], ")"),
        range = c(xmin, xmax),
        zeroline = FALSE
      ),
      yaxis = list(
        title = "",
        range = c(0.2, 2.4),
        showticklabels = FALSE,
        zeroline = FALSE
      ),
      showlegend = FALSE,
      margin = list(t = 48, b = 40)
    )
  p <- plotly::add_trace(
    p,
    x = c(xmin, xmax),
    y = c(1, 1),
    type = "scatter",
    mode = "lines",
    line = list(color = "#bdc3c7", width = 6),
    hoverinfo = "skip",
    showlegend = FALSE
  )
  p <- plotly::add_trace(
    p,
    x = c(a0, b0),
    y = c(1, 1),
    type = "scatter",
    mode = "lines",
    line = list(color = "#e74c3c", width = 14),
    name = "Insertion",
    hovertemplate = paste0(
      "Insertion<br>chr", chr, ":",
      format(a0, big.mark = ","), "-",
      format(b0, big.mark = ","),
      "<extra></extra>"
    ),
    showlegend = FALSE
  )
  if (!is.null(overlap) && nrow(overlap)) {
    syms <- unique(overlap$symbol)
    syms <- syms[!is.na(syms) & nzchar(syms)]
    n <- min(length(syms), 25L)
    if (n > 0L) {
      syms <- head(syms[order(syms)], n)
      yvals <- seq(1.5, 1.5 + 0.08 * (n - 1L), length.out = n)
      p <- plotly::add_trace(
        p,
        x = rep((a0 + b0) / 2, n),
        y = yvals,
        type = "scatter",
        mode = "markers+text",
        text = syms,
        textposition = "middle right",
        marker = list(size = 8, color = "#2980b9"),
        hovertemplate = "%{text}<extra>overlapping feature</extra>",
        showlegend = FALSE
      )
    }
  }
  p
}

ui <- page_navbar(
  title = "CreMAP",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2c3e50"),
  header = tags$div(
    class = "border-bottom px-3 py-2",
    tags$div(
      class = "d-flex flex-wrap align-items-center justify-content-between gap-3",
      tags$div(
        class = "text-muted small",
        "Mouse & human spleen single-cell RNA-seq — pick a species below, then explore expression, ",
        tags$a(href = "https://www.mousemine.org/mousemine", "MouseMine"),
        " (MGI) Cre drivers, integration-locus browser, cell types, and differential expression."
      ),
      tags$div(
        class = "d-flex align-items-center gap-2",
        tags$span(class = "small fw-semibold text-nowrap", "Species"),
        radioButtons(
          "view_species",
          label = NULL,
          choices = c("Mouse" = "mouse", "Human" = "human"),
          selected = "mouse",
          inline = TRUE
        )
      )
    )
  ),
  nav_panel(
    title = "Overview & data",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        h5("Data source"),
        radioButtons(
          "data_mode",
          NULL,
          choices = c(
            "Built-in default (reference spleen; synthetic fallback)" = "demo",
            "Seurat RDS (precomputed object)" = "rds",
            "Cell Ranger 10x matrix (MTX + barcodes + features)" = "cellranger"
          ),
          selected = "demo"
        ),
        helpText(
          "If ", tags$code("reference/mouse_spleen_reference.rds"), " and ",
          tags$code("reference/human_spleen_reference.rds"), " exist (run ",
          tags$code("Rscript scripts/download_reference_data.R"), " once), the app loads them ",
          "directly—no ExperimentHub at runtime. Otherwise demo mode shows synthetic data. ",
          if (cremap_remote_references_enabled()) {
            "Remote references are enabled, so Tabula Muris Senis and the He atlas load in the background."
          } else {
            "Remote atlas loading is disabled on this host."
          }
        ),
        conditionalPanel(
          condition = "input.data_mode == 'rds'",
          radioButtons(
            "rds_species",
            "Assign this object to",
            choices = c("Mouse" = "mouse", "Human" = "human"),
            inline = TRUE
          ),
          fileInput("rds_file", "Seurat RDS", accept = c(".rds", ".RDS")),
          textInput("rds_meta_col", "Metadata column for identities", value = default_meta_col),
          helpText("RDS must contain RNA assay and a UMAP named ", tags$code("umap"), ".")
        ),
        conditionalPanel(
          condition = "input.data_mode == 'cellranger'",
          radioButtons(
            "cr_species",
            "Species for this matrix",
            choices = c("Mouse" = "mouse", "Human" = "human"),
            inline = TRUE
          ),
          radioButtons(
            "cr_mode",
            "How to provide data",
            choices = c(
              "Folder path on this computer" = "path",
              "Upload a ZIP of the matrix folder" = "zip"
            ),
            selected = "path"
          ),
          conditionalPanel(
            condition = "input.data_mode == 'cellranger' && input.cr_mode == 'path'",
            textInput(
              "cr_path",
              "Path to matrix folder (or any parent)",
              value = "",
              placeholder = "/path/to/filtered_feature_bc_matrix"
            )
          ),
          conditionalPanel(
            condition = "input.data_mode == 'cellranger' && input.cr_mode == 'zip'",
            fileInput(
              "cr_zip",
              "ZIP file",
              accept = c(".zip", "application/zip"),
              buttonLabel = "Browse…"
            ),
            helpText(
              "ZIP should expand to a folder that contains ",
              tags$code("matrix.mtx"), " or ", tags$code("matrix.mtx.gz"),
              ", barcodes, and features (Cell Ranger layout). ",
              "Maximum upload: ", format(cremap_max_upload_mb), " MB."
            )
          ),
          hr(),
          numericInput(
            "cr_max_cells",
            "Max cells (0 = use all; subsample for speed)",
            value = 20000L,
            min = 0L,
            max = 1e7L,
            step = 1000L
          ),
          fluidRow(
            column(
              6,
              numericInput("cr_resolution", "Cluster resolution", value = 0.5, min = 0.1, max = 2, step = 0.1)
            ),
            column(
              6,
              numericInput("cr_npc", "PCA dimensions (max)", value = 50L, min = 5L, max = 100L, step = 5L)
            )
          ),
          helpText(
            "Expected files in the matrix directory (names may end in ",
            tags$code(".gz"), "): ",
            tags$code("matrix.mtx"), ", ",
            tags$code("barcodes.tsv"), ", ",
            tags$code("features.tsv"), " (or ", tags$code("genes.tsv"), " for older Cell Ranger). ",
            "The app runs normalization, PCA, UMAP, and clustering; ",
            tags$code("cell_type"), " is set to cluster IDs until you annotate."
          ),
          actionButton("load_cellranger", "Load Cell Ranger matrix", class = "btn-primary w-100")
        ),
        hr(),
        actionButton("reload_data", "Reload default spleen data", class = "btn-outline-primary w-100"),
        hr(),
        h5("Display"),
        selectInput(
          "color_meta",
          "Color UMAP by metadata column",
          choices = default_meta_col,
          selected = default_meta_col
        ),
        sliderInput(
          "pt_size",
          "UMAP point size",
          min = 3,
          max = 5,
          value = 4,
          step = 0.25
        ),
        sliderInput(
          "pt_alpha",
          "UMAP point opacity",
          min = 0.1,
          max = 1,
          value = 0.85,
          step = 0.05
        )
      ),
      card(
        card_header("Dataset summary"),
        verbatimTextOutput("data_summary")
      )
    )
  ),
  nav_panel(
    title = "Integration locus",
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        textInput(
          "gb_query",
          "MGI allele ID or symbol",
          value = "",
          placeholder = "MGI:2176173 or Tg(Nes-cre)1Kln"
        ),
        numericInput(
          "gb_ucsc_pad",
          "UCSC window padding (bp each side)",
          value = 75000L,
          min = 5000L,
          max = 500000L,
          step = 5000L
        ),
        actionButton("gb_fetch", "Load from MouseMine", class = "btn-primary w-100"),
        checkboxInput(
          "gb_auto_load",
          "Auto-load when a Cre allele is selected (MGI Cre drivers tab)",
          value = TRUE
        ),
        hr(),
        helpText(
          "Many transgenic Cre lines have published ",
          tags$strong("confounded integration sites"),
          " (insertion disrupts endogenous loci). MGI reports genomic intervals for the ",
          "transgene feature and overlapping genes at that locus."
        ),
        helpText(
          "Enter an ",
          tags$code("MGI:#####"),
          " ID or allele symbol (partial symbol match supported). Requires network access."
        ),
        tags$div(
          class = "small border rounded px-2 py-2 bg-light",
          uiOutput("gb_quick_links")
        )
      ),
      layout_columns(
        col_widths = c(12, 12, 6, 6),
        gap = "0.75rem",
        card(
          card_header("Locus summary"),
          uiOutput("gb_summary")
        ),
        card(
          card_header("UCSC Genome Browser"),
          uiOutput("gb_ucsc_frame"),
          class = "mb-0"
        ),
        card(
          card_header("Integration interval(s)"),
          DTOutput("gb_loci_table")
        ),
        card(
          card_header("Overlapping genome features (MGI)"),
          DTOutput("gb_overlap_table")
        ),
        card(
          full_screen = TRUE,
          card_header("Locus map (insertion vs overlapping genes)"),
          plotlyOutput("gb_locus_plot", height = "320px")
        )
      )
    )
  ),
  nav_panel(
    title = "Gene expression",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        selectizeInput(
          "gene_mouse",
          "Mouse gene (RNA rownames)",
          choices = NULL,
          options = list(placeholder = "e.g. Cd19", maxOptions = 40)
        ),
        selectizeInput(
          "gene_human",
          "Human gene (RNA rownames)",
          choices = NULL,
          options = list(placeholder = "e.g. CD19", maxOptions = 40)
        ),
        checkboxInput("sync_human_gene", "Auto-pick human gene from mouse symbol", value = TRUE),
        helpText(
          "Plots follow the species selected in the header. UMAP point size and opacity: ",
          strong("Overview & data"), " → Display."
        ),
        hr(),
        checkboxInput("gene_log", "Log1p color scale", value = TRUE)
      ),
      conditionalPanel(
        condition = "input.view_species == 'mouse'",
        card(
          full_screen = TRUE,
          card_header("Mouse — UMAP (gene expression)"),
          plotlyOutput("plot_gene_umap_mouse", height = "480px")
        )
      ),
      conditionalPanel(
        condition = "input.view_species == 'human'",
        card(
          full_screen = TRUE,
          card_header("Human — UMAP (gene expression)"),
          plotlyOutput("plot_gene_umap_human", height = "480px")
        )
      ),
      conditionalPanel(
        condition = "input.view_species == 'mouse'",
        card(
          card_header("Mouse — mean expression by cell type"),
          plotlyOutput("plot_gene_bar_mouse", height = "300px")
        )
      ),
      conditionalPanel(
        condition = "input.view_species == 'human'",
        card(
          card_header("Human — mean expression by cell type"),
          plotlyOutput("plot_gene_bar_human", height = "300px")
        )
      )
    )
  ),
  nav_panel(
    title = "MGI Cre drivers",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        textInput("cre_search", "Filter allele symbol (contains)", value = ""),
        numericInput("cre_limit", "Max alleles to fetch", value = 400, min = 50, max = 2000, step = 50),
        actionButton("cre_fetch", "Query MouseMine", class = "btn-primary w-100"),
        hr(),
        helpText(
          "Queries recombinase alleles from ",
          tags$a(href = "https://www.informatics.jax.org/", "MGI"),
          " via ",
          tags$a(href = "https://www.mousemine.org/mousemine", "MouseMine"),
          ". Driver gene symbols are joined when annotated in MGI."
        ),
        helpText(
          tags$strong("Driver bar chart:"), " uses ",
          tags$code("id_col"),
          " from the ",
          strong("Cell types"),
          " tab; cell types are alphabetical on the y-axis."
        ),
        helpText(
          "After selecting a row, open ",
          strong("Integration locus"),
          " to inspect insertion coordinates, overlapping genes, and a UCSC genome-browser view."
        ),
        hr(),
        tags$p(class = "small fw-semibold mb-1", "Driver vs datasets"),
        tags$div(
          class = "small border rounded px-2 py-1 bg-light",
          style = "max-height: 12rem; overflow-y: auto;",
          verbatimTextOutput("cre_driver_status")
        )
      ),
      layout_columns(
        col_widths = c(5, 7),
        gap = "0.75rem",
        card(
          card_header(class = "py-2 fs-6", "Cre alleles (select a row)"),
          DTOutput("cre_table")
        ),
        tags$div(
          class = "cre-mgi-bar-col",
          style = "min-width:0;min-height:0;display:flex;flex-direction:column;max-height:65vh;",
          conditionalPanel(
            condition = "input.view_species == 'mouse'",
            tags$div(
              class = "card border shadow-none",
              style = "min-height:0;display:flex;flex-direction:column;flex:1 1 auto;overflow:hidden;",
              card_header(class = "py-2 fs-6", "Mean driver expression by cell type"),
              tags$div(
                style = "min-height:0;flex:1 1 auto;overflow:auto;width:100%;",
                uiOutput("ui_cre_bar_mouse")
              )
            )
          ),
          conditionalPanel(
            condition = "input.view_species == 'human'",
            tags$div(
              class = "card border shadow-none",
              style = "min-height:0;display:flex;flex-direction:column;flex:1 1 auto;overflow:hidden;",
              card_header(class = "py-2 fs-6", "Mean driver expression by cell type"),
              tags$div(
                style = "min-height:0;flex:1 1 auto;overflow:auto;width:100%;",
                uiOutput("ui_cre_bar_human")
              )
            )
          )
        )
      ),
      conditionalPanel(
        condition = "input.view_species == 'mouse'",
        card(
          full_screen = TRUE,
          card_header(class = "py-2", "Mouse — UMAP (driver expression)"),
          plotlyOutput("plot_cre_umap_mouse", height = "78vh")
        )
      ),
      conditionalPanel(
        condition = "input.view_species == 'human'",
        card(
          full_screen = TRUE,
          card_header(class = "py-2", "Human — UMAP (driver expression)"),
          plotlyOutput("plot_cre_umap_human", height = "78vh")
        )
      )
    )
  ),
  nav_panel(
    title = "Cell types",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        selectInput("id_col", "Identity column for analysis", choices = default_meta_col),
        hr(),
        h5("Marker-based labeling"),
        helpText(
          "Scores mouse and human (when loaded) with splenocyte marker panels; writes ",
          tags$code("cell_type_suggested"), " per species."
        ),
        helpText(
          "With reference data, identities come from the atlas metadata (cell ontology / author labels). ",
          "The synthetic fallback uses the same granular scheme for mouse and human, e.g. ",
          tags$code("B_follicular"), ", ", tags$code("B_marginal_zone"), ", ",
          tags$code("Plasma_cell"), ", ", tags$code("T_CD4_naive"), ", ",
          tags$code("T_CD8_effector"), ", ", tags$code("Treg"), ", ",
          tags$code("cDC1"), ", ", tags$code("pDC"), ", ",
          tags$code("Neutrophil"), ", ", tags$code("Macrophage_RP"), ", ",
          tags$code("Stromal_FRC"), ", and more."
        ),
        actionButton("apply_suggest", "Suggest labels from markers (mouse & human)", class = "btn-warning w-100"),
        hr(),
        h5("Manual override (mouse only)"),
        helpText("Edits apply to mouse ", tags$code("cell_type_manual"), " only."),
        selectizeInput("cell_pick", "Pick cell (barcode)", choices = NULL),
        textInput("manual_label", "New label", value = ""),
        actionButton("apply_manual", "Set label for selected cell", class = "btn-secondary w-100")
      ),
      conditionalPanel(
        condition = "input.view_species == 'mouse'",
        card(
          full_screen = TRUE,
          card_header("Mouse — UMAP by identity"),
          plotlyOutput("plot_id_umap_mouse", height = "520px")
        )
      ),
      conditionalPanel(
        condition = "input.view_species == 'human'",
        card(
          full_screen = TRUE,
          card_header("Human — UMAP by identity"),
          plotlyOutput("plot_id_umap_human", height = "520px")
        )
      )
    )
  ),
  nav_panel(
    title = "Differential expression",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        selectInput("de_ident_col", "Identity column", choices = default_meta_col),
        selectInput("de_ref", "Reference cell type (ident.1)", choices = NULL),
        selectizeInput(
          "de_other",
          "Comparison cell types (ident.2, pooled)",
          choices = NULL,
          multiple = TRUE
        ),
        selectInput(
          "de_test",
          "DE test",
          choices = c("Wilcoxon rank-sum" = "wilcox"),
          selected = "wilcox"
        ),
        numericInput("de_min_pct", "min.pct", value = 0.1, min = 0, max = 1, step = 0.05),
        actionButton("run_de", "Run FindMarkers", class = "btn-primary w-100"),
        helpText(
          "Compares all cells in the reference type vs the union of selected comparison types on the ",
          "dataset for the species selected in the header (Mouse / Human)."
        )
      ),
      card(
        card_header("DE results"),
        DTOutput("de_table")
      ),
      card(
        full_screen = TRUE,
        card_header("Volcano (Wilcox logFC vs -log10 p)"),
        plotlyOutput("de_volcano", height = "440px")
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    obj_mouse = NULL,
    obj_human = NULL,
    cre_df = NULL,
    cre_sel = NULL,
    gb_loci = NULL,
    gb_overlap = NULL,
    gb_status = "",
    gb_loaded_id = NULL,
    last_de = NULL,
    load_diag = "",
    demo_ref_scheduled = FALSE
  )

  active_seurat <- reactive({
    req(input$view_species)
    if (identical(input$view_species, "mouse")) {
      rv$obj_mouse
    } else {
      rv$obj_human
    }
  })

  cre_bar_mouse_bundle <- reactive({
    obj <- rv$obj_mouse
    sel <- rv$cre_sel
    if (is.null(obj) || is.null(sel) || !nrow(sel)) {
      return(NULL)
    }
    dg <- sel$driver_gene_symbol[[1]]
    if (!nzchar(dg) || is.na(dg) || !dg %in% rownames(obj)) {
      return(NULL)
    }
    colm <- input$id_col
    if (!colm %in% colnames(obj@meta.data)) {
      colm <- colnames(obj@meta.data)[[1L]]
    }
    vx <- as.numeric(Seurat::FetchData(obj, vars = dg)[[1]])
    df <- data.frame(ct = obj@meta.data[[colm]], expr = vx, stringsAsFactors = FALSE) |>
      dplyr::group_by(ct) |>
      dplyr::summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")
    df <- cre_bar_df_fixed_levels(df)
    list(df = df, gene = dg, colm = colm)
  })

  cre_bar_human_bundle <- reactive({
    obj <- rv$obj_human
    sel <- rv$cre_sel
    if (is.null(obj) || is.null(sel) || !nrow(sel)) {
      return(NULL)
    }
    dg <- sel$driver_gene_symbol[[1]]
    if (!nzchar(dg) || is.na(dg)) {
      return(NULL)
    }
    rh <- rownames(obj)
    gplot <- if (dg %in% rh) {
      dg
    } else {
      hg <- toupper(dg)
      if (hg %in% rh) {
        hg
      } else {
        idx <- match(tolower(dg), tolower(rh))
        if (is.na(idx)) {
          return(NULL)
        }
        rh[[idx]]
      }
    }
    colm <- input$id_col
    if (!colm %in% colnames(obj@meta.data)) {
      colm <- colnames(obj@meta.data)[[1L]]
    }
    vx <- as.numeric(Seurat::FetchData(obj, vars = gplot)[[1]])
    df <- data.frame(ct = obj@meta.data[[colm]], expr = vx, stringsAsFactors = FALSE) |>
      dplyr::group_by(ct) |>
      dplyr::summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")
    df <- cre_bar_df_fixed_levels(df)
    list(df = df, gene = gplot, colm = colm)
  })

  output$ui_cre_bar_mouse <- renderUI({
    bd <- cre_bar_mouse_bundle()
    h <- if (!is.null(bd)) cre_bar_plot_height_px(nrow(bd$df)) else 280L
    plotlyOutput("plot_cre_bar_mouse", height = paste0(h, "px"), width = "100%")
  })

  output$ui_cre_bar_human <- renderUI({
    bd <- cre_bar_human_bundle()
    h <- if (!is.null(bd)) cre_bar_plot_height_px(nrow(bd$df)) else 280L
    plotlyOutput("plot_cre_bar_human", height = paste0(h, "px"), width = "100%")
  })

  observeEvent(input$data_mode, ignoreInit = TRUE, {
    rv$obj_mouse <- NULL
    rv$obj_human <- NULL
    if (!identical(input$data_mode, "demo")) {
      rv$demo_ref_scheduled <- FALSE
    }
  })

  observe({
    req(input$data_mode)
    if (!identical(input$data_mode, "demo")) {
      return(invisible(NULL))
    }
    if (isTRUE(rv$demo_ref_scheduled)) {
      return(invisible(NULL))
    }
    if (cremap_local_references_ready()) {
      p <- cremap_reference_paths()
      mm <- read_reference_seurat_file(p$mouse, "mouse")
      hh <- read_reference_seurat_file(p$human, "human")
      if (!is.null(mm) && !is.null(hh)) {
        on.exit(try(removeNotification("boot"), silent = TRUE), add = TRUE)
        showNotification(
          "Loading reference data from reference/*.rds (local files)…",
          id = "boot",
          duration = NULL
        )
        rv$demo_ref_scheduled <- TRUE
        rv$obj_mouse <- mm
        rv$obj_human <- hh
        rv$load_diag <- paste0(
          "Mouse: ", default_data_kind(rv$obj_mouse), "; ",
          "Human: ", default_data_kind(rv$obj_human), "."
        )
        gm <- sort(rownames(rv$obj_mouse))
        gh <- sort(rownames(rv$obj_human))
        sel_m <- if ("Cd19" %in% gm) "Cd19" else gm[[1L]]
        sel_h <- if ("CD19" %in% gh) "CD19" else gh[[1L]]
        updateSelectizeInput(session, "gene_mouse", choices = gm, selected = sel_m, server = TRUE)
        updateSelectizeInput(session, "gene_human", choices = gh, selected = sel_h, server = TRUE)
        removeNotification("boot")
        showNotification(
          paste0(
            "Loaded mouse and human reference from local RDS (n=",
            ncol(rv$obj_mouse), " / ", ncol(rv$obj_human), ")."
          ),
          type = "message",
          duration = 12
        )
        return(invisible(NULL))
      }
    }
    if (!is.null(rv$obj_mouse) && !is.null(rv$obj_human)) {
      return(invisible(NULL))
    }
    on.exit(try(removeNotification("boot"), silent = TRUE), add = TRUE)
    showNotification(
      "Preparing demo data (synthetic first, then reference atlases load in the background)…",
      id = "boot",
      duration = NULL
    )
    rv$demo_ref_scheduled <- TRUE
    rv$load_diag <- "Synthetic placeholder; reference atlases loading in background if packages/network allow."
    ok <- tryCatch(
      {
        rv$obj_mouse <- build_demo_spleno_seurat()
        rv$obj_human <- build_demo_human_spleen_seurat()
        TRUE
      },
      error = function(e) {
        rv$load_diag <- paste0("Synthetic build failed: ", conditionMessage(e))
        FALSE
      }
    )
    if (!ok) {
      rv$demo_ref_scheduled <- FALSE
      showNotification(rv$load_diag, type = "error", duration = NULL)
      return(invisible(NULL))
    }
    gm <- sort(rownames(rv$obj_mouse))
    gh <- sort(rownames(rv$obj_human))
    sel_m <- if ("Cd19" %in% gm) "Cd19" else gm[[1L]]
    sel_h <- if ("CD19" %in% gh) "CD19" else gh[[1L]]
    updateSelectizeInput(session, "gene_mouse", choices = gm, selected = sel_m, server = TRUE)
    updateSelectizeInput(session, "gene_human", choices = gh, selected = sel_h, server = TRUE)
    removeNotification("boot")
    showNotification(
      paste(
        "Interactive plots use synthetic data immediately.",
        "Tabula Muris Senis (mouse) and He organ atlas (human) load next without freezing this tab;",
        "watch the R console for download progress (often 10–30+ minutes on first use)."
      ),
      type = "message",
      duration = 12
    )
    sess <- session
    later::later(
      function() {
        if (tryCatch(isTRUE(sess$isClosed()), error = function(e) TRUE)) {
          return(invisible(NULL))
        }
        shiny::withReactiveDomain(sess, {
          shiny::isolate({
            shiny::showNotification(
              "Loading reference atlases in the background (R console shows hub/network activity)…",
              id = "refbg",
              duration = NULL
            )
            ref_err <- tryCatch(
              {
                shiny::withProgress(
                  message = "Reference atlases (background)",
                  session = sess,
                  value = 0,
                  {
                    shiny::setProgress(0.1, detail = "Mouse: Tabula Muris Senis")
                    m <- load_mouse_spleen_reference()
                    if (!is.null(m)) {
                      rv$obj_mouse <- m
                    }
                    shiny::setProgress(0.55, detail = "Human: He organ atlas")
                    h <- load_human_spleen_reference()
                    if (!is.null(h)) {
                      rv$obj_human <- h
                    }
                    rv$load_diag <- paste0(
                      "Mouse: ", default_data_kind(rv$obj_mouse), "; ",
                      "Human: ", default_data_kind(rv$obj_human), "."
                    )
                    NULL
                  }
                )
                NULL
              },
              error = function(e) {
                em <- conditionMessage(e)
                message("CreMAP: background reference load failed: ", em)
                paste0("Reference upgrade error: ", em)
              }
            )
            gm2 <- sort(rownames(rv$obj_mouse))
            gh2 <- sort(rownames(rv$obj_human))
            sel_m2 <- if ("Cd19" %in% gm2) "Cd19" else gm2[[1L]]
            sel_h2 <- if ("CD19" %in% gh2) "CD19" else gh2[[1L]]
            shiny::updateSelectizeInput(sess, "gene_mouse", choices = gm2, selected = sel_m2, server = TRUE)
            shiny::updateSelectizeInput(sess, "gene_human", choices = gh2, selected = sel_h2, server = TRUE)
            shiny::removeNotification("refbg")
            if (is.null(ref_err)) {
              shiny::showNotification(
                paste("Reference atlases applied.", rv$load_diag),
                type = "message",
                duration = 20
              )
            } else {
              shiny::showNotification(
                paste(ref_err, "— keeping synthetic data."),
                type = "warning",
                duration = 22
              )
            }
          })
        })
      },
      delay = 0.2
    )
  })

  observeEvent(input$reload_data, {
    if (input$data_mode == "demo") {
      rv$load_diag <- ""
      showNotification(
        "Reloading default data in the background (R console shows progress)…",
        id = "reld",
        duration = NULL
      )
      sess <- session
      later::later(
        function() {
          if (tryCatch(isTRUE(sess$isClosed()), error = function(e) TRUE)) {
            return(invisible(NULL))
          }
          shiny::withReactiveDomain(sess, {
            shiny::isolate({
              rel_err <- tryCatch(
                {
                  shiny::withProgress(
                    message = "Reloading default spleen data…",
                    session = sess,
                    value = 0.1,
                    {
                      shiny::setProgress(0.15, detail = "Mouse")
                      m <- load_mouse_spleen_reference()
                      if (is.null(m)) {
                        m <- build_demo_spleno_seurat()
                      }
                      shiny::setProgress(0.55, detail = "Human")
                      h <- load_human_spleen_reference()
                      if (is.null(h)) {
                        h <- build_demo_human_spleen_seurat()
                      }
                      rv$obj_mouse <- m
                      rv$obj_human <- h
                      rv$load_diag <- paste0(
                        "Mouse: ", default_data_kind(rv$obj_mouse), "; ",
                        "Human: ", default_data_kind(rv$obj_human), "."
                      )
                      NULL
                    }
                  )
                  NULL
                },
                error = function(e) {
                  em <- conditionMessage(e)
                  rv$load_diag <- paste0("Error: ", em)
                  m2 <- tryCatch(build_demo_spleno_seurat(), error = function(e2) NULL)
                  h2 <- tryCatch(build_demo_human_spleen_seurat(), error = function(e2) NULL)
                  if (!is.null(m2) && !is.null(h2)) {
                    rv$obj_mouse <- m2
                    rv$obj_human <- h2
                    rv$load_diag <- paste0(rv$load_diag, " Synthetic fallback.")
                  }
                  em
                }
              )
              if (!is.null(rv$obj_mouse) && !is.null(rv$obj_human)) {
                gm <- sort(rownames(rv$obj_mouse))
                gh <- sort(rownames(rv$obj_human))
                sel_m <- if ("Cd19" %in% gm) "Cd19" else gm[[1L]]
                sel_h <- if ("CD19" %in% gh) "CD19" else gh[[1L]]
                shiny::updateSelectizeInput(
                  sess,
                  "gene_mouse",
                  choices = gm,
                  selected = sel_m,
                  server = TRUE
                )
                shiny::updateSelectizeInput(
                  sess,
                  "gene_human",
                  choices = gh,
                  selected = sel_h,
                  server = TRUE
                )
              }
              shiny::removeNotification("reld")
              if (is.null(rel_err)) {
                shiny::showNotification(paste("Default data reloaded.", rv$load_diag), type = "message")
              } else {
                shiny::showNotification(paste("Reload issue:", rv$load_diag), type = "warning", duration = 20)
              }
            })
          })
        },
        delay = 0.05
      )
    } else if (input$data_mode == "cellranger") {
      showNotification("Use “Load Cell Ranger matrix” for 10x data.", type = "message")
    }
  })

  observeEvent(input$load_cellranger, {
    req(input$data_mode == "cellranger")
    maxc <- input$cr_max_cells
    if (!is.numeric(maxc) || maxc < 1L) {
      maxc <- NULL
    } else {
      maxc <- as.integer(maxc)
    }
    res <- tryCatch(
      {
        withProgress(message = "Cell Ranger → Seurat…", value = 0, {
          setProgress(0.05, detail = "Reading matrix")
          obj <- if (identical(input$cr_mode, "zip")) {
            zp <- input$cr_zip
            dp <- if (!is.null(zp)) zp$datapath else ""
            if (is.null(zp) || !nzchar(dp) || !file.exists(dp)) {
              stop("Choose a valid ZIP file.")
            }
            seurat_from_cellranger_zip(
              dp,
              max_cells = maxc,
              resolution = input$cr_resolution,
              npcs = as.integer(input$cr_npc),
              verbose = FALSE
            )
          } else {
            p <- trimws(input$cr_path)
            if (!nzchar(p)) {
              stop("Enter the path to the matrix folder (or a parent directory).")
            }
            if (!dir.exists(p)) {
              stop("That directory does not exist on this machine.")
            }
            mtx <- find_cellranger_matrix_dir(p)
            if (is.na(mtx)) {
              stop("Could not find matrix.mtx or matrix.mtx.gz under that path.")
            }
            setProgress(0.15, detail = basename(mtx))
            seurat_from_cellranger_mtx(
              mtx,
              max_cells = maxc,
              resolution = input$cr_resolution,
              npcs = as.integer(input$cr_npc),
              verbose = FALSE
            )
          }
          setProgress(1, detail = "Done")
          obj
        })
      },
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
        NULL
      }
    )
    if (!is.null(res)) {
      rn <- sort(rownames(res))
      sel_mouse <- if ("Cd19" %in% rn) "Cd19" else rn[[1]]
      sel_human <- if ("CD19" %in% rn) "CD19" else rn[[1]]
      if (identical(input$cr_species, "human")) {
        rv$obj_human <- res
        updateSelectizeInput(session, "gene_human", choices = rn, selected = sel_human, server = TRUE)
      } else {
        rv$obj_mouse <- res
        updateSelectizeInput(session, "gene_mouse", choices = rn, selected = sel_mouse, server = TRUE)
      }
      showNotification(
        paste0("Loaded ", ncol(res), " cells (", input$cr_species, ") from Cell Ranger."),
        type = "message"
      )
    }
  })

  observeEvent(input$rds_file, {
    req(input$data_mode == "rds")
    req(input$rds_file$datapath)
    tryCatch(
      {
        x <- readRDS(input$rds_file$datapath)
        msg <- validate_user_seurat(x)
        if (!is.na(msg)) {
          showNotification(msg, type = "error")
        } else {
          mc <- input$rds_meta_col
          if (!nzchar(mc) || !mc %in% colnames(x@meta.data)) {
            showNotification("Metadata column not found; keeping object Idents as in file.", type = "warning")
          } else {
            Seurat::Idents(x) <- x@meta.data[[mc]]
          }
          if (identical(input$rds_species, "human")) {
            rv$obj_human <- x
            gh <- sort(rownames(x))
            selh <- if ("CD19" %in% gh) "CD19" else gh[[1]]
            updateSelectizeInput(session, "gene_human", choices = gh, selected = selh, server = TRUE)
          } else {
            rv$obj_mouse <- x
            gm <- sort(rownames(x))
            selm <- if ("Cd19" %in% gm) "Cd19" else gm[[1]]
            updateSelectizeInput(session, "gene_mouse", choices = gm, selected = selm, server = TRUE)
          }
          showNotification(paste("Loaded user Seurat object (", input$rds_species, ")."), type = "message")
        }
        NULL
      },
      error = function(e) {
        showNotification(conditionMessage(e), type = "error")
      }
    )
  })

  observeEvent(input$gene_mouse, {
    req(isTRUE(input$sync_human_gene))
    hm <- rv$obj_human
    req(hm)
    gm <- input$gene_mouse
    if (is.null(gm) || !nzchar(gm)) {
      return(invisible(NULL))
    }
    rh <- rownames(hm)
    cand <- toupper(gm)
    pick <- if (cand %in% rh) {
      cand
    } else {
      idx <- match(tolower(gm), tolower(rh))
      if (is.na(idx)) {
        return(invisible(NULL))
      }
      rh[[idx]]
    }
    updateSelectizeInput(
      session,
      "gene_human",
      choices = sort(rh),
      selected = pick,
      server = TRUE
    )
  }, ignoreInit = TRUE)

  observeEvent(input$sync_human_gene, {
    req(isTRUE(input$sync_human_gene))
    hm <- isolate(rv$obj_human)
    if (is.null(hm)) {
      return(invisible(NULL))
    }
    gm <- isolate(input$gene_mouse)
    if (is.null(gm) || !nzchar(gm)) {
      return(invisible(NULL))
    }
    rh <- rownames(hm)
    cand <- toupper(gm)
    pick <- if (cand %in% rh) {
      cand
    } else {
      idx <- match(tolower(gm), tolower(rh))
      if (is.na(idx)) {
        return(invisible(NULL))
      }
      rh[[idx]]
    }
    updateSelectizeInput(
      session,
      "gene_human",
      choices = sort(rownames(hm)),
      selected = pick,
      server = TRUE
    )
  }, ignoreInit = TRUE)

  observe({
    obj <- active_seurat()
    if (is.null(obj)) {
      return(invisible(NULL))
    }
    meta_cols <- colnames(obj@meta.data)
    cur_cm <- isolate(input$color_meta)
    sel_cm <- if (!is.null(cur_cm) && cur_cm %in% meta_cols) cur_cm else meta_cols[[1]]
    updateSelectInput(session, "color_meta", choices = meta_cols, selected = sel_cm)
    sel_id <- if (default_meta_col %in% meta_cols) default_meta_col else meta_cols[[1]]
    updateSelectInput(session, "id_col", choices = meta_cols, selected = sel_id)
    updateSelectInput(session, "de_ident_col", choices = meta_cols, selected = sel_id)
  })

  observe({
    obj <- active_seurat()
    if (is.null(obj)) {
      return(invisible(NULL))
    }
    col <- input$de_ident_col
    req(col %in% colnames(obj@meta.data))
    labs <- sort(unique(as.character(obj@meta.data[[col]])))
    updateSelectInput(session, "de_ref", choices = labs, selected = labs[[1]])
    updateSelectizeInput(session, "de_other", choices = labs, server = TRUE)
  })

  observe({
    obj <- rv$obj_mouse
    req(obj)
    updateSelectizeInput(
      session,
      "cell_pick",
      choices = colnames(obj),
      server = TRUE
    )
  })

  output$data_summary <- renderPrint({
    ld <- rv$load_diag
    if (length(ld) && nzchar(as.character(ld[[1L]]))) {
      cat("Last load note:\n", as.character(ld[[1L]]), "\n\n", sep = "")
    }
    if (!is.null(rv$obj_mouse)) {
      cat("=== Mouse ===\n")
      cat("Default source:", default_data_kind(rv$obj_mouse), "\n")
      cat("Cells:", ncol(rv$obj_mouse), "  Genes:", nrow(rv$obj_mouse), "\n")
      cat("Reductions:", paste(names(rv$obj_mouse@reductions), collapse = ", "), "\n")
      cat("Metadata columns:\n")
      print(colnames(rv$obj_mouse@meta.data))
    } else {
      cat("Mouse: (not loaded)\n")
    }
    cat("\n")
    if (!is.null(rv$obj_human)) {
      cat("=== Human ===\n")
      cat("Default source:", default_data_kind(rv$obj_human), "\n")
      cat("Cells:", ncol(rv$obj_human), "  Genes:", nrow(rv$obj_human), "\n")
      cat("Reductions:", paste(names(rv$obj_human@reductions), collapse = ", "), "\n")
      cat("Metadata columns:\n")
      print(colnames(rv$obj_human@meta.data))
    } else {
      cat("Human: (not loaded)\n")
    }
  })

  output$plot_gene_umap_mouse <- renderPlotly({
    obj <- rv$obj_mouse
    if (is.null(obj)) {
      return(plotly_no_data("No mouse dataset loaded"))
    }
    g <- input$gene_mouse
    if (is.null(g) || !nzchar(g)) {
      return(plotly_no_data("Select a mouse gene"))
    }
    colm <- input$color_meta
    if (!colm %in% colnames(obj@meta.data)) {
      colm <- if (default_meta_col %in% colnames(obj@meta.data)) {
        default_meta_col
      } else {
        colnames(obj@meta.data)[[1]]
      }
    }
    df <- umap_df(obj, color_col = colm, expr_gene = g)
    z <- df[[g]]
    if (isTRUE(input$gene_log)) {
      z <- log1p(z)
    }
    df$expr_plot <- z
    has_meta <- nzchar(colm) && colm %in% colnames(df)
    df$tip <- if (has_meta) paste(df$cell, df[[colm]], sep = "<br>") else df$cell
    plot_ly(
      df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      color = ~expr_plot,
      text = ~tip,
      type = "scattergl",
      mode = "markers",
      marker = list(size = input$pt_size, opacity = input$pt_alpha)
    ) |>
      layout(
        title = list(text = g, font = list(size = 12)),
        xaxis = list(title = "UMAP_1"),
        yaxis = list(title = "UMAP_2"),
        legend = list(title = list(text = if (isTRUE(input$gene_log)) "log1p expr" else "expr"))
      )
  })

  output$plot_gene_umap_human <- renderPlotly({
    obj <- rv$obj_human
    if (is.null(obj)) {
      return(plotly_no_data("No human dataset loaded"))
    }
    g <- input$gene_human
    if (is.null(g) || !nzchar(g)) {
      return(plotly_no_data("Select a human gene"))
    }
    colm <- input$color_meta
    if (!colm %in% colnames(obj@meta.data)) {
      colm <- if (default_meta_col %in% colnames(obj@meta.data)) {
        default_meta_col
      } else {
        colnames(obj@meta.data)[[1]]
      }
    }
    df <- umap_df(obj, color_col = colm, expr_gene = g)
    z <- df[[g]]
    if (isTRUE(input$gene_log)) {
      z <- log1p(z)
    }
    df$expr_plot <- z
    has_meta <- nzchar(colm) && colm %in% colnames(df)
    df$tip <- if (has_meta) paste(df$cell, df[[colm]], sep = "<br>") else df$cell
    plot_ly(
      df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      color = ~expr_plot,
      text = ~tip,
      type = "scattergl",
      mode = "markers",
      marker = list(size = input$pt_size, opacity = input$pt_alpha)
    ) |>
      layout(
        title = list(text = g, font = list(size = 12)),
        xaxis = list(title = "UMAP_1"),
        yaxis = list(title = "UMAP_2"),
        legend = list(title = list(text = if (isTRUE(input$gene_log)) "log1p expr" else "expr"))
      )
  })

  output$plot_gene_bar_mouse <- renderPlotly({
    obj <- rv$obj_mouse
    if (is.null(obj)) {
      return(plotly_no_data("No mouse data"))
    }
    g <- input$gene_mouse
    if (is.null(g) || !nzchar(g) || !g %in% rownames(obj)) {
      return(plotly_empty())
    }
    colm <- input$id_col
    if (!colm %in% colnames(obj@meta.data)) {
      colm <- colnames(obj@meta.data)[[1]]
    }
    vx <- as.numeric(Seurat::FetchData(obj, vars = g)[[1]])
    df <- data.frame(
      ct = obj@meta.data[[colm]],
      expr = vx,
      stringsAsFactors = FALSE
    ) |>
      dplyr::group_by(ct) |>
      dplyr::summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")
    plot_ly(df, x = ~ct, y = ~mean_expr, type = "bar") |>
      layout(title = g, xaxis = list(title = ""), yaxis = list(title = "Mean expr"))
  })

  output$plot_gene_bar_human <- renderPlotly({
    obj <- rv$obj_human
    if (is.null(obj)) {
      return(plotly_no_data("No human data"))
    }
    g <- input$gene_human
    if (is.null(g) || !nzchar(g) || !g %in% rownames(obj)) {
      return(plotly_empty())
    }
    colm <- input$id_col
    if (!colm %in% colnames(obj@meta.data)) {
      colm <- colnames(obj@meta.data)[[1]]
    }
    vx <- as.numeric(Seurat::FetchData(obj, vars = g)[[1]])
    df <- data.frame(
      ct = obj@meta.data[[colm]],
      expr = vx,
      stringsAsFactors = FALSE
    ) |>
      dplyr::group_by(ct) |>
      dplyr::summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")
    plot_ly(df, x = ~ct, y = ~mean_expr, type = "bar") |>
      layout(title = g, xaxis = list(title = ""), yaxis = list(title = "Mean expr"))
  })

  observeEvent(input$cre_fetch, {
    withProgress(message = "MouseMine query…", {
      df <- fetch_cre_alleles(search = input$cre_search, limit = input$cre_limit)
      rv$cre_df <- df
      if (!nrow(df)) {
        showNotification("No alleles returned (network or filter).", type = "warning")
      } else {
        showNotification(paste("Fetched", nrow(df), "alleles."), type = "message")
      }
    })
  })

  output$cre_table <- renderDT({
    df <- rv$cre_df
    validate(need(!is.null(df), 'Click "Query MouseMine" to load Cre alleles.'))
    ids <- trimws(as.character(df$mgi_id))
    urls <- mgi_allele_summary_url(df$mgi_id)
    has_link <- !is.na(urls)
    mgi_html <- ifelse(
      has_link,
      sprintf(
        "<a href=\"%s\" target=\"_blank\" rel=\"noopener noreferrer\">%s</a>",
        htmltools::htmlEscape(urls, attribute = TRUE),
        htmltools::htmlEscape(ids, attribute = FALSE)
      ),
      ifelse(is.na(ids), "", htmltools::htmlEscape(ids, attribute = FALSE))
    )
    esc_txt <- function(x) {
      x <- as.character(x)
      x[is.na(x)] <- ""
      htmltools::htmlEscape(x, attribute = FALSE)
    }
    disp <- data.frame(
      mgi_html,
      esc_txt(df$driver_gene_symbol),
      esc_txt(df$symbol),
      esc_txt(df$name),
      stringsAsFactors = FALSE
    )
    colnames(disp) <- c(
      "MGI allele reference number",
      "Gene name",
      "Model name",
      "Description"
    )
    datatable(
      disp,
      selection = "single",
      filter = "top",
      rownames = FALSE,
      class = "display compact nowrap",
      escape = FALSE,
      options = list(
        pageLength = 5,
        lengthMenu = c(5, 10, 20),
        scrollX = TRUE,
        scrollY = "168px",
        scrollCollapse = TRUE,
        searching = TRUE,
        info = FALSE,
        dom = "lftip"
      )
    )
  })

  observeEvent(input$cre_table_rows_selected, {
    df <- rv$cre_df
    req(nrow(df))
    i <- input$cre_table_rows_selected
    req(length(i) == 1L)
    rv$cre_sel <- df[i, , drop = FALSE]
  })

  load_gb_context <- function(query_raw) {
    raw <- trimws(as.character(query_raw))
    if (!length(raw) || is.na(raw) || !nzchar(raw)) {
      rv$gb_loci <- NULL
      rv$gb_overlap <- NULL
      rv$gb_loaded_id <- NULL
      rv$gb_status <- "Enter an MGI allele ID (MGI:#####) or allele symbol."
      showNotification(rv$gb_status, type = "warning")
      return(invisible(FALSE))
    }
    mgi_id <- normalize_gb_query_id(raw)
    if (is.na(mgi_id) || !nzchar(mgi_id)) {
      rv$gb_loci <- NULL
      rv$gb_overlap <- NULL
      rv$gb_loaded_id <- NULL
      rv$gb_status <- paste0(
        "No allele found for \"",
        raw,
        "\". Try an MGI:##### ID or a symbol substring (e.g. Nes-cre)."
      )
      showNotification(rv$gb_status, type = "warning")
      return(invisible(FALSE))
    }
    loci <- NULL
    overlap <- NULL
    withProgress(message = "MouseMine — integration locus…", {
      loci <- fetch_allele_integration_sites(mgi_id)
      overlap <- fetch_allele_overlapping_features(mgi_id)
    })
    rv$gb_loci <- if (!is.null(loci) && nrow(loci)) loci else NULL
    rv$gb_overlap <- if (!is.null(overlap) && nrow(overlap)) overlap else NULL
    rv$gb_loaded_id <- mgi_id
    if (is.null(rv$gb_loci)) {
      rv$gb_status <- paste0(
        "Allele ",
        mgi_id,
        " found, but MouseMine has no mapped insertion coordinates. ",
        "Some lines lack genomic intervals in MGI."
      )
      showNotification(rv$gb_status, type = "warning")
    } else {
      n_ov <- if (is.null(rv$gb_overlap)) 0L else nrow(rv$gb_overlap)
      sym <- rv$gb_loci$symbol[[1]]
      rv$gb_status <- paste0(
        "Loaded ",
        mgi_id,
        if (nzchar(sym)) paste0(" (", sym, ")") else "",
        " — ",
        nrow(rv$gb_loci),
        " interval(s), ",
        n_ov,
        " overlapping feature(s)."
      )
      showNotification(rv$gb_status, type = "message")
    }
    invisible(!is.null(rv$gb_loci))
  }

  observeEvent(input$gb_fetch, {
    load_gb_context(input$gb_query)
  })

  observeEvent(rv$cre_sel, {
    sel <- rv$cre_sel
    if (is.null(sel) || !nrow(sel)) {
      return()
    }
    mid <- trimws(as.character(sel$mgi_id[[1]]))
    if (!nzchar(mid) || is.na(mid)) {
      return()
    }
    updateTextInput(session, "gb_query", value = mid)
    if (isTRUE(input$gb_auto_load)) {
      load_gb_context(mid)
    }
  }, ignoreInit = TRUE)

  output$gb_summary <- renderUI({
    loci <- rv$gb_loci
    status <- rv$gb_status
    if (is.null(loci) || !nrow(loci)) {
      return(tags$p(class = "text-muted mb-0", status %||% "Load an allele to view integration data."))
    }
    loc <- loci[1, , drop = FALSE]
    mid <- loc$mgi_id[[1]]
    sym <- loc$symbol[[1]]
    dg <- loc$driver_gene_symbol[[1]]
    url <- mgi_allele_summary_url(mid)
    chr <- loc$chromosome[[1]]
    a0 <- loc$start[[1]]
    b0 <- loc$end[[1]]
    asm <- loc$assembly[[1]]
    strand <- loc$strand[[1]]
    if (identical(strand, "-1")) {
      strand_lab <- "reverse (−)"
    } else if (identical(strand, "1")) {
      strand_lab <- "forward (+)"
    } else {
      strand_lab <- if (nzchar(strand)) strand else "—"
    }
    size_bp <- b0 - a0 + 1L
    overlap <- rv$gb_overlap
    confound <- character()
    if (nzchar(dg) && !is.na(dg) && !is.null(overlap) && nrow(overlap)) {
      if (dg %in% overlap$symbol) {
        confound <- c(
          confound,
          tags$div(
            class = "alert alert-warning py-2 small mb-2",
            tags$strong("Driver gene at insertion locus: "),
            "MGI reports ",
            tags$code(dg),
            " among overlapping features. Expression driven by this line may reflect ",
            "insertional effects on the endogenous locus as well as Cre activity."
          )
        )
      }
    }
    if (!is.null(overlap) && nrow(overlap) >= 3L) {
      confound <- c(
        confound,
        tags$div(
          class = "alert alert-info py-2 small mb-2",
          "Multiple genes overlap this insertion interval (",
          nrow(overlap),
          " in MGI). Review the table and UCSC view for nearby endogenous genes."
        )
      )
    }
    tags$div(
      confound,
      tags$dl(
        class = "row small mb-0",
        tags$dt(class = "col-sm-3", "MGI allele"),
        tags$dd(
          class = "col-sm-9",
          if (!is.na(url)) {
            tags$a(href = url, target = "_blank", rel = "noopener noreferrer", mid)
          } else {
            mid
          }
        ),
        tags$dt(class = "col-sm-3", "Symbol"),
        tags$dd(class = "col-sm-9", sym),
        tags$dt(class = "col-sm-3", "Type"),
        tags$dd(class = "col-sm-9", loc$allele_type[[1]] %||% "—"),
        tags$dt(class = "col-sm-3", "Driver (MGI)"),
        tags$dd(class = "col-sm-9", if (nzchar(dg) && !is.na(dg)) dg else "—"),
        tags$dt(class = "col-sm-3", "Assembly"),
        tags$dd(class = "col-sm-9", asm),
        tags$dt(class = "col-sm-3", "Insertion"),
        tags$dd(
          class = "col-sm-9",
          sprintf("chr%s:%s–%s (%s strand; %s bp)", chr, format(a0, big.mark = ","), format(b0, big.mark = ","), strand_lab, format(size_bp, big.mark = ","))
        ),
        tags$dt(class = "col-sm-3", "Description"),
        tags$dd(class = "col-sm-9", loc$name[[1]] %||% "—")
      ),
      tags$p(class = "text-muted small mt-2 mb-0", status)
    )
  })

  output$gb_quick_links <- renderUI({
    mid <- rv$gb_loaded_id
    loci <- rv$gb_loci
    if (is.null(mid) || !nzchar(mid)) {
      return(tags$p(class = "text-muted small mb-0", "Quick links appear after loading."))
    }
    url_mgi <- mgi_allele_summary_url(mid)
    links <- list(
      if (!is.na(url_mgi)) {
        tags$a(href = url_mgi, target = "_blank", rel = "noopener noreferrer", "MGI allele report")
      }
    )
    if (!is.null(loci) && nrow(loci)) {
      loc <- loci[1, , drop = FALSE]
      u <- ucsc_browser_url(
        loc$chromosome[[1]],
        loc$start[[1]],
        loc$end[[1]],
        loc$assembly[[1]],
        pad_bp = input$gb_ucsc_pad
      )
      if (!is.na(u)) {
        links <- c(links, list(tags$a(href = u, target = "_blank", rel = "noopener noreferrer", "Open UCSC in new tab")))
      }
    }
    tags$div(class = "d-flex flex-column gap-1", links)
  })

  output$gb_ucsc_frame <- renderUI({
    loci <- rv$gb_loci
    if (is.null(loci) || !nrow(loci)) {
      return(
        tags$p(
          class = "text-muted small mb-0 px-2 py-3",
          rv$gb_status %||% "Load an allele with mapped coordinates to embed UCSC."
        )
      )
    }
    loc <- loci[1, , drop = FALSE]
    u <- ucsc_browser_url(
      loc$chromosome[[1]],
      loc$start[[1]],
      loc$end[[1]],
      loc$assembly[[1]],
      pad_bp = input$gb_ucsc_pad
    )
    db <- assembly_to_ucsc_db(loc$assembly[[1]])
    if (is.na(u)) {
      return(tags$p(class = "text-warning small mb-0", "UCSC embed not available for assembly ", loc$assembly[[1]], "."))
    }
    tags$div(
      tags$p(class = "small text-muted mb-1", "Embedded ", tags$code(db), " view (GRCm39 → mm39)."),
      tags$iframe(
        src = u,
        width = "100%",
        height = "420",
        style = "border:1px solid #dee2e6;border-radius:4px;",
        title = "UCSC Genome Browser"
      )
    )
  })

  output$gb_loci_table <- renderDT({
    loci <- rv$gb_loci
    validate(need(!is.null(loci) && nrow(loci), rv$gb_status %||% "Load an allele first."))
    disp <- loci |>
      dplyr::transmute(
        Chromosome = .data$chromosome,
        Start = .data$start,
        End = .data$end,
        `Size (bp)` = .data$end - .data$start + 1L,
        Strand = .data$strand,
        Assembly = .data$assembly,
        `Feature ID` = .data$feature_id
      )
    datatable(disp, rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  output$gb_overlap_table <- renderDT({
    overlap <- rv$gb_overlap
    loci <- rv$gb_loci
    validate(need(!is.null(loci), "Load an allele with coordinates first."))
    if (is.null(overlap) || !nrow(overlap)) {
      validate(need(FALSE, "No overlapping features returned from MouseMine for this insertion."))
    }
    dg <- loci$driver_gene_symbol[[1]]
    disp <- overlap |>
      dplyr::mutate(
        `Driver locus?` = ifelse(!is.na(dg) & nzchar(dg) & .data$symbol == dg, "yes", "")
      ) |>
      dplyr::select(Symbol = .data$symbol, Name = .data$name, `MGI feature` = .data$feature_id, `Driver locus?`)
    datatable(
      disp,
      rownames = FALSE,
      options = list(pageLength = 12, scrollX = TRUE)
    ) |>
      DT::formatStyle(
        columns = "Driver locus?",
        backgroundColor = DT::styleEqual("yes", "#fff3cd")
      )
  })

  output$gb_locus_plot <- renderPlotly({
    loci <- rv$gb_loci
    sym <- if (!is.null(loci) && nrow(loci)) loci$symbol[[1]] else ""
    plot_gb_locus(loci, rv$gb_overlap, symbol = sym %||% "")
  })

  output$cre_driver_status <- renderText({
    om <- rv$obj_mouse
    oh <- rv$obj_human
    sel <- rv$cre_sel
    if (is.null(sel) || !nrow(sel)) {
      return("Select a Cre allele from the table.")
    }
    dg <- sel$driver_gene_symbol[[1]]
    if (!nzchar(dg) || is.na(dg)) {
      return(paste0("Allele: ", sel$symbol[[1]], "\nNo driver gene symbol in MouseMine for this row."))
    }
    lines <- paste0("Driver gene symbol from MGI: ", dg, "\n")
    if (is.null(om)) {
      lines <- paste0(lines, "Mouse: (not loaded)\n")
    } else if (dg %in% rownames(om)) {
      lines <- paste0(lines, "Mouse: present in RNA assay rownames.\n")
    } else {
      lines <- paste0(lines, "Mouse: not in RNA assay rownames.\n")
    }
    hg <- toupper(dg)
    if (is.null(oh)) {
      lines <- paste0(lines, "Human: (not loaded)")
    } else if (hg %in% rownames(oh)) {
      lines <- paste0(lines, "Human: present as ", hg, " in RNA assay rownames.")
    } else {
      idx <- match(tolower(dg), tolower(rownames(oh)))
      if (!is.na(idx)) {
        lines <- paste0(lines, "Human: present as ", rownames(oh)[[idx]], " in RNA assay rownames.")
      } else {
        lines <- paste0(lines, "Human: no matching symbol in RNA assay rownames.")
      }
    }
    lines
  })

  output$plot_cre_bar_mouse <- renderPlotly({
    bd <- cre_bar_mouse_bundle()
    if (is.null(bd)) {
      obj <- rv$obj_mouse
      sel <- rv$cre_sel
      if (is.null(obj) || is.null(sel) || !nrow(sel)) {
        return(plotly_no_data("No mouse data or Cre selection"))
      }
      return(plotly_empty())
    }
    df <- bd$df
    h_px <- cre_bar_plot_height_px(nrow(df))
    n <- nrow(df)
    tick_sz <- max(7L, min(11L, as.integer(220 / sqrt(n + 4))))
    pad_l <- max(72L, min(400L, as.integer(5.5 * max(nchar(as.character(df$ct)), na.rm = TRUE))))
    ylev <- levels(df$ct)
    plot_ly(
      df,
      x = ~mean_expr,
      y = ~ct,
      type = "bar",
      orientation = "h",
      hovertemplate = "%{y}<br>mean: %{x:.3f}<extra></extra>"
    ) |>
      layout(
        title = list(
          text = paste0("Mouse — ", bd$gene, "\nMean by ", bd$colm),
          font = list(size = 13),
          x = 0,
          xanchor = "left",
          pad = list(t = 2, b = 0)
        ),
        height = h_px,
        autosize = TRUE,
        margin = list(l = pad_l, r = 20, t = 52, b = 24),
        xaxis = list(
          title = list(text = "Mean expression", font = list(size = 11)),
          automargin = TRUE
        ),
        yaxis = list(
          title = "",
          automargin = TRUE,
          tickfont = list(size = tick_sz),
          categoryorder = "array",
          categoryarray = ylev
        ),
        bargap = 0.12
      ) |>
      plotly::config(responsive = TRUE, displayModeBar = TRUE)
  })

  output$plot_cre_bar_human <- renderPlotly({
    bd <- cre_bar_human_bundle()
    if (is.null(bd)) {
      obj <- rv$obj_human
      sel <- rv$cre_sel
      if (is.null(obj) || is.null(sel) || !nrow(sel)) {
        return(plotly_no_data("No human data or Cre selection"))
      }
      return(plotly_empty())
    }
    df <- bd$df
    h_px <- cre_bar_plot_height_px(nrow(df))
    n <- nrow(df)
    tick_sz <- max(7L, min(11L, as.integer(220 / sqrt(n + 4))))
    pad_l <- max(72L, min(400L, as.integer(5.5 * max(nchar(as.character(df$ct)), na.rm = TRUE))))
    ylev <- levels(df$ct)
    plot_ly(
      df,
      x = ~mean_expr,
      y = ~ct,
      type = "bar",
      orientation = "h",
      hovertemplate = "%{y}<br>mean: %{x:.3f}<extra></extra>"
    ) |>
      layout(
        title = list(
          text = paste0("Human — ", bd$gene, "\nMean by ", bd$colm),
          font = list(size = 13),
          x = 0,
          xanchor = "left",
          pad = list(t = 2, b = 0)
        ),
        height = h_px,
        autosize = TRUE,
        margin = list(l = pad_l, r = 20, t = 52, b = 24),
        xaxis = list(
          title = list(text = "Mean expression", font = list(size = 11)),
          automargin = TRUE
        ),
        yaxis = list(
          title = "",
          automargin = TRUE,
          tickfont = list(size = tick_sz),
          categoryorder = "array",
          categoryarray = ylev
        ),
        bargap = 0.12
      ) |>
      plotly::config(responsive = TRUE, displayModeBar = TRUE)
  })

  output$plot_cre_umap_mouse <- renderPlotly({
    obj <- rv$obj_mouse
    sel <- rv$cre_sel
    if (is.null(obj) || is.null(sel) || !nrow(sel)) {
      return(plotly_no_data("No mouse data or Cre selection"))
    }
    dg <- sel$driver_gene_symbol[[1]]
    if (!nzchar(dg) || is.na(dg) || !dg %in% rownames(obj)) {
      return(plotly_empty())
    }
    df <- umap_df(obj, expr_gene = dg)
    df$expr_cre <- log1p(df[[dg]])
    plot_ly(
      df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      color = ~expr_cre,
      text = ~cell,
      type = "scattergl",
      mode = "markers",
      marker = list(size = input$pt_size, opacity = input$pt_alpha)
    ) |>
      layout(title = paste("Mouse:", dg), xaxis = list(title = "UMAP_1"), yaxis = list(title = "UMAP_2"))
  })

  output$plot_cre_umap_human <- renderPlotly({
    obj <- rv$obj_human
    sel <- rv$cre_sel
    if (is.null(obj) || is.null(sel) || !nrow(sel)) {
      return(plotly_no_data("No human data or Cre selection"))
    }
    dg <- sel$driver_gene_symbol[[1]]
    if (!nzchar(dg) || is.na(dg)) {
      return(plotly_empty())
    }
    rh <- rownames(obj)
    gplot <- if (dg %in% rh) {
      dg
    } else {
      hg <- toupper(dg)
      if (hg %in% rh) {
        hg
      } else {
        idx <- match(tolower(dg), tolower(rh))
        if (is.na(idx)) {
          return(plotly_empty())
        }
        rh[[idx]]
      }
    }
    df <- umap_df(obj, expr_gene = gplot)
    df$expr_cre <- log1p(df[[gplot]])
    plot_ly(
      df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      color = ~expr_cre,
      text = ~cell,
      type = "scattergl",
      mode = "markers",
      marker = list(size = input$pt_size, opacity = input$pt_alpha)
    ) |>
      layout(title = paste("Human:", gplot), xaxis = list(title = "UMAP_1"), yaxis = list(title = "UMAP_2"))
  })

  output$plot_id_umap_mouse <- renderPlotly({
    obj <- rv$obj_mouse
    if (is.null(obj)) {
      return(plotly_no_data("No mouse dataset"))
    }
    colm <- input$id_col
    if (!colm %in% colnames(obj@meta.data)) {
      return(plotly_no_data("Identity column missing in mouse"))
    }
    df <- umap_df(obj, color_col = colm)
    df$color_lab <- as.character(df[[colm]])
    df$tip_id <- paste(df$cell, df$color_lab, sep = "<br>")
    plot_ly(
      df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      color = ~color_lab,
      text = ~tip_id,
      type = "scattergl",
      mode = "markers",
      marker = list(size = input$pt_size, opacity = input$pt_alpha)
    ) |>
      layout(xaxis = list(title = "UMAP_1"), yaxis = list(title = "UMAP_2"))
  })

  output$plot_id_umap_human <- renderPlotly({
    obj <- rv$obj_human
    if (is.null(obj)) {
      return(plotly_no_data("No human dataset"))
    }
    colm <- input$id_col
    if (!colm %in% colnames(obj@meta.data)) {
      colm <- colnames(obj@meta.data)[[1]]
    }
    df <- umap_df(obj, color_col = colm)
    df$color_lab <- as.character(df[[colm]])
    df$tip_id <- paste(df$cell, df$color_lab, sep = "<br>")
    plot_ly(
      df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      color = ~color_lab,
      text = ~tip_id,
      type = "scattergl",
      mode = "markers",
      marker = list(size = input$pt_size, opacity = input$pt_alpha)
    ) |>
      layout(xaxis = list(title = "UMAP_1"), yaxis = list(title = "UMAP_2"))
  })

  observeEvent(input$apply_suggest, {
    nms <- character()
    if (!is.null(rv$obj_mouse)) {
      x <- rv$obj_mouse
      sug <- suggest_cell_types_from_markers(x, SPLENO_MARKERS)
      x$cell_type_suggested <- sug
      Seurat::Idents(x) <- sug
      rv$obj_mouse <- x
      nms <- c(nms, "mouse")
    }
    if (!is.null(rv$obj_human)) {
      x <- rv$obj_human
      sug <- suggest_cell_types_from_markers(x, HUMAN_SPLEEN_MARKERS)
      x$cell_type_suggested <- sug
      Seurat::Idents(x) <- sug
      rv$obj_human <- x
      nms <- c(nms, "human")
    }
    if (!length(nms)) {
      showNotification("No datasets loaded.", type = "warning")
    } else {
      showNotification(paste("Suggested labels for:", paste(nms, collapse = ", ")), type = "message")
    }
  })

  observeEvent(input$apply_manual, {
    obj <- rv$obj_mouse
    req(obj)
    cell <- input$cell_pick
    lab <- trimws(input$manual_label)
    req(nzchar(cell), nzchar(lab))
    idc <- input$id_col
    req(idc %in% colnames(obj@meta.data))
    if (!"cell_type_manual" %in% colnames(obj@meta.data)) {
      obj$cell_type_manual <- as.character(obj@meta.data[[idc]])
    }
    obj$cell_type_manual <- as.character(obj$cell_type_manual)
    obj$cell_type_manual[[cell]] <- lab
    Seurat::Idents(obj) <- obj$cell_type_manual
    rv$obj_mouse <- obj
    showNotification(paste("Updated mouse cell", cell), type = "message")
  })

  observeEvent(input$run_de, {
    obj <- active_seurat()
    req(obj)
    colm <- input$de_ident_col
    req(colm %in% colnames(obj@meta.data))
    ref <- input$de_ref
    others <- input$de_other
    others <- setdiff(others, ref)
    req(length(others) >= 1L)

    withProgress(message = "FindMarkers…", {
      Seurat::Idents(obj) <- obj@meta.data[[colm]]
      markers <- Seurat::FindMarkers(
        obj,
        ident.1 = ref,
        ident.2 = others,
        test.use = input$de_test,
        min.pct = input$de_min_pct,
        verbose = FALSE
      )
      markers$gene <- rownames(markers)
      rv$last_de <- markers
    })
    showNotification("DE finished.", type = "message")
  })

  output$de_table <- renderDT({
    m <- rv$last_de
    validate(need(!is.null(m), "Run FindMarkers to populate results."))
    datatable(
      m |> dplyr::arrange(p_val) |> dplyr::select(gene, dplyr::everything()),
      filter = "top",
      options = list(pageLength = 15, scrollX = TRUE)
    )
  })

  output$de_volcano <- renderPlotly({
    m <- rv$last_de
    validate(need(!is.null(m), "Run FindMarkers first."))
    fc_col <- if ("avg_log2FC" %in% colnames(m)) "avg_log2FC" else "avg_logFC"
    validate(need(fc_col %in% colnames(m) && "p_val" %in% colnames(m), "Unexpected DE table columns."))
    df <- m |>
      dplyr::mutate(
        neglog10p = -log10(pmax(.data$p_val, 1e-300)),
        gene = rownames(m),
        fc = .data[[fc_col]]
      )
    plot_ly(
      df,
      x = ~fc,
      y = ~neglog10p,
      text = ~gene,
      type = "scattergl",
      mode = "markers",
      marker = list(size = 6, opacity = 0.65, color = "#34495e")
    ) |>
      layout(
        xaxis = list(title = fc_col),
        yaxis = list(title = "-log10 p")
      )
  })
}

shinyApp(ui, server)
