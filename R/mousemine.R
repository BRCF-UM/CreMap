#' MouseMine (MGI) helpers — Cre recombinase alleles
#'
#' Data source: https://www.mousemine.org/mousemine (Mouse Genome Informatics).

MOUSEMINE_BASE <- "https://www.mousemine.org/mousemine/service"

#' Escape XML special characters for path-query fragments
xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&apos;", x, fixed = TRUE)
  x
}

#' Build InterMine XML query string for recombinase alleles
cre_allele_query_xml <- function(search = "", limit = 500L) {
  search <- trimws(search)
  base_view <- paste(
    "Allele.primaryIdentifier",
    "Allele.symbol",
    "Allele.name",
    "Allele.drivenBy.feature.symbol",
    sep = " "
  )
  if (nzchar(search)) {
    pat <- xml_escape(search)
    constraints <- sprintf(
      "<constraint path=\"Allele.isRecombinase\" op=\"=\" value=\"true\" code=\"A\"/>%s%s%s",
      "<constraint path=\"Allele.symbol\" op=\"CONTAINS\" value=\"",
      pat,
      "\" code=\"B\"/>"
    )
    logic <- 'constraintLogic="A and B"'
  } else {
    constraints <- "<constraint path=\"Allele.isRecombinase\" op=\"=\" value=\"true\"/>"
    logic <- ""
  }
  sprintf(
    "<query name=\"\" model=\"genomic\" view=\"%s\" sortOrder=\"Allele.symbol asc\" %s>%s</query>",
    base_view,
    logic,
    constraints
  )
}

#' Build MGI allele summary page URL from MouseMine primary identifier
#'
#' @param primary_ids Character vector of MGI allele IDs (e.g. `MGI:#####`).
#' @return Character vector of URLs; `NA_character_` where id is missing or empty.
#' @noRd
mgi_allele_summary_url <- function(primary_ids) {
  ids <- trimws(as.character(primary_ids))
  out <- rep(NA_character_, length(ids))
  ok <- !is.na(ids) & nzchar(ids)
  out[ok] <- paste0("https://www.informatics.jax.org/allele/", ids[ok])
  out
}

#' Fetch Cre alleles from MouseMine
#'
#' @param search Optional substring filter on allele symbol.
#' @param limit Max rows (capped at 2000).
#' @return A data.frame with columns `mgi_id` (MGI allele ID), `symbol` (allele
#'   symbol / model shorthand), `name` (long description), `driver_gene_symbol`
#'   (driver gene when annotated). The Shiny table shows these in display order:
#'   MGI reference, gene name, model name, description.
#' @export
fetch_cre_alleles <- function(search = "", limit = 500L) {
  limit <- min(max(as.integer(limit), 1L), 2000L)
  q <- cre_allele_query_xml(search, limit)
  req <- httr2::request(MOUSEMINE_BASE) |>
    httr2::req_url_path_append("query/results") |>
    httr2::req_url_query(query = q, format = "json", size = limit) |>
    httr2::req_timeout(120) |>
    httr2::req_error(is_error = function(resp) FALSE)

  resp <- tryCatch(httr2::req_perform(req), error = function(e) NULL)
  if (is.null(resp) || httr2::resp_status(resp) >= 400) {
    return(data.frame(
      mgi_id = character(),
      symbol = character(),
      name = character(),
      driver_gene_symbol = character(),
      stringsAsFactors = FALSE
    ))
  }

  body <- httr2::resp_body_json(resp, simplifyVector = FALSE)
  cn <- c("mgi_id", "symbol", "name", "driver_gene_symbol")
  mousemine_results_df(body, cn)
}

#' Run a MouseMine InterMine XML query and return parsed JSON body
#'
#' @param query_xml InterMine query XML string.
#' @param limit Max rows.
#' @return Parsed list from MouseMine, or `NULL` on transport / HTTP failure.
#' @noRd
mousemine_query_json <- function(query_xml, limit = 500L) {
  limit <- min(max(as.integer(limit), 1L), 2000L)
  req <- httr2::request(MOUSEMINE_BASE) |>
    httr2::req_url_path_append("query/results") |>
    httr2::req_url_query(query = query_xml, format = "json", size = limit) |>
    httr2::req_timeout(120) |>
    httr2::req_error(is_error = function(resp) FALSE)
  resp <- tryCatch(httr2::req_perform(req), error = function(e) NULL)
  if (is.null(resp) || httr2::resp_status(resp) >= 400) {
    return(NULL)
  }
  httr2::resp_body_json(resp, simplifyVector = FALSE)
}

#' Map MGI assembly name to a UCSC `db` parameter
#'
#' @param assembly Assembly string from MouseMine (e.g. `GRCm39`).
#' @return UCSC database name or `NULL` if unknown.
#' @export
assembly_to_ucsc_db <- function(assembly) {
  a <- toupper(trimws(as.character(assembly)))
  a <- a[!is.na(a) & nzchar(a)]
  if (!length(a)) {
    return(NULL)
  }
  a <- a[[1L]]
  switch(
    a,
    "GRCM39" = "mm39",
    "GRCM38" = "mm10",
    "NCBIM37" = "mm9",
    "MGSCV37" = "mm9",
    NULL
  )
}

#' Format MGI chromosome id for UCSC (e.g. `12` -> `chr12`)
#'
#' @param chr_id Chromosome primary identifier from MouseMine.
#' @return UCSC-style chromosome name.
#' @export
mgi_chr_ucsc <- function(chr_id) {
  x <- trimws(as.character(chr_id))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) {
    return(NA_character_)
  }
  x <- x[[1L]]
  if (grepl("^chr", x, ignore.case = TRUE)) {
    return(x)
  }
  paste0("chr", x)
}

#' Build a UCSC Genome Browser URL centered on a locus
#'
#' @param chr UCSC chromosome (e.g. `chr12`).
#' @param start Start coordinate (1-based).
#' @param end End coordinate (1-based).
#' @param assembly MGI assembly label (e.g. `GRCm39`).
#' @param pad_bp Padding added on each side of the interval.
#' @return Character URL or `NA_character_` when coordinates are invalid.
#' @export
ucsc_browser_url <- function(chr, start, end, assembly = "GRCm39", pad_bp = 50000L) {
  db <- assembly_to_ucsc_db(assembly)
  if (is.null(db)) {
    return(NA_character_)
  }
  chr <- mgi_chr_ucsc(chr)
  start <- suppressWarnings(as.integer(start))
  end <- suppressWarnings(as.integer(end))
  pad_bp <- max(0L, as.integer(pad_bp))
  if (is.na(chr) || !nzchar(chr) || is.na(start) || is.na(end)) {
    return(NA_character_)
  }
  if (end < start) {
    tmp <- start
    start <- end
    end <- tmp
  }
  win_start <- max(1L, start - pad_bp)
  win_end <- end + pad_bp
  pos <- sprintf("%s:%d-%d", chr, win_start, win_end)
  paste0(
    "https://genome.ucsc.edu/cgi-bin/hgTracks?db=",
    utils::URLencode(db, reserved = TRUE),
    "&position=",
    utils::URLencode(pos, reserved = TRUE),
    "&knownGene=pack&refSeqComposite=full&ncbiRefSeqCurated=full&ignoreCookie=1"
  )
}

allele_genomic_query_xml <- function(mgi_id, limit = 50L) {
  mid <- xml_escape(trimws(as.character(mgi_id)))
  view <- paste(
    "Allele.primaryIdentifier",
    "Allele.symbol",
    "Allele.name",
    "Allele.alleleType",
    "Allele.drivenBy.feature.symbol",
    "Allele.feature.primaryIdentifier",
    "Allele.feature.locations.start",
    "Allele.feature.locations.end",
    "Allele.feature.locations.strand",
    "Allele.feature.locations.locatedOn.primaryIdentifier",
    "Allele.feature.locations.assembly",
    sep = " "
  )
  sprintf(
    paste0(
      "<query name=\"\" model=\"genomic\" view=\"%s\">",
      "<constraint path=\"Allele.primaryIdentifier\" op=\"=\" value=\"%s\"/>",
      "</query>"
    ),
    view,
    mid
  )
}

allele_overlap_query_xml <- function(mgi_id, limit = 500L) {
  mid <- xml_escape(trimws(as.character(mgi_id)))
  view <- paste(
    "Allele.primaryIdentifier",
    "Allele.symbol",
    "Allele.feature.overlappingFeatures.primaryIdentifier",
    "Allele.feature.overlappingFeatures.symbol",
    "Allele.feature.overlappingFeatures.name",
    sep = " "
  )
  sprintf(
    paste0(
      "<query name=\"\" model=\"genomic\" view=\"%s\">",
      "<constraint path=\"Allele.primaryIdentifier\" op=\"=\" value=\"%s\"/>",
      "</query>"
    ),
    view,
    mid
  )
}

#' Parse MouseMine tabular `results` into a data.frame
#'
#' MouseMine JSON with `simplifyVector = TRUE` collapses a single result row into a
#' flat atomic vector; this normalizes all shapes (vector, matrix, list-of-rows).
#' @noRd
mousemine_results_df <- function(body, col_names) {
  ncn <- length(col_names)
  empty <- data.frame(matrix(ncol = ncn, nrow = 0), stringsAsFactors = FALSE)
  colnames(empty) <- col_names
  if (is.null(body) || !isTRUE(body$wasSuccessful) || is.null(body$results) || !length(body$results)) {
    return(empty)
  }
  res <- body$results
  mat <- NULL
  if (is.matrix(res)) {
    mat <- res
  } else if (is.data.frame(res)) {
    mat <- as.matrix(res)
  } else if (is.list(res) && length(res) > 0L && is.list(res[[1L]])) {
    mat <- do.call(
      rbind,
      lapply(res, function(row) {
        v <- unlist(row, use.names = FALSE)
        if (length(v) < ncn) {
          v <- c(v, rep(NA, ncn - length(v)))
        }
        if (length(v) > ncn) {
          v <- v[seq_len(ncn)]
        }
        v
      })
    )
  } else {
    v <- unlist(res, use.names = FALSE)
    if (length(v) == ncn) {
      mat <- matrix(v, nrow = 1L)
    } else if (length(v) > ncn && length(v) %% ncn == 0L) {
      mat <- matrix(v, ncol = ncn, byrow = TRUE)
    }
  }
  if (is.null(mat) || !nrow(mat)) {
    return(empty)
  }
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(out) <- col_names
  out
}

#' Fetch transgene / allele integration coordinates from MouseMine
#'
#' @param mgi_id MGI allele primary identifier (e.g. `MGI:2176173`).
#' @return Data frame of integration intervals and metadata; zero rows if none.
#' @export
fetch_allele_integration_sites <- function(mgi_id) {
  mgi_id <- trimws(as.character(mgi_id))
  if (!nzchar(mgi_id) || is.na(mgi_id)) {
    return(data.frame(
      mgi_id = character(),
      symbol = character(),
      name = character(),
      allele_type = character(),
      driver_gene_symbol = character(),
      feature_id = character(),
      start = integer(),
      end = integer(),
      strand = character(),
      chromosome = character(),
      assembly = character(),
      stringsAsFactors = FALSE
    ))
  }
  body <- mousemine_query_json(allele_genomic_query_xml(mgi_id), limit = 50L)
  cn <- c(
    "mgi_id", "symbol", "name", "allele_type", "driver_gene_symbol",
    "feature_id", "start", "end", "strand", "chromosome", "assembly"
  )
  out <- mousemine_results_df(body, cn)
  if (nrow(out)) {
    out$start <- suppressWarnings(as.integer(out$start))
    out$end <- suppressWarnings(as.integer(out$end))
    chr <- as.character(out$chromosome)
    out <- out[
      !is.na(out$start) &
        !is.na(out$end) &
        !is.na(chr) &
        nzchar(chr),
      ,
      drop = FALSE
    ]
    rownames(out) <- NULL
  }
  out
}

#' Fetch genes and features overlapping the allele insertion feature
#'
#' @param mgi_id MGI allele primary identifier.
#' @return Data frame with `feature_id`, `symbol`, `name` (deduplicated).
#' @export
fetch_allele_overlapping_features <- function(mgi_id) {
  mgi_id <- trimws(as.character(mgi_id))
  if (!nzchar(mgi_id) || is.na(mgi_id)) {
    return(data.frame(
      mgi_id = character(),
      symbol = character(),
      feature_id = character(),
      name = character(),
      stringsAsFactors = FALSE
    ))
  }
  body <- mousemine_query_json(allele_overlap_query_xml(mgi_id), limit = 500L)
  cn <- c("mgi_id", "allele_symbol", "feature_id", "symbol", "name")
  out <- mousemine_results_df(body, cn)
  if (!nrow(out)) {
    return(data.frame(
      mgi_id = character(),
      symbol = character(),
      feature_id = character(),
      name = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- out[, c("mgi_id", "symbol", "feature_id", "name"), drop = FALSE]
  out <- out[!is.na(out$symbol) & nzchar(out$symbol), , drop = FALSE]
  out <- unique(out)
  rownames(out) <- NULL
  out
}

#' Resolve an allele MGI id from symbol via MouseMine
#'
#' @param symbol Allele symbol (exact match).
#' @return First matching `MGI:#####` id or `NA_character_`.
#' @export
resolve_allele_mgi_id <- function(symbol) {
  sym <- trimws(as.character(symbol))
  if (!nzchar(sym) || is.na(sym)) {
    return(NA_character_)
  }
  pat <- xml_escape(sym)
  view <- "Allele.primaryIdentifier Allele.symbol"
  query_symbol <- function(op) {
    q <- sprintf(
      "<query name=\"\" model=\"genomic\" view=\"%s\" sortOrder=\"Allele.symbol asc\"><constraint path=\"Allele.symbol\" op=\"%s\" value=\"%s\"/></query>",
      view,
      op,
      pat
    )
    mousemine_results_df(mousemine_query_json(q, limit = 20L), c("mgi_id", "symbol"))
  }
  out <- query_symbol("=")
  if (!nrow(out)) {
    out <- query_symbol("CONTAINS")
  }
  if (!nrow(out)) {
    return(NA_character_)
  }
  trimws(out$mgi_id[[1L]])
}
