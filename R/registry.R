.gp_recipe_raw_n <- function(value) {
  value <- .gp_scalar(value)
  if (!nzchar(value)) return(NA_integer_)
  parsed <- try(jsonlite::fromJSON(value, simplifyVector = FALSE), silent = TRUE)
  if (inherits(parsed, "try-error")) return(NA_integer_)
  as.integer(parsed$n_genes_raw %||% parsed$raw_gene_count %||% NA_integer_)[1]
}

.run_feature_set_overlap_qc <- function(readers, writers, params, context) {
  registry <- readers$feature_sets(); sets <- .gp_registry_sets(registry)
  .gp_assert(length(sets) > 0L, "No feature sets match requested filters.")
  min_size <- as.integer(params$min_set_size %||% 1L)
  max_size <- as.numeric(params$max_set_size %||% Inf)
  min_fraction <- as.numeric(params$min_detected_fraction %||% 0)
  universe <- as.character(unlist(registry$feature_universe, use.names = FALSE))
  rows <- vector("list", length(sets)); signatures <- character(length(sets))
  for (i in seq_along(sets)) {
    item <- sets[[i]]; mapped <- intersect(item$member_ids, universe)
    raw_n <- .gp_recipe_raw_n(item$recipe_json)
    if (is.na(raw_n)) raw_n <- as.integer(item$n_features)
    fraction <- if (raw_n > 0L) length(mapped) / raw_n else 0
    reasons <- character()
    if (length(mapped) < min_size) reasons <- c(reasons, "below_min_set_size")
    if (length(mapped) > max_size) reasons <- c(reasons, "above_max_set_size")
    if (fraction < min_fraction) reasons <- c(reasons, "below_min_detected_fraction")
    signatures[[i]] <- paste(sort(unique(mapped), method = "radix"), collapse = "\r")
    rows[[i]] <- data.frame(namespace = .gp_scalar(item$namespace),
      feature_set = item$feature_set, set_name = item$feature_set,
      matrix_id = registry$matrix_id, n_genes_raw = raw_n,
      n_genes_mapped_to_bundle = as.integer(item$n_features),
      n_genes_detected_in_matrix = length(mapped), detected_fraction = fraction,
      screening_status = if (length(reasons)) "exclude" else "include",
      reason = if (length(reasons)) paste(reasons, collapse = ";") else "ok",
      source_id = .gp_scalar(item$source_id, NA_character_),
      source_version = .gp_scalar(item$source_version, NA_character_),
      stringsAsFactors = FALSE)
  }
  qc <- do.call(rbind, rows)
  if (isTRUE(params$deduplicate_sets %||% TRUE)) {
    dup <- duplicated(signatures); qc$duplicate_signature <- dup | duplicated(signatures, fromLast = TRUE)
    qc$duplicate_group <- match(signatures, unique(signatures))
    qc$screening_status[dup] <- "exclude"
    qc$reason[dup] <- ifelse(qc$reason[dup] == "ok", "duplicate_set_after_dedup",
      paste(qc$reason[dup], "duplicate_set_after_dedup", sep = ";"))
  } else { qc$duplicate_signature <- FALSE; qc$duplicate_group <- seq_len(nrow(qc)) }
  namespaces <- sort(unique(qc$namespace), method = "radix")
  summary <- do.call(rbind, lapply(namespaces, function(namespace) {
    sub <- qc[qc$namespace == namespace, , drop = FALSE]
    data.frame(namespace = namespace, matrix_id = registry$matrix_id,
      n_sets = nrow(sub), n_include = sum(sub$screening_status == "include"),
      n_exclude = sum(sub$screening_status == "exclude"),
      median_detected_fraction = stats::median(sub$detected_fraction, na.rm = TRUE),
      median_detected_genes = stats::median(sub$n_genes_detected_in_matrix, na.rm = TRUE),
      min_set_size = min_size, max_set_size = if (is.infinite(max_size)) NA_real_ else max_size,
      min_detected_fraction = min_fraction,
      deduplicate_sets = isTRUE(params$deduplicate_sets %||% TRUE), stringsAsFactors = FALSE)
  }))
  meta <- list(matrix_id = registry$matrix_id, namespaces = namespaces,
    op = "feature_set_overlap_qc", min_set_size = min_size,
    max_set_size = if (is.infinite(max_size)) "Inf" else max_size,
    min_detected_fraction = min_fraction,
    deduplicate_sets = isTRUE(params$deduplicate_sets %||% TRUE))
  .gp_write(writers, "overlap_qc", scworkbench::wb_prim_obs_table(qc,
    "feature_set_overlap_qc", "feature_set", meta))
  .gp_write(writers, "namespace_summary", scworkbench::wb_prim_obs_table(summary,
    "feature_set_overlap_qc_summary", "namespace", meta))
  invisible(list(status = "ok"))
}

.run_pathway_enrich <- function(readers, writers, params, context) {
  matrix <- readers$matrix(); input <- readers$feature_table()$table
  registry <- readers$feature_sets(); sets <- .gp_registry_sets(registry)
  feature_col <- params$feature_col %||% "feature_id"
  .gp_assert(feature_col %in% names(input), "feature_col not found.")
  q <- input; q[[feature_col]] <- as.character(q[[feature_col]])
  q <- q[nzchar(q[[feature_col]]), , drop = FALSE]
  pcol <- params$p_adj_col %||% "p_val_adj"
  if (!is.null(pcol) && nzchar(pcol) && pcol %in% names(q)) {
    pv <- as.numeric(q[[pcol]]); q <- q[!is.na(pv) & pv <= as.numeric(params$p_adj_max %||% .05), , drop = FALSE]
  }
  lcol <- params$lfc_col %||% "avg_log1p_fc"; direction <- params$direction %||% "any"
  if (!is.null(lcol) && nzchar(lcol) && lcol %in% names(q)) {
    value <- as.numeric(q[[lcol]])
    keep <- if (identical(direction, "up")) value >= as.numeric(params$lfc_min %||% .25) else
      if (identical(direction, "down")) value <= -as.numeric(params$lfc_min %||% .25) else
        abs(value) >= as.numeric(params$lfc_min_abs %||% 0)
    q <- q[!is.na(keep) & keep, , drop = FALSE]
  }
  universe <- as.character(unlist(registry$feature_universe, use.names = FALSE))
  query <- intersect(unique(q[[feature_col]]), universe)
  .gp_assert(length(query) > 0L, "No query features overlap matrix universe.")
  rows <- list(); n <- 0L
  for (item in sets) {
    members <- unique(intersect(item$member_ids, universe)); size <- length(members)
    if (size < as.integer(params$min_set_size %||% 5L) || size > as.integer(params$max_set_size %||% 5000L)) next
    overlap <- intersect(query, members); k <- length(overlap); n <- n + 1L
    rows[[n]] <- data.frame(feature_id = item$feature_set, pathway = item$feature_set,
      namespace = params$namespace, matrix_id = matrix$matrix_id,
      set_size = size, query_size = length(query), universe_size = length(universe),
      overlap = k, enrichment_ratio = if (k) (k / length(query)) / (size / length(universe)) else 0,
      p_val = if (k) stats::phyper(k - 1L, size, length(universe) - size,
        length(query), lower.tail = FALSE) else 1,
      overlap_genes = paste(head(overlap, 200L), collapse = ","), stringsAsFactors = FALSE)
  }
  if (length(rows)) {
    out <- do.call(rbind, rows); out$p_val_adj <- stats::p.adjust(out$p_val,
      method = params$p_adjust_method %||% "BH")
    out <- out[order(out$p_val_adj, -out$enrichment_ratio, -out$overlap), , drop = FALSE]
    rownames(out) <- NULL
  } else out <- data.frame(feature_id = character(), pathway = character(), namespace = character(),
    matrix_id = character(), set_size = integer(), query_size = integer(), universe_size = integer(),
    overlap = integer(), enrichment_ratio = numeric(), p_val = numeric(), p_val_adj = numeric(),
    overlap_genes = character(), stringsAsFactors = FALSE)
  .gp_write(writers, "pathways", scworkbench::wb_prim_obs_table(out,
    "pathway_enrich", "pathway", list(matrix_id = matrix$matrix_id,
      namespace = params$namespace, query_size = length(query), universe_size = length(universe))))
  invisible(list(status = "ok"))
}

