.gp_equal_bins <- function(x, nbin) {
  n <- length(x); nbin <- max(2L, min(as.integer(nbin), n))
  ord <- order(as.numeric(x) + seq_along(x) * .Machine$double.eps,
    decreasing = FALSE, na.last = TRUE)
  out <- integer(n); out[ord] <- ceiling(seq_len(n) * nbin / n)
  names(out) <- names(x); out
}

.run_addmodulescore <- function(readers, writers, params, context) {
  .gp_require("Matrix", "module scoring")
  cells <- as.character(context$effective_cells)
  groups <- .gp_groups(params$features, params$name %||% "Module")
  probe <- readers$matrix(cells = cells)
  universe <- as.character(probe$feature_id)
  groups <- lapply(groups, function(x) intersect(x, universe))
  groups <- groups[lengths(groups) > 0L]
  .gp_assert(length(groups) > 0L, "None of the module features exist in the matrix.")
  pool <- intersect(unique(as.character(params$pool %||% universe)), universe)
  .gp_assert(length(pool) >= 2L, "Control pool is too small after filtering.")
  pool_matrix <- readers$matrix(features = pool, cells = cells)$values
  avg <- Matrix::rowMeans(pool_matrix); names(avg) <- pool
  bins <- .gp_equal_bins(avg, params$nbin %||% 24L)
  bin_map <- split(names(bins), bins)
  seed <- params$seed %||% 1L
  if (!is.null(seed)) set.seed(as.integer(seed))
  keys <- names(groups); values <- vector("list", length(keys))
  for (i in seq_along(keys)) {
    genes <- groups[[i]]
    controls <- unique(unlist(lapply(genes, function(g) {
      candidates <- setdiff(bin_map[[as.character(bins[[g]])]], genes)
      if (!length(candidates)) return(character())
      sample(candidates, min(as.integer(params$ctrl %||% 100L), length(candidates)))
    }), use.names = FALSE))
    .gp_assert(length(controls) > 0L,
      paste0("No control genes available for module '", keys[[i]], "'."))
    mod <- readers$matrix(features = genes, cells = cells)$values
    ctl <- readers$matrix(features = controls, cells = cells)$values
    values[[i]] <- as.numeric(Matrix::colMeans(mod) - Matrix::colMeans(ctl))
  }
  out <- scworkbench::wb_prim_cell_vector_family(rep(keys, each = length(cells)),
    rep(cells, times = length(keys)), unlist(values, use.names = FALSE),
    value_type = "continuous")
  .gp_write(writers, "module_scores", out)
  invisible(list(status = "ok"))
}

.gp_ucell_score <- function(x, signatures, max_rank, ncores, backend) {
  .gp_require("UCell", "UCell signature scoring")
  .gp_require("BiocParallel", "UCell parallel backend selection")
  x <- methods::as(x, "dgCMatrix")
  bp <- if (identical(backend, "serial")) BiocParallel::SerialParam() else NULL
  scored <- UCell::ScoreSignatures_UCell(x, features = signatures,
    maxRank = as.integer(max_rank), ncores = as.integer(ncores), BPPARAM = bp,
    name = "")
  if (is.matrix(scored)) scored <- as.data.frame(scored, stringsAsFactors = FALSE)
  scored
}

.gp_ucell_summary <- function(meta, family, matrix_id, policy) {
  rows <- lapply(names(meta), function(key) {
    vals <- family$value[family$key == key]
    q <- stats::quantile(vals, c(.25, .5, .75), na.rm = TRUE,
      names = FALSE, type = 7)
    m <- meta[[key]]
    data.frame(entity_id = key, signature = key,
      n_features_input = length(m$input), n_features_matched = length(m$matched),
      n_features_missing = length(m$missing),
      matched_features = paste(m$matched, collapse = ";"),
      missing_features = paste(m$missing, collapse = ";"),
      score_min = min(vals, na.rm = TRUE), score_q25 = q[[1]],
      score_median = q[[2]], score_mean = mean(vals, na.rm = TRUE),
      score_q75 = q[[3]], score_max = max(vals, na.rm = TRUE),
      n_cells = length(vals), stringsAsFactors = FALSE)
  })
  scworkbench::wb_prim_obs_table(do.call(rbind, rows),
    table_type = "ucell_signature_summary", entity_kind = "signature",
    meta = list(matrix_id = matrix_id, missing_features_policy = policy))
}

.run_ucell <- function(readers, writers, params, context) {
  .gp_require("Matrix", "UCell signature scoring")
  capability <- context$typed_slice_reads$matrix_chunks
  .gp_assert(is.list(capability), "UCell requires its declared typed-slice grant.")
  cells <- as.character(capability$cells)
  universe <- as.character(capability$features)
  groups <- .gp_groups(params$features, params$name %||% "Signature")
  policy <- params$missing_features_policy %||% "drop"
  meta <- lapply(groups, function(x) {
    input <- unique(as.character(x)); matched <- intersect(input, universe)
    missing <- setdiff(input, universe)
    if (identical(policy, "error") && length(missing))
      stop("UCell signature contains missing features.", call. = FALSE)
    .gp_assert(length(matched) > 0L, "UCell signature has no matched features.")
    list(input = input, matched = matched, missing = missing)
  })
  groups_use <- lapply(meta, `[[`, "matched")
  chunk_size <- as.integer(params$chunk_size %||% 2000L)
  chunks <- scworkbench::wb_partition_axis(cells, chunk_size)
  matrix_id <- NULL
  pieces <- lapply(chunks, function(chunk_cells) {
    primitive <- readers$matrix_chunks(features = universe, cells = chunk_cells)
    if (is.null(matrix_id)) matrix_id <<- primitive$matrix_id
    score <- .gp_ucell_score(primitive$values, groups_use,
      params$max_rank %||% 1500L, params$ncores %||% 1L,
      params$parallel_backend %||% "auto")
    returned <- as.character(rownames(score))
    .gp_assert(!anyDuplicated(returned) && setequal(returned, chunk_cells),
      "UCell returned invalid cell identities.")
    score[match(chunk_cells, returned), names(groups), drop = FALSE]
  })
  merged <- do.call(rbind, pieces)
  merged <- merged[match(cells, rownames(merged)), names(groups), drop = FALSE]
  family <- scworkbench::wb_prim_cell_vector_family(
    rep(names(groups), each = length(cells)), rep(cells, times = length(groups)),
    unlist(lapply(names(groups), function(key) as.numeric(merged[[key]])),
      use.names = FALSE), value_type = "continuous")
  .gp_write(writers, "signature_scores", family)
  .gp_write(writers, "signature_summary", .gp_ucell_summary(meta, family,
    matrix_id, policy))
  invisible(list(status = "ok"))
}
