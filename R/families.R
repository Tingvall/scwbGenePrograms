.run_family_reduce <- function(readers, writers, params, context) {
  sm <- .gp_signal(readers$signal(), params$signal_keys,
    context$effective_cells)
  mat <- sm$matrix; cells <- sm$cells
  unknown <- as.character(params$unknown_label %||% "unknown")[1]
  ties <- params$ties %||% "unknown"; tol <- as.numeric(params$tie_tol %||% 1e-8)
  label <- rep(unknown, length(cells)); top <- second <- margin <- rep(NA_real_, length(cells))
  reason <- rep("none", length(cells))
  for (i in seq_along(cells)) {
    vals <- mat[i, ]; ok <- is.finite(vals)
    if (!any(ok)) next
    ordered <- order(vals[ok], decreasing = TRUE, na.last = TRUE)
    vv <- vals[ok][ordered]; kk <- names(vals[ok])[ordered]
    top[[i]] <- vv[[1]]; second[[i]] <- if (length(vv) > 1L) vv[[2]] else NA_real_
    margin[[i]] <- if (length(vv) > 1L) vv[[1]] - vv[[2]] else NA_real_
    tied <- abs(vv - vv[[1]]) <= tol
    if (sum(tied) > 1L && identical(ties, "unknown")) { reason[[i]] <- "tie"; next }
    winner <- if (sum(tied) > 1L) sort(kk[tied], method = "radix")[[1]] else kk[[1]]
    if (!is.null(params$min_score) && vv[[1]] < as.numeric(params$min_score)) {
      reason[[i]] <- "min_score"; next
    }
    if (!is.null(params$min_margin) && is.finite(margin[[i]]) &&
        margin[[i]] < as.numeric(params$min_margin)) {
      reason[[i]] <- "min_margin"; next
    }
    label[[i]] <- winner
  }
  primitives <- list(
    label = scworkbench::wb_prim_cell_vector(stats::setNames(label, cells), "categorical"),
    top_score = scworkbench::wb_prim_cell_vector(stats::setNames(top, cells), "continuous"),
    runner_up_score = scworkbench::wb_prim_cell_vector(stats::setNames(second, cells), "continuous"),
    score_margin = scworkbench::wb_prim_cell_vector(stats::setNames(margin, cells), "continuous"),
    abstain_reason = scworkbench::wb_prim_cell_vector(stats::setNames(reason, cells), "categorical"))
  decision <- data.frame(cell_id = cells, label = label, top_score = top,
    runner_up_score = second, score_margin = margin, abstain_reason = reason,
    stringsAsFactors = FALSE)
  summary <- do.call(rbind, lapply(split(decision, decision$label), function(x)
    data.frame(label = x$label[[1]], n_cells = nrow(x), pct_cells = nrow(x) / length(cells),
      mean_top_score = mean(x$top_score, na.rm = TRUE),
      mean_margin = mean(x$score_margin, na.rm = TRUE),
      n_abstained = sum(x$abstain_reason != "none"), stringsAsFactors = FALSE)))
  rownames(summary) <- NULL
  summary <- summary[order(summary$n_cells, decreasing = TRUE, summary$label), , drop = FALSE]
  meta <- list(rule = params$rule %||% "argmax", unknown_label = unknown,
    signal_keys = sm$keys)
  primitives$decision_table <- scworkbench::wb_prim_obs_table(decision,
    "cell_vector_family_reduce_decision_table", "cell", meta)
  primitives$label_summary <- scworkbench::wb_prim_obs_table(summary,
    "cell_vector_family_reduce_summary", "label", meta)
  for (name in names(primitives)) .gp_write(writers, name, primitives[[name]])
  invisible(list(status = "ok"))
}

.run_family_calibrate <- function(readers, writers, params, context) {
  sm <- .gp_signal(readers$signal(), params$signal_keys, context$effective_cells)
  method <- params$method %||% "percentile"; ties <- params$ties_method %||% "average"
  out <- sm$matrix; summaries <- vector("list", length(sm$keys))
  for (j in seq_along(sm$keys)) {
    raw <- sm$matrix[, j]; ok <- is.finite(raw); fitted <- raw[ok]
    vals <- rep(NA_real_, length(raw)); center <- scale <- NA_real_
    if (length(fitted)) {
      if (identical(method, "percentile")) {
        vals[ok] <- if (length(fitted) == 1L) 1 else
          (rank(fitted, ties.method = ties) - 1) / (length(fitted) - 1)
      } else if (identical(method, "zscore")) {
        center <- if (isTRUE(params$center %||% TRUE)) mean(fitted) else 0
        scale <- if (isTRUE(params$scale %||% TRUE)) stats::sd(fitted) else 1
        if (!is.finite(scale) || scale == 0) scale <- 1
        vals[ok] <- (fitted - center) / scale
      } else {
        center <- stats::median(fitted); scale <- stats::mad(fitted, center, constant = 1)
        if (!is.finite(scale) || scale == 0) scale <- 1
        vals[ok] <- (fitted - center) / scale
      }
    }
    out[, j] <- vals
    rq <- stats::quantile(raw, c(.25, .5, .75), na.rm = TRUE, names = FALSE)
    cq <- stats::quantile(vals, c(.25, .5, .75), na.rm = TRUE, names = FALSE)
    summaries[[j]] <- data.frame(entity_id = sm$keys[[j]], key = sm$keys[[j]],
      n_cells = length(raw), n_finite = sum(ok), raw_min = min(raw, na.rm = TRUE),
      raw_q25 = rq[[1]], raw_median = rq[[2]], raw_mean = mean(raw, na.rm = TRUE),
      raw_q75 = rq[[3]], raw_max = max(raw, na.rm = TRUE),
      calibrated_min = min(vals, na.rm = TRUE), calibrated_q25 = cq[[1]],
      calibrated_median = cq[[2]], calibrated_mean = mean(vals, na.rm = TRUE),
      calibrated_q75 = cq[[3]], calibrated_max = max(vals, na.rm = TRUE),
      calibration_center = center, calibration_scale = scale, stringsAsFactors = FALSE)
  }
  .gp_write(writers, "calibrated_signal", .gp_family(out))
  .gp_write(writers, "calibration_summary", scworkbench::wb_prim_obs_table(
    do.call(rbind, summaries), "cell_vector_family_calibration_summary", "signal_key",
    list(method = method, signal_keys = sm$keys)))
  invisible(list(status = "ok"))
}

.gp_distribution_row <- function(values, id, prefix = list()) {
  finite <- values[is.finite(values)]
  q <- if (length(finite)) stats::quantile(finite, c(.25, .5, .75),
    names = FALSE) else rep(NA_real_, 3L)
  c(prefix, list(entity_id = id, n_finite = length(finite),
    pct_finite = length(finite) / max(1L, length(values)),
    score_min = if (length(finite)) min(finite) else NA_real_, score_q25 = q[[1]],
    score_median = q[[2]], score_mean = if (length(finite)) mean(finite) else NA_real_,
    score_q75 = q[[3]], score_max = if (length(finite)) max(finite) else NA_real_,
    score_sd = if (length(finite) > 1L) stats::sd(finite) else NA_real_,
    dynamic_range = if (length(finite)) diff(range(finite)) else NA_real_))
}

.run_family_qc <- function(readers, writers, params, context) {
  sm <- .gp_signal(readers$signal(), params$signal_keys, context$effective_cells)
  mat <- sm$matrix; tol <- as.numeric(params$tie_tol %||% 1e-8)
  methods <- as.character(unlist(params$correlation_methods %||%
    c("pearson", "spearman"), use.names = FALSE))
  winner <- tied <- stats::setNames(integer(length(sm$keys)), sm$keys)
  for (i in seq_len(nrow(mat))) {
    vals <- mat[i, ]; vals <- vals[is.finite(vals)]; if (!length(vals)) next
    keys <- names(vals)[abs(vals - max(vals)) <= tol]
    if (length(keys) == 1L) winner[keys] <- winner[keys] + 1L else tied[keys] <- tied[keys] + 1L
  }
  key_rows <- lapply(sm$keys, function(key) {
    x <- .gp_distribution_row(mat[, key], key, list(key = key, n_cells_total = nrow(mat)))
    x$top1_count <- as.integer(winner[[key]]); x$top1_pct <- winner[[key]] / nrow(mat)
    x$tied_top_count <- as.integer(tied[[key]])
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    x[, c("entity_id", "key", "n_cells_total", setdiff(names(x),
      c("entity_id", "key", "n_cells_total"))), drop = FALSE]
  })
  pair_rows <- list(); p <- 0L
  if (length(sm$keys) > 1L) for (i in seq_len(length(sm$keys) - 1L))
    for (j in seq.int(i + 1L, length(sm$keys))) {
      a <- mat[, i]; b <- mat[, j]; keep <- is.finite(a) & is.finite(b); p <- p + 1L
      pair_rows[[p]] <- data.frame(key_a = sm$keys[[i]], key_b = sm$keys[[j]],
        n_overlap = sum(keep),
        pearson_cor = if (sum(keep) > 1L && "pearson" %in% methods) stats::cor(a[keep], b[keep]) else NA_real_,
        spearman_cor = if (sum(keep) > 1L && "spearman" %in% methods) stats::cor(a[keep], b[keep], method = "spearman") else NA_real_,
        stringsAsFactors = FALSE)
    }
  pair <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame(
    key_a = character(), key_b = character(), n_overlap = integer(),
    pearson_cor = numeric(), spearman_cor = numeric(), stringsAsFactors = FALSE)
  group <- data.frame(entity_id = character(), group_id = character(), key = character(),
    n_cells_group = integer(), n_finite = integer(), pct_finite = numeric(),
    score_min = numeric(), score_q25 = numeric(), score_median = numeric(),
    score_mean = numeric(), score_q75 = numeric(), score_max = numeric(),
    score_sd = numeric(), dynamic_range = numeric(), top1_count = integer(),
    top1_pct = numeric(), tied_top_count = integer(), stringsAsFactors = FALSE)
  grouping_levels <- character()
  if ("grouping" %in% names(readers)) {
    g <- readers$grouping(); values <- as.character(g$value); names(values) <- g$cell_id
    values <- values[sm$cells]; grouping_levels <- unique(values[!is.na(values) & nzchar(values)])
    rows <- list(); n <- 0L
    for (level in grouping_levels) for (key in sm$keys) {
      ids <- names(values)[values == level]; vals <- mat[ids, key]
      d <- .gp_distribution_row(vals, paste0(level, "::", key),
        list(group_id = level, key = key, n_cells_group = length(vals)))
      best <- apply(mat[ids, , drop = FALSE], 1L, function(v) {
        if (!any(is.finite(v))) return(character())
        names(v)[abs(v - max(v, na.rm = TRUE)) <= tol]
      })
      d$top1_count <- sum(vapply(best, function(x) length(x) == 1L && key %in% x, logical(1)))
      d$top1_pct <- d$top1_count / max(1L, length(ids))
      d$tied_top_count <- sum(vapply(best, function(x) length(x) > 1L && key %in% x, logical(1)))
      n <- n + 1L; rows[[n]] <- as.data.frame(d, stringsAsFactors = FALSE)
    }
    if (length(rows)) {
      group <- do.call(rbind, rows)
      group <- group[, c("entity_id", "group_id", "key", "n_cells_group",
        setdiff(names(group), c("entity_id", "group_id", "key",
          "n_cells_group"))), drop = FALSE]
    }
  }
  .gp_write(writers, "key_summary", scworkbench::wb_prim_obs_table(do.call(rbind, key_rows),
    "cell_vector_family_qc_key_summary", "signal_key",
    list(signal_keys = sm$keys, tie_tol = tol)))
  .gp_write(writers, "pairwise_summary", scworkbench::wb_prim_obs_table(pair,
    "cell_vector_family_qc_pairwise_summary", "signal_key_pair",
    list(signal_keys = sm$keys, correlation_methods = methods)))
  .gp_write(writers, "group_summary", scworkbench::wb_prim_obs_table(group,
    "cell_vector_family_qc_group_summary", "signal_key_group",
    list(signal_keys = sm$keys, grouping_levels = grouping_levels,
      has_grouping = "grouping" %in% names(readers))))
  invisible(list(status = "ok"))
}

.run_family_graph_smooth <- function(readers, writers, params, context) {
  .gp_require("RANN", "cell-vector-family graph smoothing")
  sm <- .gp_signal(readers$signal(), params$signal_keys, context$effective_cells)
  space <- readers$space(); coords <- as.matrix(space$values)
  rownames(coords) <- as.character(space$cell_id)
  dims <- params$dims
  if (!is.null(dims)) coords <- coords[, as.integer(dims), drop = FALSE]
  coords <- coords[match(sm$cells, rownames(coords)), , drop = FALSE]
  .gp_assert(!anyNA(match(sm$cells, rownames(coords))), "Space is missing signal cells.")
  finite <- stats::complete.cases(coords); coords <- coords[finite, , drop = FALSE]
  raw <- sm$matrix[finite, , drop = FALSE]
  .gp_assert(nrow(coords) > 1L, "Need at least two cells with finite coordinates.")
  k <- min(as.integer(params$k %||% 30L), nrow(coords) - !isTRUE(params$include_self))
  alpha <- as.numeric(params$alpha %||% .5); qsize <- as.integer(params$query_chunk_size %||% 512L)
  include_self <- isTRUE(params$include_self); policy <- params$na_policy %||% "ignore"
  out <- matrix(NA_real_, nrow(raw), ncol(raw), dimnames = dimnames(raw)); counts <- integer(nrow(raw))
  chunks <- split(seq_len(nrow(raw)), ceiling(seq_len(nrow(raw)) / qsize))
  for (chunk in chunks) {
    nn <- RANN::nn2(coords, coords[chunk, , drop = FALSE],
      k = min(nrow(coords), k + if (include_self) 0L else 1L))$nn.idx
    for (j in seq_along(chunk)) {
      i <- chunk[[j]]; neighbors <- as.integer(nn[j, ])
      if (!include_self) neighbors <- neighbors[neighbors != i]
      neighbors <- head(neighbors, k); counts[[i]] <- length(neighbors)
      values <- raw[neighbors, , drop = FALSE]
      mean_neighbor <- if (identical(policy, "zero")) {
        values[!is.finite(values)] <- 0; colMeans(values)
      } else if (identical(policy, "propagate")) {
        x <- colMeans(values); x[colSums(!is.finite(values)) > 0L] <- NA_real_; x
      } else suppressWarnings(colMeans(values, na.rm = TRUE))
      base <- raw[i, ]
      if (identical(policy, "zero")) { base[!is.finite(base)] <- 0; mean_neighbor[!is.finite(mean_neighbor)] <- 0 }
      blended <- alpha * base + (1 - alpha) * mean_neighbor
      if (!identical(policy, "zero")) {
        blended[is.finite(base) & !is.finite(mean_neighbor)] <- base[is.finite(base) & !is.finite(mean_neighbor)]
        blended[!is.finite(base) & is.finite(mean_neighbor)] <- mean_neighbor[!is.finite(base) & is.finite(mean_neighbor)]
      }
      out[i, ] <- blended
    }
  }
  ancestry <- if ("source_ancestry" %in% names(readers)) readers$source_ancestry() else list()
  source_artifact <- params$source_signal_artifact_id %||% ancestry$artifact_id
  source_run <- params$source_run_id %||% ancestry$source_run_id
  source_snapshot <- params$source_snapshot_id %||% ancestry$source_snapshot_id
  meta <- list(source_signal_artifact_id = .gp_scalar(source_artifact, NA_character_),
    source_run_id = .gp_scalar(source_run, NA_character_),
    source_snapshot_id = .gp_scalar(source_snapshot, NA_character_),
    source_space_id = space$space_id, neighbor_artifact_id = NA_character_,
    selected_keys = sm$keys, selected_cells = NA_character_,
    n_cells_input_signal = nrow(sm$matrix), n_cells_smoothed = nrow(out),
    dims = if (is.null(dims)) NULL else as.integer(dims), k = as.integer(k), alpha = alpha,
    smoothing_method = "knn_neighbor_mean_blend",
    smoothing_formula = "smoothed_score = alpha * raw_score + (1 - alpha) * neighbor_mean_score",
    distance_metric = "euclidean", query_chunk_size = qsize,
    include_self = include_self, na_policy = policy,
    scworkbench_version = as.character(utils::packageVersion("scworkbench")))
  summary <- do.call(rbind, lapply(sm$keys, function(key) data.frame(key = key,
    n_cells = nrow(raw), n_cells_smoothed = sum(is.finite(out[, key])),
    n_missing_before = sum(!is.finite(raw[, key])), n_missing_after = sum(!is.finite(out[, key])),
    raw_mean = mean(raw[, key], na.rm = TRUE), smoothed_mean = mean(out[, key], na.rm = TRUE),
    raw_sd = stats::sd(raw[, key], na.rm = TRUE), smoothed_sd = stats::sd(out[, key], na.rm = TRUE),
    avg_neighbor_count = mean(counts), k = as.integer(k), alpha = alpha,
    smoothing_method = "knn_neighbor_mean_blend", stringsAsFactors = FALSE)))
  .gp_write(writers, "smoothed_signal", .gp_family(out,
    as.character(readers$signal()$value_type %||% "continuous")[1]))
  .gp_write(writers, "smoothing_summary", scworkbench::wb_prim_obs_table(summary,
    "cell_vector_family_graph_smooth_summary", "signal_key", meta))
  invisible(list(status = "ok"))
}
