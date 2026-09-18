.gp_winner_token <- function(x) {
  out <- gsub("[^A-Za-z0-9_]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out); out[!nzchar(out)] <- "unknown"; out
}

.run_cell_vector_winners <- function(readers, writers, params, context) {
  sm <- .gp_signal(readers$signal(), params$signal_keys, context$effective_cells)
  mat <- sm$matrix; cells <- sm$cells; keys <- sm$keys
  tie_method <- params$tie_method %||% "first"
  .gp_assert(tie_method %in% c("first", "last"), "Unsupported tie_method.")
  ok <- rowSums(is.finite(mat)) > 0L; winner_i <- rep(NA_integer_, nrow(mat))
  if (any(ok)) {
    x <- mat[ok, , drop = FALSE]; x[!is.finite(x)] <- -Inf
    if (identical(tie_method, "last")) {
      rev_i <- max.col(x[, rev(seq_len(ncol(x))), drop = FALSE], ties.method = "first")
      winner_i[ok] <- ncol(x) - rev_i + 1L
    } else winner_i[ok] <- max.col(x, ties.method = "first")
  }
  runner_i <- winner_i
  if (ncol(mat) > 1L && any(ok)) runner_i[ok] <- apply(mat[ok, , drop = FALSE], 1L,
    function(v) order(v, decreasing = TRUE, na.last = TRUE)[[2]])
  winner <- ifelse(is.na(winner_i), NA_character_, keys[winner_i])
  runner <- ifelse(is.na(runner_i), NA_character_, keys[runner_i])
  row <- seq_len(nrow(mat)); win_score <- run_score <- rep(NA_real_, nrow(mat))
  win_score[ok] <- mat[cbind(row[ok], winner_i[ok])]
  run_score[ok] <- mat[cbind(row[ok], runner_i[ok])]
  margin <- win_score - run_score
  tie_count <- vapply(row, function(i) if (is.finite(win_score[[i]]))
    sum(mat[i, ] == win_score[[i]], na.rm = TRUE) else 0L, integer(1))
  percentile <- rep(NA_real_, nrow(mat))
  for (key in keys) {
    pct <- rank(mat[, key], ties.method = "max", na.last = "keep") / nrow(mat)
    hit <- !is.na(winner) & winner == key; percentile[hit] <- pct[hit]
  }
  confidence <- margin * percentile
  min_winners <- as.integer(params$min_winners %||% 10L)
  counts <- table(winner[!is.na(winner) & nzchar(winner)])
  winner_count <- as.integer(counts[winner]); winner_count[is.na(winner_count)] <- 0L
  passes <- winner_count >= min_winners
  pruned <- ifelse(passes, winner, NA_character_)
  set_name <- ifelse(passes, paste0("set:winner_", .gp_winner_token(winner)), NA_character_)
  status <- ifelse(passes, "kept", "dropped_min_winners")
  status[is.na(winner) | !nzchar(winner)] <- "no_finite_scores"
  table <- data.frame(cell_id = cells, winner = ifelse(is.na(winner), "", winner),
    winner_score = win_score, runner_up = ifelse(is.na(runner), "", runner),
    runner_up_score = run_score, margin = margin, winner_percentile = percentile,
    confidence_score = confidence,
    confidence_rank_low = as.integer(rank(-confidence, ties.method = "max", na.last = "keep")),
    confidence_rank_high = as.integer(rank(-confidence, ties.method = "min", na.last = "keep")),
    tie_count = tie_count, winner_count = winner_count, passes_pruning = passes,
    pruning_status = status, pruned_winner = ifelse(is.na(pruned), "", pruned),
    pruned_cell_set = ifelse(is.na(set_name), "", set_name), stringsAsFactors = FALSE)
  meta <- list(source_score_artifact_id = .gp_scalar(params$source_score_artifact_id, NA_character_),
    source_run_id = .gp_scalar(params$source_run_id, NA_character_),
    source_snapshot_id = .gp_scalar(params$source_snapshot_id, NA_character_),
    vector_order = keys, winner_method = "argmax", tie_method = tie_method,
    min_winners = min_winners, confidence_method = "margin * winner_percentile",
    confidence_rank_method = "rank(-confidence_score), low=max, high=min",
    scworkbench_version = as.character(utils::packageVersion("scworkbench")))
  .gp_write(writers, "winners", scworkbench::wb_prim_obs_table(table,
    "cell_vector_winners", "cell", meta))
  invisible(list(status = "ok"))
}
.run_winner_cell_sets <- function(readers, writers, params, context) {
  winners <- readers$winners(); .gp_kind(winners, "obs_table")
  .gp_assert(identical(as.character(winners$table_type), "cell_vector_winners"),
    "winners must have table_type='cell_vector_winners'.")
  table <- winners$table; winner_col <- params$winner_col %||% "winner"
  passes_col <- params$passes_col %||% "passes_pruning"
  .gp_assert(all(c("cell_id", winner_col) %in% names(table)),
    "winner table is missing required columns.")
  keep <- if (passes_col %in% names(table)) as.logical(table[[passes_col]]) else rep(TRUE, nrow(table))
  keep[is.na(keep)] <- FALSE
  keys <- sort(unique(as.character(table[[winner_col]][keep])), method = "radix")
  keys <- keys[!is.na(keys) & nzchar(keys)]
  prefix <- as.character(params$cell_set_prefix %||% "winner_")[1]
  sets <- lapply(keys, function(key) {
    cells <- as.character(table$cell_id[keep & as.character(table[[winner_col]]) == key])
    list(name = paste0("set:", prefix, .gp_winner_token(key)), cells = cells,
      description = paste0("Winner cells for ", key, " from cell_vector_winners"),
      selection_rule = paste0(winner_col, " == '", key, "' && ", passes_col, " == TRUE"))
  })
  rows <- lapply(seq_along(keys), function(i) data.frame(
    cell_set = sets[[i]]$name, winner = keys[[i]], n_cells = length(sets[[i]]$cells),
    selection_rule = sets[[i]]$selection_rule,
    source_winner_artifact_id = .gp_scalar(params$source_winner_artifact_id),
    source_score_artifact_id = .gp_scalar(winners$meta$source_score_artifact_id),
    source_run_id = .gp_scalar(winners$meta$source_run_id),
    source_snapshot_id = .gp_scalar(winners$meta$source_snapshot_id),
    stringsAsFactors = FALSE))
  manifest <- if (length(rows)) do.call(rbind, rows) else data.frame(
    cell_set = character(), winner = character(), n_cells = integer(),
    selection_rule = character(), source_winner_artifact_id = character(),
    source_score_artifact_id = character(), source_run_id = character(),
    source_snapshot_id = character(), stringsAsFactors = FALSE)
  meta <- list(source_winner_artifact_id = .gp_scalar(params$source_winner_artifact_id),
    source_score_artifact_id = .gp_scalar(winners$meta$source_score_artifact_id),
    source_run_id = .gp_scalar(winners$meta$source_run_id),
    source_snapshot_id = .gp_scalar(winners$meta$source_snapshot_id),
    cell_set_prefix = prefix, winner_col = winner_col, passes_col = passes_col,
    replace = isTRUE(params$replace %||% TRUE))
  .gp_write(writers, "manifest", scworkbench::wb_prim_obs_table(manifest,
    "winner_cell_set_manifest", "cell_set", meta))
  writers$request_cell_set_family_admission("cell_set_family", list(
    kind = "cell_set_family_admission", sets = sets,
    mode = if (isTRUE(params$replace %||% TRUE)) "replace" else "add",
    source = list(label = as.character(params$source %||% "cell_vector_winners")[1],
      created_by = if (is.null(params$created_by)) NULL else as.character(params$created_by)[1],
      input_artifact_id = if (is.null(params$source_winner_artifact_id)) NULL else as.character(params$source_winner_artifact_id)[1],
      input_run_id = if (is.null(winners$meta$source_run_id) || is.na(winners$meta$source_run_id)) NULL else as.character(winners$meta$source_run_id)[1],
      input_snapshot_id = if (is.null(winners$meta$source_snapshot_id) || is.na(winners$meta$source_snapshot_id)) NULL else as.character(winners$meta$source_snapshot_id)[1])))
  invisible(list(status = "ok"))
}
