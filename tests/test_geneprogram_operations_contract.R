test_geneprogram_addmodulescore_contract <- function() {
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  result <- geneprogram_execute(f, "addmodulescore", list(matrix = f$matrix),
    list(features = list(A = f$features[1:2], B = f$features[3:4]), pool = f$features,
      nbin = 4L, ctrl = 2L, name = "Module", seed = 7L))
  scores <- geneprogram_output(result, "module_scores")
  stopifnot(identical(scores$key, rep(c("A", "B"), each = 8L)),
    identical(scores$cell_id, rep(f$cells, 2L)), all(is.finite(scores$value)))
  TRUE
}

test_geneprogram_ucell_typed_slice_contract <- function() {
  if (!requireNamespace("UCell", quietly = TRUE)) return(structure(TRUE,
    skip = "UCell is not installed"))
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  parameters <- list(features = list(A = f$features[1:2], B = f$features[3:4]),
    name = "Signature", max_rank = 100L, ncores = 1L,
    parallel_backend = "serial", missing_features_policy = "drop")
  one <- geneprogram_execute(f, "ucell", list(matrix = f$matrix),
    c(parameters, list(chunk_size = 8L)))
  backed_fixture <- geneprogram_fixture()
  on.exit(wb_close_bundle(backed_fixture$bundle), add = TRUE)
  backed <- structure(list(kind = "matrix", matrix_id = "RNA",
    dimensions = as.integer(dim(backed_fixture$values)),
    cell_id = backed_fixture$cells, feature_id = backed_fixture$features,
    storage_profile = "matrix_backed_chunked_v1"), class = "wb_primitive")
  several <- geneprogram_execute(backed_fixture, "ucell", list(matrix = backed),
    c(parameters, list(chunk_size = 3L)))
  a <- geneprogram_output(one, "signature_scores")
  b <- geneprogram_output(several, "signature_scores")
  stopifnot(identical(a$key, b$key), identical(a$cell_id, b$cell_id),
    isTRUE(all.equal(a$value, b$value, tolerance = 1e-12)))
  usage <- if (is.null(several$receipt$typed_slice_usage))
    several$receipt$typed_slice_reads else several$receipt$typed_slice_usage
  stopifnot(length(usage) == 1L)
  TRUE
}

test_geneprogram_reduce_contract <- function() {
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  result <- geneprogram_execute(f, "cell_vector_family_reduce", list(signal = f$signal),
    list(rule = "argmax", min_score = NULL, min_margin = .1,
      unknown_label = "unknown", ties = "unknown", tie_tol = 1e-8))
  stopifnot(identical(geneprogram_output(result, "label")$cell_id, f$cells),
    identical(geneprogram_output(result, "decision_table")$table$cell_id, f$cells))
  TRUE
}

test_geneprogram_calibrate_contract <- function() {
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  result <- geneprogram_execute(f, "cell_vector_family_calibrate", list(signal = f$signal),
    list(signal_keys = NULL, method = "percentile", center = TRUE, scale = TRUE,
      ties_method = "average", na_rm = TRUE))
  out <- geneprogram_output(result, "calibrated_signal")
  stopifnot(identical(out$key, f$signal$key), identical(out$cell_id, f$signal$cell_id),
    all(out$value >= 0 & out$value <= 1))
  TRUE
}

test_geneprogram_qc_contract <- function() {
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  result <- geneprogram_execute(f, "cell_vector_family_qc",
    list(signal = f$signal, grouping = f$grouping),
    list(signal_keys = NULL, correlation_methods = c("pearson", "spearman"), tie_tol = 1e-8))
  stopifnot(nrow(geneprogram_output(result, "key_summary")$table) == 2L,
    nrow(geneprogram_output(result, "pairwise_summary")$table) == 1L,
    nrow(geneprogram_output(result, "group_summary")$table) == 4L)
  TRUE
}

test_geneprogram_graph_smooth_contract <- function() {
  if (!requireNamespace("RANN", quietly = TRUE)) return(structure(TRUE,
    skip = "RANN is not installed"))
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  result <- geneprogram_execute(f, "cell_vector_family_graph_smooth",
    list(signal = f$signal, space = f$space), list(signal_keys = NULL, dims = NULL,
      k = 2L, alpha = .5, distance_metric = "euclidean", query_chunk_size = 3L,
      include_self = FALSE, na_policy = "ignore", source_signal_artifact_id = NULL,
      source_run_id = NULL, source_snapshot_id = NULL))
  out <- geneprogram_output(result, "smoothed_signal")
  stopifnot(identical(out$key, f$signal$key), identical(out$cell_id, f$signal$cell_id),
    all(is.finite(out$value)))
  TRUE
}

test_geneprogram_winners_contract <- function() {
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  result <- geneprogram_execute(f, "cell_vector_winners", list(signal = f$signal),
    list(signal_keys = NULL, min_winners = 1L, tie_method = "first",
      source_score_artifact_id = NULL, source_run_id = NULL, source_snapshot_id = NULL))
  winners <- geneprogram_output(result, "winners")
  stopifnot(nrow(winners$table) == 8L, all(winners$table$winner %in% c("A", "B")))
  TRUE
}

test_geneprogram_winner_cell_sets_contract <- function() {
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  winners_result <- geneprogram_execute(f, "cell_vector_winners", list(signal = f$signal),
    list(signal_keys = NULL, min_winners = 1L, tie_method = "first",
      source_score_artifact_id = NULL, source_run_id = NULL, source_snapshot_id = NULL))
  winners <- geneprogram_output(winners_result, "winners")
  result <- geneprogram_execute(f, "winner_cell_sets", list(winners = winners),
    list(source_winner_artifact_id = NULL, cell_set_prefix = "winner_",
      winner_col = "winner", passes_col = "passes_pruning",
      source = "cell_vector_winners", created_by = NULL, replace = TRUE))
  manifest <- geneprogram_output(result, "manifest")
  stopifnot(identical(manifest$table$cell_set, c("set:winner_A", "set:winner_B")),
    length(result$receipt$registry_requests) == 1L)
  TRUE
}

test_geneprogram_winner_cell_sets_empty_family_contract <- function() {
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  empty <- wb_prim_obs_table(data.frame(cell_id = f$cells,
    winner = rep(NA_character_, length(f$cells)),
    passes_pruning = rep(FALSE, length(f$cells)), stringsAsFactors = FALSE),
    "cell_vector_winners", "cell")
  run <- wb_session_run_op(wb_session_new(), f$bundle, "winner_cell_sets",
    inputs = list(winners = empty, cell_set_prefix = "winner_",
      winner_col = "winner", passes_col = "passes_pruning",
      replace = TRUE, min_cells = 1L))
  stopifnot(nrow(run$run$outputs$manifest$table) == 0L)
  promotion <- wb_session_promote_current(f$bundle, run$session,
    outputs = "manifest", created_by = "contract")
  stopifnot(all(promotion$promoted$status == "ok"),
    !any(startsWith(as.character(wb_cell_sets(f$bundle)), "set:winner_")))
  run$session <- wb_session_discard(run$session)
  TRUE
}

test_geneprogram_overlap_qc_contract <- function() {
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  result <- geneprogram_execute(f, "feature_set_overlap_qc", list(matrix = f$matrix),
    list(namespace = "programs", namespaces = NULL, feature_sets = NULL,
      min_set_size = 1L, max_set_size = NULL, min_detected_fraction = 0,
      deduplicate_sets = TRUE))
  qc <- geneprogram_output(result, "overlap_qc")
  stopifnot(nrow(qc$table) == 2L, all(qc$table$screening_status == "include"))
  TRUE
}

test_geneprogram_pathway_enrich_contract <- function() {
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle))
  result <- geneprogram_execute(f, "pathway_enrich",
    list(matrix = f$matrix, feature_table = f$feature_table),
    list(namespace = "programs", feature_col = "feature_id", p_adj_col = "p_val_adj",
      p_adj_max = .05, lfc_col = "avg_log1p_fc", direction = "any",
      lfc_min_abs = 0, lfc_min = .25, min_set_size = 1L,
      max_set_size = 100L, p_adjust_method = "BH"))
  out <- geneprogram_output(result, "pathways")
  stopifnot(nrow(out$table) == 2L, all(out$table$namespace == "programs"))
  TRUE
}

test_geneprogram_retained_parameter_declarations_contract <- function() {
  path <- system.file("scworkbench", "operation-pack.json",
    package = "scwbGenePrograms")
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  operations <- stats::setNames(manifest$operations,
    vapply(manifest$operations, `[[`, "", "name"))
  expected <- list(
    addmodulescore = c("module_definitions"),
    ucell = c("signature_definitions", "matrix_id", "n_cells"),
    cell_vector_family_reduce = c("effective_signal_keys", "n_cells"),
    cell_vector_family_calibrate = c("effective_signal_keys", "n_keys",
      "n_cells_per_key"),
    cell_vector_family_qc = c("effective_signal_keys", "n_keys", "n_cells",
      "has_grouping", "grouping_levels"),
    cell_vector_family_graph_smooth = c("effective_signal_keys", "n_keys",
      "n_cells", "source_space_id", "smoothing_method"),
    cell_vector_winners = c("effective_signal_keys", "n_cells",
      "winner_method", "confidence_method"),
    winner_cell_sets = c("n_cell_sets"),
    feature_set_overlap_qc = c("matrix_id", "effective_namespaces",
      "operation_marker", "effective_max_set_size"),
    pathway_enrich = c("matrix_id", "n_query", "n_universe", "n_pathways"))
  for (name in names(expected)) {
    retained <- operations[[name]]$retained_parameters
    if (is.null(retained)) retained <- list()
    observed <- vapply(retained, `[[`, "", "name")
    stopifnot(identical(observed, expected[[name]]))
  }
  TRUE
}

test_geneprogram_bounded_family_axis_fail_closed_contract <- function() {
  registry <- getFromNamespace(".wb_operation_pack_registry", "scworkbench")
  name <- "scwbgeneprograms::addmodulescore"
  registration <- registry$packs$scwbgeneprograms$registration
  original <- registration$runners[[name]]
  on.exit(registry$packs$scwbgeneprograms$registration$runners[[name]] <- original,
    add = TRUE)
  cases <- c("missing", "extra", "reordered", "duplicate")
  for (case in cases) {
    f <- geneprogram_fixture()
    faulty <- local({
      mode <- case
      function(readers, writers, params, context) {
        cells <- as.character(context$effective_cells)
        keys <- names(params$features)
        primitive <- scworkbench::wb_prim_cell_vector_family(
          rep(keys, each = length(cells)), rep(cells, length(keys)),
          seq_len(length(cells) * length(keys)), "continuous")
        if (mode == "missing") for (field in c("key", "cell_id", "value"))
          primitive[[field]] <- primitive[[field]][-1L]
        if (mode == "extra") {
          for (field in c("key", "cell_id", "value"))
            primitive[[field]] <- c(primitive[[field]], primitive[[field]][1L])
        }
        if (mode == "reordered") {
          order <- rev(seq_along(primitive$cell_id))
          for (field in c("key", "cell_id", "value"))
            primitive[[field]] <- primitive[[field]][order]
        }
        if (mode == "duplicate") primitive$cell_id[2L] <- primitive$cell_id[1L]
        writers$write_terminal_primitive("module_scores", primitive)
        invisible(list(status = "ok"))
      }
    })
    registry$packs$scwbgeneprograms$registration$runners[[name]] <- faulty
    before <- wb_art_list(f$bundle)
    error <- tryCatch(geneprogram_execute(f, "addmodulescore",
      list(matrix = f$matrix), list(features = list(A = f$features[1:2],
        B = f$features[3:4]), pool = f$features, nbin = 4L, ctrl = 2L,
        name = "Module", seed = 7L)), error = identity)
    stopifnot(inherits(error, "wb_operation_pack_error"),
      identical(before, wb_art_list(f$bundle)))
    wb_close_bundle(f$bundle)
  }
  TRUE
}

test_geneprogram_winner_family_atomicity_and_interruption_contract <- function() {
  run_winner <- function(fixture, value, replace) wb_session_run_op(
    wb_session_new(), fixture$bundle, "winner_cell_sets", inputs = list(
      winners = value, source_winner_artifact_id = NULL,
      cell_set_prefix = "winner_", winner_col = "winner",
      passes_col = "passes_pruning", source = "cell_vector_winners",
      created_by = NULL, replace = replace, min_cells = 1L))
  execute_winner <- function(fixture, value, replace) geneprogram_execute(
    fixture, "winner_cell_sets", list(winners = value), list(
      source_winner_artifact_id = NULL, cell_set_prefix = "winner_",
      winner_col = "winner", passes_col = "passes_pruning",
      source = "cell_vector_winners", created_by = NULL, replace = replace))
  f <- geneprogram_fixture(); on.exit(wb_close_bundle(f$bundle), add = TRUE)
  winners <- wb_prim_obs_table(data.frame(cell_id = f$cells,
    winner = rep(c("A", "B"), each = 4L), passes_pruning = TRUE,
    stringsAsFactors = FALSE), "cell_vector_winners", "cell")
  wb_cell_set_add(f$bundle, "set:winner_A", f$cells[[8L]])
  before_a <- wb_cells(f$bundle, "set:winner_A")
  old <- getOption("scworkbench.operation_pack_v2_test_hook")
  on.exit(options(scworkbench.operation_pack_v2_test_hook = old), add = TRUE)
  options(scworkbench.operation_pack_v2_test_hook = function(event, ...)
    if (identical(event, "before_publication")) stop("injected interruption"))
  interrupted <- tryCatch(execute_winner(f, winners, TRUE), error = identity)
  options(scworkbench.operation_pack_v2_test_hook = old)
  stopifnot(inherits(interrupted, "wb_operation_pack_error"),
    identical(before_a, wb_cells(f$bundle, "set:winner_A")),
    !"set:winner_B" %in% as.character(wb_cell_sets(f$bundle)))

  retry <- run_winner(f, winners, TRUE)
  retry_promotion <- wb_session_promote_current(f$bundle, retry$session,
    outputs = "manifest", created_by = "contract")
  stopifnot(identical(wb_cells(f$bundle, "set:winner_A"), f$cells[1:4]),
    identical(wb_cells(f$bundle, "set:winner_B"), f$cells[5:8]),
    all(retry_promotion$promoted$status == "ok"))
  retry$session <- wb_session_discard(retry$session)

  collision <- geneprogram_fixture()
  on.exit(wb_close_bundle(collision$bundle), add = TRUE)
  collision_winners <- wb_prim_obs_table(data.frame(cell_id = collision$cells,
    winner = rep(c("A", "B"), each = 4L), passes_pruning = TRUE,
    stringsAsFactors = FALSE), "cell_vector_winners", "cell")
  wb_cell_set_add(collision$bundle, "set:winner_A", collision$cells[[8L]])
  added <- run_winner(collision, collision_winners, FALSE)
  added_promotion <- wb_session_promote_current(collision$bundle, added$session,
    outputs = "manifest", created_by = "contract")
  stopifnot(identical(wb_cells(collision$bundle, "set:winner_A"),
      c(collision$cells[1:4], collision$cells[[8L]])),
    identical(wb_cells(collision$bundle, "set:winner_B"), collision$cells[5:8]),
    all(added_promotion$promoted$status == "ok"))
  added$session <- wb_session_discard(added$session)

  partial <- geneprogram_fixture()
  on.exit(wb_close_bundle(partial$bundle), add = TRUE)
  invalid_winners <- wb_prim_obs_table(data.frame(
    cell_id = c(partial$cells[1:4], ".not_in_bundle."),
    winner = c(rep("A", 4L), "B"), passes_pruning = TRUE,
    stringsAsFactors = FALSE), "cell_vector_winners", "cell")
  invalid <- tryCatch(execute_winner(partial, invalid_winners, TRUE), error = identity)
  stopifnot(inherits(invalid, "error"),
    !any(c("set:winner_A", "set:winner_B") %in% as.character(wb_cell_sets(partial$bundle))))
  TRUE
}
