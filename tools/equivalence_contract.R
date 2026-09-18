`%||%` <- function(x, y) if (is.null(x)) y else x

.gp_equivalence_abort <- function(code, method, detail) {
  condition <- structure(list(message = paste0(code, " [", method, "]: ",
    detail), call = NULL, code = code, method = method),
    class = c("geneprogram_equivalence_error", "error", "condition"))
  stop(condition)
}

.gp_json <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE,
  null = "null", digits = NA))
.gp_equal <- function(x, y) identical(.gp_json(x), .gp_json(y))
.gp_chars <- function(x) as.character(unlist(x %||% list(), use.names = FALSE))
.gp_scalar_chr <- function(x) as.character(unlist(x %||% NA_character_,
  use.names = FALSE))[[1L]]
.gp_scalar_num <- function(x) as.numeric(unlist(x %||% NA_real_,
  use.names = FALSE))[[1L]]
.gp_scalar_int <- function(x) as.integer(unlist(x %||% NA_integer_,
  use.names = FALSE))[[1L]]
.gp_scalar_lgl <- function(x) isTRUE(unlist(x %||% FALSE,
  use.names = FALSE)[[1L]])

.gp_parameter_fields <- list(
  addmodulescore = list(
    old = c("matrix_id", "n_modules", "module_names", "n_cells", "nbin",
      "ctrl", "seed"),
    new = c("module_definitions", "features", "pool", "nbin",
      "ctrl", "name", "seed")),
  ucell = list(
    old = c("matrix_id", "n_signatures", "signature_names", "n_cells",
      "max_rank", "ncores", "parallel_backend", "missing_features_policy",
      "chunking", "chunk_sizes"),
    new = c("signature_definitions", "matrix_id", "n_cells",
      "features", "name", "max_rank", "ncores", "parallel_backend",
      "chunk_size", "missing_features_policy")),
  cell_vector_family_reduce = list(
    old = c("rule", "signal_keys", "n_keys", "n_cells", "n_unknown",
      "min_score", "min_margin", "unknown_label", "ties", "tie_tol"),
    new = c("effective_signal_keys", "n_cells", "rule", "min_margin",
      "unknown_label", "ties", "tie_tol"),
    diagnostic_only = "n_unknown"),
  cell_vector_family_calibrate = list(
    old = c("method", "signal_keys", "n_keys", "n_rows", "center", "scale",
      "ties_method", "na_rm"),
    new = c("effective_signal_keys", "n_keys", "n_cells_per_key", "method",
      "center", "scale", "ties_method", "na_rm")),
  cell_vector_family_qc = list(
    old = c("signal_keys", "n_keys", "n_cells", "has_grouping", "n_groups",
      "correlation_methods", "tie_tol"),
    new = c("effective_signal_keys", "n_keys", "n_cells", "has_grouping",
      "grouping_levels", "correlation_methods", "tie_tol")),
  cell_vector_family_graph_smooth = list(
    old = c("signal_keys", "n_keys", "n_cells", "source_space_id", "dims",
      "k", "alpha", "smoothing_method", "distance_metric",
      "query_chunk_size", "include_self", "na_policy"),
    new = c("effective_signal_keys", "n_keys", "n_cells", "source_space_id",
      "smoothing_method", "k", "alpha", "distance_metric",
      "query_chunk_size", "include_self", "na_policy",
      "source_signal_artifact_id", "source_run_id", "source_snapshot_id")),
  cell_vector_winners = list(
    old = c("signal_keys", "n_keys", "n_cells", "min_winners", "tie_method",
      "winner_method", "confidence_method", "n_winners_kept",
      "n_winners_dropped"),
    new = c("effective_signal_keys", "n_cells", "winner_method",
      "confidence_method", "min_winners", "tie_method",
      "source_score_artifact_id", "source_run_id", "source_snapshot_id"),
    diagnostic_only = c("n_winners_kept", "n_winners_dropped")),
  winner_cell_sets = list(
    old = c("cell_set_prefix", "winner_col", "passes_col", "replace",
      "n_cell_sets", "source_winner_artifact_id", "source"),
    new = c("n_cell_sets", "source_winner_artifact_id", "cell_set_prefix",
      "winner_col", "passes_col", "source", "replace")),
  feature_set_overlap_qc = list(
    old = c("matrix_id", "namespaces", "op", "min_set_size", "max_set_size",
      "min_detected_fraction", "deduplicate_sets"),
    new = c("matrix_id", "effective_namespaces", "operation_marker",
      "effective_max_set_size", "namespace", "min_set_size",
      "min_detected_fraction", "deduplicate_sets")),
  pathway_enrich = list(
    old = c("matrix_id", "namespace", "n_query", "n_universe", "n_pathways",
      "direction", "p_adj_max", "min_set_size", "max_set_size",
      "p_adjust_method"),
    new = c("matrix_id", "n_query", "n_universe", "n_pathways", "namespace",
      "feature_col", "p_adj_col", "p_adj_max", "lfc_col", "direction",
      "lfc_min_abs", "lfc_min", "min_set_size", "max_set_size",
      "p_adjust_method")))

.gp_ucell_chunking <- function(record) {
  usage <- record$transition$receipt$typed_slice_usage[[1L]]
  list(enabled = as.integer(usage$call_count) > 1L, axis = "cells",
    chunk_size = as.integer(usage$requested_chunk_size),
    n_chunks = as.integer(usage$call_count))
}

.gp_effective_parameters <- function(method, record, side) {
  p <- record$run$params
  sym <- record$run$symbolic_inputs
  exec <- record$run$execution_inputs
  old <- identical(side, "old")
  keys <- function() if (old) .gp_chars(p$signal_keys) else
    .gp_chars(p$effective_signal_keys)
  switch(method,
    addmodulescore = list(
      matrix_id = if (old) .gp_scalar_chr(p$matrix_id) else .gp_scalar_chr(sym$matrix),
      n_modules = if (old) .gp_scalar_int(p$n_modules) else
        length(names(p$module_definitions)),
      module_names = if (old) .gp_chars(p$module_names) else names(p$module_definitions),
      n_cells = if (old) .gp_scalar_int(p$n_cells) else .gp_scalar_int(exec$effective_n),
      nbin = .gp_scalar_int(p$nbin), ctrl = .gp_scalar_int(p$ctrl),
      seed = .gp_scalar_int(p$seed)),
    ucell = list(
      matrix_id = .gp_scalar_chr(p$matrix_id),
      n_signatures = if (old) .gp_scalar_int(p$n_signatures) else
        length(names(p$signature_definitions)),
      signature_names = if (old) .gp_chars(p$signature_names) else
        names(p$signature_definitions),
      n_cells = .gp_scalar_int(p$n_cells), max_rank = .gp_scalar_int(p$max_rank),
      ncores = .gp_scalar_int(p$ncores),
      parallel_backend = .gp_scalar_chr(p$parallel_backend),
      missing_features_policy = .gp_scalar_chr(p$missing_features_policy),
      chunking = if (old) p$chunking else .gp_ucell_chunking(record),
      chunk_sizes = if (old) .gp_chars(p$chunk_sizes) else as.character(vapply(
        record$transition$receipt$typed_slice_usage[[1L]]$calls,
        function(x) as.integer(x$cell_count), integer(1)))),
    cell_vector_family_reduce = list(
      rule = .gp_scalar_chr(p$rule), signal_keys = keys(),
      n_keys = if (old) .gp_scalar_int(p$n_keys) else length(keys()),
      n_cells = .gp_scalar_int(p$n_cells), min_score = p$min_score %||% NULL,
      min_margin = .gp_scalar_num(p$min_margin),
      unknown_label = .gp_scalar_chr(p$unknown_label), ties = .gp_scalar_chr(p$ties),
      tie_tol = .gp_scalar_num(p$tie_tol)),
    cell_vector_family_calibrate = list(
      method = .gp_scalar_chr(p$method), signal_keys = keys(),
      n_keys = .gp_scalar_int(p$n_keys),
      n_rows = if (old) .gp_scalar_int(p$n_rows) else
        .gp_scalar_int(p$n_keys) * .gp_scalar_int(p$n_cells_per_key),
      center = .gp_scalar_lgl(p$center), scale = .gp_scalar_lgl(p$scale),
      ties_method = .gp_scalar_chr(p$ties_method), na_rm = .gp_scalar_lgl(p$na_rm)),
    cell_vector_family_qc = list(
      signal_keys = keys(), n_keys = .gp_scalar_int(p$n_keys),
      n_cells = .gp_scalar_int(p$n_cells),
      has_grouping = .gp_scalar_lgl(p$has_grouping),
      n_groups = if (old) .gp_scalar_int(p$n_groups) else
        length(.gp_chars(p$grouping_levels)),
      correlation_methods = .gp_chars(p$correlation_methods),
      tie_tol = .gp_scalar_num(p$tie_tol)),
    cell_vector_family_graph_smooth = list(
      signal_keys = keys(), n_keys = .gp_scalar_int(p$n_keys),
      n_cells = .gp_scalar_int(p$n_cells),
      source_space_id = .gp_scalar_chr(p$source_space_id),
      dims = p$dims %||% NULL, k = .gp_scalar_int(p$k),
      alpha = .gp_scalar_num(p$alpha),
      smoothing_method = .gp_scalar_chr(p$smoothing_method),
      distance_metric = .gp_scalar_chr(p$distance_metric),
      query_chunk_size = .gp_scalar_int(p$query_chunk_size),
      include_self = .gp_scalar_lgl(p$include_self),
      na_policy = .gp_scalar_chr(p$na_policy)),
    cell_vector_winners = list(
      signal_keys = keys(),
      n_keys = if (old) .gp_scalar_int(p$n_keys) else length(keys()),
      n_cells = .gp_scalar_int(p$n_cells),
      min_winners = .gp_scalar_int(p$min_winners),
      tie_method = .gp_scalar_chr(p$tie_method),
      winner_method = .gp_scalar_chr(p$winner_method),
      confidence_method = .gp_scalar_chr(p$confidence_method)),
    winner_cell_sets = list(
      cell_set_prefix = .gp_scalar_chr(p$cell_set_prefix),
      winner_col = .gp_scalar_chr(p$winner_col),
      passes_col = .gp_scalar_chr(p$passes_col), replace = .gp_scalar_lgl(p$replace),
      n_cell_sets = .gp_scalar_int(p$n_cell_sets),
      source_winner_artifact_id = .gp_scalar_chr(p$source_winner_artifact_id),
      source = .gp_scalar_chr(p$source)),
    feature_set_overlap_qc = list(
      matrix_id = .gp_scalar_chr(p$matrix_id),
      namespaces = if (old) .gp_chars(p$namespaces) else
        .gp_chars(p$effective_namespaces),
      op = if (old) .gp_scalar_chr(p$op) else .gp_scalar_chr(p$operation_marker),
      min_set_size = .gp_scalar_int(p$min_set_size),
      max_set_size = if (old) .gp_scalar_chr(p$max_set_size) else
        .gp_scalar_chr(p$effective_max_set_size),
      min_detected_fraction = .gp_scalar_num(p$min_detected_fraction),
      deduplicate_sets = .gp_scalar_lgl(p$deduplicate_sets)),
    pathway_enrich = list(
      matrix_id = .gp_scalar_chr(p$matrix_id), namespace = .gp_scalar_chr(p$namespace),
      n_query = .gp_scalar_int(p$n_query), n_universe = .gp_scalar_int(p$n_universe),
      n_pathways = .gp_scalar_int(p$n_pathways), direction = .gp_scalar_chr(p$direction),
      p_adj_max = .gp_scalar_num(p$p_adj_max),
      min_set_size = .gp_scalar_int(p$min_set_size),
      max_set_size = .gp_scalar_int(p$max_set_size),
      p_adjust_method = .gp_scalar_chr(p$p_adjust_method)))
}

.gp_input_rows <- function(x) {
  x <- x[, c("input_name", "primitive_kind", "content_identity_scheme",
    "content_sha256"), drop = FALSE]
  x[order(x$input_name, method = "radix"), , drop = FALSE]
}

.gp_binding_rows <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) return(data.frame())
  x <- x[as.character(x$direction) == "input", , drop = FALSE]
  rownames(x) <- NULL
  x
}

.gp_assert_self_identity <- function(method, value, identity, label) {
  if (!requireNamespace("scworkbench", quietly = TRUE)) .gp_equivalence_abort(
    "EQUIVALENCE_RECEIPT_MISSING", method,
    "scworkbench is required to authenticate transition evidence")
  candidate <- value
  candidate$.equivalence_identity <- identity
  assert <- getFromNamespace(".wb_operation_pack_v3_assert_self_identity", "scworkbench")
  valid <- tryCatch({
    assert(candidate, ".equivalence_identity", "EQUIVALENCE_RECEIPT_MISMATCH",
      label)
    TRUE
  }, error = function(error) FALSE)
  if (!valid) .gp_equivalence_abort("EQUIVALENCE_RECEIPT_MISMATCH", method,
    paste0(label, " self-identity is invalid"))
  TRUE
}

.gp_assert_receipt <- function(method, record) {
  transition <- record$transition
  if (!is.list(transition) || !is.list(transition$request) ||
      !is.list(transition$receipt)) .gp_equivalence_abort(
        "EQUIVALENCE_RECEIPT_MISSING", method,
        "authenticated published transition record is absent")
  request <- transition$request
  receipt <- transition$receipt
  if (!isTRUE(transition$authenticated_reopen) ||
      !isTRUE(transition$payloads_loaded)) .gp_equivalence_abort(
        "EQUIVALENCE_RECEIPT_MISSING", method,
        "transition was not reopened with payload validation")
  .gp_assert_self_identity(method, request, transition$request_identity,
    "published transition request")
  .gp_assert_self_identity(method, receipt, transition$receipt_identity,
    "published transition receipt")
  qualified <- paste0("scwbgeneprograms::", method)
  if (!identical(as.character(receipt$status), "ok") ||
      !identical(as.character(request$operation_authentication$name), qualified) ||
      !identical(as.character(receipt$transition_id),
        as.character(transition$transition_id))) .gp_equivalence_abort(
          "EQUIVALENCE_RECEIPT_MISMATCH", method,
          "transition status, operation authentication, or identity disagrees")
  observed_outputs <- sort(vapply(receipt$outputs %||% list(), `[[`, "", "name"),
    method = "radix")
  expected_outputs <- sort(as.character(record$outputs$output_name), method = "radix")
  if (!identical(observed_outputs, expected_outputs)) .gp_equivalence_abort(
    "EQUIVALENCE_RECEIPT_MISMATCH", method,
    "receipt output set differs from the published run output set")
  for (item in receipt$outputs %||% list()) {
    descriptor <- item$descriptor
    required <- c(descriptor$content_identity$scheme,
      descriptor$content_identity$sha256,
      descriptor$payload_manifest_identity$scheme,
      descriptor$payload_manifest_identity$sha256,
      descriptor$observed_descriptor_identity$sha256)
    if (length(required) != 5L || anyNA(required) || any(!nzchar(required)))
      .gp_equivalence_abort("EQUIVALENCE_RECEIPT_MISMATCH", method,
        paste0("output receipt identity is incomplete for ", item$name))
    descriptor_identity <- descriptor$observed_descriptor_identity
    descriptor$observed_descriptor_identity <- NULL
    .gp_assert_self_identity(method, descriptor, descriptor_identity,
      paste0("output descriptor ", item$name))
  }
  request_inputs <- sort(vapply(request$inputs %||% list(), `[[`, "", "name"),
    method = "radix")
  identity_inputs <- sort(as.character(record$input_identities$input_name),
    method = "radix")
  if (!identical(request_inputs, identity_inputs)) .gp_equivalence_abort(
    "EQUIVALENCE_INPUT_BINDING_MISMATCH", method,
    "request input descriptors do not cover the authoritative input ledger")
  if (identical(method, "ucell")) {
    usage <- receipt$typed_slice_usage %||% list()
    if (length(usage) != 1L || !identical(as.character(usage[[1L]]$coverage_status),
      "exact")) .gp_equivalence_abort("EQUIVALENCE_RECEIPT_MISMATCH", method,
        "typed-slice receipt is absent or does not prove exact coverage")
  }
  if (identical(method, "cell_vector_family_graph_smooth")) {
    requested <- request$provenance_reads %||% list()
    observed <- receipt$provenance_reads %||% list()
    if (length(requested) != 1L || length(observed) != 1L ||
        !identical(requested[[1L]]$descriptor_identity$sha256,
          observed[[1L]]$descriptor_identity$sha256))
      .gp_equivalence_abort("EQUIVALENCE_INPUT_BINDING_MISMATCH", method,
        "authenticated source-ancestry request and observation disagree")
  }
  if (method %in% c("feature_set_overlap_qc", "pathway_enrich")) {
    requested <- request$registry_reads %||% list()
    reads <- receipt$registry_reads %||% list()
    if (length(requested) != 1L || length(reads) != 1L ||
        !nzchar(as.character(reads[[1L]]$descriptor$content_identity$sha256 %||% "")))
      .gp_equivalence_abort("EQUIVALENCE_INPUT_BINDING_MISMATCH", method,
        "authenticated feature-set registry observation is absent")
    request_descriptor <- requested[[1L]]$descriptor
    receipt_descriptor <- reads[[1L]]$descriptor
    request_identity <- request_descriptor$descriptor_identity
    receipt_identity <- receipt_descriptor$observed_descriptor_identity
    request_descriptor$descriptor_identity <- NULL
    receipt_descriptor$observed_descriptor_identity <- NULL
    .gp_assert_self_identity(method, request_descriptor, request_identity,
      "registry request descriptor")
    .gp_assert_self_identity(method, receipt_descriptor, receipt_identity,
      "registry receipt descriptor")
    for (field in c("content_identity", "source_identity",
        "ancestry_identity", "payload_manifest_identity")) {
      if (!identical(request_descriptor[[field]], receipt_descriptor[[field]]))
        .gp_equivalence_abort("EQUIVALENCE_INPUT_BINDING_MISMATCH", method,
          paste0("registry request and receipt differ at ", field))
    }
  }
  if (identical(method, "winner_cell_sets") &&
      length(receipt$registry_requests %||% list()) != 1L)
    .gp_equivalence_abort("EQUIVALENCE_RECEIPT_MISMATCH", method,
      "cell-set-family registry receipt is absent")
  TRUE
}

compare_geneprogram_equivalence <- function(old, new, output = NULL) {
  expected_schema <- "phase3b_geneprograms_equivalence_snapshot.v2"
  if (!identical(old$schema_version, expected_schema) ||
      !identical(new$schema_version, expected_schema)) stop("snapshot schema mismatch")
  methods <- names(.gp_parameter_fields)
  if (!identical(names(old$methods), methods) ||
      !identical(names(new$methods), methods)) stop("method ledger mismatch")
  rows <- vector("list", length(methods))
  for (i in seq_along(methods)) {
    method <- methods[[i]]
    a <- old$methods[[method]]
    b <- new$methods[[method]]
    fields <- .gp_parameter_fields[[method]]
    missing_old <- setdiff(fields$old, names(a$run$params))
    missing_new <- setdiff(fields$new, names(b$run$params))
    extra_old <- setdiff(names(a$run$params), fields$old)
    extra_new <- setdiff(names(b$run$params), fields$new)
    if (length(missing_old) || length(missing_new) || length(extra_old) ||
        length(extra_new)) .gp_equivalence_abort(
      "EQUIVALENCE_PARAMETER_FIELD_UNMAPPED", method,
      paste(c(paste0("missing-old:", missing_old),
        paste0("missing-new:", missing_new), paste0("extra-old:", extra_old),
        paste0("extra-new:", extra_new)), collapse = ","))
    if (!.gp_equal(a$run$symbolic_inputs, b$run$symbolic_inputs))
      .gp_equivalence_abort("EQUIVALENCE_PARAMETER_MISMATCH", method,
        "persisted symbolic request inputs differ")
    if (!.gp_equal(.gp_effective_parameters(method, a, "old"),
        .gp_effective_parameters(method, b, "new")))
      .gp_equivalence_abort("EQUIVALENCE_PARAMETER_MISMATCH", method,
        "mapped completed-run effective parameters differ")
    if (!identical(.gp_input_rows(a$input_identities),
        .gp_input_rows(b$input_identities)) ||
        !identical(.gp_binding_rows(a$entity_bindings),
          .gp_binding_rows(b$entity_bindings)))
      .gp_equivalence_abort("EQUIVALENCE_INPUT_BINDING_MISMATCH", method,
        "scheme-qualified input identities or persisted entity bindings differ")
    if (!.gp_equal(a$registry_binding, b$registry_binding))
      .gp_equivalence_abort("EQUIVALENCE_INPUT_BINDING_MISMATCH", method,
        "authenticated registry selection, content, or source identity differs")
    if (!isTRUE(a$promotion$atomic) || !isTRUE(b$promotion$atomic))
      .gp_equivalence_abort("EQUIVALENCE_PROMOTION_ATOMICITY_MISMATCH", method,
        "promotion did not publish the complete declared terminal set")
    ao <- a$outputs[order(a$outputs$output_name, method = "radix"), , drop = FALSE]
    bo <- b$outputs[order(b$outputs$output_name, method = "radix"), , drop = FALSE]
    if (!identical(as.character(ao$output_name), as.character(bo$output_name)))
      .gp_equivalence_abort("EQUIVALENCE_OUTPUT_IDENTITY_MISMATCH", method,
        "terminal output sets differ")
    for (j in seq_len(nrow(ao))) {
      if (!identical(as.character(ao$replay_equivalence_policy[[j]]),
          "exact_canonical_value") ||
          !identical(as.character(bo$replay_equivalence_policy[[j]]),
            "exact_canonical_value"))
        .gp_equivalence_abort("EQUIVALENCE_TERMINAL_POLICY_MISMATCH", method,
          paste0("terminal policy differs for ", ao$output_name[[j]]))
      if (!identical(as.character(ao$primitive_kind[[j]]),
          as.character(bo$primitive_kind[[j]])) ||
          !identical(as.character(ao$content_identity_scheme[[j]]),
            as.character(bo$content_identity_scheme[[j]])) ||
          !identical(as.character(ao$content_sha256[[j]]),
            as.character(bo$content_sha256[[j]])) ||
          !identical(as.character(ao$replay_equivalence_reference[[j]]),
            as.character(ao$content_sha256[[j]])) ||
          !identical(as.character(bo$replay_equivalence_reference[[j]]),
            as.character(bo$content_sha256[[j]])) ||
          !identical(as.character(ao$verification_level[[j]]), "content") ||
          !identical(as.character(bo$verification_level[[j]]), "content"))
        .gp_equivalence_abort("EQUIVALENCE_OUTPUT_IDENTITY_MISMATCH", method,
          paste0("canonical kind/content/payload identity differs for ",
            ao$output_name[[j]]))
    }
    .gp_assert_receipt(method, b)
    rows[[i]] <- data.frame(method = method,
      output_count = nrow(ao), canonical_identity = TRUE,
      effective_parameters = TRUE, symbolic_inputs = TRUE,
      input_identities = TRUE, entity_bindings = TRUE,
      receipt_authenticated = TRUE, replay_policy = TRUE,
      promotion_atomicity = TRUE,
      diagnostic_only_fields = paste(fields$diagnostic_only %||% character(),
        collapse = ";"), stringsAsFactors = FALSE)
  }
  rows <- do.call(rbind, rows)
  if (!is.null(output)) utils::write.csv(rows, output, row.names = FALSE,
    quote = TRUE)
  rows
}
