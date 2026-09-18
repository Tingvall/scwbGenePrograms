#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
option <- function(name, required = TRUE, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) {
    if (required) stop("exactly one ", prefix, " argument is required")
    return(default)
  }
  if (length(hit) != 1L) stop("exactly one ", prefix, " argument is required")
  sub(prefix, "", hit, fixed = TRUE)
}

mode <- option("mode")
output <- normalizePath(option("output"), mustWork = FALSE)
fixture <- normalizePath(option("fixture"), mustWork = FALSE)
if (!mode %in% c("old-core", "v3-pack")) stop("unsupported mode")

if (mode == "old-core") {
  devtools::load_all(option("core-source"), quiet = TRUE)
} else {
  library(scworkbench)
  library(scwbGenePrograms)
  manifest_path <- system.file("scworkbench", "operation-pack.json",
    package = "scwbGenePrograms")
  wb_operation_pack_activate(manifest_path, lib.loc = .libPaths())
}

copy_directory <- function(from, to) {
  if (!dir.exists(from)) stop("fixture bundle does not exist: ", from)
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(from, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (length(entries) && !all(file.copy(entries, to, recursive = TRUE,
    copy.mode = TRUE, copy.date = TRUE))) stop("failed to copy fixture bundle")
  invisible(to)
}

make_fixture <- function(path) {
  b <- wb_init_bundle(path, overwrite = TRUE)
  on.exit(wb_close_bundle(b), add = TRUE)
  cells <- paste0("cell_", sprintf("%02d", seq_len(8)))
  features <- paste0("gene_", seq_len(8))
  values <- matrix(c(
    9,8,7,6,1,0,1,0, 8,7,6,5,0,1,0,1,
    1,0,1,0,9,8,7,6, 0,1,0,1,8,7,6,5,
    2,3,2,3,2,3,2,3, 3,2,3,2,3,2,3,2,
    4,4,3,3,2,2,1,1, 1,1,2,2,3,3,4,4),
    nrow = 8L, byrow = TRUE, dimnames = list(features, cells))
  wb_define_cells(b, cells)
  wb_attach_matrix(b, "RNA", Matrix::Matrix(values, sparse = TRUE),
    feature_type = "gene")
  wb_feature_set_replace(b, "program_a", "RNA", features[1:3],
    namespace = "programs", source_id = "fixture",
    recipe = list(n_genes_raw = 4L))
  wb_feature_set_replace(b, "program_b", "RNA", features[3:6],
    namespace = "programs", source_id = "fixture",
    recipe = list(n_genes_raw = 4L))
  coords <- cbind(x = seq_along(cells), y = c(1,2,1,2,1,2,1,2))
  rownames(coords) <- cells
  wb_add_space(b, "fixture_space", coords, overwrite = TRUE)

  signal <- wb_prim_cell_vector_family(rep(c("A", "B"), each = length(cells)),
    rep(cells, 2L), c(seq(.1, .8, .1), rev(seq(.1, .8, .1))),
    "continuous")
  signal_run <- wb_session_run_op(wb_session_new(), b,
    "cell_vector_family_calibrate", inputs = list(signal = signal,
      signal_keys = NULL, method = "percentile", center = TRUE, scale = TRUE,
      ties_method = "average", na_rm = TRUE, min_cells = 1L))
  signal_promotion <- wb_session_promote_current(b, signal_run$session,
    outputs = "calibrated_signal", name_prefix = "equivalence_source_signal",
    created_by = "phase3b-equivalence")
  signal_run$session <- wb_session_discard(signal_run$session)

  grouping_run <- wb_session_run_op(wb_session_new(), b, "cell_suffix_label",
    inputs = list(sep = "_", mapping = NULL, default_label = "UNKNOWN",
      min_cells = 1L))
  grouping_promotion <- wb_session_promote_current(b, grouping_run$session,
    outputs = "labels", name_prefix = "equivalence_source_grouping",
    created_by = "phase3b-equivalence")
  grouping_run$session <- wb_session_discard(grouping_run$session)

  feature_run <- wb_session_run_op(wb_session_new(), b, "hvg_select",
    inputs = list(matrix = "RNA", n_features = 8L, method = "variance",
      min_mean = NULL, max_mean = NULL, min_dispersion = NULL,
      batch_key = NULL, exclude_features = NULL, feature_type = NULL,
      seed = 1L, chunk_size = 4L, min_cells = 1L))
  feature_promotion <- wb_session_promote_current(b, feature_run$session,
    outputs = "selected_features", name_prefix = "equivalence_source_features",
    created_by = "phase3b-equivalence")
  feature_run$session <- wb_session_discard(feature_run$session)

  signal_artifact <- as.character(signal_promotion$promoted$artifact_id[1])
  winners_run <- wb_session_run_op(wb_session_new(), b,
    "cell_vector_winners", inputs = list(
      signal = wb_ref_artifact(signal_artifact, "calibrated_signal"),
      signal_keys = NULL,
      min_winners = 1L, tie_method = "first",
      source_score_artifact_id = signal_artifact,
      source_run_id = signal_run$run$run_id,
      source_snapshot_id = signal_promotion$provenance_snapshot_id,
      min_cells = 1L))
  winners_promotion <- wb_session_promote_current(b, winners_run$session,
    outputs = "winners", name_prefix = "equivalence_source_winners",
    created_by = "phase3b-equivalence")
  winners_run$session <- wb_session_discard(winners_run$session)

  manifest <- list(
    schema_version = "phase3b_geneprograms_equivalence_fixture.v2",
    cells = cells, features = features,
    signal = list(artifact_id = signal_artifact,
      run_id = signal_run$run$run_id,
      snapshot_id = signal_promotion$provenance_snapshot_id),
    grouping = list(
      artifact_id = as.character(grouping_promotion$promoted$artifact_id[1]),
      run_id = grouping_run$run$run_id,
      snapshot_id = grouping_promotion$provenance_snapshot_id),
    feature_table = list(
      artifact_id = as.character(feature_promotion$promoted$artifact_id[1]),
      run_id = feature_run$run$run_id,
      snapshot_id = feature_promotion$provenance_snapshot_id),
    winners = list(
      artifact_id = as.character(winners_promotion$promoted$artifact_id[1]),
      run_id = winners_run$run$run_id,
      snapshot_id = winners_promotion$provenance_snapshot_id))
  saveRDS(manifest, paste0(path, ".manifest.rds"), version = 3L)
  invisible(manifest)
}

if (!dir.exists(fixture) || !file.exists(paste0(fixture, ".manifest.rds"))) {
  if (mode != "old-core") stop("only old-core mode may create the fixture")
  dir.create(dirname(fixture), recursive = TRUE, showWarnings = FALSE)
  make_fixture(fixture)
}
fixture_manifest_path <- paste0(fixture, ".manifest.rds")
if (!file.exists(fixture_manifest_path)) stop("missing fixture manifest")
fixture_manifest <- readRDS(fixture_manifest_path)
if (!identical(fixture_manifest$schema_version,
  "phase3b_geneprograms_equivalence_fixture.v2")) stop("fixture schema mismatch")

artifact_outputs <- c(signal = "calibrated_signal", grouping = "labels",
  feature_table = "selected_features", winners = "winners")
artifact_ref <- function(name) wb_ref_artifact(
  fixture_manifest[[name]]$artifact_id, artifact_outputs[[name]])
space_primitive <- function() {
  coords <- cbind(x = seq_along(fixture_manifest$cells),
    y = c(1,2,1,2,1,2,1,2))
  rownames(coords) <- fixture_manifest$cells
  structure(list(kind = "space", space_id = "fixture_space", values = coords,
    dimensions = as.integer(dim(coords)), cell_id = fixture_manifest$cells,
    storage_profile = "space_dense_v1"), class = "wb_primitive")
}
specs <- list(
  addmodulescore = list(params = list(
    features = list(A = fixture_manifest$features[1:2],
      B = fixture_manifest$features[3:4]), pool = fixture_manifest$features,
    nbin = 4L, ctrl = 2L, name = "Module", seed = 7L)),
  ucell = list(params = list(
    features = list(A = fixture_manifest$features[1:2],
      B = fixture_manifest$features[3:4]), name = "Signature",
    max_rank = 100L, ncores = 1L, parallel_backend = "serial",
    chunk_size = 3L, missing_features_policy = "drop")),
  cell_vector_family_reduce = list(params = list(
    signal_keys = NULL, rule = "argmax", min_score = NULL, min_margin = .1,
    unknown_label = "unknown", ties = "unknown", tie_tol = 1e-8)),
  cell_vector_family_calibrate = list(params = list(
    signal_keys = NULL, method = "percentile", center = TRUE, scale = TRUE,
    ties_method = "average", na_rm = TRUE)),
  cell_vector_family_qc = list(params = list(
    signal_keys = NULL, correlation_methods = c("pearson", "spearman"),
    tie_tol = 1e-8)),
  cell_vector_family_graph_smooth = list(params = list(
    signal_keys = NULL, dims = NULL, k = 2L, alpha = .5,
    distance_metric = "euclidean", query_chunk_size = 3L,
    include_self = FALSE, na_policy = "ignore",
    source_signal_artifact_id = fixture_manifest$signal$artifact_id,
    source_run_id = fixture_manifest$signal$run_id,
    source_snapshot_id = fixture_manifest$signal$snapshot_id)),
  cell_vector_winners = list(params = list(
    signal_keys = NULL, min_winners = 1L, tie_method = "first",
    source_score_artifact_id = fixture_manifest$signal$artifact_id,
    source_run_id = fixture_manifest$signal$run_id,
    source_snapshot_id = fixture_manifest$signal$snapshot_id)),
  winner_cell_sets = list(params = list(
    source_winner_artifact_id = fixture_manifest$winners$artifact_id,
    cell_set_prefix = "winner_", winner_col = "winner",
    passes_col = "passes_pruning", source = "cell_vector_winners",
    created_by = NULL, replace = TRUE)),
  feature_set_overlap_qc = list(params = list(
    namespace = "programs", namespaces = NULL, feature_sets = NULL,
    min_set_size = 1L, max_set_size = NULL, min_detected_fraction = 0,
    deduplicate_sets = TRUE)),
  pathway_enrich = list(params = list(namespace = "programs",
    feature_col = "feature_id", p_adj_col = "not_present", p_adj_max = .05,
    lfc_col = "not_present", direction = "any", lfc_min_abs = 0,
    lfc_min = .25, min_set_size = 1L, max_set_size = 100L,
    p_adjust_method = "BH")))

method_inputs <- function(name) switch(name,
  addmodulescore = list(matrix = "RNA"),
  ucell = list(matrix = "RNA"),
  cell_vector_family_reduce = list(signal = artifact_ref("signal")),
  cell_vector_family_calibrate = list(signal = artifact_ref("signal")),
  cell_vector_family_qc = list(signal = artifact_ref("signal"),
    grouping = artifact_ref("grouping")),
  cell_vector_family_graph_smooth = list(signal = artifact_ref("signal"),
    space = space_primitive()),
  cell_vector_winners = list(signal = artifact_ref("signal")),
  winner_cell_sets = list(winners = artifact_ref("winners")),
  feature_set_overlap_qc = list(matrix = "RNA"),
  pathway_enrich = list(matrix = "RNA",
    feature_table = artifact_ref("feature_table")))

resolved_inputs <- function(bundle, name) {
  load <- function(key) wb_art_load(bundle,
    fixture_manifest[[key]]$artifact_id)
  switch(name,
    addmodulescore = list(matrix = wb_prim_matrix("RNA")),
    ucell = list(matrix = wb_prim_matrix("RNA")),
    cell_vector_family_reduce = list(signal = load("signal")),
    cell_vector_family_calibrate = list(signal = load("signal")),
    cell_vector_family_qc = list(signal = load("signal"),
      grouping = load("grouping")),
    cell_vector_family_graph_smooth = list(signal = load("signal"),
      space = space_primitive()),
    cell_vector_winners = list(signal = load("signal")),
    winner_cell_sets = list(winners = load("winners")),
    feature_set_overlap_qc = list(matrix = wb_prim_matrix("RNA")),
    pathway_enrich = list(matrix = wb_prim_matrix("RNA"),
      feature_table = load("feature_table")))
}

execute_old <- function(bundle, name, resolved, params) {
  if (identical(name, "winner_cell_sets")) {
    value <- wb_session_run_op(wb_session_new(), bundle, name,
      inputs = c(method_inputs(name), params, list(min_cells = 1L)),
      notes = "phase3b old-core equivalence")
    promotion <- wb_session_promote_current(bundle, value$session,
      name_prefix = paste0("old_", name),
      created_by = "phase3b-equivalence", outputs = names(value$run$outputs))
    value$session <- wb_session_discard(value$session)
    return(list(run = value$run, result = value$result,
      promotion = promotion))
  }
  value <- switch(name,
    addmodulescore = wb_op_addmodulescore(bundle, "RNA", params$features,
      params$pool, fixture_manifest$cells, params$nbin, params$ctrl,
      params$name, params$seed),
    ucell = wb_op_ucell(bundle, "RNA", params$features,
      fixture_manifest$cells, params$name, params$max_rank, params$ncores,
      params$parallel_backend, params$chunk_size,
      params$missing_features_policy),
    cell_vector_family_reduce = wb_op_cell_vector_family_reduce(resolved$signal,
      params$signal_keys, params$rule, params$min_score, params$min_margin,
      params$unknown_label, params$ties, params$tie_tol),
    cell_vector_family_calibrate = wb_op_cell_vector_family_calibrate(
      resolved$signal, params$signal_keys, params$method, params$center,
      params$scale, params$ties_method, params$na_rm),
    cell_vector_family_qc = wb_op_cell_vector_family_qc(resolved$signal,
      params$signal_keys, resolved$grouping, params$correlation_methods,
      params$tie_tol),
    cell_vector_family_graph_smooth = wb_op_cell_vector_family_graph_smooth(
      bundle, resolved$signal, "fixture_space", fixture_manifest$cells,
      params$signal_keys, params$dims, params$k, params$alpha,
      params$distance_metric, params$query_chunk_size, params$include_self,
      params$na_policy,
      source_signal_artifact_id = params$source_signal_artifact_id,
      source_run_id = params$source_run_id,
      source_snapshot_id = params$source_snapshot_id),
    cell_vector_winners = wb_op_cell_vector_winners(resolved$signal,
      params$signal_keys, params$min_winners, params$tie_method,
      params$source_score_artifact_id, params$source_run_id,
      params$source_snapshot_id),
    feature_set_overlap_qc = wb_op_feature_set_overlap_qc(bundle, "RNA",
      namespace = params$namespace, namespaces = params$namespaces,
      feature_sets = params$feature_sets, min_set_size = params$min_set_size,
      max_set_size = params$max_set_size %||% Inf,
      min_detected_fraction = params$min_detected_fraction,
      deduplicate_sets = params$deduplicate_sets),
    pathway_enrich = wb_op_pathway_enrich(bundle, "RNA",
      resolved$feature_table$table, params$namespace, params$feature_col,
      params$p_adj_col, params$p_adj_max, params$lfc_col, params$direction,
      params$lfc_min_abs, params$lfc_min, params$min_set_size,
      params$max_set_size, params$p_adjust_method))
  user_inputs <- c(method_inputs(name), params, list(min_cells = 1L))
  run <- wb_new_run(name,
    inputs = list(cell_set = "all", effective_n = length(fixture_manifest$cells),
      parents = character(), provenance_refs = method_inputs(name),
      user_inputs = user_inputs,
      symbolic_inputs = getFromNamespace(".wb_symbolic_inputs_snapshot",
        "scworkbench")(name, user_inputs),
      execution_inputs = getFromNamespace(".wb_execution_inputs_snapshot",
        "scworkbench")(resolved, fixture_manifest$cells)),
    params = value$params, outputs = value$outputs,
    notes = "phase3b old-core equivalence")
  runs <- stats::setNames(list(run), run$run_id)
  promotion <- wb_promote_analysis(bundle, runs, run$run_id,
    name_prefix = paste0("old_", name), created_by = "phase3b-equivalence",
    outputs = names(run$outputs))
  list(run = run, result = value, promotion = promotion)
}

execute_new <- function(bundle, name, params) {
  inputs <- c(method_inputs(name), params, list(min_cells = 1L))
  value <- wb_session_run_op(wb_session_new(), bundle,
    paste0("scwbgeneprograms::", name), inputs = inputs,
    notes = "phase3b v3-pack equivalence")
  promotion <- wb_session_promote_current(bundle, value$session,
    name_prefix = paste0("v3_", name), created_by = "phase3b-equivalence",
    outputs = names(value$run$outputs))
  value$session <- wb_session_discard(value$session)
  list(run = value$run, result = value$result, promotion = promotion)
}

canonical_binding_rows <- function(details) {
  rows <- details$entity_bindings
  if (!is.data.frame(rows) || !nrow(rows)) return(data.frame())
  keep <- intersect(c("direction", "binding_name", "entity_kind", "entity_id",
    "relationship_type", "source", "content_hash", "content_identity_scheme"),
    names(rows))
  rows <- unique(rows[, keep, drop = FALSE])
  rows[do.call(order, c(unname(rows), list(na.last = TRUE))), , drop = FALSE]
}

input_identity_rows <- function(bundle, resolved) {
  rows <- lapply(names(resolved), function(name) {
    primitive <- resolved[[name]]
    identity <- getFromNamespace(".wb_prov_payload_identity",
      "scworkbench")(bundle, primitive)
    data.frame(input_name = name, primitive_kind = as.character(primitive$kind),
      content_identity_scheme = as.character(identity$content_identity_scheme),
      content_sha256 = as.character(identity$content_sha256),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

registry_binding <- function(bundle, name, resolved, params, transition) {
  if (!name %in% c("feature_set_overlap_qc", "pathway_enrich")) return(NULL)
  declaration <- list(name = "feature_sets", matrix_input = "matrix",
    namespace_parameter = "namespace", no_selector_policy = "all_namespaces",
    empty_match_policy = "error")
  if (identical(name, "feature_set_overlap_qc")) {
    declaration$namespaces_parameter <- "namespaces"
    declaration$feature_sets_parameter <- "feature_sets"
  }
  snapshot <- getFromNamespace(
    ".wb_operation_pack_v2_feature_set_registry_snapshot", "scworkbench")(
      bundle, declaration, resolved, params)
  scoped <- getFromNamespace(".wb_operation_pack_v2_identity_scoped",
    "scworkbench")
  identities <- list(
    selection_identity = scoped(
      "scworkbench.feature_set_registry.selection.v1",
      "normalized selector binding", snapshot$selection),
    content_identity = scoped(
      "scworkbench.feature_set_registry.snapshot.v1",
      "canonical matrix-bound registry snapshot", snapshot),
    source_identity = scoped(
      "scworkbench.feature_set_registry.source.v1",
      "ordered registry source state", list(matrix_id = snapshot$matrix_id,
        feature_universe = snapshot$feature_universe, sets = snapshot$sets)))
  if (!is.null(transition)) {
    request <- transition$request$registry_reads[[1L]]$descriptor
    receipt <- transition$receipt$registry_reads[[1L]]$descriptor
    if (!identical(request$selection_identity, identities$selection_identity) ||
        !identical(request$content_identity, identities$content_identity) ||
        !identical(request$source_identity, identities$source_identity) ||
        !identical(receipt$content_identity, identities$content_identity) ||
        !identical(receipt$source_identity, identities$source_identity)) {
      stop("published registry identity differs from the authoritative snapshot")
    }
  }
  identities
}

snapshot_one <- function(name) {
  cat("SNAPSHOT METHOD mode=", mode, " method=", name, "\n", sep = "")
  target <- tempfile(paste0("geneprograms-equivalence-", name, "-"))
  copy_directory(fixture, target)
  b <- wb_open_bundle(target)
  on.exit({
    try(wb_close_bundle(b), silent = TRUE)
    unlink(target, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  resolved <- resolved_inputs(b, name)
  executed <- if (mode == "old-core") execute_old(b, name, resolved,
    specs[[name]]$params) else execute_new(b, name, specs[[name]]$params)
  sid <- as.character(executed$promotion$provenance_snapshot_id)
  details <- wb_prov_run_details(b, executed$run$run_id, sid)
  transition_record <- NULL
  if (mode == "v3-pack") {
    transition <- executed$run$operation_pack_transition
    transition_path <- file.path(getFromNamespace(
      ".wb_operation_pack_v2_store", "scworkbench")(b),
      as.character(transition$transition_id))
    request_with_identity <- jsonlite::read_json(
      file.path(transition_path, "request.json"), simplifyVector = FALSE)
    receipt_with_identity <- jsonlite::read_json(
      file.path(transition_path, "receipt.json"), simplifyVector = FALSE)
    validated <- getFromNamespace(".wb_operation_pack_v2_record_validate",
      "scworkbench")(transition_path, load_payloads = TRUE)
    if (!isTRUE(validated$ok) ||
        !identical(as.character(validated$state$state), "reopened")) {
      stop("published transition receipt failed authenticated reopen")
    }
    transition_record <- list(transition_id = transition$transition_id,
      request = validated$request, receipt = validated$receipt,
      request_identity = request_with_identity$request_identity,
      receipt_identity = receipt_with_identity$receipt_identity,
      authenticated_reopen = TRUE, payloads_loaded = TRUE)
  }
  output_rows <- DBI::dbGetQuery(b$con,
    paste(
      "SELECT output_name, primitive_kind, content_sha256,",
      "content_identity_scheme, verification_level, verification_status,",
      "replay_equivalence_policy, replay_equivalence_reference,",
      "replay_equivalence_status FROM prov_run_outputs",
      "WHERE snapshot_id=? AND run_id=? ORDER BY output_name;"),
    params = list(sid, executed$run$run_id))
  promoted <- executed$promotion$promoted
  expected_outputs <- sort(names(executed$run$outputs), method = "radix")
  promotion_atomic <- is.data.frame(promoted) &&
    all(as.character(promoted$status) == "ok") &&
    identical(sort(as.character(promoted$output), method = "radix"),
      expected_outputs) &&
    identical(sort(as.character(output_rows$output_name), method = "radix"),
      expected_outputs)
  list(
    operation_name = as.character(details$run$op_name),
    run = list(params = details$run$params,
      symbolic_inputs = details$run$symbolic_inputs,
      execution_inputs = details$run$execution_inputs),
    input_identities = input_identity_rows(b, resolved),
    registry_binding = registry_binding(b, name, resolved,
      specs[[name]]$params, transition_record),
    entity_bindings = canonical_binding_rows(details),
    outputs = output_rows,
    promotion = list(atomic = promotion_atomic,
      output_names = expected_outputs, snapshot_id = sid),
    transition = transition_record)
}

result <- lapply(names(specs), snapshot_one)
names(result) <- names(specs)
snapshot <- list(schema_version = "phase3b_geneprograms_equivalence_snapshot.v2",
  mode = mode, fixture = fixture_manifest, methods = result)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(snapshot, output, version = 3L)
cat("SNAPSHOT OK mode=", mode, " methods=", length(result),
  " outputs=", sum(vapply(result, function(x) nrow(x$outputs), integer(1))),
  " output=", output, "\n", sep = "")
