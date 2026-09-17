args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE)

string <- function(default = NULL, choices = NULL, nullable = is.null(default)) {
  out <- list(type = if (nullable) c("string", "null") else "string")
  if (!is.null(default)) out$default <- default
  if (!is.null(choices)) out$enum <- choices
  if (!nullable) out$minLength <- 1L
  out
}
number <- function(default = NULL, minimum = NULL, maximum = NULL, nullable = is.null(default)) {
  out <- list(type = if (nullable) c("number", "null") else "number")
  if (!is.null(default)) out$default <- default
  if (!is.null(minimum)) out$minimum <- minimum
  if (!is.null(maximum)) out$maximum <- maximum
  out
}
integer <- function(default = NULL, minimum = NULL, maximum = NULL, nullable = is.null(default)) {
  out <- list(type = if (nullable) c("integer", "null") else "integer")
  if (!is.null(default)) out$default <- default
  if (!is.null(minimum)) out$minimum <- minimum
  if (!is.null(maximum)) out$maximum <- maximum
  out
}
boolean <- function(default) list(type = "boolean", default = default)
strings <- function(default = NULL, required = FALSE) {
  out <- list(type = if (required) "array" else c("array", "null"),
    uniqueItems = TRUE, items = list(type = "string", minLength = 1L))
  if (!required) out$default <- default
  if (required) out$minItems <- 1L
  out
}
integers <- function(default = NULL, required = FALSE) {
  out <- list(type = if (required) "array" else c("array", "null"),
    uniqueItems = TRUE, items = list(type = "integer", minimum = 1L))
  if (!required) out$default <- default
  if (required) out$minItems <- 1L
  out
}
groups <- function() list(type = "object", minProperties = 1L,
  additionalProperties = list(type = "array", minItems = 1L,
    uniqueItems = TRUE, items = list(type = "string", minLength = 1L)))
schema <- function(properties, required = character()) {
  out <- list(type = "object", properties = properties,
    additionalProperties = FALSE)
  if (length(required)) out$required <- unname(as.list(required))
  out
}
input <- function(name, kind, required = TRUE, alignment = if (kind == "matrix")
  "cell_and_feature_ordered" else if (kind %in% c("cell_vector", "cell_vector_family", "space"))
  "cell_ordered" else if (kind == "feature_table") "feature_ordered" else "none")
  list(name = name, kind = kind, required = required, alignment = alignment)
output <- function(name, kind, group = name, table_type = NULL, value_kind = NULL,
                   category = "terminal_artifact", profile = "terminal_primitive_v1") {
  out <- list(name = name, category = category, kind = kind, storage_profile = profile,
    allocated_id_key = paste0(name, "_id"), required = TRUE, promotion_group = group)
  if (!is.null(table_type)) out$table_type <- table_type
  if (!is.null(value_kind)) out$value_kind <- value_kind
  out
}
table_policy <- function(name) list(output = name, policy = "data_dependent_table_v1")
cell_policy <- function(name, input_name, missing = "canonical_null_allowed")
  list(output = name, policy = "authenticated_output_axes_v1", shape = "one_per_cell",
    cell_axis = list(source = "distinct_input_cells_first_occurrence_v1", input = input_name),
    row_order = "input_cell_order_v1", missing_values = missing)
family_policy <- function(name, input_name = NULL, selector = "signal_keys",
                          groups_parameter = NULL, name_parameter = NULL,
                          missing_parameter = NULL, missing = "canonical_null_allowed") {
  if (!is.null(groups_parameter)) {
    key_axis <- list(source = "resolved_parameter_feature_groups_v1",
      groups_parameter = groups_parameter, matrix_input = input_name,
      name_prefix_parameter = name_parameter, matching = "exact_feature_id_v1",
      normalization = "ordered_unique_features_v1",
      generated_names = "legacy_module_group_names_v1",
      missing_feature_policy = if (is.null(missing_parameter)) list(source = "fixed_drop_v1") else
        list(source = "parameter_v1", parameter = missing_parameter),
      empty_group_policy = "drop", zero_key_policy = "error")
    cell_source <- "effective_cells_v1"
  } else {
    key_axis <- list(source = "input_keys_with_optional_selector_v1", input = input_name,
      selector_parameter = selector, no_selector_order = "first_occurrence_v1",
      selector_order = "normalized_parameter_order_v1", unknown_key_policy = "error",
      zero_key_policy = "error")
    cell_source <- "effective_cells_v1"
  }
  list(output = name, policy = "authenticated_output_axes_v1", shape = "cell_key_cartesian",
    cell_axis = list(source = cell_source, input = input_name), key_axis = key_axis,
    row_order = "key_major_then_cell_v1", missing_values = missing)
}
dependency <- function(name, probe, requirement_id = paste0(tolower(name), "_runtime"),
                       spec = ">=0.1.0") list(requirement_id = requirement_id,
  kind = "r_package", name = name,
  version = list(scheme = "r_package_version_constraint_v1", specifier = spec),
  disposition = "required", import_probe = probe)
core <- dependency("scworkbench", "wb_prim_obs_table", spec = ">=0.1.0,<0.2.0")
op <- function(name, title, description, runner, inputs, params_schema, policies,
               outputs, dependencies = list(core), closure = NULL, replay = NULL,
               registry_reads = NULL, typed_slice_reads = NULL, provenance_reads = NULL) {
  out <- list(name = name, qualified_name = paste0("scwbgeneprograms::", name),
    title = title, description = description, runner = runner, inputs = inputs,
    params_schema = params_schema, dimension_policies = policies, outputs = outputs,
    dependencies = dependencies,
    promotion_closure = lapply(closure %||% stats::setNames(lapply(
      vapply(outputs, `[[`, "", "name"), function(x) x),
      vapply(outputs, `[[`, "", "name")), function(x) unname(as.list(sort(x, method = "radix")))),
    replay_equivalence = as.list(replay %||% stats::setNames(
      rep("exact_canonical_value", length(outputs)),
      vapply(outputs, `[[`, "", "name"))))
  if (!is.null(registry_reads)) out$registry_reads <- registry_reads
  if (!is.null(typed_slice_reads)) out$typed_slice_reads <- typed_slice_reads
  if (!is.null(provenance_reads)) out$provenance_reads <- provenance_reads
  out
}
`%||%` <- function(x, y) if (is.null(x)) y else x

signal_params <- list(signal_keys = strings(), rule = string("argmax", "argmax", FALSE),
  min_score = number(), min_margin = number(), unknown_label = string("unknown", nullable = FALSE),
  ties = string("unknown", c("unknown", "first"), FALSE), tie_tol = number(1e-8, 0, 1, FALSE))

operations <- list(
  op("addmodulescore", "Add module score",
    "Seurat-style module scoring using matched control features by expression bin.",
    ".run_addmodulescore", list(input("matrix", "matrix")),
    schema(list(features = groups(), pool = strings(), nbin = integer(24L, 2L, 200L, FALSE),
      ctrl = integer(100L, 1L, 1000L, FALSE), name = string("Module", nullable = FALSE),
      seed = integer(1L, nullable = FALSE)), "features"),
    list(family_policy("module_scores", "matrix", groups_parameter = "features",
      name_parameter = "name", missing = "forbidden")),
    list(output("module_scores", "cell_vector_family", value_kind = "continuous")),
    list(core, dependency("Matrix", "colMeans", spec = ">=1.3.0"))),

  op("ucell", "UCell Signature Scoring",
    "Rank-based per-cell signature scoring using pack-owned chunking over canonical typed slices.",
    ".run_ucell", list(input("matrix", "matrix")),
    schema(list(features = groups(), name = string("Signature", nullable = FALSE),
      max_rank = integer(1500L, 1L, 1000000L, FALSE), ncores = integer(1L, 1L, 1024L, FALSE),
      parallel_backend = string("auto", c("auto", "serial"), FALSE),
      chunk_size = integer(2000L, 1L, 1000000L, FALSE),
      missing_features_policy = string("drop", c("drop", "error"), FALSE)), "features"),
    list(family_policy("signature_scores", "matrix", groups_parameter = "features",
      name_parameter = "name", missing_parameter = "missing_features_policy", missing = "forbidden"),
      table_policy("signature_summary")),
    list(output("signature_scores", "cell_vector_family", "ucell", value_kind = "continuous"),
      output("signature_summary", "obs_table", "ucell", "ucell_signature_summary")),
    list(core, dependency("Matrix", "colMeans", spec = ">=1.3.0"),
      dependency("UCell", "ScoreSignatures_UCell", spec = ">=2.0.0"),
      dependency("BiocParallel", "SerialParam", spec = ">=1.28.0")),
    closure = list(signature_scores = c("signature_scores", "signature_summary"),
      signature_summary = c("signature_scores", "signature_summary")),
    typed_slice_reads = list(list(name = "matrix_chunks", input = "matrix",
      reader = "read_typed_slices", kind = "matrix", supported_axes = list("cells"),
      chunk_axis = "cells", chunk_size_parameter = "chunk_size",
      call_discipline = "ordered_partition_exactly_once", empty_selection_policy = "reject",
      duplicate_policy = "reject", selection = list(cells = list(policy = "effective_input_axis_v1"),
        features = list(policy = "effective_input_axis_v1"))))),

  op("cell_vector_family_reduce", "Reduce Cell Vector Family",
    "Reduce a keyed signal family to labels and deterministic diagnostics.", ".run_family_reduce",
    list(input("signal", "cell_vector_family")), schema(signal_params),
    c(lapply(c("label", "top_score", "runner_up_score", "score_margin", "abstain_reason"),
      cell_policy, input_name = "signal"), list(table_policy("decision_table"), table_policy("label_summary"))),
    list(output("label", "cell_vector", value_kind = "categorical"),
      output("top_score", "cell_vector", value_kind = "continuous"),
      output("runner_up_score", "cell_vector", value_kind = "continuous"),
      output("score_margin", "cell_vector", value_kind = "continuous"),
      output("abstain_reason", "cell_vector", value_kind = "categorical"),
      output("decision_table", "obs_table", table_type = "cell_vector_family_reduce_decision_table"),
      output("label_summary", "obs_table", table_type = "cell_vector_family_reduce_summary"))),

  op("cell_vector_family_calibrate", "Calibrate Cell Vector Family",
    "Calibrate each selected signal key independently.", ".run_family_calibrate",
    list(input("signal", "cell_vector_family")), schema(list(signal_keys = strings(),
      method = string("percentile", c("percentile", "zscore", "robust_zscore"), FALSE),
      center = boolean(TRUE), scale = boolean(TRUE),
      ties_method = string("average", c("average", "first", "random", "max", "min"), FALSE),
      na_rm = boolean(TRUE))),
    list(family_policy("calibrated_signal", "signal"), table_policy("calibration_summary")),
    list(output("calibrated_signal", "cell_vector_family", "calibration", value_kind = "continuous"),
      output("calibration_summary", "obs_table", "calibration", "cell_vector_family_calibration_summary")),
    closure = list(calibrated_signal = c("calibrated_signal", "calibration_summary"),
      calibration_summary = c("calibrated_signal", "calibration_summary"))),

  op("cell_vector_family_qc", "QC Cell Vector Family",
    "Compute per-key, pairwise, and optional grouped diagnostics.", ".run_family_qc",
    list(input("signal", "cell_vector_family"), input("grouping", "cell_vector", FALSE)),
    schema(list(signal_keys = strings(), correlation_methods = strings(c("pearson", "spearman")),
      tie_tol = number(1e-8, 0, 1, FALSE))),
    lapply(c("key_summary", "pairwise_summary", "group_summary"), table_policy),
    list(output("key_summary", "obs_table", table_type = "cell_vector_family_qc_key_summary"),
      output("pairwise_summary", "obs_table", table_type = "cell_vector_family_qc_pairwise_summary"),
      output("group_summary", "obs_table", table_type = "cell_vector_family_qc_group_summary"))),

  op("cell_vector_family_graph_smooth", "Graph Smooth Cell Vector Family",
    "Smooth signal keys over nearest-neighbor geometry with authenticated source ancestry.",
    ".run_family_graph_smooth", list(input("signal", "cell_vector_family"), input("space", "space")),
    schema(list(signal_keys = strings(), dims = integers(), k = integer(30L, 1L, 100000L, FALSE),
      alpha = number(.5, 0, 1, FALSE), distance_metric = string("euclidean", "euclidean", FALSE),
      query_chunk_size = integer(512L, 1L, 1000000L, FALSE), include_self = boolean(FALSE),
      na_policy = string("ignore", c("ignore", "zero", "propagate"), FALSE),
      source_signal_artifact_id = string(), source_run_id = string(), source_snapshot_id = string())),
    list(family_policy("smoothed_signal", "signal"), table_policy("smoothing_summary")),
    list(output("smoothed_signal", "cell_vector_family", "smoothing", value_kind = "continuous"),
      output("smoothing_summary", "obs_table", "smoothing", "cell_vector_family_graph_smooth_summary")),
    list(core, dependency("RANN", "nn2", spec = ">=2.6.0")),
    closure = list(smoothed_signal = c("smoothed_signal", "smoothing_summary"),
      smoothing_summary = c("smoothed_signal", "smoothing_summary")),
    provenance_reads = list(list(name = "source_ancestry",
      category = "authenticated_core_provenance_snapshot", kind = "artifact_ancestry",
      reader = "read_artifact_ancestry", input = "signal",
      activation = "when_resolved_from_artifact", lifecycle_policy = "require_current",
      parameter_assertions = list(artifact_id = "source_signal_artifact_id",
        run_id = "source_run_id", snapshot_id = "source_snapshot_id")))),

  op("cell_vector_winners", "Cell Vector Winner Assignment",
    "Assign deterministic winners and pruning diagnostics from a signal family.",
    ".run_cell_vector_winners", list(input("signal", "cell_vector_family")),
    schema(list(signal_keys = strings(), min_winners = integer(10L, 0L, 100000000L, FALSE),
      tie_method = string("first", c("first", "last"), FALSE),
      source_score_artifact_id = string(), source_run_id = string(), source_snapshot_id = string())),
    list(table_policy("winners")), list(output("winners", "obs_table", table_type = "cell_vector_winners"))),

  op("winner_cell_sets", "Winner Cell-Set Registration",
    "Request atomic core-owned admission of the complete winner-derived cell-set family.",
    ".run_winner_cell_sets", list(input("winners", "obs_table")),
    schema(list(source_winner_artifact_id = string(), cell_set_prefix = string("winner_", nullable = FALSE),
      winner_col = string("winner", nullable = FALSE), passes_col = string("passes_pruning", nullable = FALSE),
      source = string("cell_vector_winners", nullable = FALSE), created_by = string(), replace = boolean(TRUE))),
    list(table_policy("manifest")),
    list(output("manifest", "obs_table", "winner_family", "winner_cell_set_manifest"),
      output("cell_set_family", "cell_set_family_admission", "winner_family", category = "core_registry_request",
        profile = "cell_set_family_admission_v1")),
    closure = list(manifest = c("cell_set_family", "manifest")),
    replay = list(manifest = "exact_canonical_value",
      cell_set_family = "exact_family_membership_and_order")),

  op("feature_set_overlap_qc", "Feature-Set Matrix Overlap QC",
    "Screen authenticated feature-set-registry entries against a matrix universe.",
    ".run_feature_set_overlap_qc", list(input("matrix", "matrix")),
    schema(list(namespace = string(), namespaces = strings(), feature_sets = strings(),
      min_set_size = integer(1L, 0L, nullable = FALSE), max_set_size = number(minimum = 1),
      min_detected_fraction = number(0, 0, 1, FALSE), deduplicate_sets = boolean(TRUE))),
    list(table_policy("overlap_qc"), table_policy("namespace_summary")),
    list(output("overlap_qc", "obs_table", "overlap", "feature_set_overlap_qc"),
      output("namespace_summary", "obs_table", "overlap", "feature_set_overlap_qc_summary")),
    closure = list(overlap_qc = c("overlap_qc", "namespace_summary"),
      namespace_summary = c("overlap_qc", "namespace_summary")),
    registry_reads = list(list(name = "feature_sets", category = "authenticated_core_registry_snapshot",
      kind = "feature_set_registry", reader = "read_feature_set_registry", matrix_input = "matrix",
      namespace_parameter = "namespace", namespaces_parameter = "namespaces",
      feature_sets_parameter = "feature_sets", no_selector_policy = "all_namespaces",
      empty_match_policy = "error"))),

  op("pathway_enrich", "Pathway Enrichment",
    "Hypergeometric enrichment against an authenticated feature-set namespace snapshot.",
    ".run_pathway_enrich", list(input("matrix", "matrix"), input("feature_table", "feature_table")),
    schema(list(namespace = string(nullable = FALSE), feature_col = string("feature_id", nullable = FALSE),
      p_adj_col = string("p_val_adj", nullable = FALSE), p_adj_max = number(.05, 0, 1, FALSE),
      lfc_col = string("avg_log1p_fc", nullable = FALSE),
      direction = string("any", nullable = FALSE),
      lfc_min_abs = number(0, 0, nullable = FALSE), lfc_min = number(.25, 0, nullable = FALSE),
      min_set_size = integer(5L, 1L, nullable = FALSE), max_set_size = integer(5000L, 1L, nullable = FALSE),
      p_adjust_method = string("BH", nullable = FALSE)), "namespace"),
    list(table_policy("pathways")), list(output("pathways", "obs_table", table_type = "pathway_enrich")),
    registry_reads = list(list(name = "feature_sets", category = "authenticated_core_registry_snapshot",
      kind = "feature_set_registry", reader = "read_feature_set_registry", matrix_input = "matrix",
      namespace_parameter = "namespace", no_selector_policy = "all_namespaces",
      empty_match_policy = "error")))
)

manifest <- list(schema_version = "scworkbench.operation_pack.v3",
  pack_abi_version = "3", identity_scheme = "scworkbench.operation_pack.canonical_json.v1",
  namespace = "scwbgeneprograms", package = list(name = "scwbGenePrograms", version = "0.1.0"),
  backend = "R", compatibility = list(scworkbench = list(minimum = "0.1.0",
    maximum_exclusive = "0.2.0"), pack_abi = list(minimum = "3", maximum = "3")),
  entry_point = "scwbGenePrograms::register_operation_pack", operations = operations)
public_interfaces <- jsonlite::read_json(file.path(root, "tools", "public_interfaces.json"),
  simplifyVector = FALSE)
for (i in seq_along(manifest$operations)) {
  name <- manifest$operations[[i]]$name
  manifest$operations[[i]]$public_interface <- public_interfaces[[name]]
}
path <- file.path(root, "inst", "scworkbench", "operation-pack.json")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE,
  null = "null", digits = NA), path, useBytes = TRUE)
cat(path, "\n")
