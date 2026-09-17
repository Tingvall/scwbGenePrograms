#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
option <- function(name) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) != 1L) stop("exactly one ", prefix, " argument is required")
  sub(prefix, "", hit, fixed = TRUE)
}
mode <- option("mode")
output <- option("output")
if (!mode %in% c("old-core", "v3-pack")) stop("unsupported mode")

if (mode == "old-core") {
  devtools::load_all(option("core-source"), quiet = TRUE)
} else {
  library(scworkbench)
  library(scwbGenePrograms)
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbGenePrograms")
  wb_operation_pack_activate(manifest, lib.loc = .libPaths())
}

bundle_path <- tempfile("scwbGenePrograms-equivalence-")
b <- wb_init_bundle(bundle_path, overwrite = TRUE)
on.exit({
  try(wb_close_bundle(b), silent = TRUE)
  unlink(bundle_path, recursive = TRUE, force = TRUE)
}, add = TRUE)
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
  namespace = "programs", source_id = "fixture", recipe = list(n_genes_raw = 4L))
wb_feature_set_replace(b, "program_b", "RNA", features[3:6],
  namespace = "programs", source_id = "fixture", recipe = list(n_genes_raw = 4L))
signal <- wb_prim_cell_vector_family(rep(c("A", "B"), each = length(cells)),
  rep(cells, 2L), c(seq(.1, .8, .1), rev(seq(.1, .8, .1))), "continuous")
grouping <- wb_prim_cell_vector(stats::setNames(rep(c("g1", "g2"), each = 4L), cells),
  "categorical")
coords <- cbind(x = seq_along(cells), y = c(1,2,1,2,1,2,1,2))
rownames(coords) <- cells
wb_add_space(b, "fixture_space", coords, overwrite = TRUE)
space <- structure(list(kind = "space", space_id = "fixture_space", values = coords,
  dimensions = as.integer(dim(coords)), cell_id = cells,
  storage_profile = "space_dense_v1"), class = "wb_primitive")
feature_table <- wb_prim_feature_table(data.frame(feature_id = features[1:5],
  p_val_adj = rep(.01, 5L), avg_log1p_fc = rep(1, 5L)), "de_wilcox")
matrix_primitive <- structure(list(kind = "matrix", matrix_id = "RNA",
  values = values, dimensions = as.integer(dim(values)), cell_id = cells,
  feature_id = features, storage_profile = "matrix_ordinary_v1"),
  class = "wb_primitive")

transport <- if (mode == "v3-pack")
  getFromNamespace(".wb_operation_pack_params_transport", "scworkbench") else NULL
execute <- function(name, inputs, params) {
  if (mode == "v3-pack") {
    result <- getFromNamespace(".wb_operation_pack_v2_execute", "scworkbench")(
      b, paste0("scwbgeneprograms::", name), inputs, transport(params))
    return(list(outputs = result$outputs, params = result$request$params,
      receipt = result$receipt))
  }
  value <- switch(name,
    addmodulescore = wb_op_addmodulescore(b, "RNA", params$features,
      params$pool, cells, params$nbin, params$ctrl, params$name, params$seed),
    ucell = wb_op_ucell(b, "RNA", params$features, cells, params$name,
      params$max_rank, params$ncores, params$parallel_backend,
      params$chunk_size, params$missing_features_policy),
    cell_vector_family_reduce = wb_op_cell_vector_family_reduce(inputs$signal,
      params$signal_keys, params$rule, params$min_score, params$min_margin,
      params$unknown_label, params$ties, params$tie_tol),
    cell_vector_family_calibrate = wb_op_cell_vector_family_calibrate(inputs$signal,
      params$signal_keys, params$method, params$center, params$scale,
      params$ties_method, params$na_rm),
    cell_vector_family_qc = wb_op_cell_vector_family_qc(inputs$signal,
      params$signal_keys, inputs$grouping, params$correlation_methods,
      params$tie_tol),
    cell_vector_family_graph_smooth = wb_op_cell_vector_family_graph_smooth(b,
      inputs$signal, "fixture_space", cells, params$signal_keys, params$dims,
      params$k, params$alpha, params$distance_metric, params$query_chunk_size,
      params$include_self, params$na_policy, source_signal_artifact_id = NULL,
      source_run_id = NULL, source_snapshot_id = NULL),
    cell_vector_winners = wb_op_cell_vector_winners(inputs$signal,
      params$signal_keys, params$min_winners, params$tie_method,
      source_score_artifact_id = NULL, source_run_id = NULL,
      source_snapshot_id = NULL),
    feature_set_overlap_qc = wb_op_feature_set_overlap_qc(b, "RNA",
      namespace = params$namespace, namespaces = params$namespaces,
      feature_sets = params$feature_sets, min_set_size = params$min_set_size,
      max_set_size = params$max_set_size %||% Inf,
      min_detected_fraction = params$min_detected_fraction,
      deduplicate_sets = params$deduplicate_sets),
    pathway_enrich = wb_op_pathway_enrich(b, "RNA", inputs$feature_table$table,
      params$namespace, params$feature_col, params$p_adj_col, params$p_adj_max,
      params$lfc_col, params$direction, params$lfc_min_abs, params$lfc_min,
      params$min_set_size, params$max_set_size, params$p_adjust_method))
  list(outputs = value$outputs, params = value$params, receipt = NULL)
}

specs <- list(
  addmodulescore = list(inputs = list(matrix = matrix_primitive), params = list(
    features = list(A = features[1:2], B = features[3:4]), pool = features,
    nbin = 4L, ctrl = 2L, name = "Module", seed = 7L)),
  ucell = list(inputs = list(matrix = matrix_primitive), params = list(
    features = list(A = features[1:2], B = features[3:4]), name = "Signature",
    max_rank = 100L, ncores = 1L, parallel_backend = "serial", chunk_size = 3L,
    missing_features_policy = "drop")),
  cell_vector_family_reduce = list(inputs = list(signal = signal), params = list(
    signal_keys = NULL, rule = "argmax", min_score = NULL, min_margin = .1,
    unknown_label = "unknown", ties = "unknown", tie_tol = 1e-8)),
  cell_vector_family_calibrate = list(inputs = list(signal = signal), params = list(
    signal_keys = NULL, method = "percentile", center = TRUE, scale = TRUE,
    ties_method = "average", na_rm = TRUE)),
  cell_vector_family_qc = list(inputs = list(signal = signal, grouping = grouping),
    params = list(signal_keys = NULL, correlation_methods = c("pearson", "spearman"),
      tie_tol = 1e-8)),
  cell_vector_family_graph_smooth = list(inputs = list(signal = signal, space = space),
    params = list(signal_keys = NULL, dims = NULL, k = 2L, alpha = .5,
      distance_metric = "euclidean", query_chunk_size = 3L, include_self = FALSE,
      na_policy = "ignore", source_signal_artifact_id = NULL,
      source_run_id = NULL, source_snapshot_id = NULL)),
  cell_vector_winners = list(inputs = list(signal = signal), params = list(
    signal_keys = NULL, min_winners = 1L, tie_method = "first",
    source_score_artifact_id = NULL, source_run_id = NULL,
    source_snapshot_id = NULL)),
  feature_set_overlap_qc = list(inputs = list(matrix = matrix_primitive), params = list(
    namespace = "programs", namespaces = NULL, feature_sets = NULL,
    min_set_size = 1L, max_set_size = NULL, min_detected_fraction = 0,
    deduplicate_sets = TRUE)),
  pathway_enrich = list(inputs = list(matrix = matrix_primitive,
    feature_table = feature_table), params = list(namespace = "programs",
      feature_col = "feature_id", p_adj_col = "p_val_adj", p_adj_max = .05,
      lfc_col = "avg_log1p_fc", direction = "any", lfc_min_abs = 0,
      lfc_min = .25, min_set_size = 1L, max_set_size = 100L,
      p_adjust_method = "BH")))

result <- lapply(names(specs), function(name) execute(name,
  specs[[name]]$inputs, specs[[name]]$params))
names(result) <- names(specs)
winners <- result$cell_vector_winners$outputs$winners
if (mode == "old-core") {
  manifest <- wb_register_winner_cell_sets(b, winners,
    cell_set_prefix = "winner_", replace = TRUE)
  result$winner_cell_sets <- list(outputs = list(manifest = manifest),
    params = list(source_winner_artifact_id = NULL, cell_set_prefix = "winner_",
      winner_col = "winner", passes_col = "passes_pruning",
      source = "cell_vector_winners", created_by = NULL, replace = TRUE),
    receipt = NULL)
} else {
  result$winner_cell_sets <- execute("winner_cell_sets", list(winners = winners),
    list(source_winner_artifact_id = NULL, cell_set_prefix = "winner_",
      winner_col = "winner", passes_col = "passes_pruning",
      source = "cell_vector_winners", created_by = NULL, replace = TRUE))
}
result <- result[c("addmodulescore", "ucell", "cell_vector_family_reduce",
  "cell_vector_family_calibrate", "cell_vector_family_qc",
  "cell_vector_family_graph_smooth", "cell_vector_winners", "winner_cell_sets",
  "feature_set_overlap_qc", "pathway_enrich")]
saveRDS(result, output, version = 3L)
cat("SNAPSHOT OK mode=", mode, " methods=", length(result), " output=", output, "\n", sep = "")
