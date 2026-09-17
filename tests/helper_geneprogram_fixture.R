geneprogram_activate_pack <- function() {
  library(scworkbench)
  if ("scwbgeneprograms::ucell" %in% wb_list_ops()) {
    return(invisible(TRUE))
  }
  path <- system.file("scworkbench", "operation-pack.json", package = "scwbGenePrograms")
  if (!nzchar(path)) stop("scwbGenePrograms is not installed")
  wb_operation_pack_activate(path, lib.loc = .libPaths())
  invisible(TRUE)
}

geneprogram_fixture <- function() {
  geneprogram_activate_pack()
  bundle <- wb_init_bundle(tempfile("scwbGenePrograms-test-"))
  cells <- paste0("cell_", sprintf("%02d", 1:8))
  features <- paste0("gene_", 1:8)
  wb_define_cells(bundle, cells)
  values <- matrix(c(
    9,8,7,6,1,0,1,0, 8,7,6,5,0,1,0,1,
    1,0,1,0,9,8,7,6, 0,1,0,1,8,7,6,5,
    2,3,2,3,2,3,2,3, 3,2,3,2,3,2,3,2,
    4,4,3,3,2,2,1,1, 1,1,2,2,3,3,4,4),
    nrow = 8L, byrow = TRUE, dimnames = list(features, cells))
  wb_attach_matrix(bundle, "RNA", Matrix::Matrix(values, sparse = TRUE),
    feature_type = "gene")
  wb_feature_set_replace(bundle, "program_a", "RNA", features[1:3],
    namespace = "programs", source_id = "fixture", recipe = list(n_genes_raw = 4L))
  wb_feature_set_replace(bundle, "program_b", "RNA", features[3:6],
    namespace = "programs", source_id = "fixture", recipe = list(n_genes_raw = 4L))
  matrix_primitive <- structure(list(kind = "matrix", matrix_id = "RNA",
    values = values, dimensions = as.integer(dim(values)), cell_id = cells,
    feature_id = features, storage_profile = "matrix_ordinary_v1"),
    class = "wb_primitive")
  signal <- wb_prim_cell_vector_family(rep(c("A", "B"), each = length(cells)),
    rep(cells, 2L), c(seq(.1, .8, .1), rev(seq(.1, .8, .1))), "continuous")
  grouping <- wb_prim_cell_vector(stats::setNames(rep(c("g1", "g2"), each = 4L), cells),
    "categorical")
  coords <- cbind(x = seq_along(cells), y = c(1,2,1,2,1,2,1,2))
  rownames(coords) <- cells
  space <- structure(list(kind = "space", space_id = "fixture_space", values = coords,
    dimensions = as.integer(dim(coords)), cell_id = cells,
    storage_profile = "space_dense_v1"), class = "wb_primitive")
  feature_table <- wb_prim_feature_table(data.frame(feature_id = features[1:5],
    p_val_adj = rep(.01, 5L), avg_log1p_fc = rep(1, 5L)), "de_wilcox")
  list(bundle = bundle, cells = cells, features = features, values = values,
    matrix = matrix_primitive, signal = signal, grouping = grouping,
    space = space, feature_table = feature_table)
}

geneprogram_execute <- function(fixture, name, inputs, params = list()) {
  transport <- getFromNamespace(".wb_operation_pack_params_transport", "scworkbench")
  getFromNamespace(".wb_operation_pack_v2_execute", "scworkbench")(
    fixture$bundle, paste0("scwbgeneprograms::", name), inputs,
    transport(params))
}

geneprogram_output <- function(result, name) {
  readRDS(file.path(result$path, "outputs", name, "payload.rds"))
}
