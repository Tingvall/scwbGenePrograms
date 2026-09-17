`%||%` <- function(x, y) if (is.null(x)) y else x

.gp_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.gp_require <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for %s.", package, purpose), call. = FALSE)
  }
}

.gp_kind <- function(x, kind) {
  .gp_assert(inherits(x, "wb_primitive") && identical(as.character(x$kind)[1], kind),
    paste0("Expected a ", kind, " primitive."))
  invisible(TRUE)
}

.gp_groups <- function(features, name) {
  if (is.character(features)) features <- list(module1 = features)
  .gp_assert(is.list(features) && length(features),
    "features must be a character vector or non-empty named list.")
  features <- lapply(features, function(x) unique(as.character(x)))
  if (is.null(names(features)) || any(!nzchar(names(features)))) {
    names(features) <- paste0(as.character(name)[1], seq_along(features))
  }
  features
}

.gp_signal <- function(signal, signal_keys = NULL, cells = NULL) {
  .gp_kind(signal, "cell_vector_family")
  keys <- unique(as.character(signal$key))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (!is.null(signal_keys)) keys <- intersect(as.character(signal_keys), keys)
  .gp_assert(length(keys) > 0L, "No signal_keys remain after filtering.")
  rows <- as.character(signal$key) %in% keys
  if (!is.null(cells)) rows <- rows & as.character(signal$cell_id) %in% as.character(cells)
  df <- data.frame(cell_id = as.character(signal$cell_id[rows]),
    key = as.character(signal$key[rows]), value = as.numeric(signal$value[rows]),
    stringsAsFactors = FALSE)
  .gp_assert(nrow(df) > 0L, "No rows remain after signal filtering.")
  .gp_assert(!anyDuplicated(df[c("cell_id", "key")]),
    "signal contains duplicate cell_id/key pairs; expected one value per pair.")
  cell_order <- unique(df$cell_id)
  key_order <- keys[keys %in% unique(df$key)]
  mat <- matrix(NA_real_, length(cell_order), length(key_order),
    dimnames = list(cell_order, key_order))
  mat[cbind(match(df$cell_id, cell_order), match(df$key, key_order))] <- df$value
  list(matrix = mat, cells = cell_order, keys = key_order, data = df)
}

.gp_family <- function(mat, value_type = "continuous", meta = list()) {
  keys <- colnames(mat); cells <- rownames(mat)
  scworkbench::wb_prim_cell_vector_family(
    key = rep(keys, each = length(cells)),
    cell_id = rep(cells, times = length(keys)),
    value = as.numeric(mat[cbind(rep(cells, times = length(keys)),
      rep(keys, each = length(cells)))]),
    value_type = value_type, meta = meta)
}

.gp_registry_sets <- function(snapshot) {
  .gp_assert(is.list(snapshot) && identical(snapshot$schema,
    "scworkbench.feature_set_registry.snapshot.v1"),
    "Authenticated feature-set-registry snapshot is invalid.")
  lapply(snapshot$sets %||% list(), function(item) {
    members <- vapply(item$memberships %||% list(), `[[`, character(1), "feature_id")
    item$member_ids <- members
    item
  })
}

.gp_scalar <- function(x, default = "") {
  x <- x %||% default
  if (!length(x) || is.na(x[[1]])) return(as.character(default)[1])
  as.character(x[[1]])
}

.gp_write <- function(writers, name, value) {
  writers$write_terminal_primitive(name, value)
  invisible(value)
}
