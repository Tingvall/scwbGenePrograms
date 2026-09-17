args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: extract_public_interfaces.R <core-source> <output>")
pkgload::load_all(normalizePath(args[[1]], mustWork = TRUE), quiet = TRUE)
names <- c("addmodulescore", "ucell", "cell_vector_family_reduce",
  "cell_vector_family_calibrate", "cell_vector_family_qc",
  "cell_vector_family_graph_smooth", "cell_vector_winners", "winner_cell_sets",
  "feature_set_overlap_qc", "pathway_enrich")
catalog <- getFromNamespace(".wb_op_catalog", "scworkbench")()
`%||%` <- function(x, y) if (is.null(x)) y else x
map_type <- function(type) {
  mapping <- c(string_vec_or_list = "list", logical = "bool", int_vec = "list")
  index <- match(type, names(mapping))
  if (is.na(index)) type else unname(mapping[[index]])
}
interfaces <- lapply(names, function(name) {
  op <- catalog[[name]]
  list(
    intent = unname(as.list(op$intent)),
    inputs = unname(lapply(names(op$inputs), function(input_name) {
      input <- op$inputs[[input_name]]
      list(name = input_name, label = as.character(input$label %||% input_name),
        help = as.character(input$help %||% paste("Typed", input_name, "input.")))
    })),
    parameters = unname(lapply(op$params, function(parameter) list(
      name = parameter$name, type = map_type(parameter$type), help = parameter$help))),
    outputs = unname(lapply(names(op$interface_produces), function(output_name)
      list(name = output_name,
        description = op$interface_produces[[output_name]]$description)))
  )
})
names(interfaces) <- names
interfaces$winner_cell_sets$outputs[[length(interfaces$winner_cell_sets$outputs) + 1L]] <-
  list(name = "cell_set_family",
    description = "Atomic core-owned admission request for the complete winner-derived cell-set family.")
jsonlite::write_json(interfaces, normalizePath(args[[2]], mustWork = FALSE),
  auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
