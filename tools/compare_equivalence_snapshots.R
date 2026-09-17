#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
option <- function(name) sub(paste0("^--", name, "="), "",
  args[startsWith(args, paste0("--", name, "="))])
old <- readRDS(option("old")); new <- readRDS(option("new"))
out <- option("output")
methods <- names(old)
if (!identical(methods, names(new)) || length(methods) != 10L) stop("method ledger mismatch")
canonical <- function(value) {
  attributes(value) <- attributes(value)[intersect(names(attributes(value)), "names")]
  if (is.list(value)) value <- lapply(value, canonical)
  value
}
rows <- lapply(methods, function(method) {
  old_outputs <- canonical(old[[method]]$outputs)
  new_outputs <- canonical(new[[method]]$outputs)
  ok <- isTRUE(all.equal(old_outputs, new_outputs, tolerance = 1e-12,
    check.attributes = FALSE))
  receipt <- new[[method]]$receipt
  data.frame(method = method, output_equivalent = ok,
    output_count = length(new_outputs),
    receipt_present = is.list(receipt),
    receipt_identity = as.character(receipt$receipt_identity$sha256 %||% ""),
    stringsAsFactors = FALSE)
})
rows <- do.call(rbind, rows)
utils::write.csv(rows, out, row.names = FALSE, quote = TRUE)
if (!all(rows$output_equivalent) || !all(rows$receipt_present)) {
  print(rows)
  stop("old-core/new-pack equivalence failed")
}
cat("EQUIVALENCE OK methods=10 outputs=", sum(rows$output_count), "\n", sep = "")
