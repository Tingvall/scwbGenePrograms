#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
option <- function(name) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) != 1L) stop("exactly one ", prefix, " argument is required")
  sub(prefix, "", hit, fixed = TRUE)
}
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(normalizePath(script_path)), "equivalence_contract.R"))

old <- readRDS(option("old"))
new <- readRDS(option("new"))
output <- option("output")
rows <- compare_geneprogram_equivalence(old, new, output)
cat("EQUIVALENCE OK methods=", nrow(rows), " outputs=",
  sum(rows$output_count), " parameters=effective-run-record",
  " inputs=identity-and-binding receipts=authenticated replay=declared\n", sep = "")
