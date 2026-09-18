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
invisible(compare_geneprogram_equivalence(old, new))

expect_code <- function(label, candidate, code) {
  observed <- tryCatch({
    compare_geneprogram_equivalence(old, candidate)
    NULL
  }, geneprogram_equivalence_error = function(error) error)
  if (is.null(observed) || !identical(observed$code, code)) stop(
    label, " did not fail with ", code, "; observed=",
    if (is.null(observed)) "success" else observed$code)
  cat("NEGATIVE PASS fixture=", label, " code=", code, "\n", sep = "")
}

changed_parameters <- new
changed_parameters$methods$addmodulescore$run$params$module_definitions$B <- NULL
expect_code("identical_outputs_altered_derived_module_count", changed_parameters,
  "EQUIVALENCE_PARAMETER_MISMATCH")

changed_binding <- new
changed_binding$methods$cell_vector_family_reduce$input_identities$content_sha256[[1L]] <-
  paste(rep("0", 64L), collapse = "")
expect_code("identical_outputs_altered_input_identity", changed_binding,
  "EQUIVALENCE_INPUT_BINDING_MISMATCH")

changed_policy <- new
changed_policy$methods$pathway_enrich$outputs$replay_equivalence_policy[[1L]] <-
  "not_comparable"
expect_code("identical_outputs_altered_terminal_policy", changed_policy,
  "EQUIVALENCE_TERMINAL_POLICY_MISMATCH")

cat("NEGATIVE EQUIVALENCE OK fixtures=3 intended_failures=3\n")
