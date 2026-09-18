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
cat("POSITIVE PASS fixture=unaltered_authenticated_transition methods=10 outputs=22\n")

hash <- getFromNamespace(".wb_operation_pack_sha256", "scworkbench")
zero_hash <- paste(rep("0", 64L), collapse = "")

rehash_field <- function(value, field) {
  identity <- value[[field]]
  value[[field]] <- NULL
  identity$sha256 <- hash(value)
  value[[field]] <- identity
  value
}

request_input_index <- function(candidate, method, input_name) {
  which(vapply(candidate$methods[[method]]$transition$request$inputs,
    function(x) identical(as.character(x$name), input_name), logical(1)))[[1L]]
}

receipt_output_index <- function(candidate, method, output_name) {
  which(vapply(candidate$methods[[method]]$transition$receipt$outputs,
    function(x) identical(as.character(x$name), output_name), logical(1)))[[1L]]
}

mutate_request_input <- function(candidate, method, input_name) {
  index <- request_input_index(candidate, method, input_name)
  descriptor <- candidate$methods[[method]]$transition$request$inputs[[index]]$descriptor
  descriptor$content_identity$sha256 <- zero_hash
  descriptor <- rehash_field(descriptor, "descriptor_identity")
  candidate$methods[[method]]$transition$request$inputs[[index]]$descriptor <-
    descriptor
  candidate$methods[[method]]$transition$request_identity$sha256 <- hash(
    candidate$methods[[method]]$transition$request)
  candidate
}

mutate_receipt_output <- function(candidate, method, output_name, field) {
  index <- receipt_output_index(candidate, method, output_name)
  descriptor <- candidate$methods[[method]]$transition$receipt$outputs[[index]]$descriptor
  descriptor[[field]]$sha256 <- zero_hash
  descriptor <- rehash_field(descriptor, "observed_descriptor_identity")
  candidate$methods[[method]]$transition$receipt$outputs[[index]]$descriptor <-
    descriptor
  candidate$methods[[method]]$transition$receipt_identity$sha256 <- hash(
    candidate$methods[[method]]$transition$receipt)
  candidate
}

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

for (method in c("addmodulescore", "ucell")) {
  changed <- mutate_request_input(new, method, "matrix")
  expect_code(paste0(method, "_self_consistent_request_input_content_tamper"),
    changed, "EQUIVALENCE_INPUT_BINDING_MISMATCH")
}

for (method in c("addmodulescore", "pathway_enrich")) {
  output <- if (identical(method, "addmodulescore")) "module_scores" else
    "pathways"
  changed <- mutate_receipt_output(new, method, output, "content_identity")
  expect_code(paste0(method, "_self_consistent_receipt_output_content_tamper"),
    changed, "EQUIVALENCE_OUTPUT_BINDING_MISMATCH")
  changed <- mutate_receipt_output(new, method, output,
    "payload_manifest_identity")
  expect_code(paste0(method, "_self_consistent_receipt_payload_tamper"),
    changed, "EQUIVALENCE_OUTPUT_PAYLOAD_BINDING_MISMATCH")
}

missing_input <- new
missing_input$methods$addmodulescore$transition$request$inputs <- list()
missing_input$methods$addmodulescore$transition$request_identity$sha256 <- hash(
  missing_input$methods$addmodulescore$transition$request)
expect_code("missing_request_input_descriptor", missing_input,
  "EQUIVALENCE_INPUT_BINDING_MISMATCH")

duplicate_input <- new
duplicate_input$methods$addmodulescore$transition$request$inputs[[2L]] <-
  duplicate_input$methods$addmodulescore$transition$request$inputs[[1L]]
duplicate_input$methods$addmodulescore$transition$request_identity$sha256 <- hash(
  duplicate_input$methods$addmodulescore$transition$request)
expect_code("duplicate_request_input_descriptor", duplicate_input,
  "EQUIVALENCE_INPUT_BINDING_MISMATCH")

missing_output <- new
missing_output$methods$addmodulescore$transition$receipt$outputs <- list()
missing_output$methods$addmodulescore$transition$receipt_identity$sha256 <- hash(
  missing_output$methods$addmodulescore$transition$receipt)
expect_code("missing_receipt_output_descriptor", missing_output,
  "EQUIVALENCE_OUTPUT_BINDING_MISMATCH")

duplicate_output <- new
duplicate_output$methods$addmodulescore$transition$receipt$outputs[[2L]] <-
  duplicate_output$methods$addmodulescore$transition$receipt$outputs[[1L]]
duplicate_output$methods$addmodulescore$transition$receipt_identity$sha256 <- hash(
  duplicate_output$methods$addmodulescore$transition$receipt)
expect_code("duplicate_receipt_output_descriptor", duplicate_output,
  "EQUIVALENCE_OUTPUT_BINDING_MISMATCH")

missing_payload_record <- new
missing_payload_record$methods$addmodulescore$published_payloads <- list()
expect_code("missing_authoritative_payload_record", missing_payload_record,
  "EQUIVALENCE_OUTPUT_PAYLOAD_BINDING_MISMATCH")

missing_payload_identity <- new
missing_payload_identity$methods$addmodulescore$published_payloads[[1L]]$
  payload_manifest_identity <- NULL
expect_code("absent_authoritative_payload_identity", missing_payload_identity,
  "EQUIVALENCE_OUTPUT_PAYLOAD_BINDING_MISMATCH")

cat("NEGATIVE EQUIVALENCE OK fixtures=15 intended_failures=15\n")
