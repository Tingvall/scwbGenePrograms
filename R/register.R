#' Register the scwbGenePrograms operation pack
#'
#' Installation is inert. scworkbench invokes this zero-argument entry point
#' only at its explicit operation-pack activation boundary.
#'
#' @return A validated scworkbench operation-pack registration.
#' @export
register_operation_pack <- function() {
  scworkbench::wb_operation_pack_register(
    system.file("scworkbench", "operation-pack.json", package = "scwbGenePrograms"),
    list(
      `scwbgeneprograms::addmodulescore` = .run_addmodulescore,
      `scwbgeneprograms::ucell` = .run_ucell,
      `scwbgeneprograms::cell_vector_family_reduce` = .run_family_reduce,
      `scwbgeneprograms::cell_vector_family_calibrate` = .run_family_calibrate,
      `scwbgeneprograms::cell_vector_family_qc` = .run_family_qc,
      `scwbgeneprograms::cell_vector_family_graph_smooth` = .run_family_graph_smooth,
      `scwbgeneprograms::cell_vector_winners` = .run_cell_vector_winners,
      `scwbgeneprograms::winner_cell_sets` = .run_winner_cell_sets,
      `scwbgeneprograms::feature_set_overlap_qc` = .run_feature_set_overlap_qc,
      `scwbgeneprograms::pathway_enrich` = .run_pathway_enrich
    )
  )
}
