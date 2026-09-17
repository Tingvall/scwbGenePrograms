# scwbGenePrograms

`scwbGenePrograms` is the explicitly activated scworkbench ABI-v3 operation
pack for the ten gene-program operations extracted during Core-to-Packs Phase
3B. Installation is inert. Activate it explicitly:

```r
library(scworkbench)
wb_operation_pack_load("scwbGenePrograms")
wb_list_ops()
```

It contributes `scwbgeneprograms::addmodulescore`,
`scwbgeneprograms::ucell`, `scwbgeneprograms::cell_vector_family_reduce`,
`scwbgeneprograms::cell_vector_family_calibrate`,
`scwbgeneprograms::cell_vector_family_qc`,
`scwbgeneprograms::cell_vector_family_graph_smooth`,
`scwbgeneprograms::cell_vector_winners`,
`scwbgeneprograms::winner_cell_sets`,
`scwbgeneprograms::feature_set_overlap_qc`, and
`scwbgeneprograms::pathway_enrich`. Unqualified historical names resolve only
when the pack is active; new provenance records the qualified identity.

Core authenticates every typed input and owns identifiers, staging,
publication, promotion, registry admission, provenance, replay, and recipient
validation. The pack receives only typed readers and sanctioned writers.
UCell reads authenticated bounded cell slices while the pack owns partitioning,
scoring, merge, and its final summary. Feature-set methods consume an
authenticated registry snapshot. Winner-derived cell sets are submitted as
one family request and applied atomically by core. Graph smoothing consumes
only core-attested artifact ancestry.

Promoted artifacts remain inspectable without this package. Prospective
execution and strict replay require the exact compatible pack identity.
