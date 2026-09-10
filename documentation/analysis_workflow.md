# Analysis workflow

1. Download public bulk, single-cell, and spatial datasets listed in `data_public/public_dataset_manifest.csv`.
2. Run bulk cohort preparation and the four-cohort Hedges' g REML/Hartung-Knapp synthesis.
3. Reconstruct the corrected single-cell metadata from frozen public objects, preserving raw-count GART detection (`>=1`; sensitivity `>=2`) and donor/sample identifiers. Supplementary Table S2 separates the complete 145-entry atlas map from the 79-entry ISC-like analysis subset. Historical upstream QC, doublet, ambient-RNA, and removal-count parameters that are unavailable from archived records are not inferred; no additional upstream filtering was performed during revision.
4. Fit donor-aware GART models; run donor-paired pseudobulk and enrichment analyses in 11 paired SCP259 UC donors, broader-atlas public-gene-set module and purine-enzyme analyses, pseudotime sensitivity, and five-run scTenifoldKnk workflows.
5. Summarize final donor-aware CellChat and LIANA evidence. Four predefined primary sender compartments are myeloid, fibroblast, endothelial, and T/NK; other epithelial cells excluding the ISC-like receiver are evaluated as sensitivity analyses. LIANA+ uses 100 internal permutations for internal cell-level scoring, while formal inference uses paired donor-level contrasts. Communication analyses are candidate-generating.
6. For GSE189184, apply the independent single-cell-derived program inside the locked within-section top-40% epithelial ROI and perform donor-level paired inference.
7. Generate figures exclusively from the files under `data_public/figure_source_data`.

Large public expression matrices and frozen object files are not duplicated; accessions and download links are supplied. Path-refactored public scripts retain the archived final statistical logic and parameters.
