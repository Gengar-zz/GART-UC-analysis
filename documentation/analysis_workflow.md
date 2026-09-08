# Analysis workflow

1. Download public bulk, single-cell, and spatial datasets listed in `data_public/public_dataset_manifest.csv`.
2. Run bulk cohort preparation and the four-cohort Hedges' g REML/Hartung-Knapp synthesis.
3. Reconstruct the corrected single-cell metadata from frozen public objects, preserving raw-count GART detection (`>=1`; sensitivity `>=2`) and donor/sample identifiers.
4. Fit donor-aware GART models; run donor-paired pseudobulk, public-gene-set module, purine-enzyme, pseudotime-sensitivity, and five-run scTenifoldKnk workflows.
5. Summarize final donor-aware CellChat and LIANA evidence. Communication analyses are candidate-generating.
6. For GSE189184, apply the independent single-cell-derived program inside the locked within-section top-40% epithelial ROI and perform donor-level paired inference.
7. Generate figures exclusively from the files under `data_public/figure_source_data`.

Large public expression matrices and frozen object files are not duplicated; accessions and download links are supplied. Path-refactored public scripts retain the archived final statistical logic and parameters.
