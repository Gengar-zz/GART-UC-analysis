# Exclusions and locked decisions

- The final bulk synthesis contains GSE24287, GSE36807, GSE37283, and GSE87466. GSE3629 lacks a healthy comparator and is excluded.
- In GSE37283, the 11 UC-associated-neoplasia samples are excluded; five normal controls and four quiescent non-neoplastic UC samples are retained.
- Single-cell preprocessing is dataset-specific rather than a newly rerun unified QC pipeline. The GSE214695 source workflow used Cell Ranger-filtered matrices, >100 detected genes, <65% mitochondrial transcripts, no additional hard UMI cutoff, and per-sample scDblFinder 1.8.0 with default parameters. The processed SCP259 release used >=250 detected genes, iterative mixed-lineage doublet-cluster removal with reassignment/reclustering, and source-study ambient-RNA contamination filtering.
- GART detection is defined from corrected raw counts: 0 is undetected and >=1 is detected; >=2 is sensitivity only.
- GART is excluded from the purine-biosynthesis validation score.
- Pseudotime roots were selected independently of GART. Directional transition language is not used.
- scTenifoldKnk is exploratory in silico perturbation and is not functional validation. No explicit permutation-based null procedure was used.
- Final communication evidence is donor-aware CellChat/LIANA plus the MYBL2 donor-aware expression summary. Pooled-only communication, causal-axis claims, old SCENIC activation claims, and NicheNet negative exploratory outputs are not included.
- The final spatial program is `SC_GART_NEAREST_CENTROID_PROGRAM`; the locked ROI is the within-section top 40% of the EPCAM/KRT8/KRT18/KRT19 score. ROI construction excludes GART, MYBL2, MKI67, PCNA, PAICS, OLFM4, MIF, CD74, and disease outcome.
- Earlier outcome-informed ROI and legacy spatial composite scores/results: **Superseded during revision; not used in the final analysis.** No legacy source data are distributed.
- Publication-facing Table S7 contains the de-identified, author-confirmed participant characteristics available from source records. Unavailable formal matching and other unrecoverable metadata are not inferred; no direct identifiers or private source records are distributed.
