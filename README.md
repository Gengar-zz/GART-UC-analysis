# GART expression is associated with a proliferative and purine-metabolic ISC-like epithelial state in ulcerative colitis

## Purpose

This repository contains analysis code, recovered environment records, locked parameters, and figure-level numerical source data for the final revision. Computational results are associative or hypothesis-generating unless supported by direct experiments; the repository does not establish a causal GART-mediated repair mechanism or a causal MIF/CD74–MYBL2–GART axis.

## Public datasets

- Bulk intestinal-mucosal transcriptomics: [GSE24287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24287), [GSE36807](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36807), [GSE37283](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37283), and [GSE87466](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87466).
- Single-cell RNA-seq: [GSE214695](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214695) and [SCP259](https://singlecell.broadinstitute.org/single_cell/study/SCP259).
- Spatial transcriptomics: [GSE189184](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189184).

GSE37283 retains five normal and four quiescent non-neoplastic UC samples and excludes eleven UC-associated-neoplasia samples. GSE3629 is not part of the final analysis because it has no healthy comparator.

## Repository structure

- `scripts/`: analysis and plotting scripts grouped by stage.
- `data_public/figure_source_data/`: numerical source data for main and supplementary figures.
- `data_public/supplementary_tables/`: exact publication-facing Tables S1–S9.
- `docs/` and `documentation/`: parameters, workflow, mappings, and data definitions.
- `environment/` and `session_info/`: recovered software versions and session records.
- `results_reference/`: locked statistics used for regression checks.

## Statistical units

Independent donors are the biological units for bulk, single-cell, and spatial inference; participants are the units for human tissue imaging; mice are the units for animal experiments. Cells, spots, repeated sample entries, and fields of view are not treated as independent biological replicates.

Participant- and mouse-level imaging values in the rebuilt public source-data panels are reported to two decimal places. Figure 9B retains the formatting of the previously locked public release. Statistical analyses were performed using the original unrounded archived values.

## Analysis workflow and figure reproduction

The stage order is bulk; single-cell preprocessing; GART detection/donor models; donor-paired pseudobulk, GSEA, and modules; pseudotime; scTenifoldKnk; donor-aware communication; spatial projection; figure generation; and animal endpoint verification. Panel-level input/output mappings are in `PUBLIC_SOURCE_DATA_TRACEABILITY.tsv` and `documentation/figure_to_script_map.csv`. Large expression matrices must be downloaded from the public repositories above and placed in repository-relative external-input locations described by the scripts. Figure 7–9 numerical source files reproduce de-identified participant- or mouse-level plotted values; original microscopy is not distributed.

All public scripts now use the repository-relative path contract documented in `data_external/README.md`. Large public or frozen upstream inputs are not bundled; missing inputs fail with a clear instruction. Generated intermediates and regenerated results are isolated from the locked `data_public/` directory.

## Software and random seeds

Recovered R package versions are in `environment/R_package_versions.csv`; known locked versions include R 4.5.2, Seurat 5.4.0, Harmony 1.2.4, CellChat 2.2.0.9001, LIANA+ 1.8.1, and scTenifoldKnk 1.0.3. Python package names are in `requirements.txt`; unavailable historical versions are explicitly marked rather than inferred. Recovered seeds are listed in `environment/random_seeds.csv`.

## Source-data notes

GART is excluded from the purine-biosynthesis validation module. scTenifoldKnk recurrence is a cross-run reproducibility descriptor, not permutation significance, and no explicit permutation procedure was used for gene-level testing. The primary spatial endpoint uses seven independent donors and an outcome-independent within-section top-40% epithelial ROI based on EPCAM/KRT8/KRT18/KRT19.

## Privacy and microscopy availability

Only de-identified numerical human and animal source data are included. Raw human microscopy, direct clinical identifiers, original Prism files, internal verification workbooks, and local QA files are excluded. Original microscopy is available from the corresponding author upon reasonable request, subject to institutional and ethical requirements.

## Citation and contact

Please cite the associated manuscript; machine-readable citation metadata are in `CITATION.cff`. Contact the corresponding author listed in the manuscript for data-access questions.
