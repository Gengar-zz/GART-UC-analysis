# scTenifoldKnk method provenance

The deposited R script consolidates the five archived final runs and applies the existing recurrence criteria. The original single end-to-end run-generation script was not recoverable as a complete five-run file; therefore, no call or parameter has been reconstructed by guesswork. Complete run-level gene tables are deposited in `data_public/figure_source_data/Figure4/scTenifold_5run_gene_level_complete_original_merged.csv`.

Nominal gene-level P values use the software's chi-square(1) reference for the standardized squared manifold displacement. Within-run adjusted P values use Benjamini-Hochberg FDR. No explicit permutation procedure was used. Cross-run top-100 recurrence is a reproducibility descriptor, not empirical permutation significance.
