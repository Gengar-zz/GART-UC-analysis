# External public and archived derived inputs

Large public matrices and local frozen objects are not committed. By default, scripts resolve this directory as the external-data root. Set `GART_UC_DATA_ROOT` to an alternative root if desired. Missing inputs produce an error that points back to this file. Generated intermediates go to `work/`; regenerated outputs go to `results_reproduced/`; scripts never overwrite `data_public/`.

## Expected layout

- `bulk/GSE24287/`, `bulk/GSE36807/`, and `bulk/GSE37283/`: public GEO downloads from [GSE24287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24287), [GSE36807](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36807), and [GSE37283](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37283). The repaired 21-script subset does not impose invented raw filenames for these accessions.
- `bulk/GSE87466/GSE87466_family.soft.gz`: [official GEO family SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87466/soft/GSE87466_family.soft.gz), consumed by `scripts/01_bulk/02_process_GSE87466.py`.
- `bulk/legacy_three_cohort/Fig1A_bulk_source_data.xlsx` and `bulk/legacy_three_cohort/bulk_meta_results_revised.csv`: archived derived provenance inputs used by the original bulk assembly scripts. These are not public raw data and are not required to inspect the committed final numerical source tables.
- `singlecell/GSE214695/`: public files from [GSE214695](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214695). No additional exact raw filename is imposed by the repaired scripts.
- `singlecell/SCP259/`: public files from [SCP259](https://singlecell.broadinstitute.org/single_cell/study/SCP259). No additional exact raw filename is imposed by the repaired scripts.
- `singlecell/frozen_objects/UC_GART_06c_Epithelial_Fig2AB_noGARTname_FINAL.rds`, `UC_GART_07_ISC_like_subclustered_FIXED_withDonor.rds`, and `CORRECTED_SCRNA_ANALYSIS_OBJECT.rds`: frozen upstream-derived analysis objects used by the final revision scripts.
- `singlecell/frozen_objects/RAW_LAYER_RECOVERY_AUDIT.csv`, `MODULE_GENE_SETS_AND_PARAMETERS.csv`, `PSEUDOTIME_CORRECTED_CELL_METADATA.rds`, and `PSEUDOTIME_INDEPENDENT_ROOT_SELECTION.csv`: exact archived correction inputs used by the public scripts.
- `singlecell/scTenifold_runs/Figure5_GART_virtual_KO_5seed_FINAL/`: archived five-run scTenifoldKnk outputs consumed by the consensus script.
- `spatial/GSE189184/`: public material from [GSE189184](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189184). The final scorer additionally expects `frozen_objects/GSE189184_ST_with_FINAL_Fig7_scores_STABLE.rds` and `spatial_optimization_spot_source.csv`, exact frozen upstream-derived inputs named by the authoritative script.
- `gene_sets/msigdb/h.all.v2025.1.Hs.symbols.gmt`, `c5.go.bp.v2025.1.Hs.symbols.gmt`, and `c2.cp.reactome.v2025.1.Hs.symbols.gmt`: MSigDB 2025.1.Hs files used by the donor-paired GSEA and spatial model-selection scripts. Their recovered direct URLs are recorded in `INPUT_CONTRACT_AUDIT.tsv`.
- `derived_inputs/singlecell/multiverse_singlecell_results.csv`: locked archived multiverse comparator used by donor-balanced validation.
- `derived_inputs/figure4/FigS4E_full_recurrent_gene_source_data.csv` and `FigS4F_full_programs_source_data.csv`: optional archived supplementary-render inputs; they are not required for final main panels.
- `fonts/arial.ttf` and `fonts/arialbd.ttf`: locally licensed Arial font files used by the Python Figure 4E renderer; these are not distributed by this repository.

See `INPUT_CONTRACT_AUDIT.tsv` for script-level categories, producers, requirements, and resolution status.
