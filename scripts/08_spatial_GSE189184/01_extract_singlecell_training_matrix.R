# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 19_SPATIAL_REVISED_PRIMARY_ADJUDICATION/code/02_extract_singlecell_matrix.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(Matrix))

root <- file.path(paths$work_root, "spatial_GSE189184")
obj_path <- gart_require_file(file.path(paths$external_data_root, "singlecell", "frozen_objects", "UC_GART_07_ISC_like_subclustered_FIXED_withDonor.rds"), "frozen ISC-like Seurat object")
meta_path <- gart_require_file(file.path(paths$results_root, "singlecell", "INTERMEDIATE", "ISC_corrected_cell_metadata.csv"), "upstream corrected ISC metadata")

x <- readRDS(obj_path)
assay <- attr(x, "assays")[["RNA"]]
layers <- attr(assay, "layers")
feature_map <- attr(assay, "features"); attr(feature_map, "class") <- NULL
cell_map <- attr(assay, "cells"); attr(cell_map, "class") <- NULL
layer_name <- "data.SCP_Epi.2"
mat <- layers[[layer_name]]
rownames(mat) <- rownames(feature_map)[feature_map[, layer_name]]
colnames(mat) <- rownames(cell_map)[cell_map[, layer_name]]

meta <- read.csv(meta_path, check.names=FALSE)
paired <- with(meta[meta$cohort == "SCP259", ], table(donor_id, GART_status))
paired_donors <- rownames(paired)[apply(paired, 1, function(z) all(z[c("negative", "positive")] >= 20))]
stopifnot(length(paired_donors) == 11L)
keep_meta <- meta$cohort == "SCP259" & meta$donor_id %in% paired_donors
meta <- meta[keep_meta, ]
common <- intersect(meta$cell_id_revision, colnames(mat))
stopifnot(length(common) == nrow(meta))
meta <- meta[match(common, meta$cell_id_revision), ]
mat <- mat[, common, drop=FALSE]

Matrix::writeMM(mat, file.path(root, "source_data", "singlecell_SCP259_11donor_logexpr.mtx"))
writeLines(rownames(mat), file.path(root, "source_data", "singlecell_SCP259_gene_names.txt"), useBytes=TRUE)
writeLines(colnames(mat), file.path(root, "source_data", "singlecell_SCP259_cell_names.txt"), useBytes=TRUE)
write.csv(meta, file.path(root, "source_data", "singlecell_SCP259_11donor_metadata.csv"), row.names=FALSE, quote=TRUE)
write.csv(data.frame(donor_id=paired_donors, n_negative=paired[paired_donors,"negative"], n_positive=paired[paired_donors,"positive"]), file.path(root, "source_data", "singlecell_training_donor_counts.csv"), row.names=FALSE, quote=TRUE)
cat("matrix", nrow(mat), "x", ncol(mat), "nnz", length(mat@x), "donors", length(paired_donors), "\n")
