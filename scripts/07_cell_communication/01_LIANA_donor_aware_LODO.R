# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 20_CCC_CROSS_METHOD_LOCKED_BENCHMARK/CODE/recalculate_liana_lodo_primary20.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
wd <- file.path(paths$work_root, "cell_communication")
out_dir <- file.path(paths$results_root, "cell_communication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
input_file <- gart_require_file(file.path(wd, "LIANA_per_donor_paired_differences.csv.gz"), "upstream donor-paired LIANA differences")
x <- as.data.table(read.csv(gzfile(input_file), check.names = FALSE))
x[,sender_cells:=fcase(source=="Myeloid",n_Myeloid,source=="Fibroblast",n_Fibroblast,source=="Endothelial",n_Endothelial,source=="Other epithelial",get("n_Other epithelial"),source=="T/NK",get("n_T/NK"),default=NA_real_)]
x<-x[gart_threshold==1 & expression_threshold==.10 & resource=="consensus" & ligand_complex=="MIF" & receptor_complex=="CD74_FAMILY" & sender_cells>=20 & get("n_GART-undetected ISC-like")>=20 & get("n_GART-detected ISC-like")>=20]
out<-x[,rbindlist(lapply(sort(unique(donor_id)),function(od){z<-difference[donor_id!=od];data.table(omitted_donor=od,median_difference=median(z),positive_direction=median(z)>0,n=length(z))})),by=.(method,source)]
fwrite(out, file.path(out_dir, "LIANA_MIF_CD74_LODO_primary20.csv"))
print(out[,.(runs=.N,positive_fraction=mean(positive_direction)),by=.(method,source)])
