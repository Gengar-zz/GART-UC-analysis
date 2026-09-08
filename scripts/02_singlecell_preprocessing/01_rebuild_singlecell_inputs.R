# Public-release copy: file paths were converted to repository-relative paths.
# Statistical logic, parameters, thresholds, and seeds are unchanged from the archived final script.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1]))); source(file.path(script_dir, "common.R"))
x <- read_corrected(); m <- x$metadata; m$cell_id_revision <- rownames(m)
keep <- c("cell_id_revision","cohort","condition","donor_id","sample_id","GART_count_corrected","GART_expr_corrected","GART_detected_corrected","GART_status","nCount_RNA","nFeature_RNA","S_score","G2M_score","ISC_subcluster")
write.csv(m[, keep], file.path(output_root,"INTERMEDIATE","ISC_corrected_cell_metadata.csv"), row.names=FALSE)
sample_audit <- aggregate(rep(1,nrow(m)), list(cohort=m$cohort,condition=as.character(m$condition),donor_id=m$donor_id,sample_id=m$sample_id), sum); names(sample_audit)[5] <- "n_ISClike_cells"
sample_audit$raw_GART_missing <- ave(is.na(m$GART_count_corrected), m$sample_id, FUN=sum)[match(sample_audit$sample_id,m$sample_id)]
sample_audit$included <- sample_audit$raw_GART_missing == 0
write.csv(sample_audit,file.path(output_root,"AUDIT","singlecell_sample_audit.csv"),row.names=FALSE)
donor_map <- aggregate(sample_audit$n_ISClike_cells,list(cohort=sample_audit$cohort,condition=sample_audit$condition,donor_id=sample_audit$donor_id),sum);names(donor_map)[4]<-"n_ISClike_cells"
donor_map$n_samples <- as.integer(table(sample_audit$donor_id)[donor_map$donor_id]); donor_map$sample_ids <- vapply(donor_map$donor_id,function(z)paste(sample_audit$sample_id[sample_audit$donor_id==z],collapse="|"),"")
write.csv(donor_map,file.path(output_root,"AUDIT","singlecell_donor_mapping_audit.csv"),row.names=FALSE)
file.copy(file.path(singlecell_root, "RAW_LAYER_RECOVERY_AUDIT.csv"),file.path(output_root,"AUDIT","singlecell_raw_layer_recovery_audit.csv"),overwrite=TRUE)
qc <- aggregate(cbind(nCount_RNA=m$nCount_RNA,nFeature_RNA=m$nFeature_RNA)~cohort+condition,data=m,FUN=function(v)c(n=length(v),median=median(v),q1=quantile(v,.25),q3=quantile(v,.75)))
write.csv(qc,file.path(output_root,"AUDIT","singlecell_QC_summary.csv"),row.names=FALSE)
writeLines(c("Single-cell inputs rebuilt without modifying source objects.",sprintf("ISC-like cells: %d; donors: %d; samples: %d; cohorts: %d.",nrow(m),length(unique(m$donor_id)),length(unique(m$sample_id)),length(unique(m$cohort))),sprintf("Missing verified raw GART counts: %d.",sum(is.na(m$GART_count_corrected)))),file.path(output_root,"REPORTS","01_input_rebuild_log.txt"))
