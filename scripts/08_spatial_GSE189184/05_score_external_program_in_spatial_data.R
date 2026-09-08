# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 19_SPATIAL_REVISED_PRIMARY_ADJUDICATION/code/06_unblind_score_spatial.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(Matrix))

root <- file.path(paths$work_root, "spatial_GSE189184")
obj_path <- gart_require_file(file.path(paths$external_data_root, "spatial", "GSE189184", "frozen_objects", "GSE189184_ST_with_FINAL_Fig7_scores_STABLE.rds"), "frozen GSE189184 spatial object")
spot_path <- gart_require_file(file.path(paths$external_data_root, "spatial", "GSE189184", "spatial_optimization_spot_source.csv"), "upstream spatial spot metadata")

stopifnot(file.exists(file.path(root,"reports","REVISED_PRIMARY_ADJUDICATION_LOCK.md")))
stopifnot(file.exists(file.path(root,"reports","SINGLECELL_SELECTED_MODEL_LOCK.md")))
stopifnot(file.exists(file.path(root,"reports","SPATIAL_ROI_RESOLUTION_LOCK.md")))

x <- readRDS(obj_path)
assay <- attr(x,"assays")[["Spatial"]]
layers <- attr(assay,"layers")
feature_map <- attr(assay,"features"); attr(feature_map,"class") <- NULL
cell_map <- attr(assay,"cells"); attr(cell_map,"class") <- NULL
lnames <- grep("^data\\.",names(layers),value=TRUE)
mats <- lapply(lnames,function(ln){
  z<-layers[[ln]]; rownames(z)<-rownames(feature_map)[feature_map[,ln]]; colnames(z)<-rownames(cell_map)[cell_map[,ln]]; z
})
common_genes <- Reduce(intersect,lapply(mats,rownames))
mat <- do.call(cbind,lapply(mats,function(z) z[common_genes,,drop=FALSE]))
spot <- read.csv(spot_path,check.names=FALSE)
common_spots <- intersect(colnames(mat),spot$spot_id)
stopifnot(length(common_spots)==8474L)
mat <- mat[,common_spots,drop=FALSE]
spot <- spot[match(common_spots,spot$spot_id),]

members <- read.csv(file.path(root,"source_data","candidate_signature_members.csv"),check.names=FALSE)
weights <- read.csv(file.path(root,"source_data","selected_signature_gene_list.csv"),check.names=FALSE)
sigs <- split(members$gene,members$signature)
selected_name <- "SC_GART_NEAREST_CENTROID_PROGRAM"
selected_genes <- intersect(weights$gene,rownames(mat))
selected_weights <- weights$projection_weight[match(selected_genes,weights$gene)]
validation <- c("MYBL2","PCNA","PAICS","MCM7","GINS2","OLFM4","MIF","CD74","GART")
union_genes <- unique(c(unlist(sigs),validation))
union_genes <- intersect(union_genes,rownames(mat))
idx_union <- match(union_genes,rownames(mat)); names(idx_union)<-union_genes
N <- nrow(mat)
rank_pct <- matrix(0.5,nrow=length(union_genes),ncol=ncol(mat),dimnames=list(union_genes,colnames(mat)))
for(j in seq_len(ncol(mat))){
  v <- mat[,j]
  nz <- which(v!=0)
  rr <- rank(-as.numeric(v[nz]),ties.method="average")
  pos <- match(nz,idx_union,nomatch=0L)
  keep <- pos>0
  if(any(keep)) rank_pct[pos[keep],j] <- 1-rr[keep]/N
}

score_one <- function(gs,method){
  gs<-intersect(gs,rownames(rank_pct)); r<-rank_pct[gs,,drop=FALSE]; nr<-N*(1-r); m<-length(gs)
  if(method=="mean_rank" || method=="singscore") return(colMeans(r))
  if(method=="UCell"){
    maxRank<-1500; cap<-pmin(nr,maxRank+1); u<-colSums(cap)-m*(m+1)/2; return(1-u/(m*maxRank))
  }
  if(method=="AUCell"){
    cutoff<-round(.05*N); contrib<-pmax(cutoff-nr+1,0); denom<-m*cutoff-m*(m-1)/2; return(colSums(contrib)/max(denom,1))
  }
  if(method=="ssGSEA"){
    pos<-N-nr+1; w<-pmax(pos,1)^.25; hit<-colSums(w*pos)/pmax(colSums(w),1e-12); miss<-(N*(N+1)/2-colSums(pos))/(N-m); return((hit-miss)/N)
  }
  if(method=="mean_z"){
    z<-t(scale(t(as.matrix(mat[gs,,drop=FALSE])))); z[!is.finite(z)]<-0; return(colMeans(z))
  }
}

out <- data.frame(
  spot_id=spot$spot_id,sample_id=spot$sample_id,donor_id=spot$donor_id,condition=spot$condition,
  array_row=spot$array_row,array_col=spot$array_col,pixel_row=spot$pixel_row,pixel_col=spot$pixel_col,
  GART_expression=spot$GART_expression,GART_count=spot$GART_count,
  epithelial_score=spot$ROI_A_score,stringsAsFactors=FALSE
)
out$epithelial_percentile <- NA_real_
for(sec in unique(out$sample_id)){
  ii<-which(out$sample_id==sec); out$epithelial_percentile[ii]<-rank(out$epithelial_score[ii],ties.method="average")/length(ii)
}
out$ROI_top10 <- out$epithelial_percentile>=.90
out$ROI_top20 <- out$epithelial_percentile>=.80
out$ROI_top30 <- out$epithelial_percentile>=.70
out$ROI_top40 <- out$epithelial_percentile>=.60
out$ROI_all_tissue <- TRUE
out$ROI_one_ring_top20 <- spot$ROI_neighborhood_top20
out$ROI_outcome_excluded_domain <- spot$domain_fallback

for(g in validation) out[[g]] <- if(g %in% rownames(mat)) as.numeric(mat[g,]) else NA_real_
for(sn in names(sigs)){
  for(method in c("UCell","singscore","AUCell","ssGSEA","mean_rank","mean_z")){
    out[[paste(sn,method,sep="__")]] <- score_one(sigs[[sn]],method)
  }
}
sw <- selected_weights/sum(abs(selected_weights))
out[[paste(selected_name,"nearest_centroid_similarity",sep="__")]] <- as.numeric(crossprod(sw,rank_pct[selected_genes,,drop=FALSE]))

write.csv(out,file.path(root,"source_data","all_spatial_program_scores.csv"),row.names=FALSE,quote=TRUE)
cat("unblinded scored spots",nrow(out),"selected matched genes",length(selected_genes),"of",nrow(weights),"\n")
