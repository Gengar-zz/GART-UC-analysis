# Public-release copy: file paths were converted to repository-relative paths.
# Statistical logic, parameters, thresholds, and seeds are unchanged from the archived final script.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1]))); source(file.path(script_dir, "common.R"))
suppressPackageStartupMessages({library(SeuratObject);library(Matrix)})
ep <- readRDS(file.path(singlecell_root,"UC_GART_06c_Epithelial_Fig2AB_noGARTname_FINAL.rds")); em <- ep[[]]; em$cell_id <- rownames(em)
coords <- Embeddings(ep,"umap_epi"); subtype <- as.character(em$Epi_Subtype_Fig2); subtype[is.na(subtype)|!nzchar(subtype)] <- as.character(em$Epi_Subtype_Preliminary[is.na(subtype)|!nzchar(subtype)])
set.seed(20260904); idx <- unlist(lapply(split(seq_len(nrow(em)),subtype),function(ii)sample(ii,min(length(ii),2500))))
epi_plot <- data.frame(cell_id=rownames(em)[idx],UMAP1=coords[idx,1],UMAP2=coords[idx,2],subtype=subtype[idx],cohort=em$Dataset_Final[idx],condition=em$Disease_Final[idx])
write.csv(epi_plot,file.path(output_root,"INTERMEDIATE","Fig2A_epithelial_UMAP_plot_data.csv"),row.names=FALSE)
markers <- list(`ISC-like`=c("OLFM4","LGR5","SMOC2"),`TA/proliferating`=c("MKI67","TOP2A","PCNA"),Goblet=c("MUC2","SPINK4","CLCA1"),Absorptive=c("CA1","KRT20","FABP1"),BEST4=c("BEST4","OTOP2","CA7"),Enteroendocrine=c("CHGA","CHGB","NEUROD1"),Tuft=c("POU2F3","TRPM5","AVIL"))
genes <- unique(unlist(markers)); a <- ep[["RNA"]]
z <- matrix(0, nrow=length(genes), ncol=nrow(em), dimnames=list(genes,rownames(em)))
layer_full <- function(ln){q<-methods::slot(a,"layers")[[ln]];fn<-methods::slot(a,"features")[[ln]];cn<-methods::slot(a,"cells")[[ln]];if(is.null(dim(q)))q<-matrix(q,nrow=length(fn),ncol=length(cn));if(is.null(rownames(q)))rownames(q)<-fn;if(is.null(colnames(q)))colnames(q)<-cn;q}
for(ln in grep("^data",Layers(a),value=TRUE)){q<-layer_full(ln);g<-intersect(genes,rownames(q));z[g,colnames(q)]<-as.matrix(q[g,,drop=FALSE])}
common <- colnames(z); st <- subtype[match(common,rownames(em))]
dot <- do.call(rbind,lapply(rownames(z),function(g)do.call(rbind,lapply(split(seq_along(st),st),function(ii)data.frame(gene=g,subtype=st[ii[1]],mean_expression=mean(z[g,ii]),pct_detected=mean(z[g,ii]>0)*100)))))
write.csv(dot,file.path(output_root,"RESULTS","epithelial_marker_dot_summary.csv"),row.names=FALSE)
gart <- aggregate(cbind(mean_expression=em$GART_expr,detected=as.numeric(em$GART_expr>0)),list(subtype=subtype,cohort=em$Dataset_Final),mean,na.rm=TRUE);gart$pct_detected<-100*gart$detected
write.csv(gart,file.path(output_root,"RESULTS","epithelial_GART_subtype_summary.csv"),row.names=FALSE)
rm(ep,z);gc()

isc <- readRDS(file.path(singlecell_root,"UC_GART_07_ISC_like_subclustered_FIXED_withDonor.rds")); cm <- read_corrected()$metadata; ic <- Embeddings(isc,"iscUMAP"); im <- cm[rownames(ic),,drop=FALSE]
isc_plot <- data.frame(cell_id=rownames(ic),UMAP1=ic[,1],UMAP2=ic[,2],cohort=im$cohort,condition=as.character(im$condition),donor_id=im$donor_id,sample_id=im$sample_id,GART_status=as.character(im$GART_status),GART_expr=im$GART_expr_corrected,ISC_subcluster=im$ISC_subcluster,purine_score=im$Reactome_purine_GART_excluded,TA_score=im$TA_proliferating,ISC_score=im$GO_intestinal_stem_cell_homeostasis)
write.csv(isc_plot,file.path(output_root,"INTERMEDIATE","ISC_UMAP_plot_data.csv"),row.names=FALSE)
write.csv(data.frame(metric=c("epithelial_cells","ISC_like_cells","ISC_embedding_cells","ISC_embedding_missing"),value=c(nrow(em),nrow(cm),nrow(isc_plot),nrow(cm)-nrow(isc_plot))),file.path(output_root,"AUDIT","annotation_embedding_audit.csv"),row.names=FALSE)
