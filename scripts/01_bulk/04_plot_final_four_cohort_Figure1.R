# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 20_BULK_EXTERNAL_COHORT_SCREEN/scripts/plot_fig1a_four_cohorts.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
work_root <- file.path(paths$work_root, "bulk")
outdir <- file.path(paths$results_root, "bulk", "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_old <- function(ds){
  p <- file.path(work_root,"intermediate",paste0("archived_",ds,if(ds=="GSE37283")"_retained_data.csv" else "_sample_data.csv"))
  d <- read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
  d$condition <- ifelse(d$condition=="UC","UC","Healthy")
  d[,c("dataset","GSM_ID","condition","GART_expression")]
}
dlist <- list(GSE24287=read_old("GSE24287"),GSE36807=read_old("GSE36807"),GSE37283=read_old("GSE37283"))
gnew <- read.csv(file.path(outdir,"GSE87466_sample_level.csv"),check.names=FALSE,stringsAsFactors=FALSE)
gnew <- gnew[gnew$retained_or_excluded=="retained",]
gnew$condition <- ifelse(gnew$condition=="Active UC","UC","Healthy")
dlist$GSE87466 <- gnew[,c("dataset","GSM_ID","condition","GART_expression")]

oldstats_file <- gart_require_file(file.path(paths$external_data_root, "bulk", "legacy_three_cohort", "bulk_meta_results_revised.csv"), "archived three-cohort meta-analysis result")
oldstats <- read.csv(oldstats_file, check.names = FALSE)
oldstats <- oldstats[oldstats$row_type=="cohort",]
newstats <- read.csv(file.path(outdir,"GSE87466_GART_result.csv"),check.names=FALSE)
pvals <- c(setNames(oldstats$Welch_t_P,oldstats$dataset),GSE87466=newstats$P_value)
effects <- read.csv(file.path(outdir,"four_cohort_effect_sizes.csv"),check.names=FALSE)
meta <- read.csv(file.path(outdir,"four_cohort_meta_summary.csv"),check.names=FALSE)
loo <- read.csv(file.path(outdir,"four_cohort_leave_one_out.csv"),check.names=FALSE)

burgundy <- "#8C001A"; pink <- "#F47A96"; dark <- "#202020"
pformat <- function(p) if(p<0.0001) format(p,digits=2,scientific=TRUE) else sprintf("%.4f",p)
set_font <- function(){if(.Platform$OS.type=="windows") windowsFonts(Arial=windowsFont("Arial")); par(family="Arial")}
open_tiff <- function(path,w,h){tiff(path,width=w,height=h,units="in",res=300,compression="lzw",bg="white",type="windows");set_font()}
open_pdf <- function(path,w,h){cairo_pdf(path,width=w,height=h,family="Arial",bg="white");par(family="Arial")}

cohort_panel <- function(ds, ylab=TRUE){
  d <- dlist[[ds]]; h <- d$GART_expression[d$condition=="Healthy"]; u <- d$GART_expression[d$condition=="UC"]
  vv <- c(h,u); span <- diff(range(c(0,vv))); if(span==0)span <- 1
  yl <- c(min(c(0,vv))-.06*span,max(c(0,vv))+.28*span)
  par(mar=c(3.2,if(ylab)3.65 else 2.45,2.4,.65),mgp=c(2.15,.55,0),tcl=-.22)
  plot(NA,xlim=c(.45,2.55),ylim=yl,xaxt="n",xlab="",ylab=if(ylab)"GART expression level" else "",bty="l",las=1,cex.axis=.78,cex.lab=.82)
  axis(1,at=1:2,labels=c("Health","UC"),tick=FALSE,cex.axis=.82,line=-.25)
  means <- c(mean(h),mean(u)); sds <- c(sd(h),sd(u))
  rect(c(.67,1.67),0,c(1.33,2.33),means,col=c(burgundy,pink),border=NA)
  arrows(1:2,means-sds,1:2,means+sds,angle=90,code=3,length=.06,lwd=1.1,col=dark)
  jh <- if(length(h)>1) seq(-.12,.12,length.out=length(h)) else 0
  ju <- if(length(u)>1) seq(-.12,.12,length.out=length(u)) else 0
  points(1+jh,h,pch=16,col=burgundy,cex=.72)
  points(2+ju,u,pch=15,col=pink,cex=.72)
  mtext(ds,side=3,line=.35,font=2,cex=.9)
  text(1.5,yl[2]-.08*diff(yl),paste0("P = ",pformat(pvals[ds])),cex=.78)
  text(1.5,yl[2]-.17*diff(yl),sprintf("n = %d / %d",length(h),length(u)),cex=.74)
}

draw_four <- function(){layout(matrix(1:4,nrow=1),widths=rep(1,4));for(i in seq_along(dlist))cohort_panel(names(dlist)[i],ylab=i==1)}
draw_single <- function(){cohort_panel("GSE87466",ylab=TRUE)}

open_tiff(file.path(outdir,"Fig1A_GSE87466_GART_expression.tiff"),2.15,3.35);draw_single();dev.off()
open_pdf(file.path(outdir,"Fig1A_GSE87466_GART_expression.pdf"),2.15,3.35);draw_single();dev.off()
open_tiff(file.path(outdir,"Fig1A_four_bulk_cohorts.tiff"),7.09,3.35);draw_four();dev.off()
open_pdf(file.path(outdir,"Fig1A_four_bulk_cohorts.pdf"),7.09,3.35);draw_four();dev.off()

forest_plot <- function(){
  par(mar=c(5.1,5.3,2.2,.65),mgp=c(2.55,.65,0),tcl=-.22)
  y <- c(5.2,4.2,3.2,2.2,1.05); est <- c(effects$Hedges_g,meta$pooled_Hedges_g)
  lo <- c(effects$CI_low,meta$CI_low); hi <- c(effects$CI_high,meta$CI_high)
  xr <- c(-.75,4.95)
  plot(NA,xlim=xr,ylim=c(.25,5.85),yaxt="n",xaxt="n",xlab="Hedges' g (UC - Health)",ylab="",bty="l",las=1,cex.axis=.82,cex.lab=.88)
  axis(1,at=0:3,cex.axis=.82)
  axis(2,at=y,labels=c(effects$dataset,"Pooled"),las=1,tick=FALSE,font=c(1,1,1,1,2),cex.axis=.82)
  abline(v=0,lty=2,col="#666666")
  segments(lo[1:4],y[1:4],hi[1:4],y[1:4],lwd=1.35,col=dark)
  points(est[1:4],y[1:4],pch=19,cex=.9,col=dark)
  polygon(c(lo[5],est[5],hi[5],est[5]),c(y[5],y[5]+.20,y[5],y[5]-.20),col=burgundy,border=burgundy)
  text(3.20,5.65,"g [95% CI] (weight)",pos=4,font=2,cex=.66)
  txt <- c(sprintf("%.2f [%.2f, %.2f] (%.1f%%)",effects$Hedges_g,effects$CI_low,effects$CI_high,effects$cohort_weight_percent),sprintf("%.2f [%.2f, %.2f]",meta$pooled_Hedges_g,meta$CI_low,meta$CI_high))
  text(3.20,y,txt,pos=4,cex=.63)
  title("Random-effects meta-analysis of GART expression",font.main=2,cex.main=.98)
  mtext(sprintf("REML + Hartung-Knapp: pooled g = %.3f (95%% CI %.3f to %.3f), P = %.4f",meta$pooled_Hedges_g,meta$CI_low,meta$CI_high,meta$P_value),side=1,line=3.05,cex=.67)
  mtext(sprintf("Q = %.3f (P = %.4f), I2 = %.1f%%, tau2 = %.4f",meta$Cochran_Q,meta$Q_P,meta$I2_percent,meta$tau2_REML),side=1,line=3.80,cex=.67)
}
open_tiff(file.path(outdir,"Supplementary_FigS1A_bulk_GART_meta_analysis.tiff"),7.09,4.45);forest_plot();dev.off()
open_pdf(file.path(outdir,"Supplementary_FigS1A_bulk_GART_meta_analysis.pdf"),7.09,4.45);forest_plot();dev.off()

# Exact plotting source data and statistics retained for audit.
write.csv(do.call(rbind,dlist),file.path(outdir,"Fig1A_four_bulk_cohorts_source_data.csv"),row.names=FALSE)
write.csv(loo,file.path(outdir,"four_cohort_leave_one_out.csv"),row.names=FALSE)
