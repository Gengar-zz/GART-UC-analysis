#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
ff <- grep("^--file=", args, value=TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", ff[1])))
repo <- normalizePath(file.path(script_dir, "..", ".."), mustWork=TRUE)
inp <- file.path(repo, "data_public", "figure_source_data", "Figure7", "Figure7_formal_statistical_inputs.csv")
refp <- file.path(repo, "data_public", "figure_source_data", "Figure7", "Fig7_statistical_summary_revised.csv")
d <- read.csv(inp, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
ref <- read.csv(refp, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
stopifnot(nrow(d)==752L, length(unique(paste(d$panel,d$endpoint,sep="||")))==13L)
stopifnot(all(c("panel","endpoint","anonymous_observation_id","group","value","unit") %in% names(d)))
stopifnot(!anyDuplicated(d[c("panel","endpoint","group","anonymous_observation_id")]))

fmt_p <- function(p) ifelse(is.na(p) | p < 1e-4, "<0.0001", sprintf("%.4f",p))
keys <- unique(paste(d$panel,d$endpoint,sep="||"))
rows <- lapply(keys, function(k) {
  z <- d[paste(d$panel,d$endpoint,sep="||")==k,]
  x <- z$value[z$group=="HC"]; y <- z$value[z$group=="UC"]
  is_welch <- grepl("^CD74\\+ objects per mm", z$endpoint[1]) || z$endpoint[1] %in% c("Ki67+ cells (%)","GART+Ki67+ cells (%)")
  if (is_welch) {
    tt <- t.test(y,x,var.equal=FALSE,alternative="two.sided")
    data.frame(panel=z$panel[1],endpoint=z$endpoint[1],HC_n=length(x),UC_n=length(y),
      effect_estimate=mean(y)-mean(x),CI_lower=tt$conf.int[1],CI_upper=tt$conf.int[2],
      statistic=unname(tt$statistic),df=unname(tt$parameter),P=tt$p.value,P_display=fmt_p(tt$p.value))
  } else {
    ties <- anyDuplicated(c(x,y))>0
    wt <- suppressWarnings(wilcox.test(x,y,exact=!ties,alternative="two.sided"))
    h <- if (!ties) suppressWarnings(wilcox.test(y,x,exact=TRUE,conf.int=TRUE,alternative="two.sided")) else NULL
    dif <- as.vector(outer(y,x,"-"))
    data.frame(panel=z$panel[1],endpoint=z$endpoint[1],HC_n=length(x),UC_n=length(y),
      effect_estimate=median(dif),CI_lower=if(is.null(h)) NA_real_ else h$conf.int[1],
      CI_upper=if(is.null(h)) NA_real_ else h$conf.int[2],statistic=unname(wt$statistic),df=NA_real_,
      P=if(ties) NA_real_ else wt$p.value,P_display=if(ties) "<0.0001" else fmt_p(wt$p.value))
  }
})
out <- do.call(rbind,rows)
r2 <- data.frame(panel=ref$Panel,endpoint=ref$Endpoint,HC_n=ref$HC_n,UC_n=ref$UC_n,
  effect_estimate=ref$Effect_estimate,CI_lower=ref$CI_lower,CI_upper=ref$CI_upper,
  statistic=ref$Statistic,df=ref$df,P=ref$P,P_display=ref$P_display,check.names=FALSE)
m <- merge(out,r2,by=c("panel","endpoint"),suffixes=c("_reproduced","_reference"),sort=FALSE)
eq <- function(a,b,d=2) (is.na(a)&is.na(b)) | (!is.na(a)&!is.na(b)&round(a,d)==round(b,d))
m$PASS <- m$HC_n_reproduced==m$HC_n_reference & m$UC_n_reproduced==m$UC_n_reference &
  eq(m$effect_estimate_reproduced,m$effect_estimate_reference,2) & eq(m$CI_lower_reproduced,m$CI_lower_reference,2) &
  eq(m$CI_upper_reproduced,m$CI_upper_reference,2) & eq(m$statistic_reproduced,m$statistic_reference,3) &
  eq(m$df_reproduced,m$df_reference,3) & m$P_display_reproduced==m$P_display_reference
od <- file.path(repo,"results_reproduced","imaging"); dir.create(od,recursive=TRUE,showWarnings=FALSE)
write.csv(out,file.path(od,"Figure7_statistics_reproduced_full.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(m,file.path(od,"Figure7_reproduction_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
stopifnot(nrow(m)==13L,all(m$PASS))
cat("Figure 7 current-Prism reproduction: 13/13 PASS\n")
