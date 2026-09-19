#!/usr/bin/env Rscript
a <- commandArgs(FALSE)
f <- grep("^--file=", a, value=TRUE)
sd <- dirname(normalizePath(sub("^--file=", "", f[1])))
repo <- normalizePath(file.path(sd, "..", ".."), mustWork=TRUE)
inp <- file.path(repo, "data_public", "figure_source_data", "Figure7", "Figure7_formal_statistical_inputs.csv")
refp <- file.path(repo, "data_public", "figure_source_data", "Figure7", "Fig7_statistical_summary_revised.csv")
d <- read.csv(inp, check.names=FALSE, stringsAsFactors=FALSE)
ref <- read.csv(refp, check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(nrow(d)==752, length(unique(paste(d$panel,d$endpoint,sep="||")))==13)
stopifnot(all(c("panel","endpoint","anonymous_observation_id","group","exact_formal_value","unit","source_analysis_object") %in% names(d)))
stopifnot(!anyDuplicated(d[c("panel","endpoint","group","anonymous_observation_id")]))
keys <- unique(paste(d$panel,d$endpoint,sep="||"))
rows <- lapply(keys, function(k) {
  z <- d[paste(d$panel,d$endpoint,sep="||")==k,]
  x <- z$exact_formal_value[z$group=="HC"]
  y <- z$exact_formal_value[z$group=="UC"]
  is_welch <- grepl("^CD74\\+ objects per mm", z$endpoint[1]) || z$endpoint[1] %in% c("Ki67+ cells (%)", "GART+Ki67+ cells (%)")
  if (is_welch) {
    tt <- t.test(y,x,var.equal=FALSE,alternative="two.sided")
    data.frame(panel=z$panel[1], endpoint=z$endpoint[1], HC_n=length(x), UC_n=length(y),
      HC_summary=sprintf("mean %.4f; SD %.4f",mean(x),sd(x)), UC_summary=sprintf("mean %.4f; SD %.4f",mean(y),sd(y)),
      effect_estimate=mean(y)-mean(x), CI_lower=unname(tt$conf.int[1]), CI_upper=unname(tt$conf.int[2]),
      statistic=unname(tt$statistic), df=unname(tt$parameter), P=tt$p.value, P_display=ifelse(tt$p.value<1e-4,"<0.0001",sprintf("%.4f",tt$p.value)), stringsAsFactors=FALSE)
  } else {
    has_ties <- any(duplicated(c(x,y)))
    wt <- wilcox.test(x,y,alternative="two.sided",exact=!has_ties,correct=FALSE)
    u_hc <- sum(rank(c(x,y))[seq_along(x)]) - length(x)*(length(x)+1)/2
    diffs <- as.vector(outer(y,x,"-"))
    hl <- median(diffs)
    ci <- c(NA_real_,NA_real_)
    if (!has_ties) {
      wc <- suppressWarnings(wilcox.test(y,x,alternative="two.sided",exact=TRUE,correct=FALSE,conf.int=TRUE,conf.level=.95))
      ci <- unname(wc$conf.int)
    }
    data.frame(panel=z$panel[1], endpoint=z$endpoint[1], HC_n=length(x), UC_n=length(y),
      HC_summary=sprintf("median %.4f; IQR %.4f\u2013%.4f",median(x),quantile(x,.25),quantile(x,.75)),
      UC_summary=sprintf("median %.4f; IQR %.4f\u2013%.4f",median(y),quantile(y,.25),quantile(y,.75)),
      effect_estimate=hl, CI_lower=ci[1], CI_upper=ci[2], statistic=u_hc, df=NA_real_,
      P=if(has_ties) NA_real_ else wt$p.value, P_display=ifelse(has_ties,"<0.0001",ifelse(wt$p.value<1e-4,"<0.0001",sprintf("%.4f",wt$p.value))), stringsAsFactors=FALSE)
  }
})
out <- do.call(rbind, rows)
ref2 <- data.frame(panel=ref$Panel,endpoint=ref$Endpoint,HC_n=ref$HC_n,UC_n=ref$UC_n,effect_estimate=ref$Effect_estimate,CI_lower=ref$CI_lower,CI_upper=ref$CI_upper,statistic=ref$Statistic,df=ref$df,P=ref$P,P_display=ref$P_display,stringsAsFactors=FALSE)
m <- merge(out, ref2, by=c("panel","endpoint"), suffixes=c("_reproduced","_reference"), sort=FALSE)
eq <- function(a,b,tol=1e-10) (is.na(a)&is.na(b)) | (!is.na(a)&!is.na(b)&abs(a-b)<tol)
m$n_PASS <- m$HC_n_reproduced==m$HC_n_reference & m$UC_n_reproduced==m$UC_n_reference
m$statistic_PASS <- eq(m$statistic_reproduced,m$statistic_reference)
m$effect_PASS <- eq(m$effect_estimate_reproduced,m$effect_estimate_reference)
m$CI_PASS <- eq(m$CI_lower_reproduced,m$CI_lower_reference) & eq(m$CI_upper_reproduced,m$CI_upper_reference)
m$P_PASS <- ifelse(is.na(m$P_reference),m$P_display_reproduced==m$P_display_reference,eq(m$P_reproduced,m$P_reference,1e-12))
m$exact_PASS <- m$n_PASS & m$statistic_PASS & m$effect_PASS & m$CI_PASS & m$P_PASS
od <- file.path(repo,"results_reproduced","imaging")
dir.create(od,recursive=TRUE,showWarnings=FALSE)
write.csv(out,file.path(od,"Figure7_statistics_reproduced_full.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(m,file.path(od,"Figure7_reproduction_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
stopifnot(nrow(m)==13, all(m$exact_PASS))
cat("Figure 7 exact reproduction: 13/13 PASS\n")
