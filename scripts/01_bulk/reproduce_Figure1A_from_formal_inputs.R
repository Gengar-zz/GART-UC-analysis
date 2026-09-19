#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
repo <- normalizePath(file.path(script_dir, "..", ".."), mustWork=TRUE)
input <- if (length(args) >= 1) args[1] else file.path(repo, "data_public", "bulk_frozen", "bulk_GART_formal_input_values.tsv")
output <- if (length(args) >= 2) args[2] else file.path(repo, "results_reproduced", "bulk", "Figure1A_exact_reproduction.tsv")
d <- read.delim(input, check.names=FALSE, stringsAsFactors=FALSE)
expected_p <- c(GSE24287=0.0186431019804659, GSE36807=0.00530481877805238, GSE37283=0.041140885607957, GSE87466=0.0139570795396182)
expected_g <- c(GSE24287=0.656483542432513, GSE36807=1.28679889263754, GSE37283=1.48416460764612, GSE87466=0.595855304374862)
calc <- function(z, cohort) {
  h <- z$GART_expression[z$group == "Healthy"]; u <- z$GART_expression[z$group == "UC"]
  n1 <- length(h); n2 <- length(u); df <- n1+n2-2
  sp <- sqrt(((n1-1)*var(h)+(n2-1)*var(u))/df)
  dd <- (mean(u)-mean(h))/sp
  J <- gamma(df/2)/(sqrt(df/2)*gamma((df-1)/2))
  g <- J*dd
  # Preserve the locked historical sampling-variance implementations used by
  # the legacy three-cohort export and the independently added GSE87466 route.
  effect_for_v <- if (cohort == "GSE87466") dd else g
  v <- (n1+n2)/(n1*n2)+effect_for_v^2/(2*(n1+n2))
  wt <- t.test(u,h,var.equal=FALSE)
  data.frame(n_Healthy=n1,n_UC=n2,mean_Healthy=mean(h),sd_Healthy=sd(h),mean_UC=mean(u),sd_UC=sd(u),Welch_P=wt$p.value,Hedges_g=g,variance_g=v,CI_lower=g-1.96*sqrt(v),CI_upper=g+1.96*sqrt(v))
}
cohorts <- unique(d$cohort)
res <- do.call(rbind,lapply(cohorts,function(x)cbind(cohort=x,calc(d[d$cohort==x,],x))))
res$P_exact_PASS <- abs(res$Welch_P-expected_p[res$cohort]) < 1e-12
res$g_exact_PASS <- abs(res$Hedges_g-expected_g[res$cohort]) < 1e-12
yi <- res$Hedges_g; vi <- res$variance_g; k <- length(yi)
w0<-1/vi; mu0<-sum(w0*yi)/sum(w0); Q<-sum(w0*(yi-mu0)^2); Qp<-pchisq(Q,k-1,lower.tail=FALSE); I2<-max(0,(Q-(k-1))/Q)*100
nll <- function(tau2){w<-1/(vi+tau2); mu<-sum(w*yi)/sum(w); .5*(sum(log(vi+tau2))+log(sum(w))+sum(w*(yi-mu)^2))}
tau2 <- if (Q <= k-1) 0 else optimize(nll,c(0,10),tol=.Machine$double.eps^.25)$minimum
w <- 1/(vi+tau2); pooled <- sum(w*yi)/sum(w); qhk <- sum(w*(yi-pooled)^2)/(k-1); sehk <- sqrt(qhk/sum(w))
ci <- pooled + c(-1,1)*qt(.975,k-1)*sehk; php <- 2*pt(-abs(pooled/sehk),k-1)
dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE); write.table(res,output,sep="\t",row.names=FALSE,quote=FALSE)
meta <- data.frame(method="REML with Hartung-Knapp",pooled_g=pooled,CI_lower=ci[1],CI_upper=ci[2],P=php,Q=Q,Q_P=Qp,I2_percent=I2,tau2=tau2)
write.table(meta,sub("\\.tsv$","_meta.tsv",output),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(all(res$P_exact_PASS), all(res$g_exact_PASS))
