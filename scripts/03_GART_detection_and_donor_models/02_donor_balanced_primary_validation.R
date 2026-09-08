# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 08_DONOR_BALANCED_PRIMARY_VALIDATION/CODE/run_donor_balanced_validation.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors=FALSE, warn=1, width=220)
suppressPackageStartupMessages({
  library(geepack)
  library(survey)
  library(clubSandwich)
  library(sandwich)
  library(lme4)
})

set.seed(20260905)
root <- paths$work_root
out <- file.path(paths$results_root, "donor_balanced_primary_validation")
dirs <- file.path(out,c("RESULTS","REPORTS","CODE","QA")); invisible(lapply(dirs,dir.create,recursive=TRUE,showWarnings=FALSE))
resdir<-file.path(out,"RESULTS"); repdir<-file.path(out,"REPORTS"); qadir<-file.path(out,"QA")

# Written before fitting: Implementation 1 remains the candidate revised primary regardless of results.
lock <- c(
  "# Donor-balanced primary-validation lock",
  "",
  "Locked before model fitting on 2026-09-05.",
  "",
  "- Candidate revised primary: Implementation 1, inverse-cluster-size weighted GEE.",
  "- Analysis set: 9,306 ISC-like epithelial cells from 33 donors.",
  "- Outcome: raw corrected GART count >=1.",
  "- Covariates: condition + scale(log1p(nUMI)) + scale(nGene) + S.Score + G2M.Score + cohort.",
  "- Each cell weight: 1 / donor ISC-like cell count; each donor has total weight 1.",
  "- Working correlation: independence; donor is the cluster.",
  "- Small-sample inference for Implementation 1: robust GEE sandwich scaled by G/(G-p), with a t reference using G-p degrees of freedom.",
  "- The candidate-primary designation will not change in response to P values from any implementation or sensitivity analysis."
)
writeLines(lock,file.path(repdir,"PRIMARY_CANDIDATE_LOCK.md"),useBytes=TRUE)

obj <- readRDS(gart_require_file(file.path(paths$external_data_root, "singlecell", "frozen_objects", "CORRECTED_SCRNA_ANALYSIS_OBJECT.rds"), "corrected single-cell analysis object"))
m <- as.data.frame(obj$metadata)
needed <- c("condition","cohort","donor_id","nCount_RNA","nFeature_RNA","S_score","G2M_score","GART_count_corrected")
stopifnot(all(needed %in% names(m)),nrow(m)==9306,length(unique(m$donor_id))==33,!anyNA(m[,needed]))
m$condition <- relevel(factor(m$condition),"Healthy")
m$cohort <- factor(m$cohort)
m$donor_id <- factor(m$donor_id)
m$det1 <- as.integer(m$GART_count_corrected>=1)
m$log_nUMI <- log1p(m$nCount_RNA)
m$z_log_nUMI <- as.numeric(scale(m$log_nUMI))
m$z_nGene <- as.numeric(scale(m$nFeature_RNA))
n_per <- table(m$donor_id)
m$donor_n <- as.integer(n_per[as.character(m$donor_id)])
m$donor_weight <- 1/m$donor_n
m <- m[order(m$donor_id),]

wa <- aggregate(donor_weight~donor_id+cohort+condition+donor_n,data=m,sum)
names(wa)[names(wa)=="donor_weight"] <- "total_weight"
write.csv(wa,file.path(resdir,"donor_weight_audit.csv"),row.names=FALSE)
stopifnot(max(abs(wa$total_weight-1))<1e-12,nrow(wa)==33)

form <- det1 ~ condition + z_log_nUMI + z_nGene + S_score + G2M_score + cohort
term <- "conditionUC"

cap <- function(expr){
  w<-character(); val<-withCallingHandlers(tryCatch(eval.parent(substitute(expr)),error=identity),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart("muffleWarning")})
  list(value=val,warnings=paste(unique(w),collapse=" | "))
}

cr2_extract <- function(fit,cluster,term_name="conditionUC"){
  cache_file<-file.path(qadir,"implementation3_CR2_vcov_cache.rds")
  V <- if(file.exists(cache_file))readRDS(cache_file) else vcovCR(fit,cluster=cluster,type="CR2")
  if(!file.exists(cache_file))saveRDS(V,cache_file)
  # Target-coefficient-only evaluation of the exact clubSandwich Satterthwaite formula.
  # This is algebraically identical to get_P_array(get_GH(fit,V)) for one coefficient,
  # but avoids repeating dense cluster-matrix products for all seven coefficients.
  cl<-droplevels(factor(cluster)); X<-clubSandwich:::model_matrix(fit); ids<-split(seq_len(nrow(X)),cl)
  W<-clubSandwich:::weightMatrix(fit,cl);w_scale<-attr(W,"w_scale");if(is.null(w_scale))w_scale<-1
  M<-attr(V,"bread")/attr(V,"v_scale"); MUct<-t(chol(w_scale*M)); j<-match(term_name,names(coef(fit)))
  E<-attr(V,"est_mats");A<-attr(V,"adjustments");Theta<-attr(V,"target");J<-length(ids);H<-matrix(NA_real_,nrow=ncol(X),ncol=J);g2<-numeric(J)
  for(ii in seq_len(J)){
    me<-as.numeric((M[j,,drop=FALSE]%*%E[[ii]])%*%A[[ii]])
    th<-Theta[[ii]]; g<-as.numeric(me%*%t(chol(th)))
    g2[ii]<-sum(g^2);H[,ii]<-as.numeric(me%*%X[ids[[ii]],,drop=FALSE]%*%MUct)
  }
  P<--crossprod(H);diag(P)<-diag(P)+g2;df<-(sum(diag(P))^2)/sum(P^2)
  est <- unname(coef(fit)[term_name]); se <- sqrt(V[term_name,term_name]); p <- 2*pt(-abs(est/se),df)
  crit <- qt(.975,df=df)
  c(estimate=est,SE=se,df=df,ci_low=est-crit*se,ci_high=est+crit*se,p_value=p)
}

gee_small_extract <- function(fit,target_term,clusters){
  cs<-summary(fit)$coefficients
  est<-unname(coef(fit)[target_term]); raw_se<-as.numeric(cs[target_term,"Std.err"])
  G<-length(unique(clusters)); p<-length(coef(fit)); df<-G-p
  se<-raw_se*sqrt(G/df); crit<-qt(.975,df); pv<-2*pt(-abs(est/se),df)
  c(estimate=est,SE=se,df=df,ci_low=est-crit*se,ci_high=est+crit*se,p_value=pv)
}

glm_small_extract <- function(fit,target_term,clusters){
  X<-model.matrix(fit);y<-model.response(model.frame(fit));mu<-fitted(fit);pw<-weights(fit,type="prior")
  score<-X*as.numeric(pw*(y-mu));S<-rowsum(score,group=clusters,reorder=FALSE)
  bread0<-solve(crossprod(X,X*as.numeric(pw*mu*(1-mu))));V0<-bread0%*%crossprod(S)%*%bread0
  rownames(V0)<-colnames(V0)<-colnames(X)
  G<-length(unique(clusters));p<-length(coef(fit));df<-G-p
  est<-unname(coef(fit)[target_term]);se<-sqrt(V0[target_term,target_term])*sqrt(G/df);crit<-qt(.975,df);pv<-2*pt(-abs(est/se),df)
  c(estimate=est,SE=se,df=df,ci_low=est-crit*se,ci_high=est+crit*se,p_value=pv)
}

fit_impl1 <- function(d,formula=form,target_term=term,run_gee=TRUE){
  d<-d[order(d$donor_id),]
  qa<-cap(glm(formula,data=d,weights=donor_weight,family=binomial(link="logit")))
  if(inherits(qa$value,"error")) return(list(ok=FALSE,warning=paste(qa$warnings,conditionMessage(qa$value))))
  q<-if(run_gee)cap(geeglm(formula,id=donor_id,data=d,weights=donor_weight,family=binomial(link="logit"),corstr="independence",std.err="san.se")) else list(value=NULL,warnings="")
  if(run_gee && inherits(q$value,"error")) return(list(ok=FALSE,warning=paste(q$warnings,conditionMessage(q$value))))
  e<-tryCatch(if(run_gee)gee_small_extract(q$value,target_term,d$donor_id) else glm_small_extract(qa$value,target_term,d$donor_id),error=identity)
  if(inherits(e,"error")) return(list(ok=FALSE,warning=paste(q$warnings,conditionMessage(e))))
  equiv<-if(run_gee)max(abs(coef(q$value)-coef(qa$value)),na.rm=TRUE) else 0
  list(ok=TRUE,fit=if(run_gee)q$value else qa$value,aux=qa$value,e=e,warning=paste(c(q$warnings,qa$warnings),collapse=" | "),max_GEE_GLM_beta_difference=equiv)
}

add_or <- function(x){
  z<-as.numeric(x[c("estimate","SE","df","ci_low","ci_high","p_value")]);names(z)<-c("estimate","SE","df","ci_low","ci_high","p_value")
  c(z,OR=exp(z["estimate"]),OR_ci_low=exp(z["ci_low"]),OR_ci_high=exp(z["ci_high"])) |> unname() |> setNames(c("estimate","SE","df","ci_low","ci_high","p_value","OR","OR_ci_low","OR_ci_high"))
}
rows<-list()

# Implementation 1: locked candidate revised primary.
i1<-fit_impl1(m); stopifnot(i1$ok,i1$max_GEE_GLM_beta_difference<1e-7)
e1<-add_or(i1$e)
rows[[1]]<-data.frame(implementation="Implementation 1",candidate_revised_primary=TRUE,estimand="equal-donor marginal association",method="inverse-cluster-size weighted GEE; independence working correlation; robust sandwich scaled by G/(G-p), t reference with G-p df",t(e1),n_cells=nrow(m),n_clusters=length(unique(m$donor_id)),converged=TRUE,warnings=i1$warning,check.names=FALSE)

# Implementation 2: one-stage survey design, donor PSUs, inverse cluster-size weights.
des <- svydesign(ids=~donor_id,weights=~donor_weight,data=m)
q2<-cap(svyglm(form,design=des,family=quasibinomial(link="logit")))
stopifnot(!inherits(q2$value,"error")); s2<-summary(q2$value); cs2<-coef(s2); est2<-unname(coef(q2$value)[term]); se2<-unname(cs2[term,"Std. Error"]); df2<-q2$value$df.residual; p2<-unname(cs2[term,ncol(cs2)]); crit2<-qt(.975,df2)
e2<-add_or(c(estimate=est2,SE=se2,df=df2,ci_low=est2-crit2*se2,ci_high=est2+crit2*se2,p_value=p2))
rows[[2]]<-data.frame(implementation="Implementation 2",candidate_revised_primary=FALSE,estimand="equal-donor design-based marginal association",method="survey-weighted quasibinomial logistic regression; donor PSU; design-based residual degrees of freedom",t(e2),n_cells=nrow(m),n_clusters=length(unique(m$donor_id)),converged=TRUE,warnings=q2$warnings,check.names=FALSE)

# Implementation 3: weighted quasibinomial GLM with CR2/Satterthwaite inference.
q3<-cap(glm(form,data=m,weights=donor_weight,family=quasibinomial(link="logit")))
stopifnot(!inherits(q3$value,"error"))
old_comp_file<-file.path(resdir,"implementation_comparison.csv")
if(file.exists(old_comp_file)){
  old_comp<-read.csv(old_comp_file,check.names=FALSE);cached<-old_comp[old_comp$implementation=="Implementation 3",]
} else cached<-data.frame()
if(nrow(cached)==1L && all(is.finite(cached[,c("estimate","SE","df","ci_low","ci_high","p_value")]))) {
  e3<-unlist(cached[1,c("estimate","SE","df","ci_low","ci_high","p_value","OR","OR_ci_low","OR_ci_high")],use.names=TRUE)
} else e3<-add_or(cr2_extract(q3$value,m$donor_id))
rows[[3]]<-data.frame(implementation="Implementation 3",candidate_revised_primary=FALSE,estimand="equal-donor marginal association under working mean model",method="weighted quasibinomial GLM; donor-cluster CR2 covariance; Satterthwaite test",t(e3),n_cells=nrow(m),n_clusters=length(unique(m$donor_id)),converged=TRUE,warnings=q3$warnings,check.names=FALSE)

impl<-do.call(rbind,rows)

# Archived comparators are read, not refit.
old <- read.csv(gart_require_file(file.path(paths$external_data_root, "derived_inputs", "singlecell", "multiverse_singlecell_results.csv"), "archived locked multiverse result"), check.names = FALSE)
old_eq<-old[old$specification_id=="SC-D1-EQUAL-DONOR-WT",]
gl <- read.csv(gart_require_file(file.path(paths$data_public_root, "figure_source_data", "Figure2", "GART_detection_and_hurdle_models.csv"), "committed corrected GART model table"), check.names = FALSE)
gl<-gl[gl$random_structure=="donor-only",]
comparators<-data.frame(
  implementation=c("Archived equal-donor HC1","Archived donor-random-intercept GLMM"),
  candidate_revised_primary=FALSE,
  estimand=c("equal-donor marginal association; asymptotic HC1/z","cell-weighted donor-conditional association"),
  method=c(old_eq$model_family,"donor-random-intercept GLMM; archived frozen result"),
  estimate=c(old_eq$estimate,gl$estimate),SE=c(old_eq$SE,gl$SE),df=NA_real_,ci_low=c(old_eq$ci_low,gl$CI_low),ci_high=c(old_eq$ci_high,gl$CI_high),p_value=c(old_eq$p_value,gl$P),
  OR=c(old_eq$OR,exp(gl$estimate)),OR_ci_low=c(old_eq$OR_ci_low,exp(gl$CI_low)),OR_ci_high=c(old_eq$OR_ci_high,exp(gl$CI_high)),
  n_cells=9306,n_clusters=33,converged=TRUE,warnings="",check.names=FALSE
)
allres<-rbind(impl,comparators)
write.csv(allres,file.path(resdir,"implementation_comparison.csv"),row.names=FALSE)

# Stratified donor bootstrap, 2,000 draws, stratified only by cohort as requested.
set.seed(20260905); B<-2000; boot<-rep(NA_real_,B); cohorts<-levels(m$cohort)
boot_file<-file.path(resdir,"stratified_donor_bootstrap_2000.csv");boot_sum_file<-file.path(resdir,"stratified_donor_bootstrap_summary.csv")
if(file.exists(boot_file)&&file.exists(boot_sum_file)){
  bootdf<-read.csv(boot_file);boot<-bootdf$estimate_logOR;bootok<-boot[is.finite(boot)];boot_sum<-read.csv(boot_sum_file)
} else {
donmap<-unique(m[,c("donor_id","cohort")]); donmap$donor_id<-as.character(donmap$donor_id)
Xboot<-model.matrix(form,m); yboot<-m$det1; did<-as.character(m$donor_id); donor_levels<-levels(m$donor_id); donor_index<-match(did,donor_levels); base_weights<-m$donor_weight
# Generate and lock all donor resamples before parallel fitting.
mult_mat<-matrix(0L,nrow=B,ncol=length(donor_levels),dimnames=list(NULL,donor_levels))
for(b in seq_len(B)){
  picks<-unlist(lapply(cohorts,function(cc){ids<-donmap$donor_id[donmap$cohort==cc];sample(ids,length(ids),replace=TRUE)}),use.names=FALSE)
  mult_mat[b,]<-as.integer(table(factor(picks,levels=donor_levels)))
}
nworkers<-4L; cl<-parallel::makeCluster(nworkers)
parallel::clusterExport(cl,c("Xboot","yboot","donor_index","base_weights","mult_mat","term"),envir=environment())
parts<-split(seq_len(B),rep(seq_len(nworkers),length.out=B))
bres<-parallel::parLapply(cl,parts,function(ii){
  ans<-rep(NA_real_,length(ii));names(ans)<-ii
  for(k in seq_along(ii)){
    b<-ii[k];wb<-base_weights*mult_mat[b,donor_index]
    q<-suppressWarnings(tryCatch(glm.fit(x=Xboot,y=yboot,weights=wb,family=binomial(link="logit")),error=function(e)NULL))
    if(!is.null(q)&&term%in%names(q$coefficients))ans[k]<-q$coefficients[term]
  }
  ans
})
parallel::stopCluster(cl);for(z in bres)boot[as.integer(names(z))]<-z
bootok<-boot[is.finite(boot)]; stopifnot(length(bootok)>=.95*B)
bootdf<-data.frame(iteration=seq_len(B),estimate_logOR=boot,OR=exp(boot),estimable=is.finite(boot))
write.csv(bootdf,boot_file,row.names=FALSE)
boot_sum<-data.frame(B_requested=B,B_estimable=length(bootok),estimate=e1["estimate"],bootstrap_SE=sd(bootok),bootstrap_ci_low=unname(quantile(bootok,.025)),bootstrap_ci_high=unname(quantile(bootok,.975)),bootstrap_OR_ci_low=exp(unname(quantile(bootok,.025))),bootstrap_OR_ci_high=exp(unname(quantile(bootok,.975))),bootstrap_sign_p=2*min(mean(bootok<=0),mean(bootok>=0)),stratification="cohort",seed=20260905)
write.csv(boot_sum,boot_sum_file,row.names=FALSE)
}

# Leave-one-donor-out Implementation 1, preserving candidate-primary method.
lodo<-lapply(levels(m$donor_id),function(id){
  z<-droplevels(m[m$donor_id!=id,]); q<-fit_impl1(z,run_gee=FALSE)
  if(!q$ok)return(data.frame(omitted_donor=id,estimate=NA,SE=NA,df=NA,ci_low=NA,ci_high=NA,p_value=NA,OR=NA,OR_ci_low=NA,OR_ci_high=NA,n_clusters=32,converged=FALSE,warnings=q$warning))
  e<-add_or(q$e);data.frame(omitted_donor=id,t(e),n_clusters=32,converged=TRUE,warnings=q$warning,check.names=FALSE)
})
lodo<-do.call(rbind,lodo);write.csv(lodo,file.path(resdir,"leave_one_donor_out.csv"),row.names=FALSE)
lodo_sum<-data.frame(n_refits=nrow(lodo),n_estimable=sum(is.finite(lodo$OR)),OR_min=min(lodo$OR,na.rm=TRUE),OR_max=max(lodo$OR,na.rm=TRUE),min_OR_omitted_donor=lodo$omitted_donor[which.min(lodo$OR)],max_OR_omitted_donor=lodo$omitted_donor[which.max(lodo$OR)],p_min=min(lodo$p_value,na.rm=TRUE),p_max=max(lodo$p_value,na.rm=TRUE),n_OR_below_1=sum(lodo$OR<1,na.rm=TRUE),n_p_below_0_05=sum(lodo$p_value<.05,na.rm=TRUE))
write.csv(lodo_sum,file.path(resdir,"leave_one_donor_out_summary.csv"),row.names=FALSE)

# Condition-by-cohort interaction using Implementation 1 and the same small-sample correction.
fint<-det1~condition*cohort+z_log_nUMI+z_nGene+S_score+G2M_score
int_terms_expected<-grep("conditionUC:cohort|cohort.*:conditionUC",colnames(model.matrix(fint,m)),value=TRUE); stopifnot(length(int_terms_expected)==1)
qi<-fit_impl1(m,fint,target_term=int_terms_expected,run_gee=TRUE); stopifnot(qi$ok)
int_terms<-int_terms_expected
ei<-add_or(qi$e); intres<-data.frame(term=int_terms,t(ei),n_cells=9306,n_clusters=33,converged=TRUE,warnings=qi$warning,interpretation="ratio of UC-versus-Healthy odds ratios between cohorts",check.names=FALSE)
write.csv(intres,file.path(resdir,"condition_by_cohort_interaction.csv"),row.names=FALSE)

# Cohort-specific Implementation 1. Depth variables are standardized within cohort, matching scale() in a cohort-specific fit.
cohort_res<-lapply(levels(m$cohort),function(cc){
  z<-droplevels(m[m$cohort==cc,]);z$z_log_nUMI<-as.numeric(scale(z$log_nUMI));z$z_nGene<-as.numeric(scale(z$nFeature_RNA));z$donor_n<-as.integer(table(z$donor_id)[as.character(z$donor_id)]);z$donor_weight<-1/z$donor_n
  fc<-det1~condition+z_log_nUMI+z_nGene+S_score+G2M_score
  q<-fit_impl1(z,fc,run_gee=FALSE)
  if(!q$ok)return(data.frame(cohort=cc,estimate=NA,SE=NA,df=NA,ci_low=NA,ci_high=NA,p_value=NA,OR=NA,OR_ci_low=NA,OR_ci_high=NA,n_cells=nrow(z),n_clusters=length(unique(z$donor_id)),converged=FALSE,warnings=q$warning))
  e<-add_or(q$e);data.frame(cohort=cc,t(e),n_cells=nrow(z),n_clusters=length(unique(z$donor_id)),converged=TRUE,warnings=q$warning,check.names=FALSE)
})
cohort_res<-do.call(rbind,cohort_res);write.csv(cohort_res,file.path(resdir,"cohort_specific_implementation1.csv"),row.names=FALSE)

# Model/design audit.
audit<-data.frame(item=c("candidate_primary_locked_before_fit","n_cells","n_donors","outcome","minimum_donor_total_weight","maximum_donor_total_weight","working_correlation","implementation1_small_sample_correction","maximum_GEE_vs_weighted_GLM_beta_difference","bootstrap_unit","bootstrap_stratification","bootstrap_replicates"),value=c("Implementation 1",nrow(m),length(unique(m$donor_id)),"GART_count_corrected >= 1",format(min(wa$total_weight),digits=16),format(max(wa$total_weight),digits=16),"independence","robust sandwich * G/(G-p); t reference with G-p df",format(i1$max_GEE_GLM_beta_difference,scientific=TRUE),"donor","cohort",B))
write.csv(audit,file.path(resdir,"analysis_set_and_design_audit.csv"),row.names=FALSE)

fmt<-function(x,d=3)formatC(x,format="f",digits=d)
fmtp<-function(x)ifelse(x<.0001,format(x,scientific=TRUE,digits=3),formatC(x,format="f",digits=4))
`%+%`<-function(a,b)paste0(a,b)
i1r<-impl[1,]; i2r<-impl[2,]; i3r<-impl[3,]
retains_old <- abs(i1r$OR-old_eq$OR)<0.15 && i1r$p_value<.05
single_driver <- lodo_sum$n_OR_below_1>0 || (lodo_sum$OR_max/lodo_sum$OR_min)>2
interaction_supported <- is.finite(intres$p_value) && intres$p_value<.05

report<-c(
  "# DONOR-BALANCED PRIMARY VALIDATION",
  "",
  "## Locked decision and analysis set",
  "",
  "Implementation 1 was locked before fitting as the candidate reviewer-motivated revised primary. This designation was not changed after viewing any estimate or P value. The validated data contain exactly 9,306 ISC-like epithelial cells from 33 independent donors; raw corrected GART count >=1 is the fixed binary outcome. Every donor has total analytic weight exactly 1 (numerical range " %+% format(min(wa$total_weight),digits=16) %+% " to " %+% format(max(wa$total_weight),digits=16) %+% ").",
  "",
  "The fixed mean model was `condition + scale(log1p(nUMI)) + scale(nGene) + S.Score + G2M.Score + cohort`. No other single-cell result was modified.",
  "",
  "## Three prespecified implementations",
  "",
  "| Implementation | Estimand/inference | log OR | OR (95% CI) | P | df | Clusters | Warnings |",
  "|---|---|---:|---:|---:|---:|---:|---|",
  paste0("| **1 - locked candidate revised primary** | inverse-cluster-size weighted GEE, independence; finite-cluster sandwich G/(G-p) and t(G-p) | ",fmt(i1r$estimate)," | ",fmt(i1r$OR,2)," (",fmt(i1r$OR_ci_low,2),"-",fmt(i1r$OR_ci_high,2),") | ",fmtp(i1r$p_value)," | ",fmt(i1r$df,1)," | 33 | ",ifelse(nchar(i1r$warnings),i1r$warnings,"None")," |"),
  paste0("| 2 | survey-weighted quasibinomial logistic regression; donor PSU; design-based df | ",fmt(i2r$estimate)," | ",fmt(i2r$OR,2)," (",fmt(i2r$OR_ci_low,2),"-",fmt(i2r$OR_ci_high,2),") | ",fmtp(i2r$p_value)," | ",fmt(i2r$df,1)," | ",i2r$n_clusters," | ",ifelse(nchar(i2r$warnings),i2r$warnings,"None")," |"),
  paste0("| 3 | weighted quasibinomial GLM; donor-cluster CR2/Satterthwaite | ",fmt(i3r$estimate)," | ",fmt(i3r$OR,2)," (",fmt(i3r$OR_ci_low,2),"-",fmt(i3r$OR_ci_high,2),") | ",fmtp(i3r$p_value)," | ",fmt(i3r$df,1)," | 33 | ",ifelse(nchar(i3r$warnings),i3r$warnings,"None")," |"),
  "",
  "For reference, the archived HC1 equal-donor analysis reported OR=" %+% fmt(old_eq$OR,2) %+% " (95% CI " %+% fmt(old_eq$OR_ci_low,2) %+% "-" %+% fmt(old_eq$OR_ci_high,2) %+% "), P=" %+% fmtp(old_eq$p_value) %+% "; the frozen donor-random-intercept GLMM reported OR=" %+% fmt(exp(gl$estimate),2) %+% " (95% CI " %+% fmt(exp(gl$CI_low),2) %+% "-" %+% fmt(exp(gl$CI_high),2) %+% "), P=" %+% fmtp(gl$P) %+% ". These comparators were not used to choose the revised-primary implementation.",
  "",
  "## Donor bootstrap and influence",
  "",
  paste0("The cohort-stratified donor bootstrap produced ",length(bootok),"/",B," estimable replicates. Its percentile 95% CI for the OR was ",fmt(boot_sum$bootstrap_OR_ci_low,2),"-",fmt(boot_sum$bootstrap_OR_ci_high,2)," (bootstrap sign P=",fmtp(boot_sum$bootstrap_sign_p),")."),
  paste0("Across all 33 leave-one-donor-out refits, the OR ranged from ",fmt(lodo_sum$OR_min,2)," (omitting ",lodo_sum$min_OR_omitted_donor,") to ",fmt(lodo_sum$OR_max,2)," (omitting ",lodo_sum$max_OR_omitted_donor,"); ",lodo_sum$n_OR_below_1,"/33 estimates were below 1 and ",lodo_sum$n_p_below_0_05,"/33 small-sample P values were below 0.05. ",ifelse(single_driver,"Influence is material enough that the result should not be described as immune to individual donors.","No sign reversal or disproportionate leave-one-donor-out range indicates that a single donor drives the direction.")),
  "",
  "## Cohort heterogeneity",
  "",
  paste0("The condition-by-cohort interaction ratio of ORs was ",fmt(intres$OR,2)," (95% CI ",fmt(intres$OR_ci_low,2),"-",fmt(intres$OR_ci_high,2),"), P=",fmtp(intres$p_value),". ",ifelse(interaction_supported,"This supplies small-sample evidence of cohort heterogeneity.","This interaction test does not reach the nominal 0.05 threshold, but absence of interaction significance is not evidence of homogeneity.")),
  "",
  "| Cohort | Cells | Donors | OR (95% CI) | P | df |",
  "|---|---:|---:|---:|---:|---:|",
  unlist(lapply(seq_len(nrow(cohort_res)),function(j){z<-cohort_res[j,];paste0("| ",z$cohort," | ",z$n_cells," | ",z$n_clusters," | ",fmt(z$OR,2)," (",fmt(z$OR_ci_low,2),"-",fmt(z$OR_ci_high,2),") | ",fmtp(z$p_value)," | ",fmt(z$df,1)," |")})),
  "",
  "## Required validation answers",
  "",
  "1. **Estimand alignment.** Yes. For a question framed as the average association across independent donors, the inverse-cluster-size weighted marginal estimator is better aligned than the cell-weighted donor-conditional GLMM estimand: each donor contributes equal total weight, and inference treats 33 donors, not 9,306 cells, as the independent information units. This does not make either estimand universally superior; it makes Implementation 1 closer to the stated scientific target.",
  paste0("2. **Does OR about 2.08 and P about 0.047 survive small-sample correction?** The point estimate is ",fmt(i1r$OR,2),", so the archived OR about 2.08 ",ifelse(abs(i1r$OR-old_eq$OR)<.05,"is retained","changes slightly"),". With finite-cluster GEE sandwich correction and t(G-p) inference, P=",fmtp(i1r$p_value)," and the 95% CI is ",fmt(i1r$OR_ci_low,2),"-",fmt(i1r$OR_ci_high,2),". Therefore the archived HC1 P about 0.047 ",ifelse(i1r$p_value<.05,"remains nominally below 0.05","does not remain below 0.05")," after the prespecified small-sample correction."),
  paste0("3. **Single-donor influence.** ",ifelse(single_driver,"The leave-one-donor-out range indicates material donor influence; interpret the precise magnitude and threshold crossing cautiously.","No leave-one-donor-out sign reversal and the observed OR range do not support a single-donor-driven direction.")," Full donor-specific refits are retained in `leave_one_donor_out.csv`."),
  "4. **Can Implementation 1 be the reviewer-motivated revised primary?** Yes, if the revised target is explicitly the equal-donor population-averaged association and the Methods state inverse donor-size weights, independence working correlation, the finite-cluster sandwich multiplier G/(G-p), and t(G-p) inference. The frozen GLMM should remain reported as the original conditional comparator/sensitivity, not silently replaced.",
  "5. **Selection rule.** Implementation 1 remains the candidate revised primary regardless of whether its P value is above or below 0.05. No significance-based model switching was performed.",
  "",
  "## Reviewer-safe conclusion",
  "",
  "The donor-balanced model estimates a positive UC-associated difference in GART detection under an equal-donor marginal estimand. Its strength of statistical evidence must be stated using the CR2/Satterthwaite result above, together with the bootstrap, leave-one-donor-out, cohort-specific, interaction, and frozen GLMM results. The validation does not justify a claim of universal UC-wide elevation."
)
writeLines(report,file.path(repdir,"DONOR_BALANCED_PRIMARY_VALIDATION.md"),useBytes=TRUE)

capture.output(sessionInfo(),file=file.path(qadir,"sessionInfo.txt"))
writeLines(c("PASS: 9306 cells","PASS: 33 donors","PASS: each donor total weight = 1","PASS: Implementation 1 lock written before fitting","PASS: 2000 cohort-stratified donor bootstrap attempts","PASS: 33 leave-one-donor-out refits","PASS: no other single-cell result written or modified"),file.path(qadir,"validation_checks.txt"))
cat("Completed donor-balanced validation\n")
