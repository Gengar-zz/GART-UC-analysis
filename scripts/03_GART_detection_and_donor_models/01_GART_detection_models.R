# Public-release copy: file paths were converted to repository-relative paths.
# Statistical logic, parameters, thresholds, and seeds are unchanged from the archived final script.
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1]))); source(file.path(dirname(script_dir), "02_singlecell_preprocessing", "common.R"))
suppressPackageStartupMessages({library(lme4)})
m <- read_corrected()$metadata; m$condition <- relevel(factor(m$condition),"Healthy");m$cohort<-factor(m$cohort);m$donor_id<-factor(m$donor_id);m$sample_id<-factor(m$sample_id);m$log_nUMI<-log1p(m$nCount_RNA)
results<-list(); diags<-list()
run_bin <- function(outcome,label,data=m){
  for(st in c("donor-only","donor/sample")){
    ff<-if(st=="donor-only")as.formula(paste0(outcome,"~condition+scale(log_nUMI)+scale(nFeature_RNA)+S_score+G2M_score+cohort+(1|donor_id)")) else as.formula(paste0(outcome,"~condition+scale(log_nUMI)+scale(nFeature_RNA)+S_score+G2M_score+cohort+(1|donor_id/sample_id)"))
    q<-capture_fit(glmer(ff,data=data,family=binomial(),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))); results[[length(results)+1]]<<-fit_effect(q$fit,"conditionUC",label,st,nrow(data),length(unique(data$donor_id)),q$warning)
  }
}
m$det_ge1<-as.integer(m$GART_count_corrected>=1);m$det_ge2<-as.integer(m$GART_count_corrected>=2);run_bin("det_ge1","Primary raw GART count >=1");run_bin("det_ge2","Sensitivity raw GART count >=2")
for(co in levels(m$cohort)){d<-droplevels(m[m$cohort==co,]);q<-capture_fit(glmer(det_ge1~condition+scale(log_nUMI)+scale(nFeature_RNA)+S_score+G2M_score+(1|donor_id),data=d,family=binomial(),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))));results[[length(results)+1]]<-fit_effect(q$fit,"conditionUC",paste0("Cohort-specific >=1: ",co),"donor-only",nrow(d),length(unique(d$donor_id)),q$warning)}
m$log_GART_all<-log1p(m$GART_count_corrected);q<-capture_fit(lmer(log_GART_all~condition+scale(log_nUMI)+scale(nFeature_RNA)+S_score+G2M_score+cohort+(1|donor_id),data=m,REML=FALSE));results[[length(results)+1]]<-fit_effect(q$fit,"conditionUC","Continuous log1p raw GART count","donor-only",nrow(m),length(unique(m$donor_id)),q$warning)
hp<-m[m$det_ge1==1,];q<-capture_fit(lmer(log_GART_all~condition+scale(log_nUMI)+scale(nFeature_RNA)+S_score+G2M_score+cohort+(1|donor_id),data=hp,REML=FALSE));results[[length(results)+1]]<-fit_effect(q$fit,"conditionUC","Hurdle positive-expression component","donor-only",nrow(hp),length(unique(hp$donor_id)),q$warning)
res<-do.call(rbind,results);res$primary<-res$analysis=="Primary raw GART count >=1"&res$random_structure=="donor-only";write.csv(res,file.path(output_root,"RESULTS","GART_detection_and_hurdle_models.csv"),row.names=FALSE);write.csv(res[,c("analysis","random_structure","n_cells","n_donors","isSingular","convergence_warning","donor_variance","sample_variance","residual_variance")],file.path(output_root,"RESULTS","GART_model_diagnostics.csv"),row.names=FALSE)
ds<-aggregate(m$det_ge1,list(cohort=m$cohort,condition=m$condition,donor_id=m$donor_id),mean);names(ds)[4]<-"detection_fraction";write.csv(ds,file.path(output_root,"RESULTS","GART_donor_fraction.csv"),row.names=FALSE)
tests<-do.call(rbind,lapply(split(ds,ds$cohort),function(d){uc<-d$detection_fraction[d$condition=="UC"];hc<-d$detection_fraction[d$condition=="Healthy"];tt<-t.test(uc,hc);data.frame(cohort=d$cohort[1],n_Healthy=length(hc),n_UC=length(uc),difference_UC_minus_Healthy=mean(uc)-mean(hc),CI_low=tt$conf.int[1],CI_high=tt$conf.int[2],P=tt$p.value)}));write.csv(tests,file.path(output_root,"RESULTS","GART_donor_fraction_tests.csv"),row.names=FALSE)
