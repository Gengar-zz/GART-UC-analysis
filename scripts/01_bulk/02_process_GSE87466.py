# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 20_BULK_EXTERNAL_COHORT_SCREEN/scripts/extract_analyze_gse87466.py
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from path_config import get_paths, require_file
PATHS = get_paths()
import gzip,re
import numpy as np,pandas as pd
from scipy import stats

R = PATHS.results_root / "bulk"
(R / "results").mkdir(parents=True, exist_ok=True)
p = require_file(PATHS.external_data_root / "bulk" / "GSE87466" / "GSE87466_family.soft.gz", "GSE87466 family SOFT file")
probes=['210005_PM_at','212378_PM_at','212379_PM_at','216990_PM_at','217445_PM_s_at','230766_PM_at']
rows=[]; cur=None; intable=False
with gzip.open(p,'rt',encoding='utf-8',errors='replace') as f:
    for line in f:
        line=line.rstrip('\r\n')
        if line.startswith('^SAMPLE = '):
            if cur: rows.append(cur)
            cur={'GSM_ID':line.split('=',1)[1].strip()}; intable=False
        elif cur is not None and line.startswith('!Sample_title = '): cur['original_label']=line.split('=',1)[1].strip()
        elif cur is not None and line.startswith('!Sample_characteristics_ch1 = '):
            v=line.split('=',1)[1].strip(); cur['characteristics']=cur.get('characteristics','')+(' | ' if cur.get('characteristics') else '')+v
        elif line=='!sample_table_begin': intable=True
        elif line=='!sample_table_end': intable=False
        elif intable and cur is not None:
            z=line.split('\t')
            if z[0] in probes:
                try: cur[z[0]]=float(z[1])
                except: pass
if cur: rows.append(cur)
d=pd.DataFrame(rows)
assert len(d)==108 and all(x in d for x in probes)
d['participant_ID']=d.characteristics.str.extract(r'subject:\s*([^|]+)')[0].str.strip()
d['condition']=np.where(d.characteristics.str.contains('disease: Normal',case=False,na=False),'Healthy','Active UC')
d['GART_expression']=d[probes].mean(axis=1)
d['dataset']='GSE87466'; d['activity']=np.where(d.condition=='Active UC','moderately-to-severely active','not applicable')
d['biopsy_location']='colon; exact segment not recorded in GEO sample metadata'; d['neoplasia_status']='none indicated'; d['treatment']='baseline/pre-treatment collection stated by trial context; treatment variable not present in GEO sample metadata'
d['retained_or_excluded']='retained'; d['exclusion_reason']=''; d['source_URL']='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+d.GSM_ID
h=d.loc[d.condition=='Healthy','GART_expression'].to_numpy(); u=d.loc[d.condition=='Active UC','GART_expression'].to_numpy()
tt=stats.ttest_ind(u,h,equal_var=False)
n1,n0=len(u),len(h); sp=np.sqrt(((n1-1)*u.var(ddof=1)+(n0-1)*h.var(ddof=1))/(n1+n0-2)); dc=(u.mean()-h.mean())/sp; J=1-3/(4*(n1+n0)-9); g=J*dc
vg=(n1+n0)/(n1*n0)+g*g/(2*(n1+n0-2)); se=np.sqrt(vg); ci=(g-1.96*se,g+1.96*se)
s=pd.DataFrame([{'dataset':'GSE87466','platform':'Affymetrix HT HG-U133+ PM','n_Healthy':n0,'n_Active_UC':n1,'mean_Healthy':h.mean(),'SD_Healthy':h.std(ddof=1),'median_Healthy':np.median(h),'Q1_Healthy':np.quantile(h,.25),'Q3_Healthy':np.quantile(h,.75),'mean_Active_UC':u.mean(),'SD_Active_UC':u.std(ddof=1),'median_Active_UC':np.median(u),'Q1_Active_UC':np.quantile(u,.25),'Q3_Active_UC':np.quantile(u,.75),'UC_minus_HC':u.mean()-h.mean(),'Hedges_g':g,'SE':se,'CI_low':ci[0],'CI_high':ci[1],'P_value':tt.pvalue,'test':'two-sided Welch t-test on the pre-specified mean of six unambiguously GART-annotated log2 processed probes','direction':'UC > Healthy' if u.mean()>h.mean() else 'UC < Healthy'}])
cols=['dataset','GSM_ID','participant_ID','original_label','condition','activity','biopsy_location','neoplasia_status','treatment']+probes+['GART_expression','retained_or_excluded','exclusion_reason','source_URL']
d[cols].to_csv(R/'results'/'GSE87466_sample_level.csv',index=False)
pd.DataFrame({'probe_id':probes,'gene_symbol':'GART','included_in_primary_mean':True,'cross_hybridization_warning':'none in official GPL annotation','annotation_source':'GPL13158 embedded in GSE87466 family SOFT'}).to_csv(R/'results'/'GSE87466_GART_probe_audit.csv',index=False)
s.to_csv(R/'results'/'GSE87466_GART_result.csv',index=False); print(s.to_string(index=False)); print(d.condition.value_counts())
# Pre-specified best-annotated probe sensitivity: RefSeq consensus NM_000819.
bp='212378_PM_at'; hb=d.loc[d.condition=='Healthy',bp].to_numpy(); ub=d.loc[d.condition=='Active UC',bp].to_numpy(); btt=stats.ttest_ind(ub,hb,equal_var=False)
bsp=np.sqrt(((len(ub)-1)*ub.var(ddof=1)+(len(hb)-1)*hb.var(ddof=1))/(len(ub)+len(hb)-2)); bg=(1-3/(4*(len(ub)+len(hb))-9))*(ub.mean()-hb.mean())/bsp; bse=np.sqrt((len(ub)+len(hb))/(len(ub)*len(hb))+bg*bg/(2*(len(ub)+len(hb)-2)))
pd.DataFrame([{'dataset':'GSE87466','probe_id':bp,'annotation':'RefSeq consensus NM_000819; GART','UC_minus_HC':ub.mean()-hb.mean(),'Hedges_g':bg,'CI_low':bg-1.96*bse,'CI_high':bg+1.96*bse,'P_value':btt.pvalue,'direction':'UC > Healthy' if ub.mean()>hb.mean() else 'UC < Healthy'}]).to_csv(R/'results'/'GSE87466_best_annotated_probe_sensitivity.csv',index=False)
