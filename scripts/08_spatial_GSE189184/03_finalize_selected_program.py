# Public-release copy: file paths were converted to repository-relative paths.
# Statistical logic, parameters, thresholds, and seeds are unchanged from the archived final script.
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"source_data"; REPORTS=ROOT/"reports"
genes=np.array((SRC/"singlecell_SCP259_gene_names.txt").read_text(encoding="utf-8").splitlines())
cells=np.array((SRC/"singlecell_SCP259_cell_names.txt").read_text(encoding="utf-8").splitlines())
meta=pd.read_csv(SRC/"singlecell_SCP259_11donor_metadata.csv").set_index("cell_id_revision").loc[cells].reset_index()
members=pd.read_csv(SRC/"candidate_signature_members.csv")
res=pd.read_csv(SRC/"singlecell_donor_blocked_CV_results.csv")
win=res.loc[res.selected.astype(str).str.lower().eq("true")].iloc[0]
sig=members[members.signature==win.signature].sort_values("rank")
sel_genes=sig.gene.tolist()
X=mmread(SRC/"singlecell_SCP259_11donor_logexpr.mtx").tocsc().astype(np.float32)
g2r={g:i for i,g in enumerate(genes)}
rows=np.array([g2r[g] for g in sel_genes if g in g2r])
N=X.shape[0]
R=np.empty((X.shape[1],len(rows)),dtype=np.float32)
g2pos={r:i for i,r in enumerate(rows)}
for j in range(X.shape[1]):
    lo,hi=X.indptr[j],X.indptr[j+1]
    idx=X.indices[lo:hi]; vals=X.data[lo:hi]
    rr=rankdata(-vals,method="average")
    R[j,:]=1-(N+1)/2/N
    for k,r in enumerate(idx):
        p=g2pos.get(r)
        if p is not None: R[j,p]=1-rr[k]/N
y=(meta.GART_status=="positive").to_numpy()
d=meta.donor_id.astype(str).to_numpy()
keys=np.char.add(d,np.where(y,"|1","|0")); _,inv,c=np.unique(keys,return_inverse=True,return_counts=True); w=1/c[inv]
pos=np.average(R[y],axis=0,weights=w[y]); neg=np.average(R[~y],axis=0,weights=w[~y]); weights=pos-neg
out=sig[sig.gene.isin(sel_genes)].copy(); out["projection_weight"]=[weights[sel_genes.index(g)] for g in out.gene]
out.to_csv(SRC/"selected_signature_gene_list.csv",index=False)

summary={"selected_signature":win.signature,"selected_method":win.method,"n_genes":len(sel_genes),"selection_score":float(win.SC_SELECTION_SCORE)}
(SRC/"singlecell_selection_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")

cols=["signature","method","n_genes","heldout_donor_mean_AUC","donor_direction_consistency","cell_cycle_adjusted_retention","within_phase_consistency","Brier_score","feature_stability","random_signature_specificity","abs_QC_correlation","SC_SELECTION_SCORE","selected"]
header="| "+" | ".join(cols)+" |"
sep="|"+"|".join(["---"]*len(cols))+"|"
rows_md=[]
for _,r in res[cols].iterrows():
    rows_md.append("| "+" | ".join(f"{v:.4f}" if isinstance(v,(float,np.floating)) else str(v) for v in r)+" |")
report=["# Single-cell model selection report","","Selection used only the frozen SCP259 11-paired-donor single-cell dataset. Spatial outcomes were not accessed.","",f"Winner: `{win.signature}` with `{win.method}`; SC_SELECTION_SCORE={win.SC_SELECTION_SCORE:.6f}; {len(sel_genes)} genes.","","All candidate combinations, including non-selected results, are retained below and in the CSV. Package-unavailable methods use explicitly labelled mathematical implementations; no package-native result is implied.","","The composite score followed the pre-analysis weights. Matched random signatures were evaluated with 1,000 size- and expression-decile-matched sets per fixed signature.","",header,sep,*rows_md]
(REPORTS/"SINGLECELL_MODEL_SELECTION_REPORT.md").write_text("\n".join(report),encoding="utf-8")
print(json.dumps(summary))
