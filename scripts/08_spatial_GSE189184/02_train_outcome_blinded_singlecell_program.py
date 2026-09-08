# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 19_SPATIAL_REVISED_PRIMARY_ADJUDICATION/code/03_singlecell_model_selection.py
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
from __future__ import annotations

import json
import math
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from path_config import get_paths, require_file
PATHS = get_paths()

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import rankdata, spearmanr
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import GroupKFold


ROOT = PATHS.work_root / "spatial_GSE189184"
SRC = ROOT / "source_data"
REPORTS = ROOT / "reports"
SEED_MODEL = 20260907
SEED_RANDOM = 20260908
rng = np.random.default_rng(SEED_MODEL)

DEG_PATH = require_file(PATHS.data_public_root / "figure_source_data" / "Figure3" / "SCP259_primary_pseudobulk_DEG.csv", "committed donor-paired pseudobulk DEG table")
HALLMARK_PATH = require_file(PATHS.external_data_root / "gene_sets" / "msigdb" / "h.all.v2025.1.Hs.symbols.gmt", "MSigDB Hallmark GMT")
REACTOME_PATH = require_file(PATHS.external_data_root / "gene_sets" / "msigdb" / "c2.cp.reactome.v2025.1.Hs.symbols.gmt", "MSigDB Reactome GMT")

S_GENES = {
    "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"
}
G2M_GENES = {
    "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"
}
ROI_GENES = {"EPCAM","KRT8","KRT18","KRT19"}
VALIDATION_GENES = {"MYBL2","MKI67","PCNA","PAICS","OLFM4","MIF","CD74"}


def read_gmt(path: Path):
    out = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            out[parts[0]] = set(parts[2:])
    return out


def is_qc_gene(g):
    u = g.upper()
    return u.startswith("MT-") or u.startswith("RPL") or u.startswith("RPS") or u.startswith("HBA") or u.startswith("HBB")


genes = np.array((SRC / "singlecell_SCP259_gene_names.txt").read_text(encoding="utf-8").splitlines())
cells = np.array((SRC / "singlecell_SCP259_cell_names.txt").read_text(encoding="utf-8").splitlines())
meta = pd.read_csv(require_file(SRC / "singlecell_SCP259_11donor_metadata.csv", "upstream single-cell training metadata"))
meta = meta.set_index("cell_id_revision").loc[cells].reset_index()
X = mmread(SRC / "singlecell_SCP259_11donor_logexpr.mtx").tocsr().astype(np.float32)
assert X.shape == (len(genes), len(cells))

deg = pd.read_csv(DEG_PATH)
hallmark = read_gmt(HALLMARK_PATH)
reactome = read_gmt(REACTOME_PATH)
source_sets = {
    "HALLMARK_E2F_TARGETS": hallmark["HALLMARK_E2F_TARGETS"],
    "HALLMARK_G2M_CHECKPOINT": hallmark["HALLMARK_G2M_CHECKPOINT"],
    "HALLMARK_MYC_TARGETS_V1": hallmark["HALLMARK_MYC_TARGETS_V1"],
    "HALLMARK_MYC_TARGETS_V2": hallmark["HALLMARK_MYC_TARGETS_V2"],
    "REACTOME_DNA_REPLICATION": reactome["REACTOME_DNA_REPLICATION"],
    "REACTOME_DNA_REPAIR": reactome["REACTOME_DNA_REPAIR"],
    "REACTOME_PURINE_RIBONUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHESIS": reactome["REACTOME_PURINE_RIBONUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHESIS"],
}
public_union = set().union(*source_sets.values())
expressed = set(genes)
sig_base = deg[(deg.logFC > 0) & (deg.FDR < 0.05)].copy()

def base_excluded(g):
    return g == "GART" or g in ROI_GENES or is_qc_gene(g)

sig_a = [g for g in sig_base.gene if g in public_union and not base_excluded(g)]
sig_b = [g for g in sig_a if g not in S_GENES | G2M_GENES | {"MKI67","PCNA","MYBL2"}]
sig_c_pool = sig_base[~sig_base.gene.map(base_excluded) & ~sig_base.gene.isin(S_GENES | G2M_GENES)].sort_values(["FDR","logFC"], ascending=[True,False])
sig_c = sig_c_pool.gene.head(50).tolist()
feature_pool = [g for g in sig_base.gene if not base_excluded(g) and g not in VALIDATION_GENES]

signatures = {
    "SC_GART_BROAD_PROGRAM": sig_a,
    "SC_GART_CYCLE_EXCLUDED_PROGRAM": sig_b,
    "SC_GART_STABLE_TOP_PROGRAM": sig_c,
}

audit_rows = []
for sig_name, sig_genes in signatures.items():
    selected = set(sig_genes)
    for _, row in deg.iterrows():
        g = row.gene
        pathways = [name for name, members in source_sets.items() if g in members]
        reason = "retained" if g in selected else ""
        if not reason:
            if row.logFC <= 0: reason = "excluded_nonpositive_log2FC"
            elif row.FDR >= 0.05: reason = "excluded_FDR_ge_0.05"
            elif g == "GART": reason = "excluded_GART"
            elif g in ROI_GENES: reason = "excluded_ROI_marker"
            elif is_qc_gene(g): reason = "excluded_QC_gene"
            elif sig_name == "SC_GART_BROAD_PROGRAM" and not pathways: reason = "excluded_not_in_predefined_public_union"
            elif sig_name == "SC_GART_CYCLE_EXCLUDED_PROGRAM" and g in S_GENES | G2M_GENES | {"MKI67","PCNA","MYBL2"}: reason = "excluded_cell_cycle"
            elif sig_name == "SC_GART_STABLE_TOP_PROGRAM": reason = "excluded_not_in_locked_top50"
            else: reason = "excluded_by_signature_rule"
        if g in selected or (row.logFC > 0 and row.FDR < 0.05):
            audit_rows.append({
                "signature": sig_name, "gene": g, "matched_spatial": g in expressed,
                "source_pathways": ";".join(pathways), "singlecell_log2FC": row.logFC,
                "singlecell_P": row.raw_P, "singlecell_FDR": row.FDR,
                "retained": g in selected, "reason": reason,
            })
pd.DataFrame(audit_rows).to_csv(SRC / "candidate_signature_gene_audit.csv", index=False)

# Restrict ranks/expression to all candidate and feature-pool genes.
candidate_union = sorted((set(feature_pool) | set().union(*map(set, signatures.values()))) & expressed)
gene_to_row = {g:i for i,g in enumerate(genes)}
sel_rows = np.array([gene_to_row[g] for g in candidate_union], dtype=int)
Xsel = X[sel_rows, :].T.tocsr()
n_cells, n_sel = Xsel.shape
N = X.shape[0]

# Per-cell ranks among the complete expressed universe, retaining selected genes only.
rank_sel = np.full((n_cells, n_sel), (N + 1) / 2, dtype=np.float32)
global_to_sel = np.full(N, -1, dtype=np.int32)
global_to_sel[sel_rows] = np.arange(n_sel, dtype=np.int32)
Xcsc = X.tocsc()
for j in range(n_cells):
    lo, hi = Xcsc.indptr[j], Xcsc.indptr[j+1]
    idx, vals = Xcsc.indices[lo:hi], Xcsc.data[lo:hi]
    rr = rankdata(-vals, method="average").astype(np.float32)
    pos = global_to_sel[idx]
    keep = pos >= 0
    rank_sel[j, pos[keep]] = rr[keep]

expr_sel = Xsel.toarray().astype(np.float32)
gene_means = np.asarray(X.mean(axis=1)).ravel()
gene_sqmeans = np.asarray(X.power(2).mean(axis=1)).ravel()
gene_sds = np.sqrt(np.maximum(gene_sqmeans - gene_means**2, 1e-8))

sel_index = {g:i for i,g in enumerate(candidate_union)}

def score_signature(sig_genes, method):
    idx = np.array([sel_index[g] for g in sig_genes if g in sel_index], dtype=int)
    ranks = rank_sel[:, idx]
    m = len(idx)
    if method == "UCell":
        max_rank = 1500
        capped = np.minimum(ranks, max_rank + 1)
        u = capped.sum(axis=1) - m * (m + 1) / 2
        return (1 - u / (m * max_rank)).astype(float)
    if method in {"singscore", "mean_rank"}:
        return (1 - (ranks.mean(axis=1) - 1) / (N - 1)).astype(float)
    if method == "AUCell":
        cutoff = int(round(0.05 * N))
        contrib = np.maximum(cutoff - ranks + 1, 0)
        denom = m * cutoff - m * (m - 1) / 2
        return (contrib.sum(axis=1) / max(denom, 1)).astype(float)
    if method == "ssGSEA":
        pos = N - ranks + 1
        w = np.power(np.maximum(pos, 1), 0.25)
        hit_area = (w * pos).sum(axis=1) / np.maximum(w.sum(axis=1), 1e-12)
        miss_area = (N * (N + 1) / 2 - pos.sum(axis=1)) / (N - m)
        return ((hit_area - miss_area) / N).astype(float)
    if method == "mean_z":
        rows = np.array([gene_to_row[g] for g in sig_genes if g in gene_to_row], dtype=int)
        z = (X[rows, :].T.toarray() - gene_means[rows]) / gene_sds[rows]
        return z.mean(axis=1).astype(float)
    raise ValueError(method)


def balanced_weights(donor, y):
    keys = np.char.add(np.asarray(donor, str), np.where(y==1, "|1", "|0"))
    _, inv, counts = np.unique(keys, return_inverse=True, return_counts=True)
    return 1.0 / counts[inv]


donor = meta.donor_id.astype(str).to_numpy()
y = (meta.GART_status == "positive").astype(int).to_numpy()
s_score = meta.S_score.to_numpy(float)
g2m_score = meta.G2M_score.to_numpy(float)
n_umi = meta.nCount_RNA.to_numpy(float)
n_gene = meta.nFeature_RNA.to_numpy(float)
weights = balanced_weights(donor, y)

phase = np.where((s_score <= 0) & (g2m_score <= 0), "G1", np.where(s_score >= g2m_score, "S", "G2M"))

def donor_differences(score):
    rows=[]
    for d in np.unique(donor):
        m=donor==d
        rows.append((d, score[m & (y==1)].mean()-score[m & (y==0)].mean()))
    return pd.DataFrame(rows,columns=["donor_id","difference"])

def weighted_ols_beta(score, adjusted):
    donor_levels=np.unique(donor)
    cols=[np.ones(n_cells),y.astype(float)]
    if adjusted: cols += [s_score,g2m_score]
    cols += [(donor==d).astype(float) for d in donor_levels[1:]]
    A=np.column_stack(cols)
    sw=np.sqrt(weights/weights.mean())
    beta=np.linalg.lstsq(A*sw[:,None],score*sw,rcond=None)[0]
    return beta[1]

def evaluate_score(score, signature_genes, method_name, signature_name, random_specificity):
    fold_rows=[]; pred=np.zeros(n_cells); aucs=[]
    for d in np.unique(donor):
        test=donor==d; train=~test
        lr=LogisticRegression(C=1.0,solver="lbfgs",max_iter=2000)
        lr.fit(score[train,None],y[train],sample_weight=weights[train])
        pred[test]=lr.predict_proba(score[test,None])[:,1]
        aucs.append(roc_auc_score(y[test],score[test]))
        fold_rows.append({"signature":signature_name,"method":method_name,"heldout_donor":d,"AUC":aucs[-1]})
    dif=donor_differences(score)
    donor_cons=(dif.difference>0).mean()
    unadj=weighted_ols_beta(score,False); adj=weighted_ols_beta(score,True)
    retention=max(0,min(1,adj/unadj)) if unadj>0 else 0
    phase_diffs=[]
    for d in np.unique(donor):
        for p in ("G1","S","G2M"):
            m=(donor==d)&(phase==p)
            if (m&(y==1)).sum()>=5 and (m&(y==0)).sum()>=5:
                phase_diffs.append(score[m&(y==1)].mean()-score[m&(y==0)].mean())
    within_phase=np.mean(np.array(phase_diffs)>0) if phase_diffs else np.nan
    brier=np.average((pred-y)**2,weights=weights)
    calibration=max(0,min(1,1-brier/0.25))
    # Gene feature stability: fraction of signature genes positive in >=80% donors.
    gidx=[sel_index[g] for g in signature_genes if g in sel_index]
    stable=[]
    for gi in gidx:
        gd=[]
        for d in np.unique(donor):
            m=donor==d
            gd.append(expr_sel[m&(y==1),gi].mean()-expr_sel[m&(y==0),gi].mean())
        stable.append(np.mean(np.array(gd)>0)>=0.8)
    feature_stability=float(np.mean(stable)) if stable else 0
    qc_r=max(abs(spearmanr(score,np.log1p(n_umi)).statistic),abs(spearmanr(score,n_gene).statistic))
    selection=(0.30*np.mean(aucs)+0.20*donor_cons+0.15*retention+0.10*within_phase+0.10*calibration+0.10*feature_stability+0.05*random_specificity-(0.05 if qc_r>0.30 else 0))
    return {
        "signature":signature_name,"method":method_name,"implementation":"documented_formula",
        "n_genes":len(signature_genes),"heldout_donor_mean_AUC":np.mean(aucs),"donor_direction_consistency":donor_cons,
        "cell_cycle_adjusted_retention":retention,"within_phase_consistency":within_phase,"Brier_score":brier,
        "calibration_component":calibration,"feature_stability":feature_stability,"random_signature_specificity":random_specificity,
        "abs_QC_correlation":qc_r,"QC_penalty":0.05 if qc_r>0.30 else 0,"SC_SELECTION_SCORE":selection,
        "unadjusted_effect":unadj,"cell_cycle_adjusted_effect":adj,
    }, fold_rows, pred, dif

# Matched random gene-set specificity, assessed with mean-z donor-paired effects.
eligible_random = np.array([g for g in genes if not base_excluded(g) and g not in VALIDATION_GENES])
eligible_rows = np.array([gene_to_row[g] for g in eligible_random])
means = gene_means[eligible_rows]
bins = pd.qcut(pd.Series(means).rank(method="first"), 10, labels=False).to_numpy()
bin_lookup = {g:int(b) for g,b in zip(eligible_random,bins)}
bin_members = {b:eligible_random[bins==b] for b in range(10)}

random_summary=[]
random_nulls={}
for sig_name,sig_genes in signatures.items():
    obs_genes=[g for g in sig_genes if g in gene_to_row]
    obs_rows=np.array([gene_to_row[g] for g in obs_genes])
    obs_z=((X[obs_rows,:].T.toarray()-gene_means[obs_rows])/gene_sds[obs_rows]).mean(axis=1)
    obs_eff=donor_differences(obs_z).difference.mean()
    rrng=np.random.default_rng(SEED_RANDOM+list(signatures).index(sig_name))
    random_sets=[]
    for k in range(1000):
        chosen=[]
        for g in obs_genes:
            pool=bin_members[bin_lookup.get(g,5)]
            chosen.append(rrng.choice(pool))
        random_sets.append(list(dict.fromkeys(chosen)))
    # Sparse signature-by-gene matrix, then z-standardized average scores.
    rows=[]; cols=[]; vals=[]
    for i,gs in enumerate(random_sets):
        for g in gs:
            rows.append(i); cols.append(gene_to_row[g]); vals.append(1/len(gs))
    W=sparse.csr_matrix((vals,(rows,cols)),shape=(1000,N),dtype=np.float32)
    rand_scores=(W@X).toarray().astype(np.float32)
    center=np.asarray(W@(gene_means/gene_sds)).ravel()
    # Correct the matrix term to z units.
    Wz=W@sparse.diags(1/gene_sds)
    rand_scores=(Wz@X).toarray().astype(np.float32)-center[:,None]
    rand_eff=np.zeros(1000)
    for k in range(1000): rand_eff[k]=donor_differences(rand_scores[k]).difference.mean()
    percentile=(np.sum(rand_eff<obs_eff)+0.5*np.sum(rand_eff==obs_eff))/1000
    random_nulls[sig_name]=rand_eff
    random_summary.append({"signature":sig_name,"observed_mean_z_effect":obs_eff,"empirical_percentile":percentile,"n_random":1000,"matched_on":"signature size and gene mean-expression decile"})
pd.DataFrame(random_summary).to_csv(SRC/"singlecell_matched_random_summary.csv",index=False)

results=[]; cv_rows=[]; score_store={}; pred_store={}
base_methods=["UCell","singscore","AUCell","ssGSEA","mean_rank","mean_z"]
for sig_name,sig_genes in signatures.items():
    specificity=next(x["empirical_percentile"] for x in random_summary if x["signature"]==sig_name)
    component={}
    for method in base_methods:
        sc=score_signature(sig_genes,method)
        component[method]=sc
        res,fold,pred,dif=evaluate_score(sc,sig_genes,method,sig_name,specificity)
        results.append(res); cv_rows.extend(fold); score_store[(sig_name,method)]=sc; pred_store[(sig_name,method)]=pred
    ens=np.mean(np.column_stack([rankdata(component[m])/n_cells for m in ["UCell","singscore","AUCell","ssGSEA","mean_rank"]]),axis=1)
    res,fold,pred,dif=evaluate_score(ens,sig_genes,"ensemble_rank",sig_name,specificity)
    results.append(res); cv_rows.extend(fold); score_store[(sig_name,"ensemble_rank")]=ens; pred_store[(sig_name,"ensemble_rank")]=pred

# Elastic-net: equal 100 cells per donor-state, nested donor-group CV for fixed grid.
pool_genes=[g for g in feature_pool if g in sel_index]
pidx=np.array([sel_index[g] for g in pool_genes])
Xp=expr_sel[:,pidx]
train_idx=[]
for d in np.unique(donor):
    for yy in (0,1):
        ids=np.flatnonzero((donor==d)&(y==yy)); train_idx.extend(rng.choice(ids,size=min(100,len(ids)),replace=False))
train_idx=np.array(sorted(train_idx))
enet_pred=np.zeros(n_cells); coef_folds=[]; enet_auc=[]
for d in np.unique(donor):
    outer_train=train_idx[donor[train_idx]!=d]; test=np.flatnonzero(donor==d)
    inner_groups=donor[outer_train]
    best=None
    for C in (0.03,0.1,0.3):
        for l1 in (0.2,0.5,0.8):
            vals=[]
            splitter=GroupKFold(n_splits=min(5,len(np.unique(inner_groups))))
            for tr,va in splitter.split(outer_train,y[outer_train],inner_groups):
                it,iv=outer_train[tr],outer_train[va]
                m=LogisticRegression(C=C,l1_ratio=l1,penalty="elasticnet",solver="saga",max_iter=1500,tol=1e-3,n_jobs=1)
                m.fit(Xp[it],y[it],sample_weight=balanced_weights(donor[it],y[it]))
                vals.append(roc_auc_score(y[iv],m.predict_proba(Xp[iv])[:,1]))
            key=(np.mean(vals),-C,-abs(l1-0.5))
            if best is None or key>best[0]: best=(key,C,l1)
    _,C,l1=best
    model=LogisticRegression(C=C,l1_ratio=l1,penalty="elasticnet",solver="saga",max_iter=3000,tol=1e-4,n_jobs=1)
    model.fit(Xp[outer_train],y[outer_train],sample_weight=balanced_weights(donor[outer_train],y[outer_train]))
    enet_pred[test]=model.predict_proba(Xp[test])[:,1]
    coef_folds.append(model.coef_.ravel())
    enet_auc.append(roc_auc_score(y[test],enet_pred[test]))
coef_folds=np.vstack(coef_folds)
stable_positive=(coef_folds>0).mean(axis=0)>=0.8
sig_d=[g for g,k in zip(pool_genes,stable_positive) if k]
if not sig_d:
    sig_d=[pool_genes[i] for i in np.argsort(np.median(coef_folds,axis=0))[-20:] if np.median(coef_folds[:,i])>0]
signatures["SC_GART_SPARSE_ELASTIC_NET_PROGRAM"]=sig_d
spec_d=next(x["empirical_percentile"] for x in random_summary if x["signature"]=="SC_GART_STABLE_TOP_PROGRAM")
res,fold,pred,dif=evaluate_score(enet_pred,sig_d,"elastic_net_probability","SC_GART_SPARSE_ELASTIC_NET_PROGRAM",spec_d)
res["heldout_donor_mean_AUC"]=np.mean(enet_auc); results.append(res); cv_rows.extend(fold); score_store[(res["signature"],res["method"])]=enet_pred

# Nearest-centroid, donor-held-out and rank transformed over the same feature pool.
R=1-rank_sel[:,pidx]/N
cent_pred=np.zeros(n_cells); cent_weights=[]; cent_auc=[]
for d in np.unique(donor):
    tr=donor!=d; te=~tr
    w=balanced_weights(donor[tr],y[tr])
    pos=np.average(R[tr][y[tr]==1],axis=0,weights=w[y[tr]==1])
    neg=np.average(R[tr][y[tr]==0],axis=0,weights=w[y[tr]==0])
    cw=pos-neg; cent_weights.append(cw)
    cent_pred[te]=R[te]@cw
    cent_auc.append(roc_auc_score(y[te],cent_pred[te]))
cent_weights=np.vstack(cent_weights)
stable=(cent_weights>0).mean(axis=0)>=0.8
sig_e=[g for g,k in zip(pool_genes,stable) if k]
signatures["SC_GART_NEAREST_CENTROID_PROGRAM"]=sig_e
res,fold,pred,dif=evaluate_score(cent_pred,sig_e,"nearest_centroid_similarity","SC_GART_NEAREST_CENTROID_PROGRAM",spec_d)
res["heldout_donor_mean_AUC"]=np.mean(cent_auc); results.append(res); cv_rows.extend(fold); score_store[(res["signature"],res["method"])]=cent_pred

res_df=pd.DataFrame(results).sort_values("SC_SELECTION_SCORE",ascending=False).reset_index(drop=True)
# Apply locked tie rule.
top=res_df.iloc[0]
if len(res_df)>1 and top.SC_SELECTION_SCORE-res_df.iloc[1].SC_SELECTION_SCORE<0.01:
    tied=res_df[res_df.SC_SELECTION_SCORE>=top.SC_SELECTION_SCORE-0.01].copy()
    rank_methods={"UCell":0,"singscore":0,"AUCell":0,"ssGSEA":0,"mean_rank":0,"ensemble_rank":0,"mean_z":1,"elastic_net_probability":2,"nearest_centroid_similarity":2}
    tied["tie_method"] = tied.method.map(rank_methods)
    tied=tied.sort_values(["tie_method","n_genes","cell_cycle_adjusted_retention"],ascending=[True,True,False])
    winner=tied.iloc[0]
else: winner=top

res_df["selected"]=(res_df.signature==winner.signature)&(res_df.method==winner.method)
res_df.to_csv(SRC/"singlecell_donor_blocked_CV_results.csv",index=False)
pd.DataFrame(cv_rows).to_csv(SRC/"singlecell_LODO_fold_results.csv",index=False)

sig_rows=[]
for name,gs in signatures.items():
    for rank,g in enumerate(gs,1):
        row=deg.loc[deg.gene==g].iloc[0] if (deg.gene==g).any() else None
        sig_rows.append({"signature":name,"rank":rank,"gene":g,"log2FC":None if row is None else row.logFC,"P":None if row is None else row.raw_P,"FDR":None if row is None else row.FDR,"matched_singlecell":g in expressed})
sig_df=pd.DataFrame(sig_rows)
sig_df.to_csv(SRC/"candidate_signature_members.csv",index=False)
selected_genes=signatures[winner.signature]
sig_df[sig_df.signature==winner.signature].to_csv(SRC/"selected_signature_gene_list.csv",index=False)

selected_score=score_store[(winner.signature,winner.method)]
selected_cell=meta[["cell_id_revision","donor_id","GART_status","S_score","G2M_score","nCount_RNA","nFeature_RNA"]].copy()
selected_cell["selected_program_score"]=selected_score
selected_cell.to_csv(SRC/"selected_singlecell_program_scores.csv",index=False)

method_software={
    "UCell":"formula-compatible implementation; UCell package not installed",
    "singscore":"formula-compatible up-only centered rank implementation; singscore package not installed",
    "AUCell":"formula-compatible top-5%-rank AUC implementation; AUCell package not installed",
    "ssGSEA":"closed-form single-sample running-sum implementation, tau=0.25; GSVA package not installed",
    "mean_rank":"in-house fixed rank-percentile mean",
    "mean_z":"in-house global-gene standardized mean",
    "ensemble_rank":"mean within-cell percentile rank of UCell/singscore/AUCell/ssGSEA/mean-rank",
    "elastic_net_probability":"scikit-learn 1.9.0 LogisticRegression saga elastic-net",
    "nearest_centroid_similarity":"in-house donor-held-out rank-centroid projection",
}
res_df["software_note"]=res_df.method.map(method_software)
res_df.to_csv(SRC/"singlecell_donor_blocked_CV_results.csv",index=False)

lock_lines=[
    "# Single-cell selected model lock","",
    "This lock was written after single-cell-only model selection and before any spatial GART, condition, validation-marker or spatial P-value outcome was read.","",
    f"Selected signature: `{winner.signature}`.",
    f"Selected algorithm: `{winner.method}`.",
    f"SC_SELECTION_SCORE: {winner.SC_SELECTION_SCORE:.6f}.",
    f"Signature size: {len(selected_genes)} genes.",
    "Tie rule applied exactly as pre-specified." if len(res_df)>1 and top.SC_SELECTION_SCORE-res_df.iloc[1].SC_SELECTION_SCORE<0.01 else "No <0.01 top-score tie required adjudication.",
    "",
    "Spatial outcomes used for selection: none.",
    "",
    "Selected genes:","",
    ", ".join(selected_genes),"",
    "Algorithm implementation:","",
    method_software[winner.method],"",
    "The model identity is now frozen and cannot be changed after spatial unblinding."
]
(REPORTS/"SINGLECELL_SELECTED_MODEL_LOCK.md").write_text("\n".join(lock_lines),encoding="utf-8")

report=["# Single-cell model selection report","",
        "Selection used only the frozen SCP259 11-paired-donor single-cell dataset. Spatial outcomes were not accessed.","",
        f"Winner: {winner.signature} with {winner.method}; SC_SELECTION_SCORE={winner.SC_SELECTION_SCORE:.6f}; {len(selected_genes)} genes.","",
        "All candidate combinations, including non-selected and unavailable package-native implementations, are retained in `singlecell_donor_blocked_CV_results.csv`. Formula-compatible implementations are labelled explicitly and are not presented as package-native runs.","",
        "The composite selection score followed the pre-analysis weights. Matched random signatures were evaluated by mean-z donor-paired effects using 1,000 size- and expression-decile-matched sets per fixed signature.","",
        res_df[["signature","method","n_genes","heldout_donor_mean_AUC","donor_direction_consistency","cell_cycle_adjusted_retention","within_phase_consistency","Brier_score","feature_stability","random_signature_specificity","abs_QC_correlation","SC_SELECTION_SCORE","selected"]].to_markdown(index=False)]
(REPORTS/"SINGLECELL_MODEL_SELECTION_REPORT.md").write_text("\n".join(report),encoding="utf-8")

summary={"selected_signature":winner.signature,"selected_method":winner.method,"n_genes":len(selected_genes),"selection_score":float(winner.SC_SELECTION_SCORE)}
(SRC/"singlecell_selection_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
