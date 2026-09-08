# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 20_BULK_EXTERNAL_COHORT_SCREEN/scripts/compute_four_cohort_meta.py
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from path_config import get_paths, require_file
PATHS = get_paths()
import numpy as np
import pandas as pd
from scipy import optimize, stats

ROOT = PATHS.repo_root
OUT = PATHS.results_root / "bulk" / "results"
OUT.mkdir(parents=True, exist_ok=True)
OLD = require_file(PATHS.external_data_root / "bulk" / "legacy_three_cohort" / "bulk_meta_results_revised.csv", "archived three-cohort meta-analysis result")


def reml_tau2(y, v):
    # Intercept-only REML log-likelihood, constants omitted.
    def objective(tau2):
        w = 1.0 / (v + tau2)
        mu = np.sum(w * y) / np.sum(w)
        return 0.5 * (np.sum(np.log(v + tau2)) + np.log(np.sum(w)) + np.sum(w * (y - mu) ** 2))
    upper = max(10.0, float(np.var(y, ddof=1) * 20 + np.max(v)))
    fit = optimize.minimize_scalar(objective, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-14})
    tau2 = float(fit.x)
    if objective(0.0) <= objective(tau2) + 1e-10:
        tau2 = 0.0
    return tau2


def hk_meta(frame):
    y = frame.Hedges_g.to_numpy(float)
    se = frame.SE.to_numpy(float)
    v = se ** 2
    k = len(y)
    tau2 = reml_tau2(y, v)
    w = 1.0 / (v + tau2)
    mu = float(np.sum(w * y) / np.sum(w))
    q_hk = float(np.sum(w * (y - mu) ** 2) / (k - 1))
    se_hk = float(np.sqrt(q_hk / np.sum(w)))
    crit = float(stats.t.ppf(.975, k - 1))
    ci = (mu - crit * se_hk, mu + crit * se_hk)
    p = float(2 * stats.t.sf(abs(mu / se_hk), k - 1))
    w_fe = 1.0 / v
    mu_fe = float(np.sum(w_fe * y) / np.sum(w_fe))
    Q = float(np.sum(w_fe * (y - mu_fe) ** 2))
    Qp = float(stats.chi2.sf(Q, k - 1))
    # metafor-style I2 uses tau2 relative to the typical within-study variance.
    typical_v = (k - 1) * np.sum(w_fe) / (np.sum(w_fe) ** 2 - np.sum(w_fe ** 2))
    I2 = 100 * tau2 / (tau2 + typical_v) if tau2 > 0 else 0.0
    weights = w / np.sum(w) * 100
    return {"k": k, "pooled_Hedges_g": mu, "SE_HK": se_hk, "CI_low": ci[0], "CI_high": ci[1],
            "P_value": p, "Cochran_Q": Q, "Q_P": Qp, "I2_percent": I2, "tau2_REML": tau2,
            "HK_scale": q_hk}, weights


old = pd.read_csv(OLD)
coh = old.loc[old.row_type.eq("cohort"), ["dataset", "n_Healthy", "n_UC", "Hedges_g", "SE", "CI_low", "CI_high"]].copy()
new = pd.read_csv(require_file(OUT / "GSE87466_GART_result.csv", "upstream GSE87466 result; run 02_process_GSE87466.py first"))
coh = pd.concat([coh, pd.DataFrame([{
    "dataset": "GSE87466", "n_Healthy": int(new.n_Healthy.iloc[0]), "n_UC": int(new.n_Active_UC.iloc[0]),
    "Hedges_g": new.Hedges_g.iloc[0], "SE": new.SE.iloc[0], "CI_low": new.CI_low.iloc[0], "CI_high": new.CI_high.iloc[0]
}])], ignore_index=True)

summary, weights = hk_meta(coh)
coh["cohort_weight_percent"] = weights
coh.to_csv(OUT / "four_cohort_effect_sizes.csv", index=False)
pd.DataFrame([summary]).to_csv(OUT / "four_cohort_meta_summary.csv", index=False)

loo = []
for omitted in coh.dataset:
    s, _ = hk_meta(coh.loc[~coh.dataset.eq(omitted)].reset_index(drop=True))
    s = {"omitted_dataset": omitted, **s}
    loo.append(s)
pd.DataFrame(loo).to_csv(OUT / "four_cohort_leave_one_out.csv", index=False)

# Verification: same implementation should reproduce the archived 3-cohort estimate.
check, _ = hk_meta(coh.iloc[:3].copy())
pd.DataFrame([check]).to_csv(OUT / "three_cohort_method_reproduction_check.csv", index=False)
print("Four-cohort summary")
print(pd.DataFrame([summary]).to_string(index=False))
print("\nLOO")
print(pd.DataFrame(loo)[["omitted_dataset", "pooled_Hedges_g", "CI_low", "CI_high", "P_value", "I2_percent", "tau2_REML"]].to_string(index=False))
print("\nArchived method reproduction check")
print(pd.DataFrame([check]).to_string(index=False))
