from __future__ import annotations
import itertools
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data_public" / "figure_source_data" / "Figure6"
df = pd.read_csv(SRC / "GSE189184_primary_spot_source_data.csv")
mask = df["inside_locked_primary_ROI"].astype(str).str.upper().eq("TRUE")
rows = []
for donor, d in df.loc[mask].groupby("donor_id", sort=False):
    neg = d.loc[~d.GART_detected.astype(str).str.upper().eq("TRUE"), "external_program_score"]
    pos = d.loc[d.GART_detected.astype(str).str.upper().eq("TRUE"), "external_program_score"]
    if len(neg) and len(pos):
        rows.append({"donor_id": donor, "section_id": d.section_id.iloc[0], "condition": d.condition.iloc[0],
                     "GART_undetected_mean": neg.mean(), "GART_detected_mean": pos.mean(),
                     "difference": pos.mean()-neg.mean(), "n_undetected": len(neg), "n_detected": len(pos)})
don = pd.DataFrame(rows)
assert len(don) == 7
d = don.difference.to_numpy(float)
obs = d.mean()
signs = np.array(list(itertools.product([-1, 1], repeat=len(d))))
null = (signs*d).mean(axis=1)
p_exact = np.mean(np.abs(null) >= abs(obs)-1e-12)
rng = np.random.default_rng(20260909)
boot = d[rng.integers(0, len(d), (2000, len(d)))].mean(axis=1)
ci = np.quantile(boot, [.025, .975])
print({"n_donors": 7, "effect": obs, "CI_lower": ci[0], "CI_upper": ci[1],
       "exact_signflip_P": p_exact, "positive_donors": int((d > 0).sum())})
