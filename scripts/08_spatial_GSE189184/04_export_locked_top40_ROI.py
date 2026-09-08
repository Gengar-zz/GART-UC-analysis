# Public-release copy: file paths were converted to repository-relative paths.
# Statistical logic, parameters, thresholds, and seeds are unchanged from the archived final script.
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; SRC=ROOT/"source_data"
df=pd.read_csv(SRC/"spatial_blinded_technical_data.csv")
sel=json.loads((SRC/"spatial_selection_summary.json").read_text())
df["epithelial_score"]=df[["EPCAM","KRT8","KRT18","KRT19"]].mean(axis=1)
mask=np.zeros(len(df),bool)
if sel["ROI"].startswith("epithelial_top"):
    pct=int(sel["ROI"].replace("epithelial_top",""))/100
    for sec in df.blind_section_id.unique():
        idx=np.flatnonzero(df.blind_section_id.to_numpy()==sec)
        cut=np.quantile(df.epithelial_score.to_numpy()[idx],1-pct,method="linear")
        mask[idx]=df.epithelial_score.to_numpy()[idx]>=cut
else:
    raise RuntimeError("Selected ROI membership exporter is implemented for the locked objective percentile ROI.")
out=df[["spot_id","blind_section_id"]].copy(); out["inside_locked_primary_ROI"]=mask
out.to_csv(SRC/"spatial_locked_primary_ROI_membership.csv",index=False)
assert int(mask.sum())==int(sel["n_spots"])
print(mask.sum())
