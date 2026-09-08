from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data_public" / "figure_source_data" / "Figure6"
OUT = ROOT / "outputs" / "Figure6"
OUT.mkdir(parents=True, exist_ok=True)
spots = pd.read_csv(SRC / "GSE189184_primary_spot_source_data.csv")
don = pd.read_csv(SRC / "revised_primary_donor_scores.csv")

def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

def maps(value, title, name, cmap="magma"):
    sections = list(dict.fromkeys(spots.section_id))
    fig, axes = plt.subplots(2, 4, figsize=(10, 5.2), constrained_layout=True)
    norm = plt.Normalize(np.nanquantile(spots[value], .01), np.nanquantile(spots[value], .99))
    sc = None
    for ax, sec in zip(axes.flat, sections):
        d = spots[spots.section_id == sec]
        sc = ax.scatter(d.pixel_col, d.pixel_row, c=d[value], s=3, cmap=cmap, norm=norm, linewidths=0)
        roi = d.inside_locked_primary_ROI.astype(str).str.upper().eq("TRUE")
        ax.scatter(d.loc[roi, "pixel_col"], d.loc[roi, "pixel_row"], s=5, facecolors="none", edgecolors="#E6C15A", linewidths=.2)
        ax.set_title(f"{sec} ({d.condition.iloc[0]})", fontsize=8); ax.invert_yaxis(); ax.axis("off")
    axes.flat[-1].axis("off")
    fig.colorbar(sc, ax=list(axes.flat[:-1]), fraction=.02, pad=.01, label=value.replace("_", " "))
    fig.suptitle(title, fontweight="bold")
    save(fig, name)

maps("GART_expression", "Spatial distribution of GART", "Fig6A_GART_spatial_distribution")
maps("external_program_score", "Single-cell-derived GART-associated epithelial program", "Fig6B_external_program_locked_ROI", "coolwarm")
for gene in ["MYBL2", "PCNA", "PAICS", "MCM7", "GINS2", "MIF", "CD74"]:
    maps(gene, f"Spatial distribution of {gene}", f"Fig6_{gene}_descriptive_map")

fig, ax = plt.subplots(figsize=(4.4, 4.2), constrained_layout=True)
for _, r in don.iterrows():
    ax.plot([0,1], [r.GART_undetected_mean, r.GART_detected_mean], color="#B9BDC5", lw=1)
    ax.scatter(0, r.GART_undetected_mean, color="#C9CCD2", s=34, zorder=2)
    ax.scatter(1, r.GART_detected_mean, color="#8B1E3F", s=34, zorder=2)
ax.set_xticks([0,1], ["GART-undetected", "GART-detected"]); ax.set_ylabel("External program score")
ax.set_title("Donor-level external-program difference", fontweight="bold")
ax.spines[["top","right"]].set_visible(False)
save(fig, "Fig6C_donor_level_primary")
