# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 18_FIG4E_CANDIDATE_REBUILD/CODE/plot_fig4e_candidate_rebuild.py
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from path_config import get_paths, require_file
PATHS = get_paths()
import csv
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor, white
from PIL import Image
import pypdfium2 as pdfium

ROOT = PATHS.results_root / "figures" / "Figure4"
FIG = ROOT
SRC = PATHS.data_public_root / "figure_source_data" / "Figure4"
SUPP_SRC = PATHS.external_data_root / "derived_inputs" / "figure4"
FIG.mkdir(parents=True, exist_ok=True)
(ROOT / "QA").mkdir(parents=True, exist_ok=True)

ARIAL = require_file(PATHS.external_data_root / "fonts" / "arial.ttf", "Arial TrueType font")
ARIAL_BOLD = require_file(PATHS.external_data_root / "fonts" / "arialbd.ttf", "Arial Bold TrueType font")
pdfmetrics.registerFont(TTFont("Arial", str(ARIAL)))
pdfmetrics.registerFont(TTFont("Arial-Bold", str(ARIAL_BOLD)))

BURGUNDY = HexColor("#8F1633"); PINK = HexColor("#D27A8C")
GRAY = HexColor("#858585"); INK = HexColor("#222222")
GRID = HexColor("#E7E7E7"); STEM = HexColor("#B7B7B7")

def read_csv(name):
    base = SUPP_SRC if name.startswith("FigS") else SRC
    with require_file(base / name, f"Figure 4 source table {name}").open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def category_color(cat):
    if cat == "Purine/nucleotide metabolism": return BURGUNDY
    if cat == "Replication/cell cycle": return PINK
    return GRAY

def draw_text(c, x, y, text, size=9, bold=False, color=INK, anchor="left"):
    font = "Arial-Bold" if bold else "Arial"
    c.setFont(font, size); c.setFillColor(color)
    if anchor == "right": x -= pdfmetrics.stringWidth(text, font, size)
    elif anchor == "center": x -= pdfmetrics.stringWidth(text, font, size) / 2
    c.drawString(x, y, text)

def pdf_to_tiff(pdf_path, tiff_path):
    doc = pdfium.PdfDocument(str(pdf_path)); page = doc[0]
    bitmap = page.render(scale=300/72)
    image = bitmap.to_pil().convert("RGB")
    image.save(tiff_path, format="TIFF", compression="tiff_lzw", dpi=(300,300))
    preview = image.copy(); preview.thumbnail((1800, 1800))
    preview.save(ROOT / "QA" / f"{pdf_path.stem}_preview.png", format="PNG")
    page.close(); doc.close()

def main_plot():
    rows = sorted(read_csv("Fig4E_focused_candidate_source_data.csv"), key=lambda r: float(r["median_rank"]))
    W, H = 9*72, 5.25*72
    stem = "Fig4E_recurrently_prioritized_candidate_genes_final"
    pdf = FIG / f"{stem}.pdf"; tif = FIG / f"{stem}.tiff"
    c = canvas.Canvas(str(pdf), pagesize=(W,H), pageCompression=1)
    c.setFillColor(white); c.rect(0,0,W,H,fill=1,stroke=0)
    left, plot_right, bottom, top = 103, 397, 91, 292
    category_x, recur_x = 485, 420
    draw_text(c, left, H-49, "Recurrently prioritized candidate genes", 12, True)
    draw_text(c, left, H-65, "Exploratory five-run computational perturbation", 8.3, color=HexColor("#555555"))
    draw_text(c, recur_x, top+18, "Recurrence", 8, True)
    draw_text(c, category_x, top+18, "Category", 8, True)
    def px(rank): return left + (64-rank)/64*(plot_right-left)
    n = len(rows); gap = (top-bottom)/(n-1)
    for tick in range(0,61,10):
        x=px(tick); c.setStrokeColor(GRID); c.setLineWidth(0.65); c.line(x,bottom-6,x,top+4)
        draw_text(c,x,bottom-21,str(tick),8.2,anchor="center")
    c.setStrokeColor(INK); c.setLineWidth(0.7); c.line(left,bottom-6,plot_right,bottom-6)
    for i,row in enumerate(rows):
        y=top-i*gap; rank=float(row["median_rank"]); color=category_color(row["functional_category"])
        draw_text(c,left-12,y-3,row["gene"],9.4,True,anchor="right")
        c.setStrokeColor(STEM); c.setLineWidth(2); c.line(left,y,px(rank),y)
        c.setFillColor(color); c.setStrokeColor(INK); c.setLineWidth(0.55); c.circle(px(rank),y,5.1,fill=1,stroke=1)
        draw_text(c,recur_x,y-3,f'{row["n_top100"]}/{row["n_runs"]}',8.7,True)
        draw_text(c,category_x,y-3,row["functional_category"],7.8,color=color)
    label="Median rank across five runs (lower rank indicates stronger prioritization)"
    draw_text(c,(left+plot_right)/2,42,label,8.6,anchor="center")
    c.showPage(); c.save(); pdf_to_tiff(pdf,tif)

def supplementary_plot():
    rows=sorted(read_csv("FigS4E_full_recurrent_gene_source_data.csv"),key=lambda r: float(r["median_rank"]))
    W,H=10*72,7.2*72; stem="FigS4E_full_robust_recurrent_gene_set"
    pdf=FIG/f"{stem}.pdf"; tif=FIG/f"{stem}.tiff"
    c=canvas.Canvas(str(pdf),pagesize=(W,H),pageCompression=1)
    c.setFillColor(white); c.rect(0,0,W,H,fill=1,stroke=0)
    left,plot_right,bottom,top=95,350,77,414
    colx=[385,450,500,550]
    draw_text(c,left,H-48,"Full robust recurrent non-GART gene set",12,True)
    draw_text(c,left,H-65,"Exploratory five-run computational perturbation",8.3,color=HexColor("#555555"))
    for x,h in zip(colx,["Top-100","Median","IQR","Category"]): draw_text(c,x,top+19,h,7.7,True)
    c.setStrokeColor(HexColor("#B5B5B5")); c.setLineWidth(0.7); c.line(colx[0],top+10,W-18,top+10)
    def px(rank): return left+(64-rank)/64*(plot_right-left)
    gap=(top-bottom)/(len(rows)-1)
    for tick in range(0,61,10):
        x=px(tick); c.setStrokeColor(GRID); c.setLineWidth(0.6); c.line(x,bottom-4,x,top+3)
        draw_text(c,x,bottom-18,str(tick),7.8,anchor="center")
    c.setStrokeColor(INK); c.line(left,bottom-4,plot_right,bottom-4)
    for i,row in enumerate(rows):
        y=top-i*gap; rank=float(row["median_rank"]); color=category_color(row["functional_category"])
        draw_text(c,left-10,y-2.8,row["gene"],8.5,True,anchor="right")
        c.setStrokeColor(STEM); c.setLineWidth(1.6); c.line(left,y,px(rank),y)
        c.setFillColor(color); c.setStrokeColor(INK); c.setLineWidth(0.45); c.circle(px(rank),y,4.2,fill=1,stroke=1)
        cat=row["functional_category"]
        short={"Purine/nucleotide metabolism":"Purine/nucleotide","Replication/cell cycle":"Replication/cell cycle"}.get(cat,"Other recurrent")
        vals=[f'{row["n_top100"]}/{row["n_runs"]}',f'{float(row["median_rank"]):g}',f'{float(row["rank_IQR"]):g}',short]
        for x,val in zip(colx,vals): draw_text(c,x,y-2.7,val,7.2,color=color if x==colx[-1] else INK)
    draw_text(c,(left+plot_right)/2,34,"Median rank across five runs (lower rank indicates stronger prioritization)",8.2,anchor="center")
    c.showPage(); c.save(); pdf_to_tiff(pdf,tif)

main_plot(); supplementary_plot()
