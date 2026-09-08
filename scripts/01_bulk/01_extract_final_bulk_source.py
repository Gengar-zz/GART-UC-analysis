# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 20_BULK_EXTERNAL_COHORT_SCREEN/scripts/extract_archived_fig1a_xlsx.py
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from path_config import get_paths, require_file
PATHS = get_paths()
import csv, re, zipfile
import xml.etree.ElementTree as ET

SRC = require_file(PATHS.external_data_root / "bulk" / "legacy_three_cohort" / "Fig1A_bulk_source_data.xlsx", "archived three-cohort Figure 1A workbook")
OUT = PATHS.work_root / "bulk" / "intermediate"
OUT.mkdir(exist_ok=True)
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "pr": "http://schemas.openxmlformats.org/package/2006/relationships"}


def col_num(ref):
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n - 1


with zipfile.ZipFile(SRC) as z:
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    relmap = {r.attrib["Id"]: r.attrib["Target"] for r in rels}
    sheets = []
    for s in wb.find("m:sheets", NS):
        target = relmap[s.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        sheets.append((s.attrib["name"], target))
    for name, path in sheets[:3]:
        root = ET.fromstring(z.read(path))
        rows = []
        for row in root.findall(".//m:sheetData/m:row", NS):
            vals = {}
            for c in row.findall("m:c", NS):
                idx = col_num(c.attrib["r"])
                typ = c.attrib.get("t")
                v = c.find("m:v", NS)
                if typ == "inlineStr":
                    node = c.find("m:is", NS)
                    val = "".join(t.text or "" for t in node.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")) if node is not None else ""
                elif v is None:
                    val = ""
                elif typ == "s":
                    val = shared[int(v.text)]
                else:
                    val = v.text
                vals[idx] = val
            width = max(vals.keys(), default=-1) + 1
            rows.append([vals.get(i, "") for i in range(width)])
        maxw = max(map(len, rows))
        rows = [r + [""] * (maxw - len(r)) for r in rows]
        safe = name.replace(" ", "_")
        with open(OUT / f"archived_{safe}.csv", "w", newline="", encoding="utf-8-sig") as f:
            csv.writer(f).writerows(rows)
        print(name, len(rows)-1, OUT / f"archived_{safe}.csv")
