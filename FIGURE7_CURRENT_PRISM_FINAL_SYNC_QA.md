# Figure 7 current-Prism final synchronization QA

## Release decision

**PASS.** The three task-start Prism files were copied to a read-only staging directory before extraction. All subsequent extraction and full-precision reference calculations used those frozen copies. The public 13-decimal inputs reproduce all 13 endpoints at manuscript reporting precision.

## Required checks

1. Current Fig7A Prism snapshot SHA = `5BE4820A88753C266312D41A7B2DD3A57260101B8EF3AC8447D334F69895D9EC`
2. Current Fig7B Prism snapshot SHA = `31F0B9FCC2579FE26C470174A9CAA9363DB87DBE9623EC7C3F8CC3657F9FEE6C`
3. Current Fig7C Prism snapshot SHA = `8CD357D3733E1ACF8280C77DC004877625655DA1146C584EDA760DEBC728B982`
4. Prism snapshots frozen before analysis = **YES**
5. 13 endpoints / 752 values extracted = **YES**
6. All public values use no more than 13 decimal places = **YES**
7. Fabricated trailing precision = **NO**
8. 13/13 reproduce current full-precision Prism at manuscript reporting precision = **PASS**
9. Fig7B current exact HC cross-endpoint shared-value count = **0**
10. Whole-column duplication = **NO**
11. Duplicate observation rows = **0**
12. Impossible-range values = **0**
13. GitHub old Figure7 formal CSV replaced = **YES**
14. GitHub Fig7B Prism object replaced with current snapshot = **YES**
15. Redundant old Figure7 CSV files removed = **YES**. The four `DISPLAY_ONLY_ROUNDED` CSVs and obsolete `Fig7_provenance.csv` were removed because the public reproduction script has no dependency on them.
16. Figure7.tif modified = **NO**. Figure 7 was outside this task scope and was not modified; the file was not restored or copied when absent from the submission directory.
17. Figure legends modified = **NO**
18. Manuscript only Fig7B Results numbers changed = **YES**. Only paragraph 141 changed in each manuscript version.
19. Point-by-point only Reviewer 3 latest Comment 2 changed = **YES**. Only paragraphs 193–196 changed.
20. Table S1 changed = **NO**; SHA-256 remains `A5C9993AF48A502F518523A65B1DF9FB4023BF24420BE014EB925C88C3339642`.
21. Bulk changed = **NO**
22. Data Availability changed = **NO**
23. Other Results changed = **NO**

`FIGURE7B_FIGURE_LOCKED_BY_AUTHOR = YES`

## Fig. 7B synchronized manuscript values

- Ki67-positive cells: `58.66 (50.64–66.69)` → `58.67 (50.65–66.69)` percentage points.
- GART/Ki67 double-positive cells: `41.46 (36.63–46.30)` → `41.46 (36.63–46.29)` percentage points.
- GART/MYBL2/Ki67 triple-positive cells: `37.46 (34.76–44.73)` → `32.44 (23.19–41.74)` percentage points.
- All three remain `P < 0.0001`; direction, significance, figure labels, and biological interpretation are unchanged.
- GART and MYBL2 values match the existing manuscript at its reporting precision and were not changed.

## Document integrity

- Clean manuscript: 231 paragraphs before and after; only paragraph 141 differs.
- Changes-marked manuscript: 231 paragraphs before and after; only paragraph 141 differs. Only the six changed numerical tokens are newly red; all prior red revisions are preserved.
- Point-by-point response: 197 paragraphs before and after; only paragraphs 193–196 differ.
- Automated DOCX rendering was attempted but could not be completed because LibreOffice was unavailable and the local Word background converter was not usable. DOCX package integrity, paragraph counts, paragraph-level diffs, run colors, and final-file SHA matches were verified programmatically.

## Public reproduction

Running `scripts/11_imaging_statistics/reproduce_Figure7_statistics.R` against the distributed `Figure7_formal_statistical_inputs.csv` completed successfully with:

`Figure 7 current-Prism reproduction: 13/13 PASS`
