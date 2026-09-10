# Public data dictionary

Biological replicate fields are `donor_id`/`Donor_ID` for transcriptomic analyses, `participant_id` for human tissue imaging, and `mouse_id` for animal experiments. Cells, spots, repeated biopsies, and fields of view are nested observations and are not independent biological replicates.

`GART-detected` denotes corrected raw GART count >=1; `GART-undetected` denotes raw count 0. The spatial external program is `SC_GART_NEAREST_CENTROID_PROGRAM`. The spatial ROI is the within-section top 40% EPCAM/KRT8/KRT18/KRT19 epithelial score and does not use GART, disease, the external program, or validation markers.

Figure-specific column definitions are supplied in the source-data CSV headers and the final supplementary workbooks. Human microscopy source data are de-identified participant-level numerical measurements; cross-panel participant identity is not inferred. Animal imaging source values are mouse-level summaries; fields of view are not treated as replicates.

In Supplementary Table S2, `Full_Atlas_Sample_Donor_Map` contains all 145 atlas sample entries, while `ISC_like_Analysis_Samples` contains the 79 entries contributing ISC-like cells to downstream analyses. `Preprocessing_Provenance` explicitly marks unrecoverable historical upstream parameters as unavailable. These missing parameters were not retrospectively inferred.
