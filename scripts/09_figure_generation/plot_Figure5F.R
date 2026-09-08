# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 23_FIG5_NATURE_STYLE_REBUILD/CODE/replace_fig5f_marker_summary.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

root <- file.path(paths$results_root, "figures", "Figure5")
fig_dir <- root
src_dir <- file.path(root, "SOURCE_DATA")

modules <- data.table(
  module_number = 1:3,
  module = c(
    "Inflammatory context",
    "Proliferative state",
    "Nucleotide/replication program"
  ),
  display_module = c(
    "Inflammatory context",
    "Proliferative state",
    "Nucleotide/replication\nprogram"
  ),
  markers = c("MIF; CD74", "MYBL2; Ki67; PCNA", "PAICS; MCM7; GINS2"),
  x = c(0.55, 1.75, 2.95),
  xmin = c(0.08, 1.28, 2.48),
  xmax = c(1.02, 2.22, 3.42),
  fill = c("#DDE9F0", "#F3D8DF", "#D8E9E4")
)

marker_rows <- rbindlist(lapply(seq_len(nrow(modules)), function(i) {
  genes <- strsplit(modules$markers[i], "; ", fixed = TRUE)[[1]]
  data.table(
    module_number = modules$module_number[i],
    module = modules$module[i],
    marker = genes,
    display_order = seq_along(genes),
    interpretation = "candidate tissue-validation marker",
    directionality = "none",
    causal_interpretation = FALSE
  )
}))

fwrite(marker_rows, file.path(src_dir, "Fig5F_tissue_validation_marker_summary_source_data.csv"))
fwrite(marker_rows, file.path(src_dir, "Fig5F_tissue_validation_map_nodes.csv"))
fwrite(data.table(
  from = character(), to = character(), connection_status = character(), note = character()
), file.path(src_dir, "Fig5F_tissue_validation_map_edges.csv"))

ink <- "#202124"
gray <- "#717171"
border <- "#3C4145"

p <- ggplot() +
  geom_rect(
    data = modules,
    aes(xmin = xmin, xmax = xmax, ymin = 0.32, ymax = 1.38, fill = fill),
    linewidth = 0.65, colour = border
  ) +
  geom_text(
    data = modules,
    aes(x = x, y = 1.13, label = display_module),
    family = "Arial", fontface = "bold", size = 3.35, colour = ink,
    lineheight = 0.97
  ) +
  geom_text(
    data = modules,
    aes(x = x, y = 0.73, label = gsub("; ", "\n", markers, fixed = TRUE)),
    family = "Arial", fontface = "bold", size = 4.0, colour = ink,
    lineheight = 1.22
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 3.50), ylim = c(0.18, 1.55), clip = "off") +
  labs(
    title = "Candidate tissue-validation markers associated with\nthe GART-detected ISC-like state",
    subtitle = "Computational candidates prioritized for immunofluorescence evaluation"
  ) +
  theme_void(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11.6, colour = ink,
                              hjust = 0, lineheight = 1.02),
    plot.subtitle = element_text(size = 9.4, colour = gray, hjust = 0,
                                 margin = margin(t = 3, b = 8)),
    plot.margin = margin(12, 12, 10, 12)
  )

ggsave(file.path(fig_dir, "Fig5F_tissue_validation_map.tiff"), p,
       width = 6.4, height = 3.55, units = "in", dpi = 300,
       device = "tiff", compression = "lzw", bg = "white")
ggsave(file.path(fig_dir, "Fig5F_tissue_validation_map.pdf"), p,
       width = 6.4, height = 3.55, units = "in",
       device = cairo_pdf, bg = "white")

message("Fig5F replaced with three independent marker modules; no connections drawn.")
