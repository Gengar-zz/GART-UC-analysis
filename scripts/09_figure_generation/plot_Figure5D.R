# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 23_FIG5_NATURE_STYLE_REBUILD/CODE/rebuild_fig5d_after_donor_matrix_audit.R
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
input <- gart_require_file(file.path(paths$data_public_root, "figure_source_data", "Figure3", "SCP259_primary_pseudobulk_DEG.csv"), "committed donor-paired DEG table")
out_tiff <- file.path(root, "FIGURES", "Fig5D_MYBL2_expression_rebuild.tiff")
out_pdf <- file.path(root, "FIGURES", "Fig5D_MYBL2_expression_rebuild.pdf")
out_csv <- file.path(root, "SOURCE_DATA", "Fig5D_MYBL2_expression_rebuild_source_data.csv")

frozen <- fread(input)
d <- frozen[gene == "MYBL2" & cell_threshold == 20]
stopifnot(nrow(d) == 1L)
stopifnot(abs(d$logFC - 0.526792704014918) < 1e-12)
stopifnot(abs(d$CI_low - 0.302483903129004) < 1e-12)
stopifnot(abs(d$CI_high - 0.751101504900831) < 1e-12)
stopifnot(abs(d$FDR - 0.000239919410070971) < 1e-15)
stopifnot(d$n_paired_donors == 11L)

source_data <- d[, .(
  gene,
  log2FC = logFC,
  CI_low,
  CI_high,
  nominal_P = raw_P,
  BH_FDR = FDR,
  n_paired_donors,
  cell_threshold,
  biological_unit = "donor-by-GART-state pseudobulk",
  model,
  donor_level_expression_matrix_available = FALSE,
  donor_values_reconstructed = FALSE,
  display_type = "single-row lollipop/effect plot"
)]
fwrite(source_data, out_csv)

wine <- "#7A1F35"
wine_dark <- "#551124"
rose <- "#E8C9D1"
ink <- "#202124"
gray <- "#777777"

p <- ggplot(d, aes(x = logFC, y = 1)) +
  annotate("segment", x = 0, xend = d$logFC, y = 1, yend = 1,
           linewidth = 1.55, colour = "#B98694", lineend = "round") +
  annotate("segment", x = d$CI_low, xend = d$CI_high, y = 1, yend = 1,
           linewidth = 7.5, colour = rose, lineend = "round") +
  annotate("segment", x = d$CI_low, xend = d$CI_high, y = 1, yend = 1,
           linewidth = 1.0, colour = wine_dark, lineend = "round") +
  annotate("point", x = d$logFC, y = 1, shape = 21, size = 5.4,
           stroke = 0.9, fill = wine, colour = wine_dark) +
  annotate("text", x = d$logFC, y = 1.19, label = "0.527",
           family = "Arial", fontface = "bold", size = 3.8, colour = wine_dark) +
  annotate("text", x = d$CI_low, y = 0.83, label = "0.302",
           family = "Arial", size = 3.0, colour = gray, hjust = 0.5) +
  annotate("text", x = d$CI_high, y = 0.83, label = "0.751",
           family = "Arial", size = 3.0, colour = gray, hjust = 0.5) +
  annotate("text", x = 0.90, y = 1.22,
           label = "n = 11 paired donors\nlog2FC = 0.527\n95% CI = 0.302-0.751\nBH FDR = 2.40e-4",
           family = "Arial", size = 3.35, colour = ink, hjust = 1, vjust = 1,
           lineheight = 1.12) +
  scale_x_continuous(limits = c(0, 0.94), breaks = c(0, 0.25, 0.50, 0.75),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0.72, 1.31), breaks = 1, labels = "MYBL2",
                     expand = c(0, 0)) +
  labs(
    title = "MYBL2 expression associated with the GART-detected state",
    subtitle = "Frozen donor-paired pseudobulk effect",
    x = "log2 fold change (GART-detected - GART-undetected)",
    y = NULL
  ) +
  theme_classic(base_family = "Arial", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12.2, colour = ink),
    plot.subtitle = element_text(size = 9.8, colour = gray, margin = margin(b = 8)),
    axis.title.x = element_text(size = 10.2, colour = ink, margin = margin(t = 6)),
    axis.text.x = element_text(size = 9.2, colour = ink),
    axis.text.y = element_text(size = 11.2, face = "bold", colour = ink),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.5, colour = ink),
    axis.ticks.x = element_line(linewidth = 0.45, colour = ink),
    plot.margin = margin(12, 14, 11, 12)
  )

ggsave(out_tiff, p, width = 6.8, height = 3.45, units = "in", dpi = 300,
       device = "tiff", compression = "lzw", bg = "white")
ggsave(out_pdf, p, width = 6.8, height = 3.45, units = "in",
       device = cairo_pdf, bg = "white")

message("Fig5D effect-only rebuild completed without donor-value reconstruction.")
