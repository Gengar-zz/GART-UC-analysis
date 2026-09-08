# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 16_FIG3_ORIGINAL_STYLE_REBUILD/CODE/plot_fig3g_representative_gene_effects.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

input_file <- gart_require_file(file.path(paths$data_public_root, "figure_source_data", "Figure3", "SCP259_primary_pseudobulk_DEG.csv"), "committed donor-paired DEG table")
out_root <- file.path(paths$results_root, "figures", "Figure3")
fig_dir <- file.path(out_root, "FIGURES")
src_dir <- file.path(out_root, "SOURCE_DATA")
report_dir <- file.path(out_root, "REPORTS")
qa_dir <- file.path(out_root, "QA")
code_dir <- file.path(out_root, "CODE")
invisible(lapply(c(fig_dir, src_dir, report_dir, qa_dir, code_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

if (.Platform$OS.type == "windows") suppressWarnings(windowsFonts(Arial = windowsFont("Arial")))

gene_pool <- c("PCNA", "MCM6", "MCM7", "FEN1", "TK1", "TYMS", "GINS2", "PAICS", "MKI67", "CDK4", "IMPDH2")
deg <- read.csv(input_file, check.names = FALSE)
plot_data <- deg[match(gene_pool, deg$gene), ]
stopifnot(nrow(plot_data) == length(gene_pool), all(plot_data$gene == gene_pool))
stopifnot(all(plot_data$n_paired_donors == 11L))

# The display order is fixed before drawing as descending archived log2FC.
plot_data <- plot_data[order(-plot_data$logFC, match(plot_data$gene, gene_pool)), ]
plot_data$display_order <- seq_len(nrow(plot_data))
plot_data$neg_log10_FDR <- -log10(plot_data$FDR)
plot_data$FDR_status <- ifelse(plot_data$FDR < 0.05, "FDR <0.05", "ns")
plot_data$analysis_type <- "Frozen SCP259 donor-paired edgeR quasi-likelihood pseudobulk DEG"
plot_data$effect_definition <- "log2 fold change: GART-detected minus GART-undetected"
plot_data$selection_display_rule <- paste0(
  "A priori proliferation/nucleotide-metabolism gene pool; all available genes retained, including IMPDH2 as nonsignificant; ",
  "display order fixed by descending archived log2FC; no donor-level matrix reconstructed"
)
plot_data$FDR_display <- ifelse(plot_data$FDR < 0.001,
                               formatC(plot_data$FDR, format = "e", digits = 1),
                               sprintf("%.3f", plot_data$FDR))
plot_data$point_label <- ifelse(plot_data$FDR_status == "ns",
                                paste0("FDR ", plot_data$FDR_display, " (ns)"),
                                paste0("FDR ", plot_data$FDR_display))
plot_data$gene_factor <- factor(plot_data$gene, levels = plot_data$gene)

source_out <- plot_data[c(
  "gene", "logFC", "SE", "CI_low", "CI_high", "raw_P", "FDR",
  "n_paired_donors", "cell_threshold", "model", "neg_log10_FDR",
  "FDR_status", "display_order", "analysis_type", "effect_definition",
  "selection_display_rule"
)]
names(source_out)[names(source_out) == "logFC"] <- "effect_log2FC"
names(source_out)[names(source_out) == "raw_P"] <- "P"
write.csv(source_out,
          file.path(src_dir, "Fig3G_representative_proliferation_nucleotide_gene_effects_source_data.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

p <- ggplot(plot_data, aes(gene_factor, logFC)) +
  geom_segment(aes(xend = gene_factor, y = 0, yend = logFC),
               linewidth = 0.8, color = "#D7A0B0") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.16,
                linewidth = 0.55, color = "#454545") +
  geom_point(aes(size = neg_log10_FDR, fill = FDR_status), shape = 21,
             color = "#8E1B3E", stroke = 0.85) +
  geom_text(aes(y = CI_high + 0.045, label = point_label),
            family = "Arial", size = 2.55, angle = 90, hjust = 0,
            color = "#303030") +
  annotate("label", x = Inf, y = Inf, label = "11 paired donors",
           hjust = 1.05, vjust = 1.25, family = "Arial", size = 3.0,
           fill = "white") +
  scale_fill_manual(values = c("FDR <0.05" = "#951D45", "ns" = "white"),
                    breaks = c("FDR <0.05", "ns"), name = NULL) +
  scale_size_continuous(range = c(3.4, 8.2), name = "−log10(FDR)") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.29))) +
  labs(
    title = "Representative proliferation and nucleotide-metabolism genes\nassociated with GART detection",
    x = "Genes",
    y = "log2 fold change\n(GART-detected − GART-undetected)"
  ) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, lineheight = 1.05, margin = margin(b = 10)),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8.5, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 8.5, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 8.5),
    legend.title = element_text(size = 8.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 12, 9, 9)
  )

stem <- "Fig3G_representative_proliferation_nucleotide_gene_effects"
tiff(file.path(fig_dir, paste0(stem, ".tiff")), width = 8.5, height = 5.8,
     units = "in", res = 300, compression = "lzw", type = "windows", bg = "white")
print(p)
dev.off()
cairo_pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = 8.5, height = 5.8,
          family = "Arial", bg = "white", onefile = TRUE)
print(p)
dev.off()

legend_text <- paste0(
  "# Figure 3G legend core text\n\n",
  "Representative proliferation and nucleotide-metabolism genes associated with GART detection in the frozen SCP259 donor-paired pseudobulk analysis. ",
  "Points show archived log2 fold changes for GART-detected minus GART-undetected ISC-like epithelial cells across 11 paired donors; vertical intervals show the archived 95% confidence intervals, and point size represents −log10(FDR). ",
  "Genes were selected a priori from proliferation and nucleotide metabolism programs identified in donor-paired pseudobulk analysis. Exact donor-level matrices were not reconstructed. ",
  "IMPDH2 was retained as a prespecified nucleotide-metabolism candidate and is explicitly marked nonsignificant. These associations identify a transcriptional state and do not establish GART-dependent causality.\n"
)
writeLines(legend_text,
           file.path(report_dir, "FIG3G_REPRESENTATIVE_GENE_EFFECTS_LEGEND.md"),
           useBytes = TRUE)

cat(paste(plot_data$gene, collapse = ", "), "\n")
