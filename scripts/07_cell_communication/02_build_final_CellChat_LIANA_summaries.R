# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 23_FIG5_NATURE_STYLE_REBUILD/CODE/build_fig5_nature_style.R
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
  library(cowplot)
  library(scales)
})

if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

input_root <- file.path(paths$work_root, "cell_communication", "figure5_positive_main")
out_root <- file.path(paths$results_root, "figures", "Figure5")
fig_dir <- file.path(out_root, "FIGURES")
src_dir <- file.path(out_root, "SOURCE_DATA")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

wine <- "#7A1F35"
wine_dark <- "#5B1026"
rose <- "#D9A4AE"
blue <- "#5B86A6"
gold <- "#C99A4A"
teal <- "#438B83"
ink <- "#202124"
gray <- "#777777"
light_gray <- "#E8E8E8"
sender_order <- c("Myeloid", "Fibroblast", "Endothelial", "T/NK")

theme_cell <- function(base_size = 10) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.2, colour = ink, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = gray, margin = margin(b = 6)),
      axis.title = element_text(size = base_size, colour = ink),
      axis.text = element_text(size = base_size - 1, colour = ink),
      axis.line = element_line(linewidth = 0.45, colour = ink),
      axis.ticks = element_line(linewidth = 0.4, colour = ink),
      legend.title = element_text(size = base_size - 0.5, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(9, 10, 8, 9)
    )
}

save_plot_pair <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".tiff")), plot = plot,
         width = width, height = height, units = "in", dpi = 300,
         device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot = plot,
         width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
}

fmt_fdr <- function(x) {
  ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2),
         formatC(x, format = "f", digits = 3))
}

# Frozen inputs only; no statistical model is fit in this script.
cc_all <- fread(file.path(input_root, "SOURCE_DATA", "Fig5A_CellChat_source_data.csv"))
li_all <- fread(file.path(input_root, "SOURCE_DATA", "Fig5B_LIANA_source_data.csv"))
mybl2 <- fread(file.path(input_root, "SOURCE_DATA", "Fig5D_MYBL2_paired_source_data.csv"))
genes <- fread(file.path(input_root, "SOURCE_DATA", "Fig5E_positive_genes_source_data.csv"))
input_files <- file.path(input_root, "SOURCE_DATA", c("Fig5A_CellChat_source_data.csv", "Fig5B_LIANA_source_data.csv", "Fig5D_MYBL2_paired_source_data.csv", "Fig5E_positive_genes_source_data.csv"))
invisible(lapply(input_files, gart_require_file, label = "upstream Figure 5 source table"))


# A. Donor-aware CellChat display.
cc <- cc_all[source %chin% sender_order]
cc[, source := factor(source, levels = sender_order)]
fwrite(cc, file.path(src_dir, "Fig5A_MIF_CellChat_donor_source_data.csv"))

pA <- ggplot(cc, aes(source, difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#9A9A9A") +
  geom_boxplot(width = 0.48, outlier.shape = NA, fill = alpha(rose, 0.35),
               colour = wine_dark, linewidth = 0.55) +
  geom_jitter(width = 0.13, height = 0, size = 2.15, alpha = 0.86,
              colour = wine, shape = 16) +
  annotate("text", x = 4.42, y = max(cc$difference, na.rm = TRUE) * 0.96,
           label = "P = 0.001953\nBH FDR = 0.008275\nn = 11 donors",
           hjust = 1, vjust = 1, family = "Arial", size = 3.05, colour = ink) +
  scale_y_continuous(labels = label_number(accuracy = 0.002), expand = expansion(mult = c(0.08, 0.18))) +
  labs(title = "MIF/CD74 communication",
       subtitle = "CellChat donor-paired differences",
       x = "Sender compartment",
       y = "MIF/CD74 communication score difference\n(GART-detected - GART-undetected)") +
  theme_cell(10) +
  theme(axis.text.x = element_text(angle = 24, hjust = 1), plot.margin = margin(10, 10, 8, 10))

save_plot_pair(pA, "Fig5A_MIF_CellChat_donor", 5.9, 4.45)

# B. LIANA+ donor display plus compact five-method consistency dots.
li <- li_all[method == "LIANA consensus" & source %chin% sender_order]
li[, source := factor(source, levels = sender_order)]
fwrite(li, file.path(src_dir, "Fig5B_LIANA_consensus_donor_source_data.csv"))

pBmain <- ggplot(li, aes(source, difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#9A9A9A") +
  geom_boxplot(width = 0.48, outlier.shape = NA, fill = alpha(blue, 0.24),
               colour = "#315E7C", linewidth = 0.55) +
  geom_jitter(width = 0.13, height = 0, size = 2.15, alpha = 0.86,
              colour = blue, shape = 16) +
  annotate("text", x = 4.42, y = max(li$difference, na.rm = TRUE) * 0.96,
           label = "P = 0.0009766\nBH FDR = 0.01020\nn = 11 donors",
           hjust = 1, vjust = 1, family = "Arial", size = 3.05, colour = ink) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.08, 0.18))) +
  labs(title = "LIANA+ consensus communication",
       subtitle = "Donor-paired MIF/CD74-family differences",
       x = "Sender compartment",
       y = "MIF/CD74 communication score difference\n(GART-detected - GART-undetected)") +
  theme_cell(10) +
  theme(axis.text.x = element_text(angle = 24, hjust = 1), plot.margin = margin(10, 2, 8, 10))

method_order <- c("CellPhoneDB", "Connectome", "Geometric mean", "NATMI", "LIANA consensus")
method_summary <- unique(li_all[, .(method, paired_P, BH_FDR, positive_donor_n, paired_donor_n)])
method_summary <- method_summary[method %chin% method_order]
method_summary[, method := factor(method, levels = rev(method_order))]
method_summary[, direction := "Positive"]
fwrite(method_summary, file.path(src_dir, "Fig5B_method_consistency_source_data.csv"))

pBdot <- ggplot(method_summary, aes(x = 1, y = method)) +
  geom_vline(xintercept = 1, linewidth = 0.35, colour = light_gray) +
  geom_point(size = 3.2, shape = 21, stroke = 0.65, fill = wine, colour = wine_dark) +
  annotate("text", x = 1, y = 5.75, label = "5/5 positive", family = "Arial",
           fontface = "bold", size = 3.0, colour = ink) +
  coord_cartesian(xlim = c(0.65, 1.35), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "Arial", base_size = 9) +
  theme(axis.text.y = element_text(size = 8.2, colour = ink, hjust = 1),
        plot.margin = margin(27, 7, 18, 2))

pB <- plot_grid(pBmain, pBdot, nrow = 1, rel_widths = c(0.76, 0.24), align = "h", axis = "tb")
save_plot_pair(pB, "Fig5B_LIANA_consensus_donor", 6.4, 4.45)

# C. Compact biological-context heatmap. Raw frozen donor medians are retained;
# colour is scaled within each framework solely to avoid incomparable score scales.
cc_h <- cc[, .(raw_score = median(difference)), by = source]
cc_h[, framework := "CellChat"]
li_h <- li[, .(raw_score = median(difference)), by = source]
li_h[, framework := "LIANA+ consensus"]
heat <- rbindlist(list(cc_h, li_h), use.names = TRUE)
heat[, display_scaled := {
  rr <- range(raw_score, na.rm = TRUE)
  if (diff(rr) == 0) rep(0.5, .N) else (raw_score - rr[1]) / diff(rr)
}, by = framework]
heat[, source := factor(source, levels = rev(sender_order))]
heat[, framework := factor(framework, levels = c("CellChat", "LIANA+ consensus"))]
heat[, score_label := ifelse(framework == "CellChat",
                             formatC(raw_score, format = "f", digits = 4),
                             formatC(raw_score, format = "f", digits = 3))]
fwrite(heat, file.path(src_dir, "Fig5C_MIF_CD74_context_heatmap_source_data.csv"))

pC <- ggplot(heat, aes(framework, source, fill = display_scaled)) +
  geom_tile(colour = "white", linewidth = 1.4, width = 0.96, height = 0.94) +
  geom_text(aes(label = score_label), family = "Arial", size = 3.35,
            colour = ifelse(heat$display_scaled > 0.62, "white", ink), fontface = "bold") +
  scale_fill_gradient(low = "#F4E8EB", high = wine_dark, limits = c(0, 1),
                      name = "Relative\nintensity") +
  labs(title = "MIF/CD74 communication landscape",
       subtitle = "Donor-aware paired analyses\nDirectionally consistent candidate communication",
       x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11.2, colour = ink),
        plot.subtitle = element_text(size = 9.5, colour = gray, lineheight = 1.05, margin = margin(b = 7)),
        axis.text = element_text(size = 9.3, colour = ink),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8.8, face = "bold"),
        legend.text = element_text(size = 8.2),
        plot.margin = margin(10, 10, 8, 10))

save_plot_pair(pC, "Fig5C_MIF_CD74_context_heatmap", 5.8, 4.45)

# D. Effect-only MYBL2 panel; no donor-level values are reconstructed.
fwrite(mybl2, file.path(src_dir, "Fig5D_MYBL2_effect_source_data.csv"))
pD <- ggplot(mybl2, aes(y = "MYBL2", x = estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "#8A8A8A") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.13,
                 linewidth = 1.05, colour = wine_dark) +
  geom_point(size = 4.4, shape = 21, fill = wine, colour = wine_dark, stroke = 0.75) +
  annotate("text", x = 0.03, y = 1.28,
           label = "effect = 0.527  (95% CI 0.302-0.751)\nFDR = 2.40e-4; n = 11 donors",
           hjust = 0, family = "Arial", size = 3.15, colour = ink) +
  scale_x_continuous(limits = c(-0.06, 0.93), breaks = c(0, 0.25, 0.5, 0.75),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "MYBL2 expression effect",
       subtitle = "Frozen donor-paired pseudobulk estimate",
       x = "log2 fold change\n(GART-detected - GART-undetected)", y = NULL) +
  theme_cell(10) +
  theme(axis.text.y = element_text(face = "bold", size = 10.5),
        axis.ticks.y = element_blank(), plot.margin = margin(10, 10, 8, 10))

save_plot_pair(pD, "Fig5D_MYBL2_effect", 5.8, 4.15)

# E. Six frozen positive genes; PAICS and MCM7 receive a visual outline only.
genes[, gene := factor(gene, levels = rev(unique(gene[order(display_order)])))]
genes[, emphasis := gene %chin% c("PAICS", "MCM7")]
genes[, fdr_label := paste0("FDR ", fmt_fdr(FDR))]
genes[, point_size := pmin(6.8, pmax(3.2, -log10(FDR)))]
fwrite(genes, file.path(src_dir, "Fig5E_proliferation_nucleotide_lollipop_source_data.csv"))

pE <- ggplot(genes, aes(estimate, gene)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#808080") +
  geom_segment(aes(x = 0, xend = estimate, yend = gene), linewidth = 0.8,
               colour = alpha(wine, 0.62)) +
  geom_point(aes(size = point_size, fill = emphasis), shape = 21,
             colour = wine_dark, stroke = 0.85) +
  geom_text(aes(x = estimate + 0.045, label = fdr_label), hjust = 0, vjust = 0.45,
            family = "Arial", size = 2.65, colour = ink) +
  scale_fill_manual(values = c(`FALSE` = wine, `TRUE` = gold), guide = "none") +
  scale_size_identity() +
  scale_x_continuous(limits = c(0, 1.11), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "Proliferation and nucleotide programs",
       subtitle = "Frozen donor-paired pseudobulk effects",
       x = "log2 fold change\n(GART-detected - GART-undetected)", y = NULL) +
  theme_cell(10) +
  theme(axis.text.y = element_text(
          face = ifelse(levels(genes$gene) %in% c("PAICS", "MCM7"), "bold", "plain")),
        plot.margin = margin(10, 10, 8, 10))

save_plot_pair(pE, "Fig5E_proliferation_nucleotide_lollipop", 6.2, 4.45)

# F. Association map for tissue-validation markers; dashed connectors have no arrows.
nodes <- data.table(
  node = c("GART-associated\nISC-like state", "Inflammatory context\nMIF / CD74",
           "Proliferative program\nMYBL2 / Ki67", "Nucleotide program\nPAICS / PCNA"),
  x = c(0, -1.50, 1.48, 1.48), y = c(0, 0, 0.78, -0.78),
  group = c("center", "inflammatory", "proliferative", "nucleotide")
)
edges <- data.table(x = c(-1.02, 1.02, 1.02), y = c(0, 0.52, -0.52),
                    xend = c(-0.48, 0.48, 0.48), yend = c(0, 0.17, -0.17))
fwrite(nodes, file.path(src_dir, "Fig5F_tissue_validation_map_nodes.csv"))
fwrite(edges, file.path(src_dir, "Fig5F_tissue_validation_map_edges.csv"))

pF <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
               linetype = "22", linewidth = 0.9, colour = "#7A7A7A") +
  geom_label(data = nodes, aes(x, y, label = node, fill = group),
             family = "Arial", fontface = "bold", size = 3.3,
             label.size = 0.55, label.r = unit(0.17, "lines"),
             label.padding = unit(0.42, "lines"), colour = ink) +
  scale_fill_manual(values = c(center = "#F3D8DF", inflammatory = "#DDE9F0",
                               proliferative = "#F3E3BC", nucleotide = "#D8E9E4"), guide = "none") +
  coord_fixed(xlim = c(-2.25, 2.25), ylim = c(-1.28, 1.28), clip = "off") +
  labs(title = "Candidate tissue-validation markers",
       subtitle = "Associated features of the GART-detected state") +
  theme_void(base_family = "Arial", base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11.2, colour = ink, hjust = 0),
        plot.subtitle = element_text(size = 9.5, colour = gray, hjust = 0, margin = margin(b = 5)),
        plot.margin = margin(10, 10, 8, 10))

save_plot_pair(pF, "Fig5F_tissue_validation_map", 6.2, 4.25)

# Composite preview for layout selection; panel lettering is only added here.
row1 <- plot_grid(pA, pB, labels = c("A", "B"), label_fontfamily = "Arial",
                  label_fontface = "bold", label_size = 15, label_x = -0.01, label_y = 1.015,
                  rel_widths = c(1, 1.08))
row2 <- plot_grid(pC, pD, labels = c("C", "D"), label_fontfamily = "Arial",
                  label_fontface = "bold", label_size = 15, label_x = -0.01, label_y = 1.015,
                  rel_widths = c(1, 1))
row3 <- plot_grid(pE, pF, labels = c("E", "F"), label_fontfamily = "Arial",
                  label_fontface = "bold", label_size = 15, label_x = -0.01, label_y = 1.015,
                  rel_widths = c(1, 1))
body <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1.04, 1, 1))
title_grob <- ggdraw() +
  draw_label("Donor-aware inflammatory context and proliferative programs associated with\nGART-detected ISC-like epithelial cells",
             x = 0.015, hjust = 0, y = 0.56, fontfamily = "Arial", fontface = "bold",
             size = 15, colour = ink)
preview <- plot_grid(title_grob, body, ncol = 1, rel_heights = c(0.095, 0.905))
save_plot_pair(preview, "Figure5_AF_Nature_style_preview", 12.2, 13.2)

message("Figure 5 Nature-style visualization-only rebuild complete.")
