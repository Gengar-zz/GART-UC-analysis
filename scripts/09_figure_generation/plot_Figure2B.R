# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 11_FIG2B_STYLE_REBUILD/CODE/plot_fig2b_style_rebuild.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors = FALSE)
set.seed(20260905)

obj_path <- gart_require_file(file.path(paths$external_data_root, "singlecell", "frozen_objects", "UC_GART_06c_Epithelial_Fig2AB_noGARTname_FINAL.rds"), "frozen epithelial Seurat object")
out_dir <- file.path(paths$results_root, "figures", "Figure2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

state_order <- c(
  "Low-confidence / other epithelial",
  "Rare epithelial cells",
  "BEST4+ epithelial cells",
  "TA / proliferating cells",
  "Stem/ISC-like epithelial cells",
  "Regenerative epithelial cells",
  "Goblet cells",
  "Absorptive epithelial cells"
)

marker_sets <- list(
  `Stem / ISC` = c("LGR5", "ASCL2", "SMOC2", "OLFM4"),
  `TA / proliferating` = c("MKI67", "TOP2A", "PCNA", "UBE2C"),
  `Absorptive epithelial` = c("FABP1", "CA1", "KRT20", "AQP8", "CA2"),
  `Goblet / secretory` = c("MUC2", "TFF3", "FCGBP", "SPINK4", "CLCA1"),
  `BEST4+` = c("BEST4", "OTOP2", "CA7", "GUCA2B"),
  `Rare epithelial` = c("CHGA", "CHGB", "TRPM5", "POU2F3", "NEUROD1", "AVIL"),
  `Inflammatory / stress-associated` = c("LCN2", "DUOX2", "REG3A", "CCL20"),
  `Contamination check` = c("PTPRC", "COL1A1", "VWF", "RGS5")
)
marker_group_map <- setNames(rep(names(marker_sets), lengths(marker_sets)), unlist(marker_sets, use.names = FALSE))

ep <- readRDS(obj_path)
meta <- ep[[]]
stopifnot("Epi_Subtype_Fig2" %in% colnames(meta))
subtype <- as.character(meta$Epi_Subtype_Fig2)
if (anyNA(subtype) || any(!nzchar(subtype))) stop("Frozen Epi_Subtype_Fig2 contains missing labels")
if (!setequal(unique(subtype), state_order)) stop("Frozen epithelial-state labels differ from the locked eight-state set")

assay <- ep[["RNA"]]
requested_genes <- unique(unlist(marker_sets, use.names = FALSE))
data_layers <- grep("^data", Layers(assay), value = TRUE)
if (!length(data_layers)) stop("No normalized RNA data layers found")

layer_full <- function(layer_name) {
  q <- methods::slot(assay, "layers")[[layer_name]]
  feature_names <- methods::slot(assay, "features")[[layer_name]]
  cell_names <- methods::slot(assay, "cells")[[layer_name]]
  if (is.null(dim(q))) q <- matrix(q, nrow = length(feature_names), ncol = length(cell_names))
  if (is.null(rownames(q))) rownames(q) <- feature_names
  if (is.null(colnames(q))) colnames(q) <- cell_names
  q
}

available_genes <- unique(unlist(lapply(data_layers, function(layer_name) rownames(layer_full(layer_name)))))
missing_genes <- setdiff(requested_genes, available_genes)
included_genes <- requested_genes[requested_genes %in% available_genes]
if (!length(included_genes)) stop("None of the requested markers is present")

expr <- matrix(0, nrow = length(included_genes), ncol = nrow(meta),
               dimnames = list(included_genes, rownames(meta)))
for (layer_name in data_layers) {
  q <- layer_full(layer_name)
  g <- intersect(included_genes, rownames(q))
  c <- intersect(rownames(meta), colnames(q))
  if (length(g) && length(c)) expr[g, c] <- as.matrix(q[g, c, drop = FALSE])
}

dot <- do.call(rbind, lapply(included_genes, function(gene_name) {
  do.call(rbind, lapply(state_order, function(state_name) {
    idx <- which(subtype == state_name)
    data.frame(
      marker_group = unname(marker_group_map[gene_name]),
      gene = gene_name,
      epithelial_state = state_name,
      n_cells = length(idx),
      average_expression = mean(expr[gene_name, idx]),
      percent_expressed = mean(expr[gene_name, idx] > 0) * 100,
      stringsAsFactors = FALSE
    )
  }))
}))

dot$avg_exp_scaled <- ave(dot$average_expression, dot$gene, FUN = function(v) {
  z <- as.numeric(scale(v))
  z[!is.finite(z)] <- 0
  pmax(-2.5, pmin(2.5, z))
})
dot$marker_group <- factor(dot$marker_group, levels = names(marker_sets))
dot$gene <- factor(dot$gene, levels = included_genes)
dot$epithelial_state <- factor(dot$epithelial_state, levels = rev(state_order))
dot <- dot[order(dot$marker_group, dot$gene, dot$epithelial_state), ]

write.csv(dot, file.path(out_dir, "Fig2B_faceted_epithelial_marker_dotplot_source_data.csv"), row.names = FALSE)
write.csv(data.frame(
  requested_gene = requested_genes,
  included = requested_genes %in% included_genes,
  marker_group = unname(marker_group_map[requested_genes])
), file.path(out_dir, "Fig2B_marker_availability_audit.csv"), row.names = FALSE)

p <- ggplot(dot, aes(x = gene, y = epithelial_state)) +
  geom_point(aes(size = percent_expressed, colour = avg_exp_scaled), alpha = 0.95) +
  facet_grid(cols = vars(marker_group), scales = "free_x", space = "free_x",
             labeller = labeller(marker_group = label_wrap_gen(width = 19))) +
  scale_colour_gradient2(
    low = "#3B6FB6", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    limits = c(-2.5, 2.5), breaks = c(-2, 0, 2),
    name = "Average Expression\n(scaled)"
  ) +
  scale_size_continuous(
    range = c(0.4, 5.2), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
    name = "Percent expressed"
  ) +
  labs(x = NULL, y = NULL, title = "Canonical epithelial marker expression", tag = "B") +
  theme_bw(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 8)),
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0.003, 0.995),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8, colour = "black"),
    axis.text.y = element_text(face = "bold", size = 8.2, colour = "black"),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.28),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "#B8B8B8", linewidth = 0.35),
    strip.background = element_rect(fill = "#F2F2F2", colour = "#B8B8B8", linewidth = 0.4),
    strip.text.x = element_text(face = "bold", size = 8.0, colour = "black", margin = margin(4, 2, 4, 2)),
    panel.spacing.x = unit(0.10, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 8),
    legend.box.spacing = unit(0.12, "in"),
    plot.margin = margin(8, 10, 8, 8)
  ) +
  guides(
    colour = guide_colourbar(order = 1, barheight = unit(1.45, "in"), barwidth = unit(0.15, "in")),
    size = guide_legend(order = 2, override.aes = list(colour = "#7A7A7A", alpha = 1))
  )

tiff_path <- file.path(out_dir, "Fig2B_faceted_epithelial_marker_dotplot.tiff")
pdf_path <- file.path(out_dir, "Fig2B_faceted_epithelial_marker_dotplot.pdf")
grDevices::tiff(tiff_path, width = 15.6, height = 7.5, units = "in", res = 300,
                compression = "lzw", bg = "white", type = "windows", antialias = "cleartype")
print(p)
dev.off()
grDevices::cairo_pdf(pdf_path, width = 15.6, height = 7.5, family = "Arial", bg = "white", onefile = TRUE)
print(p)
dev.off()

cat("states=", length(unique(dot$epithelial_state)), "\n", sep = "")
cat("genes=", length(unique(dot$gene)), "\n", sep = "")
cat("rows=", nrow(dot), "\n", sep = "")
cat("missing=", if (length(missing_genes)) paste(missing_genes, collapse = ",") else "none", "\n", sep = "")
