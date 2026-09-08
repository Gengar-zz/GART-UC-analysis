# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 12_FIG2C_D_STYLE_REBUILD/CODE/plot_fig2cd_style_rebuild.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors = FALSE)
set.seed(20260905)

project_root <- paths$repo_root
revision_root <- file.path(paths$results_root, "singlecell")
singlecell_root <- file.path(paths$external_data_root, "singlecell", "frozen_objects")
out_dir <- file.path(paths$results_root, "figures", "Figure2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(SeuratObject)
  library(ggplot2)
})
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

ep_path <- file.path(singlecell_root, "UC_GART_06c_Epithelial_Fig2AB_noGARTname_FINAL.rds")
isc_path <- file.path(singlecell_root, "UC_GART_07_ISC_like_subclustered_FIXED_withDonor.rds")
corrected_path <- file.path(revision_root, "INTERMEDIATE", "ISC_corrected_cell_metadata.csv")
fig2a_path <- file.path(revision_root, "INTERMEDIATE", "Fig2A_epithelial_UMAP_plot_data.csv")
invisible(lapply(c(ep_path, isc_path, corrected_path, fig2a_path), gart_require_file, label = "Figure 2C/D input"))

ep <- readRDS(ep_path)
isc <- readRDS(isc_path)
cm <- read.csv(corrected_path, check.names = FALSE)
fig2a <- read.csv(fig2a_path, check.names = FALSE)
em <- ep[[]]
ep_umap <- Embeddings(ep, "umap_epi")
isc_umap <- Embeddings(isc, "iscUMAP")

stopifnot(identical(rownames(em), rownames(ep_umap)))
stopifnot(nrow(cm) == 9306L, nrow(isc_umap) == 9306L)
stopifnot(!anyDuplicated(cm$cell_id_revision), !anyDuplicated(rownames(isc_umap)))
stopifnot(setequal(cm$cell_id_revision, rownames(isc_umap)))

# Read existing raw count layers without imputing missing values as zero.
assay <- ep[["RNA"]]
count_layers <- grep("^counts", Layers(assay), value = TRUE)
raw_gart <- setNames(rep(NA_real_, ncol(ep)), colnames(ep))
for (layer_name in count_layers) {
  q <- methods::slot(assay, "layers")[[layer_name]]
  feature_names <- methods::slot(assay, "features")[[layer_name]]
  cell_names <- methods::slot(assay, "cells")[[layer_name]]
  if (is.null(dim(q))) q <- matrix(q, nrow = length(feature_names), ncol = length(cell_names))
  if (is.null(rownames(q))) rownames(q) <- feature_names
  if (is.null(colnames(q))) colnames(q) <- cell_names
  if ("GART" %in% rownames(q)) raw_gart[colnames(q)] <- as.numeric(q["GART", ])
}
if (anyNA(raw_gart)) stop("At least one epithelial cell lacks a recoverable raw GART count")

isc_state <- "Stem/ISC-like epithelial cells"
frozen_state <- as.character(em$Epi_Subtype_Fig2)
isc_ids <- rownames(em)[frozen_state == isc_state]
if (length(isc_ids) != 9306L || !setequal(isc_ids, cm$cell_id_revision)) {
  stop("Frozen epithelial ISC membership does not match corrected 9,306-cell metadata")
}
cm_idx <- match(isc_ids, cm$cell_id_revision)
if (anyNA(cm_idx)) stop("Failed to match corrected ISC metadata")
if (any(raw_gart[isc_ids] != cm$GART_count_corrected[cm_idx])) {
  stop("Frozen raw GART counts differ from corrected ISC metadata")
}

cohort <- as.character(em$Dataset_Final)
donor_id <- ifelse(
  cohort == "GSE214695",
  as.character(em$Sample_Final),
  paste0("SCP259_", as.character(em$Subject))
)
if (anyNA(cohort) || anyNA(donor_id) || any(!nzchar(donor_id))) stop("Missing cohort or donor identifiers")

status <- rep("Other epithelial cells", nrow(em))
names(status) <- rownames(em)
status[isc_ids[raw_gart[isc_ids] == 0]] <- "GART-undetected ISC-like cells"
status[isc_ids[raw_gart[isc_ids] >= 1]] <- "GART-detected ISC-like cells"
status_levels <- c("Other epithelial cells", "GART-undetected ISC-like cells", "GART-detected ISC-like cells")

fig2c <- data.frame(
  cell_id = rownames(em),
  epithelial_state = frozen_state,
  cohort = cohort,
  donor_id = donor_id,
  UMAP_1 = ep_umap[, 1],
  UMAP_2 = ep_umap[, 2],
  raw_GART_count = as.numeric(raw_gart[rownames(em)]),
  GART_detection_status = factor(status[rownames(em)], levels = status_levels),
  stringsAsFactors = FALSE
)
write.csv(fig2c, file.path(out_dir, "Fig2C_GART_detection_source_data.csv"), row.names = FALSE)

# Reproduce the Figure 2A plotting limits (its plotted coordinate range plus base-R 4% expansion).
expand_base_r <- function(x) {
  r <- range(x, finite = TRUE)
  r + c(-1, 1) * diff(r) * 0.04
}
fig2a_xlim <- expand_base_r(fig2a$UMAP1)
fig2a_ylim <- expand_base_r(fig2a$UMAP2)

other <- fig2c[fig2c$GART_detection_status == "Other epithelial cells", ]
undet <- fig2c[fig2c$GART_detection_status == "GART-undetected ISC-like cells", ]
det <- fig2c[fig2c$GART_detection_status == "GART-detected ISC-like cells", ]
gart_pal <- c(
  "Other epithelial cells" = "#C7C7C7",
  "GART-undetected ISC-like cells" = "#79A9D1",
  "GART-detected ISC-like cells" = "#8F1537"
)

p_c <- ggplot() +
  geom_point(data = other, aes(UMAP_1, UMAP_2, colour = GART_detection_status), size = 0.24, alpha = 0.24) +
  geom_point(data = undet, aes(UMAP_1, UMAP_2, colour = GART_detection_status), size = 0.80, alpha = 0.80) +
  geom_point(data = det, aes(UMAP_1, UMAP_2, colour = GART_detection_status), size = 0.80, alpha = 0.90) +
  scale_colour_manual(values = gart_pal, breaks = status_levels, drop = FALSE, name = NULL) +
  coord_fixed(xlim = fig2a_xlim, ylim = fig2a_ylim, expand = FALSE, clip = "on") +
  labs(x = "UMAP 1", y = "UMAP 2", title = "GART detection within the ISC-like epithelial compartment") +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 7)),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5, colour = "black"),
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    legend.position = "right",
    legend.text = element_text(size = 8.5),
    legend.key.height = unit(0.22, "in"),
    plot.margin = margin(7, 8, 7, 7)
  ) +
  guides(colour = guide_legend(override.aes = list(size = c(2.0, 2.8, 2.8), alpha = 1)))

p_c_no_legend <- p_c + theme(legend.position = "none")

# Figure 2D: frozen ISC embedding and frozen cluster assignments.
isc_idx <- match(rownames(isc_umap), cm$cell_id_revision)
cluster_levels <- sort(unique(as.character(cm$ISC_subcluster[isc_idx])))
cluster_pal <- setNames(hcl.colors(length(cluster_levels), "Set 2"), cluster_levels)
fig2d <- data.frame(
  cell_id = rownames(isc_umap),
  cohort = cm$cohort[isc_idx],
  donor_id = cm$donor_id[isc_idx],
  UMAP_1 = isc_umap[, 1],
  UMAP_2 = isc_umap[, 2],
  final_ISClike_cluster = factor(as.character(cm$ISC_subcluster[isc_idx]), levels = cluster_levels),
  stringsAsFactors = FALSE
)
fig2d$cluster_color <- unname(cluster_pal[as.character(fig2d$final_ISClike_cluster)])
write.csv(fig2d, file.path(out_dir, "Fig2D_ISClike_subcluster_source_data.csv"), row.names = FALSE)

centers <- do.call(rbind, lapply(cluster_levels, function(cl) {
  d <- fig2d[fig2d$final_ISClike_cluster == cl, ]
  data.frame(final_ISClike_cluster = cl, UMAP_1 = median(d$UMAP_1), UMAP_2 = median(d$UMAP_2), n_cells = nrow(d))
}))

p_d <- ggplot(fig2d, aes(UMAP_1, UMAP_2, colour = final_ISClike_cluster)) +
  geom_point(size = 0.80, alpha = 0.90) +
  geom_label(
    data = centers,
    aes(UMAP_1, UMAP_2, label = final_ISClike_cluster),
    inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 3.8,
    colour = "black", fill = "white", alpha = 0.78,
    linewidth = 0, label.padding = unit(0.11, "lines")
  ) +
  scale_colour_manual(values = cluster_pal, breaks = cluster_levels, drop = FALSE, name = "Cluster") +
  coord_fixed(expand = TRUE) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "ISC-like epithelial subclusters") +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 7)),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5, colour = "black"),
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    legend.position = "none",
    plot.margin = margin(7, 8, 7, 7)
  )
p_d_legend <- p_d +
  theme(legend.position = "right", legend.text = element_text(size = 8.5), legend.title = element_text(face = "bold", size = 9)) +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1), ncol = 1))

save_pair <- function(plot, stem, width, height) {
  grDevices::tiff(file.path(out_dir, paste0(stem, ".tiff")), width = width, height = height,
                  units = "in", res = 300, compression = "lzw", bg = "white",
                  type = if (.Platform$OS.type == "windows") "windows" else "cairo",
                  antialias = if (.Platform$OS.type == "windows") "cleartype" else "default")
  print(plot); dev.off()
  grDevices::cairo_pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width, height = height,
                       family = "Arial", bg = "white", onefile = TRUE)
  print(plot); dev.off()
}

save_pair(p_c, "Fig2C_GART_detection_within_ISClike_compartment", 6.9, 5.5)
save_pair(p_c_no_legend, "Fig2C_GART_detection_within_ISClike_compartment_no_legend", 6.9, 5.5)
save_pair(p_d, "Fig2D_ISClike_epithelial_subclusters_clean", 5.8, 5.5)
save_pair(p_d_legend, "Fig2D_ISClike_epithelial_subclusters_with_legend", 6.6, 5.5)

qa <- data.frame(
  check = c(
    "total_epithelial_cells", "ISC_like_cells", "GART_detected_ISC_cells", "GART_undetected_ISC_cells",
    "missing_raw_GART_counts_all_epithelial", "missing_raw_GART_counts_ISC", "ISC_count_sum",
    "ISC_cluster_count", "ISC_embedding_coordinate_max_abs_difference_vs_prior_export",
    "Fig2A_x_limit_low", "Fig2A_x_limit_high", "Fig2A_y_limit_low", "Fig2A_y_limit_high"
  ),
  value = c(
    nrow(fig2c), nrow(undet) + nrow(det), nrow(det), nrow(undet), sum(is.na(fig2c$raw_GART_count)),
    sum(is.na(fig2c$raw_GART_count[fig2c$epithelial_state == isc_state])), nrow(det) + nrow(undet),
    length(cluster_levels),
    max(abs(c(isc_umap[,1] - fig2d$UMAP_1, isc_umap[,2] - fig2d$UMAP_2))),
    fig2a_xlim[1], fig2a_xlim[2], fig2a_ylim[1], fig2a_ylim[2]
  )
)
write.csv(qa, file.path(out_dir, "Fig2C_D_data_QA.csv"), row.names = FALSE)
write.csv(centers, file.path(out_dir, "Fig2D_cluster_label_positions.csv"), row.names = FALSE)

cat(sprintf("detected=%d undetected=%d clusters=%d missing_raw=%d\n", nrow(det), nrow(undet), length(cluster_levels), sum(is.na(fig2c$raw_GART_count))))
