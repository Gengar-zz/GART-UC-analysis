# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 16_FIG3_ORIGINAL_STYLE_REBUILD/CODE/plot_fig3_original_style_rebuild.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors = FALSE)
set.seed(20260905)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

input_root <- file.path(paths$results_root, "singlecell")
out_root <- file.path(paths$results_root, "figures", "Figure3")
fig_dir <- file.path(out_root, "FIGURES")
src_dir <- file.path(out_root, "SOURCE_DATA")
report_dir <- file.path(out_root, "REPORTS")
qa_dir <- file.path(out_root, "QA")
code_dir <- file.path(out_root, "CODE")
int_dir <- file.path(out_root, "INTERMEDIATE")
invisible(lapply(c(fig_dir, src_dir, report_dir, qa_dir, code_dir, int_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

if (.Platform$OS.type == "windows") suppressWarnings(windowsFonts(Arial = windowsFont("Arial")))

COL <- list(
  burgundy = "#951D45",
  pink = "#E3A1B2",
  blue = "#3F76A8",
  pale_blue = "#B8D2E5",
  grey = "#D7D7D7",
  light_grey = "#E9E9E9",
  dark = "#222222"
)

theme_panel <- theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "#5F6871", margin = margin(b = 6)),
    axis.title = element_text(size = 9.5),
    axis.text = element_text(size = 8.5, color = "black"),
    axis.line = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(7, 9, 7, 7)
  )

fmt_p <- function(x, digits = 3) {
  ifelse(is.na(x), "NA",
         ifelse(x < 1e-10, "<1×10^-10",
                ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
                       formatC(x, format = "f", digits = digits))))
}

clean_pathway <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

wrap_label <- function(x, width = 34) vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))

save_pair <- function(plot, stem, width, height) {
  tiff(file.path(fig_dir, paste0(stem, ".tiff")), width = width, height = height,
       units = "in", res = 300, compression = "lzw", type = "windows", bg = "white")
  print(plot)
  dev.off()
  cairo_pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = width, height = height,
            family = "Arial", bg = "white", onefile = TRUE)
  print(plot)
  dev.off()
}

## Frozen inputs
umap <- read.csv(file.path(input_root, "INTERMEDIATE", "ISC_UMAP_plot_data.csv"), check.names = FALSE)
meta <- read.csv(file.path(input_root, "INTERMEDIATE", "ISC_corrected_cell_metadata.csv"), check.names = FALSE)
deg <- read.csv(file.path(input_root, "RESULTS", "SCP259_primary_pseudobulk_DEG.csv"), check.names = FALSE)
gsea <- read.csv(file.path(input_root, "RESULTS", "SCP259_primary_pseudobulk_GSEA.csv"), check.names = FALSE)
paired_tests <- read.csv(file.path(input_root, "RESULTS", "primary_module_paired_donor_tests.csv"), check.names = FALSE)
paired_data <- read.csv(file.path(input_root, "RESULTS", "primary_module_paired_donor_data.csv"), check.names = FALSE)
cycle <- read.csv(file.path(input_root, "RESULTS", "module_cell_cycle_stratified.csv"), check.names = FALSE)

stopifnot(nrow(umap) == 9306L, nrow(meta) == 9306L, !anyNA(meta$GART_count_corrected))
stopifnot(!anyDuplicated(umap$cell_id), !anyDuplicated(meta$cell_id_revision))
stopifnot(all(deg$n_paired_donors == 11L), all(paired_tests$n_paired_donors == 23L))

module_keep <- c(
  "Reactome_purine_GART_excluded",
  "Reactome_DNA_replication",
  "Reactome_DNA_repair",
  "Reactome_mitotic_cell_cycle"
)
module_labels <- c(
  "Reactome_purine_GART_excluded" = "Purine biosynthesis\n(GART excluded)",
  "Reactome_DNA_replication" = "DNA replication",
  "Reactome_DNA_repair" = "DNA repair",
  "Reactome_mitotic_cell_cycle" = "Mitotic cell cycle"
)
gene_pool <- c("PCNA", "MCM6", "MCM7", "FEN1", "TK1", "TYMS", "GINS2", "PAICS", "MKI67", "CDK4", "IMPDH2")

## Fig3A
idx <- match(umap$cell_id, meta$cell_id_revision)
stopifnot(!anyNA(idx))
umap$GART_count_corrected <- meta$GART_count_corrected[idx]
umap$GART_detection_status <- ifelse(umap$GART_count_corrected >= 1, "GART-detected", "GART-undetected")
uc <- umap[umap$condition == "UC", ]
uc <- uc[order(uc$GART_detection_status == "GART-detected"), ]
uc$GART_detection_status <- factor(uc$GART_detection_status, levels = c("GART-undetected", "GART-detected"))
xlim_frozen <- range(umap$UMAP1, finite = TRUE)
ylim_frozen <- range(umap$UMAP2, finite = TRUE)

pA <- ggplot(uc, aes(UMAP1, UMAP2, color = GART_detection_status)) +
  geom_point(size = 0.8, alpha = 0.88, stroke = 0) +
  scale_color_manual(values = c("GART-undetected" = "#D3D3D3", "GART-detected" = COL$burgundy)) +
  coord_fixed(xlim = xlim_frozen, ylim = ylim_frozen, clip = "off") +
  labs(title = "GART detection in UC ISC-like epithelial cells", x = "UMAP 1", y = "UMAP 2", color = NULL) +
  theme_panel +
  theme(legend.position = "right")

fig3a_source <- data.frame(
  cell_id = uc$cell_id, cohort = uc$cohort, condition = uc$condition, donor_id = uc$donor_id,
  UMAP_1 = uc$UMAP1, UMAP_2 = uc$UMAP2,
  raw_corrected_GART_count = uc$GART_count_corrected,
  GART_detection_status = as.character(uc$GART_detection_status),
  n_donors = length(unique(uc$donor_id)),
  analysis_type = "Frozen UC ISC-like UMAP visualization",
  effect = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_, FDR = NA_real_,
  selection_display_rule = "All frozen UC ISC-like cells; raw corrected count 0 versus >=1; detected cells drawn last"
)
write.csv(fig3a_source, file.path(src_dir, "Fig3A_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
save_pair(pA, "Fig3A_UC_ISClike_by_GART_detection_clean", 6.1, 5.4)

## Fig3B
deg$plot_class <- ifelse(deg$FDR < 0.05 & deg$logFC > 0, "FDR <0.05, positive",
                         ifelse(deg$FDR < 0.05 & deg$logFC < 0, "FDR <0.05, negative", "Not FDR-significant"))
deg$minus_log10_FDR <- -log10(pmax(deg$FDR, .Machine$double.xmin))
candidate_df <- deg[deg$gene %in% gene_pool & deg$FDR < 0.05, ]
candidate_df <- candidate_df[order(candidate_df$FDR, -candidate_df$logFC), ]
label_candidates <- head(candidate_df$gene, 8)
deg$labelled_gene <- deg$gene %in% label_candidates
deg$displayed_in_volcano <- deg$gene != "GART"
deg$analysis_type <- "SCP259 paired-donor edgeR quasi-likelihood pseudobulk"
deg$effect <- deg$logFC
deg$P <- deg$raw_P
deg$effect_definition <- "log2 fold change: GART-detected minus GART-undetected"
deg$selection_display_rule <- paste0(
  "All genes except GART displayed; color by within-run FDR<0.05 and logFC direction; labels restricted to at most 8 genes from the pre-specified candidate pool, ordered by frozen FDR. ",
  "GART omitted from the volcano display because it defines the compared groups but remains in source data."
)
write.csv(deg, file.path(src_dir, "Fig3B_pseudobulk_DEG_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

deg_plot <- deg[deg$displayed_in_volcano, ]
lab_plot <- deg_plot[deg_plot$labelled_gene, ]
lab_plot <- lab_plot[order(lab_plot$minus_log10_FDR, decreasing = TRUE), ]
lab_plot$label_y <- seq(max(deg_plot$minus_log10_FDR) * 0.96,
                        max(deg_plot$minus_log10_FDR) * 0.36,
                        length.out = nrow(lab_plot))
lab_plot$label_x <- max(deg_plot$logFC, na.rm = TRUE) + 0.12
pB <- ggplot(deg_plot, aes(logFC, minus_log10_FDR, color = plot_class)) +
  geom_point(size = 1.0, alpha = 0.72, stroke = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "#727B83") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 0.4, color = "#727B83") +
  geom_segment(data = lab_plot, aes(x = logFC, y = minus_log10_FDR, xend = label_x, yend = label_y),
               inherit.aes = FALSE, color = "#8A8A8A", linewidth = 0.28) +
  geom_text(data = lab_plot, aes(x = label_x, y = label_y, label = gene),
            inherit.aes = FALSE, family = "Arial", size = 2.9,
            color = "#202020", hjust = 0) +
  annotate("label", x = Inf, y = Inf, label = "11 paired donors", hjust = 1.05, vjust = 1.25,
           family = "Arial", size = 3.0, fill = "white") +
  scale_color_manual(values = c("FDR <0.05, positive" = COL$burgundy,
                                "FDR <0.05, negative" = COL$blue,
                                "Not FDR-significant" = "#CECECE"),
                     breaks = c("FDR <0.05, positive", "FDR <0.05, negative", "Not FDR-significant")) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.26))) +
  labs(title = "SCP259 donor-paired pseudobulk differential expression",
       x = "log2 fold change (GART-detected − GART-undetected)",
       y = "−log10 FDR", color = NULL) +
  theme_panel +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 1))
save_pair(pB, "Fig3B_SCP259_paired_donor_pseudobulk_volcano", 6.6, 5.4)

## Fig3C
hallmark <- gsea[gsea$collection == "Hallmark", ]
hallmark <- hallmark[order(-hallmark$NES), ]
hallmark_sig_pos <- hallmark[hallmark$padj < 0.05 & hallmark$NES > 0, ]
hallmark_display <- head(hallmark_sig_pos, 8)
hallmark$displayed <- hallmark$pathway %in% hallmark_display$pathway
hallmark$analysis_type <- "SCP259 11-paired-donor pseudobulk ranking; fgseaMultilevel Hallmark"
hallmark$n_donors <- 11L
hallmark$effect <- hallmark$NES
hallmark$CI_low <- NA_real_
hallmark$CI_high <- NA_real_
hallmark$P <- hallmark$pval
hallmark$FDR <- hallmark$padj
hallmark$selection_display_rule <- "Pre-specified: FDR<0.05, positive NES, sorted descending by NES, maximum 8; full Hallmark results retained in this source table"
write.csv(hallmark, file.path(src_dir, "Fig3C_Hallmark_GSEA_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

hallmark_display$display_label <- wrap_label(clean_pathway(hallmark_display$pathway), 28)
hallmark_display$display_label <- factor(hallmark_display$display_label, levels = rev(hallmark_display$display_label))
hallmark_display$fdr_label <- paste0("FDR ", fmt_p(hallmark_display$padj, 3))
pC <- ggplot(hallmark_display, aes(NES, display_label)) +
  geom_segment(aes(x = 0, xend = NES, yend = display_label), linewidth = 1.2, color = "#D8A4B3") +
  geom_point(size = 3.1, color = COL$burgundy) +
  geom_text(aes(label = fdr_label), hjust = -0.08, family = "Arial", size = 2.8) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.32))) +
  labs(title = "Hallmark pathway enrichment", x = "Normalized enrichment score", y = NULL) +
  theme_panel +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_pair(pC, "Fig3C_Hallmark_GSEA_paired_pseudobulk", 6.5, 5.3)

## Fig3D: exact pre-specified canonical pathways, evaluated without post hoc substitution
reactome_fixed <- c(
  "REACTOME_DNA_REPLICATION",
  "REACTOME_DNA_REPAIR",
  "REACTOME_CELL_CYCLE_MITOTIC",
  "REACTOME_CHROMOSOME_MAINTENANCE",
  "REACTOME_METABOLISM_OF_NUCLEOTIDES",
  "REACTOME_CELL_CYCLE_CHECKPOINTS"
)
go_fixed <- c(
  "GOBP_DNA_REPLICATION",
  "GOBP_DNA_REPAIR",
  "GOBP_CHROMOSOME_SEGREGATION",
  "GOBP_CHROMOSOME_ORGANIZATION",
  "GOBP_PURINE_NUCLEOTIDE_BIOSYNTHETIC_PROCESS",
  "GOBP_MITOTIC_CELL_CYCLE_PHASE_TRANSITION"
)
gsea_rg <- gsea[gsea$collection %in% c("Reactome", "GO_BP"), ]
gsea_rg$predefined_candidate <- gsea_rg$pathway %in% c(reactome_fixed, go_fixed)
gsea_rg$displayed <- gsea_rg$predefined_candidate & gsea_rg$padj < 0.05
gsea_rg$analysis_type <- "SCP259 11-paired-donor pseudobulk ranking; fgseaMultilevel"
gsea_rg$n_donors <- 11L
gsea_rg$effect <- gsea_rg$NES
gsea_rg$CI_low <- NA_real_
gsea_rg$CI_high <- NA_real_
gsea_rg$P <- gsea_rg$pval
gsea_rg$FDR <- gsea_rg$padj
gsea_rg$selection_display_rule <- "Exact pre-specified canonical replication, repair, cell-cycle, chromosome and nucleotide/purine pathways; displayed only if frozen FDR<0.05; maximum 6 per database"
write.csv(gsea_rg, file.path(int_dir, "Fig3D_full_source.csv"), row.names = FALSE, fileEncoding = "UTF-8")
gsea_d <- gsea_rg[gsea_rg$displayed, ]
stopifnot(sum(gsea_d$collection == "Reactome") <= 6, sum(gsea_d$collection == "GO_BP") <= 6)
gsea_d$database <- ifelse(gsea_d$collection == "GO_BP", "GO Biological Process", "Reactome")
gsea_d$display_label <- wrap_label(clean_pathway(gsea_d$pathway), 35)
gsea_d$fdr_label <- paste0("FDR ", fmt_p(gsea_d$padj, 3))
gsea_d <- gsea_d[order(gsea_d$database, gsea_d$NES), ]
gsea_d$row_id <- factor(seq_len(nrow(gsea_d)), levels = seq_len(nrow(gsea_d)))
write.csv(gsea_d, file.path(int_dir, "Fig3D_displayed_source.csv"), row.names = FALSE, fileEncoding = "UTF-8")

pD <- ggplot(gsea_d, aes(NES, row_id)) +
  geom_segment(aes(x = 0, xend = NES, yend = row_id), linewidth = 1.0, color = "#D8A4B3") +
  geom_point(size = 2.8, color = COL$burgundy) +
  geom_text(aes(label = fdr_label), hjust = -0.08, family = "Arial", size = 2.55) +
  facet_grid(database ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = setNames(gsea_d$display_label, as.character(gsea_d$row_id))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.34))) +
  labs(title = "Replication and repair pathway enrichment", x = "Normalized enrichment score", y = NULL) +
  theme_panel +
  theme(
    strip.background = element_rect(fill = "#F2F2F2", color = "#BDBDBD", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 8.7),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    panel.spacing.y = unit(0.35, "lines")
  )
save_pair(pD, "Fig3D_Reactome_GO_GSEA_paired_pseudobulk", 7.2, 6.7)

## Fig3E
tests4 <- paired_tests[match(module_keep, paired_tests$module), ]
stopifnot(all(tests4$module == module_keep), all(tests4$n_paired_donors == 23L))
paired4 <- paired_data[paired_data$module %in% module_keep & is.finite(paired_data$delta), ]
stopifnot(all(table(paired4$module) == 23L))

paired_long <- rbind(
  data.frame(donor_id = paired4$donor_id, cohort = paired4$cohort, module = paired4$module,
             status = "GART-undetected", score = paired4$mean_score.negative),
  data.frame(donor_id = paired4$donor_id, cohort = paired4$cohort, module = paired4$module,
             status = "GART-detected", score = paired4$mean_score.positive)
)
paired_long$status <- factor(paired_long$status, levels = c("GART-undetected", "GART-detected"))
paired_long$module_label <- factor(module_labels[paired_long$module], levels = unname(module_labels[module_keep]))
test_idx <- match(paired_long$module, tests4$module)
paired_long$n_donors <- tests4$n_paired_donors[test_idx]
paired_long$analysis_type <- "Frozen paired-donor module summary; paired Wilcoxon used uniformly"
paired_long$effect <- tests4$mean_delta[test_idx]
paired_long$CI_low <- tests4$CI_low[test_idx]
paired_long$CI_high <- tests4$CI_high[test_idx]
paired_long$P <- tests4$paired_Wilcoxon_P[test_idx]
paired_long$FDR <- NA_real_
paired_long$selection_display_rule <- "Four pre-specified public modules; all 23 complete paired donors displayed; no cell-level points"
write.csv(paired_long, file.path(src_dir, "Fig3E_paired_donor_modules_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

annE <- data.frame(
  module_label = factor(unname(module_labels[tests4$module]), levels = unname(module_labels[module_keep])),
  label = sprintf("Δ=%.3f [%.3f, %.3f]\npaired P=%s",
                  tests4$mean_delta, tests4$CI_low, tests4$CI_high, fmt_p(tests4$paired_Wilcoxon_P, 4))
)
pE <- ggplot(paired_long, aes(status, score, group = donor_id)) +
  geom_line(color = "#BFC3C7", linewidth = 0.42, alpha = 0.8) +
  geom_point(aes(color = status), size = 1.65, alpha = 0.92) +
  geom_text(data = annE, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE,
            family = "Arial", size = 2.9, vjust = 1.12, lineheight = 0.95) +
  facet_wrap(~module_label, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("GART-undetected" = COL$pale_blue, "GART-detected" = COL$burgundy)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.36))) +
  labs(title = "Donor-level functional programs associated with GART detection",
       x = NULL, y = "Mean module score", color = NULL) +
  theme_panel +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 13)),
    strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.8),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8.2),
    legend.position = "none", panel.spacing.x = unit(0.65, "lines")
  )
save_pair(pE, "Fig3E_paired_donor_core_module_scores", 11.2, 5.2)

## Fig3F
cycle4 <- cycle[cycle$module %in% module_keep, ]
cycle4 <- cycle4[match(as.vector(outer(module_keep, c("G1", "S", "G2M"), paste, sep = "__")),
                       paste(cycle4$module, cycle4$phase, sep = "__")), ]
stopifnot(nrow(cycle4) == 12L, !anyNA(cycle4$module))
cycle4$module_label <- factor(module_labels[cycle4$module], levels = rev(unname(module_labels[module_keep])))
cycle4$phase_label <- factor(cycle4$phase, levels = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M"))
cycle4$cell_label <- sprintf("%.3f\nFDR %s", cycle4$estimate, fmt_p(cycle4$FDR, 3))
cycle4$n_donors_panel <- cycle4$n_donors
cycle4$analysis_type <- "Frozen donor-aware depth-adjusted model within each cell-cycle phase"
cycle4$effect <- cycle4$estimate
cycle4$selection_display_rule <- "All four pre-specified modules and all G1, S and G2/M strata displayed irrespective of P or FDR"
write.csv(cycle4, file.path(src_dir, "Fig3F_cell_cycle_stratified_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

limF <- max(abs(c(cycle4$estimate, cycle4$CI_low, cycle4$CI_high)), na.rm = TRUE)
pF <- ggplot(cycle4, aes(phase_label, module_label, fill = estimate)) +
  geom_tile(color = "white", linewidth = 1.0) +
  geom_text(aes(label = cell_label), family = "Arial", size = 3.0, lineheight = 0.95) +
  scale_fill_gradient2(low = "#3A78B4", mid = "white", high = "#B22645", midpoint = 0,
                       limits = c(-limF, limF), name = "Adjusted\neffect") +
  labs(title = "Cell-cycle-stratified module effects", x = "Cell-cycle phase", y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 9)),
    axis.text = element_text(size = 9, color = "black"), axis.title = element_text(size = 9.5),
    panel.grid = element_blank(), legend.title = element_text(size = 8.5), legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 8, 8, 8)
  )
save_pair(pF, "Fig3F_cell_cycle_stratified_core_module_effects", 7.5, 4.8)

## Fig3G: donor-specific pseudobulk differences were not archived; do not reconstruct from expression matrices
scp <- meta[meta$cohort == "SCP259", ]
scp$state <- ifelse(scp$GART_count_corrected >= 1, "detected", "undetected")
tab <- table(scp$donor_id, scp$state)
paired11 <- rownames(tab)[tab[, "detected"] >= 20 & tab[, "undetected"] >= 20]
stopifnot(length(paired11) == 11L)
gene_deg <- deg[match(gene_pool, deg$gene), ]
stopifnot(all(gene_deg$gene == gene_pool))
fig3g_source <- expand.grid(donor_id = paired11, gene = gene_pool, stringsAsFactors = FALSE)
gi <- match(fig3g_source$gene, gene_deg$gene)
fig3g_source$standardized_pseudobulk_difference <- NA_real_
fig3g_source$aggregate_log2FC <- gene_deg$logFC[gi]
fig3g_source$effect <- gene_deg$logFC[gi]
fig3g_source$CI_low <- gene_deg$CI_low[gi]
fig3g_source$CI_high <- gene_deg$CI_high[gi]
fig3g_source$P <- gene_deg$raw_P[gi]
fig3g_source$FDR <- gene_deg$FDR[gi]
fig3g_source$n_donors <- 11L
fig3g_source$analysis_type <- "Frozen SCP259 paired-donor pseudobulk DEG; donor-specific difference unavailable"
fig3g_source$selection_display_rule <- "Fixed 11-gene candidate pool retained without significance filtering; donor-specific values left blank because they were not archived and were not reconstructed"
fig3g_source$missing_value_reason <- "Donor-specific normalized pseudobulk differences were not present in frozen outputs; no expression-matrix re-aggregation was performed"
write.csv(fig3g_source, file.path(src_dir, "Fig3G_paired_donor_gene_heatmap_source_data.csv"), row.names = FALSE,
          na = "", fileEncoding = "UTF-8")

gene_axis <- paste0(gene_pool, "\n", ifelse(gene_deg$FDR < 0.05,
                                              paste0("FDR ", fmt_p(gene_deg$FDR, 3)), "ns"))
names(gene_axis) <- gene_pool
blank_heat <- fig3g_source
blank_heat$donor_id <- factor(blank_heat$donor_id, levels = rev(paired11))
blank_heat$gene <- factor(blank_heat$gene, levels = gene_pool)
pG <- ggplot(blank_heat, aes(gene, donor_id)) +
  geom_tile(fill = "#F3F3F3", color = "white", linewidth = 0.65) +
  annotate("label", x = 6, y = 6,
           label = "Donor-specific pseudobulk differences\nwere not archived; values were not reconstructed.",
           family = "Arial", size = 3.6, color = "#4F5962", fill = "white") +
  scale_x_discrete(labels = gene_axis) +
  labs(title = "Donor-level expression of proliferation and nucleotide-metabolism genes",
       subtitle = "Fixed candidate pool; aggregate DEG FDR shown below each gene",
       x = NULL, y = "SCP259 paired donor") +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "#5F6871", margin = margin(b = 7)),
    axis.text.x = element_text(size = 7.4, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 7.5, color = "black"), axis.title.y = element_text(size = 9),
    panel.grid = element_blank(), plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 10, 8, 8)
  )
save_pair(pG, "Fig3G_paired_donor_replication_purine_gene_heatmap", 8.4, 5.8)

## Layout previews (not final panels)
pv1_top <- plot_grid(pA, pB, pC, labels = c("A", "B", "C"), nrow = 1, rel_widths = c(0.9, 1.05, 1.05))
pv1_mid <- plot_grid(pD, pE, labels = c("D", "E"), nrow = 1, rel_widths = c(0.9, 1.35))
pv1_bot <- plot_grid(pF, pG, labels = c("F", "G"), nrow = 1, rel_widths = c(0.9, 1.1))
preview1 <- plot_grid(pv1_top, pv1_mid, pv1_bot, ncol = 1, rel_heights = c(1, 1.1, 1))

pv2_top <- plot_grid(pA, pB, labels = c("A", "B"), nrow = 1)
pv2_second <- plot_grid(pC, pD, labels = c("C", "D"), nrow = 1, rel_widths = c(0.9, 1.1))
pv2_third <- plot_grid(NULL, pE, NULL, labels = c("", "E", ""), nrow = 1, rel_widths = c(0.05, 0.9, 0.05))
pv2_bottom <- plot_grid(pF, pG, labels = c("F", "G"), nrow = 1, rel_widths = c(0.9, 1.1))
preview2 <- plot_grid(pv2_top, pv2_second, pv2_third, pv2_bottom, ncol = 1, rel_heights = c(1, 1.15, 0.85, 1))

ggsave(file.path(fig_dir, "Fig3_layout_preview_version1.png"), preview1, width = 18, height = 14, dpi = 150, bg = "white")
ggsave(file.path(fig_dir, "Fig3_layout_preview_version1.pdf"), preview1, width = 18, height = 14, device = cairo_pdf, family = "Arial", bg = "white")
ggsave(file.path(fig_dir, "Fig3_layout_preview_version2.png"), preview2, width = 14, height = 18, dpi = 150, bg = "white")
ggsave(file.path(fig_dir, "Fig3_layout_preview_version2.pdf"), preview2, width = 14, height = 18, device = cairo_pdf, family = "Arial", bg = "white")

## Reports
hallmark_names <- paste(clean_pathway(hallmark_display$pathway), collapse = "; ")
d_names <- paste(paste0(gsea_d$database, ": ", clean_pathway(gsea_d$pathway)), collapse = "; ")
summary_text <- paste0(
  "# Figure 3 original-style rebuild summary\n\n",
  "## Scope\n\n",
  "Visualization-only rebuild from frozen revised outputs. No model, differential-expression test, GSEA, module score, cell-cycle model, expression matrix, donor map, threshold, seed, P value or FDR was recomputed or changed.\n\n",
  "## Frozen analysis mapping\n\n",
  "- Fig3A: corrected raw-count GART detection on the frozen UC ISC-like UMAP; coordinates and orientation match the frozen ISC-like embedding.\n",
  "- Fig3B–D: SCP259 paired-donor pseudobulk analysis with 11 donors and the archived edgeR/fgsea results.\n",
  "- Fig3E: four public modules using the same archived paired-donor Wilcoxon framework and 23 complete paired donors.\n",
  "- Fig3F: all 12 archived module-by-phase donor-aware effects; no stratum was removed by significance.\n",
  "- Fig3G: the 11-gene pool was fixed before display. Frozen outputs contain aggregate edgeR effects but not donor-specific normalized pseudobulk differences. All 121 donor-by-gene cells were therefore retained as blank, rather than reconstructed from the expression matrix. This panel is not recommended for the main figure unless an archived donor-level matrix is supplied.\n\n",
  "## Display selections\n\n",
  "Hallmark pathways: ", hallmark_names, ".\n\n",
  "Reactome/GO pathways: ", d_names, ".\n\n",
  "## Main-panel recommendation\n\n",
  "Use Fig3A–F. Treat Fig3G as an availability audit placeholder, not a publication-ready biological heatmap, because the required frozen donor-specific values were not archived.\n"
)
writeLines(summary_text, file.path(report_dir, "FIG3_ORIGINAL_STYLE_REBUILD_SUMMARY.md"), useBytes = TRUE)

manifest_text <- paste0(
  "# Figure 3 main-versus-supplementary manifest\n\n",
  "| Item | Placement | Frozen basis | Status |\n|---|---|---|---|\n",
  "| Fig3A GART detection UMAP | Main | Corrected raw count and frozen ISC-like UMAP | Rebuilt |\n",
  "| Fig3B paired pseudobulk volcano | Main | SCP259, 11 paired donors | Rebuilt |\n",
  "| Fig3C Hallmark GSEA | Main | Frozen Hallmark FDR/NES | Rebuilt |\n",
  "| Fig3D Reactome/GO GSEA | Main | Exact pre-specified pathways passing frozen FDR<0.05 | Rebuilt |\n",
  "| Fig3E paired module scores | Main | 23 paired donors; paired Wilcoxon | Rebuilt |\n",
  "| Fig3F cell-cycle-stratified effects | Main | Frozen donor-aware phase models | Rebuilt |\n",
  "| Fig3G gene heatmap | Main candidate | Donor-specific matrix absent from frozen outputs | Blank audit placeholder; not recommended |\n",
  "| Purine-enzyme detection models | Supplementary FigS3 | Existing frozen panel `03_FIGURE_REBUILD_SINGLEPANELS/FIGURES/Fig3E_purine_enzyme_detection_models.*` | Retained unchanged |\n",
  "| Complete seven-module within-donor heatmap | Supplementary FigS3 | Existing frozen panel `03_FIGURE_REBUILD_SINGLEPANELS/FIGURES/Fig3D_within_donor_module_difference.*` | Retained unchanged |\n",
  "| Program multiverse robustness | Supplementary FigS3 | Existing locked multiverse output | Retained unchanged |\n",
  "| Hurdle positive-expression component | Supplementary FigS2 | Existing frozen sensitivity output | Retained unchanged |\n",
  "| Full Hallmark/Reactome/GO GSEA | Supplementary Figure/Table | `04_SUPPLEMENTARY_REBUILD` full outputs | Retained unchanged |\n\n",
  "The Results should continue to cite the Supplementary purine-enzyme comparison because it was requested by the reviewer.\n"
)
writeLines(manifest_text, file.path(report_dir, "FIG3_MAIN_VS_SUPPLEMENTARY_MANIFEST.md"), useBytes = TRUE)

legend_text <- paste0(
  "# Figure 3 legend core text\n\n",
  "(A) Frozen UMAP of UC ISC-like epithelial cells classified using corrected raw GART counts. GART-undetected and GART-detected cells were defined as raw counts of 0 and >=1, respectively.\n\n",
  "(B) SCP259 donor-paired pseudobulk differential expression comparing GART-detected with GART-undetected ISC-like epithelial cells across 11 paired donors. Positive log2 fold changes indicate higher expression in GART-detected cells. Colors indicate FDR<0.05 and direction. GART was omitted from the volcano display because GART detection defined the compared groups; its differential expression is not treated as independent validation.\n\n",
  "(C) Positive Hallmark enrichments passing FDR<0.05, ranked by normalized enrichment score; at most eight pathways were displayed.\n\n",
  "(D) Pre-specified Reactome and GO Biological Process replication, repair, cell-cycle, chromosome and nucleotide/purine pathways passing FDR<0.05. Complete enrichment results are provided in the Supplementary Table.\n\n",
  "(E) Donor-paired public-module scores for 23 donors with both GART-detected and GART-undetected ISC-like epithelial cells. Lines connect cells summarized within the same donor. All four modules used the same frozen paired-donor Wilcoxon analysis; displayed effects and 95% confidence intervals are archived mean detected-minus-undetected differences.\n\n",
  "(F) Frozen adjusted standardized module effects within G1, S and G2/M strata. All module-by-phase combinations are shown irrespective of significance. Models were donor-aware and adjusted for sequencing depth within each phase stratum.\n\n",
  "(G) The pre-specified proliferation and nucleotide-metabolism gene pool and aggregate pseudobulk FDR annotations are shown. Donor-specific normalized pseudobulk differences were not archived; cells are therefore left blank and no values were reconstructed.\n"
)
writeLines(legend_text, file.path(report_dir, "FIG3_LEGEND_CORE_TEXT.md"), useBytes = TRUE)

results_text <- paste0(
  "# Figure 3 results core sentences\n\n",
  "Donor-paired pseudobulk analysis of SCP259 identified higher expression of multiple proliferation and replication-associated genes in GART-detected ISC-like epithelial cells, including PCNA, MCM6, MCM7, FEN1, TK1, TYMS, GINS2, PAICS, MKI67 and CDK4 at FDR<0.05.\n\n",
  "Hallmark and Reactome enrichment analyses of the same paired-donor gene ranking showed significant positive enrichment of E2F targets, the G2/M checkpoint, MYC targets, DNA replication and DNA repair programs.\n\n",
  "Across 23 donors with both GART detection states, the frozen paired-donor analyses showed higher purine biosynthesis scores excluding GART, DNA replication scores, DNA repair scores and mitotic cell-cycle scores in GART-detected cells.\n\n",
  "Adjusted effects for these programs remained positive within G1, S and G2/M strata, although individual phase-specific uncertainty and FDR values are reported without significance-based filtering.\n\n",
  "Together, these findings indicate that GART detection marks a proliferative and nucleotide-metabolic epithelial state; they do not establish GART-dependent causality.\n\n",
  "For reviewer completeness, disease-association comparisons of PPAT, GART, PFAS, PAICS, ADSL and ATIC remain reported in the Supplementary Results and Tables.\n"
)
writeLines(results_text, file.path(report_dir, "FIG3_RESULTS_CORE_SENTENCES.md"), useBytes = TRUE)

cat(sprintf("Completed panels. Fig3A UC cells=%d donors=%d; Fig3B-D paired donors=11; Fig3E paired donors=23; Fig3G donor cells available=0/121.\n",
            nrow(uc), length(unique(uc$donor_id))))
