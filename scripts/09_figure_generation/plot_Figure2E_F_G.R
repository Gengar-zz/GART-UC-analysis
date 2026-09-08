# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 13_FIG2_EFG_STYLE_REBUILD/CODE/plot_Fig2_EFG_style_rebuild.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors = FALSE)
set.seed(20260905)

revision_root <- file.path(paths$results_root, "singlecell")
donor_balanced_root <- file.path(paths$results_root, "donor_balanced_primary_validation")
out_root <- file.path(paths$results_root, "figures", "Figure2")
fig_dir <- file.path(out_root, "FIGURES")
src_dir <- file.path(out_root, "SOURCE_DATA")
report_dir <- file.path(out_root, "REPORTS")
qa_dir <- file.path(out_root, "QA")
code_dir <- file.path(out_root, "CODE")
invisible(lapply(c(fig_dir, src_dir, report_dir, qa_dir, code_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

COL <- list(
  ink = "#1B1B1B", muted = "#6F7782", wine = "#8F1537", rose = "#D69AA8",
  blue = "#6FA7CF", deepblue = "#356D93", midgrey = "#8A9098", lightgrey = "#D2D4D6",
  pale = "#F5F6F7", teal = "#3E7E76"
)

save_pair <- function(plot, stem, width_cm, height_cm) {
  grDevices::tiff(file.path(fig_dir, paste0(stem, ".tiff")), width = width_cm, height = height_cm,
                  units = "cm", res = 300, compression = "lzw", bg = "white",
                  type = if (.Platform$OS.type == "windows") "windows" else "cairo",
                  antialias = if (.Platform$OS.type == "windows") "cleartype" else "default")
  print(plot); dev.off()
  grDevices::cairo_pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = width_cm / 2.54,
                       height = height_cm / 2.54, family = "Arial", bg = "white", onefile = TRUE)
  print(plot); dev.off()
}

fmt_p <- function(x) {
  ifelse(x < 0.0001, formatC(x, format = "e", digits = 2), sprintf("%.4f", x))
}

# -----------------------------------------------------------------------------
# Fig2E: frozen feature values on the frozen ISC-like embedding
# -----------------------------------------------------------------------------
isc_path <- file.path(revision_root, "INTERMEDIATE", "ISC_UMAP_plot_data.csv")
isc <- read.csv(isc_path, check.names = FALSE)
stopifnot(nrow(isc) == 9306L)
stopifnot(all(c("UMAP1", "UMAP2", "ISC_score", "TA_score", "GART_expr", "purine_score",
                "ISC_subcluster", "condition", "GART_status") %in% names(isc)))

expand_range <- function(x, f = 0.025) {
  r <- range(x, finite = TRUE)
  r + c(-1, 1) * diff(r) * f
}
isc_xlim <- expand_range(isc$UMAP1)
isc_ylim <- expand_range(isc$UMAP2)
feature_palette <- c("#E6EBF0", "#F7F3EA", "#E2A17E", "#A9344B", "#67001F")

feature_plot <- function(data, value_col, title, ticks = TRUE) {
  d <- data[order(data[[value_col]], na.last = TRUE), ]
  p <- ggplot(d, aes(x = UMAP1, y = UMAP2, colour = .data[[value_col]])) +
    geom_point(size = 0.80, alpha = 0.90) +
    scale_colour_gradientn(colours = feature_palette, name = NULL,
                           limits = range(data[[value_col]], finite = TRUE),
                           breaks = pretty(range(data[[value_col]], finite = TRUE), n = 3)) +
    coord_fixed(xlim = isc_xlim, ylim = isc_ylim, expand = FALSE) +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    theme_classic(base_family = "Arial", base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10.5, hjust = 0.5, margin = margin(b = 4)),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7.8, colour = "black"),
      axis.line = element_line(linewidth = 0.40, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.position = "right",
      legend.text = element_text(size = 7.2),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(4, 2, 4, 4)
    ) +
    guides(colour = guide_colourbar(barheight = unit(17, "mm"), barwidth = unit(2.2, "mm"), ticks = TRUE))
  if (!ticks) p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank())
  p
}

p_e1 <- feature_plot(isc, "ISC_score", "Stem/ISC score", TRUE)
p_e2 <- feature_plot(isc, "TA_score", "TA/proliferating score", TRUE)
p_e3 <- feature_plot(isc, "GART_expr", "GART expression", TRUE)
p_e <- plot_grid(p_e1, p_e2, p_e3, nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 1, 1))

p_en1 <- feature_plot(isc, "ISC_score", "Stem/ISC score", FALSE)
p_en2 <- feature_plot(isc, "TA_score", "TA/proliferating score", FALSE)
p_en3 <- feature_plot(isc, "GART_expr", "GART expression", FALSE)
p_e_noticks <- plot_grid(p_en1, p_en2, p_en3, nrow = 1, align = "h", axis = "tb", rel_widths = c(1, 1, 1))

save_pair(p_e, "Fig2E_state_scores_and_GART_expression_clean", 22.0, 7.5)
save_pair(p_e_noticks, "Fig2E_state_scores_and_GART_expression_clean_noticks", 22.0, 7.5)

# -----------------------------------------------------------------------------
# Fig2F: descriptive frozen-subcluster profiles, no inferential calculations
# -----------------------------------------------------------------------------
cluster_levels <- sort(unique(as.character(isc$ISC_subcluster)))
profile <- do.call(rbind, lapply(cluster_levels, function(cl) {
  d <- isc[as.character(isc$ISC_subcluster) == cl, ]
  data.frame(
    final_ISClike_cluster = cl,
    n_cells = nrow(d),
    StemISC_score_mean = mean(d$ISC_score, na.rm = TRUE),
    TA_score_mean = mean(d$TA_score, na.rm = TRUE),
    GART_mean = mean(d$GART_expr, na.rm = TRUE),
    GART_detected_fraction = mean(d$GART_status == "positive", na.rm = TRUE),
    Purine_ex_GART_mean = mean(d$purine_score, na.rm = TRUE),
    UC_fraction = mean(d$condition == "UC", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
stopifnot(nrow(profile) == length(cluster_levels), sum(profile$n_cells) == 9306L)
write.csv(profile, file.path(src_dir, "Fig2F_ISClike_subcluster_profiles_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

metric_map <- c(
  StemISC_score_mean = "Stem/ISC score",
  TA_score_mean = "TA/proliferating score",
  GART_mean = "GART mean",
  GART_detected_fraction = "GART-detected cells (%)",
  Purine_ex_GART_mean = "Purine biosynthesis\n(GART excluded)",
  UC_fraction = "UC cells (%)"
)
heat <- do.call(rbind, lapply(names(metric_map), function(metric) {
  v <- profile[[metric]]
  z <- as.numeric(scale(v))
  z[!is.finite(z)] <- 0
  data.frame(
    final_ISClike_cluster = profile$final_ISClike_cluster,
    metric = unname(metric_map[metric]),
    raw_value = v,
    scaled_mean = z,
    stringsAsFactors = FALSE
  )
}))
heat$metric <- factor(heat$metric, levels = unname(metric_map))
heat$final_ISClike_cluster <- factor(heat$final_ISClike_cluster, levels = rev(cluster_levels))
max_abs <- max(abs(heat$scaled_mean), na.rm = TRUE)
p_f <- ggplot(heat, aes(x = metric, y = final_ISClike_cluster, fill = scaled_mean)) +
  geom_tile(colour = "white", linewidth = 0.75) +
  scale_fill_gradient2(low = "#356D93", mid = "#FAFAFA", high = "#A9344B", midpoint = 0,
                       limits = c(-max_abs, max_abs), name = "Scaled mean") +
  labs(x = NULL, y = "Final ISC-like cluster", title = "Stem/ISC-like subcluster profiles") +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5, margin = margin(b = 7)),
    axis.title.y = element_text(size = 9.5, margin = margin(r = 6)),
    axis.text.x = element_text(size = 8.2, angle = 42, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y = element_text(size = 8.5, face = "bold", colour = "black"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.8),
    plot.margin = margin(6, 8, 7, 7)
  ) +
  guides(fill = guide_colourbar(barheight = unit(30, "mm"), barwidth = unit(3, "mm")))
save_pair(p_f, "Fig2F_ISClike_subcluster_profiles_simplified", 15.5, 10.5)

# -----------------------------------------------------------------------------
# Fig2G: frozen donor summaries and frozen method comparison
# -----------------------------------------------------------------------------
donor <- read.csv(file.path(revision_root, "RESULTS", "GART_donor_fraction.csv"), check.names = FALSE)
donor_tests <- read.csv(file.path(revision_root, "RESULTS", "GART_donor_fraction_tests.csv"), check.names = FALSE)
impl <- read.csv(file.path(donor_balanced_root, "RESULTS", "implementation_comparison.csv"), check.names = FALSE)
interaction <- read.csv(file.path(donor_balanced_root, "RESULTS", "condition_by_cohort_interaction.csv"), check.names = FALSE)
lodo <- read.csv(file.path(donor_balanced_root, "RESULTS", "leave_one_donor_out.csv"), check.names = FALSE)
lodo_summary <- read.csv(file.path(donor_balanced_root, "RESULTS", "leave_one_donor_out_summary.csv"), check.names = FALSE)
boot <- read.csv(file.path(donor_balanced_root, "RESULTS", "stratified_donor_bootstrap_2000.csv"), check.names = FALSE)
boot_summary <- read.csv(file.path(donor_balanced_root, "RESULTS", "stratified_donor_bootstrap_summary.csv"), check.names = FALSE)

stopifnot(nrow(donor) == 33L, nrow(lodo) == 33L, nrow(boot) == 2000L, nrow(interaction) == 1L)
write.csv(donor, file.path(src_dir, "Fig2G_donor_level_detection_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

method_order <- c("Archived equal-donor HC1", "Implementation 1", "Implementation 2", "Implementation 3", "Archived donor-random-intercept GLMM")
display_labels <- c(
  "Equal donor – asymptotic z",
  "Equal donor – finite-cluster GEE",
  "Equal donor – survey design-based",
  "Equal donor – CR2/Satterthwaite",
  "Donor-random-intercept GLMM\n(reviewer-requested reference)"
)
method_idx <- match(method_order, impl$implementation)
methods <- impl[method_idx, , drop = FALSE]
stopifnot(!anyNA(method_idx), nrow(methods) == 5L)
methods$display_label <- display_labels
methods$display_role <- c("Main-text display row", "Finite-cluster GEE", "Survey sensitivity", "CR2 sensitivity", "Reviewer-requested reference")
method_source <- methods[, c("display_label", "display_role", "implementation", "estimand", "method", "OR", "OR_ci_low", "OR_ci_high", "p_value", "n_cells", "n_clusters", "converged", "warnings")]
write.csv(method_source, file.path(src_dir, "Fig2G_method_comparison_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

donor$condition <- factor(donor$condition, levels = c("Healthy", "UC"))
donor$cohort <- factor(donor$cohort, levels = c("GSE214695", "SCP259"))
n_ann <- aggregate(donor_id ~ cohort + condition, donor, length)
names(n_ann)[3] <- "n_donors"
n_wide <- reshape(n_ann, idvar = "cohort", timevar = "condition", direction = "wide")
n_wide$label <- sprintf("Healthy n=%d; UC n=%d", n_wide$n_donors.Healthy, n_wide$n_donors.UC)
p_ann <- merge(donor_tests[, c("cohort", "P")], n_wide[, c("cohort", "label")], by = "cohort")
p_ann$cohort <- factor(p_ann$cohort, levels = levels(donor$cohort))

p_g_left <- ggplot(donor, aes(condition, detection_fraction, colour = condition, shape = condition)) +
  geom_boxplot(aes(fill = condition), width = 0.48, alpha = 0.18, outlier.shape = NA,
               colour = COL$ink, linewidth = 0.45) +
  geom_jitter(width = 0.085, height = 0, size = 1.8, alpha = 0.90, stroke = 0.45) +
  geom_text(data = p_ann, aes(x = 1.5, y = 1.115, label = paste0("P = ", fmt_p(P))),
            inherit.aes = FALSE, family = "Arial", size = 3.0) +
  geom_text(data = p_ann, aes(x = 1.5, y = 1.045, label = label),
            inherit.aes = FALSE, family = "Arial", size = 2.75, colour = COL$muted) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c(Healthy = "#6E2C3A", UC = "#D69AA8"), guide = "none") +
  scale_fill_manual(values = c(Healthy = "#6E2C3A", UC = "#D69AA8"), guide = "none") +
  scale_shape_manual(values = c(Healthy = 16, UC = 15), guide = "none") +
  coord_cartesian(ylim = c(0, 1.16), clip = "off") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.50", "0.75", "1.00"), expand = c(0, 0)) +
  labs(x = NULL, y = "GART-detected fraction in\nISC-like epithelial cells",
       caption = "Each point represents one independent donor.") +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9.5, margin = margin(b = 5)),
    axis.title.y = element_text(size = 9.2),
    axis.text = element_text(size = 8, colour = "black"),
    axis.line = element_line(linewidth = 0.42),
    axis.ticks = element_line(linewidth = 0.35),
    panel.spacing.x = unit(6, "mm"),
    plot.caption = element_text(size = 7.7, colour = COL$muted, hjust = 0.5, margin = margin(t = 5)),
    plot.margin = margin(8, 4, 4, 4)
  )

methods$display_label <- factor(methods$display_label, levels = rev(display_labels))
methods$point_colour <- c(COL$wine, COL$midgrey, COL$midgrey, COL$midgrey, COL$rose)
methods$value_text <- sprintf("OR %.2f (%.2f–%.2f), P=%s", methods$OR, methods$OR_ci_low, methods$OR_ci_high, fmt_p(methods$p_value))
p_g_right <- ggplot(methods, aes(y = display_label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = COL$muted, linewidth = 0.45) +
  geom_segment(aes(x = OR_ci_low, xend = OR_ci_high, yend = display_label), colour = COL$ink, linewidth = 0.65) +
  geom_point(aes(x = OR, fill = point_colour), shape = 21, size = 2.5, colour = COL$ink, stroke = 0.45) +
  geom_text(aes(x = 9.0, label = value_text), hjust = 0, family = "Arial", size = 2.45, colour = COL$ink) +
  scale_fill_identity() +
  scale_x_log10(limits = c(0.70, 250), breaks = c(1, 2, 4, 8), labels = c("1", "2", "4", "8")) +
  labs(x = "Odds ratio (log scale)", y = NULL, title = "Inference-method comparison",
       caption = paste0("Positive point estimates across donor-balanced methods;\n",
                        "uncertainty varied by inference method.\n",
                        "Condition × cohort interaction P = ", sprintf("%.4f", interaction$p_value), ".")) +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 9.8, hjust = 0.5, margin = margin(b = 7)),
    axis.title.x = element_text(size = 8.8),
    axis.text.x = element_text(size = 7.8, colour = "black"),
    axis.text.y = element_text(size = 7.4, colour = "black", lineheight = 0.92),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.caption = element_text(size = 7.0, colour = COL$muted, hjust = 0, lineheight = 1.05, margin = margin(t = 5)),
    plot.margin = margin(8, 4, 4, 4)
  )

g_body <- plot_grid(p_g_left, p_g_right, nrow = 1, rel_widths = c(1.25, 1.10), align = "h", axis = "tb")
p_g <- ggdraw() +
  draw_label("Donor-level association of UC with GART detection", x = 0.5, y = 0.985,
             hjust = 0.5, vjust = 1, fontfamily = "Arial", fontface = "bold", size = 12) +
  draw_plot(g_body, x = 0, y = 0, width = 1, height = 0.93)
save_pair(p_g, "Fig2G_donor_level_GART_detection_and_method_summary_clean", 24.0, 10.5)

# -----------------------------------------------------------------------------
# Supplementary FigS2B: archived method, bootstrap, and LODO sensitivity
# -----------------------------------------------------------------------------
methods_s <- methods
methods_s$display_label <- factor(as.character(methods_s$display_label), levels = rev(display_labels))
p_s_method <- p_g_right +
  scale_x_log10(limits = c(0.70, 100), breaks = c(1, 2, 4, 8), labels = c("1", "2", "4", "8")) +
  labs(title = "Archived inference-method comparison", caption = NULL) +
  theme(plot.margin = margin(5, 5, 4, 5))

p_s_boot <- ggplot(boot, aes(estimate_logOR)) +
  geom_histogram(bins = 42, fill = "#D7E4EC", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = 2, colour = COL$muted, linewidth = 0.45) +
  geom_vline(xintercept = methods$estimate[1], colour = COL$wine, linewidth = 0.8) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("2,000 donor-bootstrap replicates\nPercentile OR CI %.2f–%.2f; sign P=%.4f",
                           boot_summary$bootstrap_OR_ci_low, boot_summary$bootstrap_OR_ci_high, boot_summary$bootstrap_sign_p),
           hjust = -0.02, vjust = 1.15, family = "Arial", size = 2.7, colour = COL$muted) +
  coord_cartesian(xlim = c(-2, 4)) +
  labs(x = "Log odds ratio", y = "Replicates", title = "Cohort-stratified donor bootstrap",
       caption = sprintf("Central display −2 to 4; %d extreme replicates retained in source data.",
                         sum(boot$estimate_logOR < -2 | boot$estimate_logOR > 4))) +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9.5), axis.title = element_text(size = 8.5), axis.text = element_text(size = 7.5),
        plot.caption = element_text(size = 7.0, colour = COL$muted, hjust = 0))

lodo <- lodo[order(lodo$OR), ]
lodo$rank <- seq_len(nrow(lodo))
p_s_lodo <- ggplot(lodo, aes(OR, rank)) +
  geom_vline(xintercept = 1, linetype = 2, colour = COL$muted, linewidth = 0.45) +
  geom_segment(aes(x = 1, xend = OR, yend = rank), colour = COL$lightgrey, linewidth = 0.45) +
  geom_point(shape = 21, fill = COL$deepblue, colour = COL$ink, size = 1.45, stroke = 0.3) +
  scale_x_log10(limits = c(0.85, 3.0), breaks = c(1, 1.5, 2, 2.5), labels = c("1", "1.5", "2", "2.5")) +
  scale_y_continuous(breaks = NULL) +
  labs(x = "Odds ratio (log scale)", y = NULL, title = "Leave-one-donor-out",
       subtitle = sprintf("33/33 OR >1; range %.2f–%.2f", lodo_summary$OR_min, lodo_summary$OR_max)) +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9.5), plot.subtitle = element_text(size = 7.8, colour = COL$muted), axis.title = element_text(size = 8.5), axis.text = element_text(size = 7.5), axis.line.y = element_blank())

s_bottom <- plot_grid(p_s_boot, p_s_lodo, nrow = 1, rel_widths = c(1.15, 0.85))
s_body <- plot_grid(p_s_method, s_bottom, ncol = 1, rel_heights = c(1.05, 0.95))
p_s <- ggdraw() +
  draw_label("Extended donor-balanced sensitivity", x = 0.5, y = 0.99, hjust = 0.5, vjust = 1,
             fontfamily = "Arial", fontface = "bold", size = 12) +
  draw_plot(s_body, x = 0, y = 0, width = 1, height = 0.95)
save_pair(p_s, "FigS2B_extended_donor_balanced_sensitivity", 18.0, 15.0)

supp_method <- data.frame(record_type = "method", label = as.character(methods$display_label),
                          OR = methods$OR, CI_low = methods$OR_ci_low, CI_high = methods$OR_ci_high,
                          P = methods$p_value, n_clusters = methods$n_clusters, stringsAsFactors = FALSE)
supp_lodo <- data.frame(record_type = "leave_one_donor_out", label = paste("Omit", lodo$omitted_donor),
                        OR = lodo$OR, CI_low = lodo$OR_ci_low, CI_high = lodo$OR_ci_high,
                        P = lodo$p_value, n_clusters = lodo$n_clusters, stringsAsFactors = FALSE)
supp_boot <- data.frame(record_type = "bootstrap_replicate", label = paste0("Bootstrap ", boot$iteration),
                        OR = boot$OR, CI_low = boot_summary$bootstrap_OR_ci_low,
                        CI_high = boot_summary$bootstrap_OR_ci_high, P = boot_summary$bootstrap_sign_p,
                        n_clusters = 33, stringsAsFactors = FALSE)
write.csv(rbind(supp_method, supp_lodo, supp_boot),
          file.path(src_dir, "FigS2B_extended_donor_balanced_sensitivity_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

qa_data <- data.frame(
  check = c("ISC_cells_Fig2E", "Fig2E_common_xmin", "Fig2E_common_xmax", "Fig2E_common_ymin", "Fig2E_common_ymax",
            "Fig2F_cluster_count", "Fig2F_cell_total", "Fig2G_donor_count", "Fig2G_method_count",
            "Fig2G_interaction_P", "LODO_rows", "bootstrap_rows"),
  value = c(nrow(isc), isc_xlim[1], isc_xlim[2], isc_ylim[1], isc_ylim[2], nrow(profile), sum(profile$n_cells),
            nrow(donor), nrow(methods), interaction$p_value, nrow(lodo), nrow(boot))
)
write.csv(qa_data, file.path(qa_dir, "data_integrity_audit.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Completed: ISC=%d; clusters=%d; donors=%d; methods=%d\n", nrow(isc), nrow(profile), nrow(donor), nrow(methods)))
