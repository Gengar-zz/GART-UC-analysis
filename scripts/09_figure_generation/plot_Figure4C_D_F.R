# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 17_FIG4_CDF_FINAL_REFINEMENT/CODE/plot_fig4_cdf_final_refinement.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
options(stringsAsFactors = FALSE)

out_root <- file.path(paths$results_root, "figures", "Figure4")
fig_dir <- out_root
src_dir <- file.path(paths$data_public_root, "figure_source_data", "Figure4")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (.Platform$OS.type == "windows") suppressWarnings(windowsFonts(Arial = windowsFont("Arial")))

burgundy <- "#A50024"
undetected <- "#6F7B86"
line_grey <- "#C3C3C3"
weak <- "#BDBDBD"
soft_weak <- "#C4A5AB"
dark <- "#222222"
grid <- "#E2E2E2"

save_pair <- function(stem, width, height, code) {
  expr <- substitute(code)
  env <- parent.frame()
  tiff(file.path(fig_dir, paste0(stem, ".tiff")), width = width, height = height,
       units = "in", res = 300, compression = "lzw", bg = "white",
       type = if (.Platform$OS.type == "windows") "windows" else "cairo", pointsize = 12)
  par(family = "Arial", fg = dark, col.axis = dark, col.lab = dark, col.main = dark)
  eval(expr, env)
  dev.off()
  cairo_pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = width, height = height,
            family = "Arial", bg = "white", pointsize = 12)
  par(family = "Arial", fg = dark, col.axis = dark, col.lab = dark, col.main = dark)
  eval(expr, env)
  dev.off()
}

fdr_label <- function(x) {
  ifelse(x < 1e-4, sprintf("FDR = %.1e", x), sprintf("FDR = %.4f", x))
}

# Figure 4C: frozen paired donor values; only the prespecified Wilcoxon result is displayed.
cdat <- read.csv(file.path(src_dir, "Fig4C_source_data.csv"), check.names = FALSE)
stopifnot(nrow(cdat) == 23L)
stopifnot(length(unique(cdat$displayed_P)) == 1L)
stopifnot(length(unique(cdat$parametric_sensitivity_P)) == 1L)

save_pair("Fig4C_donor_level_pseudotime_Wilcoxon_final", 7.2, 5.0, {
  par(mar = c(4.8, 5.2, 4.8, 2.2), mgp = c(2.8, .8, 0), tcl = -.25)
  yr <- range(c(cdat$GART_undetected_median, cdat$GART_detected_median))
  pad <- diff(yr) * 0.035
  plot(c(.94, 2.06), c(yr[1] - pad, yr[2] + pad), type = "n", xaxt = "n", xlab = "",
       ylab = "Donor median state-continuum position", bty = "l")
  axis(1, at = 1:2, labels = c("GART-undetected", "GART-detected"), cex.axis = .85)
  for (i in seq_len(nrow(cdat))) {
    lines(c(1, 2), c(cdat$GART_undetected_median[i], cdat$GART_detected_median[i]),
          col = line_grey, lwd = .9)
  }
  points(rep(1, nrow(cdat)), cdat$GART_undetected_median, pch = 16, cex = 1.16, col = undetected)
  points(rep(2, nrow(cdat)), cdat$GART_detected_median, pch = 16, cex = 1.16, col = burgundy)
  mtext("Donor-level pseudotime comparison", side = 3, line = 3.1, font = 2, cex = 1.12)
  mtext("n = 23 paired donors", side = 3, line = 1.7, cex = .82)
  mtext(sprintf("Paired Wilcoxon P = %.4f", unique(cdat$displayed_P)),
        side = 3, line = .45, cex = .88, font = 2, col = burgundy)
})

# Figure 4D: frozen correlations with compact labels and visual hierarchy.
ddat <- read.csv(file.path(src_dir, "Fig4D_source_data.csv"), check.names = FALSE)
ddat <- ddat[order(ddat$display_order), ]
stopifnot(nrow(ddat) == 2L)

save_pair("Fig4D_feature_exclusion_sensitivity_compact", 6.8, 3.25, {
  par(mar = c(4.3, 12.0, 2.9, 1.6), mgp = c(2.5, .7, 0), tcl = -.25)
  y <- c(2, 1)
  cols <- c(burgundy, soft_weak)
  plot(ddat$Spearman_rho, y, type = "n", xlim = c(-0.012, 0.245), ylim = c(.55, 2.45),
       yaxt = "n", xlab = "Spearman correlation with GART expression", ylab = "", bty = "l")
  abline(v = 0, col = "#777777", lwd = .8)
  abline(v = seq(.05, .20, .05), col = grid, lwd = .55)
  axis(2, at = y, labels = ddat$label, las = 1, tick = FALSE, cex.axis = .79)
  segments(0, y, ddat$Spearman_rho, y, col = cols, lwd = 2.4)
  points(ddat$Spearman_rho, y, pch = 16, cex = 1.15, col = cols)
  text(ddat$Spearman_rho + .009, y, sprintf("rho = %.3f", ddat$Spearman_rho),
       adj = c(0, .5), cex = .80, col = dark)
  title(main = "Feature-exclusion sensitivity", font.main = 2, cex.main = 1.10)
})

# Figure 4F: only the four prespecified main-text programs; no P/FDR recalculation.
fcore <- read.csv(file.path(src_dir, "Fig4F_core_programs_source_data.csv"), check.names = FALSE)
fcore <- fcore[order(fcore$display_order, decreasing = TRUE), ]
stopifnot(nrow(fcore) == 4L)

save_pair("Fig4F_core_programs_virtual_perturbation_final", 7.6, 4.6, {
  par(mar = c(4.8, 12.6, 3.4, 2.6), mgp = c(2.6, .8, 0), tcl = -.25)
  xmax <- max(fcore$minus_log10_FDR) * 1.34
  mids <- barplot(fcore$minus_log10_FDR, names.arg = fcore$label, horiz = TRUE, las = 1,
                  xlim = c(0, xmax), col = burgundy, border = NA,
                  xlab = "-log10 FDR of consensus-rank enrichment", cex.names = .80)
  abline(v = -log10(.05), lty = 2, lwd = 1, col = "#777777")
  text(fcore$minus_log10_FDR + .06, mids, fdr_label(fcore$FDR), pos = 4, cex = .72)
  title(main = "Core programs recurrently prioritized after virtual GART perturbation",
        font.main = 2, cex.main = .96)
})

# Supplementary Figure S4F: retain every archived program, including weak and nonsignificant results.
ffull_file <- gart_require_file(file.path(paths$external_data_root, "derived_inputs", "figure4", "FigS4F_full_programs_source_data.csv"), "archived supplementary Figure S4F source table")
ffull <- read.csv(ffull_file, check.names = FALSE)
ffull <- ffull[order(ffull$display_order, decreasing = TRUE), ]
stopifnot(nrow(ffull) == 7L)

save_pair("FigS4F_full_virtual_perturbation_programs", 8.0, 5.7, {
  par(mar = c(4.8, 14.8, 3.6, 2.8), mgp = c(2.6, .8, 0), tcl = -.25)
  xmax <- max(ffull$minus_log10_FDR) * 1.38
  cols <- ifelse(ffull$FDR < .05, burgundy, weak)
  mids <- barplot(ffull$minus_log10_FDR, names.arg = ffull$label, horiz = TRUE, las = 1,
                  xlim = c(0, xmax), col = cols, border = NA,
                  xlab = "-log10 FDR of consensus-rank enrichment", cex.names = .77)
  abline(v = -log10(.05), lty = 2, lwd = 1, col = "#777777")
  text(ffull$minus_log10_FDR + .06, mids, fdr_label(ffull$FDR), pos = 4, cex = .69)
  title(main = "Exploratory virtual perturbation: full program summary",
        font.main = 2, cex.main = 1.02)
})

message("Figure 4 C/D/F refinement panels written to: ", fig_dir)
