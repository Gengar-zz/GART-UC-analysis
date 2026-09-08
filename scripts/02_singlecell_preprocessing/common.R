# Public-release portability repair.
# Scientific code restored from authoritative final-analysis source: 02_SINGLECELL_REANALYSIS_REVISED/CODE/common.R
# Only path/input/output handling, missing-input messages, and clearly corrupted redaction strings were changed.
options(stringsAsFactors = FALSE)
set.seed(20260904)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this script with Rscript so the repository root can be resolved.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
script_dir <- dirname(script_file)
source(file.path(dirname(script_dir), "path_config.R"))
paths <- gart_paths(script_file)
project_root <- paths$repo_root
revision_root <- paths$work_root
singlecell_root <- file.path(paths$external_data_root, "singlecell", "frozen_objects")
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args) && nzchar(output_args[[1]])) output_args[[1]] else file.path(paths$results_root, "singlecell")
dirs <- c("FIGURES", "RESULTS", "AUDIT", "SOURCE_DATA", "CODE", "REPORTS", "INTERMEDIATE", "QA")
for (d in dirs) dir.create(file.path(output_root, d), recursive = TRUE, showWarnings = FALSE)
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

save_pair <- function(name, width = 7.09, height = 5.2, code) {
  expr <- substitute(code); env <- parent.frame()
  tf <- file.path(output_root, "FIGURES", paste0(name, ".tiff"))
  pf <- file.path(output_root, "FIGURES", paste0(name, ".pdf"))
  tiff(tf, width = width, height = height, units = "in", res = 300, compression = "lzw", bg = "white", type = "windows")
  par(family = "Arial"); eval(expr, env); dev.off()
  cairo_pdf(pf, width = width, height = height, family = "Arial", bg = "white")
  par(family = "Arial"); eval(expr, env); dev.off()
}

ptext <- function(p) ifelse(is.na(p), "P = NA", ifelse(p < 0.0001, sprintf("P = %.2e", p), sprintf("P = %.4f", p)))
zwithin <- function(x, g) ave(x, g, FUN = function(v) { z <- as.numeric(scale(v)); z[!is.finite(z)] <- 0; z })
read_corrected <- function() readRDS(gart_require_file(file.path(singlecell_root, "CORRECTED_SCRNA_ANALYSIS_OBJECT.rds"), "corrected single-cell analysis object"))

capture_fit <- function(expr) {
  warnings <- character()
  fit <- withCallingHandlers(tryCatch(eval.parent(substitute(expr)), error = identity), warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(fit = fit, warning = paste(unique(warnings), collapse = " | "))
}

fit_effect <- function(fit, term, analysis, structure, n_cells, n_donors, warning = "") {
  if (inherits(fit, "error")) return(data.frame(analysis, random_structure = structure, term, estimate = NA, SE = NA, CI_low = NA, CI_high = NA, P = NA, OR = NA, OR_CI_low = NA, OR_CI_high = NA, n_cells, n_donors, isSingular = NA, convergence_warning = conditionMessage(fit), donor_variance = NA, sample_variance = NA, residual_variance = NA))
  cs <- coef(summary(fit)); z <- cs[term, ]; est <- unname(z[1]); se <- unname(z[2]); p <- if (ncol(cs) >= 4) unname(z[4]) else 2 * pnorm(-abs(est / se))
  vc <- as.data.frame(lme4::VarCorr(fit)); dv <- sum(vc$vcov[vc$grp == "donor_id"], na.rm = TRUE); sv <- sum(vc$vcov[grepl("sample_id", vc$grp)], na.rm = TRUE); rv <- sum(vc$vcov[vc$grp == "Residual"], na.rm = TRUE)
  conv <- paste(c(warning, fit@optinfo$conv$lme4$messages), collapse = " | ")
  is_binomial <- inherits(fit, "glmerMod") && identical(stats::family(fit)$family, "binomial")
  data.frame(analysis, random_structure = structure, term, estimate = est, SE = se, CI_low = est - 1.96 * se, CI_high = est + 1.96 * se, P = p, OR = if (is_binomial) exp(est) else NA_real_, OR_CI_low = if (is_binomial) exp(est - 1.96 * se) else NA_real_, OR_CI_high = if (is_binomial) exp(est + 1.96 * se) else NA_real_, n_cells, n_donors, isSingular = lme4::isSingular(fit), convergence_warning = conv, donor_variance = dv, sample_variance = ifelse(structure == "donor/sample", sv, NA), residual_variance = rv)
}

forest_base <- function(d, labels, xlab, title, estimate = "estimate", lo = "CI_low", hi = "CI_high", fdr = NULL, color = "#8C001A") {
  y <- rev(seq_len(nrow(d))); xr <- range(c(d[[lo]], d[[hi]], 0), na.rm = TRUE); span <- diff(xr); if (!is.finite(span) || span == 0) span <- 1
  right_extra <- if (is.null(fdr)) .18 else .60
  par(mar = c(4.3, 13, 2.4, 2.6), mgp = c(2.4, .7, 0), tcl = -.25)
  plot(d[[estimate]], y, xlim = c(xr[1] - .14 * span, xr[2] + right_extra * span), ylim = c(.4, nrow(d) + .7), yaxt = "n", xlab = xlab, ylab = "", pch = 15, cex = 1.05, col = color, bty = "l", main = title)
  axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = .68); abline(v = 0, lty = 2, col = "#777777")
  segments(d[[lo]], y, d[[hi]], y, lwd = 1.4); points(d[[estimate]], y, pch = 15, cex = 1.0, col = color)
  if (!is.null(fdr)) text(max(d[[hi]], na.rm=TRUE) + .08 * span, y, sprintf("FDR=%.3g", d[[fdr]]), pos = 4, cex = .62)
}
