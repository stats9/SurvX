#!/usr/bin/env Rscript
# =====================================================================================================
#  Figures and tables of the simulation study, prepared for a BMC journal
#  (BMC Medical Research Methodology and similar)
#
#    Rscript sim_figures.R                 # after sim_published_models.R and/or sim_survx_models.R
#    Rscript sim_figures.R out_dir=results
#
#  Figure format (BMC "Preparing figures"): width 85 mm (one column) or 170 mm (two columns), at most
#  225 mm high; vector PDF plus TIFF at 600 dpi with LZW compression (PNG previews at 300 dpi);
#  sans-serif font (Arial / Helvetica) of 7-9 pt; no title inside the figure (the legend goes in the
#  manuscript); colours from the Okabe-Ito palette, which stays readable for colour-blind readers, with
#  different point shapes so that the figures also work in grey scale.
#
#  Written to figures/ (Fig*.pdf, Fig*.tiff, Fig*.png) and tables/ (Table*.csv, Table*.tex):
#    Fig1_scenarios          true CIF of both causes for four covariate profiles, eight scenarios
#    Fig2_shr_over_time      subdistribution hazard ratio of x1 over time (why Fine-Gray is biased)
#    Fig3_<study>            % reduction of the RMSE against Fine-Gray by sample size, middle horizon
#    FigS1_<study>           the same for all three horizons (additional file)
#    Fig4_<study>            per-replication RMSE ratio to Fine-Gray, n = 500, middle horizon
#    Fig5_<study>            bias of the predicted risk, n = 1000, every horizon
#    Fig6_published_vs_survx automatic choices of study 1 and study 2 side by side (if both ran)
#    Table1_scenarios        scenario definitions with observed event and censoring percentages
#    Table2_<study>          RMSE x 100 (Monte Carlo SE) and % reduction, middle horizon
#    TableS1_<study>         every summary measure, every cell (additional file)
# =====================================================================================================

out_dir <- "results"
for (a in commandArgs(trailingOnly = TRUE)) if (startsWith(a, "out_dir=")) out_dir <- sub("^out_dir=", "", a)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("the figures need the ggplot2 package: install.packages('ggplot2')")
library(ggplot2)

here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                if (length(f)) dirname(normalizePath(f)) else getwd() })
source(file.path(here, "sim_core.R"))
fig_dir <- "figures"; tab_dir <- "tables"
dir.create(fig_dir, showWarnings = FALSE); dir.create(tab_dir, showWarnings = FALSE)
SC <- scenarios()

# ---- style -------------------------------------------------------------------------------------------
okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
theme_bmc <- function(base = 8) {
  theme_bw(base_size = base, base_family = "sans") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          strip.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.3),
          strip.text = element_text(size = base, hjust = 0, margin = margin(2, 3, 2, 3)),
          panel.border = element_rect(colour = "grey60", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey40", linewidth = 0.3), axis.text = element_text(colour = "black", size = base - 1),
          legend.position = "bottom", legend.key.width = unit(6, "mm"), legend.margin = margin(0, 0, 0, 0),
          legend.text = element_text(size = base - 0.5), legend.title = element_blank(),
          plot.margin = margin(3, 4, 2, 2))
}
save_fig <- function(p, name, width_mm = 170, height_mm = 120) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  if (capabilities("cairo")) grDevices::cairo_pdf(file.path(fig_dir, paste0(name, ".pdf")), width = w, height = h)
  else grDevices::pdf(file.path(fig_dir, paste0(name, ".pdf")), width = w, height = h, useDingbats = FALSE)
  print(p); grDevices::dev.off()
  tt <- if (capabilities("cairo")) "cairo" else getOption("bitmapType")
  grDevices::tiff(file.path(fig_dir, paste0(name, ".tiff")), width = width_mm, height = height_mm, units = "mm",
                  res = 600, compression = "lzw", type = tt)
  print(p); grDevices::dev.off()
  grDevices::png(file.path(fig_dir, paste0(name, ".png")), width = width_mm, height = height_mm, units = "mm", res = 300, type = tt)
  print(p); grDevices::dev.off()
  cat("  ", name, sprintf("(%g x %g mm)\n", width_mm, height_mm))
}
scen_labels <- vapply(SC, `[[`, "", "label")
scen_factor <- function(id) factor(unname(scen_labels[id]), levels = unname(scen_labels))

# ---- LaTeX table writer (booktabs) -------------------------------------------------------------------
write_tex <- function(df, file, caption, label, align = NULL) {
  esc <- function(x) gsub("([%&_#])", "\\\\\\1", x)
  if (is.null(align)) align <- paste0("l", strrep("r", ncol(df) - 1L))
  body <- apply(df, 1L, function(r) paste(esc(as.character(r)), collapse = " & "))
  lines <- c("\\begin{table}[ht]", "\\centering", "\\footnotesize", paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
             paste0("\\begin{tabular}{", align, "}"), "\\toprule", paste(esc(names(df)), collapse = " & "), "\\\\ \\midrule",
             paste(body, "\\\\"), "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file)
}

# =====================================================================================================
# Figure 1: true cumulative incidence in each scenario
# =====================================================================================================
cat("Figures of the design\n")
profiles <- data.frame(x1 = c(0, 1, 0, 1), x2 = c(-0.5, -0.5, 0.5, 0.5))
profiles$label <- sprintf("x1 = %d, x2 = %+.1f", profiles$x1, profiles$x2)
truth_df <- do.call(rbind, lapply(SC, function(S) {
  tmax <- if (S$kind == "direct") 1.5 * max(S$times) else 1.6 * max(S$times)
  tg <- seq(0, tmax, length.out = 200)
  do.call(rbind, lapply(1:2, function(k) do.call(rbind, lapply(seq_len(nrow(profiles)), function(i) {
    X <- as.matrix(profiles[rep(i, length(tg)), c("x1", "x2")])
    data.frame(scenario = S$id, cause = paste("cause", k), profile = profiles$label[i], time = tg,
               cif = vapply(seq_along(tg), function(j) true_cif(X[j, , drop = FALSE], tg[j], S, k), 0))
  }))))
}))
truth_df$scen <- scen_factor(truth_df$scenario)
hz <- do.call(rbind, lapply(SC, function(S) data.frame(scen = scen_factor(S$id), time = S$times)))
p1 <- ggplot(truth_df, aes(time, cif, colour = profile, linetype = cause)) +
  geom_vline(data = hz, aes(xintercept = time), colour = "grey75", linewidth = 0.3, linetype = "dotted") +
  geom_line(linewidth = 0.45) +
  facet_wrap(~ scen, ncol = 4, scales = "free_x", labeller = label_wrap_gen(width = 26)) +
  scale_colour_manual(values = okabe[1:4]) + scale_linetype_manual(values = c("solid", "22")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Time (years)", y = "True cumulative incidence") + theme_bmc() +
  guides(colour = guide_legend(nrow = 1, order = 1), linetype = guide_legend(nrow = 1, order = 2))
save_fig(p1, "Fig1_scenarios", 170, 115)

# =====================================================================================================
# Figure 2: subdistribution hazard ratio of x1 over time (x2 = 0)
# =====================================================================================================
shr_df <- do.call(rbind, lapply(SC, function(S) {
  tg <- seq(0.02, 1.5 * max(S$times), length.out = 150); h <- 1e-4
  lam <- function(x1) {
    X <- matrix(c(x1, 0), 1)
    F <- function(t) true_cif(X, t, S, 1L)
    vapply(tg, function(t) (log1p(-F(t - h)) - log1p(-F(t + h))) / (2 * h), 0)      # -d/dt log(1 - F)
  }
  data.frame(scen = scen_factor(S$id), time = tg / max(S$times), shr = lam(1) / lam(0))
}))
# S1 and S7 both have a constant ratio exp(0.5): S7 is dashed so that S1 stays visible
pal8 <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#999999")
p2 <- ggplot(shr_df, aes(time, shr, colour = scen, linetype = scen)) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) + geom_line(linewidth = 0.55) +
  scale_y_log10() + scale_colour_manual(values = pal8) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "solid", "solid", "22", "solid")) +
  labs(x = "Time / last prediction horizon", y = "Subdistribution hazard ratio of x1 (log scale)") +
  theme_bmc() + theme(legend.position = "right", legend.key.height = unit(3.2, "mm")) + guides(colour = guide_legend(ncol = 1, override.aes = list(linewidth = 0.8)))
save_fig(p2, "Fig2_shr_over_time", 170, 75)

# =====================================================================================================
# Results of each study
# =====================================================================================================
# the constrained models are shown in the additional figure only: their errors in S1-S2 (up to -275%)
# would compress the scale of the main figure
extra_models <- list(published = c("Constrained PH", "Constrained PO"), survx = character(0))
main_models <- list(
  published = c("Direct (AIC)", "Direct PH", "Direct PO", "Direct GOR", "Cause-specific Cox"),
  survx = c("Direct (AIC)", "Weighted (AIC)", "Constrained PO, G by AIC", "Direct spline PH", "SurvX (AIC)", "Cause-specific Cox"))
shapes <- c(16, 17, 15, 18, 4, 1, 3, 8)

study_figures <- function(study) {
  f_raw <- file.path(out_dir, paste0("raw_", study, ".csv.gz"))
  if (!file.exists(f_raw)) { cat("no results for study '", study, "' (", f_raw, ")\n", sep = ""); return(NULL) }
  cat("Study:", study, "\n")
  raw <- utils::read.csv(f_raw, stringsAsFactors = FALSE)
  summ <- summarise_study(raw)
  utils::write.csv(summ, file.path(tab_dir, paste0("TableS1_", study, ".csv")), row.names = FALSE)
  mods <- intersect(main_models[[study]], unique(summ$model))
  mid <- do.call(rbind, lapply(split(summ, summ$scenario), function(s) { tm <- sort(unique(s$time)); s[s$time == tm[min(2L, length(tm))], ] }))
  mods_all <- c(mods, intersect(extra_models[[study]], unique(summ$model)))
  col_m <- stats::setNames(okabe[seq_along(mods_all)], mods_all); shp_m <- stats::setNames(shapes[seq_along(mods_all)], mods_all)

  # Figure 3: % reduction of the RMSE against Fine-Gray, middle horizon
  d3 <- mid[mid$model %in% mods, ]
  d3$model <- factor(d3$model, levels = mods); d3$scen <- scen_factor(d3$scenario)
  d3$lo <- 100 * (1 - d3$rel_upper); d3$hi <- 100 * (1 - d3$rel_lower)
  pd <- position_dodge(width = 0.12)
  p3 <- ggplot(d3, aes(n, reduction_pct, colour = model, shape = model)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) +
    (if (length(unique(d3$n)) > 1L) geom_line(linewidth = 0.4, position = pd)) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.3, position = pd) +
    geom_point(size = 1.4, position = pd) +
    facet_wrap(~ scen, ncol = 4, labeller = label_wrap_gen(width = 26)) +
    scale_x_log10(breaks = sort(unique(d3$n))) + scale_colour_manual(values = col_m) + scale_shape_manual(values = shp_m) +
    coord_cartesian(ylim = c(max(-60, min(d3$lo, 0, na.rm = TRUE)), min(60, max(d3$hi, 0, na.rm = TRUE)))) +
    labs(x = "Training sample size (log scale)", y = "RMSE reduction against Fine-Gray (%)") +
    theme_bmc() + guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  save_fig(p3, paste0("Fig3_", study), 170, 125)

  # Figure S1: every horizon
  dS <- summ[summ$model %in% mods_all, ]
  dS$model <- factor(dS$model, levels = mods_all); dS$scen <- scen_factor(dS$scenario)
  dS$horizon <- ave(dS$time, dS$scenario, FUN = function(t) match(t, sort(unique(t))))
  dS$horizon <- factor(paste("horizon", dS$horizon))
  pS <- ggplot(dS, aes(n, reduction_pct, colour = model, shape = model)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) + (if (length(unique(dS$n)) > 1L) geom_line(linewidth = 0.35)) + geom_point(size = 1) +
    facet_grid(scen ~ horizon, scales = "free_y", labeller = labeller(scen = function(x) sub(":.*", "", x))) +
    scale_x_log10(breaks = sort(unique(dS$n))) + scale_colour_manual(values = col_m) + scale_shape_manual(values = shp_m) +
    labs(x = "Training sample size (log scale)", y = "RMSE reduction against Fine-Gray (%)") +
    theme_bmc(7) + guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  save_fig(pS, paste0("FigS1_", study), 170, 225)

  # Figure 4: per-replication RMSE ratio, n = 500 (or the size closest to 500), middle horizon
  n4 <- unique(raw$n)[which.min(abs(unique(raw$n) - 500))]
  tm_mid <- vapply(split(raw$time, raw$scenario), function(t) sort(unique(t))[min(2L, length(unique(t)))], 0)
  r4 <- raw[raw$n == n4 & raw$time == tm_mid[raw$scenario] & raw$ok, ]
  ref <- r4[r4$model == "Fine-Gray", ]
  r4$ratio <- sqrt(r4$mse / ref$mse[match(paste(r4$scenario, r4$rep), paste(ref$scenario, ref$rep))])
  r4 <- r4[r4$model %in% mods & is.finite(r4$ratio), ]
  r4$model <- factor(r4$model, levels = rev(mods)); r4$scen <- scen_factor(r4$scenario)
  p4 <- ggplot(r4, aes(ratio, model, fill = model)) +
    geom_vline(xintercept = 1, colour = "grey40", linewidth = 0.35) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, linewidth = 0.3, width = 0.65) +
    facet_wrap(~ scen, ncol = 2, labeller = label_wrap_gen(width = 40)) + scale_x_log10() + scale_fill_manual(values = col_m, guide = "none") +
    labs(x = sprintf("RMSE of the model / RMSE of Fine-Gray in each replication (n = %d, log scale)", n4), y = NULL) +
    theme_bmc() + theme(panel.grid.major.y = element_blank())
  save_fig(p4, paste0("Fig4_", study), 170, 200)

  # Figure 5: bias, largest n, every horizon (Fine-Gray included)
  n5 <- max(summ$n)
  d5 <- summ[summ$n == n5 & summ$model %in% c(mods, "Fine-Gray"), ]
  m5 <- c(mods, "Fine-Gray"); d5$model <- factor(d5$model, levels = m5); d5$scen <- scen_factor(d5$scenario)
  p5 <- ggplot(d5, aes(time, 100 * bias, colour = model, shape = model)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) + geom_line(linewidth = 0.4) + geom_point(size = 1.3) +
    facet_wrap(~ scen, ncol = 4, scales = "free_x", labeller = label_wrap_gen(width = 26)) +
    scale_colour_manual(values = c(col_m, "Fine-Gray" = "grey35")) + scale_shape_manual(values = c(shp_m, "Fine-Gray" = 5)) +
    labs(x = "Prediction horizon (years)", y = sprintf("Mean bias of the predicted risk x 100 (n = %d)", n5)) +
    theme_bmc() + guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  save_fig(p5, paste0("Fig5_", study), 170, 120)

  # Table 2: RMSE x 100 (MCSE) and % reduction, middle horizon, every n
  t2 <- mid[mid$model %in% c("Fine-Gray", mods), ]
  t2$cell <- ifelse(t2$model == "Fine-Gray", sprintf("%.2f (%.2f)", 100 * t2$rmse, 100 * t2$rmse_mcse),
                    sprintf("%.2f (%.2f) %+.0f%%", 100 * t2$rmse, 100 * t2$rmse_mcse, -t2$reduction_pct))
  w <- stats::reshape(t2[, c("scenario", "n", "model", "cell")], idvar = c("scenario", "n"), timevar = "model", direction = "wide")
  names(w) <- sub("^cell\\.", "", names(w))
  w <- w[order(w$scenario, w$n), c("scenario", "n", intersect(c("Fine-Gray", mods), names(w)))]
  w$horizon <- tm_mid[w$scenario]; w <- w[, c("scenario", "horizon", setdiff(names(w), c("scenario", "horizon")))]
  utils::write.csv(w, file.path(tab_dir, paste0("Table2_", study, ".csv")), row.names = FALSE)
  write_tex(w, file.path(tab_dir, paste0("Table2_", study, ".tex")),
            caption = paste("RMSE $\\times$ 100 (Monte Carlo SE) of the predicted cumulative incidence of cause 1 at the middle",
                            "horizon, and the change against Fine-Gray. Negative percentages: smaller error than Fine-Gray."),
            label = paste0("tab:rmse_", study))
  invisible(list(raw = raw, summary = summ))
}

res <- list(published = study_figures("published"), survx = study_figures("survx"))

# =====================================================================================================
# Figure 6: the automatic choices of the two studies side by side
# =====================================================================================================
if (!is.null(res$published) && !is.null(res$survx)) {
  s1 <- res$published$summary; s2 <- res$survx$summary
  d6 <- rbind(s1[s1$model %in% c("Direct (AIC)", "Cause-specific Cox"), ], s2[s2$model %in% c("Weighted (AIC)", "SurvX (AIC)"), ])
  d6 <- do.call(rbind, lapply(split(d6, d6$scenario), function(s) { tm <- sort(unique(s$time)); s[s$time == tm[min(2L, length(tm))], ] }))
  lv <- c("Direct (AIC)", "Weighted (AIC)", "SurvX (AIC)", "Cause-specific Cox")
  d6$model <- factor(d6$model, levels = lv); d6$scen <- scen_factor(d6$scenario)
  p6 <- ggplot(d6, aes(n, reduction_pct, colour = model, shape = model)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) + geom_line(linewidth = 0.4) + geom_point(size = 1.4) +
    facet_wrap(~ scen, ncol = 4, labeller = label_wrap_gen(width = 26)) + scale_x_log10(breaks = sort(unique(d6$n))) +
    scale_colour_manual(values = stats::setNames(okabe[c(1, 3, 2, 6)], lv)) + scale_shape_manual(values = stats::setNames(c(16, 15, 17, 1), lv)) +
    labs(x = "Training sample size (log scale)", y = "RMSE reduction against Fine-Gray (%)") + theme_bmc()
  save_fig(p6, "Fig6_published_vs_survx", 170, 110)
}

# =====================================================================================================
# Table 1: scenarios
# =====================================================================================================
des <- NULL
for (st in c("published", "survx")) { f <- file.path(out_dir, paste0("design_", st, ".csv")); if (file.exists(f)) { des <- utils::read.csv(f); break } }
if (is.null(des) && !is.null(res$published)) des <- summarise_design(res$published$raw)
t1 <- do.call(rbind, lapply(SC, function(S) {
  truth <- switch(S$kind, direct = sprintf("direct Gompertz, alpha = (%s)", paste(S$alpha, collapse = ", ")),
                  fg = "Fine & Gray (1999) design", csh = "exponential cause-specific hazards")
  eff <- sprintf("(%s); (%s)", paste(S$beta[[1]], collapse = ", "), paste(S$beta[[2]], collapse = ", "))
  row <- data.frame(Scenario = S$id, Truth = truth, `beta cause 1; cause 2` = eff, Censoring = sprintf("U(0, %g)", S$cmax),
                    Horizons = paste(S$times, collapse = ", "), check.names = FALSE)
  if (!is.null(des)) {
    dd <- des[des$scenario == S$id & des$n == max(des$n), ]
    if (nrow(dd)) {
      nn <- dd$n[1]
      row$`Censored %` <- sprintf("%.0f", dd$censored_pct); row$`Cause 1 %` <- sprintf("%.0f", 100 * dd$events1 / nn)
      row$`Cause 2 %` <- sprintf("%.0f", 100 * dd$events2 / nn)
    }
  }
  row
}))
utils::write.csv(t1, file.path(tab_dir, "Table1_scenarios.csv"), row.names = FALSE)
write_tex(t1, file.path(tab_dir, "Table1_scenarios.tex"),
          caption = "Simulation scenarios. x1 \\textasciitilde{} Bernoulli(0.5), x2 \\textasciitilde{} Uniform(-1, 1); percentages observed in the simulated training data.",
          label = "tab:scenarios", align = paste0("ll", strrep("l", ncol(t1) - 2L)))
writeLines(unlist(lapply(SC, function(S) c(paste0(S$label), paste("  Role:", S$role), paste("  Parameters:", S$why), ""))),
           file.path(tab_dir, "scenario_descriptions.txt"))
cat("\nfigures in", normalizePath(fig_dir), "\ntables  in", normalizePath(tab_dir), "\n")
