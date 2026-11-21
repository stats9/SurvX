# Paper 2: the figures. Run from the folder that holds paper2/results (after summarize.R):
#   Rscript paper2/figures/make_figures.R full      (or quick)
# Writes paper2/figures/fig*.pdf (vector graphics, Okabe-Ito colours).
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "full"
library(SurvX)
res_dir <- file.path("paper2", "results"); fig_dir <- file.path("paper2", "figures")
dir.create(fig_dir, showWarnings = FALSE)
rd <- function(f) utils::read.csv(file.path(res_dir, paste0(f, "_", mode, ".csv")), stringsAsFactors = FALSE)
TA <- rd("tableA"); TB <- rd("tableB"); TBc <- rd("tableB_coherence"); TC <- rd("tableC")
oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73", blue = "#0072B2", vermillion = "#D55E00",
        purple = "#CC79A7", grey = "#999999")
C_CON <- oi[["blue"]]; C_JNT <- oi[["vermillion"]]
CAUSE_COL <- unname(oi[c("vermillion", "blue", "green")])
causes <- c("Circulatory", "Neoplasms", "Other")
fig <- function(name, w, h, code) {
  grDevices::pdf(file.path(fig_dir, name), width = w, height = h, pointsize = 9)
  op <- par(mgp = c(2.2, 0.6, 0), tcl = -0.3, las = 1)
  ok <- tryCatch({ code; TRUE }, error = function(e) { message(name, " failed: ", conditionMessage(e)); FALSE })
  par(op); grDevices::dev.off()
  if (ok) cat("  ", name, "\n")
}
scen_order <- function(d) d[order(d$dgp != "constrained", -d$cmax, d$n), ]

# --------------------------------------------------------------------------------------------
# Figure 1: share of subjects whose predicted long-run risks add up to more than one
# --------------------------------------------------------------------------------------------
fig("fig1_coherence.pdf", 7, 3.4, {
  par(mar = c(4.6, 3.6, 2.4, 0.5))
  j <- scen_order(TBc[TBc$estimator == "joint, unconstrained", ])
  k <- TBc[TBc$estimator == "constrained, K = 3", ]
  k <- k[match(j$scenario, k$scenario), ]
  m <- 100 * rbind(j$excess_inf, k$excess_inf)
  allzero <- all(m[2, ] == 0)   # constrained model: zero by construction, so no bars of height 0
  if (allzero) m <- m[1, , drop = FALSE]
  b <- barplot(m, beside = TRUE, col = if (allzero) C_JNT else c(C_JNT, C_CON), border = NA, ylim = c(0, 110),
               axes = FALSE, names.arg = rep("", ncol(m)), space = if (allzero) c(0, 0.8) else c(0, 0.6))
  b <- matrix(b, nrow = nrow(m))
  axis(2, at = seq(0, 100, 25)); mtext("Long-run risks sum to > 1 (% of subjects)", 2, 2.4, las = 0, cex = 0.85)
  text(b[1, ], m[1, ] + 3, sprintf("%.0f", m[1, ]), cex = 0.75, col = C_JNT)
  if (!allzero) text(b[2, ], m[2, ] + 3, sprintf("%.0f", m[2, ]), cex = 0.75, col = C_CON)
  mid <- colMeans(b)
  mtext(paste0("n = ", j$n), 1, 0.4, at = mid, cex = 0.75)
  mtext(ifelse(j$cmax > 6, "long follow-up", "short follow-up"), 1, 1.4, at = mid, cex = 0.7, col = "grey30")
  g1 <- mid[j$dgp == "constrained"]; g2 <- mid[j$dgp != "constrained"]
  mtext("Data from the constrained model", 1, 2.8, at = mean(g1), cex = 0.8, font = 2)
  mtext("Data from cause-specific hazards", 1, 2.8, at = mean(g2), cex = 0.8, font = 2)
  abline(v = mean(c(max(g1), min(g2))), lty = 3, col = oi[["grey"]])
  if (allzero) {
    legend("top", c("joint model, unconstrained", "constrained model: 0% in every replication"),
           fill = c(C_JNT, NA), border = NA, col = c(NA, C_CON), pch = c(NA, 15), text.col = c("black", C_CON),
           bty = "n", horiz = TRUE, inset = c(0, -0.12), xpd = NA, cex = 0.9)
  } else {
    legend("top", c("joint model, unconstrained", "constrained model"), fill = c(C_JNT, C_CON), border = NA, bty = "n",
           horiz = TRUE, inset = c(0, -0.12), xpd = NA, cex = 0.9)
  }
})

# --------------------------------------------------------------------------------------------
# Figure 2: long-run risks F_k(inf | z = 0) in Study B: bias (left) and coverage (right),
# constrained against unconstrained joint model, for both data-generating models
# --------------------------------------------------------------------------------------------
fig("fig2_longrun.pdf", 7, 8, {
  layout(matrix(1:4, 2, 2, byrow = TRUE), widths = c(1.45, 1))
  for (g in c("constrained", "csh")) {
    d <- TB[TB$dgp == g & grepl("_inf_0$", TB$estimand) & TB$estimator %in% c("constrained, K = 3", "joint, unconstrained"), ]
    d$cause <- as.integer(substr(d$estimand, 2, 2))
    keys <- unique(d[order(d$cause, -d$cmax, d$n), c("cause", "cmax", "n")])
    gap <- (3 - keys$cause) * 0.8; y <- rev(seq_len(nrow(keys))) + gap
    pos <- function(est) { s <- d[d$estimator == est, ]; s[match(paste(keys$cause, keys$cmax, keys$n), paste(s$cause, s$cmax, s$n)), ] }
    con <- pos("constrained, K = 3"); jnt <- pos("joint, unconstrained")
    lab <- paste0("n = ", keys$n, ", ", ifelse(keys$cmax > 6, "long", "short"), " FU")
    ttl <- if (g == "constrained") "Data from the constrained model" else "Data from cause-specific hazards (model misspecified)"
    # bias
    par(mar = c(3.6, 7.6, 2.6, 0.4))
    lim <- max(abs(c(con$bias, jnt$bias)), 0.02, na.rm = TRUE) * 1.1
    plot(NA, xlim = c(-lim, lim), ylim = range(y) + c(-0.6, 0.6), yaxt = "n", xlab = "Bias", ylab = "")
    abline(v = 0, col = oi[["grey"]]); abline(h = y, col = "grey92")
    segments(con$bias, y + 0.15, jnt$bias, y - 0.15, col = "grey75")
    points(con$bias, y + 0.15, pch = 19, col = C_CON); points(jnt$bias, y - 0.15, pch = 17, col = C_JNT)
    axis(2, at = y, labels = lab, cex.axis = 0.8, tick = FALSE)
    for (c3 in 1:3) mtext(causes[c3], 2, 6.6, at = mean(y[keys$cause == c3]), las = 0, cex = 0.8, font = 2)
    mtext(ttl, 3, 0.8, adj = 0, font = 2, cex = 0.9)
    if (g == "constrained") legend("bottomright", c("constrained", "joint, unconstrained"), pch = c(19, 17),
                                    col = c(C_CON, C_JNT), bty = "n", cex = 0.85, bg = "white")
    # coverage
    par(mar = c(3.6, 0.6, 2.6, 0.8))
    plot(NA, xlim = c(0.2, 1), ylim = range(y) + c(-0.6, 0.6), yaxt = "n", xlab = "Coverage of 95% intervals", ylab = "")
    abline(v = 0.95, lty = 2, col = oi[["grey"]]); abline(h = y, col = "grey92")
    points(pmax(con$coverage, 0.2), y + 0.15, pch = 19, col = C_CON); points(pmax(jnt$coverage, 0.2), y - 0.15, pch = 17, col = C_JNT)
  }
})

# --------------------------------------------------------------------------------------------
# Figure 3: Study A, ratio of the mean standard error to the empirical SD, with and without the
# term for the estimated censoring distribution
# --------------------------------------------------------------------------------------------
fig("fig3_variance.pdf", 7, 3.6, {
  par(mfrow = c(1, 2), mar = c(3.6, 3.6, 2.2, 0.6))
  pairs <- list(c("weighted (IPCW variance)", "weighted (censoring ignored)", "PH and PO links"),
                c("weighted GOR (alpha estimated)", "weighted GOR (censoring ignored)", "GOR link, alpha estimated"))
  key <- function(x) paste(x$n, x$alpha, x$cmax, x$estimand)
  mm <- lapply(pairs, function(pr) {
    a <- TA[TA$estimator == pr[1], ]; b <- TA[TA$estimator == pr[2], ]
    merge(transform(a, k = key(a)), transform(b, k = key(b))[, c("k", "se_sd")], by = "k", suffixes = c("_ipcw", "_ign"))
  })
  allv <- unlist(lapply(mm, function(m) c(m$se_sd_ipcw, m$se_sd_ign)))
  rng <- if (length(allv)) range(c(allv, 0.95, 1.05), na.rm = TRUE) + c(-0.02, 0.02) else c(0.6, 1.45)
  for (i in seq_along(pairs)) {
    pr <- pairs[[i]]; m <- mm[[i]]
    plot(NA, xlim = rng, ylim = rng, xlab = "SE/SD, censoring treated as known", ylab = "SE/SD, with the censoring term", main = pr[3],
         cex.main = 0.95, font.main = 1)
    abline(0, 1, lty = 2, col = oi[["grey"]]); abline(h = 1, v = 1, col = "grey90")
    if (!nrow(m)) { text(1, 1, "no results", col = oi[["grey"]]); next }
    points(m$se_sd_ign, m$se_sd_ipcw, pch = ifelse(startsWith(m$estimand, "b_"), 1, 19),
           col = ifelse(m$cmax <= 4, C_JNT, C_CON), cex = 0.9)
    legend("topleft", c("coefficients", "incidences", "short follow-up", "long follow-up"), pch = c(1, 19, 15, 15),
           col = c("black", "black", C_JNT, C_CON), bty = "n", cex = 0.75)
  }
})

# --------------------------------------------------------------------------------------------
# Figure 4: computing time per fit (Study C)
# --------------------------------------------------------------------------------------------
fig("fig4_timing.pdf", 4.6, 3.6, {
  par(mar = c(3.6, 3.8, 0.8, 0.6))
  v <- c("weighted_PH", "weighted_PO", "joint_PH", "constrained_K3", "fine_gray")
  lab <- c("direct, weighted, PH", "direct, weighted, PO", "direct, joint, PH", "constrained, K = 3", "Fine-Gray")
  pal <- c(C_CON, oi[["sky"]], C_JNT, oi[["green"]], "black")
  vals <- unlist(TC[v]); vals <- vals[vals > 0]
  plot(NA, log = "xy", xlim = range(TC$n), ylim = range(vals) * c(0.7, 1.4), xlab = "Number of subjects", ylab = "Seconds per fit",
       xaxt = "n", yaxt = "n")
  axis(1, at = TC$n, labels = format(TC$n, big.mark = ","))
  yt <- 10^(-4:3); axis(2, at = yt, labels = format(yt, scientific = FALSE, drop0trailing = TRUE))
  for (i in seq_along(v)) { y <- TC[[v[i]]]; y[y <= 0] <- NA; lines(TC$n, y, type = "b", pch = 14 + i, col = pal[i], lwd = 1.5, lty = i) }
  legend("topleft", lab, col = pal, pch = 15:19, lty = 1:5, lwd = 1.5, bty = "n", cex = 0.75)
})

# --------------------------------------------------------------------------------------------
# Real data
# --------------------------------------------------------------------------------------------
data(flchain, package = "survival")
fl <- flchain
fl$time <- pmax(fl$futime, 1) / 365.25
ch <- as.character(fl$chapter)
fl$status <- ifelse(fl$death == 0, 0L, ifelse(ch %in% "Circulatory", 1L, ifelse(ch %in% "Neoplasms", 2L, 3L)))
fl$male <- as.integer(fl$sex == "M"); fl$flc_high <- as.integer(fl$flc.grp >= 9)
fit_group <- function(d, age0) {
  d$age5 <- (d$age - age0) / 5
  list(d = d, fit = directcif(Event(time, status) ~ age5 + male + flc_high, d, baseline = "logistic3",
                              estimation = "constrained", cause = 1, cond = "weibull"))
}
g80 <- tryCatch(fit_group(fl[fl$age >= 80, ], 85), error = function(e) NULL)
gall <- tryCatch(fit_group(fl, 65), error = function(e) NULL)
aj_step <- function(cv, k, col) {
  y <- cv[[paste0("cif.", k)]]
  lines(c(0, cv$time), c(0, y), type = "s", col = col, lwd = 1.1, lty = 2)
}
legend_lines <- function(where, labels, cols) {
  legend(where, c(labels, "model (average)", "Aalen-Johansen"), col = c(cols, "grey20", "grey20"),
         lty = c(rep(1, length(labels)), 1, 2), lwd = c(rep(2.5, length(labels)), 2, 1.1), bty = "n", cex = 0.8, seg.len = 2.5)
}

# Figure 5: average predicted CIF of each cause and the Aalen-Johansen estimate, to 30 years
fig("fig5_flchain_cif.pdf", 7, 3.6, {
  par(mfrow = c(1, 2), mar = c(3.6, 3.6, 2.2, 2.6))
  for (gg in list(list(g80, "Aged 80 or more at enrolment"), list(gall, "Whole cohort (aged 50 or more)"))) {
    g <- gg[[1]]; if (is.null(g)) next
    tg <- seq(0.05, 30, length.out = 200)
    plot(NA, xlim = c(0, 30), ylim = c(0, 1), xlab = "Years since enrolment", ylab = "Cumulative incidence",
         main = gg[[2]], cex.main = 0.95, font.main = 1)
    rect(max(g$d$time), -1, 40, 2, col = "grey95", border = NA); box()
    cv <- aj(Event(time, status) ~ 1, g$d)@curves[[1L]]
    ends <- numeric(3)
    for (k in 1:3) {
      m <- colMeans(predict(g$fit, g$d, times = tg, cause = k))
      ends[k] <- mean(predict(g$fit, g$d, times = Inf, cause = k))          # long-run risk (t -> infinity)
      lines(tg, m, col = CAUSE_COL[k], lwd = 2.5); aj_step(cv, k, CAUSE_COL[k])
    }
    yl <- pmin(ends, 0.95); o <- order(yl); for (i in 2:3) if (yl[o[i]] - yl[o[i - 1]] < 0.06) yl[o[i]] <- yl[o[i - 1]] + 0.06
    text(30.6, yl, sprintf("%.2f", ends), col = CAUSE_COL, cex = 0.8, adj = 0, xpd = NA)
    text(30.6, max(yl) + 0.08, expression(t %->% infinity), cex = 0.7, adj = 0, xpd = NA, col = "grey30")
    text(max(g$d$time) + 0.5, 0.97, "beyond follow-up", adj = c(0, 1), cex = 0.7, col = "grey40")
    legend_lines("topleft", causes, CAUSE_COL)
  }
})

# Figure 6: long-run risks by covariate profile, aged 80 or more (they add up to one)
fig("fig6_flchain_profiles.pdf", 7, 3.4, {
  if (is.null(g80)) stop("the model for the subjects aged 80 or more failed")
  par(mar = c(4.8, 3.6, 0.8, 7.2))
  nd <- expand.grid(age5 = c(0, 1), male = 0:1, flc_high = 0:1)
  lr <- sapply(1:3, function(k) predict(g80$fit, nd, times = Inf, cause = k))
  b <- barplot(t(lr), col = CAUSE_COL, border = "white", names.arg = rep("", nrow(nd)), ylab = "Long-run risk of death by cause",
               ylim = c(0, 1))
  mtext(ifelse(nd$age5 == 0, "85 y", "90 y"), 1, 0.3, at = b, cex = 0.75)
  mtext(ifelse(nd$male == 1, "male", "female"), 1, 1.2, at = b, cex = 0.75)
  mtext(ifelse(nd$flc_high == 1, "FLC high", "FLC normal"), 1, 2.1, at = b, cex = 0.75)
  cum <- t(apply(lr, 1, cumsum)); mids <- cum - lr / 2
  for (k in 1:3) text(b, mids[, k], sprintf("%.2f", lr[, k]), col = "white", cex = 0.7)
  legend(par("usr")[2], 1, rev(causes), fill = rev(CAUSE_COL), border = NA, bty = "n", cex = 0.85, xpd = NA)
})

# Figure 7: fertility, weighted direct model of each cause by education and the Aalen-Johansen estimate
fig("fig7_fertility.pdf", 7, 3.6, {
  data(fertility, package = "SurvX")
  par(mfrow = c(1, 2), mar = c(3.6, 3.6, 2.2, 0.6))
  lev <- levels(fertility$education); pal <- CAUSE_COL
  spec <- list(list(1, "mix_loglogistic", "GOR", "Live birth", 1), list(2, "mix_loglogistic", "PO", "Stillbirth or abortion", 0.12))
  for (s in spec) {
    f <- directcif(Event(time, status) ~ education, fertility, baseline = s[[2]], link = s[[3]],
                   estimation = "weighted", cause = s[[1]])
    tg <- seq(0.01, max(fertility$time), length.out = 200)
    pr <- predict(f, data.frame(education = factor(lev, levels = lev)), times = tg, cause = s[[1]])
    plot(NA, xlim = c(0, max(tg)), ylim = c(0, s[[5]]), xlab = "Years since marriage", ylab = "Cumulative incidence",
         main = s[[4]], cex.main = 0.95, font.main = 1)
    for (i in seq_along(lev)) {
      lines(tg, pr[i, ], col = pal[i], lwd = 2.5)
      d1 <- fertility[fertility$education == lev[i], ]
      aj_step(aj(Event(time, status) ~ 1, d1)@curves[[1L]], s[[1]], pal[i])
    }
    legend_lines(if (s[[1]] == 1) "bottomright" else "topleft", lev, pal)
  }
})
cat("figures written to", fig_dir, "\n")
