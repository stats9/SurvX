# Part 5. Report: summary tables (CSV and LaTeX for the paper), figures (PDF) and a readable
# report (Markdown, and HTML when rmarkdown is installed) with an overall verdict.

latex_table <- function(df, file, caption, label, digits = 3) {
  esc <- function(x) gsub("([&%_#$])", "\\\\\\1", x)
  fm <- function(v) if (is.numeric(v)) fmt_num(v, digits) else esc(as.character(v))
  body <- apply(as.data.frame(lapply(df, fm), stringsAsFactors = FALSE), 1, paste, collapse = " & ")
  al <- paste(ifelse(vapply(df, is.numeric, logical(1)), "r", "l"), collapse = "")
  writeLines(c("\\begin{table}[t!]", "\\centering", "\\small", paste0("\\begin{tabular}{", al, "}"), "\\hline",
               paste(esc(names(df)), collapse = " & "), "\\\\ \\hline", paste0(body, " \\\\"), "\\hline", "\\end{tabular}",
               paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"), "\\end{table}"), file)
}

md_table <- function(df, digits = 3) {
  fm <- function(v) if (is.numeric(v)) fmt_num(v, digits) else as.character(v)
  d <- as.data.frame(lapply(df, fm), stringsAsFactors = FALSE)
  c(paste0("| ", paste(names(d), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

# dot chart over as many pages as needed (at most 40 rows per page), with the title across the page
dot_pages <- function(x, labels, xlab, title, col, xlim, vlines = numeric(0), vlty = 2, log10x = FALSE, per_page = 40L) {
  n <- length(x); pages <- split(seq_len(n), ceiling(seq_len(n) / per_page))
  lm <- min(32, max(8, 0.5 * max(nchar(labels))))
  for (pg in seq_along(pages)) {
    i <- rev(pages[[pg]])
    op <- graphics::par(mar = c(4.5, lm, 1, 1), oma = c(0, 0, 3, 0))
    graphics::plot(x[i], seq_along(i), yaxt = "n", xlab = xlab, ylab = "", pch = 16, col = col[i], xlim = xlim,
                   ylim = c(0.5, per_page + 0.5), cex = 1.1)
    graphics::axis(2, seq_along(i), labels[i], las = 1, cex.axis = 0.7, tick = FALSE)
    graphics::abline(h = seq_along(i), col = "grey92", lwd = 0.5)
    graphics::points(x[i], seq_along(i), pch = 16, col = col[i], cex = 1.1)
    if (length(vlines)) graphics::abline(v = vlines, lty = vlty, col = "grey50")
    graphics::mtext(paste0(title, if (length(pages) > 1L) paste0("   (page ", pg, " of ", length(pages), ")") else ""),
                    side = 3, outer = TRUE, adj = 0.01, line = 1, font = 2, cex = 1.05)
    graphics::par(op)
  }
}

make_report <- function(cfg, out_dir) {
  logmsg("Part 5: report")
  tab_dir <- file.path(out_dir, "tables"); fig_dir <- file.path(out_dir, "figures")
  dir.create(tab_dir, showWarnings = FALSE); dir.create(fig_dir, showWarnings = FALSE)
  rd <- function(f) { p <- file.path(out_dir, f); if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL }
  acc <- rd("accuracy.csv"); spd <- rd("speed.csv"); sim <- rd("simulation_summary.csv"); stb <- rd("stability.csv")
  header <- c("# SurvX validation report", "",
              paste0("SurvX ", format(utils::packageVersion("SurvX")), " | R ", format(getRversion()), " | ", R.version$platform,
                     " | mode: ", cfg$mode, " | ", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
              paste0("Reference packages: ", paste(vapply(REF_PKGS, function(p) if (has(p)) paste0(p, " ", format(utils::packageVersion(p)))
                                                         else paste(p, "(missing)"), character(1)), collapse = ", ")), "")
  md <- character(0)
  txt <- character(0)                         # the plain-text results file (everything numeric)
  tx <- function(title, df = NULL, lines = NULL) {
    txt <<- c(txt, "", strrep("=", 100), title, strrep("=", 100))
    if (!is.null(lines)) txt <<- c(txt, lines)
    if (!is.null(df) && nrow(df)) {
      old <- options(width = 250); on.exit(options(old))
      txt <<- c(txt, utils::capture.output(print(df, row.names = FALSE, right = FALSE)))
    }
  }
  pdf_file <- file.path(out_dir, "SurvX_validation_figures.pdf")
  grDevices::pdf(pdf_file, width = 11, height = 8.5)
  pdf_open <- TRUE
  on.exit(if (pdf_open) grDevices::dev.off(), add = TRUE)
  verdict <- character(0)

  # ---- accuracy ----
  if (!is.null(acc)) {
    acc$group <- paste(acc$func, "|", acc$reference, "|", acc$quantity)
    s <- do.call(rbind, lapply(split(acc, acc$group), function(a) data.frame(
      area = a$area[1], func = a$func[1], reference = a$reference[1], quantity = a$quantity[1], comparisons = nrow(a),
      max_abs_diff = suppressWarnings(max(a$max_abs_diff, na.rm = TRUE)), max_rel_diff = suppressWarnings(max(a$max_rel_diff, na.rm = TRUE)),
      tolerance = a$tol[1], pass = sum(a$status == "PASS"), fail = sum(a$status == "FAIL"), error = sum(a$status == "ERROR"),
      info = sum(a$status == "INFO"), stringsAsFactors = FALSE)))
    s$max_abs_diff[!is.finite(s$max_abs_diff)] <- NA; s$max_rel_diff[!is.finite(s$max_rel_diff)] <- NA
    s <- s[order(s$area, s$func, s$reference), ]
    utils::write.csv(s, file.path(tab_dir, "accuracy_summary.csv"), row.names = FALSE)
    latex_table(s[, c("func", "reference", "quantity", "comparisons", "max_abs_diff", "max_rel_diff", "fail")], file.path(tab_dir, "accuracy_summary.tex"),
                "Agreement of SurvX with established implementations: largest absolute and relative difference over all data sets.", "tab:accuracy")
    md <- c(md, "## 1. Accuracy against established implementations", "",
            sprintf("%d comparisons: %d PASS, %d FAIL, %d ERROR, %d INFO (different valid estimators, not judged).",
                    nrow(acc), sum(acc$status == "PASS"), sum(acc$status == "FAIL"), sum(acc$status == "ERROR"), sum(acc$status == "INFO")), "",
            md_table(s[, c("func", "reference", "quantity", "comparisons", "max_abs_diff", "max_rel_diff", "tolerance", "fail", "error")]), "")
    tx("1. ACCURACY: summary by function, reference and quantity (largest difference over all data sets)", s,
       sprintf("%d comparisons: %d PASS, %d FAIL, %d ERROR, %d INFO", nrow(acc), sum(acc$status == "PASS"), sum(acc$status == "FAIL"),
               sum(acc$status == "ERROR"), sum(acc$status == "INFO")))
    bad <- acc[acc$status %in% c("FAIL", "ERROR"), c("func", "reference", "data", "quantity", "max_abs_diff", "max_rel_diff", "status", "note")]
    tx("1b. ACCURACY: every comparison that did not pass", bad)
    if (nrow(bad)) {
      md <- c(md, "### Comparisons that did not pass", "", md_table(utils::head(bad, 200)), "")
      verdict <- c(verdict, sprintf("Accuracy: %d comparisons FAIL or ERROR (see section 1).", nrow(bad)))
    }
    v <- log10(pmax(s$max_rel_diff, 1e-16))
    dot_pages(v, paste0(s$func, " vs ", sub("::.*", "", s$reference), ": ", s$quantity),
              "log10 of the largest relative difference from the reference (dashed: 1e-6 and 1e-3)",
              "Accuracy: agreement of SurvX with the reference implementation (red: FAIL or ERROR)",
              ifelse(s$fail + s$error > 0, "#D55E00", ifelse(s$tolerance == "info", "grey55", "#0072B2")), c(-16, 0), c(-6, -3))
  }

  # ---- simulation ----
  if (!is.null(sim)) {
    cols <- c("scenario", "method", "parameter", "n", "R", "converged_pct", "true", "bias", "rel_bias_pct", "emp_sd", "mean_se", "se_ratio",
              "coverage_pct", "rate_pct", "verdict")
    latex_table(sim[sim$judged, setdiff(cols, c("verdict", "rate_pct"))], file.path(tab_dir, "simulation.tex"),
                "Monte Carlo study: bias, empirical standard deviation (SD), mean standard error (SE), SE/SD and coverage of 95\\% intervals.", "tab:simulation")
    md <- c(md, "## 2. Statistical properties (Monte Carlo)", "",
            "Coverage is judged OK when within the Monte Carlo band 95% +/- 1.96 sqrt(0.95 x 0.05 / R), at least 99% of fits converged and SE/SD is within 1 +/- (1.96 sqrt(1 / (2 (R - 1))) + 0.05), its Monte Carlo error plus 5%; test size within 5% +/- its Monte Carlo band. 'reported' rows are shown for comparison and not judged.", "",
            md_table(sim[, cols]), "")
    tx("2. MONTE CARLO STUDY (bias, SD, SE, SE/SD, coverage %, test rejection %)", sim[, cols])
    chk <- sim[sim$verdict == "CHECK", ]
    if (nrow(chk)) verdict <- c(verdict, sprintf("Simulation: %d rows marked CHECK (see section 2).", nrow(chk)))
    js <- sim[!is.na(sim$coverage_pct), ]
    if (nrow(js)) {
      lab <- paste0(js$method, ": ", js$parameter, " (n=", js$n, ")")
      dot_pages(js$coverage_pct, lab, "coverage of the 95% interval (%)  -  dashed: Monte Carlo band",
                "Monte Carlo: coverage of 95% intervals (red: CHECK; grey: reported, not judged)",
                ifelse(js$verdict == "CHECK", "#D55E00", ifelse(js$judged, "#0072B2", "grey55")),
                c(min(75, js$coverage_pct, na.rm = TRUE), 100), 95 + c(-1, 1) * js$mc_band[1])
      dot_pages(js$se_ratio, lab, "mean model-based SE / empirical SD of the estimates  -  dashed: Monte Carlo band",
                "Monte Carlo: accuracy of the standard errors",
                ifelse(!js$judged, "grey55", ifelse(abs(js$se_ratio - 1) > js$se_ratio_band, "#D55E00", "#0072B2")),
                range(c(0.7, 1.3, js$se_ratio), na.rm = TRUE), c(1 - js$se_ratio_band[1], 1 + js$se_ratio_band[1]))
    }
  }

  # ---- speed ----
  if (!is.null(spd)) {
    ok <- spd[is.finite(spd$seconds) & spd$status == "ok", ]
    sn <- stats::aggregate(seconds ~ model + impl + n, ok[ok$design == "n", ], stats::median)
    wide <- stats::reshape(sn, idvar = c("model", "impl"), timevar = "n", direction = "wide")
    names(wide) <- sub("^seconds\\.", "n=", names(wide))
    wide <- wide[order(wide$model, wide$impl), ]
    utils::write.csv(wide, file.path(tab_dir, "speed_by_n.csv"), row.names = FALSE)
    latex_table(wide, file.path(tab_dir, "speed_by_n.tex"), "Median elapsed seconds of one fit (p = 5 covariates).", "tab:speed", digits = 2)
    skipped <- unique(spd[spd$status != "ok", c("model", "impl", "n", "status")])
    md <- c(md, "## 3. Speed", "", "Median elapsed seconds of one fit on the same data (p = 5).", "", md_table(wide, 2), "")
    tx("3. SPEED: median seconds of one fit, p = 5", wide)
    if (nrow(skipped)) tx("3b. SPEED: not run (too slow) or failed", skipped)
    if (nrow(skipped)) md <- c(md, "Not run or failed:", "", md_table(skipped), "")
    sp <- ok[ok$design == "p", ]
    if (nrow(sp)) {
      wp <- stats::reshape(stats::aggregate(seconds ~ model + impl + p, sp, stats::median), idvar = c("model", "impl"), timevar = "p", direction = "wide")
      names(wp) <- sub("^seconds\\.", "p=", names(wp))
      md <- c(md, "Effect of the number of covariates (n = 5000):", "", md_table(wp, 2), "")
      tx("3c. SPEED: median seconds against the number of covariates p (n = 5000)", wp)
    }
    models <- unique(sn$model)
    op <- graphics::par(mfrow = c(3, 3), mar = c(4, 4.5, 2.5, 1), oma = c(0, 0, 3, 0))
    for (m in models) {
      s <- sn[sn$model == m, ]; im <- unique(s$impl); cl <- grDevices::hcl.colors(max(2, length(im)), "Dark 3")
      graphics::plot(range(sn$n), range(pmax(s$seconds, 1e-5)), log = "xy", type = "n", xlab = "sample size n", ylab = "seconds per fit",
                     main = m, cex.main = 1)
      for (i in seq_along(im)) { z <- s[s$impl == im[i], ]; z <- z[order(z$n), ]
        graphics::lines(z$n, pmax(z$seconds, 1e-5), type = "b", pch = 16, col = cl[i], lwd = 2) }
      graphics::legend("topleft", im, col = cl[seq_along(im)], lwd = 2, bty = "n", cex = 0.75)
      graphics::mtext("Speed: median seconds of one fit against n (log scales; p = 5)", side = 3, outer = TRUE, adj = 0.01, line = 1, font = 2)
    }
    graphics::par(op)
  }

  # ---- stability ----
  if (!is.null(stb)) {
    md <- c(md, "## 4. Stability and invariances", "",
            sprintf("%d checks: %d PASS, %d FAIL.", nrow(stb), sum(stb$status == "PASS"), sum(stb$status == "FAIL")), "",
            md_table(stb[, c("kind", "case", "func", "outcome", "expected", "status", "message")]), "")
    tx("4. STABILITY AND INVARIANCES", stb[, c("kind", "case", "func", "outcome", "expected", "status", "message")])
    if (any(stb$status == "FAIL")) verdict <- c(verdict, sprintf("Stability: %d checks FAIL (see section 4).", sum(stb$status == "FAIL")))
    latex_table(stb[stb$kind == "invariance", c("case", "func", "outcome", "status")], file.path(tab_dir, "invariances.tex"),
                "Invariance checks: largest difference between equivalent fits.", "tab:invariance")
  }

  md <- c(header, "## Verdict", "", if (length(verdict)) paste("-", verdict) else "- All checks passed.", "", md)
  lg <- rec_get("log")
  tx("5. RUN LOG", lines = if (is.null(lg)) character(0) else paste(lg$time, lg$message))
  txt <- c(header[nzchar(header)], "", "VERDICT", if (length(verdict)) paste("-", verdict) else "- All checks passed.", txt,
           "", strrep("=", 100), "SESSION INFO", strrep("=", 100), utils::capture.output(sessionInfo()))
  writeLines(txt, file.path(out_dir, "SurvX_validation_results.txt"))
  grDevices::dev.off(); pdf_open <- FALSE
  writeLines(md, file.path(out_dir, "validation_report.md"))
  if (has("rmarkdown") && rmarkdown::pandoc_available())
    tryCatch(rmarkdown::render(file.path(out_dir, "validation_report.md"), output_format = rmarkdown::html_document(toc = TRUE),
                               quiet = TRUE), error = function(e) logmsg("  HTML report not made: ", conditionMessage(e)))
  logmsg("  results: ", normalizePath(file.path(out_dir, "SurvX_validation_results.txt")))
  logmsg("  figures: ", normalizePath(pdf_file))
  invisible(verdict)
}
