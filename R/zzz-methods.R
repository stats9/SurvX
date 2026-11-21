# Methods of the generics defined in zz-output.R for the parametric models (collated after it).

S7::method(estimates, ParSurvFit) <- function(object, level = 0.95, ...) .par_estimates(object, NA_integer_, level)
S7::method(estimates, ParCSCFit) <- function(object, level = 0.95, ...) {
  do.call(rbind, lapply(seq_len(object@K), function(k) .par_estimates(object@models[[k]], k, level)))
}


# plot() of a parametric survival model: fitted survival curves with the Kaplan-Meier-type
# (Nelson-Aalen) estimate as reference; hazard_plot() gives the hazard and cumulative hazard
S7::method(plot, ParSurvFit) <- function(x, ...) hazard_plot(x, type = "survival", ...)

# the tests of ph_test() as a data frame
S7::method(as.data.frame, PHTest) <- function(x, ...) x@table

# cure models (generic estimates() is defined in zz-output.R)
S7::method(estimates, CureFit) <- function(object, level = 0.95, ...) {
  est <- object@coefficients; nm <- names(est)
  pl <- ncol(object@X); L <- object@n_latency; P <- length(est)
  role <- c("baseline", rep("beta", pl), rep("shape", L - 1L - pl), "cure_baseline", rep("cure", P - L - 1L))
  sc <- ifelse(role == "beta", if (object@scale == "AFT") "log time ratio (uncured)" else "log hazard ratio (uncured)",
               ifelse(role == "cure", "log odds ratio of being cured", ifelse(role == "cure_baseline", "logit of the cured fraction", nm)))
  .est_frame(NA_integer_, nm, role, sc, est, .se_of(object), level, role %in% c("beta", "cure"))
}

# S7 >= 1.0.0: remove the placeholders that `method<-` leaves in the namespace for external generics
# (run when the package is built; older S7 versions do not have, or need, this hook)
if ("S7_on_build" %in% getNamespaceExports("S7")) getExportedValue("S7", "S7_on_build")()

# HTML report of adjusted_risk(); here because report() is defined in scan.R
S7::method(report, AdjustedRisk) <- function(x, file = "adjusted_risk_report.html", width = 8, height = 4.5, ...) {
  esc <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", s)))
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf, width = width, height = height, units = "in", res = 110)
  ok <- tryCatch({ plot(x); TRUE }, error = function(e) FALSE)
  grDevices::dev.off()
  img <- if (ok && file.exists(tf)) paste0("<img alt='adjusted risks' src='data:image/png;base64,", .b64(readBin(tf, "raw", file.info(tf)$size)), "'/>") else ""
  unlink(tf)
  r <- x@risk; cc <- x@contrast
  rt <- data.frame(method = r$method, time = format(r$time), level = r$level, risk = .num(r$risk, 3L),
                   lower = .num(r$lower, 3L), upper = .num(r$upper, 3L))
  ctab <- data.frame(method = cc$method, time = format(cc$time), contrast = paste0(cc$level, " vs ", cc$reference),
                     difference = .num(cc$difference, 3L), lower = .num(cc$lower, 3L), upper = .num(cc$upper, 3L),
                     ratio = .num(cc$ratio, 3L), `ratio lower` = .num(cc$ratio_lower, 3L), `ratio upper` = .num(cc$ratio_upper, 3L),
                     p = .pval(cc$p), check.names = FALSE)
  methods <- paste0("<p>",
    if (any(r$method == "gformula")) paste0("Standardisation (the parametric g-formula; Robins 1986): the fitted outcome model predicted, for every ",
      "subject, the risk under each treatment value with the other covariates as observed, and the predictions were averaged. ") else "",
    if (any(r$method == "ipw")) paste0("Inverse probability of treatment weighting (Cole and Hernan 2004): a ",
      if (length(x@levels) == 2L) "logistic" else "multinomial logistic", " model of the treatment (", esc(paste(deparse(x@ps_formula), collapse = " ")),
      ") gave stabilised weights", if (x@ps$trim > 0) paste0(", truncated at the ", 100 * x@ps$trim, "th and ", 100 * (1 - x@ps$trim), "th percentiles") else "",
      "; the weighted ", if (x@kind == "single") "Kaplan-Meier" else "Aalen-Johansen", " estimate was computed in each treatment group ",
      "(Neumann and Billionnet 2016 for competing risks). Balance was checked by standardised mean differences (Austin and Stuart 2015). ") else "",
    if (x@B > 0L) paste0("Percentile intervals come from ", x@ok, " bootstrap samples in which the models were refitted (Efron and Tibshirani 1993).") else "",
    "</p>")
  html <- c("<!DOCTYPE html><html><head><meta charset='utf-8'><title>SurvX adjusted risks</title>",
            "<style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;color:#222;line-height:1.5}",
            "table{border-collapse:collapse;font-size:13px;margin:1em 0}th,td{border-bottom:1px solid #ddd;padding:4px 8px;text-align:right}",
            "th:first-child,td:first-child{text-align:left}pre{background:#f6f6f6;padding:1em;overflow-x:auto}img{max-width:100%}",
            ".note{background:#fff8e1;padding:.6em 1em;border-left:4px solid #f0b400}</style></head><body>",
            "<h1>Adjusted absolute risks</h1>", paste0("<p>SurvX ", as.character(utils::packageVersion("SurvX")), ", ", format(Sys.time(), "%Y-%m-%d %H:%M"), "</p>"),
            "<h2>Summary</h2>", paste0("<pre>", esc(paste(.ar_lines(x, TRUE), collapse = "\n")), "</pre>"),
            "<h2>Figure</h2>", img,
            "<h2>Risks</h2>", .html_table(rt), "<h2>Contrasts</h2>", .html_table(ctab),
            if (!is.null(x@ps) && !is.null(x@ps$smd)) c("<h2>Balance</h2>", .html_table(data.frame(covariate = x@ps$smd$covariate, level = x@ps$smd$level,
                                                    before = .num(x@ps$smd$before, 3L), after = .num(x@ps$smd$after, 3L)))),
            "<h2>Methods</h2>", methods,
            "<div class='note'><p><b>Causal interpretation</b> requires no unmeasured confounding, positivity and a correct model; ",
            "with competing risks the contrast is a total effect on the cumulative incidence.</p></div>",
            "<h2>Code</h2>", paste0("<pre>", esc(paste(deparse(x@model_call), collapse = "\n")), "</pre>"),
            "<h2>References</h2>", paste0("<p>", esc(.scan_refs(.method_keys(x))), "</p>"), "</body></html>")
  writeLines(html, file, useBytes = TRUE)
  message("report written to ", normalizePath(file, mustWork = FALSE))
  invisible(file)
}
