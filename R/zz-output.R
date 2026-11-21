# Everything the user sees: print, summary, estimates, confint, plot and the usual extractors.
# Text is assembled as character vectors by the .render_* functions; print() shows the short form,
# summary() returns a FitSummary whose print() shows the long form.

# ----------------------------------------------------------------------------------------------
# formatting helpers
# ----------------------------------------------------------------------------------------------
.num <- function(x, d = 3L) {
  out <- formatC(as.numeric(x), format = "f", digits = d)
  out[!is.finite(as.numeric(x))] <- "-"
  out
}

.pval <- function(p) {
  out <- ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3L)))
  as.character(out)
}

.ci_text <- function(lo, hi, d = 2L) paste0(.num(lo, d), " to ", .num(hi, d))

.pct <- function(level) paste0(format(100 * level, trim = TRUE), "%")

# a data frame of strings as aligned lines: first column left-justified, the others right-justified
.table_lines <- function(df, indent = 2L) {
  if (!nrow(df)) return(character(0))
  cols <- lapply(seq_along(df), function(j) {
    v <- c(names(df)[j], as.character(df[[j]]))
    w <- max(nchar(v))
    formatC(v, width = w, flag = if (j == 1L) "-" else "")
  })
  paste0(strrep(" ", indent), do.call(paste, c(cols, sep = "  ")))
}

.zq <- function(level) stats::qnorm(1 - (1 - level) / 2)

# one row per coefficient, with Wald interval and test; expo: which rows get the exponentiated columns
.est_frame <- function(cause, term, role, scale, est, se, level, expo) {
  zq <- .zq(level)
  lo <- est - zq * se; hi <- est + zq * se
  expo <- rep_len(expo, length(est))
  data.frame(cause = cause, term = term, role = role, scale = scale,
             estimate = unname(est), se = unname(se), lower = unname(lo), upper = unname(hi),
             z = unname(est / se), p = unname(2 * stats::pnorm(-abs(est / se))),
             exp_estimate = ifelse(expo, exp(unname(est)), NA_real_),
             exp_lower = ifelse(expo, exp(unname(lo)), NA_real_),
             exp_upper = ifelse(expo, exp(unname(hi)), NA_real_),
             stringsAsFactors = FALSE, row.names = NULL)
}

# ratios: 2 decimals, 3 when the ratio is close to one (effects of continuous covariates per unit)
.ratio_digits <- function(r) ifelse(is.finite(r) & abs(log(pmax(r, 1e-300))) < 0.1, 3L, 2L)

.ratio_text <- function(r, lo, hi) {
  dg <- .ratio_digits(r)
  list(ratio = vapply(seq_along(r), function(i) .num(r[i], dg[i]), character(1)),
       ci = vapply(seq_along(r), function(i) .ci_text(lo[i], hi[i], dg[i]), character(1)))
}

# readable names of the columns of the design matrix: "sexM" becomes "sex: M vs F"
.pretty_terms <- function(term, xlevels) {
  out <- term
  for (v in names(xlevels)) {
    lev <- xlevels[[v]]
    if (length(lev) < 2L) next
    for (l in lev[-1L]) out[term == paste0(v, l)] <- paste0(v, ": ", l, " vs ", lev[1L])
  }
  out
}

# the table of covariate effects shown to the user: ratio, interval, z, p
.effect_lines <- function(e, label, level, xlevels = NULL) {
  rt <- .ratio_text(e$exp_estimate, e$exp_lower, e$exp_upper)
  df <- data.frame(term = .pretty_terms(e$term, xlevels), ratio = rt$ratio, ci = rt$ci,
                   z = .num(e$z, 2L), p = .pval(e$p), stringsAsFactors = FALSE)
  names(df) <- c("", label, paste(.pct(level), "CI"), "z", "p")
  .table_lines(df)
}

.se_of <- function(object) suppressWarnings(sqrt(diag(object@vcov)))

.status_lines <- function(converged, flags) {
  out <- character(0)
  ft <- setdiff(.flag_text(flags), c("sum of fitted CIFs exceeds 1 within follow-up", "sum of fitted CIFs exceeds 1 at infinity"))
  if (!isTRUE(converged))
    out <- c(out, "WARNING: the fit did NOT converge. Standard errors, intervals and tests should not be used.")
  if (length(ft)) out <- c(out, paste0("  reason: ", paste(ft, collapse = "; ")))
  out
}

# ----------------------------------------------------------------------------------------------
# estimates(): one tidy data frame for every model class
# ----------------------------------------------------------------------------------------------

#' Estimates as a data frame
#'
#' @description
#' All coefficients of a fitted model in one data frame, with standard errors, Wald confidence
#' limits and tests, and for regression coefficients the exponentiated estimate (hazard ratio,
#' subdistribution hazard ratio or odds ratio, see the column \code{scale}) with its limits.
#' @param object a model fitted by \code{\link{cox}}, \code{\link{fg}}, \code{\link{csc}} or \code{\link{directcif}}.
#' @param level confidence level.
#' @param ... not used.
#' @details
#' Standard errors are the square roots of the diagonal of the inverse observed information (for
#' \code{fg()} the robust sandwich variance by default, for weighted direct models the sandwich with
#' one score per subject). Limits are Wald intervals \eqn{\hat\beta \pm z_{1-(1-level)/2}\,SE} on
#' the estimation scale, and \code{exp_*} columns are their exponentials; \code{p} is the two-sided
#' Wald test of a zero coefficient.
#' @section Interpreting the table:
#' \code{estimate} and \code{se} are on the scale of estimation (log hazard ratio, log odds, log
#' shape, ...); \code{exp_estimate} with \code{exp_lower} and \code{exp_upper} is the effect on the
#' scale named in \code{scale} (hazard ratio, subdistribution hazard ratio, time ratio, odds ratio of
#' cure, ...), the one to report. Confidence limits and p-values are Wald limits and tests; for small
#' samples, compare nested models with a likelihood ratio test (\code{anova}-style: twice the
#' difference of \code{logLik()}).
#' @return
#' a data frame with the columns \code{cause, term, role, scale, estimate, se, lower, upper,
#'   z, p, exp_estimate, exp_lower, exp_upper}. \code{role} is \code{"beta"} for regression
#'   coefficients; for direct models also \code{"baseline"}, \code{"alpha"} and \code{"conditional"}
#'   (on the scale used for estimation; \code{summary()} shows them on the natural scale).
#' @references
#' Wald A (1943) Tests of statistical hypotheses concerning several parameters when the number of
#'   observations is large. Transactions of the American Mathematical Society 54:426-482.
#'
#' Casella G, Berger RL (2002) Statistical Inference, 2nd ed. Duxbury, Section 10.3.
#'
#' White H (1982) Maximum likelihood estimation of misspecified models. Econometrica 50:1-25.
#' @seealso
#' \code{\link{SurvX-methods}} (print, summary, coef, confint), \code{\link{sx_references}}.
#' @examples
#' set.seed(4)
#' d <- data.frame(time = rexp(200), status = rbinom(200, 1, 0.7), x = rnorm(200))
#' estimates(cox(Event(time, status) ~ x, d))
#' @export
estimates <- S7::new_generic("estimates", "object", function(object, level = 0.95, ...) S7::S7_dispatch())

S7::method(estimates, CoxFit) <- function(object, level = 0.95, ...) {
  .est_frame(NA_integer_, names(object@coefficients), "beta", "log hazard ratio", object@coefficients, .se_of(object), level, TRUE)
}

S7::method(estimates, FineGrayFit) <- function(object, level = 0.95, ...) {
  .est_frame(object@cause, names(object@coefficients), "beta", "log subdistribution hazard ratio",
             object@coefficients, .se_of(object), level, TRUE)
}

S7::method(estimates, CSCFit) <- function(object, level = 0.95, ...) {
  do.call(rbind, lapply(seq_len(object@K), function(k) {
    m <- object@models[[k]]
    .est_frame(k, names(m$coefficients), "beta", "log cause-specific hazard ratio", m$coefficients,
               suppressWarnings(sqrt(diag(m$vcov))), level, TRUE)
  }))
}

.link_scale <- function(link) {
  if (identical(link, "PH")) "log subdistribution hazard ratio"
  else if (identical(link, "PO")) "log odds ratio of the cumulative incidence"
  else "log effect on the generalized odds rate scale"
}

.link_of_cause <- function(object, k) if (object@estimation == "joint") object@link[k] else object@link[1L]

S7::method(estimates, DirectCIFFit) <- function(object, level = 0.95, ...) {
  ix <- object@index
  scale <- ifelse(ix$role == "beta", vapply(ix$cause, function(k) .link_scale(.link_of_cause(object, k)), character(1)),
                  ifelse(ix$role == "alpha", "alpha of the generalized odds rate link",
                         ifelse(ix$role == "share_beta", "odds ratio of the share of this cause", "estimation scale")))
  .est_frame(ix$cause, ix$term, ix$role, scale, object@coefficients, .se_of(object), level, ix$role %in% c("beta", "share_beta"))
}

# ----------------------------------------------------------------------------------------------
# extractors shared by the regression models
# ----------------------------------------------------------------------------------------------
S7::method(coef, SurvXFit) <- function(object, ...) object@coefficients

S7::method(vcov, SurvXFit) <- function(object, ...) object@vcov

S7::method(nobs, SurvXFit) <- function(object, ...) object@n

S7::method(logLik, SurvXFit) <- function(object, ...) {
  structure(object@loglik, df = length(object@coefficients), nobs = object@n, class = "logLik")
}

S7::method(formula, SurvXFit) <- function(x, ...) stats::formula(x@terms)

S7::method(formula, CSCFit) <- function(x, ...) stats::formula(x@terms)

S7::method(confint, SurvXFit) <- function(object, parm, level = 0.95, ...) {
  est <- object@coefficients; se <- .se_of(object)
  a <- (1 - level) / 2
  out <- cbind(est - .zq(level) * se, est + .zq(level) * se)
  dimnames(out) <- list(names(est), paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3L), "%"))
  if (missing(parm)) out else out[parm, , drop = FALSE]
}

S7::method(coef, CSCFit) <- function(object, ...) lapply(object@models, function(m) m$coefficients)

S7::method(vcov, CSCFit) <- function(object, ...) lapply(object@models, function(m) m$vcov)

S7::method(nobs, CSCFit) <- function(object, ...) object@n

S7::method(confint, CSCFit) <- function(object, parm, level = 0.95, ...) {
  a <- (1 - level) / 2
  lapply(object@models, function(m) {
    se <- suppressWarnings(sqrt(diag(m$vcov)))
    out <- cbind(m$coefficients - .zq(level) * se, m$coefficients + .zq(level) * se)
    dimnames(out) <- list(names(m$coefficients), paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3L), "%"))
    out
  })
}

# ----------------------------------------------------------------------------------------------
# Cox, Fine-Gray, cause-specific Cox: text
# ----------------------------------------------------------------------------------------------
# lines describing strata, (entry, time] data and the robust variance of Cox-type fits
.design_lines <- function(x) {
  out <- character(0)
  if (!is.null(x@strata))
    out <- c(out, paste0("  Strata: ", paste(x@strata_vars, collapse = ", "), " (", length(x@strata_levels),
                         " strata, each with its own baseline hazard)"))
  if (!is.null(x@entry))
    out <- c(out, "  Data: (entry, time] intervals (delayed entry or time-varying covariates; Andersen & Gill 1982)")
  if (isTRUE(x@robust))
    out <- c(out, if (!is.null(x@cluster))
      paste0("  Standard errors: robust, clustered on ", length(unique(x@cluster)), " clusters (Lin & Wei 1989)")
      else "  Standard errors: robust sandwich (Lin & Wei 1989)")
  out
}

.render_cox <- function(x, level, detail) {
  e <- estimates(x, level = level)
  lr <- 2 * (x@loglik - x@loglik_null); df <- length(x@coefficients)
  pv <- .pval(stats::pchisq(lr, df, lower.tail = FALSE))
  out <- c(paste0("Cox proportional hazards regression (", x@ties, " method for ties)"),
           paste0("  ", if (is.null(x@entry)) "Subjects" else "Rows", ": ", x@n, "   events: ", x@nevent, "   censored: ", x@n - x@nevent),
           .design_lines(x),
           "",
           "Effect of covariates: hazard ratio (HR)",
           .effect_lines(e, "HR", level, x@xlevels),
           "",
           .tt_lines(x, level),
           .reading_lines(e[!grepl("^tt\\(", e$term), , drop = FALSE], x@xlevels, "PH", "the event", level),
           "",
           sprintf("Fit: log partial likelihood %.2f | AIC %.2f | %d coefficients", x@loglik, -2 * x@loglik + 2 * df, df),
           sprintf("     likelihood ratio test against no covariates: chi-square %.2f on %d df, p %s", lr, df,
                   if (startsWith(pv, "<")) sub("<", "< ", pv) else paste("=", pv)))
  if (detail)
    out <- c(out, sprintf("     log partial likelihood without covariates: %.3f; iterations: %d", x@loglik_null, x@iterations), "",
             "HR > 1: the event occurs at a higher rate; HR < 1: at a lower rate.",
             "Numeric covariates: effect of a one-unit increase; categorical covariates: level against the reference level.",
             "Check proportional hazards with ph_test(fit); absolute risks with predict(fit, newdata, times, type = \"survival\").")
  c(out, .status_lines(x@converged, x@flags))
}

# time-varying coefficients: what the tt() rows mean and the HR at the quartiles of the event times
.tt_lines <- function(x, level) {
  if (is.null(x@tt)) return(character(0))
  te <- tt_effect(x, level = level)
  tab <- data.frame(covariate = te$term, time = format(signif(te$time, 4)), HR = sprintf("%.3f", te$hr),
                    ci = sprintf("%.3f to %.3f", te$lower, te$upper), stringsAsFactors = FALSE)
  names(tab)[4] <- paste0(round(100 * level), "% CI")
  c("Time-varying effects: tt(x) is the change of log HR of x per unit of g(t); HR of x at time t = exp(b_x + b_tt g(t)):",
    .table_lines(tab), "(tt_effect(fit, times) for other times)", "")
}

.render_fg <- function(x, level, detail) {
  e <- estimates(x, level = level)
  out <- c(paste0("Fine-Gray regression for the cumulative incidence of cause ", x@cause),
           paste0("  Subjects: ", x@n, "   events of cause ", x@cause, ": ", x@nevent,
                  "   other events: ", sum(x@status != 0L & x@status != x@cause), "   censored: ", sum(x@status == 0L)),
           if (!is.null(x@entry)) "  Delayed entry: censoring and truncation weights of Geskus (2011)",
           paste0("  Standard errors: ", if (x@variance == "model") "model-based"
                  else if (x@variance == "bootstrap") paste0("bootstrap of the subjects (", attr(x@vcov_boot, "B"), " samples)")
                  else if (is.null(x@entry)) "robust sandwich (Fine & Gray 1999)"
                  else "robust sandwich over subjects (Lin & Wei 1989), weights taken as known"),
           if (!is.null(x@entry) && x@variance == "robust")
             "  Note: with delayed entry these standard errors can be too small; variance = \"bootstrap\" accounts for the weights",
           "",
           "Effect of covariates: subdistribution hazard ratio (sHR)",
           .effect_lines(e, "sHR", level, x@xlevels), "",
           sub("^  ", "", .reading_cif(e, x@xlevels, "PH", paste("cause", x@cause), level)))
  if (detail)
    out <- c(out, "",
             sprintf("Log pseudo-likelihood: %.3f (null model: %.3f); iterations: %d", x@loglik, x@loglik_null, x@iterations),
             "sHR > 1: higher cumulative incidence of this cause at every time; sHR < 1: lower.",
             "The size of an sHR is not a ratio of risks; use predict() for absolute risks.",
             "Numeric covariates: effect of a one-unit increase; categorical covariates: level against the reference level.")
  c(out, .status_lines(x@converged, x@flags))
}

.render_csc <- function(x, level, detail) {
  e <- estimates(x, level = level)
  out <- c(paste0("Cause-specific Cox regression (", x@ties, " method for ties)"),
           paste0("  ", if (is.null(x@entry)) "Subjects" else "Rows", ": ", x@n, "   censored: ", sum(x@status == 0L)),
           .design_lines(x))
  for (k in seq_len(x@K)) {
    m <- x@models[[k]]
    ek <- e[e$cause == k, , drop = FALSE]
    out <- c(out, "", paste0("Cause ", k, " (", m$nevent, " events): cause-specific hazard ratio (csHR)"),
             .effect_lines(ek, "csHR", level, x@xlevels),
             .reading_lines(ek, x@xlevels, "PH", paste0("cause ", k, " among subjects still event-free"), level, "csHR"))
    pk <- length(m$coefficients); lr <- 2 * (m$loglik - m$loglik_null)
    pv <- .pval(stats::pchisq(lr, pk, lower.tail = FALSE))
    out <- c(out, sprintf("Fit: log partial likelihood %.2f | AIC %.2f | LR test %.2f on %d df, p %s", m$loglik,
                          -2 * m$loglik + 2 * pk, lr, pk, if (startsWith(pv, "<")) sub("<", "< ", pv) else paste("=", pv)),
             .status_lines(m$converged, m$flags))
  }
  if (detail)
    out <- c(out, "",
             "csHR: effect on the rate of one cause among subjects who are still event-free.",
             "A covariate can raise the rate of a cause and still lower its cumulative incidence,",
             "because the incidence also depends on the other causes; use predict() for absolute risks.",
             "Numeric covariates: effect of a one-unit increase; categorical covariates: level against the reference level.",
             "Check proportional hazards with ph_test(fit).")
  c(out, "", "Check: the cumulative incidences of all causes and the event-free probability add up to one (coherent by construction).")
}

# ----------------------------------------------------------------------------------------------
# direct models: text
# ----------------------------------------------------------------------------------------------
.base_label <- c(gompertz = "Gompertz", gompertz3 = "three-parameter Gompertz", weibull = "Weibull",
                 mix_weibull = "Weibull with plateau", mix_loglogistic = "log-logistic with plateau",
                 mix_lognormal = "log-normal with plateau", logistic3 = "modified logistic", spline = "restricted cubic spline")

.link_label <- function(link) {
  if (identical(link, "PH")) "proportional subdistribution hazards (PH)"
  else if (identical(link, "PO")) "proportional odds (PO)"
  else if (identical(link, "GOR")) "generalized odds rate, alpha estimated (GOR)"
  else paste0("generalized odds rate, ", link)
}

.effect_name <- function(link) {
  if (identical(link, "PH")) c("sHR", "subdistribution hazard ratio (sHR)")
  else if (identical(link, "PO")) c("OR", "odds ratio of the cumulative incidence (OR)")
  else c("exp(b)", "exp(coefficient) on the generalized odds rate scale")
}

# predicted CIF of cause k at typical covariate values (see .profile) with interval; times may contain Inf
.typical_risk <- function(x, k, times, level) {
  pr <- tryCatch(.direct_predict(x, .profile(x)$X, times, k, TRUE), error = function(e) NULL)
  if (is.null(pr)) return(NULL)
  Fv <- drop(pr$cif); se <- drop(pr$se)
  ci <- .cif_interval(Fv, se, .zq(level))
  data.frame(time = times, risk = Fv, se = se, lower = ci$lower, upper = ci$upper)
}

.has_plateau <- function(x, k) {
  if (x@estimation == "constrained") return(TRUE)
  b <- if (x@estimation == "joint") x@baseline[k] else x@baseline[1L]
  !b %in% c("weibull", "spline")
}

# baseline (and conditional) parameters on their natural scale
.natural_table <- function(e, level) {
  nat <- function(term, v) ifelse(grepl("^log_", term), exp(v), ifelse(grepl("^logit_", term), stats::plogis(v), v))
  name <- sub("^logit_", "", sub("^log_", "", e$term))
  data.frame(parameter = name, estimate = nat(e$term, e$estimate), lower = nat(e$term, e$lower), upper = nat(e$term, e$upper),
             stringsAsFactors = FALSE)
}

.render_direct <- function(x, level, detail) {
  e <- estimates(x, level = level)
  est_text <- switch(x@estimation,
    joint = "joint likelihood of all causes (Jeong & Fine 2006, 2007)",
    weighted = paste0("weighted likelihood for cause ", x@interest, " alone, censoring weights as in Fine & Gray (1999)"),
    constrained = if (x@K > 2L) paste0(x@K, " causes with long-run risks adding up to one (Shi, Cheng & Jeong 2013, ",
                                        "extended with multinomial shares of the competing causes)")
                  else "two causes with long-run risks adding up to one (Shi, Cheng & Jeong 2013)")
  cnt <- table(factor(x@status, levels = 0:x@K))
  out <- c("Direct parametric regression for the cumulative incidence function",
           paste0("  Estimation: ", est_text),
           paste0("  Subjects: ", x@n, "   censored: ", cnt[[1L]], "   ",
                  paste0("cause ", seq_len(x@K), ": ", as.integer(cnt[-1L]), collapse = "   ")))
  if (x@estimation == "weighted")
    out <- c(out, paste0("  Standard errors: ", switch(x@diagnostics$variance %||% "robust",
      ipcw = "sandwich with the term for the estimated censoring distribution",
      robust = "sandwich, censoring distribution treated as known",
      model = "inverse information of the pseudo-likelihood (censoring distribution treated as known)")))
  tables <- list()
  for (k in unique(x@index$cause)) {
    ek <- e[e$cause == k, , drop = FALSE]
    if (all(ek$role %in% c("conditional", "share", "share_beta"))) {
      oth <- setdiff(seq_len(x@K), x@interest)
      cj <- x@cond[min(match(k, oth), length(x@cond))]
      cj <- c(logistic = "modified logistic", weibull = "Weibull", loglogistic = "log-logistic", lognormal = "log-normal")[[cj]]
      if (x@K > 2L) {
        out <- c(out, "", paste0("Cause ", k, ": cumulative incidence = (1 - long-run risk of cause ", x@interest, ") x share(z) x G(t)"),
                 paste0("  G: ", cj, " distribution of the time to a cause ", k, " event; no covariates"))
        es <- ek[ek$role == "share_beta", , drop = FALSE]
        if (k == oth[1L]) out <- c(out, "  Share: reference cause of the multinomial shares")
        else if (nrow(es)) {
          es$term <- sub("^share_", "", es$term)
          out <- c(out, paste0("  Effect of covariates on the share (odds of cause ", k, " rather than cause ", oth[1L],
                               " among subjects who do not fail from cause ", x@interest, "): OR"),
                   .effect_lines(es, "OR", level, x@xlevels))
        }
      } else
        out <- c(out, "", paste0("Cause ", k, ": cumulative incidence = (1 - long-run risk of cause ", x@interest, ") x G(t)"),
                 paste0("  G: ", cj, " distribution of the time to a cause ", k, " event; no covariates"))
      ek <- ek[ek$role == "conditional", , drop = FALSE]
      pl <- .typical_risk(x, k, Inf, level)
      if (!is.null(pl)) out <- c(out, sprintf("  Long-run risk at typical covariate values: %s (%s CI %s)", .num(pl$risk, 3L), .pct(level), .ci_text(pl$lower, pl$upper, 3L)))
      if (detail) {
        nt <- .natural_table(ek, level)
        out <- c(out, "  Parameters of G", .table_lines(data.frame(` ` = nt$parameter, estimate = .num(nt$estimate, 4L),
                                                                   ci = .ci_text(nt$lower, nt$upper, 4L), check.names = FALSE), 4L))
        tables[[paste0("conditional.", k)]] <- nt
      }
      next
    }
    link <- .link_of_cause(x, k); en <- .effect_name(link)
    b <- if (x@estimation == "joint") x@baseline[k] else x@baseline[1L]
    out <- c(out, "", paste0("Cause ", k, ": baseline ", .base_label[[b]], ", link ", .link_label(link)))
    eb <- ek[ek$role == "beta", , drop = FALSE]
    if (nrow(eb)) out <- c(out, paste0("  Effect of covariates: ", en[2L]), .effect_lines(eb, en[1L], level, x@xlevels),
                           .reading_cif(eb, x@xlevels, link, paste("cause", k), level))
    ea <- ek[ek$role == "alpha", , drop = FALSE]
    if (nrow(ea))
      out <- c(out, sprintf("  alpha = %s (%s CI %s); alpha = 0 is the PH model, alpha = 1 the PO model",
                            .num(ea$estimate, 3L), .pct(level), .ci_text(ea$lower, ea$upper, 3L)))
    if (.has_plateau(x, k)) {
      pl <- .typical_risk(x, k, Inf, level)
      if (!is.null(pl) && is.finite(pl$risk)) {
        if (pl$risk < 0.9995)
          out <- c(out, sprintf("  Long-run risk (plateau) at typical covariate values: %s (%s CI %s)", .num(pl$risk, 3L), .pct(level), .ci_text(pl$lower, pl$upper, 3L)))
        else out <- c(out, "  Long-run risk at typical covariate values: no plateau below 1 (the cumulative incidence keeps rising)")
      }
    }
    if (detail) {
      ev <- x@time[x@status == k]
      if (length(ev) >= 4L) {
        tq <- signif(as.numeric(stats::quantile(ev, c(0.25, 0.5, 0.75))), 3L)
        tr <- .typical_risk(x, k, tq, level)
        if (!is.null(tr)) {
          out <- c(out, "  Predicted cumulative incidence at typical covariate values",
                   .table_lines(data.frame(time = format(tr$time, trim = TRUE), risk = .num(tr$risk, 3L), ci = .ci_text(tr$lower, tr$upper, 3L)), 4L))
          tables[[paste0("risk.", k)]] <- tr
        }
      }
      nt <- .natural_table(ek[ek$role == "baseline", , drop = FALSE], level)
      out <- c(out, "  Baseline parameters", .table_lines(data.frame(` ` = nt$parameter, estimate = .num(nt$estimate, 4L),
                                                                     ci = .ci_text(nt$lower, nt$upper, 4L), check.names = FALSE), 4L))
      tables[[paste0("baseline.", k)]] <- nt
    }
  }
  out <- c(out, "")
  if (ncol(x@X)) out <- c(out, paste0("Typical covariate values: ", .profile(x)$label, " (mean of numeric, most frequent level of categorical covariates)"))
  if (x@estimation == "weighted")
    out <- c(out, sprintf("Weighted log pseudo-likelihood: %.3f on %d parameters (compare weighted fits only)", x@loglik, length(x@coefficients)))
  else out <- c(out, sprintf("Fit: log-likelihood %.2f | AIC %.2f | BIC %.2f | %d parameters", x@loglik,
                             -2 * x@loglik + 2 * length(x@coefficients), -2 * x@loglik + log(x@n) * length(x@coefficients),
                             length(x@coefficients)))
  d <- x@diagnostics
  if (isTRUE(d$n_exceed_range > 0))
    out <- c(out, sprintf("NOTE: for %d of %d subjects the fitted risks of all causes add up to more than 1 within follow-up",
                          as.integer(d$n_exceed_range), x@n),
             sprintf("  (largest sum %.3f). The joint model assumes that some subjects never fail; consider", d$max_sum_range),
             "  estimation = \"weighted\" or \"constrained\".")
  out <- c(out, .status_lines(x@converged, x@flags))
  if (detail) {
    links <- unique(vapply(unique(x@index$cause[!x@index$role %in% c("conditional", "share", "share_beta")]), function(k) .link_of_cause(x, k), character(1)))
    out <- c(out, "", "How to read the effects")
    if ("PH" %in% links) out <- c(out, "  sHR > 1: higher cumulative incidence of this cause at every time; its size is not a ratio of risks.")
    if ("PO" %in% links) out <- c(out, "  OR: ratio of the odds F/(1 - F) of having had the event by time t; the same at every time.")
    if (any(!links %in% c("PH", "PO"))) out <- c(out, "  exp(b): effect on the generalized odds rate scale; between a hazard ratio (alpha = 0) and an odds ratio (alpha = 1).")
    out <- c(out, "  Numeric covariates: effect of a one-unit increase; categorical covariates: level against the reference level.",
             "  Absolute risks for a patient: predict(fit, newdata, times).")
  }
  list(text = out, tables = tables)
}

# ----------------------------------------------------------------------------------------------
# print and summary of the regression models
# ----------------------------------------------------------------------------------------------
.render <- function(x, level, detail) {
  r <- if (S7::S7_inherits(x, DirectCIFFit)) .render_direct(x, level, detail)
       else if (S7::S7_inherits(x, ParSurvFit)) list(text = .render_par(x, level, detail), tables = list())
       else list(text = if (S7::S7_inherits(x, CoxFit)) .render_cox(x, level, detail)
                        else if (S7::S7_inherits(x, FineGrayFit)) .render_fg(x, level, detail)
                        else .render_csc(x, level, detail), tables = list())
  r$text <- c(r$text, "", .method_lines(.method_keys(x)))
  r
}

.print_fit <- function(x, level = 0.95, ...) {
  cat(.render(x, level, FALSE)$text, sep = "\n")
  invisible(x)
}

.summary_fit <- function(object, level = 0.95, ...) {
  r <- .render(object, level, TRUE)
  FitSummary(model_call = object@model_call, coefficients = estimates(object, level = level), tables = r$tables, text = r$text)
}

S7::method(print, SurvXFit) <- .print_fit
S7::method(print, CSCFit) <- .print_fit
S7::method(summary, SurvXFit) <- .summary_fit
S7::method(summary, CSCFit) <- .summary_fit

S7::method(print, FitSummary) <- function(x, ...) {
  cat(x@text, sep = "\n")
  invisible(x)
}

S7::method(as.data.frame, FitSummary) <- function(x, ...) x@coefficients

S7::method(coef, FitSummary) <- function(object, ...) object@coefficients

# ----------------------------------------------------------------------------------------------
# Kaplan-Meier and Aalen-Johansen
# ----------------------------------------------------------------------------------------------
.default_times <- function(tmax) {
  tt <- pretty(c(0, tmax), n = 5L)
  tt[tt > 0 & tt <= tmax]
}

.at_risk <- function(x, lev, times) vapply(times, function(t) sum(x@time[x@group == lev] >= t), numeric(1))

S7::method(print, KMFit) <- function(x, ...) {
  cn <- x@counts
  df <- data.frame(group = cn$group, n = cn$n, events = cn$events, median = ifelse(is.na(cn$median), "not reached", format(cn$median, trim = TRUE)),
                   ci = paste0(ifelse(is.na(cn$median.lower), "-", format(cn$median.lower, trim = TRUE)), " to ",
                               ifelse(is.na(cn$median.upper), "not reached", format(cn$median.upper, trim = TRUE))),
                   stringsAsFactors = FALSE)
  names(df) <- c("group", "n", "events", "median survival", paste(.pct(x@level), "CI"))
  if (nrow(df) == 1L) df <- df[-1L]
  out <- c("Kaplan-Meier estimate of the survival function", .table_lines(df))
  if (!is.null(x@test))
    out <- c(out, sprintf("Log-rank test: chi-square = %.2f on %d df, p = %s", x@test$statistic, x@test$df, .pval(x@test$p.value)))
  cat(out, "", .method_lines(.method_keys(x)), sep = "\n")
  invisible(x)
}

.km_tables <- function(object, times) {
  if (is.null(times)) times <- .default_times(max(object@time))
  out <- lapply(names(object@curves), function(lev) {
    cv <- object@curves[[lev]]
    data.frame(time = times, n.risk = .at_risk(object, lev, times), surv = .step_at(cv$time, cv$surv, times, 1),
               se = .step_at(cv$time, cv$se, times, 0), lower = .step_at(cv$time, cv$lower, times, 1), upper = .step_at(cv$time, cv$upper, times, 1))
  })
  stats::setNames(out, names(object@curves))
}

S7::method(summary, KMFit) <- function(object, times = NULL, ...) {
  tabs <- .km_tables(object, times)
  text <- "Kaplan-Meier estimate of the survival function"
  for (lev in names(tabs)) {
    tb <- tabs[[lev]]
    df <- data.frame(time = format(tb$time, trim = TRUE), at.risk = tb$n.risk, survival = .num(tb$surv, 3L), ci = .ci_text(tb$lower, tb$upper, 3L),
                     stringsAsFactors = FALSE)
    names(df) <- c("time", "at risk", "survival", paste(.pct(object@level), "CI"))
    text <- c(text, if (length(tabs) > 1L) c("", paste0("Group: ", lev)), .table_lines(df))
  }
  if (!is.null(object@test))
    text <- c(text, "", sprintf("Log-rank test: chi-square = %.2f on %d df, p = %s", object@test$statistic, object@test$df, .pval(object@test$p.value)))
  FitSummary(model_call = object@model_call, coefficients = object@counts, tables = tabs, text = text)
}

S7::method(as.data.frame, KMFit) <- function(x, ...) {
  do.call(rbind, lapply(names(x@curves), function(lev) cbind(group = lev, x@curves[[lev]], stringsAsFactors = FALSE)))
}

S7::method(print, AJFit) <- function(x, ...) {
  cn <- x@counts
  if (nrow(cn) == 1L) cn <- cn[-1L]
  names(cn) <- sub("^cause\\.", "cause ", names(cn))
  out <- c("Aalen-Johansen estimate of the cumulative incidence functions", .table_lines(cn))
  if (!is.null(x@tests)) {
    tt <- data.frame(cause = x@tests$cause, chi = .num(x@tests$statistic, 2L), df = x@tests$df, p = .pval(x@tests$p.value), stringsAsFactors = FALSE)
    names(tt) <- c("cause", "chi-square", "df", "p")
    out <- c(out, "Gray's test for equal cumulative incidence across groups", .table_lines(tt))
  }
  out <- c(out, "Use summary(x, times = ...) for the cumulative incidence at chosen times.", "", .method_lines(.method_keys(x)))
  cat(out, sep = "\n")
  invisible(x)
}

.aj_tables <- function(object, times) {
  if (is.null(times)) times <- .default_times(max(object@time))
  out <- lapply(names(object@curves), function(lev) {
    cv <- object@curves[[lev]]
    tb <- data.frame(time = times, n.risk = .at_risk(object, lev, times))
    for (j in seq_len(object@K)) for (nm in c("cif", "se", "lower", "upper"))
      tb[[paste0(nm, ".", j)]] <- .step_at(cv$time, cv[[paste0(nm, ".", j)]], times, if (nm == "se") 0 else if (nm == "cif") 0 else NA_real_)
    tb
  })
  stats::setNames(out, names(object@curves))
}

S7::method(summary, AJFit) <- function(object, times = NULL, ...) {
  tabs <- .aj_tables(object, times)
  text <- "Aalen-Johansen estimate of the cumulative incidence functions"
  for (lev in names(tabs)) {
    tb <- tabs[[lev]]
    df <- data.frame(time = format(tb$time, trim = TRUE), at.risk = tb$n.risk, stringsAsFactors = FALSE)
    names(df) <- c("time", "at risk")
    for (j in seq_len(object@K)) {
      df[[paste("cause", j)]] <- .num(tb[[paste0("cif.", j)]], 3L)
      df[[paste0(.pct(object@level), " CI (", j, ")")]] <- .ci_text(tb[[paste0("lower.", j)]], tb[[paste0("upper.", j)]], 3L)
    }
    text <- c(text, if (length(tabs) > 1L) c("", paste0("Group: ", lev)), .table_lines(df))
  }
  if (!is.null(object@tests)) {
    tt <- data.frame(cause = object@tests$cause, chi = .num(object@tests$statistic, 2L), df = object@tests$df, p = .pval(object@tests$p.value), stringsAsFactors = FALSE)
    names(tt) <- c("cause", "chi-square", "df", "p")
    text <- c(text, "", "Gray's test for equal cumulative incidence across groups", .table_lines(tt))
  }
  FitSummary(model_call = object@model_call, coefficients = object@counts, tables = tabs, text = text)
}

S7::method(as.data.frame, AJFit) <- function(x, ...) {
  do.call(rbind, lapply(names(x@curves), function(lev) cbind(group = lev, x@curves[[lev]], stringsAsFactors = FALSE)))
}

# ----------------------------------------------------------------------------------------------
# model comparison
# ----------------------------------------------------------------------------------------------
S7::method(print, ModelComparison) <- function(x, ...) {
  cv <- x@folds >= 2L
  tab <- x@table
  out <- paste0("Comparison of prediction models for cause ", x@cause,
                if (cv) paste0(": ", x@folds, "-fold cross-validation") else ": apparent performance (no cross-validation)")
  col <- function(nm) if (cv) paste0(nm, "_cv") else nm
  for (t in unique(tab$time)) {
    s <- tab[tab$time == t, , drop = FALSE]
    s <- s[order(s[[col("brier")]], na.last = TRUE), , drop = FALSE]
    df <- data.frame(model = s$model, brier = .num(s[[col("brier")]], 4L), ipa = .num(100 * s[[col("ipa")]], 1L), auc = .num(s[[col("auc")]], 3L),
                     cal = .num(s[[col("calib_error")]], 3L), stringsAsFactors = FALSE)
    names(df) <- c("model (best first)", "Brier", "IPA %", "AUC", "calibration error")
    out <- c(out, "", paste0("Time horizon ", format(t, trim = TRUE)), .table_lines(df))
  }
  ft <- x@fit
  note <- rep("", nrow(ft))
  if (!is.null(ft$folds_not_converged)) note <- ifelse(ft$folds_not_converged > 0, paste0("not converged in ", ft$folds_not_converged, " fold(s)"), note)
  if (!is.null(ft$folds_failed)) note <- ifelse(ft$folds_failed > 0, paste0("failed in ", ft$folds_failed, " fold(s)"), note)
  if (!is.null(ft$converged)) note <- ifelse(!is.na(ft$converged) & !ft$converged, paste(note, "full-data fit not converged"), note)
  if (!is.null(ft$flags)) note <- ifelse(is.finite(ft$flags) & bitwAnd(as.integer(ifelse(is.finite(ft$flags), ft$flags, 0)), 15L) > 0L,
                                         trimws(paste(note, "full-data fit not converged")), note)
  df <- data.frame(model = ft$model, aic = .num(ft$aic, 1L), note = trimws(note), stringsAsFactors = FALSE)
  names(df) <- c("model", "AIC", "")
  out <- c(out, "", "Fit to the full data", .table_lines(df), "",
           "Brier: mean squared error of the predicted risk (smaller is better).",
           "IPA: gain over a prediction that ignores the covariates, in percent (larger is better; 0 = no gain).",
           "AUC: probability that a subject with the event by that time has the higher predicted risk.",
           "Calibration error: mean absolute difference between predicted and observed risk in risk groups.",
           "AIC is shown for models fitted by a full likelihood of the data: joint and constrained direct models and",
           "parametric cause-specific hazards; it can be compared among these models only.", "",
           .method_lines(.method_keys(x)))
  cat(out, sep = "\n")
  invisible(x)
}

S7::method(summary, ModelComparison) <- function(object, ...) object@table

S7::method(as.data.frame, ModelComparison) <- function(x, ...) x@table

S7::method(plot, ModelComparison) <- function(x, measure = c("ipa", "brier", "auc", "calib_error"), times = NULL,
                                               highlight = TRUE, reference = TRUE, col = "#0072B2", main = NULL,
                                               cex = 0.85, ...) {
  measure <- match.arg(measure)
  cv <- x@folds >= 2L
  nm <- if (cv) paste0(measure, "_cv") else measure
  tab <- x@table
  tt <- unique(tab$time); if (!is.null(times)) tt <- tt[tt %in% times]
  if (!length(tt)) stop("none of the requested times is in the comparison", call. = FALSE)
  smaller <- measure %in% c("brier", "calib_error")
  label <- c(ipa = "IPA (gain over the null model)", brier = "Brier score (smaller is better)",
             auc = "AUC (larger is better)", calib_error = "calibration error (smaller is better)")[[measure]]
  if (cv) label <- paste0(label, ", ", x@folds, "-fold cross-validated")
  # one common order of the models, best on top: mean rank over the time horizons
  mods <- unique(tab$model)
  rk <- sapply(tt, function(t) { s <- tab[tab$time == t, , drop = FALSE]; v <- s[[nm]][match(mods, s$model)]
                                 rank(if (smaller) v else -v, na.last = "keep") })
  ord <- order(rowMeans(matrix(rk, nrow = length(mods)), na.rm = TRUE), decreasing = TRUE)
  mods <- mods[ord]; ny <- length(mods); nt <- length(tt)
  # the labels are written once, at the left; their width is measured with the font metrics of the
  # device, so neither the names nor the axis title are cut
  din <- grDevices::dev.size("in")
  lab_in <- max(.text_inches(mods, 0.9 * cex)) + 0.15
  panel_in <- max(0.9, (din[1L] - lab_in - 0.2) / nt)
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::layout(matrix(seq_len(nt), 1L), widths = c(lab_in + panel_in, rep(panel_in, nt - 1L)))
  lab_l <- .wrap_inches(label, din[1L] - 0.4, cex)
  main_l <- if (is.null(main)) character(0) else .wrap_inches(main, din[1L] - 0.4, 1.1 * cex, font = 2L)
  graphics::par(cex = cex, oma = c(0.6 + 1.1 * length(lab_l), 0, if (length(main_l)) 0.6 + 1.3 * length(main_l) else 0.3, 0),
                mgp = c(2, 0.5, 0), tcl = -0.25)
  for (j in seq_len(nt)) {
    s <- tab[tab$time == tt[j], , drop = FALSE]; v <- s[[nm]][match(mods, s$model)]
    graphics::par(mai = c(0.35, if (j == 1L) lab_in else 0.1, 0.3, 0.1) * c(cex, 1, cex, 1))
    xr <- range(c(v, if (reference && measure == "ipa") 0), finite = TRUE); xr <- xr + c(-1, 1) * 0.08 * max(diff(xr), 1e-8)
    graphics::plot(NA, xlim = xr, ylim = c(0.5, ny + 0.5), yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "l", ...)
    at <- pretty(xr, n = 3L); at <- at[at >= xr[1L] & at <= xr[2L]]
    graphics::axis(1L, at = at, labels = format(at, trim = TRUE), cex.axis = 0.85)
    graphics::abline(h = seq_len(ny), col = "grey92", lwd = 0.8)
    if (reference && measure == "ipa") graphics::abline(v = 0, lty = 2, col = "grey55")
    best <- if (all(is.na(v))) integer(0) else if (smaller) which.min(v) else which.max(v)
    pc <- rep(col, ny); if (highlight && length(best)) pc[best] <- "#D55E00"
    graphics::points(v, seq_len(ny), pch = 19, col = pc, cex = 1.1)
    graphics::title(main = paste("time", format(tt[j], trim = TRUE)), font.main = 1, cex.main = 1, line = 0.6)
    if (j == 1L) graphics::axis(2L, at = seq_len(ny), labels = mods, las = 1L, tick = FALSE, cex.axis = 0.9, line = -0.5)
  }
  for (i in seq_along(lab_l)) graphics::mtext(lab_l[i], side = 1, outer = TRUE, line = 1.1 * (i - 1L) + 0.2, cex = cex)
  ln <- 0.2
  for (i in rev(seq_along(main_l))) { graphics::mtext(main_l[i], side = 3, outer = TRUE, line = ln, font = 2, cex = 1.1 * cex); ln <- ln + 1.3 }
  invisible(x)
}

# ----------------------------------------------------------------------------------------------
# plots
# ----------------------------------------------------------------------------------------------
.step_xy <- function(time, y, y0) list(x = c(0, time), y = c(y0, y))

S7::method(plot, KMFit) <- function(x, conf.int = length(x@curves) == 1L, risk.table = TRUE, legend = TRUE, mark.times = NULL,
                                    percent = FALSE, grid = FALSE, median = FALSE, labels = NULL, col = NULL, lty = 1, lwd = 2,
                                    xlab = "Time", ylab = "Survival probability", main = NULL, xlim = NULL, ylim = c(0, 1),
                                    cex = 0.85, ...) {
  cv <- x@curves; ng <- length(cv); keys <- if (is.null(labels)) names(cv) else rep_len(labels, ng)
  if (is.null(col)) col <- rep_len(.okabe_ito, ng)
  col <- rep_len(col, ng); lty <- rep_len(lty, ng)
  nrow_rt <- if (isTRUE(risk.table)) ng else 0L
  main_l <- if (is.null(main)) character(0) else .wrap_inches(main, grDevices::dev.size("in")[1L] - 0.6, 1.2 * cex, font = 2L)
  op <- graphics::par(mar = c(4.2 + if (nrow_rt) 1.6 + 1.1 * nrow_rt else 0, 4.5 + if (nrow_rt && ng > 1L) 2 else 0,
                              1.2 + 1.3 * max(1L, length(main_l)), 1), mgp = c(2.6, 0.7, 0), cex = cex)
  on.exit(graphics::par(op))
  xr <- if (is.null(xlim)) c(0, max(x@time)) else xlim
  graphics::plot(NA, xlim = xr, ylim = ylim, xlab = xlab, ylab = ylab, bty = "l", yaxt = "n", main = "", ...)
  for (i in rev(seq_along(main_l))) graphics::mtext(main_l[i], side = 3, line = 0.8 + 1.3 * (length(main_l) - i), font = 2, cex = 1.2 * cex)
  .prob_axis(percent)
  if (isTRUE(grid)) graphics::abline(h = graphics::axTicks(2L), col = "grey92", lwd = 0.8)
  for (g in seq_len(ng)) {
    s <- cv[[g]]; tt <- c(0, s$time)
    if (isTRUE(conf.int)) {
      pg <- .step_poly(tt, c(1, s$lower), c(1, s$upper))
      if (!is.null(pg)) graphics::polygon(pg$x, pg$y, col = grDevices::adjustcolor(col[g], alpha.f = 0.15), border = NA)
    }
    graphics::lines(tt, c(1, s$surv), type = "s", col = col[g], lty = lty[g], lwd = lwd)
    cens <- s$time[s$n.event == 0]
    if (length(cens)) graphics::points(cens, .step_at(s$time, s$surv, cens, 1), pch = 3, cex = 0.5, col = col[g])
    if (length(mark.times)) {
      v <- .step_at(s$time, s$surv, mark.times, 1)
      graphics::points(mark.times, v, pch = 19, cex = 0.6, col = col[g])
      graphics::text(mark.times, v, labels = if (percent) sprintf("%.1f%%", 100 * v) else sprintf("%.3f", v), pos = 3, cex = 0.7, col = col[g], offset = 0.35)
    }
    if (isTRUE(median) && is.finite(x@counts$median[g]))
      graphics::segments(c(xr[1L], x@counts$median[g]), c(0.5, 0.5), c(x@counts$median[g], x@counts$median[g]), c(0.5, 0), lty = 3, col = col[g])
  }
  if (length(mark.times)) graphics::abline(v = mark.times, lty = 3, col = "grey65")
  if (ng > 1L && isTRUE(legend)) graphics::legend("topright", legend = keys, col = col, lty = lty, lwd = lwd, bty = "n", cex = 0.85)
  if (nrow_rt) .risk_table(x@time, if (ng > 1L) factor(x@group, levels = names(cv), labels = keys) else NULL, xr)
  invisible(x)
}

S7::method(plot, CoxFit) <- function(x, newdata = NULL, times = NULL, xlab = "Time", ylab = "Survival probability",
                                     col = NULL, lty = 1, lwd = 2, legend = TRUE, labels = NULL, ...) {
  .no_tt(x, "plot()")
  if (is.null(times)) times <- sort(unique(c(0, x@time[x@status == 1L], max(x@time))))
  P <- stats::predict(x, newdata = newdata, times = times, type = "survival")
  if (is.null(newdata)) P <- matrix(colMeans(P), 1L)
  n <- nrow(P)
  if (is.null(col)) col <- seq_len(n)
  col <- rep_len(col, n); lty <- rep_len(lty, n)
  graphics::plot(NA, xlim = c(0, max(times)), ylim = c(0, 1), xlab = xlab, ylab = ylab, las = 1, bty = "l", ...)
  for (i in seq_len(n)) graphics::lines(times, P[i, ], type = "s", col = col[i], lty = lty[i], lwd = lwd)
  if (n > 1L && isTRUE(legend)) {
    if (is.null(labels)) labels <- if (!is.null(rownames(newdata))) rownames(newdata) else seq_len(n)
    graphics::legend("topright", legend = labels, col = col, lty = lty, lwd = lwd, bty = "n")
  }
  invisible(x)
}

# ----------------------------------------------------------------------------------------------
# further S3 generics, so that update(), model.matrix() and friends work as for lm() or coxph()
# ----------------------------------------------------------------------------------------------
S7::method(getCall, SurvXFit) <- function(x, ...) x@model_call
S7::method(getCall, CSCFit) <- function(x, ...) x@model_call
S7::method(terms, SurvXFit) <- function(x, ...) x@terms
S7::method(terms, CSCFit) <- function(x, ...) x@terms
S7::method(model.matrix, SurvXFit) <- function(object, ...) object@X
S7::method(model.matrix, CSCFit) <- function(object, ...) object@X

# the partial likelihood of cause-specific Cox models factorises over causes (Prentice et al. 1978),
# so the log partial likelihood of the whole model is the sum over causes
S7::method(logLik, CSCFit) <- function(object, ...) {
  structure(sum(vapply(object@models, function(m) m$loglik, numeric(1))),
            df = sum(vapply(object@models, function(m) length(m$coefficients), integer(1))), nobs = object@n, class = "logLik")
}

#' Methods for fitted models
#'
#' Every model of SurvX is an S7 object: its parts are read with \code{@} (\code{fit@coefficients},
#' \code{fit@vcov}, \code{fit@converged}; \code{S7::props(fit)} lists them all). The usual S3 generics
#' of R work as for \code{lm()} or \code{survival::coxph()}.
#'
#' \describe{
#'   \item{\code{print(fit)}}{covariate effects as hazard ratios, subdistribution hazard ratios or odds
#'     ratios with confidence intervals and p-values; for direct models also the long-run risk
#'     (plateau) at typical covariate values and a warning when the fit did not converge.}
#'   \item{\code{summary(fit, level = 0.95)}}{everything in \code{print} and in addition the predicted
#'     cumulative incidence at the quartiles of the event times, the baseline parameters on their
#'     natural scale and a note on how to read the effects. The result is an object whose
#'     \code{@coefficients} is the table of \code{\link{estimates}} and \code{@tables} the other tables;
#'     \code{as.data.frame()} returns the coefficient table.}
#'   \item{\code{coef}, \code{vcov}, \code{confint}, \code{logLik}, \code{AIC}, \code{BIC}, \code{nobs}}{as
#'     usual; coefficients are on the scale used for estimation (log ratios, log or logit baseline
#'     parameters). For \code{csc()} fits \code{coef}, \code{vcov} and \code{confint} return one element per cause.}
#'   \item{\code{estimates(fit)}}{all of this in one data frame, including the exponentiated effects.}
#'   \item{\code{predict(fit, newdata, times, cause)}}{absolute risks; see \code{\link{predict.SurvX}}.}
#'   \item{\code{plot(fit, ...)}}{for \code{directcif()}, \code{fg()}, \code{csc()} and \code{aj()} the same as
#'     \code{\link{cif_plot}(fit, ...)}: cumulative incidence curves of every cause with confidence bands and
#'     the number at risk; every element can be switched off. For \code{cox()} predicted survival curves.}
#'   \item{\code{formula}, \code{terms}, \code{model.matrix}, \code{getCall}, \code{update}}{as for \code{lm()};
#'     \code{update(fit, link = "PO")} refits with a changed argument.}
#'   \item{\code{print}, \code{summary(x, times)}, \code{plot}, \code{as.data.frame} for \code{km()} and \code{aj()}}{
#'     median survival with its confidence limits, survival or cumulative incidence with confidence
#'     intervals and numbers at risk at chosen times, tests between groups, curves with confidence bands.}
#'   \item{\code{print}, \code{summary}, \code{plot}, \code{as.data.frame} for \code{compare_models()} and \code{cv_compare()}}{
#'     models ranked by the Brier score at every time horizon, with a short explanation of each measure.}
#'   \item{\code{plot(x, measure = "ipa", times = NULL, highlight = TRUE, reference = TRUE, col, main, cex = 0.85)} for a comparison}{
#'     one panel per time horizon with the models in the same order in every panel (best on top, by the
#'     mean rank over the horizons); the best model at each horizon is drawn in orange
#'     (\code{highlight = FALSE} switches this off) and, for the IPA, a dashed line marks no gain over the
#'     null model (\code{reference = FALSE} switches it off). \code{measure} is one of \code{"ipa"},
#'     \code{"brier"}, \code{"auc"}, \code{"calib_error"}.}
#'   \item{\code{cex}}{every plot (\code{plot()} of \code{km()}, of a comparison, and \code{cif_plot()}) has a
#'     \code{cex} argument for the size of all text, default 0.85; titles, captions and axis titles are
#'     wrapped to the width of the device so that they never run off the page.}
#' }
#'
#' @name SurvX-methods
#' @examples
#' set.seed(5)
#' n <- 300
#' x <- rbinom(n, 1, 0.5)
#' t1 <- rexp(n, 0.4 * exp(0.6 * x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
#' d <- data.frame(time = pmin(t1, t2, cens), x = x)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#' fit <- directcif(Event(time, status) ~ x, d)
#' fit
#' summary(fit)
#' confint(fit)
#' estimates(fit)
#' fit_po <- update(fit, link = "PO")
#' AIC(fit, fit_po)
#' plot(fit, newdata = data.frame(x = 0:1), labels = c("x = 0", "x = 1"))
NULL
