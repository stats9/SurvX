# Parametric cure models. The likelihood, its derivatives, the fit and the predictions are in the
# C++ core (src/survcore/cure.hpp); this file prepares data and output.

.cure_types <- c(survival = 0L, hazard = 1L, uncured = 2L, cured = 3L, cumhaz = 4L)

#' Parametric cure models
#'
#' Mixture and non-mixture (promotion time) cure models: a fraction of the subjects never
#' experiences the event, and the cured fraction depends on covariates through a logistic model.
#' The survival of those who are not cured (the latency) is any distribution of
#' \code{\link{parsurv}}, as a PH or AFT model, with its own covariates.
#'
#' @param formula \code{Event(time, status) ~ latency covariates}; status 0/1;
#'   \code{Event(time, status, entry = )} for delayed entry.
#' @param data a data frame.
#' @param cure one-sided formula of the covariates of the cured fraction, e.g. \code{~ age + stage};
#'   \code{NULL} (default): the covariates of \code{formula}; \code{~ 1}: a common cured fraction.
#' @param baseline distribution of the latency, as in \code{\link{parsurv}}.
#' @param scale \code{"PH"} or \code{"AFT"}, as in \code{\link{parsurv}}.
#' @param type \code{"mixture"} (default) or \code{"nonmixture"}.
#' @param df,knots,cuts spline and piecewise exponential latencies, as in \code{\link{parsurv}}.
#' @param init optional starting values.
#' @param max_iter,tol iteration limit and convergence tolerance.
#' @export
cure <- function(formula, data, cure = NULL, baseline = "weibull", scale = NULL, type = c("mixture", "nonmixture"),
                 df = 3L, knots = NULL, cuts = NULL, init = NULL, max_iter = 200L, tol = 1e-10) {
  cl <- match.call(); type <- match.arg(type)
  dist <- .par_dist(baseline); scale <- .par_scale(dist, scale)
  f <- .sx_frame(formula, data, entry_ok = TRUE, what = "cure()")
  if (any(f$status > 1L)) stop("cure() is for one event (status 0/1)", call. = FALSE)
  if (any(f$time <= 0)) stop("times must be positive for parametric models", call. = FALSE)
  # covariates of the cured fraction, on the rows kept by the model frame of `formula`
  cf <- if (is.null(cure)) stats::delete.response(f$terms) else stats::terms(cure)
  if (length(attr(cf, "response")) && attr(cf, "response") != 0L) stop("cure must be a one-sided formula, e.g. ~ age", call. = FALSE)
  if (missing(data) || is.null(data)) data <- environment(formula)
  mfc <- stats::model.frame(cf, data = data, na.action = stats::na.pass)
  omit <- attr(f$mf, "na.action")
  if (!is.null(omit)) mfc <- mfc[-omit, , drop = FALSE]
  if (nrow(mfc) != length(f$time)) stop("the covariates of the cured fraction must come from the same data", call. = FALSE)
  Xc <- stats::model.matrix(cf, mfc)
  if (anyNA(Xc)) stop("missing values in the covariates of the cured fraction", call. = FALSE)
  cc <- attr(Xc, "contrasts")
  Xc <- Xc[, colnames(Xc) != "(Intercept)", drop = FALSE]
  kn <- .par_knots(dist, f$time, f$status, df, knots, cuts)
  r <- .cure_fit(f$time, f$status, f$entry, f$X, Xc, dist, scale, kn, type == "mixture", init, max_iter, tol)
  CureFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
          coefficients = r$theta, vcov = r$var, loglik = r$loglik, n = length(f$time),
          nevent = as.integer(r$nevent), iterations = as.integer(r$iterations), converged = isTRUE(r$converged),
          flags = as.integer(r$flags), time = f$time, status = f$status, X = f$X,
          dist = dist, scale = scale, dist_code = .par_dists[[dist]], aft = as.integer(scale == "AFT"),
          entry = f$entry, knots = kn, mixture = type == "mixture", Xc = Xc,
          cure_terms = stats::delete.response(cf), cure_xlevels = stats::.getXlevels(cf, mfc), cure_contrasts = cc,
          n_latency = as.integer(1L + ncol(f$X) + .nshape(dist, kn)))
}

.cure_fit <- function(time, status, entry, Xl, Xc, dist, scale, kn, mixture, init = NULL, max_iter = 200L, tol = 1e-10) {
  n <- length(time); pl <- ncol(Xl); pc <- ncol(Xc)
  P <- 1L + pl + .nshape(dist, kn) + 1L + pc
  r <- .C("SX_cure_fit", as.integer(n), as.integer(pl), as.integer(pc), .par_dists[[dist]], as.integer(scale == "AFT"),
          as.integer(length(kn)), if (length(kn)) as.double(kn) else double(1), as.integer(mixture),
          as.integer(!is.null(entry)), if (is.null(entry)) double(1) else as.double(entry),
          as.double(time), as.integer(status), .rm(Xl), .rm(Xc), as.integer(!is.null(init)),
          if (is.null(init)) double(P) else as.double(init), as.integer(max_iter), as.double(tol),
          theta = double(P), var = double(P * P), info_d = double(3), info_i = integer(4), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "cure model fit")
  nm <- c("(Intercept)", colnames(Xl), .shape_name(dist, scale, kn), "cure: (Intercept)",
          if (pc) paste0("cure: ", colnames(Xc)) else character(0))
  V <- matrix(r$var, P, P, byrow = TRUE, dimnames = list(nm, nm))
  ok <- all(is.finite(r$theta)) && all(is.finite(diag(V))) && all(diag(V) >= 0) && is.finite(r$info_d[2L])
  list(theta = stats::setNames(r$theta, nm), var = V,
       loglik = r$info_d[2L], iterations = r$info_i[1L], converged = as.logical(r$info_i[2L]) && ok,
       flags = r$info_i[3L], nevent = r$info_i[4L])
}

# covariate rows of new data for the latency and the cured fraction
.cure_newX <- function(object, newdata) {
  if (is.null(newdata)) return(list(Xl = object@X, Xc = object@Xc))
  Xl <- .sx_newX(object, newdata)
  mf <- stats::model.frame(object@cure_terms, newdata, na.action = stats::na.pass, xlev = object@cure_xlevels)
  Xc <- stats::model.matrix(object@cure_terms, mf, contrasts.arg = object@cure_contrasts)
  list(Xl = Xl, Xc = Xc[, colnames(Xc) != "(Intercept)", drop = FALSE])
}

.cure_predict_X <- function(object, Xl, Xc, at, type, se, level, average = FALSE) {
  nn <- nrow(Xl); nt <- length(at)
  usevar <- isTRUE(se) && all(is.finite(object@vcov))
  cells <- if (average) nt else nn * nt
  kn <- object@knots
  r <- .C("SX_cure_predict", object@dist_code, object@aft, as.integer(length(kn)), if (length(kn)) as.double(kn) else double(1),
          as.integer(ncol(Xl)), as.integer(ncol(Xc)), as.integer(object@mixture), as.double(object@coefficients),
          as.integer(usevar), if (usevar) as.double(object@vcov) else double(1), as.integer(nn), .rm(Xl), .rm(Xc),
          as.integer(nt), as.double(at), .cure_types[[type]], as.integer(average),
          out = double(cells), se = double(cells), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "predict")
  cn <- if (type == "cured") "cured" else paste0("t=", at)
  nr <- if (average) 1L else nn
  est <- .matrix_rm(r$out, nr, nt, list(NULL, if (type == "cured") rep("cured", nt) else cn))
  if (type == "cured") est <- est[, 1L, drop = FALSE]
  if (!isTRUE(se)) return(est)
  s <- if (usevar) .matrix_rm(r$se, nr, nt) else est * NA_real_
  if (type == "cured") s <- s[, 1L, drop = FALSE]
  lim <- .pred_limits(est, s, if (type %in% c("survival", "uncured", "cured")) "survival" else "hazard", level)
  list(est = est, se = s, lower = lim$lower, upper = lim$upper)
}

S7::method(predict, CureFit) <- function(object, newdata = NULL, times = NULL,
                                         type = c("survival", "hazard", "cumhaz", "uncured", "cured"),
                                         se = FALSE, level = 0.95, average = FALSE, ...) {
  type <- match.arg(type)
  if (type != "cured" && (is.null(times) || !length(times))) stop("times must be given", call. = FALSE)
  X <- .cure_newX(object, newdata)
  .cure_predict_X(object, X$Xl, X$Xc, if (type == "cured") 1 else times, type, se, level, average)
}

S7::method(residuals, CureFit) <- function(object, type = c("martingale", "deviance", "coxsnell"), ...) {
  type <- match.arg(type)
  n <- object@n
  H <- vapply(seq_len(n), function(i)
    .cure_predict_X(object, object@X[i, , drop = FALSE], object@Xc[i, , drop = FALSE], object@time[i], "cumhaz", FALSE, 0.95)[1L, 1L], numeric(1))
  if (!is.null(object@entry)) H <- H - vapply(seq_len(n), function(i)
    if (object@entry[i] > 0) .cure_predict_X(object, object@X[i, , drop = FALSE], object@Xc[i, , drop = FALSE], object@entry[i], "cumhaz", FALSE, 0.95)[1L, 1L] else 0,
    numeric(1))
  m <- object@status - H
  switch(type, coxsnell = H, martingale = m,
         deviance = sign(m) * sqrt(-2 * (m + ifelse(object@status == 1L, log(pmax(object@status - m, 1e-300)), 0))))
}


# the latency part as a parametric survival fit (for the baseline table)
.cure_latency <- function(x) {
  L <- x@n_latency; i <- seq_len(L)
  ParSurvFit(call = x@call, terms = x@terms, xlevels = x@xlevels, contrasts = x@contrasts,
             coefficients = x@coefficients[i], vcov = x@vcov[i, i, drop = FALSE], loglik = x@loglik, n = x@n,
             nevent = x@nevent, iterations = x@iterations, converged = x@converged, flags = x@flags,
             time = x@time, status = x@status, X = x@X, dist = x@dist, scale = x@scale, dist_code = x@dist_code,
             aft = x@aft, entry = x@entry, loglik_null = NA_real_, event_label = "the event", knots = x@knots)
}

.reading_cure <- function(e, xlevels, level) {
  e <- e[is.finite(e$p), , drop = FALSE]
  if (!nrow(e)) return(character(0))
  sig <- e[e$p < 0.05, , drop = FALSE]
  if (!nrow(sig)) return(c("Reading:", "  No covariate effect on the cured fraction reached p < 0.05."))
  lab <- .pretty_terms(sub("^cure: ", "", sig$term), xlevels)
  is_level <- grepl(": ", lab, fixed = TRUE) & grepl(" vs ", lab, fixed = TRUE)
  who <- paste0(ifelse(is_level, sub(": ", " (", lab, fixed = TRUE), paste0(lab, " (per one-unit increase")), ")")
  out <- "Reading (effects with p < 0.05):"
  for (i in seq_len(nrow(sig))) {
    r <- sig$exp_estimate[i]
    rt <- .ratio_text(r, sig$exp_lower[i], sig$exp_upper[i])
    out <- c(out, paste0("  - ", who[i], ": ", if (r > 1) "higher" else "lower", " odds of being cured (OR ", rt$ratio,
                         ", ", .pct(level), " CI ", rt$ci, ")."))
  }
  out
}

.render_cure <- function(x, level, detail) {
  e <- estimates(x, level = level)
  lat <- .cure_latency(x)
  out <- c(paste0(if (x@mixture) "Mixture cure model" else "Non-mixture (promotion time) cure model",
                  ": latency ", .par_model_label(x@dist, x@scale)),
           paste0("  Subjects: ", x@n, "   events: ", x@nevent, "   censored: ", x@n - x@nevent,
                  if (!is.null(x@entry)) "   delayed entry: yes" else ""),
           paste0("  ", if (x@mixture) "S(t) = pi + (1 - pi) S_u(t)" else "S(t) = pi ^ F_u(t)",
                  ", pi = cured fraction (logistic model), S_u = survival of the uncured"))
  # cured fraction
  pr <- tryCatch(.cure_predict_X(x, matrix(colMeans(x@X), 1L), matrix(colMeans(x@Xc), 1L), 1, "cured", TRUE, level),
                 error = function(err) NULL)
  out <- c(out, "", "Cured fraction")
  if (!is.null(pr))
    out <- c(out, sprintf("  At the mean covariate values: %s (%s CI %s)", .num(pr$est[1L, 1L], 3L), .pct(level),
                          .ci_text(pr$lower[1L, 1L], pr$upper[1L, 1L], 3L)))
  ec <- e[e$role == "cure", , drop = FALSE]
  if (nrow(ec)) {
    ec$term <- sub("^cure: ", "", ec$term)
    out <- c(out, "  Effect of covariates: odds ratio of being cured (OR)", .effect_lines(ec, "OR", level, x@cure_xlevels),
             .reading_cure(e[e$role == "cure", , drop = FALSE], x@cure_xlevels, level))
  }
  # latency
  eb <- e[e$role == "beta", , drop = FALSE]
  if (nrow(eb)) {
    lab <- if (x@scale == "AFT") c("TR", "time ratio (TR)") else c("HR", "hazard ratio (HR)")
    out <- c(out, "", paste0("Latency (subjects who are not cured): ", lab[2L]), .effect_lines(eb, lab[1L], level, x@xlevels),
             .reading_lines(e, x@xlevels, x@scale, "the event among the uncured", level))
  }
  if (detail)
    out <- c(out, "", "Latency distribution (natural scale)",
             .table_lines(stats::setNames(.natural_baseline(lat, level), c("", "estimate", paste(.pct(level), "CI")))))
  P <- length(x@coefficients)
  out <- c(out, "", sprintf("Fit: log-likelihood %.2f | AIC %.2f | BIC %.2f | %d parameters", x@loglik,
                            -2 * x@loglik + 2 * P, -2 * x@loglik + log(x@n) * P, P))
  if (detail)
    out <- c(out, "",
             "OR > 1: a higher chance of being cured; the cured fraction is identified only if follow-up is long enough",
             "  for the survival curve to level off (check plot(fit)).",
             if (x@scale == "AFT") "TR > 1 (latency): among the uncured the event happens later."
             else "HR > 1 (latency): among the uncured the event occurs at a higher rate.",
             "Predictions: predict(fit, newdata, times, type = \"survival\" | \"uncured\" | \"hazard\"), and type = \"cured\".")
  c(out, .status_lines(x@converged, x@flags))
}

S7::method(print, CureFit) <- function(x, level = 0.95, ...) {
  cat(.render_cure(x, level, FALSE), "", .method_lines(.method_keys(x)), sep = "\n")
  invisible(x)
}

S7::method(summary, CureFit) <- function(object, level = 0.95, ...) {
  FitSummary(call = object@call, coefficients = estimates(object, level = level), tables = list(),
             text = c(.render_cure(object, level, TRUE), "", .method_lines(.method_keys(object))))
}

# fitted marginal survival (average over the subjects) with the Kaplan-Meier estimate, and the
# average cured fraction as a dashed line
S7::method(plot, CureFit) <- function(x, times = NULL, col = "#0072B2", xlab = "Time", ylab = "Survival probability",
                                      main = NULL, cex = 0.85, ...) {
  tg <- if (is.null(times)) seq(0, max(x@time), length.out = 201L)[-1L] else times
  S <- .cure_predict_X(x, x@X, x@Xc, tg, "survival", FALSE, 0.95, average = TRUE)[1L, ]
  pc <- mean(.cure_predict_X(x, x@X, x@Xc, 1, "cured", FALSE, 0.95))
  op <- graphics::par(cex = cex, mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op))
  inc <- .na_increments(x@time, x@status, 1L, x@entry)
  km_t <- c(0, inc$time); km_s <- c(1, cumprod(1 - inc$d / inc$Y))
  graphics::plot(km_t, km_s, type = "s", col = "grey40", ylim = c(0, 1), xlim = c(0, max(tg)), xlab = xlab, ylab = ylab,
                 main = if (is.null(main)) "Fitted survival (average over subjects) and Kaplan-Meier" else main, las = 1, bty = "l", ...)
  graphics::lines(c(0, tg), c(1, S), col = col, lwd = 2)
  graphics::abline(h = pc, lty = 2, col = col)
  graphics::legend("topright", c("Kaplan-Meier", "cure model", "cured fraction"), col = c("grey40", col, col),
                   lty = c(1, 1, 2), lwd = c(1, 2, 1), bty = "n")
  invisible(x)
}
