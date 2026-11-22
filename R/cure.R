# Parametric cure models. The likelihood, its derivatives, the fit and the predictions are in the
# C++ core (src/survcore/cure.hpp); this file prepares data and output.

.cure_types <- c(survival = 0L, hazard = 1L, uncured = 2L, cured = 3L, cumhaz = 4L)

#' Cure models: survival with a cured fraction
#'
#' @description
#' Fits a parametric cure model: a fraction of the population will never have the event (they are
#' "cured"), and the others have it after a time described by a parametric distribution (the
#' \emph{latency}). Covariates can act on the cured fraction and on the latency separately. Mixture
#' (Boag 1949; Berkson and Gage 1952; Farewell 1982) and non-mixture or promotion-time models
#' (Yakovlev and Tsodikov 1996; Chen, Ibrahim and Sinha 1999) are available, with delayed entry.
#'
#' \strong{Use it when} the Kaplan-Meier curve levels off at a clearly positive value with long
#' follow-up (many cancers after curative treatment, recurrence after surgery), and the questions are
#' "what proportion is cured?" and "when do the others fail?".
#' \strong{Use something else when} there is no plateau, or follow-up is too short to see one: the
#' cured fraction is then not identified and the model only extrapolates. Use \code{\link{parsurv}}
#' or \code{\link{cox}}.
#' @param formula \code{Event(time, status) ~ latency covariates}, status 0/1;
#'   \code{Event(time, status, entry = entry)} for delayed entry. The right-hand side gives the
#'   covariates of the latency (the timing of the event among those who will have it).
#' @param data a data frame. Rows with missing values are removed.
#' @param cure one-sided formula with the covariates of the cured fraction, for example
#'   \code{~ age + stage}. \code{NULL} (default) uses the covariates of \code{formula}; \code{~ 1}
#'   gives one cured fraction for everybody. The two sets may differ, overlap or coincide.
#' @param baseline distribution of the latency: \code{"weibull"} (default), \code{"exponential"},
#'   \code{"gompertz"}, \code{"lognormal"}, \code{"loglogistic"}, \code{"gengamma"}, \code{"rp"} or
#'   \code{"pwexp"}, as in \code{\link{parsurv}}.
#' @param scale \code{"PH"} or \code{"AFT"} for the latency covariates, as in \code{\link{parsurv}}.
#' @param type \code{"mixture"} (default): the population is a mixture of cured subjects and
#'   subjects who will fail. \code{"nonmixture"}: the promotion-time model, in which the cumulative
#'   hazard of the whole population is bounded. See Details.
#' @param df, knots, cuts spline and piecewise exponential latencies, as in \code{\link{parsurv}}.
#' @param init optional starting values, in the order of \code{coef(fit)}.
#' @param max_iter, tol iteration limit and convergence tolerance.
#' @details
#' \subsection{The models}{
#' With \eqn{\pi(z) = \mathrm{expit}(\gamma^\top z)} the cured fraction and \eqn{S_u(t \mid x)} the
#' survival of the uncured (the latency distribution),
#' \describe{
#'   \item{mixture}{\eqn{S(t \mid x, z) = \pi(z) + \{1 - \pi(z)\} S_u(t \mid x)}.}
#'   \item{non-mixture}{\eqn{S(t \mid x, z) = \pi(z)^{F_u(t \mid x)}} with \eqn{F_u = 1 - S_u}: the
#'     cumulative hazard \eqn{-\log\pi(z) F_u(t \mid x)} is bounded and the hazard is proportional to
#'     \eqn{-\log \pi(z)}, as in the promotion time model. In this parametrisation (Lambert et al.
#'     2007, and \pkg{flexsurvcure}) \eqn{\pi} is again the cured fraction.}
#' }
#' In both, \eqn{S(t) \to \pi(z)} as \eqn{t \to \infty}.
#' }
#' \subsection{Estimation}{
#' Maximum likelihood by Newton-Raphson with exact derivatives (those of \code{\link{parsurv}} for
#' the latency; numerical for the generalized gamma); standard errors from the inverse observed
#' information. With delayed entry every contribution is divided by \eqn{S(entry \mid x, z)}.
#' Starting values: the latency fitted to the event times alone, and the logit of the Kaplan-Meier
#' estimate at the last time for the cured fraction.
#' }
#' \subsection{Identifiability}{
#' A cured fraction can be estimated only when follow-up is long enough for the survival curve to
#' level off, and when the latency distribution is plausible beyond the last events: the estimate of
#' \eqn{\pi} rests on the plateau (Maller and Zhou 1996; Amico and Van Keilegom 2018).
#' \code{plot(fit)} compares the fitted population survival with the Kaplan-Meier estimate and marks
#' the average cured fraction.
#' }
#' @section Interpreting the coefficients:
#' \describe{
#'   \item{Cured fraction (\code{cure: ...})}{log odds ratios: \eqn{e^\gamma} is the \strong{odds
#'     ratio of being cured}. OR 2.0 for a treatment means twice the odds of cure; OR < 1 means less
#'     cure. \code{predict(fit, newdata, type = "cured")} gives the cured fraction itself.}
#'   \item{Latency (the covariates of \code{formula})}{hazard ratios (PH) or time ratios (AFT)
#'     \emph{among the uncured} (mixture model): they say whether those who will have the event have
#'     it earlier or later, not whether they have it.}
#'   \item{Non-mixture model}{the latency coefficients describe the timing distribution \eqn{F_u};
#'     the cure coefficients change both the cured fraction and, proportionally, the hazard of the
#'     whole population.}
#' }
#' A covariate can increase the chance of cure and at the same time shorten the time to the event
#' among the uncured; report both parts. A sentence for a report: \emph{"In a Weibull mixture cure
#' model, the estimated cured fraction was 45\% (95\% CI 38 to 52\%) for stage I and 21\% (15 to
#' 28\%) for stage II (odds ratio of cure 0.32); among patients who relapsed, relapse occurred later
#' with treatment (HR 0.70)."}
#' @section Assumptions and how to check them:
#' \describe{
#'   \item{A plateau exists}{\code{plot(km(...))}: the Kaplan-Meier curve levels off with many
#'     subjects still at risk; \code{plot(fit)} shows whether the fitted curve follows it.}
#'   \item{Latency distribution}{compare baselines with \code{AIC()}; the cured fraction can depend
#'     strongly on the distribution when the plateau is short, so report the sensitivity.}
#' }
#' @section Messages and common problems:
#' \describe{
#'   \item{\code{WARNING: the fit did NOT converge}}{often no plateau (the cured fraction drifts to 0)
#'     or too few events. Check the Kaplan-Meier curve; try a simpler latency or \code{cure = ~ 1}.}
#'   \item{A cured fraction near 0 or 1 with a huge standard error}{the cured fraction is not
#'     identified by these data.}
#' }
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}; 62 data sets): the maximised log-likelihood agrees with
#' \code{flexsurvcure::flexsurvcure} to a relative \eqn{3 \times 10^{-10}} (mixture) and \eqn{3 \times
#' 10^{-8}} (non-mixture), and predicted survival within the tolerance of its optimiser (\eqn{2 \times
#' 10^{-4}}); the log-likelihood also agrees with independent code to \eqn{10^{-12}}. In a Monte Carlo
#' study (1000 replications, \eqn{n} = 200 to 1000) the 95\% intervals of the cure and latency
#' parameters covered 93.1 to 96.0\%. On \eqn{n = 10000}, a Weibull mixture cure fit takes about
#' 0.03 s (\code{flexsurvcure}: 2.9 s).
#' @return
#' An S7 object of class \code{CureFit}. \code{coef(fit)}: the latency coefficients (as in
#' \code{\link{parsurv}}) followed by those of the cured fraction, named \code{"cure: ..."}.
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(fit)}, \code{summary(fit)} \tab latency (HR or TR) and cure (odds ratio) tables,
#'     average cured fraction, \code{Reading:} lines \cr
#'   \code{\link{estimates}(fit)}, \code{coef}, \code{vcov}, \code{confint}, \code{logLik},
#'     \code{AIC} \tab extractors \cr
#'   \code{predict(fit, newdata, times, type)} \tab \code{"survival"} (population),
#'     \code{"hazard"}, \code{"cumhaz"}, \code{"uncured"} (survival of the uncured) and
#'     \code{"cured"} (the cured fraction; no \code{times}); \code{se = TRUE} adds delta-method
#'     standard errors \cr
#'   \code{plot(fit)} \tab fitted population survival against Kaplan-Meier, with the cured fraction \cr
#'   \code{residuals(fit, type)} \tab Cox-Snell, martingale, deviance \cr
#' }
#' @references
#' \emph{Mixture cure model:}
#' Boag JW (1949). Maximum likelihood estimates of the proportion of patients cured by cancer
#'   therapy. Journal of the Royal Statistical Society, Series B 11:15-44.
#' Berkson J, Gage RP (1952). Survival curve for cancer patients following treatment. Journal of the
#'   American Statistical Association 47:501-515.
#' Farewell VT (1982). The use of mixture models for the analysis of survival data with long-term
#'   survivors. Biometrics 38:1041-1046.
#'
#' \emph{Non-mixture (promotion time) model:}
#' Yakovlev AY, Tsodikov AD (1996). Stochastic Models of Tumor Latency and Their Biostatistical
#'   Applications. World Scientific.
#' Chen MH, Ibrahim JG, Sinha D (1999). A new Bayesian model for survival data with a surviving
#'   fraction. Journal of the American Statistical Association 94:909-919.
#' Lambert PC, Thompson JR, Weston CL, Dickman PW (2007). Estimating and modeling the cure fraction
#'   in population-based cancer survival analysis. Biostatistics 8:576-594.
#'
#' \emph{Identifiability and review:}
#' Maller RA, Zhou X (1996). Survival Analysis with Long-Term Survivors. Wiley.
#' Amico M, Van Keilegom I (2018). Cure models in survival analysis. Annual Review of Statistics and
#'   Its Application 5:311-342.
#'
#' All references: \code{sx_references("cure")}.
#' @seealso
#' \code{\link{parsurv}} (the latency distributions), \code{\link{km}} (look for the plateau),
#' \code{\link{predict.SurvX}}.
#' @examples
#' set.seed(6)
#' n <- 400
#' d <- data.frame(age = rnorm(n), stage = rbinom(n, 1, 0.5))
#' cured <- runif(n) < plogis(0.5 - 1.2 * d$stage)                 # stage 1: less cure
#' t <- ifelse(cured, Inf, rweibull(n, 1.3, 2 * exp(-0.3 * d$age)))
#' cens <- runif(n, 0, 10)
#' d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
#'
#' ## 1. look for a plateau
#' plot(km(Event(time, status) ~ stage, d))
#'
#' ## 2. mixture cure model: latency depends on age, cure on stage
#' fit <- cure(Event(time, status) ~ age, d, cure = ~ stage)
#' fit
#'
#' ## 3. cured fraction by stage, and survival of the uncured
#' nd <- data.frame(age = 0, stage = 0:1)
#' predict(fit, nd, type = "cured", se = TRUE)
#' predict(fit, nd, times = c(1, 3), type = "uncured")
#' plot(fit)
#'
#' ## 4. non-mixture model for comparison
#' AIC(fit, cure(Event(time, status) ~ age, d, cure = ~ stage, type = "nonmixture"))
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
  .sx_check_collinear(Xc, "covariates of the cured fraction")
  kn <- .par_knots(dist, f$time, f$status, df, knots, cuts)
  r <- .cure_fit(f$time, f$status, f$entry, f$X, Xc, dist, scale, kn, type == "mixture", init, max_iter, tol)
  CureFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
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
  kn <- object@knots; he <- !is.null(object@entry); n <- object@n
  r <- .C("SX_cure_residuals", object@dist_code, object@aft, as.integer(length(kn)), if (length(kn)) as.double(kn) else double(1),
          as.integer(ncol(object@X)), as.integer(ncol(object@Xc)), as.integer(object@mixture), as.double(object@coefficients),
          as.integer(n), .rm(object@X), .rm(object@Xc), as.double(object@time), as.integer(he), if (he) as.double(object@entry) else double(1),
          as.integer(object@status), match(type, c("coxsnell", "martingale", "deviance")) - 1L, out = double(n), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "residuals")
  r$out
}


# the latency part as a parametric survival fit (for the baseline table)
.cure_latency <- function(x) {
  L <- x@n_latency; i <- seq_len(L)
  ParSurvFit(model_call = x@model_call, terms = x@terms, xlevels = x@xlevels, contrasts = x@contrasts,
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
  pr <- tryCatch(.cure_predict_X(x, matrix(.sx_colmeans(x@X), 1L), matrix(.sx_colmeans(x@Xc), 1L), 1, "cured", TRUE, level),
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
  ic <- .sx_ic(x@loglik, P, x@n)
  out <- c(out, "", sprintf("Fit: log-likelihood %.2f | AIC %.2f | BIC %.2f | %d parameters", x@loglik, ic$aic, ic$bic, P))
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
  FitSummary(model_call = object@model_call, coefficients = estimates(object, level = level), tables = list(),
             text = c(.render_cure(object, level, TRUE), "", .method_lines(.method_keys(object))))
}

# fitted marginal survival (average over the subjects) with the Kaplan-Meier estimate, and the
# average cured fraction as a dashed line
S7::method(plot, CureFit) <- function(x, times = NULL, col = "#0072B2", xlab = "Time", ylab = "Survival probability",
                                      main = NULL, cex = 0.85, ...) {
  tg <- if (is.null(times)) .sx_seq(0, max(x@time), 201L)[-1L] else times
  S <- .cure_predict_X(x, x@X, x@Xc, tg, "survival", FALSE, 0.95, average = TRUE)[1L, ]
  pc <- .sx_reduce("mean", .cure_predict_X(x, x@X, x@Xc, 1, "cured", FALSE, 0.95))
  op <- graphics::par(cex = cex, mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op))
  inc <- .na_increments(x@time, x@status, 1L, x@entry)
  km_t <- c(0, inc$time); km_s <- c(1, .sx_np_curve_inc(inc$d, inc$Y, "km"))
  graphics::plot(km_t, km_s, type = "s", col = "grey40", ylim = c(0, 1), xlim = c(0, max(tg)), xlab = xlab, ylab = ylab,
                 main = if (is.null(main)) "Fitted survival (average over subjects) and Kaplan-Meier" else main, las = 1, bty = "l", ...)
  graphics::lines(c(0, tg), c(1, S), col = col, lwd = 2)
  graphics::abline(h = pc, lty = 2, col = col)
  graphics::legend("topright", c("Kaplan-Meier", "cure model", "cured fraction"), col = c("grey40", col, col),
                   lty = c(1, 1, 2), lwd = c(1, 2, 1), bty = "n")
  invisible(x)
}
