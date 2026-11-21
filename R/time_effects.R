#' Time-varying effects: detect them, choose their form and refit the model
#'
#' @description
#' Checks whether the effect of each covariate of a Cox, cause-specific Cox or Fine-Gray model
#' changes over follow-up (a violation of proportional hazards), compares several forms for how it
#' changes, refits the model with the chosen forms and reports the hazard ratio as a function of
#' time, in a table, in plain words and in figures. With \code{report()} the whole analysis goes into
#' one HTML file with the methods and references.
#'
#' \strong{Use it} after \code{\link{ph_test}} has flagged a covariate, or as a routine check of a Cox
#' or Fine-Gray model whose hazard ratios you intend to report: it shows how the effect changes, not
#' only that it does.
#'
#' \strong{It assists; it does not decide for you.} The screening tests, the comparison of the forms
#' and the final model are reported separately so that each step can be checked; you can fix the form
#' of any covariate with \code{fix}. Intervals of a model chosen from the data do not account for the
#' choice unless \code{bootstrap} is used.
#' @param fit a model fitted by \code{\link{cox}} (without \code{tt()} terms), by \code{\link{csc}} with
#'   the Cox baseline, or by \code{\link{fg}}. Strata, delayed entry, (start, stop] rows, case weights
#'   and clusters of the fit are used.
#' @param forms the forms of the log hazard ratio \eqn{\beta(t)} to compare: a character vector of
#'   built-in forms, or a list that can also contain named functions of time (your own forms, see
#'   Details). \code{"constant"} is always included.
#' @param terms the terms of the model to examine (as in the formula, for example \code{"age"} or
#'   \code{"education"}); default all. The other terms keep a constant effect.
#' @param fix a named list (or named character vector) \code{term = form} that sets the form of a
#'   term instead of choosing it, for example \code{fix = list(treatment = "piecewise")}. Use it to
#'   fit a pre-specified model or for a sensitivity analysis.
#' @param screen which terms are given a chance of a time-varying form: \code{"holm"} (default)
#'   those whose test of a constant effect has a Holm-adjusted p-value below \code{alpha};
#'   \code{"unadjusted"} those with p below \code{alpha}; \code{"all"} every term (the choice then
#'   rests on \code{criterion} alone).
#' @param alpha level of the screening tests.
#' @param criterion how the form is chosen among the candidates of a term: \code{"aic"} (default) or
#'   \code{"bic"} (penalty \eqn{\log} of the number of events; favours simpler forms).
#' @param transform time scale of the screening test, as in \code{\link{ph_test}}: \code{"km"}
#'   (default; one minus the Kaplan-Meier estimate, for a Fine-Gray model one minus the weighted
#'   product-limit estimate of the subdistribution), \code{"rank"}, \code{"identity"} or \code{"log"}.
#' @param df degrees of freedom of the spline form.
#' @param breaks times at which a step in the effect is allowed (form \code{"piecewise"}). Default
#'   \code{NULL}: one change point chosen from the data (see Details).
#' @param times times at which the hazard ratios are reported. Default: the 10th, 25th, 50th, 75th
#'   and 90th percentiles of the event times (two significant digits).
#' @param level confidence level.
#' @param bootstrap number of bootstrap samples of the subjects (clusters, if the fit has them) in
#'   which the whole procedure (screening, choice of the forms, final fit) is repeated. \code{0}
#'   (default): none. 200 or more gives percentile intervals that include the choice of the forms
#'   and the frequency with which each form is chosen.
#' @param seed seed of the bootstrap; the state of the random number generator of the session is
#'   restored afterwards.
#' @param max_rows largest number of rows of the data split at every event time. Beyond it, follow-up
#'   is split at quantiles of the event times instead, and the result says so.
#' @details
#' \subsection{Model}{
#' For a covariate column \eqn{x} the log hazard ratio (cause-specific hazard ratio for \code{csc},
#' subdistribution hazard ratio for \code{fg}) is \eqn{\beta(t) = \beta + \gamma^\top B(t)}, with
#' \eqn{B} one of the forms:
#' \describe{
#'   \item{\code{"constant"}}{\eqn{B = 0}: proportional hazards.}
#'   \item{\code{"log"}}{\eqn{B(t) = \log(t / t_m)}, \eqn{t_m} the median event time, so that
#'     \eqn{e^\beta} is the hazard ratio at \eqn{t_m}. The commonest pattern: an effect that fades
#'     or builds up gradually.}
#'   \item{\code{"linear"}}{\eqn{B(t) = t - t_m}; not in the default set.}
#'   \item{\code{"piecewise"}}{\eqn{B(t) = I(t > \tau)}: one hazard ratio before \eqn{\tau} and one
#'     after. Without \code{breaks}, \eqn{\tau} is chosen by profile likelihood among the 20th to 80th
#'     percentiles (in steps of 10) of the event times and counted as one extra parameter in the
#'     criterion and in the likelihood ratio test, which is then approximate. With \code{breaks}, one
#'     step at every break.}
#'   \item{\code{"spline"}}{a natural cubic spline in \eqn{\log t} with \code{df} degrees of freedom,
#'     knots at quantiles of the log event times (Hess 1994; Gray 1992): any smooth pattern.}
#'   \item{your own}{a named function of time returning a vector, or a matrix with one column per
#'     basis function, for example \code{forms = list("constant", "log", sqrt = function(t) sqrt(t))}.}
#' }
#' A factor gets one time-varying term per level (each with the same form). Every model is fitted by
#' partial likelihood on the data split at every event time, as \code{cox(..., tt = )} does
#' (Therneau and Grambsch 2000, Chapter 6), so the time functions are evaluated exactly at the event
#' times. For \code{fg} the same is done on the weighted counting-process data of Geskus (2011),
#' which reproduce the Fine-Gray estimator, with robust standard errors clustered on the subject
#' (Lin and Wei 1989).
#' }
#' \subsection{Procedure}{
#' For each analysis (the Cox model, every cause of a cause-specific model, or the Fine-Gray model):
#' \enumerate{
#'   \item \strong{Screening.} The score test of Grambsch and Therneau (1994) of a constant effect,
#'     for every term (all its columns together) and globally; the same test as \code{\link{ph_test}}.
#'     Holm's (1979) adjustment controls the chance of flagging any term when all effects are
#'     constant.
#'   \item \strong{Comparison.} For every term in \code{terms}, each form is fitted with the other
#'     terms constant, and compared with the constant form by AIC (or BIC) and by the likelihood
#'     ratio test. For terms that passed the screening (or are in \code{fix}), the form with the
#'     smallest criterion is chosen; for the others the effect stays constant. A term can pass the
#'     screening and still keep the constant form when no form improves the criterion.
#'   \item \strong{Final model.} All chosen forms are fitted together. If that model does not converge,
#'     the term whose form gains least is made constant and the model refitted (reported in a note).
#'   \item \strong{Effects.} \eqn{e^{\beta(t)}} with a Wald interval at \code{times}, and a curve over
#'     follow-up for the figures.
#'   \item \strong{Bootstrap} (optional). Steps 1 to 4 in every bootstrap sample.
#' }
#' }
#' \subsection{Combining with the rest of the package}{
#' \code{predict()} gives survival (\code{cox}) or cumulative incidence (\code{fg}; \code{csc}, from
#' the cause-specific hazards with product-limit overall survival) with the time-varying effects;
#' these predictions can be compared with other models through \code{\link{brier_score}} or
#' \code{\link{compare_models}} with a function that calls \code{time_effects()}.
#' \code{\link{sx_references}(x)} lists the methods used. A form found here can be refitted as an
#' ordinary model: \code{cox(Event(time, status) ~ x + tt(x) + z, d, tt = function(t) log(t))}.
#' }
#' @section Interpreting the result:
#' \describe{
#'   \item{Test of a constant effect}{a small p-value means that the hazard ratio of the term changes
#'     over time. With large data sets small, unimportant changes are significant: judge their size
#'     in the table of hazard ratios over time and in \code{plot(x)}.}
#'   \item{Comparison of forms}{\code{LR} and \code{p}: likelihood ratio test of the form against a
#'     constant effect. \code{dAIC}: difference from the smallest AIC of the term; differences below
#'     2 are small (Burnham and Anderson 2002). Forms within 2 of each other describe the data
#'     equally well; prefer the simpler or the one that suits the question (a step for an effect
#'     that is expected to stop, such as surgery).}
#'   \item{Hazard ratio over time}{the quantity to report. The reading says whether the effect weakens
#'     (moves towards 1), strengthens, or changes direction, and near which time.}
#'   \item{Stability}{with \code{bootstrap}: how often each form is chosen. A form chosen in, say,
#'     40\% of the samples is not well determined; report the hazard ratios at a few times (which are
#'     usually stable) rather than the form.}
#'   \item{Fine-Gray model}{the time-varying sHR describes the effect on the subdistribution hazard;
#'     it is not a ratio of risks. \code{predict()} gives the cumulative incidence.}
#' }
#' A sentence for a report: \emph{"The effect of treatment was not constant over time (Grambsch-Therneau
#' test, Holm-adjusted p = 0.002). A log-linear change in the log hazard ratio fitted best (AIC 6.1 lower
#' than a constant effect); the hazard ratio was 0.52 (95\% CI 0.38 to 0.71) at 6 months and 0.91 (0.67
#' to 1.24) at 3 years. The form was chosen in 81\% of 200 bootstrap samples."}
#' @section Assumptions and how to check them:
#' \describe{
#'   \item{The model apart from time}{the forms describe how effects change with time; the functional
#'     form of continuous covariates and the absence of interactions are still assumptions.}
#'   \item{Enough events}{every form spends parameters; with fewer than about 10 events per parameter
#'     the comparison is unreliable and the screening test has little power.}
#'   \item{Independent censoring}{as for the model that was fitted; for \code{fg} the censoring
#'     distribution has no covariates.}
#' }
#' @section Messages and common problems:
#' \describe{
#'   \item{\code{... was kept constant}}{the model with all chosen forms did not converge (often two
#'     flexible forms for correlated covariates); use \code{fix} or fewer forms.}
#'   \item{\code{follow-up was split at ... quantiles}}{the data split at every event time would have
#'     had more than \code{max_rows} rows; the forms are then evaluated at the quantiles, which is
#'     an approximation. Increase \code{max_rows} if memory allows.}
#'   \item{Long computing times}{each term and form is a separate fit on the split data, and the
#'     bootstrap repeats everything; examine fewer \code{terms} or use fewer \code{forms}.}
#' }
#' @section Validation:
#' The screening test is the test of \code{\link{ph_test}} and agrees with it to \eqn{10^{-8}}; the
#' fits of the log form agree with \code{cox(..., tt = )}, which agrees with
#' \code{survival::coxph(tt = )}; with every effect constant the Fine-Gray fit reproduces
#' \code{\link{fg}}, and the predictions reproduce those of \code{cox} and \code{fg}. These checks
#' are in the test suite of the package.
#' @return
#' An S7 object of class \code{TimeEffects}. \code{x@units} has one element per analysis with
#' \code{screening} (one row per term: \code{chisq}, \code{df}, \code{p}, \code{p_holm},
#' \code{considered}), \code{comparison} (one row per term and form: \code{loglik}, \code{k} parameters,
#' \code{ic}, \code{delta}, \code{lr}, \code{lr_df}, \code{p_lr}, \code{converged}, \code{chosen}),
#' \code{selected} (the form of every term), \code{final} (coefficients, variance, likelihood ratio
#' test of all time-varying terms), \code{effects} (hazard ratios at \code{times}: \code{column},
#' \code{term}, \code{form}, \code{time}, \code{log_hr}, \code{se}, \code{hr}, \code{lower},
#' \code{upper}, \code{p}) and \code{notes}. \code{x@bootstrap}: selection frequencies and percentile
#' intervals.
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(x)} \tab tests, the comparison for the terms with a time-varying form, hazard ratios
#'     over time and a reading \cr
#'   \code{summary(x)} \tab also the comparison for every screened term and the constant effects \cr
#'   \code{plot(x, term = NULL, cause = NULL, residuals = TRUE, scale = "ratio")} \tab \eqn{e^{\beta(t)}}
#'     with its band, the constant fit and the scaled Schoenfeld residuals \cr
#'   \code{predict(x, newdata, times, cause = 1)} \tab survival or cumulative incidence \cr
#'   \code{as.data.frame(x, what)} \tab \code{"effects"}, \code{"screening"}, \code{"forms"},
#'     \code{"coefficients"} or \code{"stability"} \cr
#'   \code{report(x, file)} \tab a self-contained HTML report \cr
#'   \code{\link{sx_references}(x)} \tab the references of the methods \cr
#' }
#' @references
#' \emph{Test of proportional hazards and time-varying coefficients}
#'
#' Grambsch PM, Therneau TM (1994) Proportional hazards tests and diagnostics based on weighted
#' residuals. Biometrika 81:515-526.
#'
#' Therneau TM, Grambsch PM (2000) Modeling Survival Data: Extending the Cox Model. Springer,
#' Chapter 6.
#'
#' Hess KR (1994) Assessing time-by-covariate interactions in proportional hazards regression models
#' using cubic spline functions. Stat Med 13:1045-1062.
#'
#' Gray RJ (1992) Flexible methods for analyzing survival data using splines, with applications to
#' breast cancer prognosis. J Am Stat Assoc 87:942-951.
#'
#' Buchholz A, Sauerbrei W (2011) Comparison of procedures to assess non-linear and time-varying
#' effects in multivariable models for survival data. Biom J 53:308-331.
#'
#' \emph{Multiple testing and model choice}
#'
#' Holm S (1979) A simple sequentially rejective multiple test procedure. Scand J Stat 6:65-70.
#'
#' Akaike H (1974) A new look at the statistical model identification. IEEE Trans Automat Contr
#' 19:716-723.
#'
#' Schwarz G (1978) Estimating the dimension of a model. Ann Stat 6:461-464.
#'
#' Burnham KP, Anderson DR (2002) Model Selection and Multimodel Inference, 2nd ed. Springer.
#'
#' \emph{Fine-Gray model and variance}
#'
#' Fine JP, Gray RJ (1999) A proportional hazards model for the subdistribution of a competing risk.
#' J Am Stat Assoc 94:496-509.
#'
#' Geskus RB (2011) Cause-specific cumulative incidence estimation and the Fine and Gray model under
#' both left truncation and right censoring. Biometrics 67:39-49.
#'
#' Lin DY, Wei LJ (1989) The robust inference for the Cox proportional hazards model. J Am Stat Assoc
#' 84:1074-1078.
#'
#' Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap. Chapman and Hall.
#'
#' All references: \code{sx_references("time_effects")}.
#' @seealso
#' \code{\link{ph_test}} (the test and residual plot), \code{\link{cox}} (\code{tt()} and
#' \code{strata()}), \code{\link{tt_effect}}, \code{\link{fg}}, \code{\link{csc}},
#' \code{\link{model_scan}}, \code{\link{rmst}} (a summary that needs no proportional hazards).
#' @examples
#' set.seed(4)
#' n <- 400
#' x <- rbinom(n, 1, 0.5); z <- rnorm(n)
#' t <- rweibull(n, ifelse(x == 1, 0.7, 1.4), 3 * exp(-0.3 * z)); cens <- runif(n, 0, 6)
#' d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = z)
#'
#' fit <- cox(Event(time, status) ~ x + z, d)
#' te <- time_effects(fit)
#' te                                   # x: time-varying, z: constant
#' plot(te)
#' as.data.frame(te)                    # hazard ratios over time
#' predict(te, data.frame(x = 0:1, z = 0), times = c(1, 2, 4))
#'
#' ## a pre-specified form, and your own function of time
#' time_effects(fit, fix = list(x = "piecewise"), breaks = 1)
#' time_effects(fit, forms = list("constant", "log", sqrt = function(t) sqrt(t)))
#'
#' \donttest{
#' ## stability of the choice and intervals that include it
#' te_b <- time_effects(fit, forms = c("constant", "log"), bootstrap = 50)
#' summary(te_b)
#' report(te_b, file = tempfile(fileext = ".html"))
#'
#' ## Fine-Gray model: time-varying subdistribution hazard ratios
#' data(fertility)
#' time_effects(fg(Event(time, status) ~ education, fertility, cause = 1))
#' }
#' @export
time_effects <- function(fit, forms = c("constant", "log", "piecewise", "spline"), terms = NULL, fix = NULL,
                         screen = c("holm", "unadjusted", "all"), alpha = 0.05, criterion = c("aic", "bic"),
                         transform = c("km", "rank", "identity", "log"), df = 3L, breaks = NULL, times = NULL,
                         level = 0.95, bootstrap = 0L, seed = 1L, max_rows = 5e6) {
  cl <- match.call()
  screen <- match.arg(screen); criterion <- match.arg(criterion); transform <- match.arg(transform)
  if (!is.numeric(alpha) || length(alpha) != 1L || !(alpha > 0 && alpha < 1)) stop("alpha must be a number between 0 and 1", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !(level > 0 && level < 1)) stop("level must be a number between 0 and 1", call. = FALSE)
  df <- as.integer(df)
  if (length(df) != 1L || is.na(df) || df < 1L || df > 8L) stop("df must be an integer from 1 to 8", call. = FALSE)
  bootstrap <- as.integer(bootstrap)
  if (length(bootstrap) != 1L || is.na(bootstrap) || bootstrap < 0L) stop("bootstrap must be 0 or a number of resamples", call. = FALSE)
  if (bootstrap > 0L && bootstrap < 20L) stop("bootstrap: use at least 20 resamples (200 or more for intervals)", call. = FALSE)
  if (!is.null(times) && (!is.numeric(times) || any(!is.finite(times)) || any(times <= 0)))
    stop("times must be positive numbers", call. = FALSE)
  if (!is.null(breaks) && (!is.numeric(breaks) || any(!is.finite(breaks)) || any(breaks <= 0)))
    stop("breaks must be positive numbers (times at which the effect may change)", call. = FALSE)
  raw <- .te_raw(fit)
  if (!is.null(fix)) fix <- as.list(fix)
  spec <- list(forms = .te_form_specs(forms), fix = fix, screen = screen, alpha = alpha, criterion = criterion,
               transform = transform, df = df, breaks = breaks, level = level, max_rows = max_rows)
  spec$terms <- .te_terms_arg(terms, raw)
  .te_check_fix(fix, raw, spec$forms)
  units <- lapply(raw$units, function(k) .te_unit(raw, k, spec, times = times, keep = TRUE))
  names(units) <- vapply(units, `[[`, character(1), "label")
  boot <- NULL
  if (bootstrap > 0L) boot <- .te_bootstrap(raw, spec, units, bootstrap, seed)
  TimeEffects(model_call = cl, fit = fit, kind = raw$kind, units = units, settings = spec[setdiff(names(spec), "forms")],
              forms = vapply(spec$forms, `[[`, character(1), "name"), bootstrap = boot, level = level)
}

# ----------------------------------------------------------------------------------------------
# input: the fitted model
# ----------------------------------------------------------------------------------------------
.te_raw <- function(fit) {
  if (S7::S7_inherits(fit, CoxFit)) {
    if (!is.null(fit@tt)) stop("fit already has tt() terms; give time_effects() the model without them", call. = FALSE)
    kind <- "cox"; units <- NA_integer_; ties <- fit@ties; w <- fit@weights; robust <- isTRUE(fit@robust)
  } else if (S7::S7_inherits(fit, CSCFit)) {
    kind <- "csc"; units <- seq_len(fit@K); ties <- fit@ties; w <- NULL; robust <- isTRUE(fit@robust)
  } else if (S7::S7_inherits(fit, FineGrayFit)) {
    kind <- "fg"; units <- fit@cause; ties <- "breslow"; w <- NULL; robust <- TRUE
  } else {
    stop("time_effects() takes a model fitted by cox(), by csc() with a Cox baseline, or by fg()", call. = FALSE)
  }
  X <- fit@X
  asg <- attr(X, "assign")
  labs <- attr(fit@terms, "term.labels")
  if (is.null(asg) || length(asg) != ncol(X) || is.null(labs) || any(asg < 1L) || any(asg > length(labs))) {
    asg <- seq_len(ncol(X)); labs <- colnames(X)
  }
  term_of <- labs[asg]
  if (any(!is.finite(fit@time)) || any(fit@time <= 0))
    stop("time_effects() needs positive event and censoring times", call. = FALSE)
  strata <- if (kind == "fg") NULL else fit@strata
  cluster <- if (kind == "fg") NULL else fit@cluster
  list(kind = kind, fit = fit, units = units, time = fit@time, status = fit@status, X = X, term_of = term_of,
       entry = fit@entry, w = w, strata = strata, cluster = cluster, ties = ties,
       robust = robust, n = length(fit@time))
}

.te_terms_arg <- function(terms, raw) {
  all <- unique(raw$term_of)
  if (is.null(terms)) return(all)
  bad <- setdiff(terms, all)
  if (length(bad)) stop("terms not in the model: ", paste(bad, collapse = ", "), "; available: ", paste(all, collapse = ", "), call. = FALSE)
  terms
}

.te_check_fix <- function(fix, raw, forms) {
  if (is.null(fix)) return(invisible(TRUE))
  if (!is.list(fix) && !is.character(fix)) stop("fix must be a named list or character vector: term = form", call. = FALSE)
  fix <- as.list(fix)
  if (is.null(names(fix)) || any(names(fix) == "")) stop("fix must be named by term, for example fix = list(age = \"log\")", call. = FALSE)
  bad <- setdiff(names(fix), unique(raw$term_of))
  if (length(bad)) stop("fix: terms not in the model: ", paste(bad, collapse = ", "), call. = FALSE)
  fn <- vapply(forms, `[[`, character(1), "name")
  badf <- setdiff(unlist(fix), fn)
  if (length(badf)) stop("fix: forms not among 'forms': ", paste(badf, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

# ----------------------------------------------------------------------------------------------
# forms of a time-varying effect: the log hazard ratio of column x is beta + gamma' B(t)
# ----------------------------------------------------------------------------------------------
.te_builtin <- c("constant", "log", "linear", "piecewise", "spline")

.te_form_specs <- function(forms) {
  if (is.function(forms)) stop("forms: give a list, for example list(\"constant\", sqrt = function(t) sqrt(t))", call. = FALSE)
  f <- as.list(forms)
  nm <- names(f) %||% rep("", length(f))
  out <- list()
  for (i in seq_along(f)) {
    v <- f[[i]]
    if (is.character(v) && length(v) == 1L) {
      if (!v %in% .te_builtin) stop("unknown form '", v, "'; built-in forms: ", paste(.te_builtin, collapse = ", "),
                                    "; your own forms are functions of time, for example sqrt = function(t) sqrt(t)", call. = FALSE)
      out[[length(out) + 1L]] <- list(name = v, type = v, fun = NULL)
    } else if (is.function(v)) {
      if (!nzchar(nm[i])) stop("forms: a function of time needs a name, for example list(sqrt = function(t) sqrt(t))", call. = FALSE)
      if (nm[i] %in% .te_builtin) stop("forms: '", nm[i], "' is the name of a built-in form; choose another name", call. = FALSE)
      out[[length(out) + 1L]] <- list(name = nm[i], type = "user", fun = v)
    } else stop("forms: every element must be the name of a built-in form or a named function of time", call. = FALSE)
  }
  nms <- vapply(out, `[[`, character(1), "name")
  if (anyDuplicated(nms)) stop("forms: names must be unique", call. = FALSE)
  if (!"constant" %in% nms) out <- c(list(list(name = "constant", type = "constant", fun = NULL)), out)
  out
}

# turn a form specification into a basis for the event times et of one unit: list with the basis
# function B(t) (a matrix with q columns), q, the number of parameters spent on choosing it
# (a change point estimated from the data), a label and a description for the reader
.te_basis <- function(fs, et, df, breaks, brk = NULL) {
  tpos <- function(t) pmax(t, min(et) / 2)
  tm <- stats::median(et)
  switch(fs$type,
    constant = list(name = "constant", q = 0L, extra = 0L, B = NULL, label = "constant",
                    describe = "constant over time (proportional hazards)"),
    log = list(name = fs$name, q = 1L, extra = 0L, B = function(t) cbind(log(tpos(t) / tm)), tmed = tm,
               label = "log-linear in time", describe = paste0("log hazard ratio changes linearly in log(time); centred at t = ", format(signif(tm, 3L)))),
    linear = list(name = fs$name, q = 1L, extra = 0L, B = function(t) cbind(t - tm), tmed = tm,
                  label = "linear in time", describe = paste0("log hazard ratio changes linearly in time; centred at t = ", format(signif(tm, 3L)))),
    piecewise = {
      b <- brk %||% breaks
      b <- sort(unique(b))
      list(name = fs$name, q = length(b), extra = if (is.null(breaks)) 1L else 0L,
           B = function(t) matrix(vapply(b, function(s) as.numeric(t > s), numeric(length(t))), nrow = length(t)),
           breaks = b, label = if (length(b) == 1L) paste0("step at t = ", format(signif(b, 3L))) else "steps",
           describe = paste0("one hazard ratio before and one after t = ", paste(format(signif(b, 3L)), collapse = ", "),
                             if (is.null(breaks)) " (change point chosen from the data)" else ""))
    },
    spline = {
      lt <- log(et)
      kn <- if (df > 1L) stats::quantile(lt, seq_len(df - 1L) / df, names = FALSE, type = 7) else numeric(0)
      bk <- range(lt)
      list(name = fs$name, q = df, extra = 0L,
           B = function(t) { b <- splines::ns(log(tpos(t)), knots = kn, Boundary.knots = bk, intercept = FALSE)
                             matrix(as.numeric(b), nrow = length(t)) },
           knots = kn, boundary = bk, label = paste0("spline in log time (", df, " df)"),
           describe = paste0("smooth change: natural cubic spline in log(time) with ", df, " degrees of freedom"))
    },
    user = {
      g <- fs$fun
      test <- g(et)
      m <- if (is.null(dim(test))) matrix(test, ncol = 1L) else as.matrix(test)
      if (nrow(m) != length(et) || any(!is.finite(m)))
        stop("form '", fs$name, "': the function must return one finite value (or row) per time; log(t) needs t > 0", call. = FALSE)
      list(name = fs$name, q = ncol(m), extra = 0L,
           B = function(t) { r <- g(t); if (is.null(dim(r))) matrix(r, ncol = 1L) else as.matrix(r) },
           label = paste0("user form '", fs$name, "'"), describe = paste0("your function of time '", fs$name, "'"))
    })
}

# ----------------------------------------------------------------------------------------------
# one analysis unit: the Cox model, one cause of a cause-specific model, or the Fine-Gray model
# ----------------------------------------------------------------------------------------------
.te_rows <- function(raw, k) {
  n <- raw$n
  if (raw$kind == "fg") {
    e <- .fg_expand(raw$time, raw$status, raw$entry, k)
    return(list(start = e$start, stop = e$stop, status = as.integer(e$status), w = e$w, strata = NULL, id = e$id,
                X = raw$X[e$id, , drop = FALSE], cluster = e$id, has_entry = TRUE))
  }
  st <- if (raw$kind == "csc") as.integer(raw$status == k) else as.integer(raw$status)
  list(start = raw$entry %||% rep(0, n), stop = raw$time, status = st, w = raw$w, strata = raw$strata, id = seq_len(n),
       X = raw$X, cluster = raw$cluster, has_entry = !is.null(raw$entry))
}

.te_cox <- function(start, stop, status, X, w, ties, strata, entry_used) {
  r <- tryCatch(.cox_fit(stop, status, X, w, ties, 50L, 1e-9, if (entry_used) start else NULL, strata),
                error = function(e) NULL)
  if (is.null(r) || !isTRUE(r$converged) || !is.finite(r$loglik) || !all(is.finite(r$coefficients))) return(NULL)
  r
}

# cut every row (start, stop] at the event times et strictly inside it (as cox() does for tt())
.te_split <- function(start, stop, status, et) {
  lo <- findInterval(start, et)
  hi <- findInterval(stop, et, left.open = TRUE)
  rows <- pmax(hi - lo, 0L) + 1L
  idx <- rep.int(seq_along(stop), rows)
  r <- sequence(rows)
  ci <- lo[idx] + r
  last <- r == rows[idx]
  list(idx = idx, stop = ifelse(last, stop[idx], et[pmin(ci, length(et))]),
       start = ifelse(r == 1L, start[idx], et[pmax(ci - 1L, 1L)]), status = ifelse(last, status[idx], 0L))
}

.te_split_size <- function(start, stop, et) sum(as.numeric(pmax(findInterval(stop, et, left.open = TRUE) - findInterval(start, et), 0L) + 1L))

# function of time used by the screening test at the event rows of the risk-set sums
.te_gfun <- function(transform, time, status, entry, w, weighted) {
  force(transform)
  function(rows, tev) {
    if (transform == "km" && weighted) {
      et <- sort(unique(time[status == 1L]))
      d <- vapply(et, function(s) sum(w[status == 1L & time == s]), numeric(1))
      Y <- vapply(et, function(s) sum(w[time >= s & (if (is.null(entry)) TRUE else entry < s)]), numeric(1))
      S <- cumprod(1 - d / Y)
      return(1 - c(1, S)[match(tev, et)])
    }
    .time_transform(time, status, tev, transform, entry, rows)
  }
}

# score test of a constant effect (Grambsch & Therneau 1994), per column, per term and globally
.te_zph <- function(time, status, X, beta, vcov, gfun, w, entry, strata, term_of) {
  rs <- .cox_risksets(time, status, X, beta, w, entry, strata)
  sch <- X[rs$rows, , drop = FALSE] - rs$xbar
  g0 <- gfun(rs$rows, rs$time)
  g <- g0 - mean(g0)
  p <- ncol(X); nd <- length(rs$rows)
  U <- colSums(sch * g)
  Ibb <- apply(rs$V, c(2L, 3L), sum)
  Ibg <- apply(rs$V * g, c(2L, 3L), sum)
  Igg <- apply(rs$V * g^2, c(2L, 3L), sum)
  A <- Igg - t(Ibg) %*% solve(Ibb, Ibg)
  chi_col <- vapply(seq_len(p), function(j) U[j]^2 / A[j, j], numeric(1))
  terms <- unique(term_of)
  chi_term <- vapply(terms, function(tm) {
    j <- which(term_of == tm)
    drop(t(U[j]) %*% solve(A[j, j, drop = FALSE], U[j]))
  }, numeric(1))
  df_term <- vapply(terms, function(tm) sum(term_of == tm), integer(1))
  glob <- drop(t(U) %*% solve(A, U))
  scaled <- sch %*% vcov * nd + matrix(beta, nd, p, byrow = TRUE)
  colnames(scaled) <- colnames(X)
  list(columns = data.frame(column = colnames(X), term = term_of, chisq = chi_col, df = 1L,
                            p = stats::pchisq(chi_col, 1L, lower.tail = FALSE), stringsAsFactors = FALSE),
       terms = data.frame(term = terms, chisq = unname(chi_term), df = unname(df_term),
                          p = stats::pchisq(unname(chi_term), unname(df_term), lower.tail = FALSE), stringsAsFactors = FALSE),
       global = c(chisq = glob, df = p, p = stats::pchisq(glob, p, lower.tail = FALSE)),
       resid = list(time = rs$time, g = g0, y = scaled))
}

# design of the split data with the chosen form of every column: base columns, then for every column
# with a time-varying form the products x * B(t); `map` says where each column's coefficients are
.te_design <- function(Xs, stop, term_of, sel, bases) {
  cols <- list(Xs); map <- list(); off <- ncol(Xs)
  for (j in seq_len(ncol(Xs))) {
    tm <- term_of[j]; f <- sel[[tm]] %||% "constant"
    b <- bases[[tm]][[f]] %||% NULL
    if (is.null(b) || b$q == 0L) { map[[j]] <- list(col = j, form = "constant", tv = integer(0), basis = NULL); next }
    Bm <- b$B(stop)
    cols[[length(cols) + 1L]] <- Xs[, j] * Bm
    map[[j]] <- list(col = j, form = f, tv = off + seq_len(ncol(Bm)), basis = b)
    off <- off + ncol(Bm)
  }
  X <- do.call(cbind, cols)
  nm <- colnames(Xs)
  for (j in seq_along(map)) if (length(map[[j]]$tv))
    nm <- c(nm, paste0(colnames(Xs)[j], ":", map[[j]]$form, if (length(map[[j]]$tv) > 1L) seq_along(map[[j]]$tv) else ""))
  colnames(X) <- nm
  list(X = X, map = map)
}

.te_ic <- function(ll, k, d, criterion) if (criterion == "aic") -2 * ll + 2 * k else -2 * ll + log(d) * k

.te_unit <- function(raw, k, spec, times = NULL, keep = TRUE) {
  u <- .te_rows(raw, k)
  label <- switch(raw$kind, cox = "the event", csc = paste("cause", k), fg = paste("cause", k, "(subdistribution)"))
  term_of <- raw$term_of
  terms <- unique(term_of)
  w1 <- u$w %||% rep(1, length(u$stop))
  nev <- sum(u$status == 1L & w1 > 0)
  if (nev < 5L) stop(label, ": too few events (", nev, ")", call. = FALSE)
  # constant model on the unit's rows
  r0 <- .te_cox(u$start, u$stop, u$status, u$X, u$w, raw$ties, u$strata, u$has_entry)
  if (is.null(r0)) stop(label, ": the model with constant effects did not converge", call. = FALSE)
  # screening
  gf <- .te_gfun(spec$transform, u$stop, u$status, if (u$has_entry) u$start else NULL, w1,
                 weighted = raw$kind == "fg")
  zp <- tryCatch(.te_zph(u$stop, u$status, u$X, r0$coefficients, r0$vcov, gf, u$w, if (u$has_entry) u$start else NULL,
                         u$strata, term_of), error = function(e) NULL)
  scr <- if (is.null(zp)) data.frame(term = terms, chisq = NA_real_, df = NA_integer_, p = NA_real_) else zp$terms
  scr$p_holm <- if (all(is.na(scr$p))) NA_real_ else stats::p.adjust(scr$p, "holm")
  scr$considered <- switch(spec$screen,
                           holm = !is.na(scr$p_holm) & scr$p_holm < spec$alpha,
                           unadjusted = !is.na(scr$p) & scr$p < spec$alpha,
                           all = rep(TRUE, nrow(scr)))
  scr$considered <- scr$considered & scr$term %in% spec$terms
  fixed <- spec$fix %||% list()
  scr$considered[scr$term %in% names(fixed)] <- TRUE
  # split at the event times (a grid of event-time quantiles if that would be too large)
  et <- sort(unique(u$stop[u$status == 1L & w1 > 0]))
  grid <- FALSE
  if (.te_split_size(u$start, u$stop, et) > spec$max_rows) {
    m <- 400L
    repeat {
      eg <- unique(stats::quantile(et, seq(0, 1, length.out = m), type = 1, names = FALSE))
      if (.te_split_size(u$start, u$stop, eg) <= spec$max_rows || m <= 25L) break
      m <- m %/% 2L
    }
    grid <- TRUE; cut <- eg
  } else cut <- et
  sp <- .te_split(u$start, u$stop, u$status, cut)
  S <- list(start = sp$start, stop = sp$stop, status = sp$status, X = u$X[sp$idx, , drop = FALSE],
            w = if (is.null(u$w)) NULL else u$w[sp$idx], strata = if (is.null(u$strata)) NULL else u$strata[sp$idx],
            cluster = if (is.null(u$cluster)) (if (raw$robust) u$id[sp$idx] else NULL) else u$cluster[sp$idx])
  fitS <- function(Xd) .te_cox(S$start, S$stop, S$status, Xd, S$w, raw$ties, S$strata, TRUE)
  base <- fitS(S$X)
  if (is.null(base)) stop(label, ": the model could not be fitted on the split data", call. = FALSE)
  # bases of every form for every term (piecewise without breaks: change point by profile likelihood)
  bases <- list(); comp <- list()
  for (tm in terms) {
    bases[[tm]] <- list()
    cand <- scr$considered[scr$term == tm]
    if (!cand && (!keep || !tm %in% spec$terms)) next
    for (fs in spec$forms) {
      if (fs$type == "constant") {
        bases[[tm]][["constant"]] <- .te_basis(fs, et, spec$df, spec$breaks)
        comp[[length(comp) + 1L]] <- data.frame(term = tm, form = "constant", label = "constant", loglik = base$loglik,
                                                npar = length(base$coefficients), extra = 0L, converged = TRUE, stringsAsFactors = FALSE)
        next
      }
      if (fs$type == "piecewise" && is.null(spec$breaks)) {
        cands <- unique(stats::quantile(et, seq(0.2, 0.8, by = 0.1), type = 1, names = FALSE))
        best <- NULL
        for (b0 in cands) {
          bb <- .te_basis(fs, et, spec$df, NULL, brk = b0)
          dsg <- .te_design(S$X, S$stop, term_of, stats::setNames(list(fs$name), tm), stats::setNames(list(stats::setNames(list(bb), fs$name)), tm))
          r <- fitS(dsg$X)
          if (!is.null(r) && (is.null(best) || r$loglik > best$r$loglik)) best <- list(r = r, b = bb)
        }
        if (is.null(best)) { comp[[length(comp) + 1L]] <- data.frame(term = tm, form = fs$name, label = "step", loglik = NA_real_,
                                                                      npar = NA_integer_, extra = 1L, converged = FALSE, stringsAsFactors = FALSE); next }
        bases[[tm]][[fs$name]] <- best$b
        comp[[length(comp) + 1L]] <- data.frame(term = tm, form = fs$name, label = best$b$label, loglik = best$r$loglik,
                                                npar = length(best$r$coefficients), extra = 1L, converged = TRUE, stringsAsFactors = FALSE)
        next
      }
      bb <- tryCatch(.te_basis(fs, et, spec$df, spec$breaks), error = function(e) e)
      if (inherits(bb, "error")) stop(conditionMessage(bb), call. = FALSE)
      dsg <- .te_design(S$X, S$stop, term_of, stats::setNames(list(fs$name), tm), stats::setNames(list(stats::setNames(list(bb), fs$name)), tm))
      r <- fitS(dsg$X)
      bases[[tm]][[fs$name]] <- bb
      comp[[length(comp) + 1L]] <- data.frame(term = tm, form = fs$name, label = bb$label, loglik = if (is.null(r)) NA_real_ else r$loglik,
                                              npar = if (is.null(r)) NA_integer_ else length(r$coefficients), extra = bb$extra,
                                              converged = !is.null(r), stringsAsFactors = FALSE)
    }
  }
  comp <- if (length(comp)) do.call(rbind, comp) else
    data.frame(term = character(0), form = character(0), label = character(0), loglik = numeric(0), npar = integer(0),
               extra = integer(0), converged = logical(0), stringsAsFactors = FALSE)
  comp$k <- comp$npar + comp$extra
  comp$ic <- .te_ic(comp$loglik, comp$k, nev, spec$criterion)
  ll0 <- base$loglik
  comp$lr <- 2 * (comp$loglik - ll0)
  comp$lr_df <- comp$k - length(base$coefficients)
  comp$p_lr <- ifelse(comp$lr_df > 0 & is.finite(comp$lr), stats::pchisq(pmax(comp$lr, 0), pmax(comp$lr_df, 1), lower.tail = FALSE), NA_real_)
  comp$delta <- NA_real_
  sel <- list()
  for (tm in terms) {
    i <- which(comp$term == tm)
    if (!length(i)) { sel[[tm]] <- "constant"; next }
    comp$delta[i] <- comp$ic[i] - min(comp$ic[i], na.rm = TRUE)
    if (!is.null(fixed[[tm]])) { sel[[tm]] <- fixed[[tm]]; next }
    if (!scr$considered[scr$term == tm]) { sel[[tm]] <- "constant"; next }
    ok <- i[comp$converged[i] & is.finite(comp$ic[i])]
    sel[[tm]] <- if (length(ok)) comp$form[ok][which.min(comp$ic[ok])] else "constant"
  }
  comp$chosen <- comp$form == vapply(comp$term, function(tm) sel[[tm]], character(1))
  # final model with every chosen form; terms whose form makes the joint model fail fall back to constant
  notes <- character(0)
  repeat {
    dsg <- .te_design(S$X, S$stop, term_of, sel, bases)
    fin <- fitS(dsg$X)
    if (!is.null(fin)) break
    tv <- names(sel)[vapply(sel, function(f) f != "constant", logical(1)) & !names(sel) %in% names(fixed)]
    if (!length(tv)) stop(label, ": the final model did not converge", call. = FALSE)
    gain <- vapply(tv, function(tm) {
      i <- comp$term == tm
      comp$ic[i & comp$form == "constant"][1L] - comp$ic[i & comp$form == sel[[tm]]][1L]
    }, numeric(1))
    drop1 <- tv[which.min(gain)]
    notes <- c(notes, paste0("the joint model with all chosen forms did not converge; '", drop1, "' was kept constant"))
    sel[[drop1]] <- "constant"
  }
  comp$chosen <- comp$form == vapply(comp$term, function(tm) sel[[tm]], character(1))
  vm <- fin$vcov
  if (!is.null(S$cluster)) fin$vcov <- .cox_robust(fin, S$stop, S$status, dsg$X, S$w, raw$ties, S$start, S$strata, S$cluster)
  any_tv <- any(vapply(sel, function(f) f != "constant", logical(1)))
  lr_all <- 2 * (fin$loglik - ll0); lr_df <- length(fin$coefficients) - length(base$coefficients)
  if (raw$kind == "fg") notes <- c(notes, paste0("the Fine-Gray log-likelihood is a weighted pseudo-likelihood: ", toupper(spec$criterion),
                                                  " and likelihood ratio tests are guides, not exact tests"))
  if (grid) notes <- c(notes, paste0("follow-up was split at ", length(cut), " quantiles of the event times instead of every event time ",
                                    "(more than max_rows rows); the forms are evaluated at these points"))
  # effects at chosen times and on a grid for plots
  tq <- sort(unique(times %||% signif(stats::quantile(et, c(0.1, 0.25, 0.5, 0.75, 0.9), type = 1, names = FALSE), 2L)))
  tg <- seq(min(et), max(et), length.out = 101L)
  eff <- .te_effects(fin, fin$vcov, dsg$map, term_of, tq, spec$level)
  curve <- .te_effects(fin, fin$vcov, dsg$map, term_of, tg, spec$level)
  out <- list(label = label, cause = k, kind = raw$kind, n = raw$n, nevent = nev, term_of = term_of, screening = scr,
              screening_columns = if (is.null(zp)) NULL else zp$columns, global = if (is.null(zp)) NULL else zp$global,
              comparison = comp, selected = unlist(sel), final = list(coefficients = fin$coefficients, vcov = fin$vcov,
              vcov_model = vm, loglik = fin$loglik, loglik_constant = ll0, lr = lr_all, lr_df = lr_df,
              p = if (lr_df > 0) stats::pchisq(max(lr_all, 0), lr_df, lower.tail = FALSE) else NA_real_,
              robust = !is.null(S$cluster), map = dsg$map),
              effects = eff, curve = curve, times = tq, event_times = et, any_tv = any_tv, notes = notes,
              constant = list(coefficients = base$coefficients, vcov = base$vcov), grid = grid)
  if (keep) {
    out$resid <- if (is.null(zp)) NULL else zp$resid
    out$split <- list(start = S$start, stop = S$stop, status = S$status, X = dsg$X, w = S$w, strata = S$strata, ties = raw$ties)
  }
  out
}

# log hazard ratio of every column at times tt: beta + gamma' B(t), with the delta-method interval
.te_effects <- function(fin, V, map, term_of, tt, level) {
  b <- fin$coefficients; z <- .zq(level)
  out <- lapply(map, function(m) {
    j <- m$col
    if (!length(m$tv)) {
      est <- rep(b[[j]], length(tt)); se <- rep(sqrt(V[j, j]), length(tt))
    } else {
      Bm <- m$basis$B(tt)
      a <- cbind(1, Bm)
      ix <- c(j, m$tv)
      est <- drop(a %*% b[ix])
      se <- sqrt(pmax(rowSums((a %*% V[ix, ix, drop = FALSE]) * a), 0))
    }
    data.frame(column = names(b)[j], term = term_of[j], form = m$form, time = tt, log_hr = est, se = se,
               hr = exp(est), lower = exp(est - z * se), upper = exp(est + z * se),
               p = 2 * stats::pnorm(-abs(est / se)), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# ----------------------------------------------------------------------------------------------
# bootstrap of the whole procedure: screening, choice of the forms and final fit
# ----------------------------------------------------------------------------------------------
.te_bootstrap <- function(raw, spec, units, B, seed) {
  cl <- raw$cluster %||% seq_len(raw$n)
  ids <- unique(cl)
  rows_of <- split(seq_along(cl), match(cl, ids))
  res <- .with_seed(seed, lapply(seq_len(B), function(b) {
    pick <- sample.int(length(ids), length(ids), replace = TRUE)
    rows <- unlist(rows_of[pick], use.names = FALSE)
    rb <- raw
    rb$time <- raw$time[rows]; rb$status <- raw$status[rows]; rb$X <- raw$X[rows, , drop = FALSE]
    rb$entry <- if (is.null(raw$entry)) NULL else raw$entry[rows]
    rb$w <- if (is.null(raw$w)) NULL else raw$w[rows]
    rb$strata <- if (is.null(raw$strata)) NULL else raw$strata[rows]
    rb$cluster <- if (is.null(raw$cluster)) NULL else rep(seq_along(pick), lengths(rows_of[pick]))
    rb$n <- length(rows)
    lapply(seq_along(units), function(q) {
      r <- tryCatch(.te_unit(rb, raw$units[q], spec, times = units[[q]]$times, keep = FALSE), error = function(e) NULL)
      if (is.null(r)) return(NULL)
      list(selected = r$selected, effects = r$effects)
    })
  }))
  z <- 1 - (1 - spec$level) / 2
  lapply(seq_along(units), function(q) {
    rr <- Filter(Negate(is.null), lapply(res, `[[`, q))
    nb <- length(rr)
    if (nb < max(10L, B / 2)) return(list(B = B, ok = nb, freq = NULL, ci = NULL))
    terms <- names(units[[q]]$selected)
    forms <- unique(unlist(lapply(rr, `[[`, "selected")))
    freq <- do.call(rbind, lapply(terms, function(tm) {
      s <- vapply(rr, function(r) r$selected[[tm]], character(1))
      data.frame(term = tm, form = forms, share = vapply(forms, function(f) mean(s == f), numeric(1)), stringsAsFactors = FALSE)
    }))
    freq <- freq[freq$share > 0, , drop = FALSE]
    e0 <- units[[q]]$effects
    key <- paste(e0$column, e0$time)
    L <- vapply(rr, function(r) r$effects$log_hr[match(key, paste(r$effects$column, r$effects$time))], numeric(nrow(e0)))
    L <- matrix(L, nrow = nrow(e0))
    ci <- data.frame(column = e0$column, time = e0$time,
                     lower = exp(apply(L, 1L, stats::quantile, probs = 1 - z, na.rm = TRUE, names = FALSE)),
                     upper = exp(apply(L, 1L, stats::quantile, probs = z, na.rm = TRUE, names = FALSE)), stringsAsFactors = FALSE)
    list(B = B, ok = nb, freq = freq, ci = ci)
  })
}

# ----------------------------------------------------------------------------------------------
# reading the result
# ----------------------------------------------------------------------------------------------
.te_ratio_name <- function(kind) switch(kind, cox = "HR", csc = "cause-specific HR", fg = "sHR")

.te_pretty <- function(x, s) .pretty_terms(s, x@fit@xlevels)

# one sentence per time-varying column
.te_reading <- function(x, u) {
  rn <- .te_ratio_name(x@kind)
  out <- character(0)
  for (m in u$final$map) {
    if (!length(m$tv)) next
    col <- names(u$final$coefficients)[m$col]
    e <- u$effects[u$effects$column == col, , drop = FALSE]
    cv <- u$curve[u$curve$column == col, , drop = FALSE]
    first <- e[1L, ]; last <- e[nrow(e), ]
    clear <- function(r) r$lower > 1 || r$upper < 1          # interval excludes 1
    shape <- if (sign(first$log_hr) != sign(last$log_hr) && clear(first) && clear(last)) {
      inr <- cv[cv$time >= first$time & cv$time <= last$time, , drop = FALSE]
      cross <- which(diff(sign(inr$log_hr)) != 0)
      paste0("changes direction (from ", if (first$log_hr > 0) "above" else "below", " 1 to ", if (last$log_hr > 0) "above" else "below",
             " 1)", if (length(cross)) paste0(" near t = ", format(signif(inr$time[cross[1L] + 1L], 2L))) else "")
    } else if (abs(last$log_hr - first$log_hr) < 0.05) "hardly changes between these times"
      else if (abs(last$log_hr) < abs(first$log_hr)) "moves towards 1 (the effect weakens)"
      else "moves away from 1 (the effect strengthens)"
    out <- c(out, paste0("- ", .te_pretty(x, col), " (", m$basis$label, "): ", rn, " ", .num(first$hr, 2L), " (", .ci_text(first$lower, first$upper),
                         ") at t = ", format(first$time), ", ", .num(last$hr, 2L), " (", .ci_text(last$lower, last$upper), ") at t = ",
                         format(last$time), "; the ", rn, " ", shape, "."))
  }
  out
}

.te_lines <- function(x, detail = FALSE) {
  rn <- .te_ratio_name(x@kind)
  lv <- .pct(x@level)
  head <- switch(x@kind, cox = "Cox model", csc = "cause-specific Cox model", fg = "Fine-Gray model")
  s <- x@settings
  out <- c(paste0("Time-varying effects in a ", head),
           paste0("  Screening: Grambsch-Therneau score test per term (time scale ", s$transform, "), ",
                  switch(s$screen, holm = paste0("Holm-adjusted p < ", s$alpha), unadjusted = paste0("p < ", s$alpha),
                         all = "every term considered")),
           paste0("  Forms compared: ", paste(x@forms, collapse = ", "), "; chosen by ", toupper(s$criterion)), "")
  for (u in x@units) {
    if (length(x@units) > 1L) out <- c(out, paste0("== ", .te_title(u$label), " =="))
    sc <- u$screening
    tab <- data.frame(term = .te_pretty(x, sc$term), `chi-square` = .num(sc$chisq, 2L), df = sc$df, p = .pval(sc$p),
                      `p (Holm)` = .pval(sc$p_holm), considered = ifelse(sc$considered, "yes", "no"), check.names = FALSE)
    out <- c(out, paste0("Test of a constant effect (", u$nevent, " events)"), .table_lines(tab))
    if (!is.null(u$global)) out <- c(out, paste0("  global: chi-square ", .num(u$global[["chisq"]], 2L), " on ", u$global[["df"]], " df, p ",
                                                 .pval(u$global[["p"]])))
    cmp <- u$comparison[u$comparison$term %in% names(u$selected)[u$selected != "constant"] |
                          (detail & u$comparison$term %in% u$screening$term[u$screening$considered]), , drop = FALSE]
    if (nrow(cmp)) {
      ct <- data.frame(term = .te_pretty(x, cmp$term), form = cmp$label, `log-lik` = .num(cmp$loglik, 2L),
                       par = cmp$k, LR = .num(cmp$lr, 2L), p = .pval(cmp$p_lr), d = .num(cmp$delta, 1L),
                       chosen = ifelse(cmp$chosen, "<-", ""), check.names = FALSE)
      names(ct)[7L] <- paste0("d", toupper(s$criterion))
      out <- c(out, "", "Comparison of forms (LR: against a constant effect)", .table_lines(ct))
    }
    chosen <- u$selected[u$selected != "constant"]
    out <- c(out, "", if (length(chosen)) paste0("Chosen: ", paste(paste0(.te_pretty(x, names(chosen)), " = ", chosen), collapse = "; "))
                      else "Chosen: every effect constant (proportional hazards kept).")
    if (u$any_tv) {
      out <- c(out, paste0("Final model: likelihood ratio test of all time-varying terms ", .num(u$final$lr, 2L), " on ", u$final$lr_df,
                           " df, p ", .pval(u$final$p), if (u$final$robust) "; robust standard errors" else ""))
      e <- u$effects[u$effects$form != "constant", , drop = FALSE]
      if (!is.null(x@bootstrap) && !is.null(x@bootstrap[[match(u$label, names(x@units))]]$ci)) {
        bc <- x@bootstrap[[match(u$label, names(x@units))]]$ci
        i <- match(paste(e$column, e$time), paste(bc$column, bc$time))
        e$blo <- bc$lower[i]; e$bhi <- bc$upper[i]
      }
      et <- data.frame(term = .te_pretty(x, e$column), time = format(e$time), ratio = .num(e$hr, 2L),
                       ci = .ci_text(e$lower, e$upper), p = .pval(e$p), check.names = FALSE)
      names(et)[3:4] <- c(rn, paste0(lv, " CI"))
      if (!is.null(e$blo)) et[[paste0("bootstrap ", lv, " CI")]] <- .ci_text(e$blo, e$bhi)
      out <- c(out, "", paste0(rn, " over time"), .table_lines(et), "", "Reading:", .te_reading(x, u))
    }
    const <- u$effects[u$effects$form == "constant" & !duplicated(u$effects$column), , drop = FALSE]
    if (detail && nrow(const)) {
      ct <- data.frame(term = .te_pretty(x, const$column), ratio = .num(const$hr, 2L), ci = .ci_text(const$lower, const$upper),
                       p = .pval(const$p), check.names = FALSE)
      names(ct)[2:3] <- c(rn, paste0(lv, " CI"))
      out <- c(out, "", "Constant effects in the final model", .table_lines(ct))
    }
    if (!is.null(x@bootstrap)) {
      bt <- x@bootstrap[[match(u$label, names(x@units))]]
      if (!is.null(bt$freq)) {
        fq <- bt$freq[bt$freq$term %in% u$screening$term, , drop = FALSE]
        fq <- fq[order(match(fq$term, u$screening$term), -fq$share), , drop = FALSE]
        ft <- data.frame(term = .te_pretty(x, fq$term), form = fq$form, `chosen in` = paste0(round(100 * fq$share), "%"), check.names = FALSE)
        out <- c(out, "", paste0("Stability of the choice (", bt$ok, " of ", bt$B, " bootstrap samples)"), .table_lines(ft))
        top <- tapply(fq$share, fq$term, max)
        tvs <- names(u$selected)[u$selected != "constant"]
        unst <- intersect(names(top)[top < 0.5], tvs)
        cst <- vapply(tvs, function(tm) sum(fq$share[fq$term == tm & fq$form != "constant"]), numeric(1))
        if (length(tvs)) out <- c(out, paste0("  A time-varying effect was found in ", paste0(round(100 * cst), "%", collapse = ", "),
                                              " of the samples for ", paste(.te_pretty(x, tvs), collapse = ", "), "."))
        if (length(unst)) out <- c(out, paste0("  The form is not well determined for ", paste(.te_pretty(x, unst), collapse = ", "),
                                               " (no form chosen in half of the samples): report the hazard ratios at a few times, ",
                                               "with the bootstrap intervals, rather than the form."))
      } else out <- c(out, "", paste0("Bootstrap: only ", bt$ok, " of ", bt$B, " samples could be analysed; no stability results"))
    }
    if (length(u$notes)) out <- c(out, "", paste0("Note: ", u$notes))
    out <- c(out, "")
  }
  out <- c(out, "Cautions:",
           "- The forms were chosen from the same data; the intervals of the final model do not include that choice.",
           if (is.null(x@bootstrap)) "  Use bootstrap = 200 for intervals and a stability check that include it." else
             "  The bootstrap intervals include the screening and the choice of the forms.",
           "- A non-significant test is not proof of proportional hazards; with few events the test has little power.",
           "- Judge the size of the change on plot(x), not only the p-values.")
  c(out, "", .method_lines(.method_keys(x)))
}

.te_title <- function(s) paste0(toupper(substr(s, 1L, 1L)), substr(s, 2L, nchar(s)))

S7::method(print, TimeEffects) <- function(x, ...) { cat(.te_lines(x, detail = FALSE), sep = "\n"); invisible(x) }
S7::method(summary, TimeEffects) <- function(object, ...) { cat(.te_lines(object, detail = TRUE), sep = "\n"); invisible(object) }

S7::method(as.data.frame, TimeEffects) <- function(x, row.names = NULL, optional = FALSE, what = c("effects", "screening", "forms", "coefficients", "stability"), ...) {
  what <- match.arg(what)
  pick <- function(u) switch(what,
    effects = u$effects,
    screening = u$screening,
    forms = u$comparison[, c("term", "form", "label", "loglik", "k", "ic", "delta", "lr", "lr_df", "p_lr", "converged", "chosen")],
    coefficients = {
      b <- u$final$coefficients; se <- sqrt(diag(u$final$vcov))
      data.frame(term = names(b), estimate = unname(b), se = unname(se), z = unname(b / se),
                 p = 2 * stats::pnorm(-abs(unname(b / se))), stringsAsFactors = FALSE)
    },
    stability = {
      if (is.null(x@bootstrap)) stop("no bootstrap: use time_effects(..., bootstrap = 200)", call. = FALSE)
      x@bootstrap[[match(u$label, names(x@units))]]$freq
    })
  out <- do.call(rbind, lapply(x@units, function(u) {
    d <- pick(u)
    if (is.null(d) || !nrow(d)) return(NULL)
    if (length(x@units) > 1L) d <- cbind(cause = u$cause, d)
    d
  }))
  rownames(out) <- NULL
  out
}

S7::method(plot, TimeEffects) <- function(x, term = NULL, cause = NULL, residuals = TRUE, scale = c("ratio", "log"),
                                          col = "#0072B2", ...) {
  scale <- match.arg(scale)
  us <- x@units
  if (!is.null(cause)) {
    us <- us[vapply(us, function(u) identical(as.integer(u$cause), as.integer(cause)), logical(1))]
    if (!length(us)) stop("cause ", cause, " is not in the result", call. = FALSE)
  }
  panels <- list()
  for (u in us) {
    cols <- vapply(u$final$map, function(m) names(u$final$coefficients)[m$col], character(1))
    tv <- vapply(u$final$map, function(m) length(m$tv) > 0L, logical(1))
    want <- if (is.null(term)) cols[tv] else cols[cols %in% term | u$term_of %in% term]
    if (!is.null(term) && !length(want)) stop("term not in the model: ", paste(term, collapse = ", "), call. = FALSE)
    for (cn in want) panels[[length(panels) + 1L]] <- list(u = u, col = cn)
  }
  if (!length(panels)) {
    message("every effect is constant; plotting the scaled Schoenfeld residuals of the terms instead")
    for (u in us) for (cn in names(u$constant$coefficients)) panels[[length(panels) + 1L]] <- list(u = u, col = cn)
  }
  np <- length(panels)
  op <- graphics::par(mfrow = if (np == 1L) c(1L, 1L) else c(ceiling(np / min(np, 2L)), min(np, 2L)), mar = c(4, 4.2, 2.6, 1))
  on.exit(graphics::par(op))
  rn <- .te_ratio_name(x@kind)
  for (pn in panels) {
    u <- pn$u; cn <- pn$col
    cv <- u$curve[u$curve$column == cn, , drop = FALSE]
    f <- if (scale == "ratio") identity else log
    tr <- if (scale == "ratio") exp else identity
    y <- f(cv$hr); lo <- f(cv$lower); hi <- f(cv$upper)
    rs <- u$resid
    pts <- NULL
    if (residuals && !is.null(rs) && cn %in% colnames(rs$y)) {
      ry <- rs$y[, cn]; rt <- rs$time
      keep <- is.finite(ry) & rt >= min(cv$time) & rt <= max(cv$time)
      pts <- list(t = rt[keep], y = ry[keep])
      if (length(pts$t) > 1000L) { s <- seq(1L, length(pts$t), length.out = 1000L); pts <- list(t = pts$t[s], y = pts$y[s]) }
    }
    # the range of the curve and its band, widened by half on the log scale; residuals outside it are clipped
    lr <- range(log(c(cv$lower, cv$upper, exp(u$constant$coefficients[[cn]]))), na.rm = TRUE)
    lr <- lr + c(-0.5, 0.5) * max(diff(lr), 0.1)
    yl <- if (scale == "ratio") exp(lr) else lr
    graphics::plot(cv$time, y, type = "n", ylim = yl, log = if (scale == "ratio") "y" else "", xlab = "Time",
                   ylab = if (scale == "ratio") rn else paste0("log ", rn), ...)
    if (!is.null(pts)) graphics::points(pts$t, if (scale == "ratio") exp(pts$y) else pts$y, pch = 16, cex = 0.4, col = "grey75")
    graphics::polygon(c(cv$time, rev(cv$time)), c(lo, rev(hi)), col = grDevices::adjustcolor(col, 0.2), border = NA)
    graphics::lines(cv$time, y, lwd = 2, col = col)
    graphics::abline(h = tr(0), lty = 3, col = "grey40")
    b0 <- u$constant$coefficients[[cn]]
    graphics::abline(h = tr(b0), lty = 2, col = "grey30")
    m <- u$final$map[[match(cn, names(u$final$coefficients))]]
    ttl <- paste0(.te_pretty(x, cn), ": ", if (length(m$tv)) m$basis$label else "constant", if (length(x@units) > 1L) paste0(" (", u$label, ")") else "")
    graphics::title(ttl, cex.main = 0.9)
    graphics::rug(u$event_times, col = "grey60", ticksize = 0.02)
  }
  graphics::legend("topright", c("chosen form", "constant (PH) fit", if (residuals) "scaled Schoenfeld residuals"),
                   lty = c(1, 2, NA), pch = c(NA, NA, if (residuals) 16), col = c(col, "grey30", "grey75"), bty = "n", cex = 0.75)
  invisible(x)
}

# cumulative function along each row of a matrix, keeping the matrix shape
.te_rowcum <- function(M, f) {
  if (!ncol(M)) return(M)
  out <- t(apply(M, 1L, f))
  matrix(out, nrow = nrow(M), ncol = ncol(M))
}

# predictions with the time-varying effects: survival (cox), cumulative incidence (fg, csc)
S7::method(predict, TimeEffects) <- function(object, newdata = NULL, times, cause = 1L, ...) {
  if (missing(times)) stop("times are needed", call. = FALSE)
  fit <- object@fit
  Xn <- .sx_newX(fit, newdata)
  sn <- if (object@kind == "fg") NULL else .sx_newstrata(fit, newdata)
  haz <- function(u) {                               # hazard increments: rows of newdata x event times
    sp <- u$split
    b <- u$final$coefficients
    bh <- .cox_baseline_raw(b, sp$stop, sp$status, sp$X, sp$w, sp$ties, sp$start, sp$strata)
    tj <- sort(unique(bh$time))
    nb <- ncol(fit@X)
    lp0 <- drop(Xn %*% b[seq_len(nb)])
    LP <- matrix(lp0, nrow(Xn), length(tj))
    for (m in u$final$map) if (length(m$tv)) LP <- LP + outer(Xn[, m$col], rep(1, length(tj))) * matrix(drop(m$basis$B(tj) %*% b[m$tv]), nrow(Xn), length(tj), byrow = TRUE)
    st <- if (is.null(sn)) rep(0L, nrow(Xn)) else sn
    dH <- matrix(0, nrow(Xn), length(tj))
    for (s in unique(st)) {
      bs <- bh[bh$stratum == s, , drop = FALSE]
      if (!nrow(bs)) next
      inc <- diff(c(0, bs$cumhaz))
      j <- match(bs$time, tj)
      i <- which(st == s)
      dH[i, j] <- outer(rep(1, length(i)), inc) * exp(LP[i, j, drop = FALSE])
    }
    list(t = tj, dH = dH)
  }
  at <- function(tj, M) {                            # value of the cumulative matrix M (rows x tj) at times
    k <- findInterval(times, tj)
    out <- matrix(0, nrow(M), length(times))
    for (q in seq_along(times)) if (k[q] > 0L) out[, q] <- M[, k[q]]
    out
  }
  out <- switch(object@kind,
    cox = { h <- haz(object@units[[1L]]); exp(-at(h$t, .te_rowcum(h$dH, cumsum))) },
    fg = { h <- haz(object@units[[1L]]); 1 - exp(-at(h$t, .te_rowcum(h$dH, cumsum))) },
    csc = {
      hs <- lapply(object@units, haz)
      tj <- sort(unique(unlist(lapply(hs, `[[`, "t"))))
      dk <- lapply(hs, function(h) { M <- matrix(0, nrow(Xn), length(tj)); M[, match(h$t, tj)] <- h$dH; M })
      tot <- Reduce(`+`, dk)
      S <- if (isTRUE(fit@product_limit)) .te_rowcum(pmax(1 - tot, 0), cumprod) else exp(-.te_rowcum(tot, cumsum))
      Sm <- cbind(1, S[, -ncol(S), drop = FALSE])
      k <- match(as.integer(cause), vapply(object@units, `[[`, integer(1), "cause"))
      if (is.na(k)) stop("cause ", cause, " is not in the model", call. = FALSE)
      at(tj, .te_rowcum(Sm * dk[[k]], cumsum))
    })
  dimnames(out) <- list(NULL, paste0("t=", times))
  out
}

# ----------------------------------------------------------------------------------------------
# report(): one self-contained HTML file
# ----------------------------------------------------------------------------------------------
S7::method(report, TimeEffects) <- function(x, file = "time_effects_report.html", width = 7, height = 5, ...) {
  esc <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", s)))
  img <- function(u, cn) {
    tf <- tempfile(fileext = ".png")
    grDevices::png(tf, width = width, height = height, units = "in", res = 110)
    ok <- tryCatch({ plot(x, term = cn, cause = if (length(x@units) > 1L) u$cause else NULL); TRUE }, error = function(e) FALSE)
    grDevices::dev.off()
    if (!ok || !file.exists(tf)) return("")
    on.exit(unlink(tf))
    paste0("<img alt='", esc(cn), "' src='data:image/png;base64,", .b64(readBin(tf, "raw", file.info(tf)$size)), "'/>")
  }
  rn <- .te_ratio_name(x@kind)
  s <- x@settings
  figs <- unlist(lapply(x@units, function(u) {
    cols <- vapply(u$final$map, function(m) if (length(m$tv)) names(u$final$coefficients)[m$col] else NA_character_, character(1))
    vapply(cols[!is.na(cols)], function(cn) img(u, cn), character(1))
  }))
  tables <- unlist(lapply(x@units, function(u) {
    sc <- u$screening
    cmp <- u$comparison
    e <- u$effects
    c(if (length(x@units) > 1L) paste0("<h3>", esc(.te_title(u$label)), "</h3>"),
      "<h4>Test of a constant effect</h4>",
      .html_table(data.frame(term = .te_pretty(x, sc$term), `chi-square` = .num(sc$chisq, 2L), df = sc$df, p = .pval(sc$p),
                             `p (Holm)` = .pval(sc$p_holm), considered = ifelse(sc$considered, "yes", "no"), check.names = FALSE)),
      "<h4>Comparison of forms</h4>",
      .html_table(data.frame(term = .te_pretty(x, cmp$term), form = cmp$label, `log-likelihood` = .num(cmp$loglik, 2L),
                             parameters = cmp$k, `LR vs constant` = .num(cmp$lr, 2L), p = .pval(cmp$p_lr),
                             delta = .num(cmp$delta, 1L), chosen = ifelse(cmp$chosen, "yes", ""), check.names = FALSE)),
      paste0("<h4>", rn, " at chosen times (final model)</h4>"),
      .html_table(data.frame(term = .te_pretty(x, e$column), form = e$form, time = format(e$time), ratio = .num(e$hr, 3L),
                             lower = .num(e$lower, 3L), upper = .num(e$upper, 3L), p = .pval(e$p), check.names = FALSE)))
  }))
  methods <- paste0(
    "<p>For every term the proportional hazards assumption was tested with the score test of Grambsch and Therneau (1994) ",
    "(time scale: ", s$transform, "), ",
    switch(s$screen, holm = paste0("and terms with a Holm (1979) adjusted p-value below ", s$alpha, " were considered for a time-varying effect. "),
           unadjusted = paste0("and terms with a p-value below ", s$alpha, " were considered for a time-varying effect. "),
           all = "and every term was considered for a time-varying effect. "),
    "For each such term the log hazard ratio was modelled as beta + gamma'B(t) with the forms ", paste(x@forms, collapse = ", "),
    " (log-linear: B(t) = log(t / median event time); step: an indicator of t after a change point chosen by profile likelihood among ",
    "the 20th to 80th percentiles of the event times, counted as one extra parameter; spline: a natural cubic spline in log t ",
    "with ", s$df, " degrees of freedom, as in Hess (1994)), each fitted by partial likelihood on follow-up split at every event ",
    "time (Therneau and Grambsch 2000, Chapter 6), and the form with the smallest ", toupper(s$criterion), " was chosen ",
    if (s$criterion == "aic") "(Akaike 1974)" else "(Schwarz 1978)", ". The final model contains all chosen forms together. ",
    if (x@kind == "fg") "For the Fine-Gray model the same was done on the weighted counting-process data of Geskus (2011), with robust standard errors (Lin and Wei 1989). " else "",
    if (!is.null(x@bootstrap)) paste0("The whole procedure was repeated in ", x@bootstrap[[1L]]$B, " bootstrap samples of the subjects (Efron and Tibshirani 1993); ",
                                     "the percentile intervals and the frequencies with which each form was chosen include the screening and the choice. ") else "",
    "</p>")
  refs <- .scan_refs(.method_keys(x))
  code <- paste(deparse(x@model_call), collapse = "\n")
  html <- c("<!DOCTYPE html><html><head><meta charset='utf-8'><title>SurvX time-varying effects</title>",
            "<style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;color:#222;line-height:1.5}",
            "table{border-collapse:collapse;font-size:13px;margin:1em 0}th,td{border-bottom:1px solid #ddd;padding:4px 8px;text-align:right}",
            "th:first-child,td:first-child{text-align:left}pre{background:#f6f6f6;padding:1em;overflow-x:auto}",
            "img{max-width:100%}.note{background:#fff8e1;padding:.6em 1em;border-left:4px solid #f0b400}</style></head><body>",
            "<h1>Time-varying effects</h1>", paste0("<p>SurvX ", as.character(utils::packageVersion("SurvX")), ", ", format(Sys.time(), "%Y-%m-%d %H:%M"), "</p>"),
            "<h2>Summary</h2>", paste0("<pre>", esc(paste(.te_lines(x, detail = TRUE), collapse = "\n")), "</pre>"),
            "<h2>Tables</h2>", tables,
            if (length(figs)) c("<h2>Figures</h2>", "<p>Solid line and band: the chosen form with its pointwise confidence band; dashed: the constant (proportional hazards) fit; grey points: scaled Schoenfeld residuals.</p>", figs),
            "<h2>Methods</h2>", methods,
            paste0("<div class='note'><p><b>Cautions.</b> The forms were chosen from the same data, so the model-based intervals do not include that choice",
                   if (is.null(x@bootstrap)) "; rerun with <code>bootstrap = 200</code> for intervals and a stability check that do." else
                     "; the bootstrap intervals and selection frequencies do.",
                   " A non-significant test is not evidence of proportional hazards when there are few events.</p></div>"),
            "<h2>Code</h2>", paste0("<pre>", esc(code), "</pre>"),
            "<h2>References</h2>", paste0("<p>", esc(refs), "</p>"),
            "</body></html>")
  writeLines(html, file, useBytes = TRUE)
  message("report written to ", normalizePath(file, mustWork = FALSE))
  invisible(file)
}
