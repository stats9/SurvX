.base_codes <- c(gompertz = 0L, gompertz3 = 1L, weibull = 2L, mix_weibull = 3L, mix_loglogistic = 4L,
                 mix_lognormal = 5L, logistic3 = 6L, spline = 7L)

.cond_codes <- c(logistic = 0L, weibull = 1L, loglogistic = 2L, lognormal = 3L)

.cond_names <- function(cond) {
  switch(cond, logistic = c("log_b", "c"), lognormal = c("mu", "log_sigma"), c("log_lambda", "log_kappa"))
}

.parse_link <- function(lk) {
  if (is.numeric(lk)) return(list(free = 0L, value = as.double(lk), name = paste0("alpha=", lk)))
  if (identical(lk, "PH")) return(list(free = 0L, value = 0, name = "PH"))
  if (identical(lk, "PO")) return(list(free = 0L, value = 1, name = "PO"))
  if (identical(lk, "GOR")) return(list(free = 1L, value = 0.5, name = "GOR"))
  stop("link must be \"PH\", \"PO\", \"GOR\" or a number", call. = FALSE)
}

.base_names <- function(baseline, nk) {
  switch(baseline,
    gompertz = c("rho", "log_tau"),
    gompertz3 = c("rho", "log_tau", "zeta"),
    weibull = c("log_lambda", "log_kappa"),
    mix_weibull = , mix_loglogistic = c("logit_p", "log_lambda", "log_kappa"),
    mix_lognormal = c("logit_p", "mu", "log_sigma"),
    logistic3 = c("log_b", "c", "logit_p"),
    spline = paste0("gamma", seq_len(max(nk, 2L)) - 1L))
}

#' Direct parametric regression for cumulative incidence functions
#'
#' @description
#' Models the cumulative incidence of every cause directly with a parametric form, following Jeong
#' and Fine (2006, 2007):
#' \deqn{F_k(t \mid z) = 1 - \{1 + \alpha_k e^{\beta_k^\top z} u_k(t)\}^{-1/\alpha_k},}
#' a generalized odds-rate model in which \eqn{\alpha_k \to 0} gives proportional subdistribution
#' hazards (the Fine-Gray model with a parametric baseline) and \eqn{\alpha_k = 1} proportional odds
#' of the cumulative incidence. Three estimators are available: the full likelihood of all causes,
#' a weighted estimator for one cause, and a constrained model for two causes whose incidences add up
#' to one.
#'
#' \strong{Use it when} smooth, parametric cumulative incidences are wanted (extrapolation, the
#' long-term incidence \eqn{F_k(\infty)}, health-economic models), when the effect is better described
#' on the odds scale, or when all causes are to be modelled together with coherent predictions.
#' \strong{Use something else when} no assumption on the shape of the incidence is wanted:
#' \code{\link{fg}}; or the question is about rates: \code{\link{csc}}.
#' @param formula \code{Event(time, status) ~ covariates} (or \code{~ 1}); status 0 = censored,
#'   1, 2, ..., K = cause. Delayed entry is not supported.
#' @param data a data frame. Rows with missing values are removed.
#' @param baseline the baseline family, one for all causes or one per cause (see Details):
#'   \code{"gompertz"} (default), \code{"gompertz3"}, \code{"weibull"}, \code{"mix_weibull"},
#'   \code{"mix_loglogistic"}, \code{"mix_lognormal"}, \code{"logistic3"} or \code{"spline"}.
#'   Families with a plateau (Gompertz with \eqn{\rho < 0}, the mixtures, \code{"logistic3"}) estimate
#'   the long-term incidence \eqn{F_k(\infty)}; \code{"weibull"} and \code{"spline"} go to 1.
#' @param link \code{"PH"} (\eqn{\alpha = 0}, default), \code{"PO"} (\eqn{\alpha = 1}), \code{"GOR"}
#'   (\eqn{\alpha} estimated), or a number giving a fixed \eqn{\alpha}; one for all causes or one per
#'   cause.
#' @param estimation \code{"joint"} (default): all causes by the full likelihood. \code{"weighted"}:
#'   the cause of interest alone, with censoring weights as in Fine and Gray, so that no model for
#'   the other causes is needed. \code{"constrained"}: two causes whose incidences add up to one at
#'   infinity (Shi, Cheng and Jeong 2013). See Details.
#' @param cause the cause of interest for \code{"weighted"} and \code{"constrained"}.
#' @param cond \code{"constrained"} only: the distribution \eqn{G} of the time to the competing
#'   event given that it occurs: \code{"logistic"} (modified logistic, as in the paper),
#'   \code{"weibull"}, \code{"loglogistic"} or \code{"lognormal"}.
#' @param df degrees of freedom of the spline baseline.
#' @param knots optional list with the knots (log time, boundary knots first and last) of every
#'   cause, for the spline baseline.
#' @param constraint \code{"joint"} only: keep the sum of the fitted incidences below 1 within
#'   follow-up (see Details).
#' @param censoring_grid \code{"weighted"} only: coarsen the estimated censoring distribution to at
#'   most this many support points (0 = exact). With the PH link the exact fit is linear in the
#'   sample size; with other links its cost grows with the number of competing events times the
#'   number of support points (a few seconds at \eqn{n = 10000}), and a grid can shorten it.
#' @param init optional starting values on the estimation scale.
#' @param max_iter, tol iteration limit and convergence tolerance.
#' @details
#' \subsection{Estimators}{
#' \describe{
#'   \item{\code{"joint"}}{all causes jointly by the full likelihood of Jeong and Fine (2006, 2007):
#'     events of cause \eqn{k} contribute \eqn{f_k(t)}, censored subjects
#'     \eqn{1 - \sum_k F_k(t)}. Efficient when every cause is modelled correctly; appropriate when a
#'     share of the subjects never fails (the incidences level off below a total of one).}
#'   \item{\code{"weighted"}}{the cause of interest alone. Subjects who fail from a competing cause
#'     enter through the Kaplan-Meier estimate of the censoring distribution, as in Fine and Gray
#'     (1999), so no model for the other causes is needed (Geskus 2011; Lambert, Wilkes and
#'     Crowther 2017). Standard errors are from the sandwich estimator with one score per subject;
#'     the estimation of the censoring distribution is ignored. \code{logLik} is a weighted
#'     pseudo-likelihood: use it only to compare links or baselines of weighted fits.}
#'   \item{\code{"constrained"}}{two causes with asymptotes that add up to one (Shi, Cheng and Jeong
#'     2013): the cause of interest as above and \eqn{F_2(t \mid z) = \{1 - F_1(\infty \mid z)\} G(t)}.
#'     Needs a baseline with a plateau for the cause of interest (not \code{"weibull"}, not
#'     \code{"spline"}). Appropriate when every subject eventually fails from one of the two causes.}
#' }
#' }
#' \subsection{Baseline families}{
#' Baseline families (\code{baseline}) and their parameters on the scale used for estimation:
#' \describe{
#'   \item{\code{"gompertz"}}{\eqn{u(t) = \tau (e^{\rho t} - 1)/\rho}; \code{rho, log_tau}. Jeong and Fine (2006).}
#'   \item{\code{"gompertz3"}}{three-parameter Gompertz; \code{rho, log_tau, zeta}. Haile et al. (2016).}
#'   \item{\code{"weibull"}}{\eqn{u(t) = (\lambda t)^\kappa}; \code{log_lambda, log_kappa}.}
#'   \item{\code{"mix_weibull"}, \code{"mix_loglogistic"}, \code{"mix_lognormal"}}{baseline CIF
#'     \eqn{F_0(t) = p\,G(t)} with \eqn{G} a Weibull, log-logistic or log-normal distribution function;
#'     \code{logit_p} and two parameters of \eqn{G}. Cheng (2009), Larson and Dinse (1985).}
#'   \item{\code{"logistic3"}}{modified three-parameter logistic baseline CIF; \code{log_b, c, logit_p}.
#'     Cheng (2009), Shi, Cheng and Jeong (2013).}
#'   \item{\code{"spline"}}{\eqn{\log u(t)} a restricted cubic spline in \eqn{\log t} with \code{df}
#'     degrees of freedom; \code{gamma0, gamma1, ...}. Royston and Parmar (2002), Mozumder et al. (2018).}
#' }
#' For baselines given as a CIF the baseline is the CIF at \eqn{z = 0}; centre continuous covariates.
#' }
#' \subsection{Coherence and convergence}{
#' The published models do not force \eqn{\sum_k F_k \le 1}. By default the model is fitted as
#' published and \code{fit@diagnostics} reports how many subjects have a fitted sum above 1 within
#' follow-up and at infinity. \code{constraint = TRUE} keeps the sum below 1 within follow-up.
#'
#' A fit is reported as converged only at a stationary point with a positive definite information
#' matrix. Families with a plateau are not identified when the data show no plateau; such fits
#' are flagged and their standard errors should not be used.
#' }
#' @section Interpreting the coefficients:
#' \describe{
#'   \item{PH link}{\eqn{e^\beta} is a \strong{subdistribution hazard ratio}, read as in
#'     \code{\link{fg}}: its direction is that of the effect on the cumulative incidence.}
#'   \item{PO link}{\eqn{e^\beta} is an \strong{odds ratio of the cumulative incidence}: the odds
#'     \eqn{F_k(t)/\{1 - F_k(t)\}} are multiplied by \eqn{e^\beta} at every time. OR 1.8 means 80\%
#'     higher odds of having had the event by any time.}
#'   \item{GOR link}{\eqn{\alpha} between the two; the effect is on the scale
#'     \eqn{\log[\{(1 - F)^{-\alpha} - 1\}/\alpha]}. Compare links with \code{AIC()} and report the
#'     effect on the incidence itself.}
#'   \item{Baseline parameters}{describe the incidence at \eqn{z = 0}; centre continuous covariates
#'     so that this is a typical subject.}
#' }
#' The long-term incidence \code{predict(fit, newdata, times = Inf)} is often the most useful
#' summary with plateau families. A sentence for a report: \emph{"In a direct Gompertz model for
#' the cumulative incidence of relapse (proportional odds), treatment reduced the odds of relapse by
#' 35\% (OR 0.65, 95\% CI 0.48 to 0.88); the predicted long-term incidence was 28\% versus 37\%."}
#' @section Assumptions and how to check them:
#' \describe{
#'   \item{The baseline family fits}{\code{plot(fit)} and \code{cif_plot(fit, type = "unadjusted")}
#'     compare the fitted incidences with the non-parametric \code{\link{aj}} estimate; compare
#'     families and links with \code{AIC()} (joint fits; for weighted fits, the pseudo-likelihood).}
#'   \item{A plateau for plateau families}{without one in the data, \eqn{F_k(\infty)} is not
#'     identified; such fits are flagged as not converged.}
#'   \item{Joint estimation models every cause}{a wrong model for one cause can bias the other;
#'     the weighted estimator avoids this for the cause of interest.}
#' }
#' @section Validation:
#' There is no other implementation of these models in R to compare with. The log-likelihoods,
#' scores and information matrices were checked against independent code and numerical
#' differentiation, and the fast evaluation of the weighted estimator against the direct
#' evaluation of every (subject, censoring time) pair to \eqn{10^{-9}}. In a Monte Carlo study of the
#' validation suite (\url{https://github.com/stats9/SurvX}, directory \code{validation}; 1000
#' replications, \eqn{n} = 200 to 1000, PH and PO links, all three estimators) the standard errors
#' matched the empirical standard deviation (ratio 0.93 to 1.04) and the 95\% intervals covered 93.2
#' to 95.8\%. On \eqn{n = 10000}, a joint fit takes about 0.6 s and a weighted fit 0.1 s.
#' @return
#' An S7 object of class \code{DirectCIFFit}. Coefficients are named \code{"<cause>:<term>"};
#' \code{fit@diagnostics} reports how many subjects have fitted incidences adding up to more than 1.
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(fit)}, \code{summary(fit)} \tab effects per cause (sHR or OR) with intervals,
#'     baseline parameters, \code{Reading:} lines \cr
#'   \code{\link{estimates}(fit)}, \code{coef}, \code{vcov}, \code{confint}, \code{logLik},
#'     \code{AIC}, \code{BIC}, \code{nobs} \tab extractors \cr
#'   \code{predict(fit, newdata, times, cause)} \tab cumulative incidence; \code{times = Inf} gives
#'     the long-term incidence \cr
#'   \code{plot(fit)}, \code{\link{cif_plot}}, \code{\link{cif_data}} \tab cumulative incidence
#'     curves \cr
#' }
#' @references
#' Jeong J-H, Fine JP (2006) Direct parametric inference for the cumulative incidence function.
#' J R Stat Soc C 55:187-200.
#'
#' Jeong J-H, Fine JP (2007) Parametric regression on cumulative incidence function.
#' Biostatistics 8:184-196.
#'
#' Cheng Y (2009) Modeling cumulative incidences of dementia and dementia-free death using a novel
#' three-parameter logistic function. Int J Biostat 5(1), Article 29.
#'
#' Shi H, Cheng Y, Jeong J-H (2013) Constrained parametric model for simultaneous inference of two
#' cumulative incidence functions. Biom J 55:82-96.
#'
#' Haile SR, Jeong J-H, Chen X, Cheng Y (2016) A 3-parameter Gompertz distribution for survival data
#' with competing risks, with an application to breast cancer data. J Appl Stat 43:2239-2253.
#'
#' Royston P, Parmar MKB (2002) Flexible parametric proportional-hazards and proportional-odds
#' models for censored survival data. Stat Med 21:2175-2197.
#'
#' Fine JP, Gray RJ (1999) A proportional hazards model for the subdistribution of a competing risk.
#' J Am Stat Assoc 94:496-509.
#'
#' Geskus RB (2011) Cause-specific cumulative incidence estimation and the Fine and Gray model under
#' both left truncation and right censoring. Biometrics 67:39-49.
#'
#' Lambert PC, Wilkes SR, Crowther MJ (2017) Flexible parametric modelling of the cause-specific
#' cumulative incidence function. Stat Med 36:1429-1446.
#'
#' Mozumder SI, Rutherford MJ, Lambert PC (2018) Direct likelihood inference on the cause-specific
#' cumulative incidence function: a flexible parametric regression modelling approach.
#' Stat Med 37:82-97.
#'
#' All references: \code{sx_references("directcif")}.
#' @seealso
#' \code{\link{fg}} (semi-parametric Fine-Gray), \code{\link{csc}}, \code{\link{aj}},
#' \code{\link{cif_plot}}, \code{\link{predict.SurvX}}.
#' @examples
#' set.seed(2)
#' n <- 400
#' d <- data.frame(x = rbinom(n, 1, 0.5))
#' t1 <- rexp(n, 0.4 * exp(0.6 * d$x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
#' d$time <- pmin(t1, t2, cens)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#'
#' ## 1. joint Gompertz model for both causes, proportional subdistribution hazards
#' fit <- directcif(Event(time, status) ~ x, d, baseline = "gompertz", link = "PH")
#' fit
#' predict(fit, data.frame(x = 0:1), times = c(1, 2, Inf))    # Inf: long-term incidence
#'
#' ## 2. proportional odds, and comparison of the links
#' fpo <- directcif(Event(time, status) ~ x, d, link = "PO")
#' AIC(fit, fpo)
#'
#' ## 3. cause 1 alone, without a model for cause 2
#' directcif(Event(time, status) ~ x, d, estimation = "weighted", cause = 1)
#'
#' ## 4. two causes whose incidences add up to one at infinity
#' fc <- directcif(Event(time, status) ~ x, d, estimation = "constrained", cause = 1, cond = "weibull")
#' nd <- data.frame(x = 0:1)
#' predict(fc, nd, times = Inf, cause = 1) + predict(fc, nd, times = Inf, cause = 2)
#' @export
directcif <- function(formula, data, baseline = "gompertz", link = "PH",
                      estimation = c("joint", "weighted", "constrained"), cause = 1L,
                      cond = c("logistic", "weibull", "loglogistic", "lognormal"), df = 3L, knots = NULL,
                      constraint = FALSE, censoring_grid = 0L, init = NULL, max_iter = 500L, tol = 1e-10) {
  cl <- match.call()
  estimation <- match.arg(estimation); cond <- match.arg(cond)
  f <- .sx_frame(formula, data, what = "directcif()")
  .directcif_core(cl, f, baseline, link, estimation, as.integer(cause), cond, df, knots, constraint, censoring_grid, init, max_iter, tol)
}

# the fit itself, from a model frame list as returned by .sx_frame (also used to refit without covariates)
.directcif_core <- function(cl, f, baseline, link, estimation, cause, cond, df, knots, constraint, censoring_grid, init, max_iter, tol) {
  n <- length(f$time); p <- ncol(f$X); K <- max(f$status)
  if (K < 1L) stop("no events", call. = FALSE)
  bad <- setdiff(baseline, names(.base_codes))
  if (length(bad)) stop("unknown baseline: ", paste(bad, collapse = ", "), "; choose from ", paste(names(.base_codes), collapse = ", "), call. = FALSE)
  if (estimation != "joint")
    return(.directcif_single(cl, f, estimation, as.integer(cause), cond, baseline, link, df, knots, censoring_grid, init, max_iter, tol))
  baseline <- rep_len(baseline, K)
  link <- rep_len(as.list(link), K)
  alpha_free <- integer(K); alpha_value <- numeric(K); link_name <- character(K)
  for (k in seq_len(K)) {
    lk <- .parse_link(link[[k]])
    alpha_free[k] <- lk$free; alpha_value[k] <- lk$value; link_name[k] <- lk$name
  }
  code <- unname(.base_codes[baseline])
  kn <- vector("list", K)
  for (k in seq_len(K)) {
    kn[[k]] <- numeric(0)
    if (baseline[k] == "spline") {
      if (!is.null(knots)) kn[[k]] <- as.double(knots[[k]])
      else {
        r <- .C("SX_direct_knots", as.integer(n), f$time, f$status, as.integer(k), as.integer(df),
                knots = double(df + 1L), nk = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
        .sx_rc(r$rc, "directcif (knots)")
        kn[[k]] <- r$knots[seq_len(r$nk)]
      }
    }
  }
  nk <- as.integer(lengths(kn))
  nm <- character(0); index <- NULL
  for (k in seq_len(K)) {
    bn <- .base_names(baseline[k], nk[k])
    terms_k <- c(bn, if (alpha_free[k]) "alpha", colnames(f$X))
    role_k <- c(rep("baseline", length(bn)), if (alpha_free[k]) "alpha", rep("beta", p))
    nm <- c(nm, paste0(k, ":", terms_k))
    index <- rbind(index, data.frame(cause = k, term = terms_k, role = role_k, stringsAsFactors = FALSE))
  }
  P <- length(nm)
  if (!is.null(init) && length(init) != P) stop("init must have length ", P, call. = FALSE)
  r <- .C("SX_direct_fit", as.integer(n), as.integer(p), as.integer(K), f$time, f$status, .rm(f$X),
          as.integer(alpha_free), as.double(alpha_value), as.integer(code), nk, as.double(unlist(kn)),
          as.integer(isTRUE(constraint)), as.integer(!is.null(init)), if (is.null(init)) double(P) else as.double(init),
          as.integer(max_iter), as.double(tol), theta = double(P), var = double(P * P),
          info_d = double(3), info_i = integer(4), diag = double(5), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "directcif")
  dg <- list(n_exceed_range = r$diag[1L], max_sum_range = r$diag[2L], n_exceed_inf = r$diag[3L],
             max_sum_inf = r$diag[4L], first_crossing_time = r$diag[5L], grad_maxabs = r$info_d[3L])
  DirectCIFFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
               coefficients = stats::setNames(r$theta, nm), vcov = matrix(r$var, P, P, dimnames = list(nm, nm)),
               loglik = r$info_d[2L], n = n, nevent = as.integer(sum(f$status > 0L)), iterations = r$info_i[1L],
               converged = as.logical(r$info_i[2L]), flags = r$info_i[3L], time = f$time, status = f$status, X = f$X,
               K = as.integer(K), baseline = baseline, link = link_name, base_code = as.integer(code),
               alpha_free = as.integer(alpha_free), alpha_value = as.double(alpha_value), n_knots = nk, knots = kn,
               constraint = isTRUE(constraint), diagnostics = dg, index = index,
               estimation = "joint", interest = NA_integer_, cond = NA_character_, vcov_model = NULL)
}

# estimation = "weighted" or "constrained": one directly modelled cause
.directcif_single <- function(cl, f, estimation, cause, cond, baseline, link, df, knots, grid, init, max_iter, tol) {
  n <- length(f$time); p <- ncol(f$X); K <- max(f$status)
  if (length(cause) != 1L || is.na(cause) || !any(f$status == cause)) stop("no events of cause ", cause, call. = FALSE)
  pick <- function(v) if (length(v) >= cause && length(v) > 1L) v[[cause]] else v[[1L]]
  baseline <- pick(baseline)
  lk <- .parse_link(pick(as.list(link)))
  code <- unname(.base_codes[[baseline]])
  kn <- numeric(0)
  bn <- NULL
  if (estimation == "constrained") {
    if (K != 2L) stop("estimation = \"constrained\" needs exactly two causes (status 1 and 2)", call. = FALSE)
    if (baseline %in% c("weibull", "spline")) stop("estimation = \"constrained\" needs a baseline with a plateau (not weibull or spline)", call. = FALSE)
  } else if (baseline == "spline") {
    if (!is.null(knots)) kn <- as.double(if (is.list(knots)) pick(knots) else knots)
    else {
      r <- .C("SX_direct_knots", as.integer(n), f$time, f$status, cause, as.integer(df),
              knots = double(df + 1L), nk = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
      .sx_rc(r$rc, "directcif (knots)")
      kn <- r$knots[seq_len(r$nk)]
    }
  }
  nk <- length(kn)
  bn <- .base_names(baseline, nk)
  terms1 <- c(bn, if (lk$free) "alpha", colnames(f$X))
  role1 <- c(rep("baseline", length(bn)), if (lk$free) "alpha", rep("beta", p))
  nm <- paste0(cause, ":", terms1)
  index <- data.frame(cause = cause, term = terms1, role = role1, stringsAsFactors = FALSE)
  if (estimation == "constrained") {
    other <- setdiff(1:2, cause)
    gn <- .cond_names(cond)
    nm <- c(nm, paste0(other, ":", gn))
    index <- rbind(index, data.frame(cause = other, term = gn, role = "conditional", stringsAsFactors = FALSE))
  }
  P <- length(nm)
  if (!is.null(init) && length(init) != P) stop("init must have length ", P, call. = FALSE)
  ini <- if (is.null(init)) double(P) else as.double(init)
  vm <- NULL; nsup <- NA_integer_
  if (estimation == "weighted") {
    r <- .C("SX_subdist_fit", as.integer(n), as.integer(p), f$time, f$status, .rm(f$X), cause,
            as.integer(lk$free), as.double(lk$value), as.integer(code), as.integer(nk), as.double(kn), as.integer(grid),
            as.integer(!is.null(init)), ini, as.integer(max_iter), as.double(tol),
            theta = double(P), var = double(P * P), var_model = double(P * P), info_d = double(3), info_i = integer(4),
            n_support = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "directcif (weighted)")
    vm <- matrix(r$var_model, P, P, dimnames = list(nm, nm)); nsup <- r$n_support
    nev <- as.integer(sum(f$status == cause))
  } else {
    st12 <- ifelse(f$status == 0L, 0L, ifelse(f$status == cause, 1L, 2L))
    r <- .C("SX_constrained_fit", as.integer(n), as.integer(p), f$time, as.integer(st12), .rm(f$X),
            as.integer(lk$free), as.double(lk$value), as.integer(code), as.integer(.cond_codes[[cond]]),
            as.integer(!is.null(init)), ini, as.integer(max_iter), as.double(tol),
            theta = double(P), var = double(P * P), info_d = double(3), info_i = integer(4), rc = integer(1),
            NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "directcif (constrained)")
    nev <- as.integer(sum(f$status > 0L))
  }
  dg <- list(n_exceed_range = 0, max_sum_range = NA_real_, n_exceed_inf = 0, max_sum_inf = NA_real_,
             first_crossing_time = NA_real_, grad_maxabs = r$info_d[3L], n_support = nsup)
  DirectCIFFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
               coefficients = stats::setNames(r$theta, nm), vcov = matrix(r$var, P, P, dimnames = list(nm, nm)),
               loglik = r$info_d[2L], n = n, nevent = nev, iterations = r$info_i[1L],
               converged = as.logical(r$info_i[2L]), flags = r$info_i[3L], time = f$time, status = f$status, X = f$X,
               K = as.integer(K), baseline = baseline, link = lk$name, base_code = as.integer(code),
               alpha_free = as.integer(lk$free), alpha_value = as.double(lk$value), n_knots = as.integer(nk), knots = list(kn),
               constraint = FALSE, diagnostics = dg, index = index,
               estimation = estimation, interest = cause, cond = if (estimation == "constrained") cond else NA_character_,
               vcov_model = vm)
}

.direct_predict <- function(object, Xn, times, cause, se) {
  nn <- nrow(Xn); nt <- length(times); p <- ncol(object@X)
  if (is.null(cause)) cause <- if (object@estimation == "joint") 1L else object@interest
  cause <- as.integer(cause)
  usevar <- isTRUE(se) && all(is.finite(object@vcov))
  v <- if (usevar) as.double(object@vcov) else double(1)
  if (object@estimation == "joint") {
    if (cause < 1L || cause > object@K) stop("cause must be between 1 and ", object@K, call. = FALSE)
    r <- .C("SX_direct_predict", as.integer(p), object@K, object@alpha_free, as.double(object@alpha_value), object@base_code,
            object@n_knots, as.double(unlist(object@knots)), as.double(object@coefficients), as.integer(usevar), v,
            cause, as.integer(nn), .rm(Xn), as.integer(nt), as.double(times),
            cif = double(nn * nt), se = double(nn * nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  } else if (object@estimation == "weighted") {
    if (cause != object@interest) stop("this model was fitted for cause ", object@interest, " only", call. = FALSE)
    r <- .C("SX_direct_predict", as.integer(p), 1L, object@alpha_free, as.double(object@alpha_value), object@base_code,
            object@n_knots, as.double(unlist(object@knots)), as.double(object@coefficients), as.integer(usevar), v,
            1L, as.integer(nn), .rm(Xn), as.integer(nt), as.double(times),
            cif = double(nn * nt), se = double(nn * nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  } else {
    if (!cause %in% 1:2) stop("cause must be 1 or 2", call. = FALSE)
    r <- .C("SX_constrained_predict", as.integer(p), object@alpha_free, as.double(object@alpha_value), object@base_code,
            as.integer(.cond_codes[[object@cond]]), as.double(object@coefficients), as.integer(usevar), v,
            if (cause == object@interest) 1L else 2L, as.integer(nn), .rm(Xn), as.integer(nt), as.double(times),
            cif = double(nn * nt), se = double(nn * nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  }
  .sx_rc(r$rc, "directcif prediction")
  dn <- list(NULL, paste0("t=", times))
  cif <- .matrix_rm(r$cif, nn, nt, dn)
  if (!isTRUE(se)) return(cif)
  list(cif = cif, se = if (usevar) .matrix_rm(r$se, nn, nt, dn) else matrix(NA_real_, nn, nt, dimnames = dn))
}

S7::method(predict, DirectCIFFit) <- function(object, newdata = NULL, times, cause = NULL, se = FALSE, ...) {
  .direct_predict(object, .sx_newX(object, newdata), times, cause, se)
}
