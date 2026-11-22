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
#' a weighted estimator for one cause, and a constrained model for two or more causes whose
#' incidences add up to one.
#'
#' \strong{Use it when} smooth, parametric cumulative incidences are wanted (extrapolation, the
#' long-term incidence \eqn{F_k(\infty)}, health-economic models), when the effect is better described
#' on the odds scale, or when all causes are to be modelled together with coherent predictions.
#' \strong{Use something else when} no assumption on the shape of the incidence is wanted:
#' \code{\link{fg}}; or the question is about rates: \code{\link{csc}}.
#' @param formula \code{Event(time, status) ~ covariates} (or \code{~ 1}); status 0 = censored,
#'   1, 2, ..., K = cause. Delayed entry is not supported. For interval-censored times,
#'   \code{Interval(left, right, cause = cause) ~ covariates} (see the section Interval-censored
#'   data).
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
#'   the other causes is needed. \code{"constrained"}: all causes, with incidences that add up to one
#'   at infinity for every covariate value (Shi, Cheng and Jeong 2013 for two causes; extended here to
#'   any number of causes). See Details.
#' @param cause the cause of interest for \code{"weighted"} and \code{"constrained"}.
#' @param cond \code{"constrained"} only: the distribution \eqn{G_k} of the time to a competing
#'   event given that it occurs: \code{"logistic"} (modified logistic, as in Shi et al.; default),
#'   \code{"weibull"}, \code{"loglogistic"} or \code{"lognormal"}. One value for all competing
#'   causes, or one per competing cause (in the order of the cause codes, the cause of interest left
#'   out).
#' @param timing \code{"constrained"} only: covariates of the time to a competing event, through an
#'   accelerated failure time model for \eqn{G_k} (see the section Covariates of the time
#'   distributions). A one-sided formula such as \code{~ age} for every competing cause, or a list
#'   of them named by the competing causes, such as \code{list("2" = ~ age, "3" = ~ age + sex)}.
#'   Its terms must also be terms of the model formula. \code{NULL} (default): \eqn{G_k} has no
#'   covariates, the model of Shi, Cheng and Jeong (2013).
#' @param shape \code{"constrained"} only: covariates of the shape of \eqn{G_k} (the same forms as
#'   \code{timing}). With \code{timing}, the time ratio then changes with time.
#' @param df degrees of freedom of the spline baseline.
#' @param knots optional list with the knots (log time, boundary knots first and last) of every
#'   cause, for the spline baseline.
#' @param constraint \code{"joint"} only: keep the sum of the fitted incidences below 1 within
#'   follow-up (see Details).
#' @param censoring_grid \code{"weighted"} only: coarsen the estimated censoring distribution to at
#'   most this many support points (0 = exact). With the PH link the exact fit is linear in the
#'   sample size; with other links its cost grows with the number of competing events times the
#'   number of support points (a few seconds at \eqn{n = 10000}), and a grid can shorten it.
#' @param variance \code{"weighted"} only: the standard errors. \code{"ipcw"} (default): sandwich
#'   variance that includes the effect of estimating the censoring distribution (see Details).
#'   \code{"robust"}: sandwich variance with the censoring distribution treated as known.
#'   \code{"model"}: inverse information of the weighted pseudo-likelihood. All three matrices are
#'   kept in \code{fit@diagnostics} (\code{vcov_ipcw}, \code{vcov_robust}) and
#'   \code{fit@vcov_model}. Ignored for the other estimators.
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
#'     Crowther 2017). The weights depend on the Kaplan-Meier estimate of the censoring
#'     distribution, and by default (\code{variance = "ipcw"}) the standard errors account for it:
#'     the sandwich estimator uses, for every subject, the derivative of the total score with respect
#'     to a case weight of that subject with the Kaplan-Meier estimate recomputed (the infinitesimal
#'     jackknife), the discrete form of the correction term of Fine and Gray (1999, Sec. 3.2). It is
#'     computed on the exact censoring distribution also when \code{censoring_grid > 0}. The
#'     correction mostly concerns the baseline parameters, and with them the predicted and long-term
#'     incidences; for the regression coefficients it is usually small. \code{logLik} is a weighted
#'     pseudo-likelihood: use it only to compare links or baselines of weighted fits.}
#'   \item{\code{"constrained"}}{causes with asymptotes that add up to one. With two causes this is
#'     the model of Shi, Cheng and Jeong (2013): the cause of interest as above and
#'     \eqn{F_2(t \mid z) = \{1 - F_1(\infty \mid z)\} G(t)}. With \eqn{K \ge 3} causes the
#'     remaining long-run risk is shared among the competing causes by a multinomial logit,
#'     \deqn{F_k(t \mid z) = \{1 - F_1(\infty \mid z)\}\, \pi_k(z)\, G_k(t), \qquad
#'       \pi_k(z) = \frac{\exp(\gamma_k + \delta_k^\top z)}{\sum_{j} \exp(\gamma_j + \delta_j^\top z)},}
#'     with the first competing cause as reference (\eqn{\gamma = 0, \delta = 0}), so that
#'     \eqn{\sum_k F_k(\infty \mid z) = 1} for every \eqn{z}. The cause of interest keeps the direct
#'     regression model and its interpretation; the share of each competing cause can depend on the
#'     covariates, in the manner of the mixture model of Larson and Dinse (1985). By default
#'     \eqn{G_k} has no covariates; \code{timing} and \code{shape} let covariates change when the
#'     competing events occur (section Covariates of the time distributions). Needs a baseline with a
#'     plateau for the cause of interest (not \code{"weibull"}, not \code{"spline"}). Appropriate
#'     when every subject eventually fails from one of the causes.}
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
#'   \item{Shares (constrained model, three or more causes)}{\code{share_} coefficients of cause
#'     \eqn{k}: \eqn{e^{\delta}} is an \strong{odds ratio of failing from cause \eqn{k} rather than
#'     from the reference competing cause}, among subjects who do not fail from the cause of interest;
#'     \code{share_(Intercept)} is that log odds at \eqn{z = 0}. The effect of a covariate on the
#'     incidence of a competing cause combines its effect on the cause of interest (through
#'     \eqn{1 - F_1(\infty \mid z)}) and on the share; \code{predict()} gives the incidences.}
#'   \item{Time ratios (constrained model, \code{timing})}{\code{time_} coefficients of cause
#'     \eqn{k}: \eqn{e^{\eta}} is a \strong{time ratio}: among subjects who fail from cause \eqn{k},
#'     a one-unit increase of the covariate multiplies the time to the event by \eqn{e^{\eta}}
#'     (TR > 1: later events). It does not change how many fail from cause \eqn{k} in the long run.}
#'   \item{Shape effects (\code{shape})}{\code{shape_} coefficients: \eqn{e^{\zeta}} multiplies the
#'     shape parameter of \eqn{G_k} (\eqn{\kappa} for the Weibull and log-logistic forms,
#'     \eqn{\sigma} for the log-normal form, the steepness for the modified logistic form).}
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
#' The variance with the censoring term (\code{variance = "ipcw"}) agrees with an independent
#' implementation of the formula to \eqn{10^{-10}} for every baseline family and link, and with the
#' exact infinitesimal jackknife, the numerical derivative of the estimate with respect to case
#' weights with the Kaplan-Meier estimate recomputed, to the accuracy of the numerical derivative.
#' The constrained model with \eqn{K \ge 3} causes reproduces the two-cause model exactly for
#' \eqn{K = 2}; its log-likelihood agrees with independent code to \eqn{10^{-12}}, its score and
#' information with numerical differentiation to \eqn{10^{-9}}; in 300 simulated data sets with
#' three causes (\eqn{n = 800}) the bias of every parameter was below 0.01, the ratio of mean
#' standard error to empirical standard deviation 0.94 to 1.03, and the coverage of the 95\%
#' intervals 93.7 to 95.7\%.
#' @section Covariates of the time distributions:
#' In the constrained model the covariates act on the long-run risk of the cause of interest and on
#' the shares of the competing causes, but by default not on when the competing events occur. When
#' they do (older patients die of other causes sooner, for example), \code{timing} makes
#' \eqn{G_k} an accelerated failure time model,
#' \deqn{G_k(t \mid z) = G_{k0}\{t\, e^{-\eta_k^\top z}\},}
#' and \code{shape} lets the shape of \eqn{G_{k0}} depend on covariates as well. Every
#' \eqn{G_k(\cdot \mid z)} is a distribution function, so the incidences still add up to one at
#' infinity for every covariate value, and the long-run risks and shares keep their meaning; only the
#' timing of the competing events changes. The model without these covariates is nested in the model
#' with them: it is fitted first, and the fit with them starts from it, so its log-likelihood is never
#' lower; compare the two with \code{AIC()} or a likelihood ratio test. The parameters of \eqn{G_k}
#' describe \eqn{G_k} at zero values of its covariates; the fit itself uses centred covariates, which
#' does not change the estimates. Each time ratio or shape effect is a parameter of the time
#' distribution of one cause and is informed mainly by the events of that cause (and by
#' the censored subjects): a warning is given when a cause has fewer than 10 events per parameter of its time distribution. A shape
#' effect needs more events than a time ratio.
#'
#' Validation: the log-likelihood agrees with independent code written from the model definition to
#' \eqn{10^{-15}} for all four forms of \eqn{G_k}, the score and information with numerical
#' differentiation to \eqn{10^{-8}}, and the maximum and its covariance with a general-purpose
#' optimiser of the same likelihood; without these covariates the results are identical, digit for
#' digit, to those of the model without them. In simulations from the model (\eqn{n = 1000}, two to
#' four causes, time ratios and shape effects) the bias was below 0.26 standard deviations for every
#' parameter, the ratio of mean standard error to empirical standard deviation 0.83 to 1.06 (0.93 to
#' 1.10 with \eqn{n = 4000}) and the coverage of the 95\% intervals 91 to 97\%.
#' @section Interval-censored data:
#' When the time of the event is only known to lie between two visits but its type is known
#' (\code{\link{Interval}(left, right, cause = )}), the joint likelihood uses
#' \eqn{F_k(R) - F_k(L)} for an event of cause \eqn{k} in \eqn{(L, R]} (\eqn{F_k(R)} when
#' \eqn{L = 0}), \eqn{1 - \sum_k F_k(L)} for a subject event-free at the last visit \eqn{L} and the
#' density for exact times (Hudgens, Satten and Longini 2001 for the likelihood; the direct parametric
#' models of Jeong and Fine 2006, 2007 for \eqn{F_k}). Only \code{estimation = "joint"} is available
#' (the weighted and constrained estimators need exact times). The fit starts from the fit to the
#' midpoints of the intervals. Predictions, \code{estimates()} and \code{\link{cif_plot}} work as
#' usual; the Aalen-Johansen reference and the numbers at risk are not drawn, as they are not defined
#' for these data.
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
#' Larson MG, Dinse GE (1985) A mixture model for the regression analysis of competing risks data.
#' J R Stat Soc C 34:201-211.
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
#' \code{\link{Interval}}, \code{\link{turnbull}};
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
#'
#' ## 5. three causes: the competing causes share the remaining long-run risk
#' d3 <- d
#' d3$status[d3$status == 2 & runif(n) < 0.4] <- 3
#' f3 <- directcif(Event(time, status) ~ x, d3, estimation = "constrained", cause = 1,
#'                 cond = "weibull")
#' f3
#' rowSums(sapply(1:3, function(k) predict(f3, nd, times = Inf, cause = k)))
#'
#' ## 6. x also changes when the competing events occur (time ratio of cause 2 and 3 events)
#' f4 <- directcif(Event(time, status) ~ x, d3, estimation = "constrained", cause = 1,
#'                 cond = "weibull", timing = ~ x)
#' f4
#' AIC(f3, f4)
#' @export
directcif <- function(formula, data, baseline = "gompertz", link = "PH",
                      estimation = c("joint", "weighted", "constrained"), cause = 1L,
                      cond = c("logistic", "weibull", "loglogistic", "lognormal"), timing = NULL, shape = NULL,
                      df = 3L, knots = NULL,
                      constraint = FALSE, censoring_grid = 0L, variance = c("ipcw", "robust", "model"),
                      init = NULL, max_iter = 500L, tol = 1e-10) {
  cl <- match.call()
  estimation <- match.arg(estimation); variance <- match.arg(variance)
  cond <- if (missing(cond)) "logistic"
          else vapply(cond, match.arg, character(1), choices = c("logistic", "weibull", "loglogistic", "lognormal"), USE.NAMES = FALSE)
  f <- .sx_frame(formula, data, what = "directcif()", interval_ok = TRUE)
  .directcif_core(cl, f, baseline, link, estimation, as.integer(cause), cond, df, knots, constraint, censoring_grid, init, max_iter, tol,
                  variance = variance, timing = timing, shape = shape)
}

# the fit itself, from a model frame list as returned by .sx_frame (also used to refit without covariates)
.directcif_core <- function(cl, f, baseline, link, estimation, cause, cond, df, knots, constraint, censoring_grid, init, max_iter, tol,
                            variance = "ipcw", timing = NULL, shape = NULL) {
  n <- length(f$time); p <- ncol(f$X); K <- max(f$status)
  if (K < 1L) stop("no events", call. = FALSE)
  bad <- setdiff(baseline, names(.base_codes))
  if (length(bad)) stop("unknown baseline: ", paste(bad, collapse = ", "), "; choose from ", paste(names(.base_codes), collapse = ", "), call. = FALSE)
  ic <- !is.null(f$time2)
  if (ic && estimation != "joint")
    stop("interval-censored data (Interval()) are fitted with estimation = \"joint\" (the full likelihood); ",
         "the weighted and constrained estimators need exact event times", call. = FALSE)
  if ((!is.null(timing) || !is.null(shape)) && estimation != "constrained")
    stop("timing and shape (covariates of the time distributions G_k) are for estimation = \"constrained\"", call. = FALSE)
  tk <- if (ic) .ic_mid(f$time, f$time2) else f$time       # times for knots and starting values
  if (estimation != "joint")
    return(.directcif_single(cl, f, estimation, as.integer(cause), cond, baseline, link, df, knots, censoring_grid, init, max_iter, tol,
                             variance, timing = timing, shape = shape))
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
        r <- .C("SX_direct_knots", as.integer(n), as.double(tk), f$status, as.integer(k), as.integer(df),
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
  if (ic) {
    # start: the fit to the midpoints of the intervals (unless init is given)
    if (is.null(init)) {
      r0 <- .C("SX_direct_fit", as.integer(n), as.integer(p), as.integer(K), as.double(tk), f$status, .rm(f$X),
               as.integer(alpha_free), as.double(alpha_value), as.integer(code), nk, as.double(unlist(kn)),
               as.integer(isTRUE(constraint)), 0L, double(P), as.integer(max_iter), as.double(tol), theta = double(P), var = double(P * P),
               info_d = double(3), info_i = integer(4), diag = double(5), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
      if (r0$rc == 0L && all(is.finite(r0$theta)) && is.finite(r0$info_d[2L])) init <- r0$theta
    }
    r <- .C("SX_direct_fit_ic", as.integer(n), as.integer(p), as.integer(K), as.double(f$time), as.double(f$time2), f$status, .rm(f$X),
            as.integer(alpha_free), as.double(alpha_value), as.integer(code), nk, as.double(unlist(kn)),
            as.integer(isTRUE(constraint)), as.integer(!is.null(init)), if (is.null(init)) double(P) else as.double(init),
            as.integer(max_iter), as.double(tol), theta = double(P), var = double(P * P),
            info_d = double(3), info_i = integer(4), diag = double(5), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "directcif (interval-censored)")
  } else {
  r <- .C("SX_direct_fit", as.integer(n), as.integer(p), as.integer(K), f$time, f$status, .rm(f$X),
          as.integer(alpha_free), as.double(alpha_value), as.integer(code), nk, as.double(unlist(kn)),
          as.integer(isTRUE(constraint)), as.integer(!is.null(init)), if (is.null(init)) double(P) else as.double(init),
          as.integer(max_iter), as.double(tol), theta = double(P), var = double(P * P),
          info_d = double(3), info_i = integer(4), diag = double(5), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "directcif")
  }
  dg <- list(n_exceed_range = r$diag[1L], max_sum_range = r$diag[2L], n_exceed_inf = r$diag[3L],
             max_sum_inf = r$diag[4L], first_crossing_time = r$diag[5L], grad_maxabs = r$info_d[3L])
  DirectCIFFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
               coefficients = stats::setNames(r$theta, nm), vcov = matrix(r$var, P, P, dimnames = list(nm, nm)),
               loglik = r$info_d[2L], n = n, nevent = as.integer(sum(f$status > 0L)), iterations = r$info_i[1L],
               converged = as.logical(r$info_i[2L]), flags = r$info_i[3L],
               time = if (ic) ifelse(is.finite(f$time2), f$time2, f$time) else f$time, status = f$status, X = f$X,
               K = as.integer(K), baseline = baseline, link = link_name, base_code = as.integer(code),
               alpha_free = as.integer(alpha_free), alpha_value = as.double(alpha_value), n_knots = nk, knots = kn,
               constraint = isTRUE(constraint), diagnostics = dg, index = index,
               estimation = "joint", interest = NA_integer_, cond = NA_character_, vcov_model = NULL,
               interval = if (ic) list(left = f$time, right = f$time2, type = .ic_type(f$time, f$time2)) else NULL)
}

# estimation = "weighted" or "constrained": one directly modelled cause
.directcif_single <- function(cl, f, estimation, cause, cond, baseline, link, df, knots, grid, init, max_iter, tol,
                              variance = "ipcw", timing = NULL, shape = NULL) {
  n <- length(f$time); p <- ncol(f$X); K <- max(f$status)
  if (length(cause) != 1L || is.na(cause) || !any(f$status == cause)) stop("no events of cause ", cause, call. = FALSE)
  pick <- function(v) if (length(v) >= cause && length(v) > 1L) v[[cause]] else v[[1L]]
  baseline <- pick(baseline)
  lk <- .parse_link(pick(as.list(link)))
  code <- unname(.base_codes[[baseline]])
  kn <- numeric(0)
  bn <- NULL
  has_cc <- FALSE
  if (estimation == "constrained") {
    if (K < 2L) stop("estimation = \"constrained\" needs at least two causes", call. = FALSE)
    miss <- setdiff(seq_len(K), unique(f$status[f$status > 0L]))
    if (length(miss)) stop("estimation = \"constrained\": no events of cause ", paste(miss, collapse = ", "), call. = FALSE)
    if (length(cond) != 1L && length(cond) != K - 1L)
      stop("cond must have length 1 or ", K - 1L, " (one distribution per competing cause)", call. = FALSE)
    cond <- rep_len(cond, K - 1L)
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
    other <- setdiff(seq_len(K), cause)                 # competing causes; the first is the reference of the shares
    for (j in seq_along(other)) {
      gn <- .cond_names(cond[j])
      nm <- c(nm, paste0(other[j], ":", gn))
      index <- rbind(index, data.frame(cause = other[j], term = gn, role = "conditional", stringsAsFactors = FALSE))
    }
    if (K > 2L) for (j in seq_along(other)[-1L]) {
      sn <- c("share_(Intercept)", if (p) paste0("share_", colnames(f$X)))
      nm <- c(nm, paste0(other[j], ":", sn))
      index <- rbind(index, data.frame(cause = other[j], term = sn, role = c("share", rep("share_beta", p)), stringsAsFactors = FALSE))
    }
    # covariates of the G_k: time ratios of every competing cause, then shape effects (as in the C++ core)
    wuse <- .cond_use(timing, f, other, "timing")
    vuse <- .cond_use(shape, f, other, "shape")
    for (part in list(list(u = wuse, pre = "time_", role = "timing"), list(u = vuse, pre = "shape_", role = "shape")))
      for (j in seq_along(other)) {
        cols <- which(part$u[j, ] == 1L)
        if (!length(cols)) next
        tn <- paste0(part$pre, colnames(f$X)[cols])
        nm <- c(nm, paste0(other[j], ":", tn))
        index <- rbind(index, data.frame(cause = other[j], term = tn, role = part$role, stringsAsFactors = FALSE))
      }
    has_cc <- any(wuse == 1L) || any(vuse == 1L)
    if (has_cc) {
      for (j in seq_along(other)) {
        npj <- 2L + sum(wuse[j, ]) + sum(vuse[j, ])
        nej <- sum(f$status == other[j])
        if (npj > 2L && nej < 10L * npj)
          warning("cause ", other[j], ": ", nej, " events for the ", npj, " parameters of its time distribution; ",
                  "the time ratios and shape effects may be poorly estimated", call. = FALSE)
      }
    }
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
    vrob <- matrix(r$var, P, P, dimnames = list(nm, nm))
    vip <- matrix(NA_real_, P, P, dimnames = list(nm, nm))
    if (all(is.finite(r$var_model))) {
      q <- .C("SX_subdist_var_ipcw", as.integer(n), as.integer(p), f$time, f$status, .rm(f$X), cause,
              as.integer(lk$free), as.double(lk$value), as.integer(code), as.integer(nk), as.double(kn),
              as.double(r$theta), as.double(r$var_model), var = double(P * P), rc = integer(1),
              NAOK = TRUE, PACKAGE = "SurvX")
      if (q$rc == 0L) vip[] <- q$var
    }
    if (variance == "ipcw" && !all(is.finite(vip))) {
      warning("the variance with the censoring term could not be computed; the robust variance is used", call. = FALSE)
      variance <- "robust"
    }
    r$var <- switch(variance, ipcw = vip, robust = vrob, model = vm)
  } else if (has_cc) {
    stk <- ifelse(f$status == 0L, 0L, ifelse(f$status == cause, 1L, match(f$status, other) + 1L))
    r <- .C("SX_constrainedK_fit_cov", as.integer(n), as.integer(p), f$time, as.integer(stk), .rm(f$X), as.integer(K),
            as.integer(lk$free), as.double(lk$value), as.integer(code), as.integer(unname(.cond_codes[cond])),
            as.integer(t(wuse)), as.integer(t(vuse)),
            as.integer(!is.null(init)), ini, as.integer(max_iter), as.double(tol),
            theta = double(P), var = double(P * P), info_d = double(3), info_i = integer(4), rc = integer(1),
            NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "directcif (constrained, covariates of G)")
    nev <- as.integer(sum(f$status > 0L))
  } else if (K > 2L) {
    stk <- ifelse(f$status == 0L, 0L, ifelse(f$status == cause, 1L, match(f$status, other) + 1L))
    r <- .C("SX_constrainedK_fit", as.integer(n), as.integer(p), f$time, as.integer(stk), .rm(f$X), as.integer(K),
            as.integer(lk$free), as.double(lk$value), as.integer(code), as.integer(unname(.cond_codes[cond])),
            as.integer(!is.null(init)), ini, as.integer(max_iter), as.double(tol),
            theta = double(P), var = double(P * P), info_d = double(3), info_i = integer(4), rc = integer(1),
            NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "directcif (constrained)")
    nev <- as.integer(sum(f$status > 0L))
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
  if (estimation == "weighted") { dg$variance <- variance; dg$vcov_robust <- vrob; dg$vcov_ipcw <- vip }
  DirectCIFFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
               coefficients = stats::setNames(r$theta, nm), vcov = matrix(r$var, P, P, dimnames = list(nm, nm)),
               loglik = r$info_d[2L], n = n, nevent = nev, iterations = r$info_i[1L],
               converged = as.logical(r$info_i[2L]), flags = r$info_i[3L], time = f$time, status = f$status, X = f$X,
               K = as.integer(K), baseline = baseline, link = lk$name, base_code = as.integer(code),
               alpha_free = as.integer(lk$free), alpha_value = as.double(lk$value), n_knots = as.integer(nk), knots = list(kn),
               constraint = FALSE, diagnostics = dg, index = index,
               estimation = estimation, interest = cause, cond = if (estimation == "constrained") cond else NA_character_,
               vcov_model = vm,
               cond_cov = if (has_cc) list(wuse = wuse, vuse = vuse) else NULL)
}

# Covariates of the time distributions G_k of the constrained model (arguments timing and shape):
# a one-sided formula for every competing cause, or a list of them named by the competing causes.
# Their terms must be terms of the model formula, so that covariate profiles, newdata and averages
# over the sample are defined once for the whole model. Returns the (K - 1) x p flags of the
# columns of the design matrix, rows in the order of the competing causes.
.cond_use <- function(spec, f, other, what) {
  p <- ncol(f$X)
  use <- matrix(0L, length(other), p, dimnames = list(as.character(other), colnames(f$X)))
  if (is.null(spec)) return(use)
  if (p == 0L) stop(what, " needs covariates in the model formula", call. = FALSE)
  if (inherits(spec, "formula")) {
    use[, .cond_columns(spec, f, what)] <- 1L
    return(use)
  }
  if (!is.list(spec) || is.null(names(spec)) || any(!nzchar(names(spec))))
    stop(what, " must be a one-sided formula such as ~ age, or a list of them named by the competing causes, ",
         "such as list(\"2\" = ~ age, \"3\" = ~ age + sex)", call. = FALSE)
  j <- match(names(spec), as.character(other))
  if (anyNA(j))
    stop(what, ": ", paste(names(spec)[is.na(j)], collapse = ", "), " is not a competing cause (the competing causes are ",
         paste(other, collapse = ", "), ")", call. = FALSE)
  if (anyDuplicated(j)) stop(what, ": a competing cause is named twice", call. = FALSE)
  for (i in seq_along(spec)) use[j[i], .cond_columns(spec[[i]], f, what)] <- 1L
  use
}

# columns of the design matrix of the terms of a one-sided formula; every term must be a term of
# the model formula (interactions match in any order of their variables)
.cond_columns <- function(fml, f, what) {
  if (!inherits(fml, "formula") || length(fml) != 2L)
    stop(what, " must be a one-sided formula such as ~ age (or a list of them)", call. = FALSE)
  lab <- attr(stats::terms(fml), "term.labels")
  if (!length(lab)) return(integer(0))
  norm <- function(x) vapply(strsplit(gsub(" ", "", x, fixed = TRUE), ":", fixed = TRUE),
                             function(v) paste(sort(v), collapse = ":"), character(1))
  m <- match(norm(lab), norm(attr(f$terms, "term.labels")))
  if (anyNA(m))
    stop(what, ": ", paste(lab[is.na(m)], collapse = ", "), " must also be a term of the model formula", call. = FALSE)
  which(attr(f$X, "assign") %in% m)
}

# the covariates of the G_k of a fit (NULL without them, also for fits made by older versions)
.cc_of <- function(fit) tryCatch(fit@cond_cov, error = function(e) NULL)

# predictions of the constrained model with K >= 3 causes or with covariates of the G_k:
# average = 1 gives the mean over the rows of Xn
.ck_predict <- function(f, kk, Xn, times, usevar, v, average) {
  nn <- nrow(Xn); nt <- length(times); p <- ncol(Xn); m <- if (average) 1L else nn
  cc <- .cc_of(f)
  if (is.null(cc))
    .C("SX_constrainedK_predict", as.integer(p), f@K, f@alpha_free, as.double(f@alpha_value), f@base_code,
       as.integer(unname(.cond_codes[f@cond])), as.double(f@coefficients), as.integer(usevar), v,
       as.integer(kk), as.integer(average), as.integer(nn), .rm(Xn), as.integer(nt), as.double(times),
       cif = double(m * nt), se = double(m * nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  else
    .C("SX_constrainedK_predict_cov", as.integer(p), f@K, f@alpha_free, as.double(f@alpha_value), f@base_code,
       as.integer(unname(.cond_codes[f@cond])), as.integer(t(cc$wuse)), as.integer(t(cc$vuse)),
       as.double(f@coefficients), as.integer(usevar), v,
       as.integer(kk), as.integer(average), as.integer(nn), .rm(Xn), as.integer(nt), as.double(times),
       cif = double(m * nt), se = double(m * nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
}

# index of a cause in the constrained model: 1 = cause of interest, 2..K = the others in order
.constrained_index <- function(object, cause) {
  if (cause == object@interest) 1L else match(cause, setdiff(seq_len(object@K), object@interest)) + 1L
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
    if (!cause %in% seq_len(object@K)) stop("cause must be between 1 and ", object@K, call. = FALSE)
    kk <- .constrained_index(object, cause)
    if (object@K > 2L || !is.null(.cc_of(object)))
      r <- .ck_predict(object, kk, Xn, times, usevar, v, 0L)
    else
      r <- .C("SX_constrained_predict", as.integer(p), object@alpha_free, as.double(object@alpha_value), object@base_code,
              as.integer(.cond_codes[[object@cond]]), as.double(object@coefficients), as.integer(usevar), v,
              kk, as.integer(nn), .rm(Xn), as.integer(nt), as.double(times),
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
