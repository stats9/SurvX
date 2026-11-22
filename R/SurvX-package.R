#' SurvX: survival and competing risks analysis
#'
#' @description
#' SurvX fits the main non-parametric, semi-parametric and parametric methods for time-to-event
#' data with one event type or with competing risks, in a fast C++ core with one consistent
#' interface. Every result names the published methods it rests on (the \code{Methods:} line of every
#' print, and \code{\link{sx_references}}), and every function has been checked against established
#' reference packages (see the \code{Validation} section of each help page).
#' @section Which function answers my question?:
#' \tabular{ll}{
#'   \strong{Question} \tab \strong{Function} \cr
#'   How does survival look over time, and does it differ between groups? \tab \code{\link{km}},
#'     \code{\link{logrank_test}} \cr
#'   What is the risk of each cause, when other events compete? \tab \code{\link{aj}},
#'     \code{\link{gray_test}} \cr
#'   How do covariates change the event rate (hazard ratios)? \tab \code{\link{cox}} \cr
#'   Is the hazard ratio constant over time? \tab \code{\link{ph_test}} \cr
#'   How does an effect change over time? \tab \code{\link{time_effects}} \cr
#'   Adjusted absolute risks and causal contrasts \tab \code{\link{adjusted_risk}} \cr
#'   Event times known only between two visits (interval censoring) \tab \code{\link{Interval}}, \code{\link{turnbull}}, \code{\link{parsurv}}, \code{\link{directcif}} \cr
#'   Clustered subjects with shared unmeasured risk (frailty) \tab \code{\link{cox}} or \code{\link{parsurv}} with \code{+ \link{frailty}(id)} \cr
#'   How do covariates change the rate of each cause? \tab \code{\link{csc}} \cr
#'   How do covariates change the risk (cumulative incidence) of one cause? \tab \code{\link{fg}},
#'     \code{\link{directcif}} \cr
#'   Smooth, parametric survival curves, time ratios, extrapolation \tab \code{\link{parsurv}} \cr
#'   What proportion is cured, and when do the others fail? \tab \code{\link{cure}} \cr
#'   An effect measure without proportional hazards \tab \code{\link{rmst}} \cr
#'   What does the hazard look like? Which distribution fits? \tab \code{\link{hazard_np}},
#'     \code{\link{hazard_plot}}, \code{\link{dist_plot}} \cr
#'   Absolute risks for given patients \tab \code{\link{predict.SurvX}}, \code{\link{cif_plot}} \cr
#'   How good are the predictions? Which model predicts best? \tab \code{\link{brier_score}},
#'     \code{\link{compare_models}}, \code{\link{cv_compare}} \cr
#'   Residuals and model checks \tab \code{\link{residuals.SurvX}}, \code{\link{ph_test}} \cr
#' }
#' @section Getting started:
#' \enumerate{
#'   \item Code the outcome with \code{\link{Event}(time, status)}: status 0 = censored, 1 = event,
#'     or 1, 2, ... = causes with competing risks; add \code{entry =} for delayed entry.
#'   \item Describe the data: \code{km()} (one event type) or \code{aj()} (competing risks).
#'   \item Fit a model, read its print (effects with intervals, \code{Reading:} lines that state
#'     the effects in words, and the \code{Methods:} line), and check its assumptions.
#'   \item Report absolute risks with \code{predict()} or \code{cif_plot()} next to the effects.
#' }
#' All model objects work with the usual generics (\code{print}, \code{summary}, \code{coef},
#' \code{vcov}, \code{confint}, \code{logLik}, \code{AIC}, \code{predict}, \code{residuals},
#' \code{plot}); see \code{\link{SurvX-methods}}. \code{\link{estimates}} returns any table of
#' effects as a data frame.
#' @section Conventions:
#' \itemize{
#'   \item Times that differ only by rounding error are treated as tied, as in \pkg{survival}
#'     (\code{timefix} of \code{\link{Event}}).
#'   \item Tied event times use the Efron method in Cox models by default, as in \pkg{survival}.
#'   \item Confidence intervals are 95\% Wald intervals unless stated; probabilities use
#'     transformed scales (log or complementary log-log) so that the limits stay in (0, 1).
#'   \item The fit of every iterative model is checked: a model that did not converge says so in its
#'     print, with the reason.
#' }
#' @references
#' The references of each method are listed on the help page of its function, grouped by what they
#' are used for, and printed by \code{\link{sx_references}}.
#' @seealso
#' \code{vignette(package = "SurvX")} for worked examples; \url{https://github.com/stats9/SurvX}
#' for the source code and the validation suite.
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(arm = rep(c("A", "B"), each = n / 2), age = round(rnorm(n, 60, 10)))
#' t1 <- rexp(n, ifelse(d$arm == "B", 0.06, 0.1)); t2 <- rexp(n, 0.04); cens <- runif(n, 2, 10)
#' d$time <- pmin(t1, t2, cens)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#'
#' aj(Event(time, status) ~ arm, d)                      # describe
#' fg(Event(time, status) ~ arm + age, d)                # risk of cause 1
#' csc(Event(time, status) ~ arm + age, d)               # rates of both causes
#' sx_references("fg")                                   # the literature behind fg()
"_PACKAGE"
