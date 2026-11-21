\name{time_effects}
\alias{time_effects}
\alias{TimeEffects}
\title{Time-varying effects: detect them, choose their form and refit the model}
\description{
Checks whether the effect of each covariate of a Cox, cause-specific Cox or Fine-Gray model
changes over follow-up (a violation of proportional hazards), compares several forms for how it
changes, refits the model with the chosen forms and reports the hazard ratio as a function of
time, in a table, in plain words and in figures. With \code{report()} the whole analysis goes into
one HTML file with the methods and references.

\strong{Use it} after \code{\link{ph_test}} has flagged a covariate, or as a routine check of a Cox
or Fine-Gray model whose hazard ratios you intend to report: it shows how the effect changes, not
only that it does.

\strong{It assists; it does not decide for you.} The screening tests, the comparison of the forms
and the final model are reported separately so that each step can be checked; you can fix the form
of any covariate with \code{fix}. Intervals of a model chosen from the data do not account for the
choice unless \code{bootstrap} is used.
}
\usage{
time_effects(fit, forms = c("constant", "log", "piecewise", "spline"),
    terms = NULL, fix = NULL, screen = c("holm", "unadjusted", "all"),
    alpha = 0.05, criterion = c("aic", "bic"),
    transform = c("km", "rank", "identity", "log"), df = 3L, breaks = NULL,
    times = NULL, level = 0.95, bootstrap = 0L, seed = 1L, max_rows = 5e6)
}
\arguments{
\item{fit}{a model fitted by \code{\link{cox}} (without \code{tt()} terms), by \code{\link{csc}} with
  the Cox baseline, or by \code{\link{fg}}. Strata, delayed entry, (start, stop] rows, case weights
  and clusters of the fit are used.}
\item{forms}{the forms of the log hazard ratio \eqn{\beta(t)} to compare: a character vector of
  built-in forms, or a list that can also contain named functions of time (your own forms, see
  Details). \code{"constant"} is always included.}
\item{terms}{the terms of the model to examine (as in the formula, for example \code{"age"} or
  \code{"education"}); default all. The other terms keep a constant effect.}
\item{fix}{a named list (or named character vector) \code{term = form} that sets the form of a
  term instead of choosing it, for example \code{fix = list(treatment = "piecewise")}. Use it to
  fit a pre-specified model or for a sensitivity analysis.}
\item{screen}{which terms are given a chance of a time-varying form: \code{"holm"} (default)
  those whose test of a constant effect has a Holm-adjusted p-value below \code{alpha};
  \code{"unadjusted"} those with p below \code{alpha}; \code{"all"} every term (the choice then
  rests on \code{criterion} alone).}
\item{alpha}{level of the screening tests.}
\item{criterion}{how the form is chosen among the candidates of a term: \code{"aic"} (default) or
  \code{"bic"} (penalty \eqn{\log} of the number of events; favours simpler forms).}
\item{transform}{time scale of the screening test, as in \code{\link{ph_test}}: \code{"km"}
  (default; one minus the Kaplan-Meier estimate, for a Fine-Gray model one minus the weighted
  product-limit estimate of the subdistribution), \code{"rank"}, \code{"identity"} or \code{"log"}.}
\item{df}{degrees of freedom of the spline form.}
\item{breaks}{times at which a step in the effect is allowed (form \code{"piecewise"}). Default
  \code{NULL}: one change point chosen from the data (see Details).}
\item{times}{times at which the hazard ratios are reported. Default: the 10th, 25th, 50th, 75th
  and 90th percentiles of the event times (two significant digits).}
\item{level}{confidence level.}
\item{bootstrap}{number of bootstrap samples of the subjects (clusters, if the fit has them) in
  which the whole procedure (screening, choice of the forms, final fit) is repeated. \code{0}
  (default): none. 200 or more gives percentile intervals that include the choice of the forms
  and the frequency with which each form is chosen.}
\item{seed}{seed of the bootstrap; the state of the random number generator of the session is
  restored afterwards.}
\item{max_rows}{largest number of rows of the data split at every event time. Beyond it, follow-up
  is split at quantiles of the event times instead, and the result says so.}
}
\details{
\subsection{Model}{
For a covariate column \eqn{x} the log hazard ratio (cause-specific hazard ratio for \code{csc},
subdistribution hazard ratio for \code{fg}) is \eqn{\beta(t) = \beta + \gamma^\top B(t)}, with
\eqn{B} one of the forms:
\describe{
  \item{\code{"constant"}}{\eqn{B = 0}: proportional hazards.}
  \item{\code{"log"}}{\eqn{B(t) = \log(t / t_m)}, \eqn{t_m} the median event time, so that
    \eqn{e^\beta} is the hazard ratio at \eqn{t_m}. The commonest pattern: an effect that fades
    or builds up gradually.}
  \item{\code{"linear"}}{\eqn{B(t) = t - t_m}; not in the default set.}
  \item{\code{"piecewise"}}{\eqn{B(t) = I(t > \tau)}: one hazard ratio before \eqn{\tau} and one
    after. Without \code{breaks}, \eqn{\tau} is chosen by profile likelihood among the 20th to 80th
    percentiles (in steps of 10) of the event times and counted as one extra parameter in the
    criterion and in the likelihood ratio test, which is then approximate. With \code{breaks}, one
    step at every break.}
  \item{\code{"spline"}}{a natural cubic spline in \eqn{\log t} with \code{df} degrees of freedom,
    knots at quantiles of the log event times (Hess 1994; Gray 1992): any smooth pattern.}
  \item{your own}{a named function of time returning a vector, or a matrix with one column per
    basis function, for example \code{forms = list("constant", "log", sqrt = function(t) sqrt(t))}.}
}
A factor gets one time-varying term per level (each with the same form). Every model is fitted by
partial likelihood on the data split at every event time, as \code{cox(..., tt = )} does
(Therneau and Grambsch 2000, Chapter 6), so the time functions are evaluated exactly at the event
times. For \code{fg} the same is done on the weighted counting-process data of Geskus (2011),
which reproduce the Fine-Gray estimator, with robust standard errors clustered on the subject
(Lin and Wei 1989).
}
\subsection{Procedure}{
For each analysis (the Cox model, every cause of a cause-specific model, or the Fine-Gray model):
\enumerate{
  \item \strong{Screening.} The score test of Grambsch and Therneau (1994) of a constant effect,
    for every term (all its columns together) and globally; the same test as \code{\link{ph_test}}.
    Holm's (1979) adjustment controls the chance of flagging any term when all effects are
    constant.
  \item \strong{Comparison.} For every term in \code{terms}, each form is fitted with the other
    terms constant, and compared with the constant form by AIC (or BIC) and by the likelihood
    ratio test. For terms that passed the screening (or are in \code{fix}), the form with the
    smallest criterion is chosen; for the others the effect stays constant. A term can pass the
    screening and still keep the constant form when no form improves the criterion.
  \item \strong{Final model.} All chosen forms are fitted together. If that model does not converge,
    the term whose form gains least is made constant and the model refitted (reported in a note).
  \item \strong{Effects.} \eqn{e^{\beta(t)}} with a Wald interval at \code{times}, and a curve over
    follow-up for the figures.
  \item \strong{Bootstrap} (optional). Steps 1 to 4 in every bootstrap sample.
}
}
\subsection{Combining with the rest of the package}{
\code{predict()} gives survival (\code{cox}) or cumulative incidence (\code{fg}; \code{csc}, from
the cause-specific hazards with product-limit overall survival) with the time-varying effects;
these predictions can be compared with other models through \code{\link{brier_score}} or
\code{\link{compare_models}} with a function that calls \code{time_effects()}.
\code{\link{sx_references}(x)} lists the methods used. A form found here can be refitted as an
ordinary model: \code{cox(Event(time, status) ~ x + tt(x) + z, d, tt = function(t) log(t))}.
}
}
\section{Interpreting the result}{
\describe{
  \item{Test of a constant effect}{a small p-value means that the hazard ratio of the term changes
    over time. With large data sets small, unimportant changes are significant: judge their size
    in the table of hazard ratios over time and in \code{plot(x)}.}
  \item{Comparison of forms}{\code{LR} and \code{p}: likelihood ratio test of the form against a
    constant effect. \code{dAIC}: difference from the smallest AIC of the term; differences below
    2 are small (Burnham and Anderson 2002). Forms within 2 of each other describe the data
    equally well; prefer the simpler or the one that suits the question (a step for an effect
    that is expected to stop, such as surgery).}
  \item{Hazard ratio over time}{the quantity to report. The reading says whether the effect weakens
    (moves towards 1), strengthens, or changes direction, and near which time.}
  \item{Stability}{with \code{bootstrap}: how often each form is chosen. A form chosen in, say,
    40\% of the samples is not well determined; report the hazard ratios at a few times (which are
    usually stable) rather than the form.}
  \item{Fine-Gray model}{the time-varying sHR describes the effect on the subdistribution hazard;
    it is not a ratio of risks. \code{predict()} gives the cumulative incidence.}
}
A sentence for a report: \emph{"The effect of treatment was not constant over time (Grambsch-Therneau
test, Holm-adjusted p = 0.002). A log-linear change in the log hazard ratio fitted best (AIC 6.1 lower
than a constant effect); the hazard ratio was 0.52 (95\% CI 0.38 to 0.71) at 6 months and 0.91 (0.67
to 1.24) at 3 years. The form was chosen in 81\% of 200 bootstrap samples."}
}
\section{Assumptions and how to check them}{
\describe{
  \item{The model apart from time}{the forms describe how effects change with time; the functional
    form of continuous covariates and the absence of interactions are still assumptions.}
  \item{Enough events}{every form spends parameters; with fewer than about 10 events per parameter
    the comparison is unreliable and the screening test has little power.}
  \item{Independent censoring}{as for the model that was fitted; for \code{fg} the censoring
    distribution has no covariates.}
}
}
\section{Messages and common problems}{
\describe{
  \item{\code{... was kept constant}}{the model with all chosen forms did not converge (often two
    flexible forms for correlated covariates); use \code{fix} or fewer forms.}
  \item{\code{follow-up was split at ... quantiles}}{the data split at every event time would have
    had more than \code{max_rows} rows; the forms are then evaluated at the quantiles, which is
    an approximation. Increase \code{max_rows} if memory allows.}
  \item{Long computing times}{each term and form is a separate fit on the split data, and the
    bootstrap repeats everything; examine fewer \code{terms} or use fewer \code{forms}.}
}
}
\section{Validation}{
The screening test is the test of \code{\link{ph_test}} and agrees with it to \eqn{10^{-8}}; the
fits of the log form agree with \code{cox(..., tt = )}, which agrees with
\code{survival::coxph(tt = )}; with every effect constant the Fine-Gray fit reproduces
\code{\link{fg}}, and the predictions reproduce those of \code{cox} and \code{fg}. These checks
are in the test suite of the package.
}
\value{
An S7 object of class \code{TimeEffects}. \code{x@units} has one element per analysis with
\code{screening} (one row per term: \code{chisq}, \code{df}, \code{p}, \code{p_holm},
\code{considered}), \code{comparison} (one row per term and form: \code{loglik}, \code{k} parameters,
\code{ic}, \code{delta}, \code{lr}, \code{lr_df}, \code{p_lr}, \code{converged}, \code{chosen}),
\code{selected} (the form of every term), \code{final} (coefficients, variance, likelihood ratio
test of all time-varying terms), \code{effects} (hazard ratios at \code{times}: \code{column},
\code{term}, \code{form}, \code{time}, \code{log_hr}, \code{se}, \code{hr}, \code{lower},
\code{upper}, \code{p}) and \code{notes}. \code{x@bootstrap}: selection frequencies and percentile
intervals.

Methods:
\tabular{ll}{
  \code{print(x)} \tab tests, the comparison for the terms with a time-varying form, hazard ratios
    over time and a reading \cr
  \code{summary(x)} \tab also the comparison for every screened term and the constant effects \cr
  \code{plot(x, term = NULL, cause = NULL, residuals = TRUE, scale = "ratio")} \tab \eqn{e^{\beta(t)}}
    with its band, the constant fit and the scaled Schoenfeld residuals \cr
  \code{predict(x, newdata, times, cause = 1)} \tab survival or cumulative incidence \cr
  \code{as.data.frame(x, what)} \tab \code{"effects"}, \code{"screening"}, \code{"forms"},
    \code{"coefficients"} or \code{"stability"} \cr
  \code{report(x, file)} \tab a self-contained HTML report \cr
  \code{\link{sx_references}(x)} \tab the references of the methods \cr
}
}
\references{
\emph{Test of proportional hazards and time-varying coefficients}

Grambsch PM, Therneau TM (1994) Proportional hazards tests and diagnostics based on weighted
residuals. Biometrika 81:515-526.

Therneau TM, Grambsch PM (2000) Modeling Survival Data: Extending the Cox Model. Springer,
Chapter 6.

Hess KR (1994) Assessing time-by-covariate interactions in proportional hazards regression models
using cubic spline functions. Stat Med 13:1045-1062.

Gray RJ (1992) Flexible methods for analyzing survival data using splines, with applications to
breast cancer prognosis. J Am Stat Assoc 87:942-951.

Buchholz A, Sauerbrei W (2011) Comparison of procedures to assess non-linear and time-varying
effects in multivariable models for survival data. Biom J 53:308-331.

\emph{Multiple testing and model choice}

Holm S (1979) A simple sequentially rejective multiple test procedure. Scand J Stat 6:65-70.

Akaike H (1974) A new look at the statistical model identification. IEEE Trans Automat Contr
19:716-723.

Schwarz G (1978) Estimating the dimension of a model. Ann Stat 6:461-464.

Burnham KP, Anderson DR (2002) Model Selection and Multimodel Inference, 2nd ed. Springer.

\emph{Fine-Gray model and variance}

Fine JP, Gray RJ (1999) A proportional hazards model for the subdistribution of a competing risk.
J Am Stat Assoc 94:496-509.

Geskus RB (2011) Cause-specific cumulative incidence estimation and the Fine and Gray model under
both left truncation and right censoring. Biometrics 67:39-49.

Lin DY, Wei LJ (1989) The robust inference for the Cox proportional hazards model. J Am Stat Assoc
84:1074-1078.

Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap. Chapman and Hall.

All references: \code{sx_references("time_effects")}.
}
\seealso{
\code{\link{ph_test}} (the test and residual plot), \code{\link{cox}} (\code{tt()} and
\code{strata()}), \code{\link{tt_effect}}, \code{\link{fg}}, \code{\link{csc}},
\code{\link{model_scan}}, \code{\link{rmst}} (a summary that needs no proportional hazards).
}
\examples{
set.seed(4)
n <- 400
x <- rbinom(n, 1, 0.5); z <- rnorm(n)
t <- rweibull(n, ifelse(x == 1, 0.7, 1.4), 3 * exp(-0.3 * z)); cens <- runif(n, 0, 6)
d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = z)

fit <- cox(Event(time, status) ~ x + z, d)
te <- time_effects(fit)
te                                   # x: time-varying, z: constant
plot(te)
as.data.frame(te)                    # hazard ratios over time
predict(te, data.frame(x = 0:1, z = 0), times = c(1, 2, 4))

## a pre-specified form, and your own function of time
time_effects(fit, fix = list(x = "piecewise"), breaks = 1)
time_effects(fit, forms = list("constant", "log", sqrt = function(t) sqrt(t)))

\donttest{
## stability of the choice and intervals that include it
te_b <- time_effects(fit, forms = c("constant", "log"), bootstrap = 50)
summary(te_b)
report(te_b, file = tempfile(fileext = ".html"))

## Fine-Gray model: time-varying subdistribution hazard ratios
data(fertility)
time_effects(fg(Event(time, status) ~ education, fertility, cause = 1))
}
}
