\name{Event}
\alias{Event}
\title{The event-time outcome of every SurvX model}
\description{
Builds the outcome (left-hand side of a model formula) from the follow-up time, the status and,
optionally, an entry time. Every model and estimator of SurvX uses it, for one event type and for
competing risks alike, for example \code{cox(Event(time, status) ~ age, data)} or
\code{aj(Event(time, status) ~ group, data)}.
}
\usage{
Event(time, status, entry = NULL, timefix = TRUE)
}
\arguments{
\item{time}{follow-up time: the time of the event, or the last time the subject was seen without
  the event (censoring). Non-negative numbers, in any unit (days, months, years); estimates such as
  hazards and medians are then in that unit. Measure it from a common origin (diagnosis,
  randomisation, surgery).}
\item{status}{what happened at \code{time}:
  \tabular{ll}{
    \code{0} \tab censored: the event had not happened when follow-up ended \cr
    \code{1} \tab the event (one event type) \cr
    \code{1, 2, ..., K} \tab competing risks: the type of event that happened first \cr
  }
  A logical (\code{TRUE} = event) is read as 1/0. A factor is read with its \emph{first} level as
  censoring and the following levels as causes 1, 2, ...; check the order with \code{levels()}.}
\item{entry}{optional entry time for delayed entry (left truncation): the subject is under
  observation, and in the risk sets, only from \code{entry} to \code{time}. Typical uses: age as the
  time scale with entry at the age of recruitment, or registries that include patients some time
  after diagnosis. It must be smaller than \code{time}. Also the start of (start, stop] rows for
  time-varying covariates in \code{\link{cox}} and \code{\link{csc}}. Supported by \code{cox},
  \code{csc}, \code{fg}, \code{parsurv}, \code{cure}, \code{hazard_np} and \code{dist_plot}; the
  other functions stop with a message.}
\item{timefix}{\code{TRUE} (default): times that differ only by floating-point rounding (relative
  difference below about \eqn{1.5 \times 10^{-8}}) are treated as tied, exactly as
  \pkg{survival} does in \code{coxph()} and \code{survfit()}. Such near-ties arise when times are
  computed (differences of dates divided by 365.25, for example). \code{FALSE} keeps the times as
  they are, as \pkg{cmprsk} and \pkg{prodlim} do.}
}
\details{
Censoring must be \emph{non-informative}: the reason follow-up ended must not depend on the risk of
the event (given the covariates in the model). With competing risks, a competing event is
\emph{not} censoring: code it with its own number, so that the methods for competing risks
(\code{\link{aj}}, \code{\link{fg}}, \code{\link{csc}}, \code{\link{directcif}}) treat it
correctly. Treating competing events as censored and using \code{\link{km}} overestimates the
probability of the event (Putter, Fiocco and Geskus 2007).

Rows with a missing value of \code{time} or \code{status} are removed by the model functions.
The time-fixing rule is that of \code{survival::aeqSurv}: the distinct values of \code{time} and
\code{entry} are sorted, and a value is merged into the previous one when the gap is at most
\eqn{\sqrt{\epsilon}} (\eqn{\epsilon} the machine precision), in absolute value or relative to the
mean absolute time.
}
\value{
A numeric matrix of class \code{"SurvXEvent"} with the columns \code{time}, \code{status} and,
with delayed entry, \code{entry}. It prints as \code{5}, \code{8+} (censored) or \code{3 (2)}
(cause 2), and as \code{(2, 8]} with an entry time.
}
\section{Messages}{
\describe{
  \item{\code{negative times}}{times must be zero or positive.}
  \item{\code{status must be 0, 1, 2, ...}}{non-integer or negative status codes.}
  \item{\code{entry must be smaller than time for every subject}}{a subject with zero or negative
    observation time; correct the data or drop the row.}
}
}
\references{
\emph{Censoring and left truncation:}
Klein JP, Moeschberger ML (2003). Survival Analysis: Techniques for Censored and Truncated Data,
  2nd ed. Springer, Sections 3.2-3.4.

\emph{Competing risks are not censoring:}
Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: competing risks and multi-state
  models. Statistics in Medicine 26:2389-2430.

\emph{Time fixing:}
Therneau TM. Roundoff error and tied times. Vignette of the R package \pkg{survival}
  (\code{vignette("tiedtimes", package = "survival")}).
}
\seealso{
\code{\link{cox}}, \code{\link{km}}, \code{\link{aj}}, \code{\link{fg}}, \code{\link{csc}},
\code{\link{parsurv}}.
}
\examples{
## one event type
Event(c(5, 8, 3), c(1, 0, 1))

## competing risks: 0 censored, 1 relapse, 2 death without relapse
Event(c(5, 8, 3, 6), c(1, 0, 2, 1))

## a factor: the first level is censoring
st <- factor(c("censored", "relapse", "death"), levels = c("censored", "relapse", "death"))
Event(c(4, 2, 7), st)

## delayed entry: observed from age at recruitment to age at event or censoring
Event(time = c(71, 68, 80), status = c(1, 0, 1), entry = c(65, 60, 70))

## times that differ only by rounding are tied by default
t <- c(1/3, 1 - 2/3, 0.5)
Event(t, c(1, 1, 0))
unclass(Event(t, c(1, 1, 0), timefix = FALSE))[, "time"] - 1/3
}
