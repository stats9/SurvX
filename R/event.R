#' The event-time outcome of every SurvX model
#'
#' @description
#' Builds the outcome (left-hand side of a model formula) from the follow-up time, the status and,
#' optionally, an entry time. Every model and estimator of SurvX uses it, for one event type and for
#' competing risks alike, for example \code{cox(Event(time, status) ~ age, data)} or
#' \code{aj(Event(time, status) ~ group, data)}.
#' @param time follow-up time: the time of the event, or the last time the subject was seen without
#'   the event (censoring). Non-negative numbers, in any unit (days, months, years); estimates such as
#'   hazards and medians are then in that unit. Measure it from a common origin (diagnosis,
#'   randomisation, surgery).
#' @param status what happened at \code{time}:
#'   \tabular{ll}{
#'   \code{0} \tab censored: the event had not happened when follow-up ended \cr
#'   \code{1} \tab the event (one event type) \cr
#'   \code{1, 2, ..., K} \tab competing risks: the type of event that happened first \cr
#'   }
#'   A logical (\code{TRUE} = event) is read as 1/0. A factor is read with its \emph{first} level as
#'   censoring and the following levels as causes 1, 2, ...; check the order with \code{levels()}.
#' @param entry optional entry time for delayed entry (left truncation): the subject is under
#'   observation, and in the risk sets, only from \code{entry} to \code{time}. Typical uses: age as the
#'   time scale with entry at the age of recruitment, or registries that include patients some time
#'   after diagnosis. It must be smaller than \code{time}. Also the start of (start, stop] rows for
#'   time-varying covariates in \code{\link{cox}} and \code{\link{csc}}. Supported by \code{cox},
#'   \code{csc}, \code{fg}, \code{parsurv}, \code{cure}, \code{hazard_np} and \code{dist_plot}; the
#'   other functions stop with a message.
#' @param timefix \code{TRUE} (default): times that differ only by floating-point rounding (relative
#'   difference below about \eqn{1.5 \times 10^{-8}}) are treated as tied, exactly as
#'   \pkg{survival} does in \code{coxph()} and \code{survfit()}. Such near-ties arise when times are
#'   computed (differences of dates divided by 365.25, for example). \code{FALSE} keeps the times as
#'   they are, as \pkg{cmprsk} and \pkg{prodlim} do.
#' @details
#' Censoring must be \emph{non-informative}: the reason follow-up ended must not depend on the risk of
#' the event (given the covariates in the model). With competing risks, a competing event is
#' \emph{not} censoring: code it with its own number, so that the methods for competing risks
#' (\code{\link{aj}}, \code{\link{fg}}, \code{\link{csc}}, \code{\link{directcif}}) treat it
#' correctly. Treating competing events as censored and using \code{\link{km}} overestimates the
#' probability of the event (Putter, Fiocco and Geskus 2007).
#'
#' Rows with a missing value of \code{time} or \code{status} are removed by the model functions.
#' The time-fixing rule is that of \code{survival::aeqSurv}: the distinct values of \code{time} and
#' \code{entry} are sorted, and a value is merged into the previous one when the gap is at most
#' \eqn{\sqrt{\epsilon}} (\eqn{\epsilon} the machine precision), in absolute value or relative to the
#' mean absolute time.
#' @section Messages:
#' \describe{
#'   \item{\code{negative times}}{times must be zero or positive.}
#'   \item{\code{status must be 0, 1, 2, ...}}{non-integer or negative status codes.}
#'   \item{\code{entry must be smaller than time for every subject}}{a subject with zero or negative
#'     observation time; correct the data or drop the row.}
#' }
#' @return
#' A numeric matrix of class \code{"SurvXEvent"} with the columns \code{time}, \code{status} and,
#' with delayed entry, \code{entry}. It prints as \code{5}, \code{8+} (censored) or \code{3 (2)}
#' (cause 2), and as \code{(2, 8]} with an entry time.
#' @references
#' \emph{Censoring and left truncation:}
#' Klein JP, Moeschberger ML (2003). Survival Analysis: Techniques for Censored and Truncated Data,
#'   2nd ed. Springer, Sections 3.2-3.4.
#'
#' \emph{Competing risks are not censoring:}
#' Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: competing risks and multi-state
#'   models. Statistics in Medicine 26:2389-2430.
#'
#' \emph{Time fixing:}
#' Therneau TM. Roundoff error and tied times. Vignette of the R package \pkg{survival}
#'   (\code{vignette("tiedtimes", package = "survival")}).
#' @seealso
#' \code{\link{cox}}, \code{\link{km}}, \code{\link{aj}}, \code{\link{fg}}, \code{\link{csc}},
#' \code{\link{parsurv}}.
#' @examples
#' ## one event type
#' Event(c(5, 8, 3), c(1, 0, 1))
#'
#' ## competing risks: 0 censored, 1 relapse, 2 death without relapse
#' Event(c(5, 8, 3, 6), c(1, 0, 2, 1))
#'
#' ## a factor: the first level is censoring
#' st <- factor(c("censored", "relapse", "death"), levels = c("censored", "relapse", "death"))
#' Event(c(4, 2, 7), st)
#'
#' ## delayed entry: observed from age at recruitment to age at event or censoring
#' Event(time = c(71, 68, 80), status = c(1, 0, 1), entry = c(65, 60, 70))
#'
#' ## times that differ only by rounding are tied by default
#' t <- c(1/3, 1 - 2/3, 0.5)
#' Event(t, c(1, 1, 0))
#' unclass(Event(t, c(1, 1, 0), timefix = FALSE))[, "time"] - 1/3
#' @export
Event <- function(time, status, entry = NULL, timefix = TRUE) {
  if (missing(status)) status <- rep(1L, length(time))
  if (is.factor(status)) status <- as.integer(status) - 1L
  if (is.logical(status)) status <- as.integer(status)
  if (!is.numeric(time) || !is.numeric(status)) stop("time and status must be numeric", call. = FALSE)
  if (length(time) != length(status)) stop("time and status differ in length", call. = FALSE)
  ok <- !is.na(time) & !is.na(status)
  if (any(time[ok] < 0)) stop("negative times", call. = FALSE)
  if (any(status[ok] < 0 | status[ok] != round(status[ok]))) stop("status must be 0, 1, 2, ...", call. = FALSE)
  y <- cbind(time = as.numeric(time), status = as.numeric(status))
  if (!is.null(entry)) {
    if (!is.numeric(entry) || length(entry) != length(time)) stop("entry must be numeric and as long as time", call. = FALSE)
    ok <- ok & !is.na(entry)
    if (any(entry[ok] < 0)) stop("negative entry times", call. = FALSE)
    if (any(entry[ok] >= time[ok])) stop("entry must be smaller than time for every subject", call. = FALSE)
    y <- cbind(y, entry = as.numeric(entry))
  }
  if (isTRUE(timefix)) y <- .timefix(y)
  class(y) <- "SurvXEvent"
  y
}

#' @export
"[.SurvXEvent" <- function(x, i, j, drop = FALSE) {
  if (missing(j)) {
    y <- unclass(x)[i, , drop = FALSE]
    class(y) <- "SurvXEvent"
    y
  } else {
    unclass(x)[i, j, drop = drop]
  }
}

#' @export
is.na.SurvXEvent <- function(x) {
  as.vector(rowSums(is.na(unclass(x))) > 0)
}

#' @export
print.SurvXEvent <- function(x, ...) {
  y <- unclass(x)
  s <- ifelse(y[, 2L] == 0, "+", ifelse(y[, 2L] == 1 & max(y[, 2L], na.rm = TRUE) <= 1, "", paste0(" (", y[, 2L], ")")))
  tt <- format(y[, 1L])
  if (ncol(y) >= 3L) tt <- paste0("(", format(y[, 3L]), ", ", trimws(tt), "]")
  print(paste0(tt, s), quote = FALSE, ...)
  invisible(x)
}


# Times that differ only by floating-point rounding are made equal, as survival::aeqSurv() does for
# coxph() and survfit() (timefix = TRUE): the distinct values of time and entry are sorted, and a
# value is merged into the previous one when the gap is at most sqrt(.Machine$double.eps) in
# absolute value or relative to the mean absolute time (vignette "Roundoff error and tied times" of the
# survival package).
.timefix <- function(y, tol = sqrt(.Machine$double.eps)) {
  m <- unclass(y)
  cols <- intersect(c("time", "entry"), colnames(m))
  v <- m[, cols]
  u <- sort(unique(v[is.finite(v)]))
  if (length(u) < 2L) return(y)
  du <- diff(u)
  tied <- du <= tol | du / mean(abs(u)) <= tol
  if (!any(tied)) return(y)
  cuts <- u[c(TRUE, !tied)]
  for (cn in cols) {
    ok <- is.finite(m[, cn])
    m[ok, cn] <- cuts[findInterval(m[ok, cn], cuts)]
  }
  if ("entry" %in% cols) {
    ok <- is.finite(m[, "entry"]) & is.finite(m[, "time"])
    if (any(m[ok, "entry"] >= m[ok, "time"]))
      stop("entry and time differ only by rounding error for some subject (an interval of length 0); ",
           "correct the data or use Event(..., timefix = FALSE)", call. = FALSE)
  }
  class(m) <- class(y)
  m
}

# the same rule for a plain vector of times (functions that take time and status vectors)
.timefix_vec <- function(time) {
  y <- cbind(time = as.numeric(time))
  as.vector(unclass(.timefix(y))[, "time"])
}

#' Interval-censored event times
#'
#' @description
#' Builds the outcome of a model formula when the event time is not seen exactly but is known to lie
#' in an interval: the subject was event-free at the visit at time \code{left} and had had the event
#' by the visit at time \code{right}. Typical sources: periodic examinations (tumour recurrence or
#' progression seen at scheduled scans, HIV or HPV infection found at a test, dental caries, onset of
#' a chronic condition between two surveys). Use it on the left-hand side of
#' \code{\link{turnbull}} (non-parametric), \code{\link{parsurv}} (parametric regression) and
#' \code{\link{directcif}} (competing risks), for example
#' \code{parsurv(Interval(last_neg, first_pos) ~ age, data)}.
#' @param left the last time the subject was known to be event-free. \code{0} or \code{NA}: the
#'   event was already there at the first visit (left censoring, the event lies in
#'   \code{(0, right]}).
#' @param right the first time the event was known to have happened. \code{NA} or \code{Inf}: no
#'   event was seen (right censoring at \code{left}). \code{right = left}: the event time is known
#'   exactly.
#' @param cause competing risks only: the type of event (1, 2, ...) of the subjects with an event
#'   (finite \code{right}); ignored for right-censored subjects. A factor is read with its levels as
#'   causes 1, 2, ... (no censoring level: censoring is given by \code{right}).
#' @details
#' Every subject contributes the probability of the interval in which the event is known to lie:
#' \tabular{lll}{
#'   \code{left < right < Inf} \tab interval-censored \tab \eqn{S(L) - S(R)} \cr
#'   \code{left = 0, right < Inf} \tab left-censored \tab \eqn{1 - S(R)} \cr
#'   \code{right = Inf} \tab right-censored \tab \eqn{S(L)} \cr
#'   \code{left = right} \tab exact \tab density \eqn{f(L)} \cr
#' }
#' (with competing risks \eqn{F_k(R) - F_k(L)}, \eqn{1 - \sum_k F_k(L)} and \eqn{f_k(L)}). Imputing
#' the midpoint or the right end point of the interval and then treating the times as exact biases the
#' estimates and understates their uncertainty (Lindsey and Ryan 1998; Law and Brookmeyer 1992). The
#' inspection process must be non-informative: the visit times may differ between subjects but must
#' not depend on the event itself (given the covariates). Rows with \code{left = 0} and
#' \code{right = Inf} carry no information and are dropped by the model functions.
#' @section Messages:
#' \describe{
#'   \item{\code{right must be at least left}}{an interval with its end before its start.}
#'   \item{\code{cause must be 1, 2, ...}}{a subject with an event and no (or a zero) cause.}
#' }
#' @return
#' A numeric matrix of class \code{"SurvXInterval"} with the columns \code{left}, \code{right}
#' (\code{Inf} for right censoring) and \code{status} (0 right-censored, otherwise the cause). It
#' prints as \code{(2, 5]}, \code{3} (exact), \code{4+} (right-censored) and \code{(0, 5]}
#' (left-censored).
#' @references
#' Turnbull BW (1976). The empirical distribution function with arbitrarily grouped, censored and
#'   truncated data. Journal of the Royal Statistical Society B 38:290-295.
#'
#' Lindsey JC, Ryan LM (1998). Tutorial in biostatistics: methods for interval-censored data.
#'   Statistics in Medicine 17:219-238.
#'
#' Law CG, Brookmeyer R (1992). Effects of mid-point imputation on the analysis of doubly censored
#'   data. Statistics in Medicine 11:1569-1578.
#'
#' Sun J (2006). The Statistical Analysis of Interval-censored Failure Time Data. Springer.
#' @seealso
#' \code{\link{turnbull}}, \code{\link{parsurv}}, \code{\link{directcif}}, \code{\link{Event}}.
#' @examples
#' ## event between visits, exact, right-censored, left-censored
#' Interval(left = c(2, 3, 4, 0), right = c(5, 3, NA, 6))
#'
#' ## competing risks: the type of event is known
#' Interval(c(2, 1, 4), c(5, 3, Inf), cause = c(1, 2, 0))
#' @export
Interval <- function(left, right, cause = NULL) {
  if (!is.numeric(left) || !is.numeric(right)) stop("left and right must be numeric", call. = FALSE)
  if (length(left) != length(right)) stop("left and right differ in length", call. = FALSE)
  L <- as.numeric(left); R <- as.numeric(right)
  L[is.na(L)] <- 0; R[is.na(R)] <- Inf
  if (any(L < 0) || any(R < 0)) stop("negative times", call. = FALSE)
  if (any(!is.finite(L))) stop("left must be finite (use right = Inf for right censoring)", call. = FALSE)
  if (any(R < L)) stop("right must be at least left for every subject", call. = FALSE)
  if (any(R == 0)) stop("an event at time 0 (right = 0): times must be positive", call. = FALSE)
  ev <- is.finite(R)
  st <- as.numeric(ev)
  if (!is.null(cause)) {
    if (length(cause) != length(L)) stop("cause must be as long as left", call. = FALSE)
    if (is.factor(cause)) cause <- as.integer(cause)
    if (is.logical(cause)) cause <- as.integer(cause)
    if (!is.numeric(cause)) stop("cause must be numeric or a factor", call. = FALSE)
    ce <- cause[ev]
    if (any(is.na(ce) | ce < 1 | ce != round(ce))) stop("cause must be 1, 2, ... for every subject with an event (finite right)", call. = FALSE)
    st[ev] <- ce
  }
  y <- cbind(left = L, right = R, status = st)
  class(y) <- "SurvXInterval"
  y
}

#' @export
"[.SurvXInterval" <- function(x, i, j, drop = FALSE) {
  if (missing(j)) {
    y <- unclass(x)[i, , drop = FALSE]
    class(y) <- "SurvXInterval"
    y
  } else {
    unclass(x)[i, j, drop = drop]
  }
}

#' @export
is.na.SurvXInterval <- function(x) {
  as.vector(rowSums(is.na(unclass(x))) > 0)
}

#' @export
format.SurvXInterval <- function(x, ...) {
  y <- unclass(x)
  L <- y[, 1L]; R <- y[, 2L]; s <- y[, 3L]
  fl <- format(L, trim = TRUE, ...); fr <- format(R, trim = TRUE, ...)
  out <- ifelse(!is.finite(R), paste0(fl, "+"), ifelse(L == R, fl, paste0("(", fl, ", ", fr, "]")))
  if (max(s, 0) > 1) out <- ifelse(s > 0, paste0(out, " (", s, ")"), out)
  out
}

#' @export
print.SurvXInterval <- function(x, ...) {
  print(format(x), quote = FALSE)
  invisible(x)
}

#' @export
as.character.SurvXInterval <- function(x, ...) format(x)

# the observation type of interval-censored rows: 0 right-censored, 1 exact, 2 interval, 3 left-censored
.ic_type <- function(L, R) ifelse(!is.finite(R), 0L, ifelse(L == R, 1L, ifelse(L == 0, 3L, 2L)))

# a single time per subject for starting values and knots: the midpoint of the interval
.ic_mid <- function(L, R) ifelse(!is.finite(R), L, ifelse(L == 0, R / 2, (L + R) / 2))

.is_ic <- function(fit) {
  v <- tryCatch(S7::prop(fit, "interval"), error = function(e) NULL)
  !is.null(v)
}

.no_ic <- function(fit, what, instead = NULL) {
  if (.is_ic(fit))
    stop(what, " is not available for interval-censored data (Interval())",
         if (!is.null(instead)) paste0("; ", instead), call. = FALSE)
  invisible(TRUE)
}
