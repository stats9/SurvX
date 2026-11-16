#' Event-time response
#'
#' @param time follow-up time (non-negative): the time of the event or of censoring.
#' @param status 0 = censored; 1 = event, or 1, 2, ..., K = cause of failure with competing risks.
#'   A logical is read as 0/1. A factor is read with its first level as censoring.
#' @param entry optional entry time for delayed entry (left truncation): a subject is under
#'   observation from \code{entry} to \code{time} and contributes to the risk sets only in that
#'   interval. Must be smaller than \code{time}. Supported by \code{\link{parsurv}} and by
#'   \code{\link{csc}} with a parametric baseline; the other functions stop with a message.
#' @return a matrix of class \code{"SurvXEvent"} with the columns \code{time}, \code{status} and, with
#'   delayed entry, \code{entry}, for the left-hand side of a model formula.
#' @references Klein JP, Moeschberger ML (2003) Survival Analysis: Techniques for Censored and
#'   Truncated Data, 2nd ed. Springer, Section 3.4 (left truncation).
#' @examples
#' Event(c(5, 8, 3), c(1, 0, 2))
#' Event(c(5, 8, 3), c(1, 0, 1), entry = c(0, 2, 1))
#' @export
Event <- function(time, status, entry = NULL) {
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
