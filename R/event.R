#' Event-time response
#'
#' @param time follow-up time (non-negative).
#' @param status 0 = censored; 1 = event, or 1, 2, ..., K = cause of failure with competing risks.
#'   A logical is read as 0/1. A factor is read with its first level as censoring.
#' @return a two-column matrix of class \code{"Event"} for the left-hand side of a model formula.
#' @export
Event <- function(time, status) {
  if (missing(status)) status <- rep(1L, length(time))
  if (is.factor(status)) status <- as.integer(status) - 1L
  if (is.logical(status)) status <- as.integer(status)
  if (!is.numeric(time) || !is.numeric(status)) stop("time and status must be numeric", call. = FALSE)
  if (length(time) != length(status)) stop("time and status differ in length", call. = FALSE)
  ok <- !is.na(time) & !is.na(status)
  if (any(time[ok] < 0)) stop("negative times", call. = FALSE)
  if (any(status[ok] < 0 | status[ok] != round(status[ok]))) stop("status must be 0, 1, 2, ...", call. = FALSE)
  y <- cbind(time = as.numeric(time), status = as.numeric(status))
  class(y) <- "Event"
  y
}

#' @export
"[.Event" <- function(x, i, j, drop = FALSE) {
  if (missing(j)) {
    y <- unclass(x)[i, , drop = FALSE]
    class(y) <- "Event"
    y
  } else {
    unclass(x)[i, j, drop = drop]
  }
}

#' @export
is.na.Event <- function(x) {
  as.vector(rowSums(is.na(unclass(x))) > 0)
}

#' @export
print.Event <- function(x, ...) {
  y <- unclass(x)
  s <- ifelse(y[, 2L] == 0, "+", ifelse(y[, 2L] == 1 & max(y[, 2L], na.rm = TRUE) <= 1, "", paste0(" (", y[, 2L], ")")))
  print(paste0(format(y[, 1L]), s), quote = FALSE, ...)
  invisible(x)
}
