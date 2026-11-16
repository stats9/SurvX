# Internal helpers.

.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}

.sx_rc <- function(rc, what) {
  if (rc == 0L) return(invisible(TRUE))
  msg <- switch(as.character(rc),
    "-1" = "invalid dimensions or a missing input",
    "-2" = "times must be finite and non-negative (strictly positive for baselines defined on log time)",
    "-3" = "no events of the requested type",
    paste("internal error code", rc))
  stop(what, ": ", msg, call. = FALSE)
}

# value at `times` of a right-continuous step function with jumps at x (increasing) to values y
.step_at <- function(x, y, times, left = 0) {
  if (!length(x)) return(rep(left, length(times)))
  c(left, y)[findInterval(times, x) + 1L]
}

# row-major copy of a matrix
.rm <- function(X) as.double(t(X))

# delayed entry is supported by the parametric models only; everything else stops with a clear message
.no_entry <- function(y, what) {
  if (ncol(unclass(y)) >= 3L)
    stop(what, " does not support delayed entry (left truncation) yet; use parsurv() or csc() with a ",
         "parametric baseline, or remove 'entry' from Event()", call. = FALSE)
}

.sx_frame <- function(formula, data, entry_ok = FALSE, what = "this function") {
  if (missing(data) || is.null(data)) data <- environment(formula)
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y <- mf[[1L]]
  if (!inherits(y, "Event")) stop("the left-hand side of the formula must be Event(time, status)", call. = FALSE)
  if (!entry_ok) .no_entry(y, what)
  tt <- attr(mf, "terms")
  X <- stats::model.matrix(tt, mf)
  ctr <- attr(X, "contrasts")
  asg <- attr(X, "assign")
  keep <- colnames(X) != "(Intercept)"
  X <- X[, keep, drop = FALSE]
  attr(X, "assign") <- asg[keep]                      # term of every column, used for covariate profiles
  y <- unclass(y)
  list(time = as.double(y[, 1L]), status = as.integer(round(y[, 2L])), X = X, terms = tt,
       xlevels = stats::.getXlevels(tt, mf), contrasts = ctr, mf = mf,
       entry = if (ncol(y) >= 3L) as.double(y[, 3L]) else NULL)
}

.sx_newX <- function(object, newdata) {
  if (is.null(newdata)) return(object@X)
  tt <- stats::delete.response(object@terms)
  mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass, xlev = object@xlevels)
  X <- stats::model.matrix(tt, mf, contrasts.arg = object@contrasts)
  X[, colnames(X) != "(Intercept)", drop = FALSE]
}

# one grouping factor from the right-hand side of a formula such as Event(time, status) ~ arm
.sx_groups <- function(formula, data, what = "this function", entry_ok = FALSE) {
  if (missing(data) || is.null(data)) data <- environment(formula)
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y <- mf[[1L]]
  if (!inherits(y, "Event")) stop("the left-hand side of the formula must be Event(time, status)", call. = FALSE)
  if (!entry_ok) .no_entry(y, what)
  y <- unclass(y)
  g <- if (ncol(mf) > 1L) interaction(mf[-1L], drop = TRUE, sep = ", ") else factor(rep("All patients", nrow(mf)))
  list(time = as.double(y[, 1L]), status = as.integer(round(y[, 2L])), group = droplevels(as.factor(g)),
       entry = if (ncol(y) >= 3L) as.double(y[, 3L]) else NULL)
}

.matrix_rm <- function(v, nrow, ncol, dn = NULL) {
  matrix(v, nrow = nrow, ncol = ncol, byrow = TRUE, dimnames = dn)
}

.flag_text <- function(flags) {
  bits <- c("1" = "iteration limit reached", "2" = "information matrix not positive definite",
            "4" = "non-finite log-likelihood", "8" = "a parameter is diverging or drifts on a flat likelihood",
            "16" = "sum of fitted CIFs exceeds 1 within follow-up", "32" = "sum of fitted CIFs exceeds 1 at infinity",
            "64" = "constrained fit stopped on the boundary")
  on <- bitwAnd(as.integer(flags), as.integer(names(bits))) > 0L
  unname(bits[on])
}


# confidence limits of a cumulative incidence F with standard error se on the complementary log-log scale:
# F^exp(-/+ z se / (F log F)); NA where F is 0 or 1
.cif_interval <- function(Fv, se, zq) {
  ok <- is.finite(Fv) & is.finite(se) & Fv > 0 & Fv < 1
  e <- suppressWarnings(zq * se / (Fv * log(Fv)))
  lower <- suppressWarnings(Fv^exp(-e)); upper <- suppressWarnings(Fv^exp(e))
  lower[!ok] <- NA_real_; upper[!ok] <- NA_real_
  list(lower = lower, upper = upper)
}
