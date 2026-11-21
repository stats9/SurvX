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
    "-6" = "every entry (start) time must be smaller than its time (stop)",
    paste("internal error code", rc))
  stop(what, ": ", msg, call. = FALSE)
}

# value at `times` of a right-continuous step function with jumps at x (increasing) to values y
.step_at <- function(x, y, times, left = 0) {
  if (!length(x)) return(rep(left, length(times)))
  c(left, y)[findInterval(times, x) + 1L]
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# row-major copy of a matrix
.rm <- function(X) as.double(t(X))

# delayed entry is supported by the parametric models only; everything else stops with a clear message
.no_entry <- function(y, what) {
  if (ncol(unclass(y)) >= 3L)
    stop(what, " does not support delayed entry (left truncation) yet; use cox(), csc() or parsurv(), ",
         "or remove 'entry' from Event()", call. = FALSE)
}

.sx_frame <- function(formula, data, entry_ok = FALSE, what = "this function", specials_ok = FALSE, tt_ok = FALSE) {
  if (missing(data) || is.null(data)) data <- environment(formula)
  # strata(), cluster() and tt() are read from the formula by SurvX itself (see .sx_specials_env), so
  # they work without survival and are never mistaken for ordinary covariates
  tt0 <- stats::terms(formula, specials = c("strata", "cluster", "tt"), data = data)
  sp <- attr(tt0, "specials")
  has_sp <- !is.null(sp$strata) || !is.null(sp$cluster)
  if (has_sp && !specials_ok)
    stop("strata() and cluster() are available in cox() and csc() with a Cox baseline, not in ", what, call. = FALSE)
  if (!is.null(sp$tt) && !tt_ok)
    stop("tt() (time-varying coefficients) is available in cox(), not in ", what, call. = FALSE)
  environment(tt0) <- .sx_specials_env(environment(formula))
  mf <- stats::model.frame(tt0, data = data, na.action = stats::na.omit)
  y <- mf[[1L]]
  if (!inherits(y, "SurvXEvent")) stop("the left-hand side of the formula must be Event(time, status)", call. = FALSE)
  if (!entry_ok) .no_entry(y, what)
  tt <- attr(mf, "terms")
  strata <- NULL; strata_levels <- NULL; strata_vars <- NULL; cluster <- NULL; cluster_var <- NULL
  if (has_sp) {
    fac <- attr(tt, "factors")
    vars <- rownames(fac)
    drop_terms <- integer(0)
    for (nm in c("strata", "cluster")) {
      ix <- sp[[nm]]
      if (is.null(ix)) next
      for (i in ix) {
        tm <- which(fac[i, ] > 0)
        if (any(colSums(fac[, tm, drop = FALSE] > 0) > 1L))
          stop(nm, "() cannot be part of an interaction", call. = FALSE)
        drop_terms <- c(drop_terms, tm)
      }
      if (nm == "strata") {
        sv <- lapply(vars[ix], function(v) mf[[v]])
        g <- if (length(sv) == 1L) sv[[1L]] else interaction(sv, drop = TRUE, sep = ", ")
        g <- droplevels(as.factor(g))
        strata_levels <- levels(g); strata <- as.integer(g) - 1L; strata_vars <- vars[ix]
      } else {
        if (length(ix) > 1L) stop("only one cluster() term is allowed", call. = FALSE)
        cluster <- mf[[vars[ix]]]; cluster_var <- vars[ix]
      }
    }
    tt <- if (length(unique(drop_terms)) == ncol(fac)) stats::terms(stats::update(formula, . ~ 1))
          else stats::drop.terms(tt, unique(drop_terms), keep.response = TRUE)
    environment(tt) <- environment(tt0)
  }
  X <- stats::model.matrix(tt, mf)
  ctr <- attr(X, "contrasts")
  asg <- attr(X, "assign")
  keep <- colnames(X) != "(Intercept)"
  X <- X[, keep, drop = FALSE]
  attr(X, "assign") <- asg[keep]                      # term of every column, used for covariate profiles
  tt_cols <- integer(0)
  if (!is.null(sp$tt)) {                              # columns of the tt() terms: x * g(t), built by cox()
    lab <- attr(tt, "term.labels")
    tt_terms <- which(grepl("^tt\\(", lab))
    if (any(grepl("tt\\(", lab) & grepl(":", lab, fixed = TRUE)))
      stop("tt() cannot be part of an interaction; write tt(x) on its own", call. = FALSE)
    tt_cols <- which(attr(X, "assign") %in% tt_terms)
  }
  y <- unclass(y)
  list(time = as.double(y[, 1L]), status = as.integer(round(y[, 2L])), X = X, terms = tt,
       xlevels = stats::.getXlevels(tt, mf), contrasts = ctr, mf = mf,
       entry = if (ncol(y) >= 3L) as.double(y[, 3L]) else NULL,
       strata = strata, strata_levels = strata_levels, strata_vars = strata_vars,
       cluster = cluster, cluster_var = cluster_var, tt_cols = tt_cols)
}

# environment in which strata() and cluster() mean what SurvX needs: strata() combines its
# arguments into one factor, cluster() returns its argument unchanged
.sx_specials_env <- function(parent) {
  e <- new.env(parent = parent)
  e$strata <- function(..., sep = ", ") {
    a <- list(...)
    if (length(a) == 1L) as.factor(a[[1L]]) else interaction(a, drop = TRUE, sep = sep)
  }
  e$cluster <- function(x) x
  e$tt <- function(x) {
    if (!is.numeric(x) && !is.logical(x))
      stop("tt() needs a numeric covariate; for a factor, create 0/1 indicator variables", call. = FALSE)
    as.numeric(x)
  }
  e
}

# stratum codes (0-based) of new data; NULL for unstratified fits
.sx_newstrata <- function(object, newdata) {
  if (is.null(object@strata)) return(NULL)
  if (is.null(newdata)) return(object@strata)
  env <- .sx_specials_env(environment(object@terms))
  sv <- lapply(object@strata_vars, function(v) eval(str2lang(v), newdata, env))
  g <- if (length(sv) == 1L) as.character(sv[[1L]]) else as.character(interaction(sv, sep = ", "))
  code <- match(g, object@strata_levels) - 1L
  if (anyNA(code)) stop("newdata has strata that are not in the data: ", paste(unique(g[is.na(code)]), collapse = ", "), call. = FALSE)
  code
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
  if (!inherits(y, "SurvXEvent")) stop("the left-hand side of the formula must be Event(time, status)", call. = FALSE)
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

# methods that need the covariates at a single time point do not apply to models with tt() terms
.no_tt <- function(object, what) {
  if (S7::S7_inherits(object, CoxFit) && !is.null(object@tt))
    stop(what, " is not available for a model with time-varying coefficients (tt()); ",
         "use estimates(fit) and tt_effect(fit, times) for the effect at chosen times", call. = FALSE)
  invisible(TRUE)
}

# evaluate expr with the given seed and restore the state of the random number generator afterwards,
# so that calling a SurvX function never changes the random numbers the user draws next
.with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
            if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

.rng_save <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
}
.rng_restore <- function(state) {
  if (is.null(state)) { if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()) }
  else assign(".Random.seed", state, envir = globalenv())
}
