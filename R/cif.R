# Cumulative incidence curves of any model, for one or several models, every cause, chosen or
# typical covariate values, the average over the sample or no covariates; drawn with base graphics
# or ggplot2. cif_data() computes the curves, cif_plot() draws them.

utils::globalVariables(c("time", "cif", "lower", "upper", "curve", "cause_label", "step"))

.okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442")

# ----------------------------------------------------------------------------------------------
# covariate profiles
# ----------------------------------------------------------------------------------------------

# typical covariate values: the mean of numeric covariates and the most frequent level of categorical
# covariates (treatment contrasts); other terms (interactions, splines, other contrasts) at the mean of
# their columns. Returns the design row and a description.
.profile <- function(object) {
  X <- object@X; p <- ncol(X)
  if (!p) return(list(X = matrix(0, 1L, 0L), label = "no covariates"))
  asg <- attr(X, "assign"); tt <- object@terms
  row <- colMeans(X)
  if (is.null(asg)) return(list(X = matrix(row, 1L, dimnames = list(NULL, colnames(X))), label = "mean of every covariate"))
  labs <- attr(tt, "term.labels"); dc <- attr(tt, "dataClasses"); ctr <- object@contrasts
  parts <- character(0); other <- FALSE
  for (j in unique(asg)) {
    cols <- which(asg == j); lab <- labs[j]
    cls <- if (!is.null(dc) && lab %in% names(dc)) dc[[lab]] else NA_character_
    treat <- is.null(ctr[[lab]]) || identical(ctr[[lab]], "contr.treatment")
    if (!is.na(cls) && cls %in% c("factor", "character", "logical", "ordered") && treat) {
      fr <- colMeans(X[, cols, drop = FALSE]); ref_fr <- 1 - sum(fr)
      levs <- object@xlevels[[lab]]
      ref <- if (length(levs)) levs[1L] else if (cls == "logical") "FALSE" else "reference"
      row[cols] <- 0
      if (max(fr) > ref_fr) {
        w <- cols[which.max(fr)]; row[w] <- 1
        lev <- substring(colnames(X)[w], nchar(lab) + 1L)
      } else lev <- ref
      parts <- c(parts, paste0(lab, " = ", lev))
    } else if (length(cols) == 1L && !is.na(cls) && cls == "numeric") {
      parts <- c(parts, paste0(lab, " = ", format(signif(row[cols], 3L), trim = TRUE)))
    } else other <- TRUE
  }
  if (other) parts <- c(parts, "other terms at their mean")
  list(X = matrix(row, 1L, dimnames = list(NULL, colnames(X))), label = paste(parts, collapse = ", "))
}

# labels of the rows of newdata: row names when set, else the values of the variables that differ
.newdata_labels <- function(object, newdata) {
  rn <- rownames(newdata)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(newdata))))) return(rn)
  vars <- intersect(all.vars(stats::delete.response(object@terms)), names(newdata))
  if (!length(vars)) return(paste("profile", seq_len(nrow(newdata))))
  vary <- vars[vapply(vars, function(v) length(unique(newdata[[v]])) > 1L, logical(1))]
  if (!length(vary)) vary <- vars
  apply(newdata[, vary, drop = FALSE], 1L, function(r) paste0(vary, " = ", trimws(format(r)), collapse = ", "))
}

# ----------------------------------------------------------------------------------------------
# refits: without covariates, or the weighted direct model for another cause
# ----------------------------------------------------------------------------------------------
.link_arg <- function(fit, k) {
  if (fit@alpha_free[k]) return("GOR")
  a <- fit@alpha_value[k]
  if (a == 0) "PH" else if (a == 1) "PO" else a
}

.frame_of <- function(fit, without_covariates) {
  list(time = fit@time, status = fit@status,
       X = if (without_covariates) matrix(0, fit@n, 0L) else fit@X,
       terms = if (without_covariates) stats::terms(stats::as.formula("~ 1")) else fit@terms,
       xlevels = if (without_covariates) list() else fit@xlevels,
       contrasts = if (without_covariates) NULL else fit@contrasts)
}

.refit_direct <- function(fit, without_covariates = FALSE, cause = NULL) {
  f <- .frame_of(fit, without_covariates)
  est <- fit@estimation
  if (est == "joint") {
    kn <- if (any(fit@baseline == "spline")) fit@knots else NULL
    links <- lapply(seq_len(fit@K), function(k) .link_arg(fit, k))
    return(.directcif_core(fit@call, f, fit@baseline, links, "joint", 1L, "logistic", 3L, kn, fit@constraint, 0L, NULL, 500L, 1e-10))
  }
  k <- if (is.null(cause)) fit@interest else cause
  kn <- if (fit@baseline[1L] == "spline" && k == fit@interest) fit@knots[[1L]] else NULL
  .directcif_core(fit@call, f, fit@baseline[1L], .link_arg(fit, 1L), est, as.integer(k),
                  if (is.na(fit@cond)) "logistic" else fit@cond, 3L, kn, FALSE, 0L, NULL, 500L, 1e-10)
}

# ----------------------------------------------------------------------------------------------
# curves of one model
# ----------------------------------------------------------------------------------------------
.grid_smooth <- function(tmax) seq(0, tmax, length.out = 201L)

.grid_step <- function(time, status, tmax) {
  g <- sort(unique(c(0, time[status > 0L], tmax)))
  if (length(g) > 600L) g <- unique(c(0, stats::quantile(g, seq(0, 1, length.out = 600L), type = 1L, names = FALSE), tmax))
  g
}

.interval <- function(Fv, se, lo_p, up_p, level, how) {
  if (how == "percentile" && !is.null(lo_p)) return(list(lower = lo_p, upper = up_p))
  .cif_interval(Fv, se, .zq(level))
}

.curve_frame <- function(model, k, profile, time, Fv, se, ci, step, method) {
  data.frame(model = model, cause = as.integer(k), profile = profile, time = time, cif = as.numeric(Fv),
             se = as.numeric(se), lower = as.numeric(ci$lower), upper = as.numeric(ci$upper), step = step,
             method = method, stringsAsFactors = FALSE)
}

# rows: list(X = design rows, labels = labels) for profiles or newdata; NULL for type = "marginal"
.curves_direct <- function(fit, name, causes, rows, times, conf.int, level) {
  out <- list()
  for (k in causes) {
    f <- fit; kk <- k
    if (f@estimation == "weighted" && k != f@interest) {
      f <- tryCatch(.refit_direct(fit, ncol(fit@X) == 0L, cause = k), error = function(e) NULL)
      if (is.null(f)) { warning(name, ": the weighted model could not be fitted for cause ", k, call. = FALSE); next }
    }
    if (f@estimation == "joint" && k > f@K) next
    usevar <- conf.int && all(is.finite(f@vcov))
    if (is.null(rows)) {                               # average over the subjects
      Xn <- f@X; nn <- nrow(Xn); nt <- length(times); p <- ncol(Xn)
      v <- if (usevar) as.double(f@vcov) else double(1)
      r <- if (f@estimation == "constrained")
        .C("SX_constrained_predict_mean", as.integer(p), f@alpha_free, as.double(f@alpha_value), f@base_code,
           as.integer(.cond_codes[[f@cond]]), as.double(f@coefficients), as.integer(usevar), v,
           if (kk == f@interest) 1L else 2L, as.integer(nn), .rm(Xn), as.integer(nt), as.double(times),
           cif = double(nt), se = double(nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
      else
        .C("SX_direct_predict_mean", as.integer(p), if (f@estimation == "joint") f@K else 1L, f@alpha_free, as.double(f@alpha_value),
           f@base_code, f@n_knots, as.double(unlist(f@knots)), as.double(f@coefficients), as.integer(usevar), v,
           if (f@estimation == "joint") as.integer(kk) else 1L, as.integer(nn), .rm(Xn), as.integer(nt), as.double(times),
           cif = double(nt), se = double(nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
      .sx_rc(r$rc, "cif_data")
      se <- if (usevar) r$se else rep(NA_real_, nt)
      out[[length(out) + 1L]] <- .curve_frame(name, k, "average over the sample", times, r$cif, se, .cif_interval(r$cif, se, .zq(level)), FALSE, "delta method")
    } else {
      for (i in seq_len(nrow(rows$X))) {
        pr <- .direct_predict(f, rows$X[i, , drop = FALSE], times, kk, TRUE)
        se <- if (usevar) drop(pr$se) else rep(NA_real_, length(times))
        Fv <- drop(pr$cif)
        out[[length(out) + 1L]] <- .curve_frame(name, k, rows$labels[i], times, Fv, se, .cif_interval(Fv, se, .zq(level)), FALSE, "delta method")
      }
    }
  }
  do.call(rbind, out)
}

.curves_boot <- function(fit, name, causes, rows, times, conf.int, level, B, seed, interval, adjust) {
  model <- if (S7::S7_inherits(fit, FineGrayFit)) 0L else 1L
  ties <- if (model == 1L && fit@ties == "breslow") 0L else 1L
  pl <- if (model == 1L) as.integer(fit@product_limit) else 1L
  # without covariates the cause-specific Cox model reduces to Nelson-Aalen hazards; with Breslow
  # increments and the product-limit survival the CIF is then exactly the Aalen-Johansen estimate
  if (model == 1L && adjust == 0L) { ties <- 0L; pl <- 1L }
  K <- max(fit@status); n <- fit@n; p <- ncol(fit@X)
  Xn <- if (is.null(rows)) fit@X else rows$X
  labels <- if (is.null(rows)) "average over the sample" else rows$labels
  avg <- as.integer(is.null(rows)); nn <- nrow(Xn); nt <- length(times); cells <- (if (avg) 1L else nn) * nt
  Bu <- if (conf.int) as.integer(B) else 0L
  out <- list()
  for (k in causes) {
    if (k > K) next
    r <- .C("SX_cif_bootstrap", model, as.integer(n), as.integer(p), fit@time, fit@status, .rm(fit@X), as.integer(K), as.integer(k),
            ties, pl, as.integer(adjust), as.integer(nn), .rm(Xn), avg, as.integer(nt), as.double(times), Bu, as.integer(seed), as.double(level),
            est = double(cells), se = double(cells), lower = double(cells), upper = double(cells), n_ok = integer(1), rc = integer(1),
            NAOK = TRUE, PACKAGE = "SurvX")
    if (r$rc == -4L) { warning(name, ": the model could not be fitted for cause ", k, call. = FALSE); next }
    .sx_rc(r$rc, "cif_data")
    if (Bu > 0L && r$n_ok < Bu) warning(name, ", cause ", k, ": ", Bu - r$n_ok, " of ", Bu, " bootstrap samples could not be fitted", call. = FALSE)
    est <- .matrix_rm(r$est, cells / nt, nt)
    se <- if (Bu > 1L) .matrix_rm(r$se, cells / nt, nt) else matrix(NA_real_, cells / nt, nt)
    lo <- .matrix_rm(r$lower, cells / nt, nt); up <- .matrix_rm(r$upper, cells / nt, nt)
    for (i in seq_len(nrow(est))) {
      ci <- if (Bu > 1L) .interval(est[i, ], se[i, ], lo[i, ], up[i, ], level, interval) else list(lower = rep(NA_real_, nt), upper = rep(NA_real_, nt))
      out[[length(out) + 1L]] <- .curve_frame(name, k, labels[i], times, est[i, ], se[i, ], ci, TRUE,
                                              if (Bu > 1L) paste0("bootstrap (", r$n_ok, " samples)") else "none")
    }
  }
  do.call(rbind, out)
}

.aj_frame <- function(time, status, level) {
  n <- length(time); K <- max(status)
  r <- .C("SX_aj", as.integer(n), as.double(time), as.integer(status), as.integer(K), time = double(n), surv = double(n),
          cif = double(n * K), se = double(n * K), m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "Aalen-Johansen")
  m <- r$m
  list(time = r$time[seq_len(m)], cif = .matrix_rm(r$cif[seq_len(m * K)], m, K), se = .matrix_rm(r$se[seq_len(m * K)], m, K), K = K)
}

.curves_aj <- function(curves, name, causes, times, level, K) {
  out <- list()
  for (g in names(curves)) for (k in causes) {
    if (k > K) next
    cv <- curves[[g]]
    Fv <- .step_at(cv$time, cv[[paste0("cif.", k)]], times, 0); se <- .step_at(cv$time, cv[[paste0("se.", k)]], times, 0)
    out[[length(out) + 1L]] <- .curve_frame(name, k, g, times, Fv, se, .cif_interval(Fv, se, .zq(level)), TRUE, "delta method")
  }
  do.call(rbind, out)
}

# ----------------------------------------------------------------------------------------------
# cif_data()
# ----------------------------------------------------------------------------------------------

#' Cumulative incidence curves of fitted models
#'
#' Computes the cumulative incidence function of every cause for one or several models, at chosen
#' covariate values, at typical covariate values, as the average over the sample, or from the same
#' models fitted without covariates. \code{cif_plot()} draws the result; \code{cif_data()} returns it
#' as a data frame for tables or other graphics.
#'
#' @details
#' \describe{
#'   \item{\code{type = "profile"} (default)}{one curve at typical covariate values: the mean of every
#'     numeric covariate and the most frequent level of every categorical covariate.}
#'   \item{\code{newdata}}{one curve for every row; \code{type} is then ignored.}
#'   \item{\code{type = "marginal"}}{the average of the predicted curves of all subjects.}
#'   \item{\code{type = "unadjusted"}}{the same model fitted without covariates (for a direct model the
#'     parametric curve of Jeong and Fine 2006; for Fine-Gray and cause-specific Cox the models with all
#'     coefficients zero).}
#' }
#' Curves of all causes are computed: a Fine-Gray model is refitted for every other cause and a
#' weighted direct model is refitted for every other cause, with the same covariates.
#'
#' Pointwise confidence intervals: delta method for direct models and the Aalen-Johansen estimator,
#' nonparametric bootstrap of the subjects for Fine-Gray and cause-specific Cox models (Efron and
#' Tibshirani 1993); both on the complementary log-log scale, \eqn{F^{\exp(\pm z\,se/(F \log F))}},
#' so that the limits stay in [0, 1]. \code{interval = "percentile"} gives bootstrap percentile limits instead.
#'
#' @param x a model fitted by \code{\link{directcif}}, \code{\link{fg}} or \code{\link{csc}}, an
#'   \code{\link{aj}} estimate, or a named list of these for a comparison.
#' @param newdata data frame of covariate values, one curve per row.
#' @param cause causes to show; default all.
#' @param type \code{"profile"}, \code{"marginal"} or \code{"unadjusted"}; see Details.
#' @param times time points; default a fine grid (direct models) or the event times (step functions).
#' @param conf.int compute pointwise confidence intervals.
#' @param level confidence level.
#' @param B number of bootstrap samples for Fine-Gray and cause-specific Cox models.
#' @param seed seed of the bootstrap (the same seed gives the same intervals).
#' @param interval \code{"cloglog"} or \code{"percentile"} for bootstrap intervals.
#' @return a data frame with the columns \code{model, cause, profile, time, cif, se, lower, upper, step,
#'   method}; attribute \code{"description"} says which covariate values were used.
#' @references
#' Jeong J-H, Fine JP (2006) J R Stat Soc C 55:187-200.
#' Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap. Chapman and Hall, New York.
#' Kalbfleisch JD, Prentice RL (2002) The Statistical Analysis of Failure Time Data, 2nd ed. Wiley,
#'   Section 8.2 (cumulative incidence functions).
#' @examples
#' set.seed(6)
#' n <- 300
#' x <- rbinom(n, 1, 0.5)
#' t1 <- rexp(n, 0.4 * exp(0.6 * x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
#' d <- data.frame(time = pmin(t1, t2, cens), x = x)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#' fit <- directcif(Event(time, status) ~ x, d)
#' head(cif_data(fit))
#' cif_plot(fit)
#' cif_plot(fit, newdata = data.frame(x = 0:1), cause = 1, risk.table = FALSE)
#' cif_plot(list("direct PH" = fit, "Fine-Gray" = fg(Event(time, status) ~ x, d)), B = 50)
#' cif_plot(fit, layout = "stacked", conf.int = FALSE)
#' @export
cif_data <- function(x, newdata = NULL, cause = NULL, type = c("profile", "marginal", "unadjusted"), times = NULL,
                     conf.int = TRUE, level = 0.95, B = 200L, seed = 1L, interval = c("cloglog", "percentile")) {
  type <- match.arg(type); interval <- match.arg(interval)
  models <- if (is.list(x) && !S7::S7_inherits(x)) x else list(x)
  if (is.null(names(models))) names(models) <- if (length(models) == 1L) .model_label(models[[1L]]) else paste("model", seq_along(models))
  bad <- names(models) == ""
  names(models)[bad] <- vapply(models[bad], .model_label, character(1))
  ref <- models[[1L]]
  K <- max(ref@status)
  causes <- if (is.null(cause)) seq_len(K) else as.integer(cause)
  tmax <- max(ref@time)
  if (!is.null(newdata)) type <- "newdata"
  out <- list(); desc <- character(0)
  for (nm in names(models)) {
    m <- models[[nm]]
    if (S7::S7_inherits(m, AJFit)) {
      tg <- if (is.null(times)) .grid_step(m@time, m@status, tmax) else times
      out[[nm]] <- .curves_aj(m@curves, nm, causes, tg, level, m@K)
      next
    }
    is_direct <- S7::S7_inherits(m, DirectCIFFit)
    if (!is_direct && !S7::S7_inherits(m, FineGrayFit) && !S7::S7_inherits(m, CSCFit))
      stop("cif_data() takes models fitted by directcif(), fg(), csc() or an aj() estimate", call. = FALSE)
    tg <- if (!is.null(times)) times else if (is_direct) .grid_smooth(tmax) else .grid_step(m@time, m@status, tmax)
    rows <- NULL; adjust <- 1L; fit <- m
    if (type == "newdata") {
      rows <- list(X = .sx_newX(m, newdata), labels = .newdata_labels(m, newdata))
    } else if (type == "profile") {
      pf <- .profile(m); rows <- list(X = pf$X, labels = "typical covariate values"); desc[nm] <- pf$label
    } else if (type == "unadjusted") {
      if (is_direct) {
        fit <- tryCatch(.refit_direct(m, TRUE), error = function(e) stop(nm, ": the model could not be fitted without covariates: ", conditionMessage(e), call. = FALSE))
        rows <- list(X = matrix(0, 1L, 0L), labels = "without covariates")
      } else {
        rows <- list(X = matrix(0, 1L, ncol(m@X)), labels = "without covariates"); adjust <- 0L
      }
    }
    out[[nm]] <- if (is_direct) .curves_direct(fit, nm, causes, rows, tg, conf.int, level)
                 else .curves_boot(fit, nm, causes, rows, tg, conf.int, level, B, seed, interval, adjust)
  }
  d <- do.call(rbind, out)
  rownames(d) <- NULL
  attr(d, "description") <- if (type == "profile" && !length(desc)) NULL else switch(type,
    profile = if (length(unique(desc)) <= 1L) paste0("Typical covariate values: ", desc[1L], " (mean of numeric, most frequent level of categorical covariates)")
              else paste0("Typical covariate values (mean of numeric, most frequent level of categorical covariates): ", paste0(names(desc), ": ", desc, collapse = "; ")),
    marginal = "Average of the predicted curves of all subjects",
    unadjusted = "Models fitted without covariates",
    newdata = "Covariate values given in newdata")
  attr(d, "level") <- level
  attr(d, "risk_time") <- ref@time
  attr(d, "groups") <- if (S7::S7_inherits(ref, AJFit)) ref@group else NULL
  d
}

.model_label <- function(m) {
  if (S7::S7_inherits(m, AJFit)) return("Aalen-Johansen")
  if (S7::S7_inherits(m, FineGrayFit)) return("Fine-Gray")
  if (S7::S7_inherits(m, CSCFit)) return("cause-specific Cox")
  if (S7::S7_inherits(m, DirectCIFFit)) return(paste0("direct ", m@baseline[1L], " ", m@link[1L], if (m@estimation != "joint") paste0(", ", m@estimation)))
  "model"
}

# ----------------------------------------------------------------------------------------------
# cif_plot()
# ----------------------------------------------------------------------------------------------

#' @rdname cif_data
#' @param layout \code{"panels"} (one panel per cause), \code{"overlay"} (all causes in one panel) or
#'   \code{"stacked"} (the cumulative incidences of the causes stacked on each other, with the
#'   event-free probability on top; one panel per model and covariate profile).
#' @param risk.table show the number of subjects at risk under the time axis.
#' @param legend show the legend.
#' @param caption show which covariate values were used.
#' @param reference add the Aalen-Johansen estimate of the whole sample for comparison.
#' @param mark.times time points at which the values of the curves are written on the plot.
#' @param percent label the probability axis in percent.
#' @param grid draw light horizontal grid lines.
#' @param labels labels of the curves (model and profile combinations), in the order of the legend.
#' @param cause.names names of the causes, e.g. \code{c("relapse", "death")}.
#' @param col,lty,lwd colours, line types and width of the curves.
#' @param xlab,ylab,main,xlim,ylim axis labels, title and limits.
#' @param cex size of all text, relative to the device font (default 0.85). With
#'   \code{engine = "ggplot2"} the theme's base size is \code{12 * cex} points.
#' @param engine \code{"base"} (default, no dependencies) or \code{"ggplot2"} (returns a ggplot object;
#'   the table of numbers at risk is not drawn).
#' @param ... passed to \code{plot.default} (base engine).
#' @export
cif_plot <- function(x, newdata = NULL, cause = NULL, type = c("profile", "marginal", "unadjusted"),
                     layout = c("panels", "overlay", "stacked"), times = NULL, conf.int = TRUE, level = 0.95,
                     B = 200L, seed = 1L, interval = c("cloglog", "percentile"), risk.table = TRUE, legend = TRUE,
                     caption = TRUE, reference = FALSE, mark.times = NULL, percent = FALSE, grid = FALSE,
                     labels = NULL, cause.names = NULL, col = NULL, lty = NULL, lwd = 2, xlab = "Time", ylab = NULL,
                     main = NULL, xlim = NULL, ylim = NULL, cex = 0.85, engine = c("base", "ggplot2"), ...) {
  type <- match.arg(type); layout <- match.arg(layout); engine <- match.arg(engine); interval <- match.arg(interval)
  if (layout == "stacked") conf.int <- FALSE
  d <- cif_data(x, newdata = newdata, cause = cause, type = type, times = times, conf.int = conf.int, level = level,
                B = B, seed = seed, interval = interval)
  if (isTRUE(reference)) {
    ref <- if (is.list(x) && !S7::S7_inherits(x)) x[[1L]] else x
    a <- .aj_frame(ref@time, ref@status, level)
    cv <- data.frame(time = a$time, a$cif, a$se); names(cv) <- c("time", paste0("cif.", seq_len(a$K)), paste0("se.", seq_len(a$K)))
    tg <- .grid_step(ref@time, ref@status, max(ref@time))
    ra <- .curves_aj(list("all subjects" = cv), "Aalen-Johansen", unique(d$cause), tg, level, a$K)
    if (!conf.int) ra$lower <- ra$upper <- NA_real_
    keep <- list(description = attr(d, "description"), level = attr(d, "level"), risk_time = attr(d, "risk_time"), groups = attr(d, "groups"))
    d <- rbind(d, ra)
    for (a_ in names(keep)) if (!is.null(keep[[a_]])) attr(d, a_) <- keep[[a_]]
  }
  if (!nrow(d)) stop("nothing to plot", call. = FALSE)
  # one key per curve: model and profile (only what distinguishes the curves)
  # (the Aalen-Johansen reference curve does not count as a profile, so that a comparison of models at one
  # profile is labelled by the model names only)
  nmod <- length(unique(d$model)); nprof <- length(unique(d$profile[d$model != "Aalen-Johansen"]))
  d$curve <- if (nmod > 1L && nprof > 1L) paste0(d$model, ": ", d$profile) else if (nmod > 1L) d$model else d$profile
  keys <- unique(d$curve)
  if (!is.null(labels)) { if (length(labels) != length(keys)) stop("labels must have ", length(keys), " elements", call. = FALSE); d$curve <- labels[match(d$curve, keys)]; keys <- labels }
  causes <- sort(unique(d$cause))
  cn <- if (is.null(cause.names)) paste("Cause", causes) else rep_len(as.character(cause.names), max(causes))[causes]
  d$cause_label <- factor(cn[match(d$cause, causes)], levels = cn)
  cap <- if (isTRUE(caption)) attr(d, "description") else NULL
  if (engine == "ggplot2") return(.cif_ggplot(d, keys, layout, conf.int, percent, col, xlab, ylab, main, cap, xlim, ylim, lwd, cex))
  .cif_base(d, keys, causes, cn, layout, conf.int, risk.table, legend, mark.times, percent, grid, col, lty, lwd,
            xlab, ylab, main, cap, xlim, ylim, cex = cex, ...)
  invisible(d)
}

.step_poly <- function(t, lo, hi) {
  ok <- is.finite(lo) & is.finite(hi)
  if (sum(ok) < 2L) return(NULL)
  t <- t[ok]; lo <- lo[ok]; hi <- hi[ok]; n <- length(t)
  xs <- c(rep(t, each = 2L)[-1L], t[n]); yl <- rep(lo, each = 2L); yh <- rep(hi, each = 2L)
  list(x = c(xs, rev(xs)), y = c(yl, rev(yh)))
}

.prob_axis <- function(percent) {
  at <- graphics::axTicks(2L)
  graphics::axis(2L, at = at, labels = if (percent) paste0(round(100 * at), "%") else format(at), las = 1L)
}

.cif_base <- function(d, keys, causes, cn, layout, conf.int, risk.table, legend, mark.times, percent, grid, col, lty, lwd,
                      xlab, ylab, main, cap, xlim, ylim, cex = 0.85, ...) {
  nk <- length(keys)
  if (is.null(col)) col <- rep_len(.okabe_ito, nk)
  col <- rep_len(col, nk)
  lty <- rep_len(if (is.null(lty)) 1L else lty, nk)
  rt <- attr(d, "risk_time")
  groups <- attr(d, "groups")
  nrow_rt <- if (isTRUE(risk.table)) (if (!is.null(groups) && nlevels(groups) > 1L) nlevels(groups) else 1L) else 0L
  panels <- if (layout == "panels") lapply(causes, function(k) d[d$cause == k, , drop = FALSE])
            else if (layout == "overlay") list(d)
            else lapply(keys, function(k) d[d$curve == k, , drop = FALSE])
  titles <- if (layout == "panels") cn else if (layout == "overlay") "" else keys
  np <- length(panels)
  nc <- if (np <= 3L) np else ceiling(sqrt(np)); nr <- ceiling(np / nc)
  # long titles and captions are wrapped to the width of the device, so that no text runs off the page
  win <- grDevices::dev.size("in")[1L] - 0.4
  main_l <- if (!is.null(main)) .wrap_inches(main, win, 1.1 * cex, font = 2L) else character(0)
  cap_l <- if (!is.null(cap)) .wrap_inches(cap, win, 0.8 * cex) else character(0)
  top <- (if (length(main_l)) 0.4 + 1.25 * length(main_l) else 0) + (if (length(cap_l)) 0.3 + 1.0 * length(cap_l) else 0)
  op <- graphics::par(mfrow = c(nr, nc), mar = c(4.2 + if (nrow_rt) 1.6 + 1.1 * nrow_rt else 0, 4.5, 2.2, 1),
                      oma = c(0, 0, top, 0), mgp = c(2.6, 0.7, 0))
  op2 <- graphics::par(cex = cex)            # after mfrow, which resets cex; all text and margins scale with it
  on.exit(graphics::par(c(op2, op)))
  xr <- if (is.null(xlim)) c(0, max(d$time)) else xlim
  for (j in seq_len(np)) {
    s <- panels[[j]]
    if (layout == "stacked") {
      .stacked_panel(s, causes, cn, xr, percent, grid, xlab, ylab, titles[j], legend && j == 1L, ...)
    } else {
      top <- max(c(s$cif, if (conf.int) s$upper), na.rm = TRUE)
      yr <- if (is.null(ylim)) c(0, min(1, max(0.05, 1.06 * top))) else ylim
      yl <- if (!is.null(ylab)) ylab else "Cumulative incidence"
      graphics::plot(NA, xlim = xr, ylim = yr, xlab = xlab, ylab = yl, bty = "l", yaxt = "n", main = titles[j], font.main = 1, ...)
      .prob_axis(percent)
      if (isTRUE(grid)) graphics::abline(h = graphics::axTicks(2L), col = "grey92", lwd = 0.8)
      overlay_cause <- layout == "overlay" && length(causes) > 1L
      for (q in seq_len(nk)) for (k in causes) {
        c_ <- s[s$curve == keys[q] & s$cause == k, , drop = FALSE]
        if (!nrow(c_)) next
        cc <- if (overlay_cause && nk == 1L) .okabe_ito[match(k, causes)] else col[q]
        ll <- if (overlay_cause && nk > 1L) match(k, causes) else lty[q]
        if (conf.int && any(is.finite(c_$lower))) {
          pg <- if (c_$step[1L]) .step_poly(c_$time, c_$lower, c_$upper)
                else { ok <- is.finite(c_$lower) & is.finite(c_$upper); if (sum(ok) > 1L) list(x = c(c_$time[ok], rev(c_$time[ok])), y = c(c_$lower[ok], rev(c_$upper[ok]))) }
          if (!is.null(pg)) graphics::polygon(pg$x, pg$y, col = grDevices::adjustcolor(cc, alpha.f = 0.15), border = NA)
        }
        graphics::lines(c_$time, c_$cif, type = if (c_$step[1L]) "s" else "l", col = cc, lty = ll, lwd = lwd)
        if (length(mark.times)) {
          v <- vapply(mark.times, function(t) { i <- max(which(c_$time <= t), 0L); if (i) c_$cif[i] else 0 }, numeric(1))
          graphics::points(mark.times, v, pch = 19, cex = 0.6, col = cc)
          # with several curves the labels are stacked (one text line per curve) so that they do not overlap
          graphics::text(mark.times, v, labels = if (percent) sprintf("%.1f%%", 100 * v) else sprintf("%.3f", v), pos = 3, cex = 0.7,
                         col = cc, offset = 0.35 + 0.85 * (q - 1L))
        }
      }
      if (length(mark.times)) graphics::abline(v = mark.times, lty = 3, col = "grey65")
      if (isTRUE(legend) && j == 1L) {
        if (overlay_cause && nk == 1L) graphics::legend("topleft", legend = cn, col = .okabe_ito[seq_along(causes)], lty = 1, lwd = lwd, bty = "n", cex = 0.85)
        else if (overlay_cause) graphics::legend("topleft", legend = c(keys, cn), col = c(col, rep("grey30", length(causes))),
                                                 lty = c(rep(1, nk), seq_along(causes)), lwd = lwd, bty = "n", cex = 0.85)
        else if (nk > 1L) graphics::legend("topleft", legend = keys, col = col, lty = lty, lwd = lwd, bty = "n", cex = 0.85)
      }
    }
    if (nrow_rt) .risk_table(rt, groups, xr)
  }
  .outer_lines(main_l, cap_l, cex)
}

# width in inches of each string at character expansion cex, measured with the font metrics of a
# throw-away pdf device at the point size of the current device (no plot is needed on the current one)
.text_inches <- function(txt, cex = 1, ps = NULL, font = 1L) {
  if (!length(txt)) return(numeric(0))
  if (is.null(ps)) ps <- graphics::par("ps")
  cur <- grDevices::dev.cur()
  grDevices::pdf(NULL, pointsize = ps); on.exit({ grDevices::dev.off(); if (cur > 1L) grDevices::dev.set(cur) })
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  graphics::strwidth(txt, units = "inches", cex = cex, font = font)
}

# wrap a text so that every line fits in width_in inches at character expansion cex
.wrap_inches <- function(txt, width_in, cex, ps = NULL, font = 1L) {
  out <- character(0)
  for (z in txt) {
    n <- nchar(z); lines <- z
    while (n > 15L && max(.text_inches(lines, cex, ps, font)) > 0.95 * width_in) { n <- floor(n * 0.9); lines <- strwrap(z, width = n) }
    out <- c(out, lines)
  }
  out
}

# write the wrapped title and caption lines in the outer top margin, the caption below the title
.outer_lines <- function(main_l, cap_l, cex = 1) {
  # mtext() in the outer margin does not scale its text by par("cex"), so the size is given in full
  # here, and it is the same size at which .wrap_inches() measured the lines
  ln <- 0.3
  for (i in rev(seq_along(cap_l))) { graphics::mtext(cap_l[i], side = 3, outer = TRUE, line = ln, cex = 0.8 * cex, col = "grey30"); ln <- ln + 1.0 }
  if (length(cap_l)) ln <- ln + 0.2
  for (i in rev(seq_along(main_l))) { graphics::mtext(main_l[i], side = 3, outer = TRUE, line = ln, font = 2, cex = 1.1 * cex); ln <- ln + 1.25 }
}

.risk_table <- function(rt, groups, xr) {
  at <- graphics::axTicks(1L); at <- at[at >= xr[1L] & at <= xr[2L]]
  usr <- graphics::par("usr"); left <- usr[1L] - 0.02 * diff(usr[1:2])
  graphics::mtext("Number at risk", side = 1, line = 3.9, at = usr[1L], adj = 0, cex = 0.72, font = 2)
  if (!is.null(groups) && nlevels(groups) > 1L) {
    for (g in seq_len(nlevels(groups))) {
      lev <- levels(groups)[g]; tt <- rt[groups == lev]
      graphics::mtext(vapply(at, function(a) sum(tt >= a), numeric(1)), side = 1, line = 3.9 + 1.1 * g, at = at, cex = 0.72)
      graphics::mtext(lev, side = 1, line = 3.9 + 1.1 * g, at = left, adj = 1, cex = 0.72)
    }
  } else graphics::mtext(vapply(at, function(a) sum(rt >= a), numeric(1)), side = 1, line = 5.0, at = at, cex = 0.72)
}

.stacked_panel <- function(s, causes, cn, xr, percent, grid, xlab, ylab, title, legend, ...) {
  tg <- sort(unique(s$time))
  M <- sapply(causes, function(k) { c_ <- s[s$cause == k, , drop = FALSE]; .step_at(c_$time, c_$cif, tg, 0) })
  if (is.null(dim(M))) M <- matrix(M, ncol = length(causes))
  cum <- t(apply(M, 1L, cumsum)); if (length(causes) == 1L) cum <- t(cum)
  over <- which(cum[, ncol(cum)] > 1 + 1e-9)
  cum <- pmin(cum, 1)
  graphics::plot(NA, xlim = xr, ylim = c(0, 1), xlab = xlab, ylab = if (is.null(ylab)) "Probability" else ylab, bty = "l",
                 yaxt = "n", main = title, font.main = 1, ...)
  .prob_axis(percent)
  stepped <- isTRUE(s$step[1L])
  sx <- function(v) if (stepped) list(x = c(rep(tg, each = 2L)[-1L], tg[length(tg)]), y = rep(v, each = 2L)) else list(x = tg, y = v)
  low <- rep(0, length(tg)); fills <- .okabe_ito[seq_along(causes)]
  for (j in seq_along(causes)) {
    a <- sx(cum[, j]); b <- sx(low)
    graphics::polygon(c(a$x, rev(b$x)), c(a$y, rev(b$y)), col = grDevices::adjustcolor(fills[j], alpha.f = 0.55), border = fills[j])
    low <- cum[, j]
  }
  if (isTRUE(grid)) graphics::abline(h = graphics::axTicks(2L), col = "grey92", lwd = 0.8)
  # the note goes inside the plotting region, top right, so that it never collides with the title
  if (length(over)) graphics::mtext(sprintf("sum of the risks > 1 from t = %s", format(signif(tg[over[1L]], 3L))),
                                    side = 3, line = -1.1, adj = 0.98, cex = 0.7, col = "#D55E00")
  if (isTRUE(legend)) graphics::legend("topleft", legend = c(rev(cn), "no event yet"), fill = c(rev(grDevices::adjustcolor(fills, alpha.f = 0.55)), "white"),
                                       border = c(rev(fills), "grey50"), bty = "n", cex = 0.85)
}

.cif_ggplot <- function(d, keys, layout, conf.int, percent, col, xlab, ylab, main, cap, xlim, ylim, lwd, cex = 0.85) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("engine = \"ggplot2\" needs the ggplot2 package", call. = FALSE)
  d$curve <- factor(d$curve, levels = keys)
  pal <- if (is.null(col)) rep_len(.okabe_ito, length(keys)) else rep_len(col, length(keys))
  yl <- if (is.null(ylab)) (if (layout == "stacked") "Probability" else "Cumulative incidence") else ylab
  if (layout == "stacked") {
    # the same stacking as the base engine: cause 1 at the bottom, the running sum capped at 1
    st <- do.call(rbind, lapply(levels(d$curve), function(k) {
      s <- d[d$curve == k, , drop = FALSE]; cs <- levels(s$cause_label)[levels(s$cause_label) %in% s$cause_label]
      tg <- sort(unique(s$time))
      M <- sapply(cs, function(cl) { c_ <- s[s$cause_label == cl, , drop = FALSE]; .step_at(c_$time, c_$cif, tg, 0) })
      M <- matrix(M, nrow = length(tg)); cum <- pmin(t(apply(M, 1L, cumsum)), 1)
      if (length(cs) == 1L) cum <- matrix(cum, ncol = 1L)
      if (isTRUE(s$step[1L]) && length(tg) > 1L) {          # step curves: draw the stairs, not straight lines
        n <- length(tg); r <- c(rep(seq_len(n - 1L), each = 2L), n)
        tg <- c(tg[1L], rep(tg[-1L], each = 2L)); cum <- cum[r, , drop = FALSE]; M <- M[r, , drop = FALSE]
      }
      low <- cbind(0, cum[, -ncol(cum), drop = FALSE])
      data.frame(curve = k, cause_label = rep(cs, each = length(tg)), time = rep(tg, length(cs)),
                 ymin = as.vector(low), ymax = as.vector(cum), over = rep(rowSums(M) > 1 + 1e-9, length(cs)))
    }))
    st$curve <- factor(st$curve, levels = levels(d$curve)); st$cause_label <- factor(st$cause_label, levels = levels(d$cause_label))
    p <- ggplot2::ggplot(st, ggplot2::aes(x = time, ymin = ymin, ymax = ymax, fill = cause_label)) +
      ggplot2::geom_ribbon(alpha = 0.7, colour = NA) +
      ggplot2::facet_wrap(~ curve) + ggplot2::scale_fill_manual(values = .okabe_ito, name = NULL) +
      ggplot2::coord_cartesian(ylim = c(0, 1))
    if (any(st$over)) {
      t1 <- tapply(st$time[st$over], st$curve[st$over], min)
      note <- paste0("Sum of the risks of all causes > 1 from t = ", paste(format(signif(t1[is.finite(t1)], 3L)), collapse = ", "),
                     " (shown capped at 1).")
      cap <- c(cap, note)
    }
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = time, y = cif, colour = curve, fill = curve, group = interaction(curve, cause_label)))
    if (conf.int) {
      p <- p + ggplot2::geom_ribbon(data = d[is.finite(d$lower) & !d$step, , drop = FALSE], ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA)
      ds <- d[is.finite(d$lower) & d$step, , drop = FALSE]
      if (nrow(ds)) p <- p + ggplot2::geom_rect(data = transform(ds, xend = stats::ave(ds$time, ds$curve, ds$cause, FUN = function(v) c(v[-1L], v[length(v)]))),
                                                ggplot2::aes(xmin = time, xmax = xend, ymin = lower, ymax = upper), alpha = 0.15, colour = NA)
    }
    p <- p + ggplot2::geom_line(data = d[!d$step, , drop = FALSE], linewidth = lwd / 3) +
      ggplot2::geom_step(data = d[d$step, , drop = FALSE], linewidth = lwd / 3) +
      ggplot2::scale_colour_manual(values = pal, name = NULL) + ggplot2::scale_fill_manual(values = pal, name = NULL, guide = "none")
    if (layout == "panels") p <- p + ggplot2::facet_wrap(~ cause_label, scales = "free_y")
    else p <- p + ggplot2::aes(linetype = cause_label) + ggplot2::labs(linetype = NULL)
  }
  # ggplot2 does not wrap titles: wrap them here for a figure about 7 inches wide; the title and the
  # subtitle start at the left edge of the whole plot, not of the panel
  # the width is that of the open device (where the plot is usually printed), else 7 inches
  win <- (if (grDevices::dev.cur() > 1L) grDevices::dev.size("in")[1L] else 7) - 0.5
  wrapc <- function(z, pt, face) if (is.null(z) || !length(z)) NULL else
    paste(.wrap_inches(z, win, pt / 12, ps = 12, font = if (face == "bold") 2L else 1L), collapse = "\n")
  p <- p + ggplot2::labs(x = xlab, y = yl, title = wrapc(main, 13 * cex, "bold"), subtitle = wrapc(cap, 10 * cex, "plain")) + ggplot2::theme_minimal(base_size = 12 * cex) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(), plot.title.position = "plot",
                   plot.title = ggplot2::element_text(face = "bold", size = 13 * cex),
                   plot.subtitle = ggplot2::element_text(size = 10 * cex, colour = "grey30", lineheight = 1.1),
                   plot.margin = ggplot2::margin(8, 12, 6, 8))
  if (percent) p <- p + ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 * v), "%"))
  if (!is.null(xlim) || !is.null(ylim)) p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  if (length(keys) == 1L && layout != "stacked") p <- p + ggplot2::theme(legend.position = "none")
  p
}

# plot() of the regression models and of aj() draws their cumulative incidence curves
S7::method(plot, DirectCIFFit) <- function(x, ...) cif_plot(x, ...)
S7::method(plot, FineGrayFit) <- function(x, ...) cif_plot(x, ...)
S7::method(plot, CSCFit) <- function(x, ...) cif_plot(x, ...)
S7::method(plot, AJFit) <- function(x, ...) cif_plot(x, ...)

# ggplot2::autoplot(), registered when ggplot2 is loaded
autoplot <- S7::new_external_generic("ggplot2", "autoplot", "object")
S7::method(autoplot, DirectCIFFit) <- function(object, ...) cif_plot(object, ..., engine = "ggplot2")
S7::method(autoplot, FineGrayFit) <- function(object, ...) cif_plot(object, ..., engine = "ggplot2")
S7::method(autoplot, CSCFit) <- function(object, ...) cif_plot(object, ..., engine = "ggplot2")
S7::method(autoplot, AJFit) <- function(object, ...) cif_plot(object, ..., engine = "ggplot2")
