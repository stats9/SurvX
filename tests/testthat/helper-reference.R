# The R implementations that were replaced by the C++ core (src/capi_support.cpp and later
# modules), kept as references: test-core-parity.R checks that the core gives the same numbers.

ref_cif_interval <- function(Fv, se, zq) {
  ok <- is.finite(Fv) & is.finite(se) & Fv > 0 & Fv < 1
  e <- suppressWarnings(zq * se / (Fv * log(Fv)))
  lower <- suppressWarnings(Fv^exp(-e)); upper <- suppressWarnings(Fv^exp(e))
  lower[!ok] <- NA_real_; upper[!ok] <- NA_real_
  list(lower = lower, upper = upper)
}

ref_step_at <- function(x, y, times, left = 0) {
  if (!length(x)) return(rep(left, length(times)))
  c(left, y)[findInterval(times, x) + 1L]
}

ref_frailty_lrt <- function(ll, ll0, boundary) {
  st <- if (boundary) 0 else max(0, 2 * (ll - ll0))
  list(statistic = st, p.value = 0.5 * stats::pchisq(st, 1, lower.tail = FALSE))
}

ref_frailty_frame <- function(x, level = 0.95) {
  a <- (1 - level) / 2
  ev <- as.vector(tapply(x@status, factor(x@cluster_index, seq_along(x@cluster)), sum))
  nn <- tabulate(x@cluster_index, length(x@cluster))
  if (x@boundary)
    return(data.frame(cluster = x@cluster, n = nn, events = ev, frailty = 1, lower = 1, upper = 1, stringsAsFactors = FALSE))
  data.frame(cluster = x@cluster, n = nn, events = ev, frailty = x@shape / x@rate,
             lower = stats::qgamma(a, x@shape, x@rate), upper = stats::qgamma(1 - a, x@shape, x@rate),
             stringsAsFactors = FALSE)
}

ref_breslow_cumhaz <- function(bh, X, beta, times) {
  H0 <- vapply(times, function(t) sum(bh$hazard[bh$time <= t]), numeric(1))
  outer(exp(drop(X %*% beta)), H0)
}

ref_na_frame <- function(inc, level) {
  H <- cumsum(inc$d / inc$Y); se <- sqrt(cumsum(inc$d / inc$Y^2))
  zq <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(time = inc$time, n.risk = inc$Y, n.event = inc$d, cumhaz = H, se = se,
             lower = H * exp(-zq * se / H), upper = H * exp(zq * se / H))
}

ref_na_increments <- function(time, status, k, entry = NULL) {
  ev <- status == k
  et <- sort(unique(time[ev]))
  if (!length(et)) return(NULL)
  d <- tabulate(match(time[ev], et), length(et))
  Y <- vapply(et, function(s) sum(time >= s & (if (is.null(entry)) TRUE else entry < s)), numeric(1))
  list(time = et, d = as.numeric(d), Y = Y)
}

ref_km_limits <- function(S, se, zq) {
  w <- suppressWarnings(exp(zq * se / S))
  list(lower = ifelse(S > 0, S / w, NA_real_), upper = ifelse(S > 0, pmin(1, S * w), NA_real_))
}

# adjusted_risk(): propensity model, weights, balance, weighted Aalen-Johansen
# multinomial (binary: logistic) model of the treatment and stabilised weights
ref_ar_ps <- function(d, treatment, lev, psf, trim) {
  mf <- stats::model.frame(psf, d, na.action = stats::na.pass)
  X <- stats::model.matrix(psf, mf)
  a <- factor(as.character(d[[treatment]]), levels = lev)
  P <- ref_ar_mlogit(X, a)
  pa <- P[cbind(seq_len(nrow(P)), as.integer(a))]
  marg <- as.numeric(table(a)) / length(a)
  w <- marg[as.integer(a)] / pa
  if (trim > 0) {
    q <- stats::quantile(w, c(trim, 1 - trim), names = FALSE)
    w <- pmin(pmax(w, q[1L]), q[2L])
  }
  Xc <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  smd <- if (ncol(Xc)) {
    ref <- a == lev[1L]
    do.call(rbind, lapply(lev[-1L], function(l) {
      o <- a == l
      s <- sqrt((apply(Xc[ref, , drop = FALSE], 2L, stats::var) + apply(Xc[o, , drop = FALSE], 2L, stats::var)) / 2)
      s[!is.finite(s) | s == 0] <- NA_real_
      wm <- function(i) colSums(Xc[i, , drop = FALSE] * w[i]) / sum(w[i])
      data.frame(covariate = colnames(Xc), level = l,
                 before = (colMeans(Xc[o, , drop = FALSE]) - colMeans(Xc[ref, , drop = FALSE])) / s,
                 after = (wm(o) - wm(ref)) / s, stringsAsFactors = FALSE, row.names = NULL)
    }))
  } else NULL
  ess <- vapply(lev, function(l) { i <- a == l; sum(w[i])^2 / sum(w[i]^2) }, numeric(1))
  diag <- list(ps_range = t(vapply(lev, function(l) range(P[, l]), numeric(2))),
               weights = data.frame(level = lev, n = as.integer(table(a)), ess = unname(ess),
                                    max_weight = vapply(lev, function(l) max(w[a == l]), numeric(1)), row.names = NULL),
               smd = smd, trim = trim)
  list(w = w, diag = diag)
}

# multinomial logistic regression by Newton-Raphson; returns the n x K matrix of probabilities
ref_ar_mlogit <- function(X, a) {
  K <- nlevels(a); n <- nrow(X); p <- ncol(X)
  Y <- matrix(0, n, K); Y[cbind(seq_len(n), as.integer(a))] <- 1
  b <- matrix(0, p, K - 1L)
  prob <- function(b) { e <- cbind(0, X %*% b); e <- exp(e - apply(e, 1L, max)); e / rowSums(e) }
  ll <- function(P) sum(Y * log(pmax(P, 1e-300)))
  P <- prob(b); l0 <- ll(P)
  for (it in seq_len(100L)) {
    g <- as.vector(crossprod(X, Y[, -1L, drop = FALSE] - P[, -1L, drop = FALSE]))
    H <- matrix(0, p * (K - 1L), p * (K - 1L))
    for (j in seq_len(K - 1L)) for (k in seq_len(K - 1L)) {
      wjk <- P[, j + 1L] * ((j == k) - P[, k + 1L])
      H[(j - 1L) * p + seq_len(p), (k - 1L) * p + seq_len(p)] <- crossprod(X * wjk, X)
    }
    step <- tryCatch(solve(H + diag(1e-8, nrow(H)), g), error = function(e) NULL)
    if (is.null(step)) break
    s <- 1
    repeat {
      bn <- b + matrix(step * s, p, K - 1L)
      Pn <- prob(bn); ln <- ll(Pn)
      if (is.finite(ln) && ln >= l0 - 1e-10) break
      s <- s / 2
      if (s < 1e-8) break
    }
    b <- bn; P <- Pn
    if (abs(ln - l0) < 1e-10 * (abs(l0) + 1)) break
    l0 <- ln
  }
  colnames(P) <- levels(a)
  P
}

# weighted Aalen-Johansen estimate of the cumulative incidence of cause k (one event type: 1 - KM)
ref_ar_wcif <- function(time, status, entry, w, k, tt) {
  et <- sort(unique(time[status > 0L]))
  if (!length(et)) return(rep(0, length(tt)))
  Y <- vapply(et, function(s) sum(w[time >= s & (if (is.null(entry)) TRUE else entry < s)]), numeric(1))
  dall <- vapply(et, function(s) sum(w[time == s & status > 0L]), numeric(1))
  dk <- vapply(et, function(s) sum(w[time == s & status == k]), numeric(1))
  ok <- Y > 0
  haz_all <- ifelse(ok, dall / Y, 0); haz_k <- ifelse(ok, dk / Y, 0)
  S <- cumprod(1 - haz_all)
  Sm <- c(1, S[-length(S)])
  F <- cumsum(Sm * haz_k)
  ref_step_at(et, F, tt, 0)
}

ref_ar_tables <- function(est, boot, lev, times, tt, level) {
  z <- stats::qnorm(1 - (1 - level) / 2); a <- (1 - level) / 2
  ms <- intersect(c("gformula", "ipw"), names(est))
  have_b <- !is.null(boot) && boot$ok >= max(10L, boot$B / 2)
  pick <- function(m) if (have_b) lapply(boot$res, `[[`, m) else NULL
  risk <- list(); con <- list(); curve <- list()
  for (m in ms) {
    R <- est[[m]]                                   # levels x tt
    bs <- pick(m)
    arr <- if (have_b) simplify2array(bs) else NULL # levels x tt x B
    for (i in seq_along(lev)) {
      sd <- if (have_b) apply(arr[i, , , drop = FALSE], 2L, stats::sd) else rep(NA_real_, length(tt))
      lo <- if (have_b) apply(arr[i, , , drop = FALSE], 2L, stats::quantile, probs = a, names = FALSE) else rep(NA_real_, length(tt))
      hi <- if (have_b) apply(arr[i, , , drop = FALSE], 2L, stats::quantile, probs = 1 - a, names = FALSE) else rep(NA_real_, length(tt))
      df <- data.frame(method = m, level = lev[i], time = tt, risk = R[i, ], se = sd, lower = lo, upper = hi, stringsAsFactors = FALSE)
      curve[[length(curve) + 1L]] <- df
      risk[[length(risk) + 1L]] <- df[df$time %in% times, , drop = FALSE]
    }
    for (i in seq_along(lev)[-1L]) {
      j <- tt %in% times
      dd <- R[i, j] - R[1L, j]; rr <- R[i, j] / R[1L, j]
      if (have_b) {
        D <- matrix(arr[i, j, ] - arr[1L, j, ], nrow = sum(j))
        Q <- matrix(log(arr[i, j, ] / arr[1L, j, ]), nrow = sum(j))
        sdd <- apply(D, 1L, stats::sd); sdl <- apply(Q, 1L, function(v) stats::sd(v[is.finite(v)]))
        dlo <- apply(D, 1L, stats::quantile, probs = a, names = FALSE); dhi <- apply(D, 1L, stats::quantile, probs = 1 - a, names = FALSE)
        rlo <- exp(apply(Q, 1L, function(v) stats::quantile(v[is.finite(v)], a, names = FALSE)))
        rhi <- exp(apply(Q, 1L, function(v) stats::quantile(v[is.finite(v)], 1 - a, names = FALSE)))
        p <- 2 * stats::pnorm(-abs(dd / sdd))
      } else { sdd <- dlo <- dhi <- rlo <- rhi <- p <- rep(NA_real_, sum(j)) }
      con[[length(con) + 1L]] <- data.frame(method = m, level = lev[i], reference = lev[1L], time = tt[j], difference = dd, se = sdd,
                                            lower = dlo, upper = dhi, ratio = rr, ratio_lower = rlo, ratio_upper = rhi, p = p,
                                            stringsAsFactors = FALSE)
    }
  }
  list(risk = do.call(rbind, risk), contrast = do.call(rbind, con), curve = do.call(rbind, curve))
}

# turnbull(): survival, curve, test, tables and medians
# survival P(T > t) (left = FALSE) or P(T >= t) (left = TRUE): estimate with the mass of an innermost
# interval spread uniformly over it, and the lower and upper limits over all spreads
ref_tb_surv <- function(tb, t, left = FALSE) {
  lo <- tb$lo; hi <- tb$hi; p <- tb$p
  pt <- lo == hi
  est <- lower <- upper <- numeric(length(t))
  for (m in seq_along(t)) {
    s <- t[m]
    sure <- if (left) ifelse(pt, lo >= s, lo >= s) else ifelse(pt, lo > s, lo >= s)
    amb <- !pt & lo < s & hi > s
    lower[m] <- sum(p[sure])
    upper[m] <- lower[m] + sum(p[amb])
    if (any(amb & !is.finite(hi) & p > 0)) est[m] <- NA_real_
    else est[m] <- lower[m] + sum(p[amb] * (hi[amb] - s) / (hi[amb] - lo[amb]))
  }
  list(estimate = est, lower = lower, upper = upper)
}

# the curve as a polyline: at every end point the left limit and the value
ref_tb_steps <- function(tb) {
  g <- sort(unique(c(0, tb$lo, tb$hi[is.finite(tb$hi)])))
  s1 <- ref_tb_surv(tb, g, left = TRUE)$estimate
  s2 <- ref_tb_surv(tb, g)$estimate
  time <- rep(g, each = 2L); surv <- as.vector(rbind(s1, s2))
  ok <- is.finite(surv)
  data.frame(time = time[ok], surv = surv[ok])
}

# weighted logrank test, Finkelstein (1986) scores, permutation central limit variance (Fay & Shaw 2010)
ref_tb_test <- function(L, R, g, max_iter, tol) {
  tb <- SurvX:::.tb_fit(L, R, NULL, max_iter, tol)
  cs <- c(rev(cumsum(rev(tb$p))), 0)            # cs[j] = sum_{l >= j} p_l
  SL <- cs[tb$a]; SR <- cs[tb$b + 1L]
  xlx <- function(s) ifelse(s > 0, s * log(s), 0)
  sc <- (xlx(SL) - xlx(SR)) / (SL - SR)
  N <- length(sc); ng <- as.vector(table(g))
  U <- as.vector(tapply(sc, g, sum)) - ng * mean(sc)
  vc <- sum((sc - mean(sc))^2) / (N - 1)
  V <- vc * (diag(ng, length(ng)) - outer(ng, ng) / N)
  G <- length(ng)
  stat <- if (G == 2L) U[1L]^2 / V[1L, 1L] else drop(crossprod(U[-G], solve(V[-G, -G, drop = FALSE], U[-G])))
  list(statistic = stat, df = G - 1L, p.value = stats::pchisq(stat, G - 1L, lower.tail = FALSE),
       table = data.frame(group = levels(g), n = ng, score = U, sd = sqrt(diag(V))))
}

ref_tb_table <- function(x, times) {
  out <- NULL
  a <- (1 - x@level) / 2
  for (lv in x@groups) {
    f <- x@fits[[lv]]
    s <- ref_tb_surv(f, times)
    lo <- hi <- rep(NA_real_, length(times))
    if (!is.null(x@boot)) {
      bm <- x@boot[[lv]]
      for (m in seq_along(times)) {
        j <- which.min(abs(x@times - times[m]))
        if (abs(x@times[j] - times[m]) < 1e-8 * max(1, times[m])) v <- bm[, j]
        else {
          # between two grid times: the bootstrap curves interpolated linearly
          j1 <- max(which(x@times < times[m]), 1L); j2 <- min(which(x@times > times[m]), length(x@times))
          if (j1 == j2 || !is.finite(times[m]) || times[m] > max(x@times)) next
          f <- (times[m] - x@times[j1]) / (x@times[j2] - x@times[j1])
          v <- (1 - f) * bm[, j1] + f * bm[, j2]
        }
        v <- v[is.finite(v)]
        if (length(v) >= 10L) {
          q <- stats::quantile(v, c(a, 1 - a), names = FALSE); lo[m] <- q[1L]; hi[m] <- q[2L]
        }
      }
    }
    out <- rbind(out, data.frame(group = lv, time = times, survival = s$estimate, lower = lo, upper = hi,
                                 bound_low = s$lower, bound_high = s$upper, stringsAsFactors = FALSE))
  }
  out
}

# time at which a curve given on a grid first drops to p (linear interpolation)
ref_tb_quant <- function(time, surv, p = 0.5) {
  ok <- is.finite(surv)
  time <- time[ok]; surv <- surv[ok]
  i <- which(surv <= p)
  if (!length(i)) return(NA_real_)
  i <- i[1L]
  if (i == 1L) return(time[1L])
  if (surv[i - 1L] == surv[i]) return(time[i])
  time[i - 1L] + (surv[i - 1L] - p) / (surv[i - 1L] - surv[i]) * (time[i] - time[i - 1L])
}

ref_tb_medians <- function(x) {
  a <- (1 - x@level) / 2
  do.call(rbind, lapply(x@groups, function(lv) {
    st <- ref_tb_steps(x@fits[[lv]])
    md <- ref_tb_quant(st$time, st$surv)
    ci <- c(NA_real_, NA_real_)
    if (!is.null(x@boot)) {
      bm <- apply(x@boot[[lv]], 1L, function(s) ref_tb_quant(x@times, s))
      bm <- bm[is.finite(bm)]
      if (length(bm) >= 0.5 * x@B && length(bm) >= 10L) ci <- stats::quantile(bm, c(a, 1 - a), names = FALSE)
    }
    data.frame(group = lv, median = md, lower = ci[1L], upper = ci[2L], stringsAsFactors = FALSE)
  }))
}

# model_scan(): curves, measures and ranking
ref_scan_integrate <- function(grid, y) {
  if (any(!is.finite(y))) return(NA_real_)
  x <- c(0, grid); v <- c(0, y)
  sum(diff(x) * (v[-1L] + v[-length(v)]) / 2) / max(grid)
}
ref_scan_np <- function(time, status, cause, at) {
  d <- data.frame(time = time, status = status)
  if (max(status) <= 1L) {
    cv <- km(Event(time, status) ~ 1, d)@curves[[1L]]
    1 - ref_step_at(cv$time, cv$surv, at, 1)
  } else {
    cv <- aj(Event(time, status) ~ 1, d)@curves[[1L]]
    ref_step_at(cv$time, cv[[paste0("cif.", cause)]], at, 0)
  }
}
ref_measures <- function(time, status, risk, times, cause, groups = 10L) {
  t(vapply(seq_along(times), function(j) {
    b <- tryCatch(brier_score(time, status, risk[, j], times[j], cause), error = function(e) c(NA_real_, NA_real_, NA_real_))
    a <- tryCatch(auc_t(time, status, risk[, j], times[j], cause), error = function(e) NA_real_)
    ce <- tryCatch(attr(calibration(time, status, risk[, j], times[j], cause, groups), "error"), error = function(e) NA_real_)
    c(brier = b[[1L]], brier_null = b[[2L]], ipa = b[[3L]], auc = a, calib = ce)
  }, numeric(5)))
}
ref_rank <- function(tab, foldibs) {
  ok <- is.finite(tab$ibs_cv) & tab$fitted
  tab$rank <- NA_integer_; tab$diff_best <- NA_real_; tab$se_diff <- NA_real_; tab$equivalent <- NA
  if (any(ok)) {
    best <- which(ok)[which.min(tab$ibs_cv[ok])]
    for (m in which(ok)) {
      dk <- foldibs[m, ] - foldibs[best, ]
      dk <- dk[is.finite(dk)]
      tab$diff_best[m] <- tab$ibs_cv[m] - tab$ibs_cv[best]
      tab$se_diff[m] <- if (m == best) 0 else if (length(dk) >= 2L) stats::sd(dk) / sqrt(length(dk)) else NA_real_
      tab$equivalent[m] <- m == best || (is.finite(tab$se_diff[m]) && tab$diff_best[m] <= tab$se_diff[m])
    }
    eqv <- ok & tab$equivalent %in% TRUE
    tab$rank[eqv] <- 1L
    rest <- which(ok & !eqv)
    tab$rank[rest[order(tab$ibs_cv[rest])]] <- sum(eqv) + seq_along(rest)
  }
  tab
}

# ph_test() and time_effects(): score tests of proportional hazards
# risk-set sums of a Cox model (Breslow): for every event row, the weighted mean x-bar(t) and covariance
# V(t) of the covariates over the subjects still at risk
ref_cox_risksets <- function(time, status, X, beta, w = NULL, entry = NULL, strata = NULL) {
  if (!is.null(entry) || !is.null(strata)) return(ref_cox_risksets_general(time, status, X, beta, w, entry, strata))
  n <- length(time); p <- ncol(X)
  if (is.null(w)) w <- rep(1, n)
  r <- w * exp(drop(X %*% beta))
  o <- order(time)
  ts <- time[o]
  rc <- function(v) rev(cumsum(rev(v)))                 # sum over j >= i in time order
  S0 <- rc(r[o])
  S1 <- apply(X[o, , drop = FALSE] * r[o], 2L, rc)
  if (!is.matrix(S1)) S1 <- matrix(S1, ncol = p)
  ev <- which(status == 1L)
  ev <- ev[order(time[ev])]
  first <- match(time[ev], ts)                           # first position of each event time: everybody with time >= t
  s0 <- S0[first]; xbar <- S1[first, , drop = FALSE] / s0
  # second moments
  S2 <- array(0, c(length(ev), p, p))
  for (a in seq_len(p)) for (b in seq_len(a)) {
    v <- rc(X[o, a] * X[o, b] * r[o])[first] / s0 - xbar[, a] * xbar[, b]
    S2[, a, b] <- v; S2[, b, a] <- v
  }
  list(rows = ev, time = time[ev], xbar = xbar, V = S2, w = w[ev])
}

# the same with (entry, time] intervals and strata: risk set of an event at t in stratum s is
# {j in s: entry_j < t <= time_j}; O(n x events), used only when entry or strata are present
ref_cox_risksets_general <- function(time, status, X, beta, w, entry, strata) {
  n <- length(time); p <- ncol(X)
  if (is.null(w)) w <- rep(1, n)
  if (is.null(entry)) entry <- rep(-Inf, n)
  if (is.null(strata)) strata <- integer(n)
  r <- w * exp(drop(X %*% beta) - max(drop(X %*% beta)))
  ev <- which(status == 1L)
  ev <- ev[order(strata[ev], time[ev])]
  xbar <- matrix(0, length(ev), p); V <- array(0, c(length(ev), p, p))
  for (q in seq_along(ev)) {
    k <- ev[q]; t <- time[k]
    at <- which(strata == strata[k] & entry < t & time >= t)
    ra <- r[at]; Xa <- X[at, , drop = FALSE]; s0 <- sum(ra)
    m <- colSums(Xa * ra) / s0
    xbar[q, ] <- m
    V[q, , ] <- crossprod(Xa * sqrt(ra)) / s0 - tcrossprod(m)
  }
  list(rows = ev, time = time[ev], xbar = xbar, V = V, w = w[ev])
}

# ----------------------------------------------------------------------------------------------
# test of proportional hazards
# ----------------------------------------------------------------------------------------------
ref_time_transform <- function(time, status, t_event, transform, entry = NULL, rows = NULL) {
  switch(transform,
         identity = t_event,
         log = log(t_event),
         # rank of the event time among all observed (stop) times, censored included, as in
         # survival::cox.zph (version 3.0 and later)
         rank = if (is.null(rows)) rank(t_event) else rank(time)[rows],
         km = {                                                # 1 - Kaplan-Meier just before each event time
           inc <- ref_na_increments(time, status, 1L, entry)
           S <- cumprod(1 - inc$d / inc$Y)
           1 - c(1, S)[match(t_event, inc$time)]
         })
}

ref_zph_core <- function(time, status, X, beta, vcov, transform, w = NULL, entry = NULL, strata = NULL) {
  rs <- ref_cox_risksets(time, status, X, beta, w, entry, strata)
  sch <- X[rs$rows, , drop = FALSE] - rs$xbar
  g <- ref_time_transform(time, status, rs$time, transform, entry, rs$rows)
  g <- g - mean(g)
  p <- ncol(X); nd <- length(rs$rows)
  U <- colSums(sch * g)
  Ibb <- apply(rs$V, c(2L, 3L), sum)
  Ibg <- apply(rs$V * g, c(2L, 3L), sum)
  Igg <- apply(rs$V * g^2, c(2L, 3L), sum)
  Ibb_inv <- solve(Ibb)
  A <- Igg - t(Ibg) %*% Ibb_inv %*% Ibg                 # efficient information for the time-varying terms
  chi <- vapply(seq_len(p), function(k) U[k]^2 / A[k, k], numeric(1))
  glob <- drop(t(U) %*% solve(A, U))
  tab <- data.frame(chisq = c(chi, glob), df = c(rep(1L, p), p),
                    p = stats::pchisq(c(chi, glob), c(rep(1L, p), p), lower.tail = FALSE),
                    row.names = c(colnames(X), "GLOBAL"))
  scaled <- sch %*% vcov * nd + matrix(beta, nd, p, byrow = TRUE)
  colnames(scaled) <- colnames(X)
  list(table = tab, x = g + mean(ref_time_transform(time, status, rs$time, transform, entry, rs$rows)), y = scaled, time = rs$time,
       transform = transform, beta = beta, vcov = vcov)
}

# function of time used by the screening test at the event rows of the risk-set sums
ref_te_gfun <- function(transform, time, status, entry, w, weighted) {
  force(transform)
  function(rows, tev) {
    if (transform == "km" && weighted) {
      et <- sort(unique(time[status == 1L]))
      d <- vapply(et, function(s) sum(w[status == 1L & time == s]), numeric(1))
      Y <- vapply(et, function(s) sum(w[time >= s & (if (is.null(entry)) TRUE else entry < s)]), numeric(1))
      S <- cumprod(1 - d / Y)
      return(1 - c(1, S)[match(tev, et)])
    }
    ref_time_transform(time, status, tev, transform, entry, rows)
  }
}

# score test of a constant effect (Grambsch & Therneau 1994), per column, per term and globally
ref_te_zph <- function(time, status, X, beta, vcov, gfun, w, entry, strata, term_of) {
  rs <- ref_cox_risksets(time, status, X, beta, w, entry, strata)
  sch <- X[rs$rows, , drop = FALSE] - rs$xbar
  g0 <- gfun(rs$rows, rs$time)
  g <- g0 - mean(g0)
  p <- ncol(X); nd <- length(rs$rows)
  U <- colSums(sch * g)
  Ibb <- apply(rs$V, c(2L, 3L), sum)
  Ibg <- apply(rs$V * g, c(2L, 3L), sum)
  Igg <- apply(rs$V * g^2, c(2L, 3L), sum)
  A <- Igg - t(Ibg) %*% solve(Ibb, Ibg)
  chi_col <- vapply(seq_len(p), function(j) U[j]^2 / A[j, j], numeric(1))
  terms <- unique(term_of)
  chi_term <- vapply(terms, function(tm) {
    j <- which(term_of == tm)
    drop(t(U[j]) %*% solve(A[j, j, drop = FALSE], U[j]))
  }, numeric(1))
  df_term <- vapply(terms, function(tm) sum(term_of == tm), integer(1))
  glob <- drop(t(U) %*% solve(A, U))
  scaled <- sch %*% vcov * nd + matrix(beta, nd, p, byrow = TRUE)
  colnames(scaled) <- colnames(X)
  list(columns = data.frame(column = colnames(X), term = term_of, chisq = chi_col, df = 1L,
                            p = stats::pchisq(chi_col, 1L, lower.tail = FALSE), stringsAsFactors = FALSE),
       terms = data.frame(term = terms, chisq = unname(chi_term), df = unname(df_term),
                          p = stats::pchisq(unname(chi_term), unname(df_term), lower.tail = FALSE), stringsAsFactors = FALSE),
       global = c(chisq = glob, df = p, p = stats::pchisq(glob, p, lower.tail = FALSE)),
       resid = list(time = rs$time, g = g0, y = scaled))
}

# time_effects(): bases, splitting, effects
ref_te_basis <- function(fs, et, df, breaks, brk = NULL) {
  tpos <- function(t) pmax(t, min(et) / 2)
  tm <- stats::median(et)
  switch(fs$type,
    constant = list(name = "constant", q = 0L, extra = 0L, B = NULL, label = "constant",
                    describe = "constant over time (proportional hazards)"),
    log = list(name = fs$name, q = 1L, extra = 0L, B = function(t) cbind(log(tpos(t) / tm)), tmed = tm,
               label = "log-linear in time", describe = paste0("log hazard ratio changes linearly in log(time); centred at t = ", format(signif(tm, 3L)))),
    linear = list(name = fs$name, q = 1L, extra = 0L, B = function(t) cbind(t - tm), tmed = tm,
                  label = "linear in time", describe = paste0("log hazard ratio changes linearly in time; centred at t = ", format(signif(tm, 3L)))),
    piecewise = {
      b <- brk %||% breaks
      b <- sort(unique(b))
      list(name = fs$name, q = length(b), extra = if (is.null(breaks)) 1L else 0L,
           B = function(t) matrix(vapply(b, function(s) as.numeric(t > s), numeric(length(t))), nrow = length(t)),
           breaks = b, label = if (length(b) == 1L) paste0("step at t = ", format(signif(b, 3L))) else "steps",
           describe = paste0("one hazard ratio before and one after t = ", paste(format(signif(b, 3L)), collapse = ", "),
                             if (is.null(breaks)) " (change point chosen from the data)" else ""))
    },
    spline = {
      lt <- log(et)
      kn <- if (df > 1L) stats::quantile(lt, seq_len(df - 1L) / df, names = FALSE, type = 7) else numeric(0)
      bk <- range(lt)
      list(name = fs$name, q = df, extra = 0L,
           B = function(t) { b <- splines::ns(log(tpos(t)), knots = kn, Boundary.knots = bk, intercept = FALSE)
                             matrix(as.numeric(b), nrow = length(t)) },
           knots = kn, boundary = bk, label = paste0("spline in log time (", df, " df)"),
           describe = paste0("smooth change: natural cubic spline in log(time) with ", df, " degrees of freedom"))
    },
    user = {
      g <- fs$fun
      test <- g(et)
      m <- if (is.null(dim(test))) matrix(test, ncol = 1L) else as.matrix(test)
      if (nrow(m) != length(et) || any(!is.finite(m)))
        stop("form '", fs$name, "': the function must return one finite value (or row) per time; log(t) needs t > 0", call. = FALSE)
      list(name = fs$name, q = ncol(m), extra = 0L,
           B = function(t) { r <- g(t); if (is.null(dim(r))) matrix(r, ncol = 1L) else as.matrix(r) },
           label = paste0("user form '", fs$name, "'"), describe = paste0("your function of time '", fs$name, "'"))
    })
}


# cut every row (start, stop] at the event times et strictly inside it (as cox() does for tt())
ref_te_split <- function(start, stop, status, et) {
  lo <- findInterval(start, et)
  hi <- findInterval(stop, et, left.open = TRUE)
  rows <- pmax(hi - lo, 0L) + 1L
  idx <- rep.int(seq_along(stop), rows)
  r <- sequence(rows)
  ci <- lo[idx] + r
  last <- r == rows[idx]
  list(idx = idx, stop = ifelse(last, stop[idx], et[pmin(ci, length(et))]),
       start = ifelse(r == 1L, start[idx], et[pmax(ci - 1L, 1L)]), status = ifelse(last, status[idx], 0L))
}

ref_te_split_size <- function(start, stop, et) sum(as.numeric(pmax(findInterval(stop, et, left.open = TRUE) - findInterval(start, et), 0L) + 1L))


ref_te_effects <- function(fin, V, map, term_of, tt, level) {
  b <- fin$coefficients; z <- stats::qnorm(1 - (1 - level) / 2)
  out <- lapply(map, function(m) {
    j <- m$col
    if (!length(m$tv)) {
      est <- rep(b[[j]], length(tt)); se <- rep(sqrt(V[j, j]), length(tt))
    } else {
      Bm <- m$basis$B(tt)
      a <- cbind(1, Bm)
      ix <- c(j, m$tv)
      est <- drop(a %*% b[ix])
      se <- sqrt(pmax(rowSums((a %*% V[ix, ix, drop = FALSE]) * a), 0))
    }
    data.frame(column = names(b)[j], term = term_of[j], form = m$form, time = tt, log_hr = est, se = se,
               hr = exp(est), lower = exp(est - z * se), upper = exp(est + z * se),
               p = 2 * stats::pnorm(-abs(est / se)), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}


# Fine-Gray expansion, timefix, robust variance
ref_fg_expand <- function(time, status, entry, cause) {
  n <- length(time)
  if (is.null(entry)) entry <- rep(0, n)
  se <- sort(entry); st <- sort(time)
  ct <- sort(unique(time[status == 0L]))
  Y <- findInterval(ct, se, left.open = TRUE) - findInterval(ct, st, left.open = TRUE)     # #{entry < c <= time}
  dc <- tabulate(match(time[status == 0L], ct), length(ct))
  cumG <- c(1, cumprod(1 - dc / Y))
  lt <- sort(unique(entry[entry > 0]))
  R <- findInterval(lt, se) - findInterval(lt, st)                                          # #{entry <= l < time}
  el <- tabulate(match(entry[entry > 0], lt), length(lt))
  Hrev <- c(rev(cumprod(rev(1 - el / R))), 1)                                               # prod over l_k, l_k+1, ...
  G_right <- function(a) cumG[findInterval(a, ct) + 1L]                                     # G(a)
  G_left <- function(a) cumG[findInterval(a, ct, left.open = TRUE) + 1L]                    # G(a-)
  H_after <- function(a) Hrev[findInterval(a, lt) + 1L]                                     # prod over l > a
  H_left <- function(a) Hrev[findInterval(a, lt, left.open = TRUE) + 1L]                    # H(a-): prod over l >= a
  jumps <- sort(unique(c(ct, lt)))
  tmax <- max(c(0, time[status == cause]))
  base <- data.frame(id = seq_len(n), start = entry, stop = time, status = as.integer(status == cause), w = 1)
  comp <- which(status != 0L & status != cause & time < tmax)
  ext <- lapply(comp, function(i) {
    Ti <- time[i]
    cuts <- c(Ti, jumps[jumps > Ti & jumps < tmax])
    stops <- c(cuts[-1L], tmax)
    den <- G_left(Ti) * H_left(Ti)
    w <- if (den > 0) G_right(cuts) * H_after(cuts) / den else rep(0, length(cuts))
    data.frame(id = i, start = cuts, stop = stops, status = 0L, w = w)
  })
  d <- do.call(rbind, c(list(base), ext))
  d[d$w > 0 & d$stop > d$start, , drop = FALSE]
}

ref_timefix <- function(y, tol = sqrt(.Machine$double.eps)) {
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

ref_cox_robust <- function(fit, time, status, X, weights, ties, entry, strata, cluster) {
  U <- SurvX:::.cox_score(fit$coefficients, time, status, X, weights, ties, entry, strata)
  if (!is.null(cluster)) U <- rowsum(U, cluster, reorder = FALSE)
  V <- fit$vcov
  if (!all(is.finite(V))) return(V)
  out <- V %*% crossprod(U) %*% V
  dimnames(out) <- dimnames(V)
  out
}

# cause-specific Cox: cumulative incidence with strata
ref_csc_cif_general <- function(object, Xn, sn, times, cause) {
  K <- object@K
  bh <- lapply(seq_len(K), function(k) SurvX:::.csc_baseline(object, k))
  B <- lapply(object@models, function(m) m$coefficients)
  out <- matrix(0, nrow(Xn), length(times), dimnames = list(NULL, paste0("t=", times)))
  for (i in seq_len(nrow(Xn))) {
    s <- if (is.null(sn)) 0L else sn[i]
    inc <- lapply(seq_len(K), function(k) {
      b <- bh[[k]][bh[[k]]$stratum == s, , drop = FALSE]
      data.frame(time = b$time, dH = diff(c(0, b$cumhaz)) * exp(sum(Xn[i, ] * B[[k]])))
    })
    grid <- sort(unique(unlist(lapply(inc, `[[`, "time"))))
    if (!length(grid)) next
    dH <- vapply(inc, function(d) { v <- numeric(length(grid)); v[match(d$time, grid)] <- d$dH; v }, numeric(length(grid)))
    dH <- matrix(dH, length(grid), K)
    tot <- rowSums(dH)
    Sstep <- if (object@product_limit) pmax(0, 1 - tot) else exp(-tot)
    Sminus <- c(1, cumprod(Sstep))[seq_along(grid)]
    F <- cumsum(Sminus * dH[, cause])
    out[i, ] <- SurvX:::.step_at(grid, F, times, 0)
  }
  out
}

# tt(): splitting of follow-up
ref_tt_split <- function(f, weights, tt) {
  time <- f$time; status <- f$status; n <- length(time)
  entry <- f$entry SurvX:::`%||%` rep(0, n)
  et <- sort(unique(time[status == 1L]))
  lo <- findInterval(entry, et)                      # event times <= entry
  hi <- findInterval(time, et, left.open = TRUE)     # event times < time
  k <- pmax(hi - lo, 0L)                             # cut points strictly inside (entry, time)
  rows <- k + 1L
  N <- sum(as.numeric(rows))
  if (N > 5e7) stop("tt(): splitting at every event time would create ", format(N, big.mark = ","),
                    " rows; use fewer distinct times (round them) or split follow-up into periods yourself", call. = FALSE)
  idx <- rep.int(seq_len(n), rows)
  r <- sequence(rows)
  ci <- lo[idx] + r                                  # index in et of the end of a piece that is a cut
  last <- r == rows[idx]
  stop_t <- ifelse(last, time[idx], et[pmin(ci, length(et))])
  start_t <- ifelse(r == 1L, entry[idx], et[pmax(ci - 1L, 1L)])
  X <- f$X[idx, , drop = FALSE]
  funs <- if (is.function(tt)) rep(list(tt), length(f$tt_cols)) else tt
  if (!is.list(funs) || length(funs) != length(f$tt_cols) || !all(vapply(funs, is.function, logical(1))))
    stop("tt must be a function of time, or a list with one function per tt() term", call. = FALSE)
  for (j in seq_along(f$tt_cols)) {
    g <- funs[[j]](stop_t)
    if (length(g) != length(stop_t) || any(!is.finite(g)))
      stop("tt function ", j, " must return a finite value for every time (log(t) needs t > 0)", call. = FALSE)
    X[, f$tt_cols[j]] <- X[, f$tt_cols[j]] * g
  }
  list(time = stop_t, status = ifelse(last, status[idx], 0L), X = X, entry = start_t,
       weights = if (is.null(weights)) NULL else weights[idx],
       strata = if (is.null(f$strata)) NULL else f$strata[idx],
       cluster = if (!is.null(f$cluster)) f$cluster[idx] else idx)
}

# residuals of parametric and cure models
ref_par_resid <- function(m, type) {
  SurvX:::.no_ic(m, "residuals()", "compare the fit with the Turnbull estimate: hazard_plot(fit, type = \"survival\")")
  Xn <- m@X; n <- m@n
  Ht <- vapply(seq_len(n), function(i) SurvX:::.par_predict_X(m, Xn[i, , drop = FALSE], m@time[i], "cumhaz", FALSE, 0.95)[1L, 1L], numeric(1))
  if (!is.null(m@entry)) Ht <- Ht - vapply(seq_len(n), function(i)
    if (m@entry[i] > 0) SurvX:::.par_predict_X(m, Xn[i, , drop = FALSE], m@entry[i], "cumhaz", FALSE, 0.95)[1L, 1L] else 0, numeric(1))
  mres <- m@status - Ht
  switch(type, coxsnell = Ht, martingale = mres,
         deviance = sign(mres) * sqrt(-2 * (mres + ifelse(m@status == 1L, log(pmax(m@status - mres, 1e-300)), 0))),
         stop("parametric models have martingale, deviance and Cox-Snell residuals", call. = FALSE))
}

ref_cure_resid <- function(object, type = "martingale") {
    n <- object@n
  H <- vapply(seq_len(n), function(i)
    SurvX:::.cure_predict_X(object, object@X[i, , drop = FALSE], object@Xc[i, , drop = FALSE], object@time[i], "cumhaz", FALSE, 0.95)[1L, 1L], numeric(1))
  if (!is.null(object@entry)) H <- H - vapply(seq_len(n), function(i)
    if (object@entry[i] > 0) SurvX:::.cure_predict_X(object, object@X[i, , drop = FALSE], object@Xc[i, , drop = FALSE], object@entry[i], "cumhaz", FALSE, 0.95)[1L, 1L] else 0,
    numeric(1))
  m <- object@status - H
  switch(type, coxsnell = H, martingale = m,
         deviance = sign(m) * sqrt(-2 * (m + ifelse(object@status == 1L, log(pmax(object@status - m, 1e-300)), 0))))
}

# Fine-Gray prediction and log-scale limits
ref_fg_predict <- function(object, newdata = NULL, times) {
  if (!is.null(object@entry)) {
    e <- object@expanded
    bh <- SurvX:::.cox_baseline_raw(object@coefficients, e$stop, e$status, object@X[e$id, , drop = FALSE], e$w, "breslow", e$start, NULL)
    H <- SurvX:::.step_at(bh$time, bh$cumhaz, times, 0)
  } else {
    n <- object@n; p <- length(object@coefficients)
    r <- .C("SX_fg_baseline", as.integer(n), as.integer(p), object@time, object@status, SurvX:::.rm(object@X), object@cause,
            as.double(object@coefficients), time = double(n), cumhaz = double(n), m = integer(1), rc = integer(1),
            NAOK = TRUE, PACKAGE = "SurvX")
    SurvX:::.sx_rc(r$rc, "baseline subdistribution hazard")
    k <- seq_len(r$m)
    H <- SurvX:::.step_at(r$time[k], r$cumhaz[k], times, 0)
  }
  lp <- drop(SurvX:::.sx_newX(object, newdata) %*% object@coefficients)
  out <- 1 - exp(-outer(exp(lp), H))
  dimnames(out) <- list(NULL, paste0("t=", times))
  out
}

ref_pred_limits <- function(est, se, type, level) {
  zq <- SurvX:::.zq(level)
  if (type %in% c("survival", "cif")) {
    lim <- SurvX:::.cif_interval(est, se, zq)
    return(list(lower = array(lim$lower, dim(est), dimnames(est)), upper = array(lim$upper, dim(est), dimnames(est))))
  }
  e <- suppressWarnings(zq * se / est)
  lo <- suppressWarnings(est * exp(-e)); up <- suppressWarnings(est * exp(e))
  bad <- !is.finite(est) | !is.finite(se) | !(est > 0)
  lo[bad] <- NA_real_; up[bad] <- NA_real_
  list(lower = lo, upper = up)
}
