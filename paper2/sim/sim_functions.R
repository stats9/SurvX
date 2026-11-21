# Paper 2: data-generating models, true values and the estimators compared.
# Used by run_simulation.R. Every function is self-contained base R plus SurvX.

# ------------------------------------------------------------------------------------------
# Study A: one cause of interest with a plateau; the competing cause is not of the direct form
# ------------------------------------------------------------------------------------------
# F_1(t | z) = 1 - {1 + alpha E u(t)}^(-1/alpha), E = exp(b1 z1 + b2 z2),
# u(t) = tau (exp(rho t) - 1) / rho with rho < 0 (plateau). alpha = 0: 1 - exp(-E u(t)).
# With probability 1 - F_1(inf | z) the subject fails from cause 2 at a log-normal time whose
# mean depends on z1 (so the joint Gompertz model is misspecified for cause 2). Uniform
# censoring on (0, cmax).
A_par <- list(rho = -0.6, tau = 0.4, beta = c(z1 = 0.5, z2 = -0.3), t0 = 1.5)

A_u <- function(t, p = A_par) p$tau * expm1(p$rho * t) / p$rho
A_uinf <- function(p = A_par) -p$tau / p$rho

A_F1 <- function(t, z1, z2, alpha, p = A_par) {
  E <- exp(p$beta[["z1"]] * z1 + p$beta[["z2"]] * z2)
  u <- if (is.infinite(t)) A_uinf(p) else A_u(t, p)
  if (alpha == 0) -expm1(-E * u) else 1 - (1 + alpha * E * u)^(-1 / alpha)
}

A_gen <- function(n, alpha, cmax, p = A_par) {
  z1 <- rbinom(n, 1, 0.5); z2 <- rnorm(n)
  E <- exp(p$beta[["z1"]] * z1 + p$beta[["z2"]] * z2)
  p1 <- vapply(seq_len(n), function(i) A_F1(Inf, z1[i], z2[i], alpha, p), 0)
  is1 <- runif(n) < p1
  t <- numeric(n)
  v <- runif(n) * p1                                     # F_1 value, below the plateau
  u <- if (alpha == 0) -log1p(-v) / E else ((1 - v)^(-alpha) - 1) / (alpha * E)
  t1 <- log1p(p$rho * u / p$tau) / p$rho
  t[is1] <- t1[is1]
  t[!is1] <- rlnorm(sum(!is1), meanlog = 0.2 + 0.4 * z1[!is1], sdlog = 0.7)
  cens <- runif(n, 0, cmax)
  data.frame(time = pmax(pmin(t, cens), 1e-4), status = ifelse(t <= cens, ifelse(is1, 1L, 2L), 0L), z1 = z1, z2 = z2)
}

# true values of the estimands of study A
A_truth <- function(alpha, p = A_par) {
  c(b_z1 = p$beta[["z1"]], b_z2 = p$beta[["z2"]],
    Finf_0 = A_F1(Inf, 0, 0, alpha, p), Finf_1 = A_F1(Inf, 1, 0, alpha, p),
    Ft0_0 = A_F1(p$t0, 0, 0, alpha, p), Ft0_1 = A_F1(p$t0, 1, 0, alpha, p))
}

# Wald interval on the log(-log(1 - F)) scale for a probability
cif_ci <- function(F, se, z = qnorm(0.975)) {
  g <- log(-log1p(-F)); sg <- se / ((1 - F) * (-log1p(-F)))
  cbind(lower = -expm1(-exp(g - z * sg)), upper = -expm1(-exp(g + z * sg)))
}

# one fitted model -> estimates and standard errors of the six estimands
A_extract <- function(fit, t0, link_name) {
  cf <- coef(fit); V <- vcov(fit)
  i1 <- which(endsWith(names(cf), ":z1") & grepl("^1:", names(cf)))
  i2 <- which(endsWith(names(cf), ":z2") & grepl("^1:", names(cf)))
  nd <- data.frame(z1 = c(0, 1), z2 = 0)
  pr <- predict(fit, nd, times = c(Inf, t0), cause = 1, se = TRUE)
  est <- c(cf[i1], cf[i2], pr$cif[, 1], pr$cif[, 2])
  se <- c(sqrt(V[i1, i1]), sqrt(V[i2, i2]), pr$se[, 1], pr$se[, 2])
  names(est) <- names(se) <- c("b_z1", "b_z2", "Finf_0", "Finf_1", "Ft0_0", "Ft0_1")
  list(est = est, se = se)
}

# all estimators of study A on one data set; returns a long data frame
A_fit_all <- function(d, alpha, p = A_par) {
  link <- if (alpha == 0) "PH" else "PO"
  out <- list()
  add <- function(name, r) {
    if (is.null(r)) return()
    out[[length(out) + 1L]] <<- data.frame(estimator = name, estimand = names(r$est), est = unname(r$est),
                                           se = unname(r$se), stringsAsFactors = FALSE)
  }
  safe <- function(expr) tryCatch(expr, error = function(e) NULL)
  fw <- safe(directcif(Event(time, status) ~ z1 + z2, d, link = link, estimation = "weighted", cause = 1))
  if (!is.null(fw) && fw@converged) {
    add("weighted (IPCW variance)", safe(A_extract(fw, p$t0, link)))
    fr <- safe(directcif(Event(time, status) ~ z1 + z2, d, link = link, estimation = "weighted", cause = 1,
                         variance = "robust", init = coef(fw)))
    if (!is.null(fr)) add("weighted (censoring ignored)", safe(A_extract(fr, p$t0, link)))
  }
  fgor <- safe(directcif(Event(time, status) ~ z1 + z2, d, link = "GOR", estimation = "weighted", cause = 1))
  if (!is.null(fgor) && fgor@converged) {
    r <- safe(A_extract(fgor, p$t0, "GOR"))
    if (!is.null(r)) { r$est[1:2] <- NA; r$se[1:2] <- NA; add("weighted GOR (alpha estimated)", r) }
    fgr <- safe(directcif(Event(time, status) ~ z1 + z2, d, link = "GOR", estimation = "weighted", cause = 1,
                          variance = "robust", init = coef(fgor)))
    r2 <- if (is.null(fgr)) NULL else safe(A_extract(fgr, p$t0, "GOR"))
    if (!is.null(r2)) { r2$est[1:2] <- NA; r2$se[1:2] <- NA; add("weighted GOR (censoring ignored)", r2) }
  }
  fj <- safe(directcif(Event(time, status) ~ z1 + z2, d, link = link))
  if (!is.null(fj) && fj@converged) add("joint (cause 2 misspecified)", safe(A_extract(fj, p$t0, link)))
  if (alpha == 0) {
    ff <- safe(fg(Event(time, status) ~ z1 + z2, d, cause = 1))
    if (!is.null(ff)) {
      cf <- coef(ff); se <- sqrt(diag(vcov(ff)))
      add("Fine-Gray (semi-parametric)", list(est = c(b_z1 = cf[["z1"]], b_z2 = cf[["z2"]]),
                                              se = c(b_z1 = se[["z1"]], b_z2 = se[["z2"]])))
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# ------------------------------------------------------------------------------------------
# Study B: three causes, every subject eventually fails
# ------------------------------------------------------------------------------------------
# DGP "constrained": the K-cause constrained model (cause 1 Gompertz PH with a plateau; the
#   remaining risk shared between causes 2 and 3 by a logit in z1; Weibull G_k).
# DGP "csh": constant cause-specific hazards h_k(z) = lambda_k exp(g_k' z) (a Cox model for
#   every cause is correct; F_k(inf | z) = h_k / sum h).
B_par <- list(rho = -0.8, tau = 0.5, beta = c(z1 = 0.5, z2 = -0.3), gam = -0.3, delta = c(z1 = 0.6, z2 = 0),
              G = list(c(shape = 1.5, scale = 2), c(shape = 0.8, scale = 3)),
              lambda = c(0.12, 0.08, 0.05), g = rbind(c(0.5, -0.3), c(0.3, 0), c(-0.4, 0.3)), t0 = 2)

B_cif_constrained <- function(t, z1, z2, k, p = B_par) {
  E <- exp(p$beta[["z1"]] * z1 + p$beta[["z2"]] * z2)
  u <- if (is.infinite(t)) -p$tau / p$rho else p$tau * expm1(p$rho * t) / p$rho
  F1 <- -expm1(-E * u); p1 <- -expm1(E * p$tau / p$rho)
  if (k == 1L) return(F1)
  l3 <- p$gam + p$delta[["z1"]] * z1 + p$delta[["z2"]] * z2
  pi <- c(1, exp(l3)) / (1 + exp(l3))
  G <- if (is.infinite(t)) 1 else stats::pweibull(t, p$G[[k - 1L]][["shape"]], p$G[[k - 1L]][["scale"]])
  (1 - p1) * pi[k - 1L] * G
}

B_cif_csh <- function(t, z1, z2, k, p = B_par) {
  h <- p$lambda * exp(p$g[, 1] * z1 + p$g[, 2] * z2); H <- sum(h)
  (h[k] / H) * (if (is.infinite(t)) 1 else -expm1(-H * t))
}

B_gen <- function(n, dgp, cmax, p = B_par) {
  z1 <- rbinom(n, 1, 0.5); z2 <- rnorm(n)
  t <- numeric(n); cause <- integer(n)
  if (dgp == "constrained") {
    E <- exp(p$beta[["z1"]] * z1 + p$beta[["z2"]] * z2)
    p1 <- -expm1(E * p$tau / p$rho)
    l3 <- p$gam + p$delta[["z1"]] * z1 + p$delta[["z2"]] * z2
    pi3 <- exp(l3) / (1 + exp(l3))
    w <- runif(n)
    cause <- ifelse(w < p1, 1L, ifelse(runif(n) < pi3, 3L, 2L))
    v <- runif(n) * p1
    t1 <- log1p(p$rho * (-log1p(-v) / E) / p$tau) / p$rho
    t[cause == 1L] <- t1[cause == 1L]
    for (k in 2:3) {
      s <- cause == k
      t[s] <- stats::rweibull(sum(s), p$G[[k - 1L]][["shape"]], p$G[[k - 1L]][["scale"]])
    }
  } else {
    h <- sapply(1:3, function(k) p$lambda[k] * exp(p$g[k, 1] * z1 + p$g[k, 2] * z2))
    H <- rowSums(h)
    t <- stats::rexp(n, H)
    cause <- apply(h / H, 1, function(pr) sample.int(3L, 1L, prob = pr))
  }
  cens <- runif(n, 0, cmax)
  data.frame(time = pmax(pmin(t, cens), 1e-4), status = ifelse(t <= cens, cause, 0L), z1 = z1, z2 = z2)
}

B_truth <- function(dgp, p = B_par) {
  f <- if (dgp == "constrained") B_cif_constrained else B_cif_csh
  out <- c()
  for (k in 1:3) for (tt in c(Inf, p$t0)) for (zz in 0:1)
    out[paste0("F", k, "_", if (is.infinite(tt)) "inf" else "t0", "_", zz)] <- f(tt, zz, 0, k, p)
  out
}

# estimates of the 12 incidences (3 causes x {inf, t0} x z1 in {0, 1}, z2 = 0) and the share of
# the sample whose predicted incidences add up to more than 1 (at infinity and at the end of
# follow-up)
B_fit_all <- function(d, p = B_par) {
  safe <- function(expr) tryCatch(expr, error = function(e) NULL)
  nd <- data.frame(z1 = c(0, 1), z2 = 0)
  tend <- max(d$time)
  res <- list(); coh <- list()
  pred_direct <- function(fit, name) {
    est <- c(); se <- c()
    for (k in 1:3) {
      pr <- safe(predict(fit, nd, times = c(Inf, p$t0), cause = k, se = TRUE))
      if (is.null(pr)) return(NULL)
      for (j in 1:2) for (zz in 1:2) {
        nm <- paste0("F", k, "_", c("inf", "t0")[j], "_", zz - 1L)
        est[nm] <- pr$cif[zz, j]; se[nm] <- pr$se[zz, j]
      }
    }
    res[[name]] <<- data.frame(estimator = name, estimand = names(est), est = unname(est), se = unname(se[names(est)]),
                               stringsAsFactors = FALSE)
    tot_inf <- Reduce(`+`, lapply(1:3, function(k) predict(fit, d, times = Inf, cause = k)))
    tot_end <- Reduce(`+`, lapply(1:3, function(k) predict(fit, d, times = tend, cause = k)))
    coh[[name]] <<- data.frame(estimator = name, excess_inf = mean(tot_inf > 1 + 1e-8),
                               excess_end = mean(tot_end > 1 + 1e-8), max_inf = max(tot_inf), stringsAsFactors = FALSE)
  }
  fc <- safe(directcif(Event(time, status) ~ z1 + z2, d, estimation = "constrained", cause = 1, cond = "weibull"))
  if (!is.null(fc) && fc@converged) pred_direct(fc, "constrained, K = 3")
  fj <- safe(directcif(Event(time, status) ~ z1 + z2, d))
  if (!is.null(fj) && fj@converged) pred_direct(fj, "joint, unconstrained")
  fcs <- safe(csc(Event(time, status) ~ z1 + z2, d))
  if (!is.null(fcs)) {
    est <- c()
    for (k in 1:3) {
      pr <- predict(fcs, nd, times = p$t0, cause = k)
      for (zz in 1:2) est[paste0("F", k, "_t0_", zz - 1L)] <- pr[zz, 1]
    }
    res[["csc"]] <- data.frame(estimator = "cause-specific Cox", estimand = names(est), est = unname(est), se = NA_real_,
                               stringsAsFactors = FALSE)
  }
  if (!length(res)) return(NULL)
  list(est = do.call(rbind, res), coherence = if (length(coh)) do.call(rbind, coh) else NULL)
}
