# SurvX validation suite: settings, recording, data generators and helpers shared by all parts.
# Sourced by run_validation.R; every part writes its records through the rec_* functions below.

suppressPackageStartupMessages(library(SurvX))
# survival is attached so that its formula specials (strata) work inside the reference calls; SurvX
# reads strata() and cluster() itself and is not affected
if (requireNamespace("survival", quietly = TRUE)) suppressPackageStartupMessages(library(survival))

# ---- settings ----------------------------------------------------------------------------------------
# mode "quick": every check once, small replications (about 15-30 minutes);
# mode "full": the study for the paper (several hours; use cores).
make_config <- function(args) {
  mode <- if ("full" %in% args) "full" else "quick"
  cfg <- list(
    mode = mode,
    parts = 1:4,
    cores = max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L), na.rm = TRUE),
    seed = 20261116L,
    # simulation: replications and sample sizes
    sim_R = if (mode == "full") 1000L else 100L,
    sim_n = if (mode == "full") c(200L, 500L, 1000L) else 300L,
    # speed: sample sizes, replications per size, and the time limit after which a method is
    # not run at larger n
    speed_n = if (mode == "full") c(500L, 2000L, 10000L, 50000L, 200000L) else c(500L, 2000L, 10000L),
    speed_p = if (mode == "full") c(2L, 10L, 30L) else c(2L, 10L),
    speed_reps = if (mode == "full") 10L else 3L,
    speed_max_seconds = if (mode == "full") 300 else 60,
    # accuracy: replications of each simulated design
    acc_reps = if (mode == "full") 20L else 3L
  )
  for (a in grep("=", args, value = TRUE)) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    val <- strsplit(kv[2], ",", fixed = TRUE)[[1]]
    num <- suppressWarnings(as.numeric(val))
    cfg[[kv[1]]] <- if (all(is.finite(num))) (if (all(num == round(num))) as.integer(num) else num) else val
  }
  cfg
}

has <- function(p) requireNamespace(p, quietly = TRUE)
REF_PKGS <- c("survival", "cmprsk", "riskRegression", "prodlim", "flexsurv", "flexsurvcure", "survRM2", "timeROC",
              "icenReg", "interval", "frailtyEM", "parfm")

# ---- recording ---------------------------------------------------------------------------------------
.REC <- new.env()
rec_init <- function() for (k in c("accuracy", "speed", "simulation", "stability", "log")) assign(k, list(), envir = .REC)
rec_add <- function(table, row) {
  lst <- get(table, envir = .REC)
  lst[[length(lst) + 1L]] <- row
  assign(table, lst, envir = .REC)
}
rec_get <- function(table) {
  lst <- get(table, envir = .REC)
  if (!length(lst)) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(r) { for (c in setdiff(cols, names(r))) r[[c]] <- NA; r[cols] }))
}
logmsg <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(...))
  cat(m, "\n", sep = "")
  rec_add("log", data.frame(time = format(Sys.time()), message = paste0(...), stringsAsFactors = FALSE))
}

# value of expr, or an object of class "failed" with the message; warnings are collected
attempt <- function(expr) {
  warns <- character(0)
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(list(message = conditionMessage(e)), class = "failed")),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  attr(val, "warnings") <- warns
  val
}
failed <- function(x) inherits(x, "failed")
fail_msg <- function(x) if (failed(x)) x$message else ""

# ---- accuracy records --------------------------------------------------------------------------------
# One comparison of a SurvX quantity with the same quantity from a reference implementation.
# PASS when the largest absolute difference is at most tol_abs, or the largest relative difference
# (relative to max(|reference|, 1e-8)) is at most tol_rel. "exact" quantities are computed by the same
# formula in both packages (tolerance 1e-6); "optimiser" quantities come from a reference that
# maximises the likelihood with a general-purpose optimiser (tolerance 1e-3 relative).
# "info" quantities are computed by different (both valid) estimators; their difference is reported
# but not judged.
TOL <- list(exact = c(abs = 1e-6, rel = 1e-6), optimiser = c(abs = 1e-3, rel = 1e-3), loose = c(abs = 1e-2, rel = 1e-2),
            info = c(abs = Inf, rel = Inf))
check <- function(area, func, reference, data, quantity, ours, ref, tol = "exact", note = "") {
  row <- data.frame(area = area, func = func, reference = reference, data = data, quantity = quantity,
                    n_values = NA_integer_, max_abs_diff = NA_real_, max_rel_diff = NA_real_, tol = tol,
                    status = "ERROR", note = note, stringsAsFactors = FALSE)
  if (failed(ours) || failed(ref)) {
    row$note <- paste0(if (failed(ours)) paste("SurvX:", fail_msg(ours)) else "", if (failed(ref)) paste(" reference:", fail_msg(ref)) else "")
    rec_add("accuracy", row); return(invisible(row))
  }
  a <- suppressWarnings(as.numeric(unlist(ours))); b <- suppressWarnings(as.numeric(unlist(ref)))
  if (length(a) != length(b)) {
    row$note <- paste0("different lengths: SurvX ", length(a), ", reference ", length(b)); rec_add("accuracy", row); return(invisible(row))
  }
  both_na <- is.na(a) & is.na(b)
  ok <- is.finite(a) & is.finite(b)
  mismatch <- sum(xor(is.finite(a), is.finite(b)) & !both_na)
  row$n_values <- sum(ok)
  if (any(ok)) {
    d <- abs(a[ok] - b[ok])
    row$max_abs_diff <- max(d); row$max_rel_diff <- max(d / pmax(abs(b[ok]), 1e-8))
    tl <- TOL[[tol]]
    row$status <- if (tol == "info") "INFO" else if (mismatch > 0) "FAIL"
                  else if (row$max_abs_diff <= tl[["abs"]] || row$max_rel_diff <= tl[["rel"]]) "PASS" else "FAIL"
    if (mismatch > 0) row$note <- paste(row$note, mismatch, "values finite in one package only")
  } else row$note <- paste(row$note, "no finite values to compare")
  rec_add("accuracy", row)
  invisible(row)
}

# ---- data generators ---------------------------------------------------------------------------------
# One event. Weibull PH: S(t | x) = exp{-(t / scale)^shape exp(x'b)}; options for tied times, delayed
# entry, strata (baseline scale by stratum) and clusters of size 2 with a gamma frailty.
gen_surv <- function(n, p = 2L, beta = NULL, shape = 1.3, scale = 2, cens_max = 6, ties = 0L, entry = FALSE,
                     strata = FALSE, frailty = 0, plateau = 0) {
  X <- cbind(x1 = rbinom(n, 1, 0.5), matrix(rnorm(n * (p - 1L)), n, p - 1L, dimnames = list(NULL, paste0("x", seq_len(p - 1L) + 1L))))
  if (is.null(beta)) beta <- c(0.5, seq(-0.4, 0.4, length.out = p - 1L))
  lp <- drop(X %*% beta)
  g <- if (strata) sample(c("A", "B", "C"), n, TRUE) else rep("A", n)
  sc <- scale * c(A = 1, B = 0.6, C = 1.5)[g]
  id <- rep(seq_len(ceiling(n / 2)), each = 2L)[seq_len(n)]
  if (frailty > 0) lp <- lp + log(rgamma(max(id), 1 / frailty, 1 / frailty))[id]
  t <- sc * (-log(runif(n)) / exp(lp))^(1 / shape)
  if (plateau > 0) t[runif(n) < plateau] <- Inf
  cens <- runif(n, 0, cens_max)
  d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), X, g = g, id = id)
  if (ties > 0) d$time <- ceiling(d$time / max(d$time) * ties) * max(d$time) / ties
  if (entry) d$entry <- runif(n, 0, 0.5) * d$time * (runif(n) < 0.7)
  d
}

# Two causes with exponential (or Weibull) cause-specific hazards.
gen_cr <- function(n, p = 2L, ties = 0L, entry = FALSE, strata = FALSE) {
  X <- cbind(x1 = rbinom(n, 1, 0.5), matrix(rnorm(n * (p - 1L)), n, p - 1L, dimnames = list(NULL, paste0("x", seq_len(p - 1L) + 1L))))
  b1 <- c(0.6, seq(0.4, -0.4, length.out = p - 1L)); b2 <- c(-0.3, seq(-0.2, 0.3, length.out = p - 1L))
  g <- if (strata) sample(c("A", "B"), n, TRUE) else rep("A", n)
  m <- c(A = 1, B = 1.8)[g]
  t1 <- rweibull(n, 1.2, 1 / (0.35 * m * exp(drop(X %*% b1)))^(1 / 1.2))
  t2 <- rexp(n, 0.25 * exp(drop(X %*% b2)))
  cens <- runif(n, 0, 5)
  time <- pmin(t1, t2, cens)
  d <- data.frame(time = time, status = ifelse(time == cens, 0L, ifelse(t1 < t2, 1L, 2L)), X, g = g,
                  id = rep(seq_len(ceiling(n / 2)), each = 2L)[seq_len(n)])
  if (ties > 0) d$time <- ceiling(d$time / max(d$time) * ties) * max(d$time) / ties
  if (entry) d$entry <- runif(n, 0, 0.4) * d$time * (runif(n) < 0.7)
  d$dead <- as.integer(d$status > 0)
  d
}

# ---- real data (survival package) --------------------------------------------------------------------
real_data <- function() {
  if (!has("survival")) return(list())
  out <- list()
  l <- survival::lung
  out$lung <- stats::na.omit(data.frame(time = l$time, status = as.integer(l$status == 2), age = l$age,
                                        sex = factor(l$sex, 1:2, c("male", "female")), ecog = l$ph.ecog))
  v <- survival::veteran
  out$veteran <- data.frame(time = v$time, status = v$status, trt = factor(v$trt), karno = v$karno,
                            celltype = v$celltype, age = v$age)
  p <- survival::pbc[!is.na(survival::pbc$trt), ]
  out$pbc_cr <- stats::na.omit(data.frame(time = p$time, status = c(0L, 2L, 1L)[p$status + 1L],   # 1 death, 2 transplant
                                          age = p$age, bili = log(p$bili), albumin = p$albumin, edema = p$edema, trt = factor(p$trt)))
  m <- survival::mgus2
  m$time <- ifelse(m$pstat == 0, m$futime, m$ptime)
  m$status <- ifelse(m$pstat == 0, 2L * m$death, 1L)
  out$mgus2_cr <- stats::na.omit(data.frame(time = m$time, status = m$status, sex = m$sex, age = m$age, hgb = m$hgb,
                                            mspike = m$mspike))
  out$mgus2_cr <- out$mgus2_cr[out$mgus2_cr$time > 0, ]
  cc <- survival::colon[survival::colon$etype == 2, ]
  out$colon <- stats::na.omit(data.frame(time = cc$time, status = cc$status, rx = cc$rx, age = cc$age, nodes = cc$nodes,
                                         sex = factor(cc$sex)))
  out
}

# ---- small helpers -----------------------------------------------------------------------------------
rhs <- function(vars) paste(vars, collapse = " + ")
fml <- function(lhs, vars) stats::as.formula(paste(lhs, "~", rhs(vars)))
event_quantiles <- function(d, probs = c(0.25, 0.5, 0.75)) unname(stats::quantile(d$time[d$status > 0], probs, type = 1))

# step function value at times (right-continuous)
step_at <- function(x, y, at, y0 = 1) { i <- findInterval(at, x); c(y0, y)[i + 1L] }

# elapsed seconds of one call; fast calls are repeated until 0.5 s have passed and the total time is
# divided by the number of calls (the clock of Windows has a resolution of about 0.01 s)
time_it <- function(expr_fun, min_total = 0.5, max_inner = 200L) {
  t0 <- proc.time()[["elapsed"]]; v <- attempt(expr_fun()); t1 <- proc.time()[["elapsed"]] - t0
  if (failed(v)) return(list(seconds = NA_real_, value = v))
  if (t1 >= min_total) return(list(seconds = t1, value = v))
  k <- 0L; t0 <- proc.time()[["elapsed"]]
  repeat { expr_fun(); k <- k + 1L; el <- proc.time()[["elapsed"]] - t0; if (el >= min_total || k >= max_inner) break }
  list(seconds = el / k, value = v)
}

# parallel lapply that also works on Windows (PSOCK cluster; the workers load SurvX and WORKER_FILES)
plapply <- function(X, FUN, cores, export = character(0)) {
  if (cores <= 1L || length(X) < 2L) return(lapply(X, FUN))
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(files) { suppressPackageStartupMessages(library(SurvX)); for (f in files) source(f); NULL }, WORKER_FILES)
  if (length(export)) parallel::clusterExport(cl, export, envir = globalenv())
  parallel::parLapplyLB(cl, X, FUN)
}

fmt_num <- function(x, digits = 3) ifelse(is.na(x), "", formatC(x, format = "g", digits = digits))
