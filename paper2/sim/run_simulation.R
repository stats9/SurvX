# Paper 2: simulation study. Run from the repository root:
#
#   Rscript paper2/sim/run_simulation.R quick     # a few replications per scenario, to test (minutes)
#   Rscript paper2/sim/run_simulation.R full      # 1000 replications per scenario (hours)
#
# Optional second argument: number of worker processes (default: all cores but one).
# Results: paper2/results/sim_<mode>.rds, then run paper2/sim/summarize.R.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "quick"
stopifnot(mode %in% c("quick", "full"))
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)
R <- if (mode == "quick") 20L else 1000L

library(SurvX)
here <- "paper2"
source(file.path(here, "sim", "sim_functions.R"))
dir.create(file.path(here, "results"), showWarnings = FALSE)
cat("SurvX", as.character(utils::packageVersion("SurvX")), "| mode", mode, "| replications", R, "| workers", cores, "\n")

scen_A <- expand.grid(n = c(300L, 1000L), alpha = c(0, 1), cmax = c(10, 4), KEEP.OUT.ATTRS = FALSE)
scen_B <- expand.grid(n = c(500L, 1500L), dgp = c("constrained", "csh"), cmax = c(12, 5), KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)

one_A <- function(r, s) {
  set.seed(100000L * s + r)
  sc <- scen_A[s, ]
  d <- A_gen(sc$n, sc$alpha, sc$cmax)
  f <- A_fit_all(d, sc$alpha)
  if (is.null(f)) return(NULL)
  cbind(scenario = s, rep = r, f, censored = mean(d$status == 0L))
}
one_B <- function(r, s) {
  set.seed(200000L * s + r)
  sc <- scen_B[s, ]
  d <- B_gen(sc$n, sc$dgp, sc$cmax)
  f <- B_fit_all(d)
  if (is.null(f)) return(NULL)
  list(est = cbind(scenario = s, rep = r, f$est, censored = mean(d$status == 0L)),
       coherence = if (!is.null(f$coherence)) cbind(scenario = s, rep = r, f$coherence) else NULL)
}

cl <- parallel::makeCluster(cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterEvalQ(cl, { library(SurvX); NULL })
parallel::clusterExport(cl, c("scen_A", "scen_B", "one_A", "one_B", ls(pattern = "^(A_|B_|cif_ci)")))

run <- function(fun, S, label) {
  out <- vector("list", S)
  for (s in seq_len(S)) {
    t0 <- proc.time()[["elapsed"]]
    out[[s]] <- parallel::parLapply(cl, seq_len(R), fun, s = s)
    cat(sprintf("%s scenario %d/%d: %.0f s\n", label, s, S, proc.time()[["elapsed"]] - t0))
  }
  out
}
resA <- run(one_A, nrow(scen_A), "Study A")
resB <- run(one_B, nrow(scen_B), "Study B")

A <- do.call(rbind, unlist(resA, recursive = FALSE))
Bl <- unlist(resB, recursive = FALSE)
B <- do.call(rbind, lapply(Bl, `[[`, "est"))
Bc <- do.call(rbind, lapply(Bl, `[[`, "coherence"))

# Study C: computing time on one data set per size (median of 3 runs)
timing <- function(expr) median(replicate(3L, system.time(expr)[["elapsed"]]))
C <- NULL
for (n in if (mode == "quick") c(1000L, 5000L) else c(1000L, 5000L, 20000L, 50000L)) {
  set.seed(n)
  dA <- A_gen(n, 0, 10); dB <- B_gen(n, "constrained", 12)
  C <- rbind(C, data.frame(n = n,
    weighted_PH = timing(directcif(Event(time, status) ~ z1 + z2, dA, estimation = "weighted")),
    weighted_PO = timing(directcif(Event(time, status) ~ z1 + z2, dA, link = "PO", estimation = "weighted")),
    joint_PH = timing(directcif(Event(time, status) ~ z1 + z2, dA)),
    fine_gray = timing(fg(Event(time, status) ~ z1 + z2, dA)),
    constrained_K3 = timing(directcif(Event(time, status) ~ z1 + z2, dB, estimation = "constrained", cond = "weibull"))))
  cat("Study C n =", n, "done\n")
}

truthA <- lapply(seq_len(nrow(scen_A)), function(s) A_truth(scen_A$alpha[s]))
truthB <- lapply(seq_len(nrow(scen_B)), function(s) B_truth(scen_B$dgp[s]))
saveRDS(list(mode = mode, R = R, version = as.character(utils::packageVersion("SurvX")), session = utils::sessionInfo(),
             scen_A = scen_A, scen_B = scen_B, truth_A = truthA, truth_B = truthB, A = A, B = B, B_coherence = Bc, C = C),
        file.path(here, "results", paste0("sim_", mode, ".rds")))
cat("saved", file.path(here, "results", paste0("sim_", mode, ".rds")), "\n")
