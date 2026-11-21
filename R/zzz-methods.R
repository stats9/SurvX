# Methods of the generics defined in zz-output.R for the parametric models (collated after it).

S7::method(estimates, ParSurvFit) <- function(object, level = 0.95, ...) .par_estimates(object, NA_integer_, level)
S7::method(estimates, ParCSCFit) <- function(object, level = 0.95, ...) {
  do.call(rbind, lapply(seq_len(object@K), function(k) .par_estimates(object@models[[k]], k, level)))
}


# plot() of a parametric survival model: fitted survival curves with the Kaplan-Meier-type
# (Nelson-Aalen) estimate as reference; hazard_plot() gives the hazard and cumulative hazard
S7::method(plot, ParSurvFit) <- function(x, ...) hazard_plot(x, type = "survival", ...)

# the tests of ph_test() as a data frame
S7::method(as.data.frame, PHTest) <- function(x, ...) x@table

# cure models (generic estimates() is defined in zz-output.R)
S7::method(estimates, CureFit) <- function(object, level = 0.95, ...) {
  est <- object@coefficients; nm <- names(est)
  pl <- ncol(object@X); L <- object@n_latency; P <- length(est)
  role <- c("baseline", rep("beta", pl), rep("shape", L - 1L - pl), "cure_baseline", rep("cure", P - L - 1L))
  sc <- ifelse(role == "beta", if (object@scale == "AFT") "log time ratio (uncured)" else "log hazard ratio (uncured)",
               ifelse(role == "cure", "log odds ratio of being cured", ifelse(role == "cure_baseline", "logit of the cured fraction", nm)))
  .est_frame(NA_integer_, nm, role, sc, est, .se_of(object), level, role %in% c("beta", "cure"))
}

# S7 >= 1.0.0: remove the placeholders that `method<-` leaves in the namespace for external generics
# (run when the package is built; older S7 versions do not have, or need, this hook)
if ("S7_on_build" %in% getNamespaceExports("S7")) getExportedValue("S7", "S7_on_build")()
