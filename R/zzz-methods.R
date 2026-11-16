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
