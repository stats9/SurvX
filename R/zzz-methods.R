# Methods of the generics defined in zz-output.R for the parametric models (collated after it).

S7::method(estimates, ParSurvFit) <- function(object, level = 0.95, ...) .par_estimates(object, NA_integer_, level)
S7::method(estimates, ParCSCFit) <- function(object, level = 0.95, ...) {
  do.call(rbind, lapply(seq_len(object@K), function(k) .par_estimates(object@models[[k]], k, level)))
}

