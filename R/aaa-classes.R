# S7 classes of SurvX. Properties are read with `@`.

#' @noRd
SurvXFit <- S7::new_class("SurvXFit", package = "SurvX", abstract = TRUE,
  properties = list(
    call = S7::class_any,
    terms = S7::class_any,
    xlevels = S7::class_any,
    contrasts = S7::class_any,
    coefficients = S7::class_numeric,
    vcov = S7::class_any,
    loglik = S7::class_numeric,
    n = S7::class_integer,
    nevent = S7::class_integer,
    iterations = S7::class_integer,
    converged = S7::class_logical,
    flags = S7::class_integer,
    time = S7::class_numeric,
    status = S7::class_integer,
    X = S7::class_any
  ))

CoxFit <- S7::new_class("CoxFit", parent = SurvXFit, package = "SurvX",
  properties = list(
    ties = S7::class_character,
    weights = S7::class_any,
    loglik_null = S7::class_numeric
  ))

FineGrayFit <- S7::new_class("FineGrayFit", parent = SurvXFit, package = "SurvX",
  properties = list(
    cause = S7::class_integer,
    variance = S7::class_character,
    vcov_model = S7::class_any,
    vcov_robust = S7::class_any,
    loglik_null = S7::class_numeric
  ))

DirectCIFFit <- S7::new_class("DirectCIFFit", parent = SurvXFit, package = "SurvX",
  properties = list(
    K = S7::class_integer,
    baseline = S7::class_character,
    link = S7::class_character,
    base_code = S7::class_integer,
    alpha_free = S7::class_integer,
    alpha_value = S7::class_numeric,
    n_knots = S7::class_integer,
    knots = S7::class_list,
    constraint = S7::class_logical,
    diagnostics = S7::class_list,
    index = S7::class_any,
    estimation = S7::class_character,
    interest = S7::class_integer,
    cond = S7::class_character,
    vcov_model = S7::class_any
  ))

CSCFit <- S7::new_class("CSCFit", package = "SurvX",
  properties = list(
    call = S7::class_any,
    terms = S7::class_any,
    xlevels = S7::class_any,
    contrasts = S7::class_any,
    models = S7::class_list,
    K = S7::class_integer,
    ties = S7::class_character,
    product_limit = S7::class_logical,
    n = S7::class_integer,
    time = S7::class_numeric,
    status = S7::class_integer,
    X = S7::class_any
  ))

KMFit <- S7::new_class("KMFit", package = "SurvX",
  properties = list(
    call = S7::class_any,
    curves = S7::class_list,
    counts = S7::class_any,
    test = S7::class_any,
    level = S7::class_numeric,
    time = S7::class_numeric,
    status = S7::class_integer,
    group = S7::class_any
  ))

AJFit <- S7::new_class("AJFit", package = "SurvX",
  properties = list(
    call = S7::class_any,
    curves = S7::class_list,
    K = S7::class_integer,
    counts = S7::class_any,
    tests = S7::class_any,
    level = S7::class_numeric,
    time = S7::class_numeric,
    status = S7::class_integer,
    group = S7::class_any
  ))

# what summary() returns: a table that can be used in further computations and the text shown by print()
FitSummary <- S7::new_class("FitSummary", package = "SurvX",
  properties = list(
    call = S7::class_any,
    coefficients = S7::class_any,
    tables = S7::class_list,
    text = S7::class_character
  ))

ModelComparison <- S7::new_class("ModelComparison", package = "SurvX",
  properties = list(
    call = S7::class_any,
    table = S7::class_any,
    fit = S7::class_any,
    folds = S7::class_integer,
    cause = S7::class_integer
  ))
