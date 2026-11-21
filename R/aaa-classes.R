# S7 classes of SurvX. Properties are read with `@`.

#' @noRd
SurvXFit <- S7::new_class("SurvXFit", package = "SurvX", abstract = TRUE,
  properties = list(
    model_call = S7::class_any,
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
    loglik_null = S7::class_numeric,
    entry = S7::class_any,            # start of the (entry, time] interval, or NULL
    strata = S7::class_any,           # 0-based stratum codes, or NULL
    strata_levels = S7::class_any,
    strata_vars = S7::class_any,
    cluster = S7::class_any,          # cluster identifiers, or NULL
    robust = S7::class_logical,       # vcov is the robust (sandwich) variance
    vcov_model = S7::class_any,       # inverse information
    tt = S7::class_any                # time-varying coefficients: list(fun, cols), or NULL
  ))

FineGrayFit <- S7::new_class("FineGrayFit", parent = SurvXFit, package = "SurvX",
  properties = list(
    cause = S7::class_integer,
    variance = S7::class_character,
    vcov_model = S7::class_any,
    vcov_robust = S7::class_any,
    loglik_null = S7::class_numeric,
    entry = S7::class_any,            # delayed entry (Geskus 2011), or NULL
    expanded = S7::class_any,         # weighted counting-process rows of the Geskus fit
    vcov_boot = S7::class_any         # bootstrap variance (variance = "bootstrap"), or NULL
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
    model_call = S7::class_any,
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
    X = S7::class_any,
    entry = S7::class_any,
    strata = S7::class_any,
    strata_levels = S7::class_any,
    strata_vars = S7::class_any,
    cluster = S7::class_any,
    robust = S7::class_logical
  ))

KMFit <- S7::new_class("KMFit", package = "SurvX",
  properties = list(
    model_call = S7::class_any,
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
    model_call = S7::class_any,
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
    model_call = S7::class_any,
    coefficients = S7::class_any,
    tables = S7::class_list,
    text = S7::class_character
  ))

ModelComparison <- S7::new_class("ModelComparison", package = "SurvX",
  properties = list(
    model_call = S7::class_any,
    table = S7::class_any,
    fit = S7::class_any,
    folds = S7::class_integer,
    cause = S7::class_integer
  ))

# result of model_scan(): one row per candidate model in `table`
ModelScan <- S7::new_class("ModelScan", package = "SurvX",
  properties = list(
    model_call = S7::class_any,
    type = S7::class_character,
    K = S7::class_integer,
    cause = S7::class_integer,
    goal = S7::class_character,
    depth = S7::class_character,
    times = S7::class_numeric,
    grid = S7::class_numeric,
    folds = S7::class_integer,
    seed = S7::class_integer,
    n = S7::class_integer,
    counts = S7::class_integer,
    plateau = S7::class_logical,
    table = S7::class_any,
    per_time = S7::class_any,
    brier = S7::class_any,
    effects = S7::class_any,
    fits = S7::class_list,
    nonparametric = S7::class_any,
    data = S7::class_any,
    time = S7::class_numeric,
    status = S7::class_integer
  ))

# parametric survival regression: coefficients = (Intercept), covariates, shape parameter
ParSurvFit <- S7::new_class("ParSurvFit", parent = SurvXFit, package = "SurvX",
  properties = list(
    dist = S7::class_character,
    scale = S7::class_character,
    dist_code = S7::class_integer,
    aft = S7::class_integer,
    entry = S7::class_any,
    loglik_null = S7::class_numeric,
    event_label = S7::class_character,
    knots = S7::class_numeric
  ))

# parametric cause-specific hazards: one ParSurvFit per cause
ParCSCFit <- S7::new_class("ParCSCFit", package = "SurvX",
  properties = list(
    model_call = S7::class_any,
    terms = S7::class_any,
    xlevels = S7::class_any,
    contrasts = S7::class_any,
    models = S7::class_list,
    K = S7::class_integer,
    dist = S7::class_character,
    scale = S7::class_character,
    n = S7::class_integer,
    time = S7::class_numeric,
    status = S7::class_integer,
    entry = S7::class_any,
    X = S7::class_any
  ))

# non-parametric hazard: Nelson-Aalen cumulative hazard and kernel-smoothed hazard by group and cause
HazardFit <- S7::new_class("HazardFit", package = "SurvX",
  properties = list(
    model_call = S7::class_any,
    cumhaz = S7::class_list,
    hazard = S7::class_list,
    bandwidth = S7::class_any,
    K = S7::class_integer,
    level = S7::class_numeric,
    time = S7::class_numeric,
    status = S7::class_integer,
    group = S7::class_any
  ))

# result of adjusted_risk(): risks by treatment value and their contrasts, by method
AdjustedRisk <- S7::new_class("AdjustedRisk", package = "SurvX",
  properties = list(
    model_call = S7::class_any,
    fit = S7::class_any,
    treatment = S7::class_character,
    levels = S7::class_character,
    times = S7::class_numeric,
    cause = S7::class_integer,
    kind = S7::class_character,
    method = S7::class_character,
    risk = S7::class_any,
    contrast = S7::class_any,
    curve = S7::class_any,
    ps = S7::class_any,
    B = S7::class_integer,
    ok = S7::class_integer,
    level = S7::class_numeric,
    n = S7::class_integer,
    ps_formula = S7::class_any
  ))

# result of time_effects(): one element of `units` per analysis (the Cox model, each cause of a
# cause-specific model, or the Fine-Gray model)
TimeEffects <- S7::new_class("TimeEffects", package = "SurvX",
  properties = list(
    model_call = S7::class_any,
    fit = S7::class_any,
    kind = S7::class_character,
    units = S7::class_list,
    settings = S7::class_list,
    forms = S7::class_character,
    bootstrap = S7::class_any,
    level = S7::class_numeric
  ))

# result of ph_test(): one element of `tests` per cause (one for a Cox model)
PHTest <- S7::new_class("PHTest", package = "SurvX",
  properties = list(
    tests = S7::class_list,
    transform = S7::class_character,
    table = S7::class_any
  ))

# parametric cure model: latency coefficients first (as ParSurvFit), then the logistic model of the
# cured fraction
CureFit <- S7::new_class("CureFit", parent = SurvXFit, package = "SurvX",
  properties = list(
    dist = S7::class_character,
    scale = S7::class_character,
    dist_code = S7::class_integer,
    aft = S7::class_integer,
    entry = S7::class_any,
    knots = S7::class_numeric,
    mixture = S7::class_logical,
    Xc = S7::class_any,
    cure_terms = S7::class_any,
    cure_xlevels = S7::class_any,
    cure_contrasts = S7::class_any,
    n_latency = S7::class_integer
  ))
