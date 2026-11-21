# The scientific sources of every method in SurvX, in one place.
# Every printed result ends with a "Methods:" line built from these entries, and sx_references()
# prints the full references of the methods behind a result.

.bib <- list(
  # Kaplan-Meier, Aalen-Johansen, tests, RMST
  kaplan_meier1958 = c("Kaplan & Meier (1958)", "Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete observations. Journal of the American Statistical Association 53:457-481."),
  greenwood1926 = c("Greenwood (1926)", "Greenwood M (1926). The natural duration of cancer. Reports on Public Health and Medical Subjects 33:1-26. HMSO, London."),
  brookmeyer1982 = c("Brookmeyer & Crowley (1982)", "Brookmeyer R, Crowley J (1982). A confidence interval for the median survival time. Biometrics 38:29-41."),
  aalen_johansen1978 = c("Aalen & Johansen (1978)", "Aalen OO, Johansen S (1978). An empirical transition matrix for non-homogeneous Markov chains based on censored observations. Scandinavian Journal of Statistics 5:141-150."),
  mantel1966 = c("Mantel (1966)", "Mantel N (1966). Evaluation of survival data and two new rank order statistics arising in its consideration. Cancer Chemotherapy Reports 50:163-170."),
  peto1972 = c("Peto & Peto (1972)", "Peto R, Peto J (1972). Asymptotically efficient rank invariant test procedures. Journal of the Royal Statistical Society, Series A 135:185-207."),
  harrington1982 = c("Harrington & Fleming (1982)", "Harrington DP, Fleming TR (1982). A class of rank test procedures for censored survival data. Biometrika 69:553-566."),
  gray1988 = c("Gray (1988)", "Gray RJ (1988). A class of K-sample tests for comparing the cumulative incidence of a competing risk. Annals of Statistics 16:1141-1154."),
  irwin1949 = c("Irwin (1949)", "Irwin JO (1949). The standard error of an estimate of expectation of life, with special reference to expectation of tumourless life in experiments with mice. Journal of Hygiene 47:188-189."),
  royston_parmar2013 = c("Royston & Parmar (2013)", "Royston P, Parmar MKB (2013). Restricted mean survival time: an alternative to the hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. BMC Medical Research Methodology 13:152."),
  # hazards
  nelson1972 = c("Nelson (1972)", "Nelson W (1972). Theory and applications of hazard plotting for censored failure data. Technometrics 14:945-966."),
  aalen1978 = c("Aalen (1978)", "Aalen OO (1978). Nonparametric inference for a family of counting processes. Annals of Statistics 6:701-726."),
  ramlau1983 = c("Ramlau-Hansen (1983)", "Ramlau-Hansen H (1983). Smoothing counting process intensities by means of kernel functions. Annals of Statistics 11:453-466."),
  muller_wang1994 = c("Mueller & Wang (1994)", "Mueller HG, Wang JL (1994). Hazard rate estimation under random censoring with varying kernels and bandwidths. Biometrics 50:61-76."),
  patil1993 = c("Patil (1993)", "Patil PN (1993). On the least squares cross-validation bandwidth in hazard rate estimation. Annals of Statistics 21:1792-1810."),
  hall_marron1991 = c("Hall & Marron (1991)", "Hall P, Marron JS (1991). Local minima in cross-validation functions. Journal of the Royal Statistical Society, Series B 53:245-252."),
  # Cox
  cox1972 = c("Cox (1972)", "Cox DR (1972). Regression models and life-tables. Journal of the Royal Statistical Society, Series B 34:187-220."),
  kawaguchi2021 = c("Kawaguchi et al. (2021)", "Kawaguchi ES, Shen JI, Suchard MA, Li G (2021). Scalable algorithms for large competing risks data. Journal of Computational and Graphical Statistics 30:685-693."),
  latouche2013 = c("Latouche et al. (2013)", "Latouche A, Allignol A, Beyersmann J, Labopin M, Fine JP (2013). A competing risks analysis should report results on all cause-specific hazards and cumulative incidence functions. Journal of Clinical Epidemiology 66:648-653."),
  austin_fine2017 = c("Austin & Fine (2017)", "Austin PC, Fine JP (2017). Practical recommendations for reporting Fine-Gray model analyses for competing risk data. Statistics in Medicine 36:4391-4400."),
  marubini1995 = c("Marubini & Valsecchi (1995)", "Marubini E, Valsecchi MG (1995). Analysing Survival Data from Clinical Trials and Observational Studies. Wiley, Chichester."),
  putter2007 = c("Putter et al. (2007)", "Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: competing risks and multi-state models. Statistics in Medicine 26:2389-2430."),
  andersen2012 = c("Andersen et al. (2012)", "Andersen PK, Geskus RB, de Witte T, Putter H (2012). Competing risks in epidemiology: possibilities and pitfalls. International Journal of Epidemiology 41:861-870."),
  stone1974 = c("Stone (1974)", "Stone M (1974). Cross-validatory choice and assessment of statistical predictions. Journal of the Royal Statistical Society, Series B 36:111-147."),
  breslow1972 = c("Breslow (1972)", "Breslow NE (1972). Discussion of the paper by D. R. Cox. Journal of the Royal Statistical Society, Series B 34:216-217."),
  breslow1974 = c("Breslow (1974)", "Breslow NE (1974). Covariance analysis of censored survival data. Biometrics 30:89-99."),
  andersen_gill1982 = c("Andersen & Gill (1982)", "Andersen PK, Gill RD (1982). Cox's regression model for counting processes: a large sample study. Annals of Statistics 10:1100-1120."),
  lin_wei1989 = c("Lin & Wei (1989)", "Lin DY, Wei LJ (1989). The robust inference for the Cox proportional hazards model. Journal of the American Statistical Association 84:1074-1078."),
  efron1977 = c("Efron (1977)", "Efron B (1977). The efficiency of Cox's likelihood function for censored data. Journal of the American Statistical Association 72:557-565."),
  # competing risks
  prentice1978 = c("Prentice et al. (1978)", "Prentice RL, Kalbfleisch JD, Peterson AV, Flournoy N, Farewell VT, Breslow NE (1978). The analysis of failure times in the presence of competing risks. Biometrics 34:541-554."),
  benichou1990 = c("Benichou & Gail (1990)", "Benichou J, Gail MH (1990). Estimates of absolute cause-specific risk in cohort studies. Biometrics 46:813-826."),
  fine_gray1999 = c("Fine & Gray (1999)", "Fine JP, Gray RJ (1999). A proportional hazards model for the subdistribution of a competing risk. Journal of the American Statistical Association 94:496-509."),
  jeong_fine2006 = c("Jeong & Fine (2006)", "Jeong JH, Fine J (2006). Direct parametric inference for the cumulative incidence function. Journal of the Royal Statistical Society, Series C 55:187-200."),
  jeong_fine2007 = c("Jeong & Fine (2007)", "Jeong JH, Fine JP (2007). Parametric regression on cumulative incidence function. Biostatistics 8:184-196."),
  dabrowska1988 = c("Dabrowska & Doksum (1988)", "Dabrowska DM, Doksum KA (1988). Estimation and testing in a two-sample generalized odds-rate model. Journal of the American Statistical Association 83:744-749."),
  geskus2011 = c("Geskus (2011)", "Geskus RB (2011). Cause-specific cumulative incidence estimation and the Fine and Gray model under both left truncation and right censoring. Biometrics 67:39-49."),
  lambert2017 = c("Lambert et al. (2017)", "Lambert PC, Wilkes SR, Crowther MJ (2017). Flexible parametric modelling of the cause-specific cumulative incidence function. Statistics in Medicine 36:1429-1446."),
  shi2013 = c("Shi, Cheng & Jeong (2013)", "Shi H, Cheng Y, Jeong JH (2013). Constrained parametric model for simultaneous inference of two cumulative incidence functions. Biometrical Journal 55:82-96."),
  cheng2009 = c("Cheng (2009)", "Cheng Y (2009). Modeling cumulative incidences of dementia and dementia-free death using a novel three-parameter logistic function. International Journal of Biostatistics 5(1), Article 29."),
  haile2016 = c("Haile et al. (2016)", "Haile SR, Jeong JH, Chen X, Cheng Y (2016). A 3-parameter Gompertz distribution for survival data with competing risks, with an application to breast cancer data. Journal of Applied Statistics 43:2239-2253."),
  larson1985 = c("Larson & Dinse (1985)", "Larson MG, Dinse GE (1985). A mixture model for the regression analysis of competing risks data. Journal of the Royal Statistical Society, Series C 34:201-211."),
  hinchliffe2013 = c("Hinchliffe & Lambert (2013)", "Hinchliffe SR, Lambert PC (2013). Flexible parametric modelling of cause-specific hazards to estimate cumulative incidence functions. BMC Medical Research Methodology 13:13."),
  # parametric models
  kalbfleisch2002 = c("Kalbfleisch & Prentice (2002)", "Kalbfleisch JD, Prentice RL (2002). The Statistical Analysis of Failure Time Data, 2nd ed. Wiley, Hoboken."),
  lawless2003 = c("Lawless (2003)", "Lawless JF (2003). Statistical Models and Methods for Lifetime Data, 2nd ed. Wiley, Hoboken."),
  klein2003 = c("Klein & Moeschberger (2003)", "Klein JP, Moeschberger ML (2003). Survival Analysis: Techniques for Censored and Truncated Data, 2nd ed. Springer, New York."),
  royston_parmar2002 = c("Royston & Parmar (2002)", "Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards and proportional-odds models for censored survival data, with application to prognostic modelling and estimation of treatment effects. Statistics in Medicine 21:2175-2197."),
  friedman1982 = c("Friedman (1982)", "Friedman M (1982). Piecewise exponential models for survival data with covariates. Annals of Statistics 10:101-113."),
  prentice1974 = c("Prentice (1974)", "Prentice RL (1974). A log gamma model and its maximum likelihood estimation. Biometrika 61:539-544."),
  cox2007 = c("Cox et al. (2007)", "Cox C, Chu H, Schneider MF, Munoz A (2007). Parametric survival analysis and taxonomy of hazard functions for the generalized gamma distribution. Statistics in Medicine 26:4352-4374."),
  # cure models
  boag1949 = c("Boag (1949)", "Boag JW (1949). Maximum likelihood estimates of the proportion of patients cured by cancer therapy. Journal of the Royal Statistical Society, Series B 11:15-44."),
  berkson1952 = c("Berkson & Gage (1952)", "Berkson J, Gage RP (1952). Survival curve for cancer patients following treatment. Journal of the American Statistical Association 47:501-515."),
  farewell1982 = c("Farewell (1982)", "Farewell VT (1982). The use of mixture models for the analysis of survival data with long-term survivors. Biometrics 38:1041-1046."),
  yakovlev1996 = c("Yakovlev & Tsodikov (1996)", "Yakovlev AY, Tsodikov AD (1996). Stochastic Models of Tumor Latency and Their Biostatistical Applications. World Scientific, Singapore."),
  chen1999 = c("Chen, Ibrahim & Sinha (1999)", "Chen MH, Ibrahim JG, Sinha D (1999). A new Bayesian model for survival data with a surviving fraction. Journal of the American Statistical Association 94:909-919."),
  lambert2007 = c("Lambert et al. (2007)", "Lambert PC, Thompson JR, Weston CL, Dickman PW (2007). Estimating and modeling the cure fraction in population-based cancer survival analysis. Biostatistics 8:576-594."),
  maller1996 = c("Maller & Zhou (1996)", "Maller RA, Zhou X (1996). Survival Analysis with Long-Term Survivors. Wiley, Chichester."),
  amico2018 = c("Amico & Van Keilegom (2018)", "Amico M, Van Keilegom I (2018). Cure models in survival analysis. Annual Review of Statistics and Its Application 5:311-342."),
  # diagnostics
  cox_snell1968 = c("Cox & Snell (1968)", "Cox DR, Snell EJ (1968). A general definition of residuals. Journal of the Royal Statistical Society, Series B 30:248-275."),
  schoenfeld1982 = c("Schoenfeld (1982)", "Schoenfeld D (1982). Partial residuals for the proportional hazards regression model. Biometrika 69:239-241."),
  therneau1990 = c("Therneau et al. (1990)", "Therneau TM, Grambsch PM, Fleming TR (1990). Martingale-based residuals for survival models. Biometrika 77:147-160."),
  grambsch1994 = c("Grambsch & Therneau (1994)", "Grambsch PM, Therneau TM (1994). Proportional hazards tests and diagnostics based on weighted residuals. Biometrika 81:515-526."),
  therneau2000 = c("Therneau & Grambsch (2000)", "Therneau TM, Grambsch PM (2000). Modeling Survival Data: Extending the Cox Model. Springer, New York."),
  collett2015 = c("Collett (2015)", "Collett D (2015). Modelling Survival Data in Medical Research, 3rd ed. Chapman and Hall/CRC, Boca Raton."),
  # prediction performance
  graf1999 = c("Graf et al. (1999)", "Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). Assessment and comparison of prognostic classification schemes for survival data. Statistics in Medicine 18:2529-2545."),
  gerds2006 = c("Gerds & Schumacher (2006)", "Gerds TA, Schumacher M (2006). Consistent estimation of the expected Brier score in general survival models with right-censored event times. Biometrical Journal 48:1029-1040."),
  gerds2007 = c("Gerds & Schumacher (2007)", "Gerds TA, Schumacher M (2007). Efron-type measures of prediction error for survival analysis. Biometrics 63:1283-1287."),
  kattan2018 = c("Kattan & Gerds (2018)", "Kattan MW, Gerds TA (2018). The index of prediction accuracy: an intuitive measure useful for evaluating risk prediction models. Diagnostic and Prognostic Research 2:7."),
  blanche2013 = c("Blanche et al. (2013)", "Blanche P, Dartigues JF, Jacqmin-Gadda H (2013). Estimating and comparing time-dependent areas under receiver operating characteristic curves for censored event times with competing risks. Statistics in Medicine 32:5381-5397."),
  gerds2014 = c("Gerds et al. (2014)", "Gerds TA, Andersen PK, Kattan MW (2014). Calibration plots for risk prediction models in the presence of competing risks. Statistics in Medicine 33:3191-3203."),
  hastie2009 = c("Hastie et al. (2009)", "Hastie T, Tibshirani R, Friedman J (2009). The Elements of Statistical Learning, 2nd ed. Springer, New York, Section 7.10."),
  # general inference
  akaike1974 = c("Akaike (1974)", "Akaike H (1974). A new look at the statistical model identification. IEEE Transactions on Automatic Control 19:716-723."),
  schwarz1978 = c("Schwarz (1978)", "Schwarz G (1978). Estimating the dimension of a model. Annals of Statistics 6:461-464."),
  wilks1938 = c("Wilks (1938)", "Wilks SS (1938). The large-sample distribution of the likelihood ratio for testing composite hypotheses. Annals of Mathematical Statistics 9:60-62."),
  white1982 = c("White (1982)", "White H (1982). Maximum likelihood estimation of misspecified models. Econometrica 50:1-25."),
  wald1943 = c("Wald (1943)", "Wald A (1943). Tests of statistical hypotheses concerning several parameters when the number of observations is large. Transactions of the American Mathematical Society 54:426-482."),
  casella2002 = c("Casella & Berger (2002)", "Casella G, Berger RL (2002). Statistical Inference, 2nd ed. Duxbury, Pacific Grove, Section 5.5 (delta method)."),
  andersen1993 = c("Andersen et al. (1993)", "Andersen PK, Borgan O, Gill RD, Keiding N (1993). Statistical Models Based on Counting Processes. Springer, New York."),
  efron1993 = c("Efron & Tibshirani (1993)", "Efron B, Tibshirani RJ (1993). An Introduction to the Bootstrap. Chapman and Hall, New York."),
  levenberg1944 = c("Levenberg (1944)", "Levenberg K (1944). A method for the solution of certain non-linear problems in least squares. Quarterly of Applied Mathematics 2:164-168."),
  marquardt1963 = c("Marquardt (1963)", "Marquardt DW (1963). An algorithm for least-squares estimation of nonlinear parameters. Journal of the Society for Industrial and Applied Mathematics 11:431-441.")
)

# what each key is used for, by topic; values are character vectors of keys
.topics <- list(
  km = c("kaplan_meier1958", "greenwood1926", "brookmeyer1982", "mantel1966", "peto1972", "harrington1982"),
  aj = c("aalen_johansen1978", "kalbfleisch2002", "marubini1995", "andersen1993", "gray1988", "putter2007"),
  logrank_test = c("mantel1966", "peto1972", "harrington1982"),
  gray_test = c("gray1988"),
  rmst = c("irwin1949", "royston_parmar2013", "klein2003"),
  hazard_np = c("nelson1972", "aalen1978", "ramlau1983", "muller_wang1994", "patil1993", "hall_marron1991"),
  cox = c("cox1972", "breslow1974", "breslow1972", "efron1977", "andersen_gill1982", "lin_wei1989", "kalbfleisch2002",
          "therneau1990", "grambsch1994", "therneau2000", "wilks1938", "akaike1974"),
  fg = c("fine_gray1999", "geskus2011", "white1982", "lin_wei1989", "efron1993", "kawaguchi2021", "latouche2013", "austin_fine2017"),
  csc = c("prentice1978", "cox1972", "breslow1974", "efron1977", "benichou1990", "kalbfleisch2002", "andersen_gill1982", "lin_wei1989",
          "latouche2013", "andersen2012", "hinchliffe2013", "lawless2003"),
  parsurv = c("kalbfleisch2002", "lawless2003", "klein2003", "royston_parmar2002", "friedman1982", "prentice1974",
              "cox2007", "casella2002", "akaike1974", "schwarz1978", "levenberg1944", "marquardt1963"),
  directcif = c("jeong_fine2006", "jeong_fine2007", "dabrowska1988", "casella2002", "akaike1974",
                "levenberg1944", "marquardt1963", "fine_gray1999", "geskus2011", "lambert2017", "shi2013",
                "larson1985"),
  cure = c("boag1949", "berkson1952", "farewell1982", "yakovlev1996", "chen1999", "lambert2007", "maller1996", "amico2018",
           "kalbfleisch2002", "casella2002"),
  residuals = c("cox_snell1968", "schoenfeld1982", "therneau1990", "grambsch1994"),
  ph_test = c("grambsch1994", "therneau2000", "schoenfeld1982"),
  dist_plot = c("klein2003", "collett2015", "nelson1972", "cox_snell1968"),
  hazard_plot = c("nelson1972", "ramlau1983", "muller_wang1994", "klein2003"),
  brier_score = c("graf1999", "gerds2006", "kattan2018"),
  auc_t = c("blanche2013"),
  calibration = c("gerds2014", "aalen_johansen1978"),
  compare_models = c("stone1974", "graf1999", "gerds2006", "kattan2018", "blanche2013", "gerds2014", "hastie2009", "gerds2007", "akaike1974"),
  cv_compare = c("hastie2009", "gerds2007", "graf1999", "gerds2006", "kattan2018", "blanche2013", "gerds2014"),
  cif_plot = c("aalen_johansen1978", "benichou1990", "efron1993", "casella2002"),
  cif_data = c("aalen_johansen1978", "benichou1990", "efron1993", "casella2002"),
  estimates = c("wald1943", "casella2002", "white1982"),
  predict = c("breslow1974", "benichou1990", "casella2002", "kalbfleisch2002"),
  Event = c("klein2003", "geskus2011")
)

# keys behind one fitted object: the method itself and the options that were used
.method_keys <- function(x) {
  if (S7::S7_inherits(x, CoxFit))
    return(c("cox1972", if (x@ties == "efron") "efron1977" else "breslow1974", "breslow1974", .design_keys(x)))
  if (S7::S7_inherits(x, FineGrayFit))
    return(c("fine_gray1999", if (!is.null(x@entry)) c("geskus2011", "andersen_gill1982"),
             if (identical(x@variance, "robust")) (if (is.null(x@entry)) "white1982" else "lin_wei1989"),
             if (identical(x@variance, "bootstrap")) "efron1993"))
  if (S7::S7_inherits(x, CSCFit))
    return(c("prentice1978", "cox1972", if (x@ties == "efron") "efron1977", "breslow1974", "benichou1990", .design_keys(x)))
  if (S7::S7_inherits(x, CureFit))
    return(c(if (x@mixture) c("boag1949", "berkson1952", "farewell1982") else c("yakovlev1996", "chen1999"), "lambert2007",
             "kalbfleisch2002", if (x@dist == "rp") "royston_parmar2002", if (x@dist == "pwexp") "friedman1982",
             if (x@dist == "gengamma") c("prentice1974", "cox2007"), if (!is.null(x@entry)) "klein2003", "casella2002"))
  if (S7::S7_inherits(x, ParSurvFit) || S7::S7_inherits(x, ParCSCFit)) {
    d <- unique(x@dist)
    k <- c(if (S7::S7_inherits(x, ParCSCFit)) c("prentice1978", "hinchliffe2013"), "kalbfleisch2002", "lawless2003",
           if (any(d == "rp")) "royston_parmar2002", if (any(d == "pwexp")) "friedman1982",
           if (any(d == "gengamma")) c("prentice1974", "cox2007"),
           if (!is.null(x@entry)) "klein2003", "casella2002")
    return(k)
  }
  if (S7::S7_inherits(x, DirectCIFFit)) {
    k <- switch(x@estimation,
                joint = c("jeong_fine2006", "jeong_fine2007"),
                weighted = c("jeong_fine2007", "fine_gray1999", "geskus2011", "lambert2017", "white1982"),
                constrained = c("shi2013", "cheng2009", if (x@K > 2L) "larson1985"))
    b <- unique(x@baseline)
    k <- c(k, if (any(b == "gompertz3")) "haile2016", if (any(b == "logistic3")) "cheng2009",
           if (any(b == "spline")) "royston_parmar2002", if (any(startsWith(b, "mix_"))) "larson1985",
           if (any(!x@link %in% c("PH", "PO"))) "dabrowska1988", "casella2002")
    return(k)
  }
  if (S7::S7_inherits(x, KMFit)) return(c("kaplan_meier1958", "greenwood1926", "brookmeyer1982", if (!is.null(x@test)) "mantel1966"))
  if (S7::S7_inherits(x, AJFit)) return(c("aalen_johansen1978", "andersen1993", "kalbfleisch2002", if (!is.null(x@tests)) "gray1988"))
  if (S7::S7_inherits(x, HazardFit)) return(c("nelson1972", "aalen1978", if (!is.null(x@bandwidth)) c("ramlau1983", "muller_wang1994", "patil1993", "hall_marron1991")))
  if (S7::S7_inherits(x, ModelComparison)) return(c("graf1999", "gerds2006", "kattan2018", "blanche2013", "gerds2014", if (x@folds >= 2L) "hastie2009"))
  if (S7::S7_inherits(x, PHTest)) return(c("grambsch1994", "therneau2000"))
  character(0)
}

.design_keys <- function(x) {
  c(if (!is.null(x@strata)) "kalbfleisch2002", if (!is.null(x@entry)) "andersen_gill1982", if (isTRUE(x@robust)) "lin_wei1989",
    if (S7::S7_inherits(x, CoxFit) && !is.null(x@tt)) "therneau2000")
}

# the line printed under every result
.method_lines <- function(keys) {
  keys <- unique(keys[keys %in% names(.bib)])
  if (!length(keys)) return(character(0))
  short <- gsub(" ", "\001", vapply(.bib[keys], `[[`, character(1), 1L), fixed = TRUE)   # keep each citation on one line
  txt <- paste0("Methods: ", paste(short, collapse = "; "), ". Full references: sx_references(fit).")
  gsub("\001", " ", strwrap(txt, width = 100L, exdent = 2L), fixed = TRUE)
}

#' References of the methods behind a result
#'
#' @description
#' Every result of SurvX rests on published methods. \code{sx_references()} prints the full references
#' of the methods that produced a fitted object, or of a function of the package, so that they can be
#' cited and checked. The printed output of every model ends with the short form of the same list
#' (the \dQuote{Methods:} line).
#' @param x a fitted object (from \code{\link{cox}}, \code{\link{fg}}, \code{\link{csc}},
#'   \code{\link{directcif}}, \code{\link{parsurv}}, \code{\link{cure}}, \code{\link{km}}, \code{\link{aj}},
#'   \code{\link{hazard_np}}, \code{\link{compare_models}} or \code{\link{ph_test}}), or the name of a
#'   function of the package as a character string (for example \code{"brier_score"}); \code{NULL}
#'   lists all references of the package.
#' @details
#' For a fitted object the list depends on the options that were used: for example the Efron or
#' Breslow method for ties of a Cox model, the estimation method and baselines of a direct model, or
#' delayed entry and the distribution of a parametric model.
#' @return
#' invisibly, a character vector of full references.
#' @seealso
#' The \code{References} section of the help page of every function, which groups the references by
#' what they are used for; \code{\link{SurvX-package}} for an overview.
#' @examples
#' sx_references("ph_test")
#' set.seed(1)
#' d <- data.frame(time = rexp(50), status = rbinom(50, 1, 0.8), x = rnorm(50))
#' sx_references(cox(Event(time, status) ~ x, d))
#' @export
sx_references <- function(x = NULL) {
  keys <- if (is.null(x)) names(.bib)
          else if (is.character(x)) {
            if (!x %in% names(.topics)) stop("no references stored for '", x, "'; available: ",
                                             paste(names(.topics), collapse = ", "), call. = FALSE)
            .topics[[x]]
          } else .method_keys(x)
  keys <- unique(keys)
  if (!length(keys)) { cat("No references are stored for this object.\n"); return(invisible(character(0))) }
  full <- vapply(.bib[keys], `[[`, character(1), 2L)
  full <- full[order(full)]
  cat("References of the methods used", if (is.character(x)) paste0(" by ", x, "()") else "", ":\n\n", sep = "")
  for (r in full) cat(strwrap(r, width = 100L, exdent = 2L), "", sep = "\n")
  invisible(unname(full))
}
