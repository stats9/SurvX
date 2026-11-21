// SurvX: registration of the .C entry points (generated from rwrap.h).
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "rwrap.h"

static const R_CMethodDef c_methods[] = {
  {"SX_km", (DL_FUNC) &SX_km, 10, NULL},
  {"SX_na", (DL_FUNC) &SX_na, 8, NULL},
  {"SX_rmst", (DL_FUNC) &SX_rmst, 7, NULL},
  {"SX_aj", (DL_FUNC) &SX_aj, 10, NULL},
  {"SX_cox", (DL_FUNC) &SX_cox, 15, NULL},
  {"SX_cox_baseline", (DL_FUNC) &SX_cox_baseline, 13, NULL},
  {"SX_cox2", (DL_FUNC) &SX_cox2, 19, NULL},
  {"SX_cox2_baseline", (DL_FUNC) &SX_cox2_baseline, 18, NULL},
  {"SX_cox_score", (DL_FUNC) &SX_cox_score, 15, NULL},
  {"SX_cox_resid", (DL_FUNC) &SX_cox_resid, 16, NULL},
  {"SX_cure_fit", (DL_FUNC) &SX_cure_fit, 23, NULL},
  {"SX_cure_predict", (DL_FUNC) &SX_cure_predict, 20, NULL},
  {"SX_fg", (DL_FUNC) &SX_fg, 14, NULL},
  {"SX_fg_baseline", (DL_FUNC) &SX_fg_baseline, 11, NULL},
  {"SX_csc_cif", (DL_FUNC) &SX_csc_cif, 15, NULL},
  {"SX_logrank", (DL_FUNC) &SX_logrank, 15, NULL},
  {"SX_gray", (DL_FUNC) &SX_gray, 15, NULL},
  {"SX_direct_nparam", (DL_FUNC) &SX_direct_nparam, 6, NULL},
  {"SX_direct_knots", (DL_FUNC) &SX_direct_knots, 8, NULL},
  {"SX_direct_fit", (DL_FUNC) &SX_direct_fit, 22, NULL},
  {"SX_direct_eval", (DL_FUNC) &SX_direct_eval, 16, NULL},
  {"SX_direct_predict", (DL_FUNC) &SX_direct_predict, 18, NULL},
  {"SX_brier", (DL_FUNC) &SX_brier, 10, NULL},
  {"SX_auc", (DL_FUNC) &SX_auc, 8, NULL},
  {"SX_calibration", (DL_FUNC) &SX_calibration, 13, NULL},
  {"SX_compare", (DL_FUNC) &SX_compare, 20, NULL},
  {"SX_subdist_fit", (DL_FUNC) &SX_subdist_fit, 23, NULL},
  {"SX_subdist_var_ipcw", (DL_FUNC) &SX_subdist_var_ipcw, 15, NULL},
  {"SX_constrained_fit", (DL_FUNC) &SX_constrained_fit, 18, NULL},
  {"SX_constrained_predict", (DL_FUNC) &SX_constrained_predict, 16, NULL},
  {"SX_constrainedK_fit", (DL_FUNC) &SX_constrainedK_fit, 19, NULL},
  {"SX_constrainedK_predict", (DL_FUNC) &SX_constrainedK_predict, 18, NULL},
  {"SX_direct_predict_mean", (DL_FUNC) &SX_direct_predict_mean, 18, NULL},
  {"SX_constrained_predict_mean", (DL_FUNC) &SX_constrained_predict_mean, 16, NULL},
  {"SX_cif_bootstrap", (DL_FUNC) &SX_cif_bootstrap, 25, NULL},
  {"SX_par_fit", (DL_FUNC) &SX_par_fit, 20, NULL},
  {"SX_par_predict", (DL_FUNC) &SX_par_predict, 17, NULL},
  {"SX_parcsc_predict", (DL_FUNC) &SX_parcsc_predict, 19, NULL},
  {"SX_cv_folds", (DL_FUNC) &SX_cv_folds, 5, NULL},
  {"SX_kernel_hazard", (DL_FUNC) &SX_kernel_hazard, 12, NULL},
  {"SX_par_fit_ic", (DL_FUNC) &SX_par_fit_ic, 18, NULL},
  {"SX_direct_fit_ic", (DL_FUNC) &SX_direct_fit_ic, 23, NULL},
  {"SX_direct_eval_ic", (DL_FUNC) &SX_direct_eval_ic, 17, NULL},
  {"SX_turnbull", (DL_FUNC) &SX_turnbull, 16, NULL},
  {NULL, NULL, 0, NULL}
};

extern "C" void R_init_SurvX(DllInfo* dll) {
  R_registerRoutines(dll, c_methods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
