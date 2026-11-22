// SurvX: C API of the further performance measures (survcore/assess.hpp).
#include "survcore/survcore.h"

#include <cmath>
#include <limits>
#include <vector>

#include "survcore/assess.hpp"

using namespace survcore;

namespace {
const double kMissing = std::numeric_limits<double>::quiet_NaN();
bool finite_n(const double* x, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
  return true;
}
int check_outcome(int n, const double* time, const int* cause) {
  if (n <= 0 || !time || !cause) return -1;
  if (!finite_n(time, static_cast<std::size_t>(n))) return -2;
  for (int i = 0; i < n; ++i) if (cause[i] < 0) return -1;
  return 0;
}
}  // namespace

extern "C" {

int sc_assess(int n, const double* time, const int* cause, int interest, int nt, const double* times, int M,
              const double* risk, int groups, double z, double* out, int* ref) try {
  if (int rc = check_outcome(n, time, cause)) return rc;
  if (nt < 1 || M < 1 || !times || !risk || !out || !ref || groups < 1 || interest < 1 || !(z > 0.0)) return -1;
  if (!finite_n(times, static_cast<std::size_t>(nt))) return -2;
  assess_models(static_cast<std::size_t>(n), time, cause, interest, static_cast<std::size_t>(nt), times,
                static_cast<std::size_t>(M), risk, groups, z, out, ref);
  return 0;
} catch (...) { return -9; }

int sc_concordance(int n, const double* time, const int* cause, int interest, const double* risk, double tau, int uno,
                   double z, double* out, double* dfbeta) try {
  if (int rc = check_outcome(n, time, cause)) return rc;
  if (!risk || !out || interest < 1) return -1;
  if (!finite_n(risk, static_cast<std::size_t>(n)) || std::isnan(tau)) return -2;
  const ConcordanceResult r = concordance_index(static_cast<std::size_t>(n), time, cause, interest, risk, tau, uno != 0, dfbeta);
  out[0] = r.c; out[1] = r.se; out[2] = r.concordant; out[3] = r.discordant; out[4] = r.tied; out[5] = r.pairs;
  out[6] = out[7] = kMissing;
  if (std::isfinite(r.c) && std::isfinite(r.se) && z > 0.0) { out[6] = std::fmax(0.0, r.c - z * r.se); out[7] = std::fmin(1.0, r.c + z * r.se); }
  return 0;
} catch (...) { return -9; }

int sc_net_benefit(int n, const double* time, const int* cause, int interest, const double* pred, double t, int m,
                   const double* thr, double* nb, double* nb_all) try {
  if (int rc = check_outcome(n, time, cause)) return rc;
  if (!pred || m < 1 || !thr || !nb || !nb_all || interest < 1) return -1;
  if (!finite_n(pred, static_cast<std::size_t>(n)) || !std::isfinite(t)) return -2;
  std::vector<double> a, b;
  net_benefit(std::vector<double>(time, time + n), std::vector<int>(cause, cause + n), interest,
              std::vector<double>(pred, pred + n), t, std::vector<double>(thr, thr + m), a, b);
  for (int k = 0; k < m; ++k) { nb[k] = a[k]; nb_all[k] = b[k]; }
  return 0;
} catch (...) { return -9; }

int sc_smooth_calibration(int n, const double* time, const int* cause, int interest, const double* pred, double t,
                          double* observed, double* out) try {
  if (int rc = check_outcome(n, time, cause)) return rc;
  if (!pred || !out || interest < 1) return -1;
  if (!finite_n(pred, static_cast<std::size_t>(n)) || !std::isfinite(t)) return -2;
  std::vector<double> obs;
  const SmoothCalibration s = smooth_calibration(std::vector<double>(time, time + n), std::vector<int>(cause, cause + n),
                                                 interest, std::vector<double>(pred, pred + n), t, &obs);
  out[0] = s.ici; out[1] = s.e50; out[2] = s.e90;
  if (observed) for (int i = 0; i < n; ++i) observed[i] = s.ok ? obs[static_cast<std::size_t>(i)] : kMissing;
  return s.ok ? 0 : 1;
} catch (...) { return -9; }

int sc_mean_curve(int m, const double* grid, const double* y, double* out) try {
  if (m < 1 || !grid || !y || !out) return -1;
  *out = mean_curve(static_cast<std::size_t>(m), grid, y);
  return 0;
} catch (...) { return -9; }

}  // extern "C"
