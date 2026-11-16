// survcore: kernel-smoothed hazard rate from the Nelson-Aalen increments.
//
//   h(x) = sum_j b^-1 K_x((x - t_j) / b) dH_j,   dH_j = d_j / Y_j
//   Var h(x) = sum_j b^-2 K_x((x - t_j) / b)^2 d_j / Y_j^2
// (Ramlau-Hansen 1983, Ann Stat 11:453-466). K is the Epanechnikov kernel; near the ends of the
// follow-up, x < b or x > R - b with R the largest event time, the boundary kernels of Mueller and
// Wang (1994, Biometrics 50:61-76, eq. 2.4) keep the estimate free of the bias of the plain kernel:
//   K_q(z) = 12 (1+q)^-4 (z+1) [z(1-2q) + (3q^2-2q+1)/2],  -1 <= z <= q,  q = x/b,
// mirrored at the right end. The bandwidth minimises the least-squares cross-validation criterion
// of Patil (1993, Ann Stat 21:1792-1810) on [0, tau], where the risk set is not small; of its local
// minima the largest is taken (Hall & Marron 1991).
#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

namespace survcore {

inline double epan(double z) { return std::fabs(z) <= 1.0 ? 0.75 * (1.0 - z * z) : 0.0; }

// left boundary kernel; for the right end call with -z
inline double epan_boundary(double z, double q) {
  if (q >= 1.0) return epan(z);
  if (z < -1.0 || z > q) return 0.0;
  const double c = 12.0 / std::pow(1.0 + q, 4);
  return c * (z + 1.0) * (z * (1.0 - 2.0 * q) + (3.0 * q * q - 2.0 * q + 1.0) / 2.0);
}

// kernel weight of an event at t for the estimate at x with bandwidth b; R = right end of support
inline double kernel_weight(double x, double t, double b, double R, bool boundary) {
  const double z = (x - t) / b;
  if (!boundary) return epan(z);
  if (x < b) return epan_boundary(z, x / b);
  if (x > R - b) return epan_boundary(-z, (R - x) / b);
  return epan(z);
}

// et: distinct event times, d: events, Y: number at risk (length m). Estimates at the grid x (length nx).
inline void kernel_hazard(int m, const double* et, const double* d, const double* Y, double b, bool boundary,
                          int nx, const double* x, double* out_h, double* out_se) {
  const double R = m ? et[m - 1] : 0.0;
  for (int i = 0; i < nx; ++i) {
    double h = 0.0, v = 0.0;
    for (int j = 0; j < m; ++j) {
      const double k = kernel_weight(x[i], et[j], b, R, boundary);
      if (k == 0.0) continue;
      h += k * d[j] / Y[j];
      v += k * k * d[j] / (Y[j] * Y[j]);
    }
    out_h[i] = h / b;
    if (out_se) out_se[i] = std::sqrt(v) / b;
  }
}

// End of the region where the hazard is estimated reliably: the last event time at which at least
// max(10, 10% of the initial risk set) subjects are at risk. Beyond it the estimate rests on a
// handful of subjects; the cross-validation and the default plotting range stop there, as in the
// weighted criterion of Patil (1993), whose weight function vanishes in the right tail.
inline double hazard_support_end(int m, const double* et, const double* Y) {
  const double need = std::fmax(10.0, 0.1 * Y[0]);
  double end = et[0];
  for (int j = 0; j < m; ++j) if (Y[j] >= need) end = et[j];
  return end;
}

// least-squares cross-validation criterion of Patil (1993) at bandwidth b on [0, tau]:
//   CV(b) = int_0^tau h_b(x)^2 dx - 2 sum_{t_i <= tau} h_b^(-i)(t_i) dH_i,
// where h_b^(-i) leaves out the events at t_i (the tied events at t_i other than one are kept:
// weight (d_i - 1)/d_i). The integral by the trapezoidal rule on 401 points.
inline double hazard_cv(int m, const double* et, const double* d, const double* Y, double b, bool boundary, double tau) {
  const double R = et[m - 1];
  const int ng = 401;
  std::vector<double> xg(ng), hg(ng);
  for (int i = 0; i < ng; ++i) xg[i] = tau * i / (ng - 1);
  kernel_hazard(m, et, d, Y, b, boundary, ng, xg.data(), hg.data(), nullptr);
  double integral = 0.0;
  for (int i = 1; i < ng; ++i) integral += 0.5 * (hg[i] * hg[i] + hg[i - 1] * hg[i - 1]) * (xg[i] - xg[i - 1]);
  double cross = 0.0;
  int lo = 0;
  for (int i = 0; i < m && et[i] <= tau; ++i) {
    while (lo < m && et[lo] < et[i] - b) ++lo;
    double h = 0.0;
    for (int j = lo; j < m && et[j] <= et[i] + b; ++j) {
      const double w = j == i ? (d[i] - 1.0) / d[i] : 1.0;
      if (w == 0.0) continue;
      h += w * kernel_weight(et[i], et[j], b, R, boundary) * d[j] / Y[j];
    }
    cross += h / b * d[i] / Y[i];
  }
  return integral - 2.0 * cross;
}

// Bandwidth: 40 values log-spaced from tau/50 to tau/2; the largest local minimiser of the CV
// criterion, which is far more stable than the global minimiser (Hall & Marron 1991, J R Stat Soc B
// 53:245-252).
inline double hazard_bandwidth(int m, const double* et, const double* d, const double* Y, bool boundary) {
  const double tau = hazard_support_end(m, et, Y);
  if (!(tau > 0.0) || m < 3) return tau > 0.0 ? tau / 4.0 : 1.0;
  const int nb = 40;
  std::vector<double> bs(nb), cv(nb);
  for (int i = 0; i < nb; ++i) { bs[i] = tau / 50.0 * std::pow(25.0, i / (nb - 1.0)); cv[i] = hazard_cv(m, et, d, Y, bs[i], boundary, tau); }
  int best = -1;
  for (int i = 0; i < nb; ++i) {
    if (!std::isfinite(cv[i])) continue;
    const bool left = i == 0 || !std::isfinite(cv[i - 1]) || cv[i] <= cv[i - 1];
    const bool right = i == nb - 1 || !std::isfinite(cv[i + 1]) || cv[i] <= cv[i + 1];
    if (left && right) best = i;
  }
  return best >= 0 ? bs[best] : tau / 4.0;
}

}  // namespace survcore
