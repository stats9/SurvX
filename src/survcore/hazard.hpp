// survcore: kernel-smoothed hazard rate from the Nelson-Aalen increments.
//
//   h(x) = sum_j b^-1 K_x((x - t_j) / b) dH_j,   dH_j = d_j / Y_j
//   Var h(x) = sum_j b^-2 K_x((x - t_j) / b)^2 d_j / Y_j^2
// (Ramlau-Hansen 1983, Ann Stat 11:453-466). K is the Epanechnikov kernel; near the ends of the
// follow-up, x < b or x > R - b with R the largest event time, the boundary kernels of Mueller and
// Wang (1994, Biometrics 50:61-76, eq. 2.4) keep the estimate free of the bias of the plain kernel:
//   K_q(z) = 12 (1+q)^-4 (z+1) [z(1-2q) + (3q^2-2q+1)/2],  -1 <= z <= q,  q = x/b,
// mirrored at the right end. The bandwidth minimises the least-squares cross-validation criterion
// of Patil (1993, Ann Stat 21:1792-1810):
//   CV(b) = int h_b(x)^2 dx - 2 sum_{i != j} K_b(t_i - t_j) dN_i dN_j / (Y_i Y_j).
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

// least-squares cross-validation criterion of Patil (1993) at bandwidth b; the integral of h^2 by the
// trapezoidal rule on 401 points over [0, R]
inline double hazard_cv(int m, const double* et, const double* d, const double* Y, double b, bool boundary) {
  const double R = et[m - 1];
  const int ng = 401;
  std::vector<double> xg(ng), hg(ng);
  for (int i = 0; i < ng; ++i) xg[i] = R * i / (ng - 1);
  kernel_hazard(m, et, d, Y, b, boundary, ng, xg.data(), hg.data(), nullptr);
  double integral = 0.0;
  for (int i = 1; i < ng; ++i) integral += 0.5 * (hg[i] * hg[i] + hg[i - 1] * hg[i - 1]) * (xg[i] - xg[i - 1]);
  double cross = 0.0;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j) {
      const double k = epan((et[i] - et[j]) / b) / b;
      if (k == 0.0) continue;
      // pairs of distinct events: d_i d_j for i != j, d_i (d_i - 1) for tied events at the same time
      const double pairs = i == j ? d[i] * (d[i] - 1.0) : d[i] * d[j];
      cross += k * pairs / (Y[i] * Y[j]);
    }
  return integral - 2.0 * cross;
}

// bandwidth on a grid of 40 values, log-spaced from R/100 to R/2, minimising the CV criterion
inline double hazard_bandwidth(int m, const double* et, const double* d, const double* Y, bool boundary) {
  const double R = et[m - 1];
  if (!(R > 0.0) || m < 3) return R > 0.0 ? R / 4.0 : 1.0;
  double best = R / 4.0, best_cv = INFINITY;
  for (int i = 0; i < 40; ++i) {
    const double b = R / 100.0 * std::pow(50.0, i / 39.0);
    const double cv = hazard_cv(m, et, d, Y, b, boundary);
    if (std::isfinite(cv) && cv < best_cv) { best_cv = cv; best = b; }
  }
  return best;
}

}  // namespace survcore
