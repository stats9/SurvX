// survcore: chi-square upper tail probability via the regularised incomplete
// gamma function Q(a, x). Series and continued-fraction evaluation follow
// Press et al. (2007), Numerical Recipes, 3rd ed., Section 6.2.
#pragma once
#include <cmath>
#include <limits>

namespace survcore {

inline double gamma_q(double a, double x) {
  if (!(x > 0.0)) return 1.0;
  if (!(a > 0.0)) return std::numeric_limits<double>::quiet_NaN();
  const double gln = std::lgamma(a);
  if (x < a + 1.0) {                       // series for P, return 1 - P
    double ap = a, sum = 1.0 / a, del = sum;
    for (int n = 0; n < 1000; ++n) {
      ap += 1.0; del *= x / ap; sum += del;
      if (std::fabs(del) < std::fabs(sum) * 1e-16) break;
    }
    return 1.0 - sum * std::exp(-x + a * std::log(x) - gln);
  }
  const double tiny = 1e-300;              // modified Lentz continued fraction for Q
  double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
  for (int i = 1; i < 1000; ++i) {
    const double an = -i * (i - a);
    b += 2.0;
    d = an * d + b; if (std::fabs(d) < tiny) d = tiny;
    c = b + an / c; if (std::fabs(c) < tiny) c = tiny;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < 1e-16) break;
  }
  return std::exp(-x + a * std::log(x) - gln) * h;
}

inline double chisq_sf(double x, double df) { return gamma_q(0.5 * df, 0.5 * x); }

}  // namespace survcore
