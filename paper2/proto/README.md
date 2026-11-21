# Reference implementations for paper 2

Independent Python code used to check the C++ core of SurvX for the methods of paper 2.
Not part of the R package (excluded by `.Rbuildignore`).

Build the C library once, from the package's `src` directory:

    g++ -O2 -std=c++17 -shared -fPIC -o /tmp/sx2.so capi.cpp capi_par.cpp -I.

(another path can be given with the environment variable `SURVX_LIB`). Then, from this directory:

| Script | Checks |
|---|---|
| `check_psi.py` | the analytic censoring term (psi) against the exact derivative of the total score with respect to case weights, the Kaplan-Meier estimate recomputed (Python only) |
| `check_cpp.py` | C++ variance with the censoring term against `wvar.py` (Gompertz; PH, PO, GOR; with and without ties) |
| `check_cpp_bases.py` | the same for every baseline family, with F(t; z) taken from the C++ prediction routine |
| `mc_quick.py` | small Monte Carlo: standard errors with and without the censoring term against the empirical SD |
| `check_ck.py` | K-cause constrained model: K = 2 against the two-cause model; score and information against finite differences (K = 3, 4) |
| `check_ck2.py` | K-cause constrained model: log-likelihood against independent code; long-term incidences add up to one; recovery in 300 simulated data sets |

`wvar.py`: weighted estimator, its variance with and without the censoring term, and the exact
infinitesimal jackknife. `ck.py`: ctypes bindings of the K-cause model. `ext.py`: ctypes bindings
of the other direct-model routines.
