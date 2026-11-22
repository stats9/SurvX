# Memory and undefined-behaviour check of the C++ core

`memcheck_capi.cpp` calls the model entry points of the C API (`src/survcore/survcore.h`) on random
data of many sizes (ties, delayed entry, strata, clusters, interval-censored times, every baseline
and family, invalid input). `memcheck_support.cpp` does the same for the summary entry points of `src/capi_support.cpp`
(distributions, quantiles, intervals, frailty, adjusted risks, Turnbull summaries, prediction
measures, tests of proportional hazards, time-varying effects, data expansions); build it with the
same commands, replacing the file name. Both need only a C++17 compiler; R is not involved. From the package
root, on Linux or macOS:

```sh
# AddressSanitizer + UndefinedBehaviorSanitizer (any report stops the run)
g++ -std=c++17 -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined \
    -fno-sanitize-recover=undefined -Isrc validation/memcheck/memcheck_capi.cpp \
    src/capi.cpp src/capi_par.cpp src/capi_support.cpp src/capi_assess.cpp -o memcheck && ./memcheck

# UndefinedBehaviorSanitizer with clang, trap mode (no sanitizer runtime needed)
clang++ -std=c++17 -O1 -g -fsanitize=undefined -fsanitize-trap=undefined -Isrc \
    validation/memcheck/memcheck_capi.cpp src/capi.cpp src/capi_par.cpp src/capi_support.cpp src/capi_assess.cpp -o memcheck_ub && ./memcheck_ub

# valgrind (fewer repetitions, it is slow)
g++ -std=c++17 -O2 -g -DREPS=9 -Isrc validation/memcheck/memcheck_capi.cpp \
    src/capi.cpp src/capi_par.cpp src/capi_support.cpp src/capi_assess.cpp -o memcheck_vg
valgrind --error-exitcode=9 --leak-check=full --track-origins=yes ./memcheck_vg
```

A clean run prints the number of calls and `invalid inputs refused 11 of 11`, with no sanitizer
report and, for valgrind, `ERROR SUMMARY: 0 errors` and no leaks. Non-zero return codes in the
count are expected: they are invalid combinations (for example an AFT form of a family that has
none) refused by the core. Last run (2026-11-23): `memcheck_capi.cpp` 26592 calls under gcc 13 ASan + UBSan and clang 18 UBSan,
no report; valgrind 6648 calls, 0 errors, all heap blocks freed. `memcheck_support.cpp` (with an
installed index source, as the R front end uses, and every toolkit entry point): 4924 calls under
ASan + UBSan and clang UBSan, 1470 under valgrind, no report. The workflow `.github/workflows/memcheck.yaml`
runs the same checks on every push.
