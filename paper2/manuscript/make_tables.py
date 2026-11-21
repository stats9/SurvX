# Builds the LaTeX tables of the manuscript from the CSV files written by
# paper2/sim/summarize.R (copied to tables/). Run from paper2/manuscript:
#   python3 make_tables.py
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
T = lambda f: os.path.join(HERE, "tables", f)


def read(f):
    with open(T(f)) as h:
        return list(csv.DictReader(h))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def fb(x):   # bias, 3 decimals with sign; -0.000 shown as 0.000
    if x != x:
        return "--"
    if abs(x) < 0.0005:
        x = 0.0
    return f"$\\phantom{{-}}{x:.3f}$" if x >= 0 else f"${x:.3f}$"


def f2(x):
    return "--" if x != x else f"{x:.2f}"


def cp(x):   # coverage in percent, bold if outside 93.6-96.4 (95% +- 1.96 Monte Carlo SE, 1000 replications)
    if x != x:
        return "--"
    s = f"{100 * x:.1f}"
    return f"\\textbf{{{s}}}" if (x < 0.936 or x > 0.964) else s


# ---------------------------------------------------------------------------------------
# Table 1: Study A
# ---------------------------------------------------------------------------------------
A = read("tableA_full.csv")
idx = {(r["n"], r["alpha"], r["cmax"], r["estimator"], r["estimand"]): r for r in A}
est_lab = [("weighted (IPCW variance)", "Weighted"),
           ("weighted GOR (alpha estimated)", "Weighted, GOR"),
           ("joint (cause 2 misspecified)", "Joint"),
           ("Fine-Gray (semi-parametric)", "Fine--Gray")]
estimands = ["b_z1", "Finf_0", "Ft0_0"]
lines = []
for alpha in ["0", "1"]:
    for cmax in ["10", "4"]:
        for n in ["300", "1000"]:
            first = True
            for e, lab in est_lab:
                if e.startswith("Fine") and alpha != "0":
                    continue
                cells = []
                for m in estimands:
                    r = idx.get((n, alpha, cmax, e, m))
                    if r is None:
                        cells += ["--", "--", "--"]
                    else:
                        cells += [fb(num(r["bias"])), f2(num(r["se_sd"])), cp(num(r["coverage"]))]
                link = "PH" if alpha == "0" else "PO"
                fu = "long" if cmax == "10" else "short"
                head = f"{link} & {fu} & {n}" if first else " & & "
                first = False
                lines.append(f"{head} & {lab} & " + " & ".join(cells) + r" \\")
            lines.append(r"\addlinespace[2pt]")
lines = lines[:-1]
with open(os.path.join(HERE, "tab_studyA.tex"), "w") as h:
    h.write("\n".join(lines) + "\n\\bottomrule\n")

# Supplementary: every estimand of study A, including the variance ignoring the censoring term
sup = []
allest = ["b_z1", "b_z2", "Finf_0", "Finf_1", "Ft0_0", "Ft0_1"]
elab = {"b_z1": r"$\beta_1$", "b_z2": r"$\beta_2$", "Finf_0": r"$F_1(\infty\mid 0)$",
        "Finf_1": r"$F_1(\infty\mid 1)$", "Ft0_0": r"$F_1(t_0\mid 0)$", "Ft0_1": r"$F_1(t_0\mid 1)$"}
ests = ["weighted (IPCW variance)", "weighted (censoring ignored)", "weighted GOR (alpha estimated)",
        "weighted GOR (censoring ignored)", "joint (cause 2 misspecified)", "Fine-Gray (semi-parametric)"]
slab = {ests[0]: "W", ests[1]: "W$_0$", ests[2]: "WG", ests[3]: "WG$_0$", ests[4]: "J", ests[5]: "FG"}
for alpha in ["0", "1"]:
    for cmax in ["10", "4"]:
        for n in ["300", "1000"]:
            for m in allest:
                row = []
                for e in ests:
                    r = idx.get((n, alpha, cmax, e, m))
                    row.append(f"{f2(num(r['se_sd']))}/{cp(num(r['coverage']))}" if r else "--")
                link = "PH" if alpha == "0" else "PO"
                fu = "long" if cmax == "10" else "short"
                sup.append(f"{link} & {fu} & {n} & {elab[m]} & " + " & ".join(row) + r" \\")
            sup.append(r"\addlinespace[2pt]")
with open(os.path.join(HERE, "tab_studyA_full.tex"), "w") as h:
    h.write("\n".join(sup[:-1]) + "\n\\bottomrule\n")

# ---------------------------------------------------------------------------------------
# Table 2: Study B, long-run risks at z = 0, and coherence
# ---------------------------------------------------------------------------------------
B = read("tableB_full.csv")
C = read("tableB_coherence_full.csv")
bidx = {(r["n"], r["dgp"], r["cmax"], r["estimator"], r["estimand"]): r for r in B}
cidx = {(r["n"], r["dgp"], r["cmax"], r["estimator"]): r for r in C}
lines = []
for dgp, dlab in [("constrained", "Constrained"), ("csh", "Cause-specific")]:
    for cmax in ["12", "5"]:
        for n in ["500", "1500"]:
            first = True
            for e, lab in [("constrained, K = 3", "Constrained"), ("joint, unconstrained", "Joint")]:
                cells = []
                for k in "123":
                    r = bidx.get((n, dgp, cmax, e, f"F{k}_inf_0"))
                    cells += [fb(num(r["bias"])), cp(num(r["coverage"]))] if r else ["--", "--"]
                c = cidx.get((n, dgp, cmax, e))
                ex = f"{100 * num(c['excess_inf']):.0f}" if c else "--"
                mx = f"{num(c['max_inf']):.2f}" if c else "--"
                fu = "long" if cmax == "12" else "short"
                head = f"{dlab} & {fu} & {n}" if first else " & & "
                first = False
                lines.append(f"{head} & {lab} & " + " & ".join(cells) + f" & {ex} & {mx}" + r" \\")
            lines.append(r"\addlinespace[2pt]")
with open(os.path.join(HERE, "tab_studyB.tex"), "w") as h:
    h.write("\n".join(lines[:-1]) + "\n\\bottomrule\n")

# Supplementary: incidences at t0 = 2, with the cause-specific Cox model
sup = []
for dgp, dlab in [("constrained", "Constrained"), ("csh", "Cause-specific")]:
    for cmax in ["12", "5"]:
        for n in ["500", "1500"]:
            first = True
            for e, lab in [("constrained, K = 3", "Constrained"), ("joint, unconstrained", "Joint"),
                           ("cause-specific Cox", "Cause-specific Cox")]:
                cells = []
                for k in "123":
                    r = bidx.get((n, dgp, cmax, e, f"F{k}_t0_0"))
                    cells += [fb(num(r["bias"])), f"{num(r['rmse']):.3f}", cp(num(r["coverage"]))] if r else ["--"] * 3
                fu = "long" if cmax == "12" else "short"
                head = f"{dlab} & {fu} & {n}" if first else " & & "
                first = False
                sup.append(f"{head} & {lab} & " + " & ".join(cells) + r" \\")
            sup.append(r"\addlinespace[2pt]")
with open(os.path.join(HERE, "tab_studyB_t0.tex"), "w") as h:
    h.write("\n".join(sup[:-1]) + "\n\\bottomrule\n")

# ---------------------------------------------------------------------------------------
# Table 3: Study C, seconds per fit
# ---------------------------------------------------------------------------------------
Cc = read("tableC_full.csv")
cols = ["weighted_PH", "weighted_PO", "joint_PH", "constrained_K3", "fine_gray"]
lines = []
for r in Cc:
    n = int(float(r["n"]))
    lines.append(f"{n:,}".replace(",", "{,}") + " & " + " & ".join(f"{num(r[c]):.3g}" for c in cols) + r" \\")
with open(os.path.join(HERE, "tab_studyC.tex"), "w") as h:
    h.write("\n".join(lines) + "\n\\bottomrule\n")
print("tables written")
