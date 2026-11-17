#!/usr/bin/env python3
"""xi/zeta test for the carrier (task 2, 2026-11-17).

Reproduces Batygin et al. 2024 (ApJL 966, L8) Sec. 3.2 bias correction
for a given model q-distribution (here: the carrier's secular footprints).

Definitions (from the paper, verified against the text):
  xi_j = CDF_{r_j}(q_j), where r_j is object j's discovery heliocentric
  distance and CDF_{r_j} is the CDF of the model's perihelion distribution
  for discovery distance r_j. The model PDF for r_j is built from the
  simulated footprints (a>100 AU, q<30 AU, i<40 deg), each weighted by
  the fraction of its orbital period spent in the vicinity of r_j
  (geometric correction; Brown 2001; Morbidelli et al. 2004).
  If the model matches the observations, xi_j ~ Uniform(0,1); the paper
  reports a KS test: P9 model p = 0.41, P9-free p = 0.0034.

  zeta = sum_j log10(xi_j)   [see note below on the printed formula]
  The paper's typeset equation renders as a PRODUCT (Pi), but the quoted
  moments (mean -7.2, sigma 1.8 for N=17) are only consistent with a SUM
  of base-10 logs: E[log10 U] = -1/ln10 = -0.4343, so the sum of 17 has
  mean -7.38 and sigma sqrt(17 * Var[log10 U]) = sqrt(17*0.1886) = 1.79.
  A product would give mean -1 and sigma ~362. We therefore use the SUM
  of log10, verify it by Monte Carlo, and report z-scores against both
  our Monte-Carlo null moments and the paper's quoted (-7.2, 1.8).
  The paper's values: zeta(P9) = -7.9, zeta(P9-free) = -16.5 (~5 sigma).

Weight function w(a, q; r_j):
  fraction of the orbital period spent with heliocentric distance in
  [r_j - 1, r_j + 1] AU (1 AU window; the width cancels in the
  normalisation for a given r_j). Computed exactly via Kepler's equation:
  for r in [q, Q], cos E = (1 - r/a)/e; the time from perihelion to E is
  (E - e sin E)/n; the fraction in [r1, r2] is
  [(E2 - e sin E2) - (E1 - e sin E1)] / pi   (both branches).
  If no footprint orbit reaches r_j (total weight 0), xi_j is undefined:
  the model predicts zero objects discoverable at r_j. We assign xi_j = 0
  (the model assigns zero probability to that observation) and flag it.
  If q_j is below every simulated q, CDF(q_j) = 0, so xi_j = 0.

Inputs:
  --results FILE : carrier snapshot results JSON (crosser_samepop_results.json
                   or crosser_distributions_results.json). Uses the
                   snapshot (a, q, i) arrays, filtered to i < 40 deg.
  --key KEY      : which run inside the file (e.g. A_100_400_S_9).
  --objects FILE : xi_objects.json (17 objects with q and r_disc).
  --label TEXT   : label for the output.
"""
import numpy as np
from scipy import stats
import json, sys, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True)
ap.add_argument("--key", required=True)
ap.add_argument("--objects", default="xi_objects.json")
ap.add_argument("--label", default="")
args = ap.parse_args()

with open(args.results) as f:
    res = json.load(f)
with open(args.objects) as f:
    objs = json.load(f)

# ---- get snapshot arrays for the chosen run ----
if "snap_a" in res.get(args.key, {}):
    d = res[args.key]
    snap_a = np.array(d["snap_a"]); snap_q = np.array(d["snap_q"])
    snap_i = np.array(d["snap_i"])
else:
    # crosser_distributions_results.json layout: top-level keys are run tags
    d = res[args.key]
    snap_a = np.array(d["snap_a"]); snap_q = np.array(d["snap_q"])
    snap_i = np.array(d["snap_i"])
mask = (snap_a > 100) & (snap_q < 30) & (snap_i < 40)
sa, sq = snap_a[mask], snap_q[mask]
print(f"[{args.label}] run {args.key}: snapshots total {len(snap_a)}, "
      f"after cuts (a>100, q<30, i<40): {len(sa)}", flush=True)
print(f"  snapshot q range [{sq.min():.2f}, {sq.max():.2f}] AU, "
      f"median {np.median(sq):.2f}", flush=True)


def period_fraction_in_window(a, q, rlo, rhi):
    """Fraction of orbital period with heliocentric distance in [rlo, rhi].
    Vectorised over arrays a, q. Exact via Kepler's equation:
    for r in [q, Q], cos E = (1 - r/a)/e; time from perihelion to E is
    (E - e sin E)/n; the fraction in [r1, r2] (both branches) is
    [(E2 - e sin E2) - (E1 - e sin E1)] / pi."""
    a = np.asarray(a, dtype=float); q = np.asarray(q, dtype=float)
    e = 1.0 - q / a
    Q = a * (1.0 + e)
    lo = np.maximum(rlo, q)
    hi = np.minimum(rhi, Q)
    w = np.zeros_like(a)
    ok = hi > lo
    if np.any(ok):
        ao, eo = a[ok], e[ok]
        def E_of_r(r, ao_, eo_):
            c = (1.0 - r / ao_) / eo_
            return np.arccos(np.clip(c, -1.0, 1.0))
        E1 = E_of_r(lo[ok], ao, eo); E2 = E_of_r(hi[ok], ao, eo)
        w[ok] = ((E2 - eo * np.sin(E2)) - (E1 - eo * np.sin(E1))) / np.pi
    return w


def xi_for_object(qj, rj, sa, sq, halfwin=1.0):
    """xi_j = weighted CDF of model q at q_j, for discovery distance r_j."""
    w = period_fraction_in_window(sa, sq, rj - halfwin, rj + halfwin)
    tot = w.sum()
    if tot <= 0:
        return 0.0, 0.0, 0  # no model orbit reaches r_j
    order = np.argsort(sq)
    qs = sq[order]; ws = w[order]
    cw = np.cumsum(ws) / tot
    # CDF at qj: weight of snapshots with q <= qj
    idx = np.searchsorted(qs, qj, side="right")
    xi = float(cw[idx - 1]) if idx > 0 else 0.0
    n_reach = int((w > 0).sum())
    return xi, float(tot), n_reach


print(f"\n{'desig':10s} {'q':>6s} {'r_disc':>7s} {'xi':>7s} {'w_tot':>8s} {'n_reach':>8s}", flush=True)
xis = []
for o in objs:
    xi, tot, n_reach = xi_for_object(o["q"], o["r_disc"], sa, sq)
    xis.append(xi)
    print(f"{o['desig']:10s} {o['q']:6.2f} {o['r_disc']:7.2f} {xi:7.4f} {tot:8.1f} {n_reach:8d}",
          flush=True)
xis = np.array(xis)
n_zero = int((xis <= 0).sum())
print(f"\n[{args.label}] xi values: {np.round(xis, 4).tolist()}", flush=True)
print(f"  xi == 0 : {n_zero} objects (model assigns zero probability)", flush=True)
print(f"  xi median {np.median(xis):.4f}, mean {xis.mean():.4f}", flush=True)

ks = stats.kstest(xis, "uniform")
print(f"  KS test for uniformity of xi: D = {ks.statistic:.4f}, p = {ks.pvalue:.4g}",
      flush=True)

# ---- zeta: sum of log10(xi), with a floor for xi == 0 ----
# Monte-Carlo null for zeta = sum of 17 log10(U), U~Uniform(0,1)
rng = np.random.default_rng(0)
mc = np.sum(np.log10(rng.uniform(0, 1, size=(200000, 17))), axis=1)
mc_mean, mc_sd = float(mc.mean()), float(mc.std())
print(f"  Monte-Carlo null for zeta (sum of 17 log10 U): "
      f"mean {mc_mean:.3f}, sigma {mc_sd:.3f} "
      f"(paper quotes -7.2, 1.8; analytic -7.383, 1.790)", flush=True)

for floor in [None, 1e-3, 1e-4]:
    if floor is None:
        if n_zero > 0:
            print(f"  zeta (no floor): -inf  ({n_zero} xi == 0)", flush=True)
            continue
        z = float(np.sum(np.log10(xis)))
    else:
        xf = np.maximum(xis, floor)
        z = float(np.sum(np.log10(xf)))
    z_mc = (z - mc_mean) / mc_sd
    z_paper = (z - (-7.2)) / 1.8
    print(f"  zeta (floor={floor}): {z:.2f}; "
          f"z vs MC null {z_mc:+.2f} sigma; z vs paper moments {z_paper:+.2f} sigma",
          flush=True)
print(f"  For comparison: paper zeta(P9) = -7.9, zeta(P9-free) = -16.5 "
      f"(mean -7.2, sigma 1.8).", flush=True)

out = dict(label=args.label, key=args.key, n_snap=int(len(sa)),
           xi=xis.tolist(), n_zero=n_zero,
           ks_D=float(ks.statistic), ks_p=float(ks.pvalue),
           mc_mean=mc_mean, mc_sd=mc_sd)
with open(f"xi_result_{args.key}.json", "w") as f:
    json.dump(out, f, indent=1)
print(f"saved xi_result_{args.key}.json", flush=True)
