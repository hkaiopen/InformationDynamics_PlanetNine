#!/usr/bin/env python3
"""Validate the Siraj estimator by reproducing Siraj et al. 2025b on their own
Table 1 sample (46 nonresonant TNOs, 80-400 AU), matched into the current MPCORB
by (a, q, i, Omega). Published: 80-200 AU i0~15 deg, node~120 deg, 2.74 sigma;
200-400 AU ~1 sigma; 80-400 AU 2.52 sigma. Also: same 46 objects' tilt in the
FULL current sample context, and the P9/Y predicted bin-averaged tilt at 200-400.
"""
import numpy as np
import pandas as pd
from scipy.special import i0
from scipy.optimize import minimize, minimize_scalar

exec(open("mean_plane_measure.py").read().split("base = df[")[0])

# Siraj Table 1: (a, q, i, Omega) -- 33 in 80-200, 13 in 200-400
T1 = [
    (83.0, 31.6, 14.0, 78.3), (86.7, 33.9, 18.4, 213.9), (87.4, 39.9, 48.8, 153.1),
    (88.5, 37.3, 20.9, 109.8), (88.8, 38.7, 19.7, 303.4), (88.8, 37.9, 25.2, 33.6),
    (89.5, 34.3, 12.8, 176.8), (90.1, 33.2, 17.2, 192.5), (91.7, 32.5, 14.1, 125.5),
    (95.5, 35.9, 4.6, 135.1), (97.3, 37.9, 24.8, 107.4), (97.9, 33.5, 10.4, 218.2),
    (98.1, 40.0, 21.3, 212.6), (98.4, 39.2, 23.4, 327.1), (100.6, 34.2, 13.9, 141.0),
    (100.9, 36.0, 4.8, 60.6), (101.3, 31.7, 18.1, 358.1), (101.5, 31.0, 25.5, 341.6),
    (104.3, 39.1, 25.4, 239.6), (107.1, 32.9, 16.7, 157.4), (108.6, 37.3, 30.5, 69.6),
    (112.1, 30.1, 19.0, 181.2), (116.8, 35.3, 14.9, 94.5), (120.8, 40.0, 31.7, 175.7),
    (121.4, 37.7, 15.5, 96.0), (122.0, 33.6, 18.7, 161.1), (126.4, 33.5, 19.8, 192.4),
    (138.7, 38.1, 40.8, 60.8), (143.3, 35.1, 8.6, 185.2), (159.5, 38.1, 15.5, 197.9),
    (164.8, 33.2, 23.4, 33.6), (192.0, 35.2, 12.0, 24.4), (199.4, 34.5, 12.8, 206.2),
    (210.7, 36.0, 14.9, 201.6), (215.5, 34.3, 30.8, 179.4), (225.1, 31.3, 6.0, 341.5),
    (233.0, 35.1, 13.2, 119.0), (252.0, 65.9, 11.0, 72.1), (262.3, 80.5, 24.1, 90.8),
    (309.4, 34.0, 12.9, 273.1), (311.2, 38.4, 8.8, 46.1), (314.3, 38.1, 6.5, 94.7),
    (335.3, 35.0, 23.5, 192.3), (358.1, 32.2, 15.8, 294.4), (358.9, 36.1, 29.5, 67.6),
    (372.5, 35.5, 42.9, 149.2),
]
base = df[(df["q"] > 30) & (df["U"] <= 6)].copy()
cand = base[(base["a"] >= 80) & (base["a"] <= 400)].copy()
matched_idx, bad = [], []
for (a, q, i, Om) in T1:
    d = np.sqrt(((cand["a"] - a) / 2.0) ** 2 + ((cand["q"] - q) / 1.0) ** 2
                + ((cand["i"] - i) / 1.0) ** 2 + ((cand["Omega"] - Om) / 3.0) ** 2)
    j = int(np.argmin(d.values))
    if d.values[j] > 5:
        bad.append((a, q, i, Om, float(d.values[j])))
    matched_idx.append(cand.index[j])
print(f"matched {len(matched_idx)} / {len(T1)} Table-1 objects; poor matches: {bad}")
sub46 = base.loc[sorted(set(matched_idx))].copy()
print(f"unique matched: {len(sub46)}; a range {sub46['a'].min():.1f}-{sub46['a'].max():.1f}")

analyze("Table-1 subset 80-200 AU", sub46[sub46["a"] < 200])
analyze("Table-1 subset 200-400 AU", sub46[sub46["a"] >= 200])
analyze("Table-1 subset 80-400 AU", sub46)

# Predicted bin-averaged tilt at 200-400 AU from the bodies model (linear forced plane)
exec(open("bodies_coupling_fit.py").read().split('print("\\n--- giants only')[0])
P9_BB21 = ("P9 BB21", 6.2, 380.0, 16.0, 100.0)
P9_SIRAJ = ("P9 Siraj", 4.4, 290.0, 7.0, 100.0)
Y_LIGHT = ("Planet Y light", 0.06, 150.0, 20.0, 120.0)
Y_MID = ("Planet Y mid", 0.1, 150.0, 25.0, 120.0)


def bin_tilt(extra, a_values):
    ts, ns = [], []
    for a in a_values:
        t, n = forced_plane(a, extra=list(extra))
        ts.append(t)
        ns.append(n)
    return float(np.mean(ts)), float(np.mean(ns))


a200_400 = sub46[sub46["a"] >= 200]["a"].values
a200_400_all = base[(base["a"] >= 200) & (base["a"] <= 400)]["a"].values
print("\n--- predicted mean tilt in the 200-400 AU bin (linear forced plane) ---")
for label, extra in [("giants only", []), ("giants + P9 BB21", [P9_BB21]),
                     ("giants + P9 Siraj (i=7)", [P9_SIRAJ]),
                     ("giants + Planet Y light", [Y_LIGHT]),
                     ("giants + Planet Y mid", [Y_MID]),
                     ("giants + P9 BB21 + Planet Y mid", [P9_BB21, Y_MID])]:
    t1, n1 = bin_tilt(extra, a200_400)
    t2, n2 = bin_tilt(extra, a200_400_all)
    print(f"{label}: Table-1 a-values tilt={t1:.1f} deg node={n1:.0f}; "
          f"current-sample a-values tilt={t2:.1f} deg node={n2:.0f}")
print("DONE")
