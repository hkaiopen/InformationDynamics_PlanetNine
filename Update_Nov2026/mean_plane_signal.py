#!/usr/bin/env python3
"""Measure the strongest astronomer-defined signals from current MPCORB.

1. Distant-KBO mean-plane offset (Siraj et al. 2025b, MNRAS 543, L27:
   warp 80-400 AU at 98%, 80-200 AU at 96%, 154 non-resonant objects,
   "Planet Y" 0.06-1 M_E, 100-200 AU). Here: proxy samples (no resonance
   classification -- stated caveat), mean pole vs invariable pole with
   bootstrap CI. A real warp has a consistent offset DIRECTION across
   independent subsamples, not just a nonzero angle.
2. Neptune-crosser raw fraction (Batygin et al. 2024, ApJL 966, L8:
   P9-free model rejected ~5sigma after bias modelling; P9 predicts 3%
   vs 0.5% without). Raw counts only -- the 5sigma test itself needs
   their bias model + simulations and is NOT reproduced here.
3. Sednoid list (q>60, a>200): the detached objects themselves.
"""
import numpy as np
import pandas as pd

colspecs = [(0, 7), (8, 13), (14, 19), (20, 25), (26, 35), (37, 46),
            (48, 57), (59, 68), (70, 79), (80, 91), (92, 103)]
names = ["desig", "H", "G", "epoch", "M", "omega", "Omega", "i", "e", "n", "a"]
df = pd.read_fwf("MPCORB-202611.DAT", colspecs=colspecs, names=names,
                 header=None, skiprows=1)
for col in ["a", "e", "omega", "Omega", "i", "H"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df.dropna(subset=["a", "e", "omega", "Omega", "i"], inplace=True)
df["q"] = df["a"] * (1 - df["e"])
print(f"parsed {len(df)} objects")

INV = np.array([np.cos(np.radians(88.42)) * np.cos(np.radians(17.58)),
                np.cos(np.radians(88.42)) * np.sin(np.radians(17.58)),
                np.sin(np.radians(88.42))])  # invariable-plane pole
ECL = np.array([0.0, 0.0, 1.0])


def poles(sub):
    Om = np.radians(sub["Omega"].values)
    inc = np.radians(sub["i"].values)
    return np.column_stack([np.sin(inc) * np.sin(Om),
                            -np.sin(inc) * np.cos(Om),
                            np.cos(inc)])


def mean_pole_stats(label, sub, n_boot=5000, seed=3):
    X = poles(sub)
    n = len(X)
    m = X.mean(axis=0)
    R = float(np.linalg.norm(m))
    mu = m / R
    lon = float(np.degrees(np.arctan2(mu[1], mu[0])) % 360)
    lat = float(np.degrees(np.arcsin(np.clip(mu[2], -1, 1))))
    off_inv = float(np.degrees(np.arccos(np.clip(mu @ INV, -1, 1))))
    off_ecl = float(np.degrees(np.arccos(np.clip(mu @ ECL, -1, 1))))
    rng = np.random.default_rng(seed)
    offs, lons = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        mb = X[idx].mean(axis=0)
        Rb = float(np.linalg.norm(mb))
        if Rb < 1e-9:
            continue
        mub = mb / Rb
        offs.append(float(np.degrees(np.arccos(np.clip(mub @ INV, -1, 1)))))
        lons.append(float(np.degrees(np.arctan2(mub[1], mub[0])) % 360))
    offs = np.array(offs)
    lo, hi = np.percentile(offs, [2.5, 97.5])
    print(f"\n--- {label} (n={n}) ---")
    print(f"mean pole: lon {lon:.1f} deg, lat {lat:.1f} deg; R = {R:.3f}")
    print(f"offset from invariable pole: {off_inv:.2f} deg "
          f"[bootstrap 95% CI {lo:.2f}-{hi:.2f}]; offset from ecliptic north: {off_ecl:.2f} deg")
    print(f"bootstrap: P(offset < 1 deg) = {np.mean(offs < 1):.3f}")
    return dict(n=n, lon=lon, lat=lat, off_inv=off_inv, ci=(lo, hi))


mean_pole_stats("ETNO a>200, q>30, H>5", df[(df["a"] > 200) & (df["q"] > 30) & (df["H"] > 5)])
mean_pole_stats("Distant 80-400 AU, q>30 (proxy; resonants included)",
                df[(df["a"] >= 80) & (df["a"] <= 400) & (df["q"] > 30)])
mean_pole_stats("Distant 80-200 AU, q>30 (proxy)",
                df[(df["a"] >= 80) & (df["a"] <= 200) & (df["q"] > 30)])
mean_pole_stats("Distant 200-400 AU, q>30 (proxy)",
                df[(df["a"] > 200) & (df["a"] <= 400) & (df["q"] > 30)])
mean_pole_stats("Distant 80-400 AU, q>30, H>5 (proxy)",
                df[(df["a"] >= 80) & (df["a"] <= 400) & (df["q"] > 30) & (df["H"] > 5)])

print("\n--- Sednoids (q>60, a>200) ---")
sed = df[(df["q"] > 60) & (df["a"] > 200)].sort_values("q", ascending=False)
for _, r in sed.iterrows():
    varpi = (r["omega"] + r["Omega"]) % 360
    print(f"{str(r['desig']).strip():>10} a={r['a']:7.1f} q={r['q']:5.1f} "
          f"i={r['i']:5.1f} omega={r['omega']:6.1f} Omega={r['Omega']:6.1f} varpi={varpi:6.1f}")

print("\n--- Neptune-crossers, raw counts (Batygin et al. 2024 population: a>100, i<40) ---")
base = df[(df["a"] > 100) & (df["i"] < 40)]
nc = base[base["q"] < 30]
print(f"a>100, i<40: n={len(base)}, q<30: n={len(nc)} ({100*len(nc)/len(base):.1f}% raw; "
      f"paper's bias-corrected expectation 3% with P9 vs 0.5% without -- raw count is NOT the test)")
print("DONE")
