#!/usr/bin/env python3
"""Extended spherical-geometry test: virtual-space S^2 applied to Planet Nine data.

Variables (directions on a sphere, not a 1-D phase angle):
  P = perihelion direction (unit vector, ecliptic coords)
  h = orbit pole (angular-momentum direction)
Statistics: Rayleigh on the sphere, von Mises-Fisher kappa (MLE),
orientation-tensor eigenvalues (bipolar: lam1 >> lam2 ~= lam3;
girdle/planar: lam1 ~= lam2 >> lam3), Monte Carlo p vs uniform,
bootstrap kappa CI, split-half mean-direction stability.
Samples from MPCORB-202611.DAT (current, 2026-11-16).
"""
import numpy as np
import pandas as pd
from scipy.optimize import brentq

colspecs = [(0, 7), (8, 13), (14, 19), (20, 25), (26, 35), (37, 46),
            (48, 57), (59, 68), (70, 79), (80, 91), (92, 103)]
names = ["desig", "H", "G", "epoch", "M", "omega", "Omega", "i", "e", "n", "a"]
df = pd.read_fwf("MPCORB-202611.DAT", colspecs=colspecs, names=names,
                 header=None, skiprows=1)
for col in ["a", "e", "omega", "Omega", "i", "H"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df.dropna(subset=["a", "e", "omega", "Omega", "i"], inplace=True)
df["q"] = df["a"] * (1 - df["e"])
print(f"parsed {len(df)} objects", flush=True)


def dirs(sub):
    om = np.radians(sub["omega"].values % 360)
    Om = np.radians(sub["Omega"].values % 360)
    inc = np.radians(sub["i"].values)
    P = np.column_stack([
        np.cos(om) * np.cos(Om) - np.sin(om) * np.sin(Om) * np.cos(inc),
        np.cos(om) * np.sin(Om) + np.sin(om) * np.cos(Om) * np.cos(inc),
        np.sin(om) * np.sin(inc),
    ])
    H = np.column_stack([
        np.sin(inc) * np.sin(Om),
        -np.sin(inc) * np.cos(Om),
        np.cos(inc),
    ])
    return P, H


def kappa_from_R(R):
    if R < 1e-9:
        return 0.0
    if R > 0.9:
        f = lambda k: 1 / np.tanh(k) - 1 / k - R
        return float(brentq(f, 1e-6, 500))
    return float(R * (3 - R * R) / (1 - R * R))


def analyse(label, X, n_mc=5000, n_boot=500, seed=11):
    n = len(X)
    Rvec = X.mean(axis=0)
    R = float(np.linalg.norm(Rvec))
    k = kappa_from_R(R)
    S = (X[:, :, None] * X[:, None, :]).mean(axis=0)
    w = np.linalg.eigvalsh(S)[::-1]
    rng = np.random.default_rng(seed)
    # Monte Carlo: uniform directions
    cnt_R = cnt_ax = 0
    for _ in range(n_mc):
        u = rng.normal(size=(n, 3))
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        if np.linalg.norm(u.mean(axis=0)) >= R:
            cnt_R += 1
        wu = np.linalg.eigvalsh((u[:, :, None] * u[:, None, :]).mean(axis=0))[::-1]
        if wu[0] - wu[1] >= w[0] - w[1]:
            cnt_ax += 1
    # bootstrap kappa CI (resample resultant lengths, convert)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        Rb = float(np.linalg.norm(X[idx].mean(axis=0)))
        boots.append(kappa_from_R(Rb))
    klo, khi = np.percentile(boots, [16, 84])
    mean_dir = Rvec / R if R > 0 else np.array([np.nan] * 3)
    lon = float(np.degrees(np.arctan2(mean_dir[1], mean_dir[0])) % 360)
    lat = float(np.degrees(np.arcsin(np.clip(mean_dir[2], -1, 1))))
    print(f"\n--- {label} (n={n}) ---")
    print(f"Rayleigh R = {R:.3f} (uniform E[R] ~ {1/np.sqrt(n):.3f}); "
          f"Fisher kappa = {k:.2f}, bootstrap 68% CI [{klo:.2f}, {khi:.2f}]")
    print(f"mean direction: lon {lon:.1f} deg, lat {lat:.1f} deg")
    print(f"MC p (Rayleigh) = {cnt_R/n_mc:.4f}")
    print(f"orientation eig = {w[0]:.3f} / {w[1]:.3f} / {w[2]:.3f}; "
          f"lam1-lam2 = {w[0]-w[1]:.3f}; MC p (bipolar, lam1-lam2) = {cnt_ax/n_mc:.4f}")
    return dict(R=R, kappa=k, klo=klo, khi=khi, lon=lon, lat=lat,
                eig=w, mean_dir=mean_dir)


def angle_between(u, v):
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v), -1, 1))))


samples = {
    "MAIN a>200, q>30, H>5": df[(df["a"] > 200) & (df["q"] > 30) & (df["H"] > 5)],
    "ALT a>150, q>30, H>5": df[(df["a"] > 150) & (df["q"] > 30) & (df["H"] > 5)],
    "ALT a>250, q>30, H>5": df[(df["a"] > 250) & (df["q"] > 30) & (df["H"] > 5)],
}
for label, sub in samples.items():
    P, H = dirs(sub)
    print(f"\n===== {label}: n={len(sub)} =====", flush=True)
    analyse("Perihelion directions P", P)
    analyse("Orbit poles h", H, n_mc=2000, n_boot=200)

# Split-half stability on the main sample (a real signal's mean direction
# should agree between independent halves; under uniformity the angle
# between half-means is random with mean ~90 deg).
main = samples["MAIN a>200, q>30, H>5"]
P, _ = dirs(main)
splits = {
    "a: 200-350 vs >350": (main["a"] <= 350, main["a"] > 350),
    "e: <=0.9 vs >0.9": (main["e"] <= 0.9, main["e"] > 0.9),
    "q: <=40 vs >40": (main["q"] <= 40, main["q"] > 40),
}
print("\n===== Split-half stability of mean perihelion direction =====")
for name, (m1, m2) in splits.items():
    X1, X2 = P[m1.values], P[m2.values]
    r1 = X1.mean(axis=0); r1 /= np.linalg.norm(r1)
    r2 = X2.mean(axis=0); r2 /= np.linalg.norm(r2)
    print(f"{name}: n={len(X1)}/{len(X2)}, R1={np.linalg.norm(X1.mean(axis=0)):.3f}, "
          f"R2={np.linalg.norm(X2.mean(axis=0)):.3f}, "
          f"angle between half-means = {angle_between(r1, r2):.1f} deg")
print("\nDONE")
