#!/usr/bin/env python3
"""Virtual-real-space-corrected test for the Planet Nine sample.

The paper tested omega (argument of perihelion) on a circle with a
two-window (N=2) test. From the ID virtual-real viewpoint the virtual
space is a 2-sphere, so the natural projected variables are DIRECTIONS
on a sphere, not a 1-D phase angle inside the orbital plane:
  - perihelion direction P (unit vector, ecliptic coords)
  - orbit pole h (angular-momentum direction, unit vector)
The spherical analogues are:
  - von Mises-Fisher distribution (not von Mises on a circle)
  - Rayleigh test on the sphere (unimodal) and the orientation-tensor
    eigenvalue test (axial / bipolar = the spherical N=2 analogue).
Data: etno_table_nov2026.csv (n=63, paper criteria, current MPCORB).
"""
import numpy as np
import pandas as pd
from scipy.optimize import brentq

t = pd.read_csv("etno_table_nov2026.csv")
om = np.radians(t["omega"].values)
Om = np.radians(t["Omega"].values)
inc = np.radians(t["i"].values)

# Perihelion direction (ecliptic Cartesian, unit vector)
P = np.column_stack([
    np.cos(om) * np.cos(Om) - np.sin(om) * np.sin(Om) * np.cos(inc),
    np.cos(om) * np.sin(Om) + np.sin(om) * np.cos(Om) * np.cos(inc),
    np.sin(om) * np.sin(inc),
])
# Orbit pole (angular momentum direction)
H = np.column_stack([
    np.sin(inc) * np.sin(Om),
    -np.sin(inc) * np.cos(Om),
    np.cos(inc),
])

def rayleigh_sphere(X):
    n = len(X)
    Rvec = X.mean(axis=0)
    R = np.linalg.norm(Rvec)
    return R, Rvec

def axial_eig(X):
    n = len(X)
    S = (X[:, :, None] * X[:, None, :]).mean(axis=0)  # scatter matrix
    w = np.linalg.eigvalsh(S)[::-1]  # descending
    return w  # uniform expectation: all 1/3

def fisher_kappa(R, n):
    # MLE: coth(k) - 1/k = R
    if R < 1e-6:
        return 0.0
    f = lambda k: 1 / np.tanh(k) - 1 / k - R if k > 0 else -R
    try:
        return brentq(f, 1e-6, 500)
    except Exception:
        return float("nan")

def mc_pvalue(X, n_sim=20000, seed=99):
    n = len(X)
    R_obs, _ = rayleigh_sphere(X)
    w_obs = axial_eig(X)
    rng = np.random.default_rng(seed)
    cnt_R = 0
    cnt_ax = 0
    for _ in range(n_sim):
        # uniform random directions on sphere
        u = rng.normal(size=(n, 3))
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        R_s, _ = rayleigh_sphere(u)
        w_s = axial_eig(u)
        if R_s >= R_obs:
            cnt_R += 1
        if w_s[0] >= w_obs[0]:
            cnt_ax += 1
    return R_obs, w_obs, cnt_R / n_sim, cnt_ax / n_sim

def ecl_lon_lat(v):
    lon = np.degrees(np.arctan2(v[1], v[0])) % 360
    lat = np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
    return lon, lat

for label, X in [("Perihelion directions P", P), ("Orbit poles h", H)]:
    n = len(X)
    R, Rvec = rayleigh_sphere(X)
    w = axial_eig(X)
    k = fisher_kappa(R, n)
    lon, lat = ecl_lon_lat(Rvec / R) if R > 0 else (float("nan"), float("nan"))
    R_obs, w_obs, p_R, p_ax = mc_pvalue(X)
    print(f"\n===== {label} (n={n}) =====")
    print(f"Rayleigh resultant R = {R:.3f}  (uniform: E[R]~{1/np.sqrt(n):.3f})")
    print(f"  Fisher kappa (MLE) = {k:.2f}")
    print(f"  mean direction: ecliptic lon {lon:.1f} deg, lat {lat:.1f} deg")
    print(f"  Monte Carlo p (unimodal, 20000 trials) = {p_R:.4f}")
    print(f"Orientation-tensor eigenvalues: {w[0]:.3f}, {w[1]:.3f}, {w[2]:.3f}  (uniform: 1/3 each)")
    print(f"  Monte Carlo p (axial/bipolar, largest eig) = {p_ax:.4f}")

# Robustness: e<0.9 subset (n=37)
sub = t[t["e"] < 0.9]
idx = sub.index
print(f"\n--- e<0.9 subset (n={len(sub)}) ---")
for label, X in [("Perihelion directions P", P[[t.index.get_loc(i) for i in idx]]),
                 ("Orbit poles h", H[[t.index.get_loc(i) for i in idx]])]:
    R, Rvec = rayleigh_sphere(X)
    w = axial_eig(X)
    _, _, p_R, p_ax = mc_pvalue(X, n_sim=10000)
    print(f"{label}: R={R:.3f}, Fisher kappa={fisher_kappa(R, len(X)):.2f}, "
          f"eig={w[0]:.3f}/{w[1]:.3f}/{w[2]:.3f}, p_Rayleigh={p_R:.4f}, p_axial={p_ax:.4f}")
print("\nDONE")
