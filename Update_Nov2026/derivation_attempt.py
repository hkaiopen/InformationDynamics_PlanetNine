#!/usr/bin/env python3
"""Derivation attempt: virtual-space S^2 eigenmodes -> ETNO directions.

Step 1 (formal): on S^2, eigenmodes are spherical harmonics Y_lm.
  - l=1 dipole: Born-projection map p ~ (1 + r.n)/2  -> Fisher distribution
    f(x) ~ exp(kappa * mu.x), axis mu = virtual Bloch vector r-hat.
  - l=2 quadrupole: antipodal (N=2) symmetry -> Watson axial distribution
    f(x) ~ exp(kappa * (mu.x)^2), axis mu = mode axis.
Step 2 (quantitative): fit both to the perihelion directions (n=63, main
  sample) by exact MLE; likelihood-ratio vs uniform with Monte Carlo p.
Step 3 (axis hypotheses): quadrupole statistic q = <(3(a.x)^2-1)/2> for
  fixed axes a = ecliptic north, galactic north, Planet Nine aphelion
  direction (BB21/BB24); MC p vs uniform.
Step 4 (numerology check, labelled as such): l=2 Laplacian eigenvalue on
  R=2 sphere = l(l+1)/R^2 = 1.5; compare with fitted kappa.
"""
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import hyp1f1

t = pd.read_csv("etno_table_nov2026.csv")
om = np.radians(t["omega"].values % 360)
Om = np.radians(t["Omega"].values % 360)
inc = np.radians(t["i"].values)
P = np.column_stack([
    np.cos(om) * np.cos(Om) - np.sin(om) * np.sin(Om) * np.cos(inc),
    np.cos(om) * np.sin(Om) + np.sin(om) * np.cos(Om) * np.cos(inc),
    np.sin(om) * np.sin(inc),
])
n = len(P)
S = (P[:, :, None] * P[:, None, :]).mean(axis=0)
w, V = np.linalg.eigh(S)  # ascending
lam1 = w[-1]
mu_w = V[:, -1]
print(f"n={n}; scatter eigenvalues {w[0]:.3f}/{w[1]:.3f}/{w[2]:.3f}; "
      f"principal axis lon/lat = "
      f"{np.degrees(np.arctan2(mu_w[1], mu_w[0]))%360:.1f} deg, "
      f"{np.degrees(np.arcsin(mu_w[2])):.1f} deg")

# ---- Watson MLE: (1/3) * 1F1(3/2;5/2;k)/1F1(1/2;3/2;k) = lam1 ----
def watson_ratio(k):
    return (1/3) * hyp1f1(1.5, 2.5, k) / hyp1f1(0.5, 1.5, k)

def watson_kappa(lam):
    if lam <= 1/3 + 1e-9:
        return 0.0
    return float(brentq(lambda k: watson_ratio(k) - lam, 1e-9, 50))

def watson_logC(k):
    return -np.log(4 * np.pi) - np.log(hyp1f1(0.5, 1.5, k))

k_w = watson_kappa(lam1)
proj = P @ mu_w
ll_w = np.sum(k_w * proj**2 + watson_logC(k_w))
ll_u = n * (-np.log(4 * np.pi))
print(f"Watson MLE: kappa = {k_w:.3f}; 2*Delta logL vs uniform = {2*(ll_w-ll_u):.2f}")

rng = np.random.default_rng(5)
cnt = 0
N_MC = 5000
for _ in range(N_MC):
    u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1, keepdims=True)
    Su = (u[:, :, None] * u[:, None, :]).mean(axis=0)
    lam_u = np.linalg.eigvalsh(Su)[-1]
    k_u = watson_kappa(lam_u)
    mu_u = np.linalg.eigh(Su)[1][:, -1]
    ll_u2 = np.sum(k_u * (u @ mu_u)**2 + watson_logC(k_u))
    if 2 * (ll_u2 - n * (-np.log(4 * np.pi))) >= 2 * (ll_w - ll_u):
        cnt += 1
print(f"Watson vs uniform: MC p = {cnt/N_MC:.4f} (5000 trials)")

# ---- Fisher dipole LR (same as Rayleigh, for the record) ----
Rbar = float(np.linalg.norm(P.mean(axis=0)))
print(f"Fisher dipole: Rbar = {Rbar:.3f}; Rayleigh MC p = 0.64 (from directional_test)")

# ---- Fixed-axis quadrupole statistic ----
def eq_to_ecl(ra_deg, dec_deg):
    eps = np.radians(23.4392811)
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    x = np.cos(dec) * np.cos(ra); y = np.cos(dec) * np.sin(ra); z = np.sin(dec)
    xe, ye, ze = x, y * np.cos(eps) - z * np.sin(eps), y * np.sin(eps) + z * np.cos(eps)
    return np.array([xe, ye, ze])

# Galactic north pole (J2000 RA 192.86, Dec +27.13) -> ecliptic
NGP_eq = np.array([np.cos(np.radians(27.13)) * np.cos(np.radians(192.86)),
                   np.cos(np.radians(27.13)) * np.sin(np.radians(192.86)),
                   np.sin(np.radians(27.13))])
eps = np.radians(23.4392811)
NGP_ecl = np.array([NGP_eq[0],
                    NGP_eq[1] * np.cos(eps) - NGP_eq[2] * np.sin(eps),
                    NGP_eq[1] * np.sin(eps) + NGP_eq[2] * np.cos(eps)])
print(f"Galactic north pole, ecliptic lon/lat = "
      f"{np.degrees(np.arctan2(NGP_ecl[1], NGP_ecl[0]))%360:.1f}, "
      f"{np.degrees(np.arcsin(NGP_ecl[2])):.1f} deg")
# Planet Nine aphelion direction, central BB24-ish value RA ~3.5h, Dec +5
P9_ecl = eq_to_ecl(52.5, 5.0)
print(f"P9 aphelion direction (RA 3.5h, Dec +5), ecliptic lon/lat = "
      f"{np.degrees(np.arctan2(P9_ecl[1], P9_ecl[0]))%360:.1f}, "
      f"{np.degrees(np.arcsin(P9_ecl[2])):.1f} deg")

axes = {
    "ecliptic north": np.array([0.0, 0.0, 1.0]),
    "galactic north": NGP_ecl,
    "P9 aphelion dir": P9_ecl,
    "P9 perihelion dir (antipode)": -P9_ecl,
}
for name, a in axes.items():
    a = a / np.linalg.norm(a)
    q = float(np.mean((3 * (P @ a) ** 2 - 1) / 2))
    cntq = 0
    for _ in range(N_MC):
        u = rng.normal(size=(n, 3)); u /= np.linalg.norm(u, axis=1, keepdims=True)
        qu = float(np.mean((3 * (u @ a) ** 2 - 1) / 2))
        if abs(qu) >= abs(q):
            cntq += 1
    print(f"quadrupole q, axis={name}: q = {q:+.3f}; MC p (two-sided) = {cntq/N_MC:.4f}")

# Angle between measured mean perihelion direction and P9 aphelion dir
m = P.mean(axis=0); m /= np.linalg.norm(m)
ang = np.degrees(np.arccos(np.clip(m @ P9_ecl, -1, 1)))
print(f"angle(mean P dir, P9 aphelion dir) = {ang:.1f} deg "
      f"(mean dir itself is noise: Rbar={Rbar:.3f})")
print(f"l=2 eigenvalue on R=2 sphere: l(l+1)/R^2 = {2*3/4:.2f}; "
      f"on R=1: {2*3:.2f}; fitted Watson kappa = {k_w:.3f} (no derivation links them)")
print("DONE")
