#!/usr/bin/env python3
"""Warp-target geometry: Siraj et al. 2025b (MNRAS 543, L27) bias-free warp
vs naive mean pole, natural axes, perturber-node constraint.

Paper numbers (verified from the MNRAS PDF):
  warp i0 ~ 15 deg, node Omega0 ~ 120 deg (large uncertainties);
  80-200 AU: 2.74 sigma; 80-400 AU: 2.52 sigma; 50-80 and 200-400 AU: ~1 sigma (no warp).
  Invariable plane: i_inv = 1.6 deg, Omega_inv = 107 deg.
  Naive orbit-normal averaging is within 1 deg of truth only ~20% of the time
  (their Monte Carlo) -- it severely underestimates the warp.
"""
import numpy as np
import pandas as pd

EPS = np.radians(23.4392811)  # obliquity J2000


def vec_from_i_Omega(i_deg, Om_deg):
    i, Om = np.radians(i_deg), np.radians(Om_deg)
    return np.array([np.sin(i) * np.sin(Om), -np.sin(i) * np.cos(Om), np.cos(i)])


def lonlat(v):
    v = v / np.linalg.norm(v)
    return float(np.degrees(np.arctan2(v[1], v[0])) % 360), float(np.degrees(np.arcsin(np.clip(v[2], -1, 1))))


def angle(a, b):
    return float(np.degrees(np.arccos(np.clip(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)), -1, 1))))


def eq_to_ecl(ra_deg, dec_deg):
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)
    xe = x
    ye = y * np.cos(EPS) - z * np.sin(EPS)
    ze = y * np.sin(EPS) + z * np.cos(EPS)
    return np.array([xe, ye, ze])


INV = vec_from_i_Omega(1.6, 107.0)
WARP = vec_from_i_Omega(15.0, 120.0)
print("Paper warp pole (i0=15, Om0=120): lon %.1f, lat %.1f" % lonlat(WARP))
print("Angle warp pole vs invariable pole: %.1f deg" % angle(WARP, INV))
# p-q check: p0=sin(i0/2) sin Om0, q0=sin(i0/2) cos Om0
p0 = np.sin(np.radians(7.5)) * np.sin(np.radians(120))
q0 = np.sin(np.radians(7.5)) * np.cos(np.radians(120))
pinv = np.sin(np.radians(0.8)) * np.sin(np.radians(107))
qinv = np.sin(np.radians(0.8)) * np.cos(np.radians(107))
print("p-q distance warp vs invariable: %.3f (paper says ~0.1)" % np.hypot(p0 - pinv, q0 - qinv))

# Natural axes
NGP = eq_to_ecl(192.8595, 27.1283)   # north galactic pole J2000
SUN = eq_to_ecl(286.13, 63.87)       # solar spin axis J2000
print("\nNorth galactic pole (ecliptic): lon %.1f, lat %.1f" % lonlat(NGP))
print("Solar spin axis (ecliptic): lon %.1f, lat %.1f" % lonlat(SUN))
for name, ax in [("galactic pole", NGP), ("solar spin axis", SUN)]:
    print(f"angle warp pole vs {name}: {angle(WARP, ax):.1f} deg")
# Is the warp TOWARD the candidate? Compare the offset direction (INV -> WARP)
# with the direction from INV toward the candidate, in the tangent plane at INV.
def tangent_dir(target):
    t = target - (target @ INV) * INV
    return t / np.linalg.norm(t)


off_dir = tangent_dir(WARP)
for name, ax in [("galactic pole", NGP), ("solar spin axis", SUN)]:
    d = tangent_dir(ax)
    print(f"angle between warp-offset direction and direction to {name}: "
          f"{float(np.degrees(np.arccos(np.clip(off_dir @ d, -1, 1)))):.1f} deg "
          f"(0 = warp points straight at it)")

# Naive means on current MPCORB, paper's bins, for the estimator comparison
colspecs = [(0, 7), (8, 13), (14, 19), (20, 25), (26, 35), (37, 46),
            (48, 57), (59, 68), (70, 79), (80, 91), (92, 103)]
names = ["desig", "H", "G", "epoch", "M", "omega", "Omega", "i", "e", "n", "a"]
df = pd.read_fwf("MPCORB-202611.DAT", colspecs=colspecs, names=names, header=None, skiprows=1)
for c in ["a", "e", "omega", "Omega", "i", "H"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df.dropna(subset=["a", "e", "omega", "Omega", "i"], inplace=True)
df["q"] = df["a"] * (1 - df["e"])


def naive_mean(label, sub):
    Om = np.radians(sub["Omega"].values)
    inc = np.radians(sub["i"].values)
    X = np.column_stack([np.sin(inc) * np.sin(Om), -np.sin(inc) * np.cos(Om), np.cos(inc)])
    m = X.mean(axis=0)
    mu = m / np.linalg.norm(m)
    lon, lat = lonlat(mu)
    i0 = 90 - lat
    Om0 = (lon + 90) % 360
    print(f"naive {label}: n={len(sub)} pole lon {lon:.1f} lat {lat:.1f} "
          f"-> (i0={i0:.1f}, Om0={Om0:.1f}); offset from invariable pole {angle(mu, INV):.2f} deg")


naive_mean("50-80 AU, q>30", df[(df["a"] >= 50) & (df["a"] < 80) & (df["q"] > 30)])
naive_mean("80-200 AU, q>30", df[(df["a"] >= 80) & (df["a"] < 200) & (df["q"] > 30)])
naive_mean("200-400 AU, q>30", df[(df["a"] >= 200) & (df["a"] <= 400) & (df["q"] > 30)])
naive_mean("80-400 AU, q>30", df[(df["a"] >= 80) & (df["a"] <= 400) & (df["q"] > 30)])
print("\nPerturber constraint (linear secular theory): forced plane lies between")
print("invariable plane and perturber plane, node forced = node perturber.")
print("Warp node Om0~120 deg -> perturber node Om_Y ~ 120 deg (large uncertainty).")
print("Warp angle from invariable pole ~13.5 deg -> perturber inclination i_Y >~ 13.5 deg")
print("(paper's N-body fits: i_Y >~ 10 deg, mass 0.06-1 M_E, a 100-200 AU).")
print("DONE")
