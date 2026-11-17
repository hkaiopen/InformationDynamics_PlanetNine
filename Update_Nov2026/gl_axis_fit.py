#!/usr/bin/env python3
"""Fit the warp axis with the generalized GL equation (dipole-sector reduction).

Framework Eq. (1): dPsi/dt = D nabla^2 Psi - k(|Psi|^2 - rho_t) Psi + iK(x) Psi + u.
- Diffusion, compression, anti-diffusion, log potential, quartic: all isotropic.
- The ONLY symmetry-breaking slot in the minimal equation is K(x):
  a direction-dependent coupling K(n) = k0 + k1 (n.d) tilts the mode axis to d.
- The mode axis is the dipole of the output density y(n) = |Psi(n)|^2,
  i.e. the mean pole -- exactly the warp observable.

Toy models (free energy F on the unit sphere, m = mode axis):
  (V) vector: F = -A (m.z)^2 - h.m ;  z = invariable pole (base axis, INPUT).
  (A) axial : F = -A (m.z)^2 - B (m.d)^2 ;  d = second axis (plane symmetry respected).
Observed warp (Siraj et al. 2025b, bias-free): pole (lon 30, lat 75),
theta = 13.4 deg from invariable pole, node Om0 ~ 120 deg, 80-200 AU only.
"""
import numpy as np
from scipy.optimize import minimize_scalar

INV = np.array([np.cos(np.radians(88.42)) * np.cos(np.radians(17.58)),
                np.cos(np.radians(88.42)) * np.sin(np.radians(17.58)),
                np.sin(np.radians(88.42))])
WARP = np.array([np.cos(np.radians(75.0)) * np.cos(np.radians(30.0)),
                 np.cos(np.radians(75.0)) * np.sin(np.radians(30.0)),
                 np.sin(np.radians(75.0))])
THETA = float(np.degrees(np.arccos(np.clip(INV @ WARP, -1, 1))))
print(f"warp tilt theta = {THETA:.2f} deg; 1/4 rad = {np.degrees(0.25):.2f} deg (flagged coincidence, no mechanism)")


def tilt_dir():
    """Unit tangent at INV pointing toward WARP, in ecliptic components."""
    t = WARP - (WARP @ INV) * INV
    return t / np.linalg.norm(t)


t = tilt_dir()
print(f"tilt direction (tangent at invariable pole, ecliptic xyz): {np.round(t, 3)}")


# --- Vector model: F = -A cos^2 th - h_z cos th - h_perp sin th, h_z = 0 ---
# Equilibrium: 2A cos th sin th - h_perp cos th = 0 -> h_perp = 2A tan(th) cos(th)
def h_perp_over_A(th_deg):
    th = np.radians(th_deg)
    return 2 * np.tan(th) * np.cos(th)


print(f"\n[V] vector model (h_z = 0): h_perp / A = {h_perp_over_A(THETA):.3f}")
print(f"    small-angle approx tan(theta) = {np.tan(np.radians(THETA)):.3f} (ratio h_perp/(2A+h_z))")
# Numerical check: grid-minimise F with A = 1, h_perp fitted
A = 1.0
h = h_perp_over_A(THETA)


def F_vec(th_deg):
    th = np.radians(th_deg)
    return -A * np.cos(th) ** 2 - h * np.sin(th)


res = minimize_scalar(F_vec, bounds=(0, 90), method="bounded")
print(f"    numerical minimum at theta = {res.x:.2f} deg (target {THETA:.2f}) -- model reproduces the warp")

# --- Axial model: F = -A cos^2 th - B cos^2(alpha - th); A sin 2th = B sin 2(alpha - th) ---
print("\n[A] axial model: B/A as a function of the second-axis angle alpha (deg from invariable pole):")
for alpha in [30, 45, 60, 90, 120, 150]:
    a = np.radians(alpha)
    th = np.radians(THETA)
    ratio = np.sin(2 * th) / np.sin(2 * (a - th))
    print(f"    alpha = {alpha:3d} deg -> B/A = {ratio:.3f}" + ("   <-- d perpendicular: B ~ A" if alpha == 90 else ""))
print("    -> one equation, two unknowns (alpha, B/A): the warp alone cannot fix both.")

# --- Concentration constraint: pole Fisher kappa (from earlier naive means) ---
for label, R in [("80-400 AU naive mean pole", 0.915), ("ETNO n=63", 0.943)]:
    kappa = R * (3 - R ** 2) / (1 - R ** 2)
    print(f"\n{label}: R-bar = {R:.3f} -> Fisher kappa ~ {kappa:.1f} (stiffness/noise ratio; second ratio)")

print("\nPer-bin anisotropy (vector model, h_z = 0):")
print("  50-80 AU : warp ~0 (1 sigma)      -> h_perp/A ~ 0")
print(f"  80-200 AU: warp theta ~13.4 deg   -> h_perp/A ~ {h_perp_over_A(THETA):.2f}")
print("  200-400 AU: warp ~0 (1 sigma)     -> h_perp/A ~ 0")
print("  -> h(a) is a fitted per-bin profile: descriptive, not explanatory.")
print("DONE")
