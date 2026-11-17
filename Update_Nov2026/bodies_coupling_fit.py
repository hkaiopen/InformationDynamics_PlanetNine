#!/usr/bin/env python3
"""Derive the mode axis from K(x) sourced by celestial bodies (visible + invisible).

Hypothesis (user's, tested not assumed): K(x) = sum over bodies of their coupling,
each body contributing its own orbital plane, weighted by secular coupling strength.
Mode axis = forced plane = weighted mean of the bodies' inclination vectors.

Secular weights (Laplace-Lagrange, test particle at semimajor axis a):
  perturber interior (a_j < a): w_j ~ m_j * alpha * b^(1)_{3/2}(alpha), alpha = a_j/a
  perturber exterior (a_j > a): w_j ~ m_j * alpha^2 * b^(1)_{3/2}(alpha), alpha = a/a_j
  (checked against quadrupole limits: interior -> 3 m_j (a_j/a)^2; exterior -> 3 m_j (a/a_j)^3)
Inclination vector: p = sin(i/2) sin(Om), q = sin(i/2) cos(Om) (Siraj et al. 2025b form).
Forced plane I(a) = sum w_j I_j / sum w_j; tilt = angle of I(a) from invariable plane.

Bodies: Jupiter, Saturn, Uranus, Neptune (J2000 mean elements);
Planet Nine BB21 nominal (m=6.2, a=380, i=16, Om=100);
Planet Y scanned (Siraj et al. 2025b N-body ranges: m 0.06-1, a 100-200, i >~10).
Observed warp to match: tilt ~13 deg (8-20 allowed), node ~120 deg, at 80-200 AU;
tilt < ~5 deg at 50-80 AU; ~0 at 200-400 AU (paper: ~1 sigma).
"""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

# --- Laplace coefficient b^(1)_{3/2}(alpha) ---
def b132(alpha):
    if alpha <= 0:
        return 0.0
    f = lambda psi: np.cos(psi) / (1 - 2 * alpha * np.cos(psi) + alpha ** 2) ** 1.5
    val, _ = quad(f, 0, np.pi, limit=200)
    return 2 * val / np.pi

# validate against series b = 3a + 45/8 a^3 + 525/64 a^5 at small alpha
for a in [0.05, 0.1, 0.2]:
    series = 3 * a + 45 / 8 * a ** 3 + 525 / 64 * a ** 5
    print(f"b(1,3/2)({a}): numeric {b132(a):.4f}, series {series:.4f}")

agrid = np.concatenate([np.linspace(0.01, 0.95, 95), [0.97, 0.99]])
bgrid = np.array([b132(a) for a in agrid])
binterp = interp1d(agrid, bgrid, kind="cubic", fill_value="extrapolate")
print(f"b(0.6) = {float(binterp(0.6)):.2f}, b(0.667) = {float(binterp(0.667)):.2f}")

INV_P = np.sin(np.radians(0.8)) * np.sin(np.radians(107.58))
INV_Q = np.sin(np.radians(0.8)) * np.cos(np.radians(107.58))


def incl_vec(i_deg, Om_deg):
    return np.array([np.sin(np.radians(i_deg / 2)) * np.sin(np.radians(Om_deg)),
                     np.sin(np.radians(i_deg / 2)) * np.cos(np.radians(Om_deg))])


def weight(m, a_j, a):
    if a_j < a:
        alpha = a_j / a
        return m * alpha * float(binterp(alpha))
    else:
        alpha = a / a_j
        return m * alpha ** 2 * float(binterp(alpha))


# (name, mass M_E, a AU, i deg, Om deg); masses in Earth masses
GIANTS = [("Jupiter", 317.8, 5.20, 1.30, 100.46),
          ("Saturn", 95.2, 9.58, 2.49, 113.67),
          ("Uranus", 14.5, 19.22, 0.77, 74.01),
          ("Neptune", 17.1, 30.07, 1.77, 131.78)]
P9 = ("Planet Nine (BB21)", 6.2, 380.0, 16.0, 100.0)


def forced_plane(a, extra=()):
    bodies = GIANTS + list(extra)
    num = np.zeros(2)
    den = 0.0
    for _, m, a_j, i, Om in bodies:
        w = weight(m, a_j, a)
        num += w * incl_vec(i, Om)
        den += w
    p, q = num / den
    tilt = 2 * np.degrees(np.arcsin(np.clip(np.hypot(p - INV_P, q - INV_Q), 0, 1)))
    node = float(np.degrees(np.arctan2(p, q)) % 360)
    return tilt, node


print("\n--- giants only (check: forced plane ~ invariable plane at all a) ---")
for a in [65, 140, 300]:
    t, n = forced_plane(a)
    print(f"a = {a} AU: tilt {t:.2f} deg, node {n:.0f} deg")

print("\n--- giants + Planet Nine BB21 (m=6.2, a=380, i=16, Om=100) ---")
for a in [65, 140, 300]:
    t, n = forced_plane(a, extra=[P9])
    print(f"a = {a} AU: tilt {t:.2f} deg, node {n:.0f} deg")

print("\n--- Planet Y scan (node fixed at 120 deg): tilt at a=65 / 140 / 300 AU ---")
print("target: tilt(65) < 5, tilt(140) in 8-20 with node ~120, tilt(300) small (paper: no warp)")
good = []
for mY in [0.06, 0.1, 0.3, 1.0]:
    for aY in [100, 150, 200]:
        for iY in [10, 20, 30, 40, 60]:
            Y = ("Planet Y", mY, float(aY), float(iY), 120.0)
            t65, _ = forced_plane(65, extra=[P9, Y])
            t140, n140 = forced_plane(140, extra=[P9, Y])
            t300, _ = forced_plane(300, extra=[P9, Y])
            flag = ""
            if t65 < 5 and 8 <= t140 <= 20 and 90 <= n140 <= 150:
                flag = "  <-- MATCHES warp (node+amplitude) at 80-200, quiet at 50-80"
                good.append((mY, aY, iY, t65, t140, n140, t300))
            print(f"m={mY:4.2f} a={aY:3d} i={iY:2d}: tilt65={t65:5.2f} tilt140={t140:5.2f} "
                  f"node140={n140:5.0f} tilt300={t300:5.2f}{flag}")
print(f"\nMatching combinations (with P9 included): {len(good)}")
for g in good:
    print(f"  m={g[0]} a={g[1]} i={g[2]}: tilt65={g[3]:.2f} tilt140={g[4]:.2f} node140={g[5]:.0f} tilt300={g[6]:.2f}")
print("DONE")
