#!/usr/bin/env python3
"""D-P3b: candidate radial scales for the ~80 AU turn-on of the coupling kernel.

Questions:
  A. Neptune resonance locations near 80 AU (a = a_N * (p/q)^(2/3)).
  B. Solar J2 nodal precession vs giants' secular nodal precession at 50-400 AU
     (is the Laplace-plane / J2 transition anywhere near 80 AU?).
  C. Free nodal rate s(a) vs Neptune's nodal eigenfrequency s8 = -0.69 "/yr
     (secular resonance location from known planets).
  D. Bodies picture: where does the perturber's weight equal the giants' weight?
     (the "turn-on" scale for light Y and heavy Y branches).
  E. Fine tilt curves for the two branches (where does each branch's rise sit?).
  F. Secular resonance with a hypothetical Planet Y (its own nodal rate s_Y,
     and where the test-particle rate equals it).
  G. Radial eigenfunctions on a ball: spherical Bessel j_2(kr) — can a single
     radial mode give a flat profile 0.4-0.5 from 80 to 400 AU?
  H. Rubin sample size needed to settle the 80-140 AU shape (11 deg vs 1 deg
     at 100 AU).

Reuses the validated Laplace coefficient from bodies_coupling_fit.py.
"""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.special import spherical_jn
from scipy.optimize import brentq

# ---------- Laplace coefficient b^(1)_{3/2} ----------
def b132(alpha):
    if alpha <= 0:
        return 0.0
    f = lambda psi: np.cos(psi) / (1 - 2*alpha*np.cos(psi) + alpha**2)**1.5
    val, _ = quad(f, 0, np.pi, limit=200)
    return 2*val/np.pi

agrid = np.concatenate([np.linspace(0.01, 0.95, 95), [0.97, 0.99]])
bgrid = np.array([b132(a) for a in agrid])
binterp = interp1d(agrid, bgrid, kind="cubic", fill_value="extrapolate")
def b(alpha):
    return float(binterp(min(max(alpha, 0.01), 0.99)))

GIANTS = [("Jupiter", 317.8, 5.20, 1.30, 100.46),
          ("Saturn", 95.2, 9.58, 2.49, 113.67),
          ("Uranus", 14.5, 19.22, 0.77, 74.01),
          ("Neptune", 17.1, 30.07, 1.77, 131.78)]
M_SUN_EARTH = 332946.0  # solar mass in Earth masses
ARCSEC_PER_RAD = 206265.0

def weight(m, a_j, a, cap=0.95):
    """Secular weight (Laplace-Lagrange), alpha capped at `cap` like the fit."""
    if a_j < a:
        alpha = min(a_j/a, cap)
        return m*alpha*b(alpha)
    else:
        alpha = min(a/a_j, cap)
        return m*alpha**2*b(alpha)

def sec_rate(a, extra=()):
    """Free nodal precession rate s(a) in arcsec/yr from giants (+ extra bodies).
    s = -(n/4) * sum_j (m_j/M_sun) * (alpha*b or alpha^2*b), n = 2*pi/a^1.5 rad/yr."""
    n = 2*np.pi/a**1.5  # rad/yr
    tot = 0.0
    for _, m, a_j, _, _ in list(GIANTS) + list(extra):
        if a_j < a:
            alpha = min(a_j/a, 0.99)
            tot += (m/M_SUN_EARTH)*alpha*b(alpha)
        else:
            alpha = min(a/a_j, 0.99)
            tot += (m/M_SUN_EARTH)*alpha**2*b(alpha)
    return -(n/4)*tot*ARCSEC_PER_RAD

# ================= A. Neptune resonances =================
print("=== A. Neptune mean-motion resonances (a_N = 30.07 AU) ===")
for p, q in [(2,1),(5,2),(3,1),(7,2),(4,1),(9,2),(5,1),(6,1)]:
    a_res = 30.07*(p/q)**(2/3)
    print(f"  {p}:{q}  a = {a_res:6.2f} AU")

# ================= B. J2 vs secular =================
print("\n=== B. Solar J2 nodal rate vs giants' secular nodal rate ===")
J2 = 2.2e-7
R_SUN_AU = 0.00465047
print("  a(AU)   J2 rate(\"/yr)   secular rate(\"/yr)   ratio J2/secular")
for a in [5, 10, 30, 50, 80, 100, 140, 200, 300, 400]:
    n = 2*np.pi/a**1.5
    j2_rate = -1.5*J2*n*(R_SUN_AU/a)**2*ARCSEC_PER_RAD
    s = sec_rate(a)
    print(f"  {a:4d}   {j2_rate: .3e}        {s: .3e}        {abs(j2_rate/s):.2e}")

# ================= C. secular resonance vs s8 =================
print("\n=== C. Free nodal rate s(a) vs Neptune eigenfrequency s8 = -0.69 \"/yr ===")
for a in [30, 35, 40, 45, 50, 80, 100, 140, 200, 300]:
    print(f"  a={a:4d} AU: s = {sec_rate(a): .4f} \"/yr   |s/s8| = {abs(sec_rate(a)/0.69):.4f}")

# ================= D. weight turn-on =================
print("\n=== D. Where does the perturber weight equal the giants' weight? ===")
def w_giants(a):
    return sum(weight(m, aj, a) for _, m, aj, _, _ in GIANTS)
for label, mY, aY in [("light Y", 0.1, 150.0), ("heavy Y", 4.0, 128.0)]:
    f = lambda a: weight(mY, aY, a) - w_giants(a)
    try:
        a_star = brentq(f, 40, 0.98*aY)
    except ValueError:
        a_star = float("nan")
    print(f"  {label} (m={mY}, a={aY}): wY = wGiants at a* = {a_star:.1f} AU")
    print(f"    a: " + "  ".join(f"{a:>5d}" for a in [50,65,80,100,120,140,180,250,300]))
    print(f"    wY/wG: " + "  ".join(f"{weight(mY,aY,a)/w_giants(a):5.3f}" for a in [50,65,80,100,120,140,180,250,300]))

# ================= E. fine tilt curves =================
print("\n=== E. Forced-plane tilt curves, fine a-grid (linear secular, alpha cap 0.95) ===")
INV_P = np.sin(np.radians(0.8))*np.sin(np.radians(107.58))
INV_Q = np.sin(np.radians(0.8))*np.cos(np.radians(107.58))
def incl_vec(i_deg, Om_deg):
    return np.array([np.sin(np.radians(i_deg/2))*np.sin(np.radians(Om_deg)),
                     np.sin(np.radians(i_deg/2))*np.cos(np.radians(Om_deg))])
def forced_plane(a, extra=()):
    bodies = GIANTS + list(extra)
    num = np.zeros(2); den = 0.0
    for _, m, a_j, i, Om in bodies:
        w = weight(m, a_j, a)
        num += w*incl_vec(i, Om); den += w
    p, q = num/den
    tilt = 2*np.degrees(np.arcsin(np.clip(np.hypot(p-INV_P, q-INV_Q), 0, 1)))
    node = float(np.degrees(np.arctan2(p, q)) % 360)
    return tilt, node
Y_LIGHT = ("Y light", 0.1, 150.0, 25.0, 120.0)
Y_HEAVY = ("Y heavy", 4.0, 128.0, 16.8, 119.0)
P9 = ("P9 BB21", 6.2, 380.0, 16.0, 100.0)
print(f"  {'a':>4} {'heavyY':>7} {'lightY':>7} {'lightY+P9':>10} {'heavyY+P9':>10}  (deg; node in parens)")
for a in [60, 70, 80, 90, 100, 110, 120, 140, 160, 200, 300]:
    th, nh = forced_plane(a, extra=[Y_HEAVY])
    tl, nl = forced_plane(a, extra=[Y_LIGHT])
    tlp, nlp = forced_plane(a, extra=[Y_LIGHT, P9])
    thp, nhp = forced_plane(a, extra=[Y_HEAVY, P9])
    print(f"  {a:4d} {th:6.2f}  {tl:6.2f}  {tlp:6.2f}/{nlp:3.0f}  {thp:6.2f}/{nhp:3.0f}")

# ================= F. secular resonance with Planet Y =================
print("\n=== F. Planet Y's own nodal rate s_Y vs test-particle rate s(a) ===")
# s_Y: Y's nodal precession from the giants (all interior to Y)
for mY, aY in [(0.1, 150.0), (4.0, 128.0)]:
    nY = 2*np.pi/aY**1.5
    tot = sum((m/M_SUN_EARTH)*(aj/aY)*b(aj/aY) for _, m, aj, _, _ in GIANTS)
    sY = -(nY/4)*tot*ARCSEC_PER_RAD
    print(f"  Y (m={mY}, a={aY}): s_Y = {sY:.5f} \"/yr")
    # test-particle rate from giants only, at a < aY: where does it equal s_Y?
    f = lambda a: sec_rate(a) - sY
    try:
        a_res = brentq(f, 30, 0.95*aY)
        print(f"    s_giants(a) = s_Y at a = {a_res:.1f} AU  (secular resonance with Y, giants only)")
    except ValueError:
        print(f"    s_giants(a) never equals s_Y in [30, {0.95*aY:.0f}] AU (no secular resonance there)")

# ================= G. radial Bessel modes on a ball =================
print("\n=== G. Radial eigenfunction j_2(kr) on a ball, Dirichlet BC ===")
xs = np.linspace(0.01, 8, 8000)
j2 = spherical_jn(2, xs)
x_peak = xs[np.argmax(j2)]
x_zero = xs[np.where((xs > 4) & (j2 < 0))[0][0]]  # first zero after the peak
# refine zero
from scipy.optimize import brentq as bq
x_zero = bq(lambda x: spherical_jn(2, x), 5.0, 6.5)
print(f"  j_2 peaks at kr = {x_peak:.3f} (j2 = {spherical_jn(2, x_peak):.4f}); first zero at kr = {x_zero:.4f}")
for R in [269.0, 500.0]:
    # R chosen so the peak sits at 140 AU (R=269) or at 260 AU (R=500)
    k = x_zero / R
    print(f"  Ball radius R = {R:.0f} AU (peak at {x_peak/k:.0f} AU, field = 0 at r >= {R:.0f} AU):")
    print(f"    r:      " + "  ".join(f"{r:>6d}" for r in [65, 80, 100, 140, 200, 300, 400]))
    prof = [spherical_jn(2, k*r) for r in [65, 80, 100, 140, 200, 300, 400]]
    ref = spherical_jn(2, k*140)
    print(f"    j2:     " + "  ".join(f"{v:6.3f}" for v in prof))
    print(f"    j2^2/140: " + "  ".join(f"{(v/ref)**2:6.3f}" for v in prof))
    print(f"    (observed rho(a) ~ 0 below 80 AU, ~0.4-0.5 flat from 80 to 400 AU)")

# ================= H. Rubin sample size =================
print("\n=== H. Rubin sample size to settle the 80-140 AU shape ===")
# current: 1-sigma ~5.6 deg at n=32 (80-200, Table 1), ~5.9 deg at n=29 (200-400)
sig0, n0 = 5.75, 30.0
for sig_target in [3.3, 2.0, 1.5]:
    n_need = n0*(sig0/sig_target)**2
    print(f"  target 1-sigma = {sig_target:.1f} deg -> n ~ {n_need:.0f} nonresonant objects per bin")
print("  (3.3 deg separates 11 deg vs 1 deg at 100 AU at 3-sigma; 2 deg maps the rise shape)")

# ================= I. washout times (free pattern, no maintainer) =================
print("\n=== I. Washout time of a free tilted pattern (differential precession) ===")
print("  spread of s(a) across each bin -> time for poles to phase-mix by 180 deg")
for lo, hi in [(50, 80), (80, 200), (200, 400)]:
    ds = abs(sec_rate(lo) - sec_rate(hi))  # arcsec/yr
    deg_per_myr = ds*1e6/3600.0
    t_mix_myr = 180.0/deg_per_myr
    print(f"  {lo:3d}-{hi:3d} AU: s({lo})={sec_rate(lo): .4f}, s({hi})={sec_rate(hi): .5f} \"/yr,"
          f" spread {ds:.4f} \"/yr = {deg_per_myr:.2f} deg/Myr -> 180-deg mix in {t_mix_myr:.1f} Myr")
print("  (solar-system age 4500 Myr; a free pattern is gone long before today)")

# ================= J. tidal-law Mercury suppression =================
print("\n=== J. Tidal (r^2) law: suppression at Mercury vs a direction-only law ===")
a_merc = 0.387
for a_ref in [140.0]:
    print(f"  (r_Mercury / {a_ref:.0f} AU)^2 = {(a_merc/a_ref)**2:.2e}")
print("  D-P3 numbers: direction-only S_min -> 0.51 \"/cy at Mercury (bound 0.01);")
print("  with a tidal law the Mercury rate is multiplied by ~7.6e-6 -> ~4e-6 \"/cy (negligible).")
print("DONE")
