#!/usr/bin/env python3
"""Neptune-crosser secular simulation with the fitted carrier kernel K.

Method (fixed by task): secular (orbit-averaged) numerical integration.
The carrier is a delocalised field: no close-encounter scattering, so its
only dynamical effect is secular. N-body would be both too expensive and
wrong (scattering is not a carrier effect).

Hamiltonian (quadrupole, relative to carrier axis d):
  H_sec = A(a) * [(2+3e^2)(3cos^2 I - 1) + 15 e^2 sin^2 I cos(2 w')]
Delaunay (G, g=w', H, h=Om'): dG/dt=-dH/dg, dg/dt=dH/dG, dH/dt=0
(Kozai constant H conserved), dh/dt=dH/dH + s(a).  a conserved.

Calibration of A(a) (key, not guessed from normalisation):
  A(a) = A300 * (rho(a)/0.45) * (a/300 AU)^(3/2),
  S(a) = S300 * (rho(a)/0.45), rho = 0.45 * 0.5 * (1+tanh((a-80)/10)).
  A300 is fixed empirically: run a reference particle at 300 AU, measure
  its secular (e-oscillation) period, scale A300 so the period equals the
  D-P3 value (4.5 Gyr for S_min=3.41 m^2/s^2; 1.21 Gyr for S_9=12.66).
  S_9's A300 is then A300(S_min) * (12.66/3.41) (period ~ 1/S), verified
  by re-measuring its period.

Units: AU, yr, solar mass. mu = 4 pi^2 AU^3/yr^2.
1 m^2/s^2 = 4.44996e-8 AU^2/yr^2.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import json, sys

MU = 4.0 * np.pi ** 2
M2S2 = 1.0 / 2.24720e7
S_MIN = 3.41 * M2S2
S_9 = 12.66 * M2S2
T_MIN = 4.5e9
T_9 = 1.21e9
ARCSEC_TO_RAD = 4.84813681109536e-6

# carrier axis d: ecliptic lon 30 deg, lat 75 deg
LON = np.radians(30.0)
LAT = np.radians(75.0)
D = np.array([np.cos(LAT) * np.cos(LON), np.cos(LAT) * np.sin(LON), np.sin(LAT)])
ZEC = np.array([0.0, 0.0, 1.0])
XPR = np.cross(D, ZEC)
XPR = XPR / np.linalg.norm(XPR)
YPR = np.cross(D, XPR)


def rho(a):
    return 0.45 * 0.5 * (1.0 + np.tanh((a - 80.0) / 10.0))


_A_ANCH = np.array([50.0, 80.0, 200.0, 400.0])
_S_ANCH = np.array([-0.203, -0.028, -0.001, -0.00009])  # arcsec/yr


def s_of_a(a):
    """Giants' nodal rate s(a): log|s| vs log a linear interp through bin edges."""
    la = np.log(a)
    lan = np.log(_A_ANCH)
    ls = np.log(-_S_ANCH)
    ls_i = np.interp(la, lan, ls)
    return -np.exp(ls_i) * ARCSEC_TO_RAD  # rad/yr


def ecl_to_carrier(inc, Om, om):
    """Ecliptic (i, Om, om) -> carrier-frame (I, w', Om'). Angles in radians."""
    h = np.array([np.sin(inc) * np.sin(Om), -np.sin(inc) * np.cos(Om), np.cos(inc)])
    n = np.array([np.cos(Om), np.sin(Om), 0.0])
    hxn = np.cross(h, n)
    p = np.cos(om) * n + np.sin(om) * hxn
    cosI = float(np.dot(h, D))
    I = np.arccos(np.clip(cosI, -1, 1))
    dxh = np.cross(D, h)
    nn = np.linalg.norm(dxh)
    nprime = XPR.copy() if nn < 1e-12 else dxh / nn
    e2 = np.cross(h, nprime)
    wp = np.arctan2(float(np.dot(e2, p)), float(np.dot(nprime, p)))
    Op = np.arctan2(float(np.dot(nprime, YPR)), float(np.dot(nprime, XPR)))
    return I, wp, Op


def make_rhs(L, Acoef, sa):
    def rhs(t, y):
        G, g, H, h = y
        Gc = max(G, 1e-9)
        u = 1.0 - (Gc / L) ** 2
        u = min(max(u, 0.0), 0.999999)
        c = H / Gc
        c2 = min(c * c, 1.0)
        s2 = 1.0 - c2
        cg = np.cos(2 * g)
        sg = np.sin(2 * g)
        dG = 30.0 * Acoef * u * s2 * sg
        dg = Acoef * (-(6.0 * Gc / L ** 2) * (3 * c2 - 1.0)
                      - (6.0 * c2 / Gc) * (2 + 3 * u)
                      - (30.0 * Gc / L ** 2) * s2 * cg
                      + (30.0 * u * c2 / Gc) * cg)
        dh = Acoef * ((6.0 * c / Gc) * (2 + 3 * u) - (30.0 * u * c / Gc) * cg) + sa
        return [dG, dg, 0.0, dh]
    return rhs


def hsec_val(L, Acoef, G, g, H):
    Gc = max(G, 1e-9)
    u = 1.0 - (Gc / L) ** 2
    u = min(max(u, 0.0), 0.999999)
    c2 = min((H / Gc) ** 2, 1.0)
    s2 = 1.0 - c2
    return Acoef * ((2 + 3 * u) * (3 * c2 - 1.0) + 15 * u * s2 * np.cos(2 * g))


def run_particle(a, e0, I0, wp0, Acoef, sa, tmax, npts):
    L = np.sqrt(MU * a)
    G0 = L * np.sqrt(1 - e0 ** 2)
    H0 = G0 * np.cos(I0)
    y0 = [G0, wp0, H0, 0.0]
    rhs = make_rhs(L, Acoef, sa)
    tev = np.linspace(0, tmax, npts)
    sol = solve_ivp(rhs, (0, tmax), y0, method="DOP853", t_eval=tev,
                    rtol=1e-9, atol=1e-11)
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    G = sol.y[0]
    e = np.sqrt(np.clip(1 - (G / L) ** 2, 0, 1))
    g = sol.y[1]
    Hh = sol.y[2]
    return tev, e, g, G, Hh, L, sol


def measure_e_period(a, e0, I0, wp0, Acoef, sa, tmax=40e9, npts=1601):
    tev, e, g, G, Hh, L, _ = run_particle(a, e0, I0, wp0, Acoef, sa, tmax, npts)
    pk, _ = find_peaks(e, prominence=0.02)
    if len(pk) < 3:
        return None, tev, e, g
    periods = np.diff(tev[pk])
    return float(np.mean(periods)), tev, e, g


# ---------------- calibration ----------------
print("=== calibration: reference particle a=300 AU, e0=0.3, I0=60deg, w'=90deg ===", flush=True)
REF = dict(a=300.0, e0=0.3, I0=np.radians(60.0), wp0=np.radians(90.0))
sa300 = s_of_a(300.0)
A300_guess = S_MIN  # first guess: A = S
T_meas, _, _, _ = measure_e_period(REF["a"], REF["e0"], REF["I0"], REF["wp0"],
                                   A300_guess, sa300)
print(f"  with A300 = S_min: measured e-period = {T_meas/1e9:.3f} Gyr "
      f"(target {T_MIN/1e9} Gyr)", flush=True)
if T_meas is None:
    raise RuntimeError("calibration: could not measure a period")
A300_MIN = A300_guess * (T_meas / T_MIN)
print(f"  A300(S_min) = {A300_MIN:.4e} AU^2/yr^2  (= {A300_MIN/M2S2:.3f} m^2/s^2; "
      f"A/S = {A300_MIN/S_MIN:.4f})", flush=True)
# verify S_min period
T_ver, _, _, _ = measure_e_period(REF["a"], REF["e0"], REF["I0"], REF["wp0"],
                                  A300_MIN, sa300)
print(f"  verify S_min: e-period = {T_ver/1e9:.3f} Gyr (target 4.5)", flush=True)
# S_9 by scaling (period ~ 1/S)
A300_9 = A300_MIN * (12.66 / 3.41)
T_ver9, _, _, _ = measure_e_period(REF["a"], REF["e0"], REF["I0"], REF["wp0"],
                                   A300_9, sa300)
print(f"  A300(S_9) = {A300_9:.4e} AU^2/yr^2 (= {A300_9/M2S2:.3f} m^2/s^2); "
      f"verify S_9: e-period = {T_ver9/1e9:.3f} Gyr (target 1.21)", flush=True)
if abs(T_ver / T_MIN - 1) > 0.05 or abs(T_ver9 / T_9 - 1) > 0.05:
    raise RuntimeError("calibration verification FAILED (>5% off)")

# ---------------- verification tests ----------------
print("\n=== verification ===", flush=True)
ok = True
# 1. high-I Kozai: large e oscillation + w' libration
tev, e, g, G, Hh, L, _ = run_particle(300.0, 0.3, np.radians(60.0), np.radians(90.0),
                                     A300_MIN, sa300, 4.5e9, 901)
e_amp = e.max() - e.min()
gw = np.degrees(g) % 360.0
# libration: g stays away from 0/180 (does not circulate); check range of unwrapped g
g_unwrap = np.unwrap(g)
g_range = g_unwrap.max() - g_unwrap.min()
print(f"  [V1] high-I (I0=60deg): e range {e.min():.3f}-{e.max():.3f} (amp {e_amp:.3f}); "
      f"g unwrapped range {np.degrees(g_range):.0f} deg", flush=True)
v1 = (e_amp > 0.3) and (g_range < np.radians(200))
print(f"       -> {'PASS' if v1 else 'FAIL'} (expect amp>0.3, libration range<200deg)", flush=True)
ok &= v1
# 2. low-I: small e amplitude
tev, e, g, G, Hh, L, _ = run_particle(300.0, 0.3, np.radians(15.0), np.radians(0.0),
                                     A300_MIN, sa300, 4.5e9, 901)
e_amp2 = e.max() - e.min()
print(f"  [V2] low-I (I0=15deg): e range {e.min():.3f}-{e.max():.3f} (amp {e_amp2:.3f})", flush=True)
v2 = e_amp2 < 0.10
print(f"       -> {'PASS' if v2 else 'FAIL'} (expect amp<0.10)", flush=True)
ok &= v2
# 3. conservation of H_sec and Kozai H
tev, e, g, G, Hh, L, _ = run_particle(300.0, 0.5, np.radians(45.0), np.radians(45.0),
                                     A300_MIN, sa300, 4.5e9, 901)
hs = np.array([hsec_val(L, A300_MIN, G[k], g[k], Hh[k]) for k in range(len(tev))])
rel_hsec = (hs.max() - hs.min()) / max(abs(hs).max(), 1e-30)
rel_H = (Hh.max() - Hh.min()) / max(abs(Hh).max(), 1e-30)
print(f"  [V3] conservation: H_sec rel variation {rel_hsec:.2e}; H rel variation {rel_H:.2e}", flush=True)
v3 = (rel_hsec < 1e-6) and (rel_H < 1e-6)
print(f"       -> {'PASS' if v3 else 'FAIL'} (expect <1e-6)", flush=True)
ok &= v3
# 4. period already verified above (V4)
v4 = (abs(T_ver / T_MIN - 1) <= 0.05) and (abs(T_ver9 / T_9 - 1) <= 0.05)
print(f"  [V4] period match: S_min {T_ver/1e9:.3f}/4.5 Gyr, S_9 {T_ver9/1e9:.3f}/1.21 Gyr "
      f"-> {'PASS' if v4 else 'FAIL'}", flush=True)
ok &= v4
if not ok:
    raise RuntimeError("VERIFICATION FAILED -- stop, do not proceed to main runs")
print("ALL VERIFICATIONS PASSED", flush=True)

# ---------------- main runs ----------------
print("\n=== main runs: N=200, a~U[100,400], q~U[30,80], i~U[0,40deg] ===", flush=True)
rng = np.random.default_rng(42)
N = 200
a_s = rng.uniform(100.0, 400.0, N)
q_s = rng.uniform(30.0, 80.0, N)
i_s = np.radians(rng.uniform(0.0, 40.0, N))
Om_s = rng.uniform(0, 2 * np.pi, N)
om_s = rng.uniform(0, 2 * np.pi, N)
e0_s = 1.0 - q_s / a_s
I0_s = np.empty(N)
wp0_s = np.empty(N)
for k in range(N):
    I0_s[k], wp0_s[k], _ = ecl_to_carrier(i_s[k], Om_s[k], om_s[k])
print(f"  population: e0 range {e0_s.min():.3f}-{e0_s.max():.3f}; "
      f"I0 (rel carrier) range {np.degrees(I0_s).min():.1f}-{np.degrees(I0_s).max():.1f} deg; "
      f"I0>39.23deg: {(I0_s > np.radians(39.23)).sum()} particles", flush=True)

TMAX = 4.5e9
NPTS = 901
results = {}
for tag, A300 in [("S=0 control", 0.0), ("S_min", A300_MIN), ("S_9", A300_9)]:
    qmins = np.empty(N)
    emaxs = np.empty(N)
    qfinal = np.empty(N)
    frac_time_cross = np.empty(N)
    for k in range(N):
        a = a_s[k]
        Acoef = A300 * (rho(a) / 0.45) * (a / 300.0) ** 1.5 if A300 > 0 else 0.0
        sa = s_of_a(a)
        tev, e, g, G, Hh, L, _ = run_particle(a, e0_s[k], I0_s[k], wp0_s[k],
                                             Acoef, sa, TMAX, NPTS)
        q = a * (1.0 - e)
        qmins[k] = q.min()
        emaxs[k] = e.max()
        qfinal[k] = q[-1]
        frac_time_cross[k] = (q < 30.0).mean()
        if (k + 1) % 50 == 0:
            print(f"  {tag}: {k+1}/{N} particles", flush=True)
    cross = qmins < 30.0
    frac = cross.mean()
    sig = np.sqrt(frac * (1 - frac) / N) if 0 < frac < 1 else 1.0 / N
    inst_final = (qfinal < 30.0).mean()
    inst_timeavg = frac_time_cross.mean()
    results[tag] = dict(frac=float(frac), sigma=float(sig), n_cross=int(cross.sum()),
                        qmins=qmins.tolist(), emaxs=emaxs.tolist(),
                        qfinal=qfinal.tolist(), frac_time_cross=frac_time_cross.tolist(),
                        inst_final=float(inst_final), inst_timeavg=float(inst_timeavg))
    print(f"  {tag}: EVER-crossed {cross.sum()}/{N} = {frac*100:.1f}% +- {sig*100:.1f}% "
          f"(1 sigma); instantaneous at t=4.5Gyr {(inst_final)*100:.1f}%; "
          f"time-averaged instantaneous {inst_timeavg*100:.1f}%; "
          f"median q_min {np.median(qmins):.1f} AU; median e_max {np.median(emaxs):.3f}", flush=True)
    if cross.any():
        idx = np.where(cross)[0]
        print(f"       crossers: a={np.round(a_s[idx],0)}, e0={np.round(e0_s[idx],2)}, "
              f"I0={np.round(np.degrees(I0_s[idx]),0)} deg, q_min={np.round(qmins[idx],1)} AU",
              flush=True)

with open("crosser_simulation_results.json", "w") as f:
    json.dump({"A300_MIN": A300_MIN, "A300_9": A300_9, "T_ver": T_ver, "T_ver9": T_ver9,
               "a": a_s.tolist(), "q0": q_s.tolist(), "e0": e0_s.tolist(),
               "I0_deg": np.degrees(I0_s).tolist(),
               "results": {k: {kk: vv for kk, vv in v.items()} for k, v in results.items()}},
              f)
print("saved crosser_simulation_results.json")
print("DONE")
