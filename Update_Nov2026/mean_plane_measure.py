#!/usr/bin/env python3
"""Measure the 200-400 AU mean plane from current MPCORB with the Siraj et al.
2025b bias-free estimator (conditional likelihood, von Mises planes), plus
Brown & Pan (2004) and naive orbit-normal averaging for comparison.
80-200 AU (known warp) and 50-80 AU (known no warp) as controls.

Per object: r_hat (heliocentric position at its own epoch) and J_hat (orbit pole)
from (a, e, i, Om, om, M) via Kepler. Siraj logL(m, gamma) =
  gamma * m . sum(J_hat) - sum log I0(gamma * sin(theta_i)),
  cos(theta_i) = m . r_hat_i.  Significance vs invariable plane: sqrt(2 dlogL).
Resonance cut (crude proxy for Siraj's SBDynT classification): exclude a within
0.5 AU of any n:1 (n=2..45) or n:2 (odd n=5..45) Neptune resonance. Reported
with and without the cut. Uncertainty parameter U<=6 as in Siraj et al.
File epoch ~2026-06-09 (packed K2669); positions at each object's own epoch.
"""
import numpy as np
import pandas as pd
from scipy.special import i0
from scipy.optimize import minimize, minimize_scalar

K_GAUSS = 0.01720209895  # AU^3/2 / day
INV = np.array([np.cos(np.radians(88.42)) * np.cos(np.radians(17.58)),
                np.cos(np.radians(88.42)) * np.sin(np.radians(17.58)),
                np.sin(np.radians(88.42))])
INV_PQ = (np.sin(np.radians(0.8)) * np.sin(np.radians(107.58)),
          np.sin(np.radians(0.8)) * np.cos(np.radians(107.58)))


def kepler_E(M, e, tol=1e-12):
    E = M if e < 0.8 else np.pi
    for _ in range(50):
        f = E - e * np.sin(E) - M
        fp = 1 - e * np.cos(E)
        d = f / fp
        E -= d
        if abs(d) < tol:
            break
    return E


def state_vectors(a, e, inc, Om, om, Mdeg):
    """Return r_hat, v_t_hat (transverse velocity unit), J_hat (ecliptic J2000)."""
    M = np.radians(Mdeg)
    E = kepler_E(M, e)
    r = a * (1 - e * np.cos(E))
    nu = np.arctan2(np.sqrt(1 - e ** 2) * np.sin(E), np.cos(E) - e)
    # perifocal position / velocity (units: AU, AU/day; only directions used)
    x_p = r * np.cos(nu)
    y_p = r * np.sin(nu)
    vx_p = -np.sqrt(K_GAUSS ** 2 * a) / r * np.sin(E) * a  # -sqrt(mu a)/r * sinE
    vy_p = np.sqrt(K_GAUSS ** 2 * a) / r * np.sqrt(1 - e ** 2) * np.cos(E) * a
    ci, si = np.cos(np.radians(inc)), np.sin(np.radians(inc))
    cO, sO = np.cos(np.radians(Om)), np.sin(np.radians(Om))
    co, so = np.cos(np.radians(om)), np.sin(np.radians(om))
    R = np.array([[cO * co - sO * so * ci, -cO * so - sO * co * ci, sO * si],
                  [sO * co + cO * so * ci, -sO * so + cO * co * ci, -cO * si],
                  [so * si, co * si, ci]])
    rvec = R @ np.array([x_p, y_p, 0.0])
    vvec = R @ np.array([vx_p, vy_p, 0.0])
    r_hat = rvec / np.linalg.norm(rvec)
    vt = vvec - (vvec @ r_hat) * r_hat
    vt_hat = vt / np.linalg.norm(vt)
    J_hat = np.array([si * sO, -si * cO, ci])
    return r_hat, vt_hat, J_hat


colspecs = [(0, 7), (8, 13), (14, 19), (20, 25), (26, 35), (37, 46),
            (48, 57), (59, 68), (70, 79), (80, 91), (92, 103), (104, 106)]
names = ["desig", "H", "G", "epoch", "M", "omega", "Omega", "i", "e", "n", "a", "U"]
df = pd.read_fwf("MPCORB-202611.DAT", colspecs=colspecs, names=names,
                 header=None, skiprows=43)
for c in ["a", "e", "omega", "Omega", "i", "M", "U"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df.dropna(subset=["a", "e", "omega", "Omega", "i", "M"], inplace=True)
df["q"] = df["a"] * (1 - df["e"])
df["U"] = df["U"].fillna(9)
print(f"parsed {len(df)} objects; U<=6: {(df['U'] <= 6).sum()}")

# resonance centres (Neptune a_N = 30.07)
res = set()
for nn in range(2, 46):
    res.add(30.07 * nn ** (2 / 3))
for nn in range(5, 46, 2):
    res.add(30.07 * (nn / 2) ** (2 / 3))
res = np.array(sorted(res))


def near_res(a):
    return np.min(np.abs(res - a)) < 0.5


def pq_to_m(p, q):
    s = np.hypot(p, q)
    if s >= 0.5:
        return None
    i0 = 2 * np.arcsin(np.clip(s, 0, 1))
    Om0 = np.arctan2(p, q)
    return np.array([np.sin(i0) * np.sin(Om0), -np.sin(i0) * np.cos(Om0), np.cos(i0)])


def m_to_pq(m):
    lon = np.degrees(np.arctan2(m[1], m[0])) % 360
    lat = np.degrees(np.arcsin(np.clip(m[2], -1, 1)))
    i0 = 90 - lat
    Om0 = (lon + 90) % 360
    return (np.sin(np.radians(i0 / 2)) * np.sin(np.radians(Om0)),
            np.sin(np.radians(i0 / 2)) * np.cos(np.radians(Om0)), i0, Om0)


def siraj_fit(R, J):
    """R: (n,3) r_hat; J: (n,3) J_hat. Returns dict with best (p,q,gamma), logLmax,
    significance vs invariable plane, and 1-sigma contour max tilt."""
    Jsum = J.sum(axis=0)

    def neglogL(x):
        p, q, lg = x
        m = pq_to_m(p, q)
        if m is None or lg > 8:
            return 1e9
        gamma = np.exp(lg)
        sinth = np.sqrt(np.clip(1 - (R @ m) ** 2, 0, 1))
        return -(gamma * (m @ Jsum) - np.sum(np.log(i0(gamma * sinth))))

    x0 = np.array([INV_PQ[0], INV_PQ[1], np.log(20.0)])
    best = minimize(neglogL, x0, method="Nelder-Mead",
                    options={"maxiter": 4000, "xatol": 1e-5, "fatol": 1e-7})
    p, q, lg = best.x
    gamma = np.exp(lg)
    logLmax = -best.fun
    m = pq_to_m(p, q)
    _, _, i0_deg, Om0 = m_to_pq(m)
    tilt = float(np.degrees(np.arccos(np.clip(m @ INV, -1, 1))))

    def logL_at_pq(pp, qq):
        mm = pq_to_m(pp, qq)
        if mm is None:
            return -np.inf
        sinth = np.sqrt(np.clip(1 - (R @ mm) ** 2, 0, 1))

        def f(lg_):
            g = np.exp(lg_)
            return -(g * (mm @ Jsum) - np.sum(np.log(i0(g * sinth))))
        r = minimize_scalar(f, bounds=(-4, 8), method="bounded")
        return -r.fun

    logL_inv = logL_at_pq(*INV_PQ)
    sig = float(np.sqrt(max(0.0, 2 * (logLmax - logL_inv))))
    # 1-sigma contour: 2 dlogL = 2.30; max tilt inside it (coarse grid)
    grid = np.linspace(-0.30, 0.30, 61)
    max_tilt_1sig = 0.0
    for pp in grid:
        for qq in grid:
            if np.hypot(pp, qq) >= 0.5:
                continue
            ll = logL_at_pq(pp, qq)
            if 2 * (logLmax - ll) <= 2.30:
                mm = pq_to_m(pp, qq)
                max_tilt_1sig = max(max_tilt_1sig,
                                    float(np.degrees(np.arccos(np.clip(mm @ INV, -1, 1)))))
    return dict(n=len(R), i0=i0_deg, Om0=Om0, gamma=gamma, logLmax=logLmax,
                sig=sig, max_tilt_1sig=max_tilt_1sig, m=m)


def brown_pan(Vt):
    Mmat = Vt.T @ Vt
    vals, vecs = np.linalg.eigh(Mmat)
    n = vecs[:, 0]
    if n[2] < 0:
        n = -n
    _, _, i0, Om0 = m_to_pq(n)
    return i0, Om0


def naive(J):
    m = J.mean(axis=0)
    m = m / np.linalg.norm(m)
    _, _, i0, Om0 = m_to_pq(m)
    return i0, Om0


def analyze(label, sub):
    R, Vt, J = [], [], []
    for _, r in sub.iterrows():
        rh, vth, Jh = state_vectors(r["a"], r["e"], r["i"], r["Omega"], r["omega"], r["M"])
        R.append(rh)
        Vt.append(vth)
        J.append(Jh)
    R, Vt, J = map(np.array, (R, Vt, J))
    s = siraj_fit(R, J)
    bp_i, bp_O = brown_pan(Vt)
    nv_i, nv_O = naive(J)
    print(f"\n--- {label} (n={s['n']}) ---")
    print(f"Siraj bias-free: i0={s['i0']:.1f} deg, Om0={s['Om0']:.0f} deg, gamma={s['gamma']:.1f}, "
          f"significance vs invariable plane = {s['sig']:.2f} sigma; "
          f"max tilt inside 1-sigma contour = {s['max_tilt_1sig']:.1f} deg")
    print(f"Brown & Pan: i0={bp_i:.1f} deg, Om0={bp_O:.0f} deg")
    print(f"Naive mean pole: i0={nv_i:.1f} deg, Om0={nv_O:.0f} deg")


base = df[(df["q"] > 30) & (df["U"] <= 6)].copy()
base["res"] = base["a"].apply(near_res)
for lo, hi, label in [(50, 80, "50-80 AU (control: no warp expected)"),
                      (80, 200, "80-200 AU (control: warp expected)"),
                      (200, 400, "200-400 AU (TARGET)"),
                      (80, 400, "80-400 AU (combined)")]:
    sub = base[(base["a"] >= lo) & (base["a"] < hi if hi < 400 else base["a"] <= 400)]
    analyze(f"{label} -- all", sub)
    analyze(f"{label} -- resonance-cut", sub[~sub["res"]])
print("DONE")
