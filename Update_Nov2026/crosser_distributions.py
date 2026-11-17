#!/usr/bin/env python3
"""Crosser distribution details: carrier (fitted K) vs P9-secular control.

Extends crosser_simulation.py (do NOT modify it): same secular integrator,
same calibration (A300 read from crosser_simulation_results.json), same
population (seed 42). Records full trajectories and computes:
  * a distribution of crossers (ever-crossed + snapshot-weighted)
  * dwell-time distribution (duty cycle + single-episode durations)
  * omega / omega' distributions at crossing (Rayleigh + axial tests)
  * bonus: q and i distributions of crossing snapshots (directly
    comparable to Batygin et al. 2024 Figs. 3 & 5)
Plus a P9-secular control: same quadrupole form, axis = P9 orbit pole
(lon 10 deg, lat 74 deg; BB21 i=16, Om=100), A300 = A300_9.
Caveat (stated in the note): P9's real crossers are also made by
scattering, which a secular model cannot capture; the control tests
only whether the AXIS direction matters, not P9's full prediction.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy import stats
import json, sys

MU = 4.0 * np.pi ** 2
M2S2 = 1.0 / 2.24720e7
ARCSEC_TO_RAD = 4.84813681109536e-6
ZEC = np.array([0.0, 0.0, 1.0])
TMAX = 4.5e9
NPTS = 1801
DT = TMAX / (NPTS - 1)  # 2.5 Myr


def axis_from_lonlat(lon_deg, lat_deg):
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon),
                     np.sin(lat)])


def frame_basis(D):
    XPR = np.cross(D, ZEC)
    XPR = XPR / np.linalg.norm(XPR)
    YPR = np.cross(D, XPR)
    E = np.column_stack([XPR, YPR, D])
    return XPR, YPR, E


def rho(a):
    return 0.45 * 0.5 * (1.0 + np.tanh((a - 80.0) / 10.0))


_A_ANCH = np.array([50.0, 80.0, 200.0, 400.0])
_S_ANCH = np.array([-0.203, -0.028, -0.001, -0.00009])


def s_of_a(a):
    la = np.log(a)
    lan = np.log(_A_ANCH)
    ls = np.log(-_S_ANCH)
    ls_i = np.interp(la, lan, ls)
    return -np.exp(ls_i) * ARCSEC_TO_RAD


def ecl_to_frame(inc, Om, om, D, XPR, YPR):
    h = np.array([np.sin(inc) * np.sin(Om), -np.sin(inc) * np.cos(Om),
                  np.cos(inc)])
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


def frame_to_ecl_arrays(I, wp, Op, E):
    """(I, wp, Op) arrays in a frame with basis matrix E (cols XPR,YPR,D)
    -> (i, Om, om) in ecliptic coordinates, degrees. Vectorised."""
    sI = np.sin(I)
    cI = np.cos(I)
    sOp = np.sin(Op)
    cOp = np.cos(Op)
    hx = sI * sOp
    hy = -sI * cOp
    hz = cI
    nx = cOp
    ny = sOp
    # h x n in frame coords = (-cI sinOp, cI cosOp, sI)
    cx = -cI * sOp
    cy = cI * cOp
    cz = sI
    cwp = np.cos(wp)
    swp = np.sin(wp)
    px = cwp * nx + swp * cx
    py = cwp * ny + swp * cy
    pz = swp * cz
    def conv(vx, vy, vz):
        return (E[0, 0] * vx + E[0, 1] * vy + E[0, 2] * vz,
                E[1, 0] * vx + E[1, 1] * vy + E[1, 2] * vz,
                E[2, 0] * vx + E[2, 1] * vy + E[2, 2] * vz)
    hex_, hey, hez = conv(hx, hy, hz)
    pex, pey, pez = conv(px, py, pz)
    i_ecl = np.arccos(np.clip(hez, -1, 1))
    si = np.sin(i_ecl)
    si_safe = np.where(si < 1e-12, 1e-12, si)
    nnx = -hey / si_safe
    nny = hex_ / si_safe
    kx = -hez * nny
    ky = hez * nnx
    kz = hex_ * nny - hey * nnx
    om_ecl = np.arctan2(pex * kx + pey * ky + pez * kz,
                        pex * nnx + pey * nny)
    Om_ecl = np.arctan2(nny, nnx)
    om_ecl = np.where(si < 0.01, np.nan, om_ecl)  # ill-defined at tiny i
    return (np.degrees(i_ecl), np.degrees(Om_ecl % (2 * np.pi)),
            np.degrees(om_ecl % (2 * np.pi)))


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
        dh = Acoef * ((6.0 * c / Gc) * (2 + 3 * u)
                      - (30.0 * u * c / Gc) * cg) + sa
        return [dG, dg, 0.0, dh]
    return rhs


def run_traj(a, e0, I0, wp0, Acoef, sa, tmax=TMAX, npts=NPTS):
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
    g = sol.y[1]
    Hh = sol.y[2]
    h = sol.y[3]
    e = np.sqrt(np.clip(1 - (G / L) ** 2, 0, 1))
    return tev, e, g, h, G, Hh, L


def rayleigh(angles_rad):
    """Rayleigh test for a unimodal direction. Returns (Rbar, z, p, n)."""
    ang = np.asarray(angles_rad, dtype=float)
    ang = ang[np.isfinite(ang)]
    n = len(ang)
    if n == 0:
        return np.nan, np.nan, np.nan, 0
    Rbar = abs(np.exp(1j * ang).mean())
    z = n * Rbar ** 2
    p = np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * n))
    return Rbar, z, max(p, 0.0), n


def find_episodes(mask):
    eps = []
    i = 0
    n = len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j + 1 < n and mask[j + 1]:
                j += 1
            eps.append((i, j))
            i = j + 1
        else:
            i += 1
    return eps


# ---------------- load calibration + population ----------------
with open("crosser_simulation_results.json") as f:
    prev = json.load(f)
A300_MIN = prev["A300_MIN"]
A300_9 = prev["A300_9"]
print(f"A300_MIN = {A300_MIN:.4e}, A300_9 = {A300_9:.4e}", flush=True)

D_CAR = axis_from_lonlat(30.0, 75.0)
D_P9 = axis_from_lonlat(10.0, 74.0)
sep = np.degrees(np.arccos(np.clip(np.dot(D_CAR, D_P9), -1, 1)))
print(f"carrier axis vs P9 pole angular separation: {sep:.2f} deg", flush=True)

rng = np.random.default_rng(42)
N = 200
a_s = rng.uniform(100.0, 400.0, N)
q_s = rng.uniform(30.0, 80.0, N)
i_s = np.radians(rng.uniform(0.0, 40.0, N))
Om_s = rng.uniform(0, 2 * np.pi, N)
om_s = rng.uniform(0, 2 * np.pi, N)
e0_s = 1.0 - q_s / a_s
assert np.allclose(a_s, prev["a"]) and np.allclose(q_s, prev["q0"]), \
    "population mismatch with saved run"
print("population reproduced (seed 42, matches saved run)", flush=True)

# ---------------- round-trip test of the inverse transform ----------------
for tag, D in [("carrier", D_CAR), ("P9", D_P9)]:
    XPR, YPR, E = frame_basis(D)
    maxerr = 0.0
    for k in range(N):
        I, wp, Op = ecl_to_frame(i_s[k], Om_s[k], om_s[k], D, XPR, YPR)
        i2, Om2, om2 = frame_to_ecl_arrays(np.array([I]), np.array([wp]),
                                           np.array([Op]), E)
        d = lambda x, y: abs((x - y + 180.0) % 360.0 - 180.0)
        maxerr = max(maxerr, d(np.degrees(i_s[k]), i2[0]),
                     d(np.degrees(Om_s[k]), Om2[0]),
                     d(np.degrees(om_s[k]), om2[0]))
    print(f"round-trip {tag} frame: max angle error {maxerr:.2e} deg", flush=True)
    assert maxerr < 1e-6, "inverse transform FAILED"

# ---------------- runs ----------------
runs = [("carrier_Smin", D_CAR, A300_MIN, "S_min"),
        ("carrier_S9", D_CAR, A300_9, "S_9"),
        ("P9sec", D_P9, A300_9, None)]
out = {}
for tag, D, A300, prev_key in runs:
    XPR, YPR, E = frame_basis(D)
    I0_s = np.empty(N)
    wp0_s = np.empty(N)
    for k in range(N):
        I0_s[k], wp0_s[k], _ = ecl_to_frame(i_s[k], Om_s[k], om_s[k], D, XPR, YPR)
    qmins = np.empty(N)
    duty = np.empty(N)
    per_part = []
    episodes = []          # dicts: a, dur_Myr, qmin, wp_qmin, om_qmin, i_qmin, censored
    snap_a, snap_q, snap_i, snap_wp, snap_om = [], [], [], [], []
    for k in range(N):
        a = a_s[k]
        Acoef = A300 * (rho(a) / 0.45) * (a / 300.0) ** 1.5
        sa = s_of_a(a)
        tev, e, g, h, G, Hh, L = run_traj(a, e0_s[k], I0_s[k], wp0_s[k], Acoef, sa)
        q = a * (1.0 - e)
        Iarr = np.arccos(np.clip(Hh / np.maximum(G, 1e-12), -1, 1))
        i_ecl, Om_ecl, om_ecl = frame_to_ecl_arrays(Iarr, g, h, E)
        wp_deg = np.degrees(g) % 360.0
        mask = q < 30.0
        qmins[k] = q.min()
        duty[k] = mask.mean()
        eps = find_episodes(mask)
        for (i0, i1) in eps:
            seg = slice(i0, i1 + 1)
            kmin = i0 + int(np.argmin(q[seg]))
            episodes.append(dict(a=float(a), dur_Myr=float((i1 - i0 + 1) * DT / 1e6),
                                 qmin=float(q[kmin]),
                                 wp_qmin=float(wp_deg[kmin]),
                                 om_qmin=(float(om_ecl[kmin])
                                          if np.isfinite(om_ecl[kmin]) else None),
                                 i_qmin=float(i_ecl[kmin]),
                                 censored=bool(i1 == NPTS - 1)))
        per_part.append(dict(a=float(a), e0=float(e0_s[k]), I0=float(np.degrees(I0_s[k])),
                             qmin=float(qmins[k]), duty=float(duty[k]),
                             n_ep=len(eps)))
        if mask.any():
            snap_a.extend([a] * int(mask.sum()))
            snap_q.extend(q[mask].tolist())
            snap_i.extend(i_ecl[mask].tolist())
            snap_wp.extend(wp_deg[mask].tolist())
            snap_om.extend([v for v in om_ecl[mask].tolist() if np.isfinite(v)])
        if (k + 1) % 50 == 0:
            print(f"  {tag}: {k+1}/{N} particles", flush=True)
    cross = qmins < 30.0
    out[tag] = dict(n_cross=int(cross.sum()), qmins=qmins.tolist(),
                    duty=duty.tolist(), per_part=per_part, episodes=episodes,
                    snap_a=snap_a, snap_q=snap_q, snap_i=snap_i,
                    snap_wp=snap_wp, snap_om=snap_om)
    # reproduction check vs the original run
    if prev_key is not None:
        qm_old = np.array(prev["results"][prev_key]["qmins"])
        print(f"  {tag}: crossers {cross.sum()} (original: "
              f"{prev['results'][prev_key]['n_cross']}); "
              f"max |qmin_new - qmin_old| = {np.max(np.abs(qmins - qm_old)):.3f} AU",
              flush=True)

with open("crosser_distributions_results.json", "w") as f:
    json.dump(out, f)
print("saved crosser_distributions_results.json", flush=True)


# ---------------- statistics ----------------
def pct(x, p):
    return float(np.percentile(x, p))

def summ(x, name):
    x = np.asarray(x, dtype=float)
    print(f"  {name}: n={len(x)} median={pct(x,50):.2f} "
          f"[Q1={pct(x,25):.2f}, Q3={pct(x,75):.2f}] range=[{x.min():.2f}, {x.max():.2f}]",
          flush=True)

print("\n================ STATISTICS ================", flush=True)
res = {}
for tag in ["carrier_Smin", "carrier_S9", "P9sec"]:
    d = out[tag]
    qmins = np.array(d["qmins"])
    duty = np.array(d["duty"])
    cross = qmins < 30.0
    a_cross = a_s[cross]
    a_non = a_s[~cross]
    eps = d["episodes"]
    dur = np.array([e["dur_Myr"] for e in eps])
    wp_ep = np.array([e["wp_qmin"] for e in eps])
    om_ep = np.array([e["om_qmin"] for e in eps if e["om_qmin"] is not None])
    snap_q = np.array(d["snap_q"])
    snap_a = np.array(d["snap_a"])
    snap_i = np.array(d["snap_i"])
    snap_wp = np.array(d["snap_wp"])
    print(f"\n--- {tag}: crossers {cross.sum()}/200 ---", flush=True)
    summ(a_cross, "a ever-crossed (AU)")
    bins = np.histogram(a_cross, bins=[100, 150, 200, 250, 300, 350, 400])[0]
    print(f"  a ever-crossed 50-AU bins [100-150..350-400]: {bins.tolist()} "
          f"(uniform expectation {len(a_cross)/6:.1f}/bin)", flush=True)
    ks = stats.ks_2samp(a_cross, a_non)
    print(f"  KS a_cross vs a_non: D={ks.statistic:.3f} p={ks.pvalue:.3f}", flush=True)
    ks_u = stats.kstest(a_cross, "uniform", args=(100, 300))
    print(f"  KS a_cross vs U[100,400]: D={ks_u.statistic:.3f} p={ks_u.pvalue:.3f}",
          flush=True)
    summ(snap_a, "a snapshot-weighted (AU)")
    summ(snap_q, "q snapshot (AU)")
    hq = np.histogram(snap_q, bins=[14, 18, 20, 22, 24, 26, 28, 30])[0]
    print(f"  q snapshot 2-AU bins [14-18..28-30]: {hq.tolist()}", flush=True)
    print(f"  snapshot q<20 AU: {(snap_q < 20).mean()*100:.1f}%; "
          f"q<25 AU: {(snap_q < 25).mean()*100:.1f}% "
          f"(flat U[16,30] would give 28.6% / 64.3%)", flush=True)
    ks_q = stats.kstest(snap_q, "uniform", args=(16, 14))
    print(f"  KS snapshot q vs U[16,30] (P9-flat): D={ks_q.statistic:.3f} "
          f"p={ks_q.pvalue:.2e}", flush=True)
    summ(snap_i, "i snapshot (deg)")
    hi = np.histogram(snap_i, bins=[0, 5, 10, 15, 20, 25, 30, 35, 40])[0]
    print(f"  i snapshot 5-deg bins [0-5..35-40]: {hi.tolist()}", flush=True)
    dc = duty[cross]
    summ(dc * 100, "duty cycle per crosser (%)")
    print(f"  episodes: n={len(eps)}; per crosser mean={len(eps)/max(cross.sum(),1):.2f}; "
          f"right-censored (touching t=4.5Gyr): "
          f"{sum(1 for e in eps if e['censored'])}", flush=True)
    summ(dur, "single-episode duration (Myr)")
    R2, z2, p2, n2 = rayleigh(np.radians(2 * wp_ep))
    frac_win = np.mean(((wp_ep >= 60) & (wp_ep <= 120)) |
                       ((wp_ep >= 240) & (wp_ep <= 300)))
    print(f"  omega' at episode q-min (n={n2}): doubled-angle Rayleigh "
          f"Rbar2={R2:.3f} z={z2:.1f} p={p2:.2e}; "
          f"fraction in [60,120]u[240,300] = {frac_win*100:.0f}% (uniform 33%)",
          flush=True)
    hw = np.histogram(snap_wp, bins=np.arange(0, 361, 30))[0]
    print(f"  omega' snapshot 30-deg bins: {hw.tolist()} "
          f"(uniform {len(snap_wp)/12:.0f}/bin)", flush=True)
    R, z, p, n = rayleigh(np.radians(om_ep))
    print(f"  omega (ecliptic) at episode q-min (n={n}): Rayleigh Rbar={R:.3f} "
          f"z={z:.2f} p={p:.3f}", flush=True)
    ho = np.histogram(np.array(d["snap_om"]), bins=np.arange(0, 361, 30))[0]
    print(f"  omega snapshot 30-deg bins: {ho.tolist()} "
          f"(uniform {len(d['snap_om'])/12:.0f}/bin)", flush=True)
    res[tag] = dict(n_cross=int(cross.sum()),
                    a_cross_med=pct(a_cross, 50), a_cross_q1=pct(a_cross, 25),
                    a_cross_q3=pct(a_cross, 75),
                    snap_q_med=pct(snap_q, 50), snap_q_lt20=float((snap_q < 20).mean()),
                    snap_q_lt25=float((snap_q < 25).mean()),
                    ks_q_D=float(ks_q.statistic), ks_q_p=float(ks_q.pvalue),
                    duty_med=pct(dc * 100, 50),
                    n_ep=len(eps), dur_med=pct(dur, 50), dur_max=float(dur.max()),
                    wp_Rbar2=float(R2), wp_p2=float(p2), wp_frac_win=float(frac_win),
                    om_Rbar=float(R), om_p=float(p))

# cross-run comparisons
print("\n--- S_min vs S_9 (carrier) ---", flush=True)
for key, arr in [("a ever-crossed", "a"), ("duty", "duty"), ("episode duration", "dur")]:
    if arr == "a":
        x = a_s[np.array(out["carrier_Smin"]["qmins"]) < 30]
        y = a_s[np.array(out["carrier_S9"]["qmins"]) < 30]
    elif arr == "duty":
        x = np.array(out["carrier_Smin"]["duty"])[np.array(out["carrier_Smin"]["qmins"]) < 30]
        y = np.array(out["carrier_S9"]["duty"])[np.array(out["carrier_S9"]["qmins"]) < 30]
    else:
        x = np.array([e["dur_Myr"] for e in out["carrier_Smin"]["episodes"]])
        y = np.array([e["dur_Myr"] for e in out["carrier_S9"]["episodes"]])
    ks = stats.ks_2samp(x, y)
    print(f"  KS {key}: D={ks.statistic:.3f} p={ks.pvalue:.3f}", flush=True)

print("\n--- P9sec vs carrier_S9 ---", flush=True)
for key, arr in [("a ever-crossed", "a"), ("duty", "duty"), ("episode duration", "dur")]:
    if arr == "a":
        x = a_s[np.array(out["P9sec"]["qmins"]) < 30]
        y = a_s[np.array(out["carrier_S9"]["qmins"]) < 30]
    elif arr == "duty":
        x = np.array(out["P9sec"]["duty"])[np.array(out["P9sec"]["qmins"]) < 30]
        y = np.array(out["carrier_S9"]["duty"])[np.array(out["carrier_S9"]["qmins"]) < 30]
    else:
        x = np.array([e["dur_Myr"] for e in out["P9sec"]["episodes"]])
        y = np.array([e["dur_Myr"] for e in out["carrier_S9"]["episodes"]])
    ks = stats.ks_2samp(x, y)
    print(f"  KS {key}: D={ks.statistic:.3f} p={ks.pvalue:.3f}", flush=True)

# P9-control episodes seen in the CARRIER frame: does the P9 secular
# crosser sample also look concentrated in omega' around the carrier axis?
print("\n--- P9sec episodes re-expressed in the carrier frame ---", flush=True)
D = D_CAR
XPRc, YPRc, Ec = frame_basis(D)
# Re-run the P9 crossers only and convert each episode's q-minimum state
# (ecliptic i/Om/om) into the carrier frame, to test whether a P9-made
# crosser sample would ALSO look concentrated in omega' around the carrier
# axis (the two axes are only ~5 deg apart, so the answer matters).
Xp, Yp, Ep = frame_basis(D_P9)
qmins_p9 = np.array(out["P9sec"]["qmins"])
wp_car_ep = []
n_skip_nan = 0
for k in range(N):
    if qmins_p9[k] >= 30.0:
        continue
    a = a_s[k]
    I0pk, wp0pk, _ = ecl_to_frame(i_s[k], Om_s[k], om_s[k], D_P9, Xp, Yp)
    Acoef = A300_9 * (rho(a) / 0.45) * (a / 300.0) ** 1.5
    tev, e, g, h, G, Hh, L = run_traj(a, e0_s[k], I0pk, wp0pk, Acoef, s_of_a(a))
    q = a * (1.0 - e)
    Iarr = np.arccos(np.clip(Hh / np.maximum(G, 1e-12), -1, 1))
    i_ecl, Om_ecl, om_ecl = frame_to_ecl_arrays(Iarr, g, h, Ep)
    mask = q < 30.0
    for (i0, i1) in find_episodes(mask):
        kmin = i0 + int(np.argmin(q[i0:i1 + 1]))
        if not np.isfinite(om_ecl[kmin]):
            n_skip_nan += 1
            continue
        Ic, wpc, _ = ecl_to_frame(np.radians(i_ecl[kmin]), np.radians(Om_ecl[kmin]),
                                  np.radians(om_ecl[kmin]), D_CAR, XPRc, YPRc)
        wp_car_ep.append(np.degrees(wpc) % 360.0)
wp_car_ep = np.array(wp_car_ep)
R2c, z2c, p2c, n2c = rayleigh(np.radians(2 * wp_car_ep))
print(f"  P9sec episodes in CARRIER frame omega' (n={n2c}): doubled-angle "
      f"Rayleigh Rbar2={R2c:.3f} z={z2c:.1f} p={p2c:.2e}", flush=True)
res["P9sec_carrierframe_wp"] = dict(Rbar2=float(R2c), p=float(p2c), n=int(n2c))

with open("crosser_distributions_summary.json", "w") as f:
    json.dump(res, f, indent=1)
print("\nDONE")
