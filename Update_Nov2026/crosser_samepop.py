#!/usr/bin/env python3
"""Same-population re-run of the carrier crosser simulation (task 1, 2026-11-17).

Question: does the carrier's q-shape (hugging 30 AU, steep drop, ~19 AU
cutoff) survive when the initial population is Batygin's scattered-disc
population instead of the uniform one used before?

Population (reconstructed cluster_2-like, NOT the actual cluster_2 data,
which is not publicly released -- see note):
  * a: dN/da ~ a^(-3/2) (Fokker-Planck steady state for Neptune scattering;
    also the a-distribution in Batygin+2024's own validation exercise).
  * q: two-component mixture --
      70% scattering: q = 30 + Exp(mean 3 AU), truncated to [30, 45]
          (Neptune-scattered objects hug q ~ 30-35 AU);
      30% detached:  q ~ Uniform[35, 80] AU
          (cluster-raised perihelia, broad; Nesvorny+2023 cluster2).
    The 70/30 split is a judgement call, stated as such.
  * i: half-normal sigma = 15 deg (Batygin+2024 validation exercise:
    "Gaussian distribution of inclinations with a dispersion of 15 deg").
  * Omega, omega: uniform [0, 360).
Two versions:
  * A (clean): a in [100, 400] AU -- inside the measured rho(a) range.
  * B (full):  a in [100, 5000] AU (Batygin's range); rho(a) = 0.45 for
    a > 400 AU (flat extrapolation -- an ASSUMPTION, stated in the note).
Everything else is inherited unchanged from crosser_simulation.py /
crosser_distributions.py: secular quadrupole Hamiltonian, calibration
A300 (read from crosser_simulation_results.json), rho(a) profile,
4.5 Gyr integration, giants' background nodal rate s(a) (log-log interp
through the [50,80,200,400] AU anchors; clamped beyond 400 AU -- s(a)
only enters dOm'/dt and does not affect (e,I,w') in the quadrupole
Hamiltonian, so this is a bookkeeping term).

Outputs (saved incrementally after each run):
  per-particle q_min, duty cycle; per-episode q_min/duration/wp';
  crossing snapshots (a, q, i_ecl) for the xi-test and q-shape.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy import stats
import json, sys, time

MU = 4.0 * np.pi ** 2
M2S2 = 1.0 / 2.24720e7
ARCSEC_TO_RAD = 4.84813681109536e-6
ZEC = np.array([0.0, 0.0, 1.0])
TMAX = 4.5e9
NPTS = 901
DT = TMAX / (NPTS - 1)
N = 1000
SEED = 123


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
    """Measured profile to 400 AU; flat 0.45 beyond (assumption, Run B)."""
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
    sI = np.sin(I); cI = np.cos(I); sOp = np.sin(Op); cOp = np.cos(Op)
    hx = sI * sOp; hy = -sI * cOp; hz = cI
    nx = cOp; ny = sOp
    cx = -cI * sOp; cy = cI * cOp; cz = sI
    cwp = np.cos(wp); swp = np.sin(wp)
    px = cwp * nx + swp * cx; py = cwp * ny + swp * cy; pz = swp * cz
    def conv(vx, vy, vz):
        return (E[0, 0] * vx + E[0, 1] * vy + E[0, 2] * vz,
                E[1, 0] * vx + E[1, 1] * vy + E[1, 2] * vz,
                E[2, 0] * vx + E[2, 1] * vy + E[2, 2] * vz)
    hex_, hey, hez = conv(hx, hy, hz)
    pex, pey, pez = conv(px, py, pz)
    i_ecl = np.arccos(np.clip(hez, -1, 1))
    si = np.sin(i_ecl)
    si_safe = np.where(si < 1e-12, 1e-12, si)
    nnx = -hey / si_safe; nny = hex_ / si_safe
    kx = -hez * nny; ky = hez * nnx; kz = hex_ * nny - hey * nnx
    om_ecl = np.arctan2(pex * kx + pey * ky + pez * kz,
                        pex * nnx + pey * nny)
    Om_ecl = np.arctan2(nny, nnx)
    om_ecl = np.where(si < 0.01, np.nan, om_ecl)
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
        cg = np.cos(2 * g); sg = np.sin(2 * g)
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
    G = sol.y[0]; g = sol.y[1]; Hh = sol.y[2]; h = sol.y[3]
    e = np.sqrt(np.clip(1 - (G / L) ** 2, 0, 1))
    return tev, e, g, h, G, Hh, L


def find_episodes(mask):
    eps = []
    i = 0; n = len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j + 1 < n and mask[j + 1]:
                j += 1
            eps.append((i, j)); i = j + 1
        else:
            i += 1
    return eps


def sample_a_powerlaw(rng, n, amin, amax, alpha=-1.5):
    """dN/da ~ a^alpha; inverse-CDF sampling."""
    # CDF: F(a) = (a^(alpha+1) - amin^(alpha+1)) / (amax^(alpha+1) - amin^(alpha+1))
    p = alpha + 1.0  # -0.5
    u = rng.uniform(0, 1, n)
    return (amin ** p + u * (amax ** p - amin ** p)) ** (1.0 / p)


def make_population(rng, n, amin, amax):
    a_s = sample_a_powerlaw(rng, n, amin, amax)
    # q mixture: 70% scattering (30 + Exp(3), truncated [30,45]);
    #            30% detached (Uniform[35,80]).
    is_scat = rng.uniform(0, 1, n) < 0.70
    q_s = np.empty(n)
    q_s[is_scat] = 30.0 + rng.exponential(3.0, int(is_scat.sum()))
    q_s[is_scat] = np.clip(q_s[is_scat], 30.0, 45.0)
    q_s[~is_scat] = rng.uniform(35.0, 80.0, int((~is_scat).sum()))
    # i: half-normal sigma=15 deg
    i_s = np.abs(rng.normal(0.0, np.radians(15.0), n))
    Om_s = rng.uniform(0, 2 * np.pi, n)
    om_s = rng.uniform(0, 2 * np.pi, n)
    return a_s, q_s, i_s, Om_s, om_s


# ---------------- load calibration ----------------
with open("crosser_simulation_results.json") as f:
    prev = json.load(f)
A300_MIN = prev["A300_MIN"]
A300_9 = prev["A300_9"]
print(f"A300_MIN = {A300_MIN:.4e}, A300_9 = {A300_9:.4e}", flush=True)

D_CAR = axis_from_lonlat(30.0, 75.0)
XPR, YPR, E = frame_basis(D_CAR)

results = {}
for version, (amin, amax) in [("A_100_400", (100.0, 400.0)),
                              ("B_100_5000", (100.0, 5000.0))]:
    rng = np.random.default_rng(SEED)
    a_s, q_s, i_s, Om_s, om_s = make_population(rng, N, amin, amax)
    e0_s = 1.0 - q_s / a_s
    assert np.all(e0_s > 0) and np.all(e0_s < 1), "bad e0"
    I0_s = np.empty(N); wp0_s = np.empty(N)
    for k in range(N):
        I0_s[k], wp0_s[k], _ = ecl_to_frame(i_s[k], Om_s[k], om_s[k],
                                            D_CAR, XPR, YPR)
    print(f"\n=== {version}: a in [{amin:.0f},{amax:.0f}], N={N} ===", flush=True)
    print(f"  a median {np.median(a_s):.0f} AU, q0 median {np.median(q_s):.1f} AU "
          f"(q0<35: {(q_s<35).mean()*100:.0f}%, q0<40: {(q_s<40).mean()*100:.0f}%), "
          f"i median {np.degrees(np.median(i_s)):.1f} deg, "
          f"e0 median {np.median(e0_s):.3f}", flush=True)
    for tag, A300 in [("S_min", A300_MIN), ("S_9", A300_9)]:
        t0 = time.time()
        qmins = np.empty(N); duty = np.empty(N)
        per_part = []
        episodes = []
        snap_a, snap_q, snap_i = [], [], []
        for k in range(N):
            a = a_s[k]
            Acoef = A300 * (rho(a) / 0.45) * (a / 300.0) ** 1.5
            sa = s_of_a(a)
            tev, e, g, h, G, Hh, L = run_traj(a, e0_s[k], I0_s[k], wp0_s[k],
                                             Acoef, sa)
            q = a * (1.0 - e)
            Iarr = np.arccos(np.clip(Hh / np.maximum(G, 1e-12), -1, 1))
            i_ecl, Om_ecl, om_ecl = frame_to_ecl_arrays(Iarr, g, h, E)
            wp_deg = np.degrees(g) % 360.0
            mask = q < 30.0
            qmins[k] = q.min(); duty[k] = mask.mean()
            eps = find_episodes(mask)
            for (i0, i1) in eps:
                kmin = i0 + int(np.argmin(q[i0:i1 + 1]))
                episodes.append(dict(a=float(a),
                                     dur_Myr=float((i1 - i0 + 1) * DT / 1e6),
                                     qmin=float(q[kmin]),
                                     wp_qmin=float(wp_deg[kmin]),
                                     i_qmin=float(i_ecl[kmin]),
                                     censored=bool(i1 == NPTS - 1)))
            per_part.append(dict(a=float(a), q0=float(q_s[k]),
                                 e0=float(e0_s[k]),
                                 I0=float(np.degrees(I0_s[k])),
                                 qmin=float(qmins[k]), duty=float(duty[k]),
                                 n_ep=len(eps)))
            if mask.any():
                snap_a.extend([a] * int(mask.sum()))
                snap_q.extend(q[mask].tolist())
                snap_i.extend(i_ecl[mask].tolist())
            if (k + 1) % 250 == 0:
                print(f"  {version}/{tag}: {k+1}/{N} ({time.time()-t0:.0f}s)",
                      flush=True)
        cross = qmins < 30.0
        frac = cross.mean()
        sig = np.sqrt(frac * (1 - frac) / N) if 0 < frac < 1 else 1.0 / N
        inst_timeavg = duty.mean()
        key = f"{version}_{tag}"
        results[key] = dict(n=N, n_cross=int(cross.sum()), frac=float(frac),
                            sigma=float(sig), inst_timeavg=float(inst_timeavg),
                            qmins=qmins.tolist(), duty=duty.tolist(),
                            per_part=per_part, episodes=episodes,
                            snap_a=snap_a, snap_q=snap_q, snap_i=snap_i,
                            a=a_s.tolist(), q0=q_s.tolist())
        # headline statistics
        sq = np.array(snap_q)
        qm_ep = np.array([e["qmin"] for e in episodes])
        ks_flat = stats.kstest(qm_ep, "uniform", args=(16, 14)) if len(qm_ep) else None
        print(f"  {key}: EVER-crossed {cross.sum()}/{N} = {frac*100:.1f}% +- {sig*100:.1f}%; "
              f"time-avg instantaneous {inst_timeavg*100:.1f}%; "
              f"snapshot q median {np.median(sq):.1f} AU, q<20 {(sq<20).mean()*100:.1f}%, "
              f"q<25 {(sq<25).mean()*100:.1f}%, min {sq.min():.1f} AU; "
              f"episodes n={len(episodes)}; "
              f"KS per-episode qmin vs flat U[16,30]: "
              f"D={ks_flat.statistic:.3f} p={ks_flat.pvalue:.2e}" if ks_flat else "",
              flush=True)
        print(f"    ({time.time()-t0:.0f}s)", flush=True)
        with open("crosser_samepop_results.json", "w") as f:
            json.dump(results, f)
        print("  saved crosser_samepop_results.json (incremental)", flush=True)

print("\nDONE")
