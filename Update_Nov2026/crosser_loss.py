#!/usr/bin/env python3
"""Path (a2): add a loss process to the carrier crosser simulation.

Question (convergence.md): can adding a loss process bring the
same-population time-averaged instantaneous crosser fraction down from
8.4-9.1% to Batygin's ~3%?

Loss model (stated assumption): a Neptune-crosser (q<30 AU) is removed
from the population by a memoryless (Poisson) process with rate 1/tau
while it is crossing -- i.e. at each 5-Myr sample step spent at q<30 it
is removed with probability p = 1 - exp(-DT/tau). tau is the
characteristic clearing time; Batygin's N-body clears crossers in
~100 Myr. tau = 50 / 100 / 200 Myr are run, plus tau = infinity
(no loss, control -- must reproduce crosser_samepop.py exactly).
Loss acts only on crossers; non-crossers are never removed. Removed
particles leave the system (no replenishment -- stated limitation: the
model then has no steady state, the crosser population decays).

Metrics (same as before, plus loss diagnostics):
  * ever-crossed fraction (unchanged by loss -- removal happens after
    first crossing; verified, not assumed)
  * time-averaged instantaneous fraction among survivors:
    (sum over sample times of alive crossers) / (sum of alive),
    full 4.5 Gyr window (headline, same window as the no-loss baseline)
    and last-1-Gyr window (Batygin's footprint window)
  * snapshot q distribution (alive snapshots only): median, q<20,
    q<25, min; KS of per-episode qmin vs flat U[16,30]
  * removed fraction, median removal time.
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
DT = TMAX / (NPTS - 1)          # 5 Myr
N = 1000
SEED = 123
TAUS_MYR = [None, 50.0, 100.0, 200.0]   # None = no loss (control)
LAST1_N = 200                    # last 1 Gyr = last 200 sample points


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
    p = alpha + 1.0  # -0.5
    u = rng.uniform(0, 1, n)
    return (amin ** p + u * (amax ** p - amin ** p)) ** (1.0 / p)


def make_population(rng, n, amin, amax):
    a_s = sample_a_powerlaw(rng, n, amin, amax)
    is_scat = rng.uniform(0, 1, n) < 0.70
    q_s = np.empty(n)
    q_s[is_scat] = 30.0 + rng.exponential(3.0, int(is_scat.sum()))
    q_s[is_scat] = np.clip(q_s[is_scat], 30.0, 45.0)
    q_s[~is_scat] = rng.uniform(35.0, 80.0, int((~is_scat).sum()))
    i_s = np.abs(rng.normal(0.0, np.radians(15.0), n))
    Om_s = rng.uniform(0, 2 * np.pi, n)
    om_s = rng.uniform(0, 2 * np.pi, n)
    return a_s, q_s, i_s, Om_s, om_s


# ---------------- load calibration ----------------
with open("crosser_simulation_results.json") as f:
    prev = json.load(f)
A300_MIN = prev["A300_MIN"]
A300_9 = prev["A300_9"]
print(f"A300_MIN = {A300_MIN:.4e}, A300_9 = {A300_9:.4e}; DT = {DT/1e6:.0f} Myr",
      flush=True)

D_CAR = axis_from_lonlat(30.0, 75.0)
XPR, YPR, E = frame_basis(D_CAR)
JIDX = np.arange(NPTS)

results = {}
snap_store = {}
for v_idx, (version, (amin, amax)) in enumerate(
        [("A_100_400", (100.0, 400.0)), ("B_100_5000", (100.0, 5000.0))]):
    rng = np.random.default_rng(SEED)
    a_s, q_s, i_s, Om_s, om_s = make_population(rng, N, amin, amax)
    e0_s = 1.0 - q_s / a_s
    I0_s = np.empty(N); wp0_s = np.empty(N)
    for k in range(N):
        I0_s[k], wp0_s[k], _ = ecl_to_frame(i_s[k], Om_s[k], om_s[k],
                                            D_CAR, XPR, YPR)
    print(f"\n=== {version}: a in [{amin:.0f},{amax:.0f}], N={N} ===", flush=True)
    for s_idx, (tag, A300) in enumerate([("S_min", A300_MIN), ("S_9", A300_9)]):
        t0 = time.time()
        # ---- integrate all particles once, store trajectories ----
        Q = np.empty((N, NPTS)); IE = np.empty((N, NPTS)); WP = np.empty((N, NPTS))
        for k in range(N):
            a = a_s[k]
            Acoef = A300 * (rho(a) / 0.45) * (a / 300.0) ** 1.5
            sa = s_of_a(a)
            tev, e, g, h, G, Hh, L = run_traj(a, e0_s[k], I0_s[k], wp0_s[k],
                                             Acoef, sa)
            Q[k] = a * (1.0 - e)
            Iarr = np.arccos(np.clip(Hh / np.maximum(G, 1e-12), -1, 1))
            i_ecl, Om_ecl, om_ecl = frame_to_ecl_arrays(Iarr, g, h, E)
            IE[k] = i_ecl
            WP[k] = np.degrees(g) % 360.0
            if (k + 1) % 500 == 0:
                print(f"  {version}/{tag}: integrated {k+1}/{N} "
                      f"({time.time()-t0:.0f}s)", flush=True)
        mask = Q < 30.0
        ever = mask.any(axis=1)
        A_mat = np.repeat(a_s[:, None], NPTS, axis=1)
        print(f"  {version}/{tag}: trajectories done ({time.time()-t0:.0f}s); "
              f"EVER-crossed (full traj) {ever.sum()}/{N} = {ever.mean()*100:.1f}%",
              flush=True)
        # ---- apply each loss tau in post-processing ----
        for t_idx, tau_myr in enumerate(TAUS_MYR):
            if tau_myr is None:
                remove_idx = np.full(N, NPTS)
                tau_label = "no_loss"
            else:
                lrng = np.random.default_rng(7000 + v_idx * 1000 + s_idx * 100 + t_idx)
                U = lrng.uniform(size=(N, NPTS))
                p = 1.0 - np.exp(-DT / (tau_myr * 1e6))
                hazard = mask & (U < p)
                has = hazard.any(axis=1)
                remove_idx = np.where(has, hazard.argmax(axis=1), NPTS)
                tau_label = f"tau{int(tau_myr)}"
            alive = JIDX[None, :] < remove_idx[:, None]
            cross_alive = mask & alive
            # headline: (sum of alive crossers) / (sum of alive), full window
            inst_full = cross_alive.sum() / alive.sum()
            inst_last1 = (cross_alive[:, -LAST1_N:].sum()
                          / alive[:, -LAST1_N:].sum())
            alive_j = alive.sum(axis=0)
            cross_j = cross_alive.sum(axis=0)
            inst_mean = float(np.mean(cross_j / np.maximum(alive_j, 1)))
            # Batygin aperture a in [100,1000]
            sel = a_s <= 1000.0
            inst_full_a1000 = (cross_alive[sel].sum() / alive[sel].sum())
            inst_last1_a1000 = (cross_alive[sel][:, -LAST1_N:].sum()
                                / alive[sel][:, -LAST1_N:].sum())
            # snapshots (alive only)
            sq = Q[cross_alive]
            sa_snap = A_mat[cross_alive]
            si_snap = IE[cross_alive]
            # episodes on truncated trajectories
            episodes = []
            for k in range(N):
                ri = int(remove_idx[k])
                if ri == 0:
                    continue
                eps = find_episodes(mask[k, :ri])
                for (i0, i1) in eps:
                    kmin = i0 + int(np.argmin(Q[k, i0:i1 + 1]))
                    episodes.append(dict(a=float(a_s[k]),
                                         dur_Myr=float((i1 - i0 + 1) * DT / 1e6),
                                         qmin=float(Q[k, kmin]),
                                         wp_qmin=float(WP[k, kmin]),
                                         i_qmin=float(IE[k, kmin]),
                                         censored=bool(i1 == ri - 1 and ri < NPTS)))
            qm_ep = np.array([e["qmin"] for e in episodes]) if episodes else np.array([])
            ks_flat = (stats.kstest(qm_ep, "uniform", args=(16, 14))
                       if len(qm_ep) else None)
            n_removed = int((remove_idx < NPTS).sum())
            rem_times = remove_idx[remove_idx < NPTS] * DT / 1e6 if n_removed else np.array([])
            key = f"{version}_{tag}_{tau_label}"
            results[key] = dict(
                n=N, version=version, tag=tag, tau_Myr=tau_myr,
                n_ever=int(ever.sum()), frac_ever=float(ever.mean()),
                inst_full=float(inst_full), inst_last1=float(inst_last1),
                inst_mean=float(inst_mean),
                inst_full_a1000=float(inst_full_a1000),
                inst_last1_a1000=float(inst_last1_a1000),
                snap_q_median=float(np.median(sq)) if len(sq) else None,
                snap_q_lt20=float((sq < 20).mean()) if len(sq) else None,
                snap_q_lt25=float((sq < 25).mean()) if len(sq) else None,
                snap_q_min=float(sq.min()) if len(sq) else None,
                n_snap=int(len(sq)), n_episodes=int(len(episodes)),
                ks_flat_D=(float(ks_flat.statistic) if ks_flat else None),
                ks_flat_p=(float(ks_flat.pvalue) if ks_flat else None),
                n_removed=n_removed, frac_removed=float(n_removed / N),
                rem_time_median_Myr=(float(np.median(rem_times)) if n_removed else None),
                rem_time_mean_Myr=(float(np.mean(rem_times)) if n_removed else None),
                remove_idx=remove_idx.tolist(),
                qmin_full=Q.min(axis=1).tolist(),
                duty_full=mask.mean(axis=1).tolist())
            snap_store[key + "_a"] = sa_snap
            snap_store[key + "_q"] = sq
            snap_store[key + "_i"] = si_snap
            print(f"  {key}: EVER {ever.mean()*100:.1f}% | inst_full "
                  f"{inst_full*100:.2f}% | inst_last1Gyr {inst_last1*100:.2f}% | "
                  f"a<=1000 inst_full {inst_full_a1000*100:.2f}% | "
                  f"snap q med {np.median(sq):.1f} AU, q<20 {(sq<20).mean()*100:.1f}%, "
                  f"q<25 {(sq<25).mean()*100:.1f}%, min {sq.min():.1f} | "
                  f"ep n={len(episodes)}"
                  + (f", KS D={ks_flat.statistic:.3f} p={ks_flat.pvalue:.2e}" if ks_flat else "")
                  + f" | removed {n_removed}/{N} ({n_removed/N*100:.1f}%)"
                  + (f", rem t med {np.median(rem_times):.0f} Myr" if n_removed else ""),
                  flush=True)
        with open("crosser_loss_results.json", "w") as f:
            json.dump(results, f)
        np.savez_compressed("crosser_loss_snapshots.npz", **snap_store)
        print("  saved crosser_loss_results.json + crosser_loss_snapshots.npz "
              "(incremental)", flush=True)

print("\nDONE")
