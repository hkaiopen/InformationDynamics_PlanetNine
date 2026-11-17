#!/usr/bin/env python3
"""D diffusion + u replenishment test for the carrier crosser simulation.

Question (user, 2026-11-17): the verdict "crossers need a real scattering
source" was K(x)-only -- the carrier simulation used ONLY the iK(x)Psi
term of the generalized GL equation. The D and u terms were never used.
Can D (diffusion) flatten q and can u (source/replenishment) restore a
steady state, so the carrier can do the crosser job after all?

GL equation (id_principles_v4.txt, eq. 1):
    dPsi/dt = D nabla^2 Psi - kappa(|Psi|^2 - rho_target) Psi + iK(x) Psi + u
  * D nabla^2 Psi: "smooths the field and connects neighbouring regions"
    (diffusion).  * u: "an external control input", output y = |Psi|^2.
  kappa (nonlinear compression) is left aside for now (task instruction).

Implementation (stated assumptions -- these are MODELLING choices, not
derivations, same status as the K(x) mapping):
  * D (diffusion): a random walk on q (perihelion distance), applied on
    top of the deterministic secular drift. At each 5-Myr sample step,
    q <- q + (secular change from the donor trajectory) + N(0, sigma_q),
    sigma_q = sqrt(2 * Dq * DT), Dq in AU^2/Myr (free parameter).
    Reflecting boundary at q = 0 (q <- |q|). The map "diffusion in the
    information field" -> "diffusion in q" is a modelling assumption;
    in Fokker-Planck language scattering IS phase-space diffusion, which
    is why D is the candidate for the scattering-like channel. Diffusion
    acts on ALL particles at ALL times (unlike loss, which acts only on
    crossers): it both spreads the q of existing crossers (flattening)
    and pushes non-crossers below q = 30 (production of new crossers).
    Feedback of the diffused q onto the secular evolution is neglected
    (second order: diffusion is a slow spreading on top of the fast
    secular e-oscillation). Only q is diffused (a and I held at the
    donor's values); this is the minimal choice targeting the observable
    that fails (q-shape).
  * u (replenishment): when a particle is removed by the loss process,
    its slot is refilled with a fresh particle with probability f
    (f = replenishment fraction, the u parameter; f = 0 is the (a2)
    no-replenishment case, f = 1 is full replenishment). The fresh
    particle is drawn from the same scattered-disc population: its
    trajectory is a donor trajectory from the same group's library
    (bootstrap resample -- the donor's secular phase at the injection
    time is a random draw from the ensemble, valid under the
    stationarity of the secular ensemble). Its q starts at the donor's
    q at the injection step and then follows drift + diffusion. The
    injection rate is therefore f x (loss rate) -- u is slaved to the
    loss rate in this implementation, a stated modelling choice.
  * Loss: same as (a2): Poisson hazard 1/tau while q < 30, tau = 100 Myr
    (Batygin's clearing timescale; physical). tau = infinity = no loss.
  * Common random numbers: one RNG seed per group, same draw sequence
    for every (Dq, config) combo, so differences between combos are not
    Monte-Carlo noise.

Test matrix per group (version x S):
  Dq in [0, 0.01, 0.03, 0.1, 0.3] AU^2/Myr
  config in {no-loss, loss tau=100 f=0, loss tau=100 f=0.5, loss tau=100 f=1}
  => D-only (no loss), D+loss (no u), u-only+loss (D=0), D+u+loss (main).
  f is inert when there is no loss (nothing to replace).

Metrics (same as before, plus replenishment diagnostics):
  * ever-crossed: initial particles only (comparable to baseline
    25.8-27.5%) and all distinct particles (initial + injected).
  * time-averaged instantaneous fraction: full 4.5 Gyr, last 1 Gyr,
    and a in [100,1000] subset (Batygin aperture).
  * snapshot q: median, q<20, q<25, min; KS of per-episode qmin vs
    flat U[16,30] (episodes found within each particle-lifetime segment,
    so an injection boundary never merges two episodes).
  * n_injected, n_removed, steady-state check (last-1-Gyr fraction).
Snapshots (a, q, i) are saved for the A_100_400_S_9 group only (the main
xi-test group) to crosser_diffusion_snapshots.npz; xi_test.py is run
afterwards on the chosen combos via a small JSON built from the .npz.
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
DT_MYR = DT / 1e6
N = 1000
SEED = 123
DQS = [0.0, 0.01, 0.03, 0.1, 0.3]     # AU^2/Myr
TAU_MYR = 100.0                       # physical loss timescale
FS = [0.0, 0.5, 1.0]                  # replenishment fraction
LAST1_N = 200                         # last 1 Gyr = last 200 sample points
SEED_BASE = 9000


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


def run_combo(Q, IE, WP, a_s, Dq, tau_myr, f_repl, seed, record_snap):
    """One (Dq, tau, f) combo: time-stepped drift + diffusion + loss + u.

    Returns dict of summary metrics; if record_snap, also snapshot arrays.
    Common random numbers: the caller passes the same seed for every combo
    of a group, and the draw order per step is fixed (loss U, repl U,
    new donor, diffusion kick), so differences between combos are not
    Monte-Carlo noise.
    """
    rng = np.random.default_rng(seed)
    sig_q = np.sqrt(2.0 * Dq * DT_MYR) if Dq > 0 else 0.0
    p_loss = (1.0 - np.exp(-DT / (tau_myr * 1e6))
              if tau_myr is not None else 0.0)
    donor = np.arange(N, dtype=np.int64)   # current donor traj per slot
    q = Q[:, 0].copy()                     # current q per slot
    alive = np.ones(N, dtype=bool)
    seg_start = np.zeros(N, dtype=int)     # start step of current particle
    seg_bounds = [[] for _ in range(N)]    # (start, end) per particle-lifetime
    Qslot = np.empty((N, NPTS))
    Aslot = np.empty((N, NPTS))
    Islot = np.empty((N, NPTS))
    Dslot = np.empty((N, NPTS), dtype=np.int64)
    n_injected = 0
    n_removed = 0
    snap_a, snap_q, snap_i = [], [], []
    for j in range(NPTS):
        # ---- loss at step j (on current q), then replenishment ----
        if tau_myr is not None:
            cross = alive & (q < 30.0)
            U_loss = rng.uniform(size=N)
            hit = cross & (U_loss < p_loss)
            if hit.any():
                n_removed += int(hit.sum())
                for k in np.where(hit)[0]:
                    seg_bounds[k].append((int(seg_start[k]), j))
                U_repl = rng.uniform(size=N)
                do_repl = hit & (U_repl < f_repl)
                new_donor = rng.integers(0, N, size=N)
                donor[do_repl] = new_donor[do_repl]
                q[do_repl] = Q[new_donor[do_repl], j]
                seg_start[do_repl] = j
                n_injected += int(do_repl.sum())
                alive[hit & ~do_repl] = False
            else:
                # keep RNG sequence aligned across combos: still draw
                rng.uniform(size=N); rng.integers(0, N, size=N)
        else:
            rng.uniform(size=N); rng.uniform(size=N); rng.integers(0, N, size=N)
        # ---- record state at step j (post-loss, post-replacement) ----
        cross_alive = alive & (q < 30.0)
        Qslot[:, j] = q
        Aslot[:, j] = a_s[donor]
        Islot[:, j] = IE[donor, j]
        Dslot[:, j] = donor
        if record_snap and cross_alive.any():
            idx = np.where(cross_alive)[0]
            snap_a.extend(a_s[donor[idx]].tolist())
            snap_q.extend(q[idx].tolist())
            snap_i.extend(IE[donor[idx], j].tolist())
        # ---- advance to j+1: secular drift + diffusion kick ----
        if j + 1 < NPTS:
            dq_sec = Q[donor, j + 1] - Q[donor, j]
            q = q + dq_sec
            kick = rng.normal(0.0, sig_q, size=N)
            q = q + kick
            q = np.abs(q)   # reflecting boundary at q = 0
    # finalize segments for surviving particles
    for k in range(N):
        if alive[k]:
            seg_bounds[k].append((int(seg_start[k]), NPTS))
    # ---- metrics ----
    # alive status per slot per step, reconstructed from segments
    alive_mask = np.zeros((N, NPTS), dtype=bool)
    for k in range(N):
        for (s, e) in seg_bounds[k]:
            alive_mask[k, s:e] = True
    cmask = alive_mask & (Qslot < 30.0)
    inst_full = float(cmask.sum() / alive_mask.sum())
    inst_last1 = float(cmask[:, -LAST1_N:].sum() / alive_mask[:, -LAST1_N:].sum())
    sel = a_s <= 1000.0
    # a per slot per step for the a<=1000 cut: use Aslot
    sel_mask = Aslot <= 1000.0
    inst_full_a1000 = float((cmask & sel_mask).sum() / (alive_mask & sel_mask).sum())
    inst_last1_a1000 = float((cmask[:, -LAST1_N:] & sel_mask[:, -LAST1_N:]).sum()
                             / (alive_mask[:, -LAST1_N:] & sel_mask[:, -LAST1_N:]).sum())
    # ever-crossed: initial particles (first segment) and all particles
    ever_init = np.zeros(N, dtype=bool)
    n_ever_all = 0
    n_particles = N + n_injected
    for k in range(N):
        segs = seg_bounds[k]
        if segs:
            s0, e0 = segs[0]
            ever_init[k] = bool(cmask[k, s0:e0].any())
        for (s, e) in segs:
            if cmask[k, s:e].any():
                n_ever_all += 1
    # episodes within segments (an injection boundary never merges two)
    episodes = []
    for k in range(N):
        for (s, e) in seg_bounds[k]:
            for (i0, i1) in find_episodes(cmask[k, s:e]):
                kmin = s + i0 + int(np.argmin(Qslot[k, s + i0:s + i1 + 1]))
                dnr = int(Dslot[k, kmin])
                episodes.append(dict(a=float(Aslot[k, kmin]),
                                     dur_Myr=float((i1 - i0 + 1) * DT_MYR),
                                     qmin=float(Qslot[k, kmin]),
                                     wp_qmin=float(WP[dnr, kmin]),
                                     i_qmin=float(Islot[k, kmin]),
                                     censored=bool(i1 == e - s - 1 and e < NPTS)))
    sq = np.array(snap_q) if snap_q else np.array([])
    qm_ep = np.array([e["qmin"] for e in episodes]) if episodes else np.array([])
    ks_flat = (stats.kstest(qm_ep, "uniform", args=(16, 14))
               if len(qm_ep) else None)
    return dict(
        Dq=float(Dq), tau_Myr=tau_myr, f_repl=float(f_repl),
        n=N, n_injected=int(n_injected), n_removed=int(n_removed),
        n_particles=int(n_particles),
        frac_ever_init=float(ever_init.mean()),
        frac_ever_all=float(n_ever_all / max(n_particles, 1)),
        inst_full=inst_full, inst_last1=inst_last1,
        inst_full_a1000=inst_full_a1000, inst_last1_a1000=inst_last1_a1000,
        snap_q_median=float(np.median(sq)) if len(sq) else None,
        snap_q_lt20=float((sq < 20).mean()) if len(sq) else None,
        snap_q_lt25=float((sq < 25).mean()) if len(sq) else None,
        snap_q_min=float(sq.min()) if len(sq) else None,
        n_snap=int(len(sq)), n_episodes=int(len(episodes)),
        ks_flat_D=(float(ks_flat.statistic) if ks_flat else None),
        ks_flat_p=(float(ks_flat.pvalue) if ks_flat else None),
        snap_a=snap_a, snap_q=snap_q, snap_i=snap_i)


# ---------------- load calibration ----------------
with open("crosser_simulation_results.json") as f:
    prev = json.load(f)
A300_MIN = prev["A300_MIN"]
A300_9 = prev["A300_9"]
print(f"A300_MIN = {A300_MIN:.4e}, A300_9 = {A300_9:.4e}; DT = {DT_MYR:.0f} Myr; "
      f"Dq grid {DQS} AU^2/Myr; tau = {TAU_MYR:.0f} Myr; f grid {FS}", flush=True)

D_CAR = axis_from_lonlat(30.0, 75.0)
XPR, YPR, E = frame_basis(D_CAR)

results = {}
snap_store = {}
configs = ([("no_loss", None, 0.0)]
           + [(f"tau100_f{f}", TAU_MYR, f) for f in FS])
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
        print(f"  {version}/{tag}: trajectories done ({time.time()-t0:.0f}s)",
              flush=True)
        seed_g = SEED_BASE + v_idx * 10 + s_idx
        for Dq in DQS:
            for c_idx, (clabel, tau_myr, f_repl) in enumerate(configs):
                key = f"{version}_{tag}_D{Dq:g}_{clabel}"
                record_snap = (version == "A_100_400")
                out = run_combo(Q, IE, WP, a_s, Dq, tau_myr, f_repl,
                                seed_g, record_snap)
                snap_a, snap_q, snap_i = out.pop("snap_a"), out.pop("snap_q"), out.pop("snap_i")
                out["version"] = version; out["tag"] = tag
                results[key] = out
                if record_snap:
                    snap_store[key + "_a"] = np.array(snap_a)
                    snap_store[key + "_q"] = np.array(snap_q)
                    snap_store[key + "_i"] = np.array(snap_i)
                print(f"  {key}: EVER_init {out['frac_ever_init']*100:.1f}% | "
                      f"EVER_all {out['frac_ever_all']*100:.1f}% | "
                      f"inst_full {out['inst_full']*100:.2f}% | "
                      f"inst_last1Gyr {out['inst_last1']*100:.2f}% | "
                      f"a<=1000 {out['inst_full_a1000']*100:.2f}% | "
                      f"snap q med {out['snap_q_median'] if out['snap_q_median'] is not None else float('nan'):.1f} AU, "
                      f"q<20 {out['snap_q_lt20']*100 if out['snap_q_lt20'] is not None else float('nan'):.1f}%, "
                      f"q<25 {out['snap_q_lt25']*100 if out['snap_q_lt25'] is not None else float('nan'):.1f}% | "
                      f"ep n={out['n_episodes']}"
                      + (f", KS D={out['ks_flat_D']:.3f} p={out['ks_flat_p']:.2e}" if out['ks_flat_D'] is not None else "")
                      + f" | injected {out['n_injected']}, removed {out['n_removed']}",
                      flush=True)
        with open("crosser_diffusion_results.json", "w") as f:
            json.dump(results, f)
        np.savez_compressed("crosser_diffusion_snapshots.npz", **snap_store)
        print("  saved crosser_diffusion_results.json + "
              "crosser_diffusion_snapshots.npz (incremental)", flush=True)

print("\nDONE")
