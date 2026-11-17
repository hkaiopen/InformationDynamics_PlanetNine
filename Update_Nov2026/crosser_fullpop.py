#!/usr/bin/env python3
"""Full-population D+u simulation: detached-object (sednoid) production rate.

Task (b) (user, 2026-11-17): the data audit (data_audit.md §6) found
diffusion opens a non-Kozai q-lifting route -- D_q = 0.3 AU^2/Myr spreads
q by ~52 AU (1 sigma) over 4.5 Gyr, and the four detached objects need
Delta q = 34.7-50.6 AU = 0.67-0.97 sigma, with the SAME D used to flatten
the crossers (no extra tuning). But the RATE was unverified:
crosser_diffusion.py only records crosser footprints (q<30), so its
snapshots contain 0% with q>60 by construction, not by physics.
This script records the FULL population and measures how many detached
objects (q>60, a>240) the D+u model produces, against the 4 observed
(Sedna, 2012 VP113, Leleakuhonua, Ammonite).

Machinery copied verbatim from crosser_diffusion.py (secular trajectories,
population, drift + diffusion + loss + replenishment, common random
numbers). Only the recording and the metrics are new: every alive slot's
(a, q, i) is used (Qslot/Aslot/Islot already hold the full population),
and detached metrics are computed per particle-lifetime segment.

Detached definition (primary): q > 60 AU and a > 240 AU -- captures all
four observed (Ammonite a=247.5 is the binding case; a>250 would miss it).
Sensitivity: q>60 & a>250, q>60 & a>200, q>50 & a>250.
Observed comparison (MPCORB-202611.DAT, counted 2026-11-17):
a>240 & q>30: 54 known objects (30 with U<=6); q>60 among them: 4
(7.4%; 4/30 = 13% among U<=6). a>200 & q>30: 67 known, 4 detached (6.0%).
Discovery bias works AGAINST high q, so 4 is a lower bound on the true
detached fraction -- the model must produce AT LEAST ~6-7% (of the
a>240, q>30 population) to be consistent.
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
DQS = [0.0, 0.1, 0.3]     # AU^2/Myr (main D=0.3, controls D=0.1, D=0)
TAU_MYR = 100.0                       # physical loss timescale
FS = [1.0]                  # replenishment fraction (main: full, f=1.0)
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
    # ---- detached metrics (full population; NEW for this task) ----
    # Primary detached definition: q > 60 AU and a > 240 AU (captures all
    # four observed; Ammonite a=247.5 is the binding case).
    DET_Q = 60.0
    DET_A = 240.0
    det_mask = alive_mask & (Qslot > DET_Q) & (Aslot > DET_A)
    # denominator matching the observed census: a > 240 & q > 30, alive
    denom_mask = alive_mask & (Aslot > DET_A) & (Qslot > 30.0)
    # per particle-lifetime segment: born detached / ever detached /
    # became detached (born q<=60, later reached q>60, a>240)
    n_seg_a240 = 0
    n_born_det = 0
    n_ever_det = 0
    n_became_det = 0
    n_born_det_init = 0
    n_ever_det_init = 0
    n_became_det_init = 0
    first_pass = []  # (first-passage time Myr, a, q at passage)
    for k in range(N):
        for si, (s, e) in enumerate(seg_bounds[k]):
            a_seg = float(Aslot[k, s])
            if a_seg <= DET_A:
                continue
            n_seg_a240 += 1
            q_birth = float(Qslot[k, s])
            born_det = q_birth > DET_Q
            ever_det = bool(det_mask[k, s:e].any())
            if born_det:
                n_born_det += 1
            if ever_det:
                n_ever_det += 1
                if not born_det:
                    n_became_det += 1
                    j0 = s + int(np.where(det_mask[k, s:e])[0][0])
                    first_pass.append((float(j0 * DT_MYR), a_seg,
                                       float(Qslot[k, j0])))
            if si == 0:  # initial particle of this slot
                if born_det:
                    n_born_det_init += 1
                if ever_det:
                    n_ever_det_init += 1
                    if not born_det:
                        n_became_det_init += 1
    # instantaneous detached fraction (denominator: a>240 & q>30, alive)
    n_denom_final = int(denom_mask[:, -1].sum())
    n_det_final = int(det_mask[:, -1].sum())
    frac_det_final = float(n_det_final / max(n_denom_final, 1))
    n_denom_last1 = int(denom_mask[:, -LAST1_N:].sum())
    n_det_last1 = int(det_mask[:, -LAST1_N:].sum())
    frac_det_last1 = float(n_det_last1 / max(n_denom_last1, 1))
    n_denom_t0 = int(denom_mask[:, 0].sum())
    n_det_t0 = int(det_mask[:, 0].sum())
    frac_det_t0 = float(n_det_t0 / max(n_denom_t0, 1))
    # per-step detached fraction time series (same denominator)
    frac_det_t = (det_mask.sum(axis=0)
                  / np.maximum(denom_mask.sum(axis=0), 1))
    # detached q distribution at final step
    det_q_final = Qslot[det_mask[:, -1], -1] if n_det_final else np.array([])
    det_a_final = Aslot[det_mask[:, -1], -1] if n_det_final else np.array([])
    # sensitivity thresholds (same denominator rule, matched a-cut)
    def frac_det_qa(qthr, athr):
        m = alive_mask & (Qslot > qthr) & (Aslot > athr)
        d = alive_mask & (Aslot > athr) & (Qslot > 30.0)
        return float(m[:, -1].sum() / max(d[:, -1].sum(), 1))
    sens = {
        'q60_a250': frac_det_qa(60.0, 250.0),
        'q60_a200': frac_det_qa(60.0, 200.0),
        'q50_a250': frac_det_qa(50.0, 250.0),
    }
    # full-population final-step arrays (alive slots) for later analysis
    alive_final = alive_mask[:, -1]
    full_a = Aslot[alive_final, -1]
    full_q = Qslot[alive_final, -1]
    full_i = Islot[alive_final, -1]
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
        snap_a=snap_a, snap_q=snap_q, snap_i=snap_i,
        # detached (NEW)
        det_n_seg_a240=int(n_seg_a240),
        det_frac_born=float(n_born_det / max(n_seg_a240, 1)),
        det_frac_ever=float(n_ever_det / max(n_seg_a240, 1)),
        det_frac_became=float(n_became_det / max(n_seg_a240, 1)),
        det_frac_born_init=float(n_born_det_init / max(N, 1)),
        det_frac_ever_init=float(n_ever_det_init / max(N, 1)),
        det_frac_became_init=float(n_became_det_init / max(N, 1)),
        det_frac_final=frac_det_final,
        det_frac_last1=frac_det_last1,
        det_frac_t0=frac_det_t0,
        det_n_final=n_det_final,
        det_n_denom_final=n_denom_final,
        det_n_last1=n_det_last1,
        det_n_denom_last1=n_denom_last1,
        det_q_final_median=(float(np.median(det_q_final))
                            if len(det_q_final) else None),
        det_q_final_min=(float(det_q_final.min())
                         if len(det_q_final) else None),
        det_q_final_max=(float(det_q_final.max())
                         if len(det_q_final) else None),
        det_a_final_median=(float(np.median(det_a_final))
                            if len(det_a_final) else None),
        det_first_pass_median=(float(np.median([f[0] for f in first_pass]))
                               if first_pass else None),
        det_first_pass_n=int(len(first_pass)),
        det_sens=sens,
        det_frac_t=frac_det_t.tolist(),
        full_a=full_a.tolist(), full_q=full_q.tolist(),
        full_i=full_i.tolist())


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
                full_a = np.array(out.pop("full_a"))
                full_q = np.array(out.pop("full_q"))
                full_i = np.array(out.pop("full_i"))
                det_frac_t = np.array(out.pop("det_frac_t"))
                out["version"] = version; out["tag"] = tag
                results[key] = out
                if record_snap:
                    snap_store[key + "_a"] = np.array(snap_a)
                    snap_store[key + "_q"] = np.array(snap_q)
                    snap_store[key + "_i"] = np.array(snap_i)
                # full-population final-step arrays + detached time series
                snap_store[key + "_full_a"] = full_a
                snap_store[key + "_full_q"] = full_q
                snap_store[key + "_full_i"] = full_i
                snap_store[key + "_det_frac_t"] = det_frac_t
                print(f"  {key}: DET born {out['det_frac_born']*100:.1f}% | "
                      f"EVER {out['det_frac_ever']*100:.1f}% | "
                      f"BECAME {out['det_frac_became']*100:.1f}% | "
                      f"final {out['det_frac_final']*100:.1f}% "
                      f"({out['det_n_final']}/{out['det_n_denom_final']}) | "
                      f"last1Gyr {out['det_frac_last1']*100:.1f}% | "
                      f"t0 {out['det_frac_t0']*100:.1f}% | "
                      f"det q med {out['det_q_final_median'] if out['det_q_final_median'] is not None else float('nan'):.1f} | "
                      f"1st-pass med {out['det_first_pass_median'] if out['det_first_pass_median'] is not None else float('nan'):.0f} Myr | "
                      f"sens q60_a250 {out['det_sens']['q60_a250']*100:.1f}% | "
                      f"cross inst_last1 {out['inst_last1']*100:.2f}%",
                      flush=True)
        with open("crosser_fullpop_results.json", "w") as f:
            json.dump(results, f)
        np.savez_compressed("crosser_fullpop_snapshots.npz", **snap_store)
        print("  saved crosser_fullpop_results.json + "
              "crosser_fullpop_snapshots.npz (incremental)", flush=True)

print("\nDONE")
