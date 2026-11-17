#!/usr/bin/env python3
"""Step 1 of the analytic route: nail down the two-body picture.
Joint Bayesian fit of the forced-plane model to the 46 Siraj Table-1 objects,
Siraj conditional von Mises likelihood, mean plane predicted per object at its
own a (no binning). 50-80 AU no-warp constraint as a soft prior (tilt at 65 AU
~ N(0, 3 deg)). Vectorised likelihood (all 46 objects at once). No GL first
principles in this step, by design. Loads table1_vectors.npz (run
make_table1_cache.py first).
"""
import numpy as np
from scipy.special import i0
from scipy.optimize import differential_evolution, minimize_scalar

z = np.load("table1_vectors.npz")
R, J, A = z["R"], z["J"], z["A"]
print(f"fit data: n={len(A)} objects, a {A.min():.1f}-{A.max():.1f} AU", flush=True)

exec(open("bodies_coupling_fit.py").read().split('print("\\n--- giants only')[0])
ALPHA_CAP = 0.95
# giants as arrays
G_M = np.array([m for _, m, _, _, _ in GIANTS])
G_A = np.array([a for _, _, a, _, _ in GIANTS])
G_PQ = np.array([incl_vec(i, Om) for _, _, _, i, Om in GIANTS])  # (4,2)
INV = np.array([np.cos(np.radians(88.42)) * np.cos(np.radians(17.58)),
                np.cos(np.radians(88.42)) * np.sin(np.radians(17.58)),
                np.sin(np.radians(88.42))])


def b_of(alpha):
    return np.interp(np.clip(alpha, 0.01, ALPHA_CAP), agrid, bgrid)


def forced_pq_vec(a_arr, extra):
    """extra: list of (m, a_j, p_j, q_j). Returns (p_f, q_f) arrays."""
    a_arr = np.asarray(a_arr, dtype=float)
    num = np.zeros((len(a_arr), 2))
    den = np.zeros(len(a_arr))
    bodies = [(G_M[k], G_A[k], G_PQ[k, 0], G_PQ[k, 1]) for k in range(4)] + extra
    for m, a_j, pj, qj in bodies:
        interior = a_j < a_arr
        alpha = np.where(interior, a_j / a_arr, a_arr / a_j)
        w = np.where(interior, m * alpha * b_of(alpha),
                     m * alpha ** 2 * b_of(alpha))
        num[:, 0] += w * pj
        num[:, 1] += w * qj
        den += w
    return num[:, 0] / den, num[:, 1] / den


def pole_from_pq_vec(p, q):
    s = np.hypot(p, q)
    i0a = 2 * np.arcsin(np.clip(s, 0, 1))
    Om0 = np.arctan2(p, q)
    return np.column_stack([np.sin(i0a) * np.sin(Om0),
                            -np.sin(i0a) * np.cos(Om0), np.cos(i0a)])


def tilt_from_pq(p, q):
    return 2 * np.degrees(np.arcsin(np.clip(
        np.hypot(p - INV_P, q - INV_Q), 0, 1)))


def extra_from_theta(th):
    lmY, aY, iY, OY, lm9, a9, i9, O9 = th
    vY = incl_vec(iY, OY)
    v9 = incl_vec(i9, O9)
    return [(10 ** lmY, aY, vY[0], vY[1]), (10 ** lm9, a9, v9[0], v9[1])]


def loglik(th, lg):
    gamma = 10 ** lg
    p, q = forced_pq_vec(A, extra_from_theta(th))
    M = pole_from_pq_vec(p, q)
    mdotR = np.einsum("ij,ij->i", M, R)
    mdotJ = np.einsum("ij,ij->i", M, J)
    sinth = np.sqrt(np.clip(1 - mdotR ** 2, 0, 1))
    return float(np.sum(gamma * mdotJ - np.log(i0(gamma * sinth))))


def logprior_soft(th):
    p, q = forced_pq_vec(np.array([65.0]), extra_from_theta(th))
    t65 = tilt_from_pq(p[0], q[0])
    return -0.5 * (t65 / 3.0) ** 2, t65


def negpost_two(x):
    return -(loglik(x[:8], x[8]) + logprior_soft(x[:8])[0])


def negpost_yonly(x):
    th = np.array([x[0], x[1], x[2], x[3], -6.0, 500.0, 0.0, 0.0])
    return -(loglik(th, x[4]) + logprior_soft(th)[0])


BOUNDS_TWO = [(np.log10(0.002), np.log10(4)), (80, 300), (5, 55), (0, 360),
              (np.log10(0.5), np.log10(10)), (200, 800), (0, 30), (0, 360),
              (np.log10(3), np.log10(100))]
BOUNDS_Y = [(np.log10(0.002), np.log10(4)), (80, 300), (5, 55), (0, 360),
            (np.log10(3), np.log10(100))]

print("\n--- DE: two-body MLE ---", flush=True)
res2 = differential_evolution(negpost_two, BOUNDS_TWO, seed=7, maxiter=200,
                              popsize=20, tol=1e-7, polish=True)
th2, lg2 = res2.x[:8], res2.x[8]
print(f"negpost={res2.fun:.2f}  logL={loglik(th2, lg2):.2f}  "
      f"tilt65={logprior_soft(th2)[1]:.2f} deg  gamma={10**lg2:.1f}", flush=True)
print(f"Planet Y: m={10**th2[0]:.3f} M_E, a={th2[1]:.0f} AU, i={th2[2]:.1f} deg, Om={th2[3]:.0f} deg", flush=True)
print(f"Planet 9: m={10**th2[4]:.2f} M_E, a={th2[5]:.0f} AU, i={th2[6]:.1f} deg, Om={th2[7]:.0f} deg", flush=True)

print("\n--- DE: Y-only MLE ---", flush=True)
resY = differential_evolution(negpost_yonly, BOUNDS_Y, seed=7, maxiter=200,
                              popsize=20, tol=1e-7, polish=True)
thY = np.array([resY.x[0], resY.x[1], resY.x[2], resY.x[3], -6.0, 500.0, 0.0, 0.0])
print(f"negpost={resY.fun:.2f}  logL={loglik(thY, resY.x[4]):.2f}  "
      f"tilt65={logprior_soft(thY)[1]:.2f} deg  gamma={10**resY.x[4]:.1f}", flush=True)
print(f"Planet Y: m={10**thY[0]:.3f} M_E, a={thY[1]:.0f} AU, i={thY[2]:.1f} deg, Om={thY[3]:.0f} deg", flush=True)
print(f"Delta logL (two-body minus Y-only) = {loglik(th2, lg2) - loglik(thY, resY.x[4]):.2f} (+3 params)", flush=True)

rn = minimize_scalar(lambda lg: -(10 ** lg * (INV @ J.sum(axis=0))
                     - np.sum(np.log(i0(10 ** lg * np.sqrt(np.clip(1 - (R @ INV) ** 2, 0, 1)))))),
                     bounds=(np.log10(3), np.log10(100)), method="bounded")
print(f"null (invariable plane): logL={-rn.fun:.2f}, gamma={10**rn.x:.1f}", flush=True)
print(f"Delta logL (two-body minus null) = {loglik(th2, lg2) + rn.fun:.2f}", flush=True)

# --- Metropolis MCMC (two-body), 4 chains ---
rng = np.random.default_rng(3)
LO = np.array([b[0] for b in BOUNDS_TWO])
HI = np.array([b[1] for b in BOUNDS_TWO])
STEP = np.array([0.12, 12.0, 4.0, 15.0, 0.12, 50.0, 2.5, 15.0, 0.08])


def logpost(x):
    if np.any(x < LO) or np.any(x > HI):
        return -np.inf
    return -negpost_two(x)


def run_chain(seed, nsteps=20000):
    r = np.random.default_rng(seed)
    x = res2.x + r.normal(0, STEP * 0.3)
    x = np.clip(x, LO, HI)
    cur = logpost(x)
    out = np.empty((nsteps, 9))
    acc = 0
    for k in range(nsteps):
        y = x + r.normal(0, STEP)
        y[3] %= 360.0
        y[7] %= 360.0
        lp = logpost(y)
        if np.log(r.uniform()) < lp - cur:
            x, cur = y, lp
            acc += 1
        out[k] = x
    return out, acc / nsteps


print("\n--- MCMC (4 chains x 20000) ---", flush=True)
chains = []
for c in range(4):
    ch, acc = run_chain(100 + c)
    chains.append(ch[3000:])
    print(f"chain {c}: acceptance {acc:.2f}", flush=True)
S = np.vstack(chains)
np.savez("two_body_mcmc.npz", S=S, th2=th2, lg2=lg2, thY=thY, lgY=resY.x[4],
         logL2=loglik(th2, lg2), logLY=loglik(thY, resY.x[4]))
names = ["log10 mY", "aY", "iY", "OmY", "log10 m9", "a9", "i9", "Om9", "log10 gamma"]
print("\nposterior (16/50/84%):", flush=True)
for j, nm in enumerate(names):
    lo, med, hi = np.percentile(S[:, j], [16, 50, 84])
    if nm.startswith("log10 m"):
        print(f"  {nm}: {10**med:.3f} M_E  [{10**lo:.3f}, {10**hi:.3f}]", flush=True)
    elif nm.startswith("log10 gamma"):
        print(f"  {nm}: {10**med:.1f}  [{10**lo:.1f}, {10**hi:.1f}]", flush=True)
    else:
        print(f"  {nm}: {med:.1f}  [{lo:.1f}, {hi:.1f}]", flush=True)

print("\n--- predicted mean-plane tilt / node (MLE and posterior) ---", flush=True)
for a in [65, 100, 140, 180, 250, 300, 350]:
    p, q = forced_pq_vec(np.array([float(a)]), extra_from_theta(th2))
    node = float(np.degrees(np.arctan2(p[0], q[0])) % 360)
    ts = []
    for s in S[::40]:
        pp, qq = forced_pq_vec(np.array([float(a)]), extra_from_theta(s[:8]))
        ts.append(tilt_from_pq(pp[0], qq[0]))
    lo, med, hi = np.percentile(ts, [16, 50, 84])
    print(f"a={a:3d} AU: MLE tilt={tilt_from_pq(p[0], q[0]):5.2f} deg node={node:5.0f} | "
          f"posterior tilt {med:.1f} [{lo:.1f}, {hi:.1f}] deg", flush=True)
for lo_a, hi_a, lab, meas in [(80, 200, "80-200", 13.3), (200, 400, "200-400", 13.9)]:
    aa = A[(A >= lo_a) & (A <= hi_a)]
    p, q = forced_pq_vec(aa, extra_from_theta(th2))
    ts = tilt_from_pq(p, q)
    print(f"bin {lab}: MLE predicted mean tilt = {np.mean(ts):.1f} deg (measured {meas} deg)", flush=True)
print("saved two_body_mcmc.npz")
print("DONE")
