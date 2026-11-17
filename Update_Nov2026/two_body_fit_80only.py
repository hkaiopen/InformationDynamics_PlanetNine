#!/usr/bin/env python3
"""TRUE MODEL BLIND TEST: fit ONLY the 80-200 AU Table-1 objects, predict the
200-400 AU bin WITHOUT using its orientations. The 200-400 orientations are
never in the likelihood; only their semi-major axes are used to evaluate the
prediction. If the 80-200-only fit predicts the 200-400 second lobe, the BaZi
(fitted on one bin) predicts another bin -- a model blind test, not a
sample-consistency test. Same machinery as two_body_fit.py (Siraj conditional
von Mises likelihood, soft prior tilt(65) ~ N(0,3deg), alpha cap 0.95).
"""
import numpy as np
from scipy.special import i0
from scipy.optimize import differential_evolution, minimize_scalar

z = np.load("table1_vectors.npz")
R_all, J_all, A_all = z["R"], z["J"], z["A"]
mask80 = (A_all >= 80) & (A_all <= 200)
mask200 = (A_all > 200) & (A_all <= 400)
R, J, A = R_all[mask80], J_all[mask80], A_all[mask80]
R200, J200, A200 = R_all[mask200], J_all[mask200], A_all[mask200]
print(f"fit data (80-200 only): n={len(A)} objects, a {A.min():.1f}-{A.max():.1f} AU", flush=True)
print(f"blind-test bin (200-400, orientations NOT used): n={len(A200)} objects", flush=True)

exec(open("bodies_coupling_fit.py").read().split('print("\\n--- giants only')[0])
ALPHA_CAP = 0.95
G_M = np.array([m for _, m, _, _, _ in GIANTS])
G_A = np.array([a for _, _, a, _, _ in GIANTS])
G_PQ = np.array([incl_vec(i, Om) for _, _, _, i, Om in GIANTS])
INV = np.array([np.cos(np.radians(88.42)) * np.cos(np.radians(17.58)),
                np.cos(np.radians(88.42)) * np.sin(np.radians(17.58)),
                np.sin(np.radians(88.42))])


def b_of(alpha):
    return np.interp(np.clip(alpha, 0.01, ALPHA_CAP), agrid, bgrid)


def forced_pq_vec(a_arr, extra):
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


def loglik(th, lg, RR, JJ, AA):
    gamma = 10 ** lg
    p, q = forced_pq_vec(AA, extra_from_theta(th))
    M = pole_from_pq_vec(p, q)
    mdotR = np.einsum("ij,ij->i", M, RR)
    mdotJ = np.einsum("ij,ij->i", M, JJ)
    sinth = np.sqrt(np.clip(1 - mdotR ** 2, 0, 1))
    return float(np.sum(gamma * mdotJ - np.log(i0(gamma * sinth))))


def logprior_soft(th):
    p, q = forced_pq_vec(np.array([65.0]), extra_from_theta(th))
    t65 = tilt_from_pq(p[0], q[0])
    return -0.5 * (t65 / 3.0) ** 2, t65


def negpost_two(x):
    return -(loglik(x[:8], x[8], R, J, A) + logprior_soft(x[:8])[0])


def negpost_yonly(x):
    th = np.array([x[0], x[1], x[2], x[3], -6.0, 500.0, 0.0, 0.0])
    return -(loglik(th, x[4], R, J, A) + logprior_soft(th)[0])


BOUNDS_TWO = [(np.log10(0.002), np.log10(4)), (80, 300), (5, 55), (0, 360),
              (np.log10(0.5), np.log10(10)), (200, 800), (0, 30), (0, 360),
              (np.log10(3), np.log10(100))]
BOUNDS_Y = [(np.log10(0.002), np.log10(4)), (80, 300), (5, 55), (0, 360),
            (np.log10(3), np.log10(100))]

print("\n--- DE: two-body MLE (80-200 only) ---", flush=True)
res2 = differential_evolution(negpost_two, BOUNDS_TWO, seed=7, maxiter=200,
                              popsize=20, tol=1e-7, polish=True)
th2, lg2 = res2.x[:8], res2.x[8]
print(f"negpost={res2.fun:.2f}  logL(80-200)={loglik(th2, lg2, R, J, A):.2f}  "
      f"tilt65={logprior_soft(th2)[1]:.2f} deg  gamma={10**lg2:.1f}", flush=True)
print(f"Planet Y: m={10**th2[0]:.3f} M_E, a={th2[1]:.0f} AU, i={th2[2]:.1f} deg, Om={th2[3]:.0f} deg", flush=True)
print(f"Planet 9: m={10**th2[4]:.2f} M_E, a={th2[5]:.0f} AU, i={th2[6]:.1f} deg, Om={th2[7]:.0f} deg", flush=True)

print("\n--- DE: Y-only MLE (80-200 only) ---", flush=True)
resY = differential_evolution(negpost_yonly, BOUNDS_Y, seed=7, maxiter=200,
                              popsize=20, tol=1e-7, polish=True)
thY = np.array([resY.x[0], resY.x[1], resY.x[2], resY.x[3], -6.0, 500.0, 0.0, 0.0])
print(f"negpost={resY.fun:.2f}  logL(80-200)={loglik(thY, resY.x[4], R, J, A):.2f}  "
      f"tilt65={logprior_soft(thY)[1]:.2f} deg  gamma={10**resY.x[4]:.1f}", flush=True)
print(f"Planet Y: m={10**thY[0]:.3f} M_E, a={thY[1]:.0f} AU, i={thY[2]:.1f} deg, Om={thY[3]:.0f} deg", flush=True)
print(f"Delta logL (two-body minus Y-only, 80-200 only) = "
      f"{loglik(th2, lg2, R, J, A) - loglik(thY, resY.x[4], R, J, A):.2f} (+3 params)", flush=True)

rn = minimize_scalar(lambda lg: -(10 ** lg * (INV @ J.sum(axis=0))
                     - np.sum(np.log(i0(10 ** lg * np.sqrt(np.clip(1 - (R @ INV) ** 2, 0, 1)))))),
                     bounds=(np.log10(3), np.log10(100)), method="bounded")
print(f"null (invariable plane, 80-200 only): logL={-rn.fun:.2f}, gamma={10**rn.x:.1f}", flush=True)
print(f"Delta logL (two-body minus null, 80-200 only) = {loglik(th2, lg2, R, J, A) + rn.fun:.2f}", flush=True)

# --- MCMC (two-body), 4 chains ---
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


print("\n--- MCMC (4 chains x 20000, 80-200 only) ---", flush=True)
chains = []
for c in range(4):
    ch, acc = run_chain(100 + c)
    chains.append(ch[3000:])
    print(f"chain {c}: acceptance {acc:.2f}", flush=True)
S = np.vstack(chains)
np.savez("two_body_80only_mcmc.npz", S=S, th2=th2, lg2=lg2, thY=thY, lgY=resY.x[4],
         logL2=loglik(th2, lg2, R, J, A), logLY=loglik(thY, resY.x[4], R, J, A))

# --- BLIND PREDICTION: 200-400 AU tilt/node from the 80-200-only fit ---
print("\n--- BLIND PREDICTION for 200-400 AU (fit never saw these orientations) ---", flush=True)
print("measured (Table-1 subset): tilt 13.9 deg, node 125 deg, 1.34 sigma", flush=True)
print("measured (current full sample n=29): tilt 11.7 deg, node 130 deg, 1.99 sigma", flush=True)
for tag, th in [("two-body MLE", th2), ("Y-only MLE", thY)]:
    p, q = forced_pq_vec(A200, extra_from_theta(th))
    ts = tilt_from_pq(p, q)
    node = float(np.degrees(np.arctan2(p.mean(), q.mean())) % 360)
    print(f"{tag}: predicted 200-400 mean tilt = {np.mean(ts):.1f} deg "
          f"(per-object {np.min(ts):.1f}-{np.max(ts):.1f}), node = {node:.0f} deg", flush=True)
    for a in [250, 300, 350]:
        pp, qq = forced_pq_vec(np.array([float(a)]), extra_from_theta(th))
        nn = float(np.degrees(np.arctan2(pp[0], qq[0])) % 360)
        print(f"   a={a}: tilt={tilt_from_pq(pp[0], qq[0]):.1f} deg node={nn:.0f}", flush=True)
# posterior predictive for the 200-400 bin mean tilt (two-body)
t200 = []
for s in S[::40]:
    p, q = forced_pq_vec(A200, extra_from_theta(s[:8]))
    t200.append(np.mean(tilt_from_pq(p, q)))
lo, med, hi = np.percentile(t200, [16, 50, 84])
print(f"two-body posterior predictive 200-400 mean tilt: {med:.1f} [{lo:.1f}, {hi:.1f}] deg", flush=True)
print("saved two_body_80only_mcmc.npz")
print("DONE")
