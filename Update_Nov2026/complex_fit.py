#!/usr/bin/env python3
"""Complex fit of the mean plane: real part + imaginary part.

Real part: four giants, Laplace secular weights, masses fixed at ephemeris
values (the visible part; explains the invariable plane / baseline).
Imaginary part: the carrier, a single-axis kernel with its OWN weight law
(not the Laplace law) -- a smooth step in a (logistic) as primary model,
per-bin free weights as control. No Planet Y / Planet Nine as real bodies.

Model (per object at semi-major axis a):
  forced (p,q) = [W_G(a)*(p_G,q_G) + w_c(a)*(p_d,q_d)] / [W_G(a) + w_c(a)]
  W_G, (p_G,q_G): giants' Laplace-weighted total weight and mean (p,q)
  w_c(a) = C * f(a): carrier weight; f(a) = logistic step (primary) or
           per-bin step (control); C = relative coupling strength (fitted,
           in Earth-mass units, i.e. relative to the giants' weights)
  (p_d,q_d) = incl_vec(i_d, Om): carrier axis (fitted)
  M(a) = pole_from_pq_vec(p_f, q_f)
Likelihood: Siraj conditional von Mises, shared gamma (same as two_body_fit.py).
Soft prior: tilt at 65 AU ~ N(0, 3 deg) (50-80 AU no-warp constraint).
"""
import numpy as np
from scipy.special import i0 as bessel_i0
from scipy.optimize import differential_evolution

z = np.load("table1_vectors.npz")
R, J, A = z["R"], z["J"], z["A"]
print(f"fit data: n={len(A)} objects, a {A.min():.1f}-{A.max():.1f} AU", flush=True)

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


def giants_pq_w(a_arr):
    """Giants alone: return (p_G, q_G, W_G) arrays."""
    a_arr = np.asarray(a_arr, dtype=float)
    num = np.zeros((len(a_arr), 2))
    den = np.zeros(len(a_arr))
    for k in range(4):
        interior = G_A[k] < a_arr
        alpha = np.where(interior, G_A[k] / a_arr, a_arr / G_A[k])
        wj = np.where(interior, G_M[k] * alpha * b_of(alpha),
                      G_M[k] * alpha ** 2 * b_of(alpha))
        num[:, 0] += wj * G_PQ[k, 0]
        num[:, 1] += wj * G_PQ[k, 1]
        den += wj
    return num[:, 0] / den, num[:, 1] / den, den


def pole_from_pq_vec(p, q):
    s = np.hypot(p, q)
    ia = 2 * np.arcsin(np.clip(s, 0, 1))
    Om0 = np.arctan2(p, q)
    return np.column_stack([np.sin(ia) * np.sin(Om0),
                            -np.sin(ia) * np.cos(Om0), np.cos(ia)])


def tilt_from_pole(M):
    return np.degrees(np.arccos(np.clip(M @ INV, -1, 1)))


def node_from_pole(M):
    # node of the plane whose pole is M: atan2 of the pole's horizontal part
    return float(np.degrees(np.arctan2(M[0], -M[1])) % 360)


def axis_ecliptic(i_d, Om):
    """Carrier axis pole: ecliptic lon/lat and angle from invariable pole."""
    p = np.array([np.sin(np.radians(i_d)) * np.sin(np.radians(Om)),
                  -np.sin(np.radians(i_d)) * np.cos(np.radians(Om)),
                  np.cos(np.radians(i_d))])
    lon = float(np.degrees(np.arctan2(p[1], p[0])) % 360)
    lat = float(np.degrees(np.arcsin(np.clip(p[2], -1, 1))))
    ang = float(np.degrees(np.arccos(np.clip(p @ INV, -1, 1))))
    return lon, lat, ang


# --- carrier weight profiles ---
def f_logistic(a_arr, a0, wdt):
    return 1.0 / (1.0 + np.exp(-(np.asarray(a_arr, dtype=float) - a0) / wdt))


def f_bin(a_arr, w1, w2, w3):
    a_arr = np.asarray(a_arr, dtype=float)
    return np.where(a_arr < 80.0, w1, np.where(a_arr < 200.0, w2, w3))


def model_pole(a_arr, wc_arr, pd, qd):
    pG, qG, WG = giants_pq_w(a_arr)
    pf = (WG * pG + wc_arr * pd) / (WG + wc_arr)
    qf = (WG * qG + wc_arr * qd) / (WG + wc_arr)
    return pole_from_pq_vec(pf, qf)


def loglik_from_pole(M, lg):
    gamma = 10 ** lg
    mdotR = np.einsum("ij,ij->i", M, R)
    mdotJ = np.einsum("ij,ij->i", M, J)
    sinth = np.sqrt(np.clip(1 - mdotR ** 2, 0, 1))
    return float(np.sum(gamma * mdotJ - np.log(bessel_i0(gamma * sinth))))


# ---------- logistic model: x = [log10 C, a0, w, i_d, Om, log10 gamma] ----------
def log_parts(x):
    lC, a0, wdt, id_, Om, lg = x
    pd, qd = incl_vec(id_, Om)
    wc = (10 ** lC) * f_logistic(A, a0, wdt)
    M = model_pole(A, wc, pd, qd)
    ll = loglik_from_pole(M, lg)
    M65 = model_pole(np.array([65.0]),
                     (10 ** lC) * f_logistic(np.array([65.0]), a0, wdt), pd, qd)
    t65 = float(tilt_from_pole(M65)[0])
    return ll, t65


def negpost_log(x):
    ll, t65 = log_parts(x)
    return -(ll - 0.5 * (t65 / 3.0) ** 2)


# ---------- per-bin control: x = [log10 W1, log10 W2, log10 W3, i_d, Om, log10 gamma] ----------
def bin_parts(x):
    lW1, lW2, lW3, id_, Om, lg = x
    pd, qd = incl_vec(id_, Om)
    wc = f_bin(A, 10 ** lW1, 10 ** lW2, 10 ** lW3)
    M = model_pole(A, wc, pd, qd)
    ll = loglik_from_pole(M, lg)
    M65 = model_pole(np.array([65.0]),
                     f_bin(np.array([65.0]), 10 ** lW1, 10 ** lW2, 10 ** lW3),
                     pd, qd)
    t65 = float(tilt_from_pole(M65)[0])
    return ll, t65


def negpost_bin(x):
    ll, t65 = bin_parts(x)
    return -(ll - 0.5 * (t65 / 3.0) ** 2)


BOUNDS_LOG = [(np.log10(0.1), np.log10(200)), (50, 200), (1, 100),
              (0, 60), (0, 360), (np.log10(3), np.log10(100))]
BOUNDS_BIN = [(np.log10(0.001), np.log10(200)), (np.log10(0.1), np.log10(200)),
              (np.log10(0.1), np.log10(200)),
              (0, 60), (0, 360), (np.log10(3), np.log10(100))]

# --- sanity: giants alone give the invariable plane ---
for a in [65, 140, 300]:
    pG, qG, WG = giants_pq_w(np.array([float(a)]))
    MG = pole_from_pq_vec(pG, qG)
    print(f"sanity giants-only: a={a} AU tilt={tilt_from_pole(MG)[0]:.2f} deg "
          f"W_G={WG[0]:.2f}", flush=True)

print("\n--- DE: complex logistic MLE ---", flush=True)
resL = differential_evolution(negpost_log, BOUNDS_LOG, seed=7, maxiter=200,
                              popsize=20, tol=1e-7, polish=True)
xL = resL.x
llL, t65L = log_parts(xL)
print(f"negpost={resL.fun:.2f}  logL={llL:.2f}  tilt65={t65L:.2f} deg  "
      f"gamma={10**xL[5]:.1f}", flush=True)
print(f"  C={10**xL[0]:.2f} (rel. coupling), a0={xL[1]:.0f} AU, w={xL[2]:.1f} AU, "
      f"i_d={xL[3]:.1f} deg, Om={xL[4]:.0f} deg", flush=True)
lon, lat, ang = axis_ecliptic(xL[3], xL[4])
print(f"  axis d: ecliptic lon={lon:.1f} lat={lat:.1f}, "
      f"{ang:.1f} deg from invariable pole", flush=True)

print("\n--- DE: per-bin control MLE ---", flush=True)
resB = differential_evolution(negpost_bin, BOUNDS_BIN, seed=7, maxiter=200,
                              popsize=20, tol=1e-7, polish=True)
xB = resB.x
llB, t65B = bin_parts(xB)
print(f"negpost={resB.fun:.2f}  logL={llB:.2f}  tilt65={t65B:.2f} deg  "
      f"gamma={10**xB[5]:.1f}", flush=True)
print(f"  W1={10**xB[0]:.3f} (50-80), W2={10**xB[1]:.2f} (80-200), "
      f"W3={10**xB[2]:.2f} (200-400), i_d={xB[3]:.1f} deg, Om={xB[4]:.0f} deg",
      flush=True)
lon, lat, ang = axis_ecliptic(xB[3], xB[4])
print(f"  axis d: ecliptic lon={lon:.1f} lat={lat:.1f}, "
      f"{ang:.1f} deg from invariable pole", flush=True)
print(f"Delta logL (logistic minus control) = {llL - llB:.2f} (same k=6)",
      flush=True)


# ---------- MCMC ----------
def run_chain(negpost, x0, bounds, step, seed, nsteps=20000):
    r = np.random.default_rng(seed)
    LO = np.array([b[0] for b in bounds])
    HI = np.array([b[1] for b in bounds])
    x = x0 + r.normal(0, step * 0.3)
    x = np.clip(x, LO, HI)
    cur = -negpost(x)
    out = np.empty((nsteps, len(x0)))
    acc = 0
    for k in range(nsteps):
        y = x + r.normal(0, step)
        y[4] %= 360.0
        if np.any(y < LO) or np.any(y > HI):
            out[k] = x
            continue
        lp = -negpost(y)
        if np.log(r.uniform()) < lp - cur:
            x, cur = y, lp
            acc += 1
        out[k] = x
    return out, acc / nsteps


print("\n--- MCMC logistic (4 chains x 20000) ---", flush=True)
STEP_L = np.array([0.15, 10.0, 10.0, 3.0, 15.0, 0.08])
chainsL = []
for c in range(4):
    ch, acc = run_chain(negpost_log, xL, BOUNDS_LOG, STEP_L, 100 + c)
    chainsL.append(ch[3000:])
    print(f"chain {c}: acceptance {acc:.2f}", flush=True)
SL = np.vstack(chainsL)

print("\n--- MCMC per-bin control (4 chains x 20000) ---", flush=True)
STEP_B = np.array([0.15, 0.15, 0.15, 3.0, 15.0, 0.08])
chainsB = []
for c in range(4):
    ch, acc = run_chain(negpost_bin, xB, BOUNDS_BIN, STEP_B, 200 + c)
    chainsB.append(ch[3000:])
    print(f"chain {c}: acceptance {acc:.2f}", flush=True)
SB = np.vstack(chainsB)

np.savez("complex_fit_mcmc.npz", SL=SL, xL=xL, logLL=llL,
         SB=SB, xB=xB, logLB=llB)
namesL = ["log10 C", "a0", "w", "i_d", "Om", "log10 gamma"]
namesB = ["log10 W1", "log10 W2", "log10 W3", "i_d", "Om", "log10 gamma"]
print("\nposterior logistic (16/50/84%):", flush=True)
for j, nm in enumerate(namesL):
    lo, med, hi = np.percentile(SL[:, j], [16, 50, 84])
    if nm.startswith("log10 C"):
        print(f"  {nm}: {10**med:.2f}  [{10**lo:.2f}, {10**hi:.2f}]", flush=True)
    elif nm.startswith("log10 gamma"):
        print(f"  {nm}: {10**med:.1f}  [{10**lo:.1f}, {10**hi:.1f}]", flush=True)
    else:
        print(f"  {nm}: {med:.1f}  [{lo:.1f}, {hi:.1f}]", flush=True)
print("\nposterior per-bin (16/50/84%):", flush=True)
for j, nm in enumerate(namesB):
    lo, med, hi = np.percentile(SB[:, j], [16, 50, 84])
    if nm.startswith("log10 W"):
        print(f"  {nm}: {10**med:.2f}  [{10**lo:.2f}, {10**hi:.2f}]", flush=True)
    elif nm.startswith("log10 gamma"):
        print(f"  {nm}: {10**med:.1f}  [{10**lo:.1f}, {10**hi:.1f}]", flush=True)
    else:
        print(f"  {nm}: {med:.1f}  [{lo:.1f}, {hi:.1f}]", flush=True)


def tilt_curve_log(x, a_list):
    lC, a0, wdt, id_, Om, lg = x
    pd, qd = incl_vec(id_, Om)
    out = []
    for a in a_list:
        wc = (10 ** lC) * f_logistic(np.array([float(a)]), a0, wdt)
        M = model_pole(np.array([float(a)]), wc, pd, qd)
        out.append((float(tilt_from_pole(M)[0]), node_from_pole(M[0])))
    return out


def tilt_curve_bin(x, a_list):
    lW1, lW2, lW3, id_, Om, lg = x
    pd, qd = incl_vec(id_, Om)
    out = []
    for a in a_list:
        wc = f_bin(np.array([float(a)]), 10 ** lW1, 10 ** lW2, 10 ** lW3)
        M = model_pole(np.array([float(a)]), wc, pd, qd)
        out.append((float(tilt_from_pole(M)[0]), node_from_pole(M[0])))
    return out


print("\n--- predicted tilt / node (MLE) ---", flush=True)
a_list = [65, 100, 140, 180, 250, 300, 350]
cl = tilt_curve_log(xL, a_list)
cb = tilt_curve_bin(xB, a_list)
for a, (tl, nl), (tb, nb) in zip(a_list, cl, cb):
    rho_l = 2 * np.sin(np.radians(tl))
    rho_b = 2 * np.sin(np.radians(tb))
    print(f"a={a:3d} AU: logistic tilt={tl:5.2f} deg node={nl:5.0f} "
          f"(rho={rho_l:.3f}) | control tilt={tb:5.2f} deg node={nb:5.0f} "
          f"(rho={rho_b:.3f})", flush=True)

for lo_a, hi_a, lab, meas in [(80, 200, "80-200", 13.3), (200, 400, "200-400", 13.9)]:
    aa = A[(A >= lo_a) & (A <= hi_a)]
    pd, qd = incl_vec(xL[3], xL[4])
    wc = (10 ** xL[0]) * f_logistic(aa, xL[1], xL[2])
    M = model_pole(aa, wc, pd, qd)
    tl = np.mean(tilt_from_pole(M))
    pd, qd = incl_vec(xB[3], xB[4])
    wc = f_bin(aa, 10 ** xB[0], 10 ** xB[1], 10 ** xB[2])
    M = model_pole(aa, wc, pd, qd)
    tb = np.mean(tilt_from_pole(M))
    print(f"bin {lab}: logistic MLE mean tilt = {tl:.1f} deg, control = {tb:.1f} "
          f"deg (measured {meas} deg)", flush=True)

# posterior median tilt curve (logistic)
print("\n--- posterior median tilt (logistic, 16/84%) ---", flush=True)
for a in a_list:
    ts = []
    for s in SL[::40]:
        pd, qd = incl_vec(s[3], s[4])
        wc = (10 ** s[0]) * f_logistic(np.array([float(a)]), s[1], s[2])
        M = model_pole(np.array([float(a)]), wc, pd, qd)
        ts.append(float(tilt_from_pole(M)[0]))
    lo, med, hi = np.percentile(ts, [16, 50, 84])
    print(f"a={a:3d} AU: {med:.1f} [{lo:.1f}, {hi:.1f}] deg", flush=True)

# --- AIC/BIC comparison (n=46) ---
n = len(A)
rows = [("null (invariable)", 76.44, 1), ("Y-only", 79.98, 5),
        ("two-body", 80.04, 9), ("complex logistic", llL, 6),
        ("complex per-bin", llB, 6)]
print("\n--- model comparison (n=46) ---", flush=True)
for name, ll, k in rows:
    print(f"{name:18s} logL={ll:6.2f} k={k}  AIC={2*k-2*ll:7.2f}  "
          f"BIC={k*np.log(n)-2*ll:7.2f}", flush=True)
print("saved complex_fit_mcmc.npz")
print("DONE")
