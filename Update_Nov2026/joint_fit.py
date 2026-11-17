#!/usr/bin/env python3
"""Joint fit: warp + crossers in ONE likelihood, fitting C (carrier
coupling amplitude, imaginary part) and m9 (Planet Nine mass, real part).

Question (user 2026-11-17): are C and m9 linked ("one thing", two faces
of the same complex entity) or independent ("two separate things")?
Test: joint posterior -- does it show a stable functional relation
C = f(m9)?  If yes -> "one thing"; if no -> "two independent"; if the
posterior is too wide to tell -> "cannot distinguish", with the reason.

Two versions of the carrier's crosser contribution (the (a2) result):
  Version N (no loss): carrier makes crossers at the saturated rate
      f_C = 9.0% (same-population), q hugging 30 (median 28.7 AU).
      The carrier's contribution is large, so the mixture q-shape and
      the total fraction DO depend on C (via the C_* threshold below).
  Version L (with loss, tau=100 Myr, the physical one): the carrier's
      crossers decay away; today (last 1 Gyr) f_C ~ 0.15%, q median
      29.3-29.8 AU. The carrier contributes ~nothing to TODAY's
      crossers, so the crosser likelihood is ~independent of C.
      (a2) showed: no-loss 9% (3x over) -> with-loss 0.7% (4x under),
      q-shape worse, no steady state. Path (a2) is dead; this version
      uses its numbers as the physical carrier contribution.

Warp likelihood: Siraj conditional von Mises on the 46 Table-1 objects
(reused from complex_fit.py; depends on C, a0, w, i_d, Om, gamma only).
Crosser likelihood: KDE likelihood of the 17 observed crossers' q values
(with Batygin geometric discovery-bias weighting, xi_test.py machinery)
under a MIXTURE q-distribution with TWO effective components:
    F(q) = [f_H * F_H(q) + f_F * F_9(q)] / (f_H + f_F)
  Hugging-30 component (F_H = carrier snapshots; version N: samepop
       A_100_400_S_9; version L: loss tau100 surviving snapshots):
       carrier (secular) + Neptune-scattering baseline (no-P9, 0.5%).
       The baseline's q-shape is approximated by the carrier's snapshots
       (both hug 30 AU; the baseline is likely even MORE concentrated at
       30, so this slightly understates the m9=0 penalty -- conservative).
  Flat component (F_9): Planet Nine scattering, flat U[16,30] + a low-q
       tail 0.15*U[2,16] (Batygin's model is flat 16-30 but must extend
       below 16: three observed crossers have q = 2.57, 10.56, 7.88 AU;
       the flat-only form would assign them ~zero probability).
  f_F(m9) = 0.005*m9  (linear in mass; 0 at m9=0, 0.025 at m9=5;
       assumption, sensitivity-tested).  NOTE: at m9=0 there is NO flat
       component -- the no-P9 baseline is hugging-30, NOT flat.  (An
       earlier draft wrongly gave the flat shape to the m9=0 baseline.)
  f_H(C): version N: f_C(C) + 0.005, f_C = f_sat*logistic((log10C -
       log10C_*)/0.3), f_sat = 0.090, C_* = threshold C where the
       carrier's secular period ~ 4.5 Gyr (i.e. S = S_min); C_* is a free
       parameter (it encodes the unknown C<->S map, D-P2a).
       version L: f_H = 0.0015 + 0.005 = 0.0065 fixed (today's decayed
       carrier value; the carrier contributes ~nothing to today's
       crossers whatever C is, so C_* is not needed and is dropped).
  Fraction likelihood (secondary, Gaussian):
       logL_frac = -0.5*((f_H + f_F - 0.03)/0.015)^2
       (3% = Batygin P9 steady state; sigma 1.5% is a judgement).

Parameters:
  version N: [log10 C, a0, w, i_d, Om, log10 gamma, log10 m9, log10 C_*]
  version L: [log10 C, a0, w, i_d, Om, log10 gamma, log10 m9]
Priors are uniform in the stated ranges (same warp bounds as
complex_fit.py; m9 log-uniform [0.5, 10] M_E as in two_body_fit.py).
"""
import numpy as np
import json, time, sys
from scipy.special import i0 as bessel_i0
from scipy.optimize import differential_evolution

T0 = time.time()
def stamp(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)

# ================= warp data + machinery (from complex_fit.py) =================
z = np.load("table1_vectors.npz")
R, J, A = z["R"], z["J"], z["A"]
stamp(f"warp data: n={len(A)} objects, a {A.min():.1f}-{A.max():.1f} AU")

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
    a_arr = np.asarray(a_arr, dtype=float)
    num = np.zeros((len(a_arr), 2)); den = np.zeros(len(a_arr))
    for k in range(4):
        interior = G_A[k] < a_arr
        alpha = np.where(interior, G_A[k] / a_arr, a_arr / G_A[k])
        wj = np.where(interior, G_M[k] * alpha * b_of(alpha),
                      G_M[k] * alpha ** 2 * b_of(alpha))
        num[:, 0] += wj * G_PQ[k, 0]; num[:, 1] += wj * G_PQ[k, 1]; den += wj
    return num[:, 0] / den, num[:, 1] / den, den

def pole_from_pq_vec(p, q):
    s = np.hypot(p, q)
    ia = 2 * np.arcsin(np.clip(s, 0, 1))
    Om0 = np.arctan2(p, q)
    return np.column_stack([np.sin(ia) * np.sin(Om0),
                            -np.sin(ia) * np.cos(Om0), np.cos(ia)])

def tilt_from_pole(M):
    return np.degrees(np.arccos(np.clip(M @ INV, -1, 1)))

def f_logistic(a_arr, a0, wdt):
    return 1.0 / (1.0 + np.exp(-(np.asarray(a_arr, dtype=float) - a0) / wdt))

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

def warp_loglik(xw):
    """xw = [log10 C, a0, w, i_d, Om, log10 gamma]; returns (logL, tilt65)."""
    lC, a0, wdt, id_, Om, lg = xw
    pd, qd = incl_vec(id_, Om)
    wc = (10 ** lC) * f_logistic(A, a0, wdt)
    M = model_pole(A, wc, pd, qd)
    ll = loglik_from_pole(M, lg)
    M65 = model_pole(np.array([65.0]),
                     (10 ** lC) * f_logistic(np.array([65.0]), a0, wdt), pd, qd)
    t65 = float(tilt_from_pole(M65)[0])
    return ll - 0.5 * (t65 / 3.0) ** 2, t65

# ================= crosser data: 17 objects =================
objs = json.load(open("xi_objects.json"))
QJ = np.array([o["q"] for o in objs])
RJ = np.array([o["r_disc"] for o in objs])
stamp(f"crosser data: {len(objs)} objects, q median {np.median(QJ):.1f} AU, "
      f"q<25: {(QJ<25).sum()}/17")

# ================= carrier snapshots =================
res = json.load(open("crosser_samepop_results.json"))
dN = res["A_100_400_S_9"]
saN = np.array(dN["snap_a"]); sqN = np.array(dN["snap_q"]); siN = np.array(dN["snap_i"])
mN = (saN > 100) & (sqN < 30) & (siN < 40)
saN, sqN = saN[mN], sqN[mN]
stamp(f"carrier version N snapshots (A_100_400_S_9, cuts a>100 q<30 i<40): "
      f"{len(saN)}, q median {np.median(sqN):.2f}, q<20 {(sqN<20).mean()*100:.1f}%")

zL = np.load("crosser_loss_snapshots.npz")
saL = zL["A_100_400_S_9_tau100_a"]; sqL = zL["A_100_400_S_9_tau100_q"]
siL = zL["A_100_400_S_9_tau100_i"]
mL = (saL > 100) & (sqL < 30) & (siL < 40)
saL, sqL = saL[mL], sqL[mL]
stamp(f"carrier version L snapshots (A_100_400_S_9 tau100 surviving, same cuts): "
      f"{len(saL)}, q median {np.median(sqL):.2f}, q<20 {(sqL<20).mean()*100:.1f}%")

# ================= synthetic P9 snapshots =================
rng = np.random.default_rng(20261117)
N9 = 20000
# a ~ a^(-3/2) on [100, 5000] (Batygin range, population form)
p = -0.5
u = rng.uniform(0, 1, N9)
sa9 = (100.0 ** p + u * (5000.0 ** p - 100.0 ** p)) ** (1.0 / p)
# q: 0.85*U[16,30] + 0.15*U[2,16] (flat + low-q tail; see header)
is_tail = rng.uniform(0, 1, N9) < 0.15
sq9 = np.empty(N9)
sq9[~is_tail] = rng.uniform(16.0, 30.0, int((~is_tail).sum()))
sq9[is_tail] = rng.uniform(2.0, 16.0, int(is_tail.sum()))
si9 = np.abs(rng.normal(0.0, 15.0, N9))
m9 = (sa9 > 100) & (sq9 < 30) & (si9 < 40)
sa9, sq9 = sa9[m9], sq9[m9]
stamp(f"P9 synthetic snapshots: {len(sa9)}, q median {np.median(sq9):.2f}, "
      f"q<16 {(sq9<16).mean()*100:.1f}%")

# ================= geometric weight + KDE precomputation =================
def period_fraction_in_window(a, q, rlo, rhi):
    """Fraction of orbital period with heliocentric distance in [rlo, rhi].
    Exact via Kepler's equation (copied from xi_test.py)."""
    a = np.asarray(a, dtype=float); q = np.asarray(q, dtype=float)
    e = 1.0 - q / a
    Q = a * (1.0 + e)
    lo = np.maximum(rlo, q); hi = np.minimum(rhi, Q)
    w = np.zeros_like(a)
    ok = hi > lo
    if np.any(ok):
        ao, eo = a[ok], e[ok]
        def E_of_r(r, ao_, eo_):
            c = (1.0 - r / ao_) / eo_
            return np.arccos(np.clip(c, -1.0, 1.0))
        E1 = E_of_r(lo[ok], ao, eo); E2 = E_of_r(hi[ok], ao, eo)
        w[ok] = ((E2 - eo * np.sin(E2)) - (E1 - eo * np.sin(E1))) / np.pi
    return w

def precompute(sa, sq, h=1.5):
    """For each of the 17 objects j: A_j = mean over snapshots of
    w_j * K_h(q - q_j), Cw_j = mean of w_j. Returns (A, Cw) arrays len 17."""
    n = len(sa)
    Aj = np.zeros(len(QJ)); Cj = np.zeros(len(QJ))
    norm = 1.0 / (h * np.sqrt(2 * np.pi))
    for j in range(len(QJ)):
        w = period_fraction_in_window(sa, sq, RJ[j] - 1.0, RJ[j] + 1.0)
        K = norm * np.exp(-0.5 * ((sq - QJ[j]) / h) ** 2)
        Aj[j] = float(np.mean(w * K))
        Cj[j] = float(np.mean(w))
    return Aj, Cj

stamp("precomputing KDE/geometry sums (h=1.5 AU) ...")
AN, CN = precompute(saN, sqN)
AL, CL = precompute(saL, sqL)
A9, C9 = precompute(sa9, sq9)
stamp(f"precompute done. Cw (mean weight) ranges: carrier-N "
      f"[{CN.min():.3g},{CN.max():.3g}], carrier-L [{CL.min():.3g},{CL.max():.3g}], "
      f"P9 [{C9.min():.3g},{C9.max():.3g}]")
np.savez("joint_fit_precomp.npz", AN=AN, CN=CN, AL=AL, CL=CL, A9=A9, C9=C9,
         QJ=QJ, RJ=RJ)

# ================= crosser likelihood =================
# Two effective components for the q-shape:
#   Hugging-30: carrier (secular) + Neptune-scattering baseline (no-P9).
#       The baseline fraction is 0.005 (Batygin no-P9 model, 0.5%) and
#       its q-shape is approximated by the carrier's snapshots (both
#       hug 30 AU; the baseline is likely even MORE concentrated at 30,
#       so this slightly UNDERSTATES the m9=0 penalty -- conservative
#       for the m9 constraint).
#   Flat: P9 scattering. f_F(m9) = 0.005*m9 (0 at m9=0, 0.025 at m9=5);
#       q ~ 0.85*U[16,30] + 0.15*U[2,16] (flat + low-q tail).
# Total fraction f_tot = f_H + f_F; target 3% (Batygin P9 steady state).
F_SAT_N = 0.090   # carrier saturated fraction, no loss (same-population)
F_BASE = 0.005    # Neptune-scattering baseline (no-P9), hugging-30
F_C_L = 0.0015    # carrier today fraction with loss tau=100 Myr
                  # (last-1-Gyr mean over 4 groups: 0.0003-0.0027)
F_SIG = 0.015     # sigma for the fraction Gaussian (judgement)

def fF_of_m9(m9):
    """P9-scattering (flat) fraction vs mass: 0 at m9=0, 0.025 at m9=5.
    Linear in mass (scattering rate ~ gravitational cross-section);
    assumption, sensitivity-tested."""
    return 0.005 * m9

def crosser_loglik(fH, m9, A_c, C_c):
    """Mixture KDE log-likelihood of the 17 q values + fraction Gaussian.
    fH: hugging-30 fraction (carrier + baseline); m9: P9 mass;
    (A_c, C_c): hugging-component precomputed sums; flat sums (A9, C9)
    are global."""
    fF = fF_of_m9(m9)
    ftot = fH + fF
    if ftot <= 0:
        return -np.inf
    # mixture density at each object: (fH*A_c + fF*A9)/(fH*C_c + fF*C9)
    num = fH * A_c + fF * A9
    den = fH * C_c + fF * C9
    with np.errstate(divide="ignore", invalid="ignore"):
        dens = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)
    if np.any(dens <= 0) or not np.all(np.isfinite(dens)):
        return -1e6  # model assigns ~zero probability to some object
    ll_q = float(np.sum(np.log(dens)))
    ll_f = -0.5 * ((ftot - 0.03) / F_SIG) ** 2
    return ll_q + ll_f

def fC_N(log10C, log10Cs):
    """Version N carrier fraction: saturated above threshold C_*, 0 below;
    smooth logistic step in log10 C, width 0.3 dex."""
    return F_SAT_N / (1.0 + 10.0 ** (-(log10C - log10Cs) / 0.3))

# ---- sanity: crosser-only likelihoods at anchor points ----
# (fH = carrier fraction + 0.005 baseline)
stamp("--- crosser-only sanity checks (fH = carrier + 0.005 baseline) ---")
for lab, fH, m9, Ac, Cc in [
        ("carrier only N (fH=0.095, m9=0)", 0.095, 0.0, AN, CN),
        ("P9 only (fH=0.005, m9=5)", 0.005, 5.0, AN, CN),
        ("no-P9 baseline only (fH=0.005, m9=0)", 0.005, 0.0, AN, CN),
        ("mixture N (fH=0.095, m9=5)", 0.095, 5.0, AN, CN),
        ("mixture L (fH=0.0065, m9=5)", 0.0065, 5.0, AL, CL),
        ("mixture L (fH=0.0065, m9=0)", 0.0065, 0.0, AL, CL)]:
    ll = crosser_loglik(fH, m9, Ac, Cc)
    stamp(f"  {lab}: crosser logL = {ll:.2f}")

# ---- crosser-only best fit over m9 (version L hugging fH=0.0065) ----
grid_m9 = np.linspace(0.0, 10.0, 101)
ll_grid = [crosser_loglik(0.0065, m, AL, CL) for m in grid_m9]
i_best = int(np.argmax(ll_grid))
stamp(f"crosser-only (version L) best m9 = {grid_m9[i_best]:.1f} M_E, "
      f"logL = {ll_grid[i_best]:.2f}; m9=5 logL = "
      f"{crosser_loglik(0.0065, 5.0, AL, CL):.2f}; m9=0 logL = "
      f"{crosser_loglik(0.0065, 0.0, AL, CL):.2f}")
ll_grid_N = [crosser_loglik(0.095, m, AN, CN) for m in grid_m9]
i_best_N = int(np.argmax(ll_grid_N))
stamp(f"crosser-only (version N, fH=0.095) best m9 = {grid_m9[i_best_N]:.1f} M_E, "
      f"logL = {ll_grid_N[i_best_N]:.2f}")

# ================= joint MCMC =================
# Parameters:
#   version N: x = [log10C, a0, w, i_d, Om, log10gamma, log10m9, log10C_*]
#   version L: x = [log10C, a0, w, i_d, Om, log10gamma, log10m9]
# Priors: uniform in the ranges below (same warp bounds as complex_fit.py;
# m9 log-uniform [0.5, 10] M_E as in two_body_fit.py).
BOUNDS_N = np.array([[-1.0, np.log10(200.0)],   # log10 C
                     [50.0, 200.0],             # a0
                     [1.0, 100.0],              # w
                     [0.0, 60.0],               # i_d
                     [0.0, 360.0],              # Om
                     [np.log10(3.0), 2.0],      # log10 gamma
                     [np.log10(0.5), 1.0],      # log10 m9
                     [-1.0, 3.0]])              # log10 C_*
BOUNDS_L = BOUNDS_N[:7]
STEP_N = np.array([0.15, 10.0, 10.0, 3.0, 15.0, 0.08, 0.15, 0.30])
STEP_L = STEP_N[:7]
X0_N = np.array([np.log10(79.52), 76.0, 1.0, 14.8, 119.0, np.log10(15.2),
                 np.log10(5.0), 1.0])
X0_L = X0_N[:7]

def joint_logpost(x, version):
    """Joint log-posterior (log-likelihood + uniform priors in bounds).
    version 'N': crosser fH = fC_N(log10C, log10C_*) + 0.005 (carrier
        snapshots N).  version 'L': fH = 0.0065 fixed (carrier snapshots L);
        the crosser likelihood is then ~independent of C."""
    if version == "N":
        B = BOUNDS_N
    else:
        B = BOUNDS_L
    if np.any(x < B[:, 0]) or np.any(x > B[:, 1]):
        return -np.inf
    lC, a0, wdt, id_, Om, lg, lm9 = x[:7]
    m9 = 10.0 ** lm9
    ll_w, _ = warp_loglik(x[:6])
    if version == "N":
        fH = fC_N(lC, x[7]) + F_BASE
        ll_c = crosser_loglik(fH, m9, AN, CN)
    else:
        ll_c = crosser_loglik(0.0065, m9, AL, CL)
    ll = ll_w + ll_c
    if not np.isfinite(ll):
        return -np.inf
    return ll

def run_mcmc(version, nchain=4, nstep=20000, seed=12345):
    """Metropolis MCMC, same style as complex_fit.py. Returns
    (samples, acc_rates, best_x, best_ll)."""
    B = BOUNDS_N if version == "N" else BOUNDS_L
    step = STEP_N if version == "N" else STEP_L
    x0 = X0_N if version == "N" else X0_L
    rng = np.random.default_rng(seed)
    nd = len(x0)
    burn, thin = 5000, 5
    nkeep = (nstep - burn) // thin
    all_s = np.zeros((nchain, nkeep, nd))
    accs = []
    best_x, best_ll = x0.copy(), joint_logpost(x0, version)
    for c in range(nchain):
        x = x0 + rng.normal(0, 0.3 * step)
        x = np.clip(x, B[:, 0] + 1e-9, B[:, 1] - 1e-9)
        ll = joint_logpost(x, version)
        acc = 0
        kept = 0
        for s in range(nstep):
            xp = x + rng.normal(0, step)
            # wrap Om into [0, 360)
            xp[4] = xp[4] % 360.0
            llp = joint_logpost(xp, version)
            if np.isfinite(llp) and np.log(rng.uniform()) < llp - ll:
                x, ll = xp, llp
                acc += 1
            if ll > best_ll:
                best_x, best_ll = x.copy(), ll
            if s >= burn and (s - burn) % thin == 0:
                all_s[c, kept] = x
                kept += 1
        accs.append(acc / nstep)
        stamp(f"  version {version} chain {c+1}: acc {acc/nstep:.2f}, "
              f"best ll {best_ll:.2f}")
    return all_s, accs, best_x, best_ll

stamp("--- joint MCMC version N (no loss) ---")
sN, accN, bestN, bestllN = run_mcmc("N")
stamp(f"version N done. best logL = {bestllN:.2f} at "
      f"log10C={bestN[0]:.2f} a0={bestN[1]:.0f} w={bestN[2]:.1f} "
      f"i_d={bestN[3]:.1f} Om={bestN[4]:.0f} log10g={bestN[5]:.2f} "
      f"m9={10**bestN[6]:.2f} log10C_*={bestN[7]:.2f}")
stamp("--- joint MCMC version L (with loss) ---")
sL, accL, bestL, bestllL = run_mcmc("L")
stamp(f"version L done. best logL = {bestllL:.2f} at "
      f"log10C={bestL[0]:.2f} a0={bestL[1]:.0f} w={bestL[2]:.1f} "
      f"i_d={bestL[3]:.1f} Om={bestL[4]:.0f} log10g={bestL[5]:.2f} "
      f"m9={10**bestL[6]:.2f}")
np.savez("joint_fit_mcmc.npz", sN=sN, sL=sL, accN=accN, accL=accL,
         bestN=bestN, bestllN=bestllN, bestL=bestL, bestllL=bestllL)

# ================= post-processing =================
def summarise(s, version):
    """Posterior summary for one version. s: (nchain, nkeep, nd)."""
    flat = s.reshape(-1, s.shape[2])
    lC = flat[:, 0]; lm9 = flat[:, 6]
    C = 10.0 ** lC; m9 = 10.0 ** lm9
    r = float(np.corrcoef(lC, lm9)[0, 1])
    # log-log regression: log10 C = a + b * log10 m9
    Areg = np.column_stack([np.ones_like(lm9), lm9])
    coef, *_ = np.linalg.lstsq(Areg, lC, rcond=None)
    resid = lC - Areg @ coef
    dof = len(lC) - 2
    s2 = float(np.sum(resid ** 2) / dof)
    cov = s2 * np.linalg.inv(Areg.T @ Areg)
    b_se = float(np.sqrt(cov[1, 1]))
    R2 = float(1.0 - np.sum(resid ** 2) / np.sum((lC - lC.mean()) ** 2))
    # ratio C/m9 in dex
    lr = lC - lm9
    # fH implied (version N only): carrier fraction at the posterior
    out = dict(version=version, r=r, slope=float(coef[1]),
               slope_se=b_se, intercept=float(coef[0]), R2=R2,
               C_med=float(np.median(C)), C_lo=float(np.percentile(C, 16)),
               C_hi=float(np.percentile(C, 84)),
               m9_med=float(np.median(m9)), m9_lo=float(np.percentile(m9, 16)),
               m9_hi=float(np.percentile(m9, 84)),
               ratio_med=float(np.median(10.0 ** lr)),
               ratio_lo=float(np.percentile(10.0 ** lr, 16)),
               ratio_hi=float(np.percentile(10.0 ** lr, 84)),
               ratio_dex_width=float(np.percentile(lr, 84) - np.percentile(lr, 16)))
    if version == "N":
        fC = fC_N(lC, flat[:, 7])
        out["fC_med"] = float(np.median(fC))
        out["fC_lo"] = float(np.percentile(fC, 16))
        out["fC_hi"] = float(np.percentile(fC, 84))
        out["frac_Cstar_above_C"] = float(np.mean(flat[:, 7] > lC))
    return out

def report_summary(out):
    stamp(f"--- posterior summary version {out['version']} ---")
    stamp(f"  Pearson r(log10C, log10m9) = {out['r']:+.3f}")
    stamp(f"  log-log fit: log10C = {out['intercept']:.2f} + "
          f"{out['slope']:+.2f}(+/-{out['slope_se']:.2f}) * log10m9, R2 = {out['R2']:.3f}")
    stamp(f"  C: median {out['C_med']:.1f} [{out['C_lo']:.1f}, {out['C_hi']:.1f}]")
    stamp(f"  m9: median {out['m9_med']:.2f} M_E [{out['m9_lo']:.2f}, {out['m9_hi']:.2f}]")
    stamp(f"  C/m9: median {out['ratio_med']:.1f} [{out['ratio_lo']:.1f}, "
          f"{out['ratio_hi']:.1f}], 16-84% width in dex = {out['ratio_dex_width']:.2f}")
    if out['version'] == 'N':
        stamp(f"  carrier fraction fC: median {out['fC_med']:.4f} "
              f"[{out['fC_lo']:.4f}, {out['fC_hi']:.4f}]")
        stamp(f"  P(C_* > C) = {out['frac_Cstar_above_C']:.2f} "
              f"(threshold above the fitted C: carrier ON for warp, OFF for crossers)")

sumN = summarise(sN, "N")
report_summary(sumN)
sumL = summarise(sL, "L")
report_summary(sumL)

# ---- max logL comparison: joint best vs sum of separate bests ----
# warp-only max (complex_fit.py): 80.15. Crosser-only max from the grid.
crosser_max_L = max(ll_grid)
crosser_max_N = max(ll_grid_N)
stamp("--- max logL comparison ---")
stamp(f"  warp-only max logL = 80.15 (complex_fit.py)")
stamp(f"  crosser-only max logL: version L {crosser_max_L:.2f}, "
      f"version N {crosser_max_N:.2f}")
stamp(f"  sum of separate maxima: version L {80.15 + crosser_max_L:.2f}, "
      f"version N {80.15 + crosser_max_N:.2f}")
stamp(f"  joint best logL: version L {bestllL:.2f}, version N {bestllN:.2f}")
stamp(f"  joint - sum of separate: version L {bestllL - (80.15 + crosser_max_L):.2f}, "
      f"version N {bestllN - (80.15 + crosser_max_N):.2f}")
json.dump({"sumN": sumN, "sumL": sumL,
           "bestllN": float(bestllN), "bestllL": float(bestllL),
           "crosser_max_L": float(crosser_max_L),
           "crosser_max_N": float(crosser_max_N)},
          open("joint_fit_summary.json", "w"), indent=1)
stamp("saved joint_fit_mcmc.npz, joint_fit_summary.json")
stamp("DONE")
