#!/usr/bin/env python3
"""N-body cross-check of the two-body picture (Siraj-style).
REBOUND/WHFast: Sun + 4 giants + Planet Y + Planet Nine (MLE params from
two_body_mle.npz, else scan defaults) + the 46 Table-1 objects as massless test
particles with Omega/omega/M randomised (warp destroyed at t=0). Integrate and
measure the mean plane per bin with the Siraj estimator at snapshots.
Control run: giants only. Usage: python3 nbody_check.py [tmax_Myr] [tag]
"""
import sys, json, time
import numpy as np
import rebound
from scipy.special import i0
from scipy.optimize import minimize

TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
TAG = sys.argv[2] if len(sys.argv) > 2 else "y_p9"
BRANCH = sys.argv[3] if len(sys.argv) > 3 else "light"
MODE = sys.argv[4] if len(sys.argv) > 4 else "planets"  # planets | control

z = np.load("table1_vectors.npz")
ELEM = z["elem"][::3]  # subsample: 16 of 46 Table-1 objects (speed)
# branch params: (Y m, a, i, Om; P9 m, a, i, Om) -- P9 None means absent
BRANCHES = {
    "heavy": (dict(m=4.0, a=128.0, i=16.8, Om=119.0), None),
    "light": (dict(m=0.1, a=150.0, i=25.0, Om=120.0),
              dict(m=6.2, a=380.0, i=16.0, Om=100.0)),
}
Y, P9 = BRANCHES[BRANCH]
print(f"branch {BRANCH}: Y={Y}; P9={P9}")

GIANTS_RB = [("Jupiter", 0.0009543, 5.20, 0.0489, 1.30, 100.46, 14.75),
             ("Saturn", 0.0002857, 9.58, 0.0565, 2.49, 113.67, 92.86),
             ("Uranus", 0.0000437, 19.22, 0.0457, 0.77, 74.01, 96.94),
             ("Neptune", 0.0000515, 30.07, 0.0113, 1.77, 131.78, 276.34)]
ME = 3.003e-6  # Earth mass in solar masses
rng = np.random.default_rng(11)

# invariable plane (Souami & Souchay 2012) and (p, q) helpers (self-contained)
INV = np.array([np.cos(np.radians(88.42)) * np.cos(np.radians(17.58)),
                np.cos(np.radians(88.42)) * np.sin(np.radians(17.58)),
                np.sin(np.radians(88.42))])
INV_PQ = (np.sin(np.radians(0.8)) * np.sin(np.radians(107.58)),
          np.sin(np.radians(0.8)) * np.cos(np.radians(107.58)))
N_SYNTH = 8  # synthetic 50-80 AU particles


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


def build(with_planets=True):
    rng_tp = np.random.default_rng(12345)  # test particles: same ICs every run
    sim = rebound.Simulation()
    sim.units = ("AU", "yr", "Msun")
    sim.add(m=1.0)
    for _, m, a, e, i, Om, om in GIANTS_RB:
        sim.add(m=m, a=a, e=e, inc=np.radians(i), Omega=np.radians(Om),
                omega=np.radians(om), M=np.radians(rng.uniform(0, 360)))
    if with_planets:
        sim.add(m=Y["m"] * ME, a=Y["a"], e=0.1, inc=np.radians(Y["i"]),
                Omega=np.radians(Y["Om"]), omega=0.0, M=np.radians(rng.uniform(0, 360)))
        if P9 is not None:
            sim.add(m=P9["m"] * ME, a=P9["a"], e=0.2, inc=np.radians(P9["i"]),
                    Omega=np.radians(P9["Om"]), omega=0.0,
                    M=np.radians(rng.uniform(0, 360)))
    for (a, e, i, Om, om, M) in ELEM:
        sim.add(m=0.0, a=a, e=e, inc=np.radians(i),
                Omega=np.radians(rng_tp.uniform(0, 360)),
                omega=np.radians(rng_tp.uniform(0, 360)),
                M=np.radians(rng_tp.uniform(0, 360)))
    # synthetic 50-80 AU particles (to test the 50-80 warp in N-body)
    for _ in range(N_SYNTH):
        a = rng_tp.uniform(50, 80)
        e = rng_tp.uniform(0.05, 0.4)
        i = rng_tp.uniform(5, 25)
        sim.add(m=0.0, a=a, e=e, inc=np.radians(i),
                Omega=np.radians(rng_tp.uniform(0, 360)),
                omega=np.radians(rng_tp.uniform(0, 360)),
                M=np.radians(rng_tp.uniform(0, 360)))
    sim.integrator = "whfast"
    sim.dt = 0.5
    sim.move_to_com()
    return sim


def snapshot(sim):
    ps = sim.particles
    n_tp = len(ELEM) + N_SYNTH
    R, J, A = [], [], []
    for k in range(n_tp):
        p = ps[-(n_tp - k)]
        r = np.array([p.x, p.y, p.z])
        v = np.array([p.vx, p.vy, p.vz])
        r_hat = r / np.linalg.norm(r)
        h = np.cross(r, v)
        J_hat = h / np.linalg.norm(h)
        R.append(r_hat)
        J.append(J_hat)
        A.append(p.a)
    return np.array(R), np.array(J), np.array(A)


def siraj_quick(R, J):
    Jsum = J.sum(axis=0)

    def neglogL(x):
        p, q, lg = x
        m = pq_to_m(p, q)
        if m is None or lg > 8:
            return 1e9
        g = np.exp(lg)
        sinth = np.sqrt(np.clip(1 - (R @ m) ** 2, 0, 1))
        return -(g * (m @ Jsum) - np.sum(np.log(i0(g * sinth))))

    best = minimize(neglogL, np.array([INV_PQ[0], INV_PQ[1], np.log(15.0)]),
                    method="Nelder-Mead", options={"maxiter": 2000})
    p, q, lg = best.x
    m = pq_to_m(p, q)
    _, _, i0_deg, Om0 = m_to_pq(m)
    tilt = float(np.degrees(np.arccos(np.clip(m @ INV, -1, 1))))
    return tilt, Om0


def run(with_planets, tmax, tag):
    sim = build(with_planets)
    out = []
    times = [0.0] + [t for t in [25.0, 50.0, 100.0, 150.0, 200.0] if t <= tmax]
    for t in times:
        if t > 0:
            t0 = time.time()
            sim.integrate(t * 1e6)
            el = time.time() - t0
            print(f"  integrated to {t:.0f} Myr in {el:.0f} s", flush=True)
        R, J, A = snapshot(sim)
        row = {"t_Myr": t}
        for lo, hi, lab in [(50, 80, "50-80"), (80, 200, "80-200"), (200, 400, "200-400")]:
            sel = (A >= lo) & (A <= hi)
            tilt, Om0 = siraj_quick(R[sel], J[sel])
            row[lab] = [round(tilt, 1), round(float(Om0), 0), int(sel.sum())]
        out.append(row)
        print(f"  t={t:.0f} Myr: " + "  ".join(
            f"{k}: tilt {v[0]} deg node {v[1]:.0f} (n={v[2]})" for k, v in row.items()
            if k != "t_Myr"), flush=True)
        with open(f"nbody_{tag}_partial.json", "w") as f:
            json.dump({"mode": "planets" if with_planets else "control",
                       "snapshots": out, "Y": Y, "P9": P9}, f, indent=1)
    return out


if __name__ == "__main__":
    t0 = time.time()
    with_planets = (MODE == "planets")
    print(f"=== run TAG={TAG} BRANCH={BRANCH} MODE={MODE}, TMAX={TMAX} Myr ===", flush=True)
    res = run(with_planets, TMAX, TAG)
    with open(f"nbody_{TAG}.json", "w") as f:
        json.dump({"mode": MODE, "snapshots": res, "Y": Y, "P9": P9,
                   "wall_s": time.time() - t0}, f, indent=1)
    print(f"saved nbody_{TAG}.json, total {time.time()-t0:.0f} s")
    print("DONE")
