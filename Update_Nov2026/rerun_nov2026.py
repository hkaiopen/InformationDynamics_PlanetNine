#!/usr/bin/env python3
"""Corrected re-run of the Planet Nine ETNO analysis on CURRENT MPCORB data.
Fixes vs the paper's code:
- Monte Carlo p-value re-maximises the window centre in every null trial
  (the paper fixed the centre found from the real data -> p=0.0107 is invalid).
- Von Mises prediction computed as the mixture probability in the UNION of
  the two windows (paper's integral as written gives 66.0%, not 83.7%).
- Bootstrap mu handled circularly (atan2 resultant), multiple MLE starts.
- Saves the full object table so the sample is reproducible.
- Also tests longitude of perihelion varpi = omega + Omega (the classic
  Batygin-Brown clustering variable) and selection-cut variants.
"""
import os, sys, gzip, requests
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import i0
from scipy.stats import circmean

WORK = os.path.dirname(os.path.abspath(__file__))
os.chdir(WORK)
URL = "https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz"

if not os.path.exists("MPCORB-202611.DAT"):
    print("Downloading MPCORB.DAT.gz ...", flush=True)
    r = requests.get(URL, timeout=600)
    print("HTTP", r.status_code, "bytes", len(r.content), flush=True)
    open("MPCORB.DAT.gz", "wb").write(r.content)
    with gzip.open("MPCORB.DAT.gz", "rb") as f_in, open("MPCORB-202611.DAT", "wb") as f_out:
        f_out.write(f_in.read())
    print("downloaded + extracted", flush=True)
else:
    print("Using local MPCORB-202611.DAT", flush=True)

with open("MPCORB-202611.DAT", "rb") as f:
    raw_head = [f.readline() for _ in range(3)]
print("First raw lines:")
for ln in raw_head:
    print(repr(ln[:120]))

colspecs = [(0, 7), (8, 13), (14, 19), (20, 25), (26, 35), (37, 46),
            (48, 57), (59, 68), (70, 79), (80, 91), (92, 103)]
names = ["desig", "H", "G", "epoch", "M", "omega", "Omega", "i", "e", "n", "a"]
df = pd.read_fwf("MPCORB-202611.DAT", colspecs=colspecs, names=names,
                 header=None, skiprows=1)
for col in ["a", "e", "omega", "Omega", "i", "H"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df.dropna(subset=["a", "e", "omega", "Omega", "i"], inplace=True)
df["q"] = df["a"] * (1 - df["e"])
df["varpi"] = (df["omega"] + df["Omega"]) % 360
print(f"Total parsed objects: {len(df)}", flush=True)

def angdiff(a, b):
    d = np.abs(a - b) % 360
    return np.minimum(d, 360 - d)

def max_frac(omegas, half=60):
    centres = np.arange(360)
    # (360, n) distances to both windows
    d1 = angdiff(omegas[None, :], centres[:, None])
    d2 = angdiff(omegas[None, :], ((centres + 180) % 360)[:, None])
    frac = np.mean((d1 < half) | (d2 < half), axis=1)
    i = int(np.argmax(frac))
    return float(centres[i]), float(frac[i])

def mc_pvalue(n, frac_obs, half=60, n_sim=10000, seed=123):
    rng = np.random.default_rng(seed)
    centres = np.arange(360)
    cnt = 0
    chunk = 500
    for s in range(0, n_sim, chunk):
        m = min(chunk, n_sim - s)
        sims = rng.uniform(0, 360, size=(m, n))
        # (m, 360, n) too big; loop centres in blocks of 60
        best = np.zeros(m)
        for c0 in range(0, 360, 60):
            cc = centres[c0:c0 + 60]
            d1 = np.minimum(np.abs(sims[:, None, :] - cc[None, :, None]) % 360,
                            360 - (np.abs(sims[:, None, :] - cc[None, :, None]) % 360))
            d2 = np.minimum(np.abs(sims[:, None, :] - ((cc + 180) % 360)[None, :, None]) % 360,
                            360 - (np.abs(sims[:, None, :] - ((cc + 180) % 360)[None, :, None]) % 360))
            f = np.mean((d1 < half) | (d2 < half), axis=2)
            best = np.maximum(best, f.max(axis=1))
        cnt += int(np.sum(best >= frac_obs - 1e-12))
    return cnt / n_sim

def vm_pdf_deg(x_deg, mu, kappa):
    # density per DEGREE (paper's code returned per-radian -> Fig.2 scale bug)
    return np.exp(kappa * np.cos(np.radians(x_deg - mu))) / (360.0 * i0(kappa))

def fit_mixture(omegas):
    def nll(p):
        mu, k = p
        if k < 0:
            return 1e10
        pdf = 0.5 * vm_pdf_deg(omegas, mu, k) + 0.5 * vm_pdf_deg(omegas, (mu + 180) % 360, k)
        return -np.sum(np.log(np.clip(pdf, 1e-300, None)))
    best = None
    for g in [0, 90]:
        res = minimize(nll, [g, 1.0], method="L-BFGS-B", bounds=[(0, 360), (0, None)])
        if best is None or res.fun < best.fun:
            best = res
    return float(best.x[0] % 360), float(best.x[1])

def mixture_union_frac(kappa, half=60):
    import mpmath as mp
    def P(a):
        aa = mp.radians(a)
        return float(mp.quad(lambda th: mp.e**(kappa * mp.cos(th)) / (2 * mp.pi * mp.besseli(0, kappa)), [-aa, aa]))
    return P(half) + (1 - P(180 - half))

def rayleigh_p(angles_deg):
    th = np.radians(angles_deg)
    R = np.abs(np.mean(np.exp(1j * th)))
    n = len(th)
    z = n * R * R
    return float(R), float(np.exp(-z) * (1 + (2 * z - z * z) / (4 * n)))

def analyse(label, sub, mc_n=4000, boot_n=200):
    om = sub["omega"].values % 360
    n = len(om)
    print(f"\n===== {label}: n={n} =====", flush=True)
    if n < 5:
        print("too few"); return
    c, f = max_frac(om)
    p = mc_pvalue(n, f, n_sim=mc_n)
    mu, k = fit_mixture(om)
    pid = k / (1 + k)
    pred60 = mixture_union_frac(k, 60)
    pred30 = mixture_union_frac(k, 30)
    R, rp = rayleigh_p(om)
    # circular bootstrap, 1000 trials
    rng = np.random.default_rng(7)
    bmu, bk = [], []
    for _ in range(boot_n):
        s = rng.choice(om, size=n, replace=True)
        try:
            m2, k2 = fit_mixture(s)
        except Exception:
            continue
        bmu.append(m2); bk.append(k2)
    bmu = np.array(bmu); bk = np.array(bk)
    mu_circ = float(np.degrees(circmean(np.radians(bmu), low=0, high=2 * np.pi)))
    klo, khi = np.percentile(bk, [16, 84])
    print(f"best centre={c:.0f} deg  fraction={f:.4f} ({int(round(f*n))}/{n})")
    print(f"corrected MC p (maximised, {mc_n} trials) = {p:.4f}")
    print(f"von Mises mixture: mu={mu:.1f} deg  kappa={k:.2f}  pID={pid:.3f}")
    print(f"bootstrap kappa 68% CI [{klo:.2f}, {khi:.2f}] -> pID [{klo/(1+klo):.3f}, {khi/(1+khi):.3f}]")
    print(f"bootstrap mu circular mean = {mu_circ:.1f} deg")
    print(f"model predicted fraction in two +-60deg windows = {pred60:.4f} (expected {pred60*n:.1f}/{n})")
    print(f"model predicted fraction in two +-30deg windows = {pred30:.4f}")
    print(f"Rayleigh (omega, unimodal): R={R:.3f} p={rp:.4f}")
    # varpi
    vp = sub["varpi"].values % 360
    Rv, rpv = rayleigh_p(vp)
    print(f"varpi: Rayleigh R={Rv:.3f} p={rpv:.4f}; resultant mean dir={float(np.degrees(np.angle(np.mean(np.exp(1j*np.radians(vp))))))%360:.1f} deg")
    return dict(label=label, n=n, centre=c, frac=f, p_mc=p, mu=mu, kappa=k,
                pid=pid, pred60=pred60, pred30=pred30, rayleigh_p=rp, varpi_rayleigh_p=rpv)

main = df[(df["a"] > 200) & (df["q"] > 30) & (df["H"] > 5)].copy()
main_sorted = main.sort_values("desig")
main_sorted[["desig", "a", "e", "q", "H", "omega", "Omega", "i", "varpi"]].to_csv(
    "etno_table_nov2026.csv", index=False)
print(f"\nMain sample (paper criteria) n={len(main)}; table saved to etno_table_nov2026.csv")
print("Designations:", ", ".join(main_sorted["desig"].astype(str).tolist()))

results = []
results.append(analyse("MAIN a>200, q>30, H>5 (paper criteria)", main, mc_n=10000, boot_n=1000))
results.append(analyse("ALT a>150, q>30, H>5", df[(df["a"] > 150) & (df["q"] > 30) & (df["H"] > 5)]))
results.append(analyse("ALT a>250, q>30, H>5", df[(df["a"] > 250) & (df["q"] > 30) & (df["H"] > 5)]))
results.append(analyse("ALT a>200, q>30, no H cut", df[(df["a"] > 200) & (df["q"] > 30)]))
results.append(analyse("ALT a>200, q>30, H>5, e<0.9 (drop extreme-e)", main[main["e"] < 0.9]))
print("\nDONE")
