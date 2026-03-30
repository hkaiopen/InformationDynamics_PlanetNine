#!/usr/bin/env python3
"""
Information Dynamics Parameter Inversion for Planet Nine
- Downloads MPCORB, filters ETNOs (a>200 AU, q>30 AU, H>5)
- Performs fixed-window clustering test (N=2 symmetry)
- Fits a mixture of two von Mises distributions (means 180° apart) to ω data
- Estimates concentration κ, information purity p_ID = κ/(1+κ), ε/γ = p/(1-p)
- Bootstrap uncertainties
- Generates polar rose plot and fitted PDF plot
"""

import os
import sys
import gzip
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import i0

# ------------------------------
# 1. Download MPCORB if needed
# ------------------------------
url = "https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz"
file_gz = "MPCORB.DAT.gz"
file_dat = "MPCORB.DAT" 
#Downloaded from https://www.minorplanetcenter.net/data/

if not os.path.exists(file_dat):
    print("Downloading MPCORB...")
    try:
        r = requests.get(url, timeout=300)
        with open(file_gz, "wb") as f:
            f.write(r.content)
        with gzip.open(file_gz, "rb") as f_in:
            with open(file_dat, "wb") as f_out:
                f_out.write(f_in.read())
        print("Downloaded and extracted.")
    except Exception as e:
        print(f"Download failed: {e}")
        sys.exit(1)
else:
    print("Using local MPCORB.DAT")

# ------------------------------
# 2. Parse fixed‑width format
# ------------------------------
colspecs = [
    (0, 7),    # designation
    (8, 13),   # H
    (14, 19),  # G
    (20, 25),  # epoch
    (26, 35),  # M
    (37, 46),  # omega
    (48, 57),  # Omega
    (59, 68),  # i
    (70, 79),  # e
    (80, 91),  # n
    (92, 103)  # a
]
names = ["desig", "H", "G", "epoch", "M", "omega", "Omega", "i", "e", "n", "a"]

df = pd.read_fwf(file_dat, colspecs=colspecs, names=names, header=None, skiprows=1)

for col in ["a", "e", "omega", "Omega", "i", "H"]:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df.dropna(subset=["a", "e", "omega", "Omega", "i"], inplace=True)
df["q"] = df["a"] * (1 - df["e"])

# ------------------------------
# 3. Filter ETNO (base criteria)
# ------------------------------
a_min = 200
q_min = 30
h_min = 5

mask = (df["a"] > a_min) & (df["q"] > q_min) & (df["H"] > h_min)
etno = df[mask].copy()
print(f"ETNO found: {len(etno)} (a > {a_min} AU, q > {q_min} AU, H > {h_min})")
if len(etno) < 5:
    print("Too few objects, aborting.")
    sys.exit(1)

etno["omega"] = etno["omega"] % 360
omegas = etno["omega"].values
n = len(omegas)

# ------------------------------
# 4. Fixed-window clustering test
# ------------------------------
FIXED_WIN = 60

def clustering_fraction(omega_p9, omegas):
    w1 = omega_p9 % 360
    w2 = (omega_p9 + 180) % 360
    cnt = 0
    for w in omegas:
        d1 = min(abs(w - w1), 360 - abs(w - w1))
        d2 = min(abs(w - w2), 360 - abs(w - w2))
        if d1 < FIXED_WIN or d2 < FIXED_WIN:
            cnt += 1
    return cnt / len(omegas)

omega_range = np.linspace(0, 360, 360)
best_frac = 0
best_omega = None
for omega in omega_range:
    frac = clustering_fraction(omega, omegas)
    if frac > best_frac:
        best_frac = frac
        best_omega = omega

print(f"\n=== Fixed-window test (half-width={FIXED_WIN}°) ===")
print(f"Best ω_P9 = {best_omega:.1f}°")
print(f"Clustering fraction = {best_frac*100:.1f}% ({int(best_frac*n)}/{n})")

np.random.seed(42)
n_sim = 10000
counts = []
for _ in range(n_sim):
    rand_omegas = np.random.uniform(0, 360, n)
    cnt = clustering_fraction(best_omega, rand_omegas) * n
    counts.append(cnt)
p_value = np.mean(np.array(counts) >= int(best_frac * n))
print(f"Monte Carlo p-value: {p_value:.4f}")

# ------------------------------
# 5. Von Mises mixture fitting
# ------------------------------
def von_mises_pdf(x, mu, kappa):
    return np.exp(kappa * np.cos(np.radians(x - mu))) / (2 * np.pi * i0(kappa))

def neg_log_likelihood(params, omegas):
    mu, kappa = params
    if kappa < 0:
        return 1e10
    pdf1 = von_mises_pdf(omegas, mu, kappa)
    pdf2 = von_mises_pdf(omegas, (mu + 180) % 360, kappa)
    pdf = 0.5 * pdf1 + 0.5 * pdf2
    pdf = np.clip(pdf, 1e-300, None)
    return -np.sum(np.log(pdf))

def fit_mixture(omegas, guess_mu=0.0, guess_kappa=1.0):
    result = minimize(neg_log_likelihood, [guess_mu, guess_kappa],
                      args=(omegas,), method='L-BFGS-B',
                      bounds=[(0, 360), (0, None)])
    mu = result.x[0] % 360
    kappa = result.x[1]
    return mu, kappa, result.fun

def bootstrap_fit(omegas, n_boot=1000):
    boot_mu = []
    boot_kappa = []
    n = len(omegas)
    for _ in range(n_boot):
        sample = np.random.choice(omegas, size=n, replace=True)
        mu, kappa, _ = fit_mixture(sample, guess_mu=0, guess_kappa=1.0)
        boot_mu.append(mu)
        boot_kappa.append(kappa)
    return np.percentile(boot_mu, [16,84]), np.percentile(boot_kappa, [16,84])

mu_opt, kappa_opt, nll = fit_mixture(omegas)
mu_ci, kappa_ci = bootstrap_fit(omegas, n_boot=500)

# Derived quantities
p_id = kappa_opt / (1 + kappa_opt)
p_id_low = kappa_ci[0] / (1 + kappa_ci[0])
p_id_high = kappa_ci[1] / (1 + kappa_ci[1])
eps_over_gamma = p_id / (1 - p_id)

print("\n=== Information Dynamics Parameter Inversion ===")
print(f"Best ω_P9 = {mu_opt:.1f}° (68% CI: [{mu_ci[0]:.1f}, {mu_ci[1]:.1f}])")
print(f"Concentration κ = {kappa_opt:.2f} (68% CI: [{kappa_ci[0]:.2f}, {kappa_ci[1]:.2f}])")
print(f"Information purity p_ID = {p_id:.3f} (68% CI: [{p_id_low:.3f}, {p_id_high:.3f}])")
print(f"Self-organization/dissipation ε/γ = {eps_over_gamma:.3f}")

# ------------------------------
# 6. Plots
# ------------------------------
# Polar rose plot
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8,8))
theta = np.radians(omegas)
ax.hist(theta, bins=36, color='skyblue', edgecolor='black', alpha=0.7)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title(f"ETNO ω distribution (n={n})")
plt.tight_layout()
plt.savefig("planet_nine_rose.png", dpi=150)
plt.close()

# Fit comparison
x_grid = np.linspace(0, 360, 500)
pdf_fit = 0.5 * von_mises_pdf(x_grid, mu_opt, kappa_opt) + \
          0.5 * von_mises_pdf(x_grid, (mu_opt+180)%360, kappa_opt)

plt.figure(figsize=(10,5))
plt.hist(omegas, bins=36, density=True, alpha=0.5, label='Data')
plt.plot(x_grid, pdf_fit, 'r-', label=f'Mixture fit (κ={kappa_opt:.2f})')
plt.xlabel('ω (degrees)')
plt.ylabel('Density')
plt.title(f'ETNO ω distribution fit with von Mises mixture\np_ID = {p_id:.3f}')
plt.legend()
plt.tight_layout()
plt.savefig("planet_nine_von_mises_fit.png", dpi=150)
plt.close()

# ------------------------------
# 7. Final report
# ------------------------------
print("\n" + "="*60)
print("FINAL REPORT – Information Dynamics Inversion for Planet Nine")
print("="*60)
print(f"Data: MPCORB.DAT (processed {pd.Timestamp.now().strftime('%Y-%m-%d')})")
print(f"ETNO criteria: a > {a_min} AU, q > {q_min} AU, H > {h_min}")
print(f"Number of ETNO: {n}")
print("\n--- Statistical clustering test ---")
print(f"Fixed window half-width: {FIXED_WIN}°")
print(f"Best ω_P9 = {best_omega:.1f}°")
print(f"Clustering fraction = {best_frac*100:.1f}% ({int(best_frac*n)}/{n})")
print(f"Monte Carlo p-value = {p_value:.4f}")
print("\n--- Information Dynamics inversion ---")
print(f"Best ω_P9 = {mu_opt:.1f}° (68% CI: [{mu_ci[0]:.1f}, {mu_ci[1]:.1f}])")
print(f"Concentration κ = {kappa_opt:.2f} (68% CI: [{kappa_ci[0]:.2f}, {kappa_ci[1]:.2f}])")
print(f"Information purity p_ID = {p_id:.3f} (68% CI: [{p_id_low:.3f}, {p_id_high:.3f}])")
print(f"Self-organization/dissipation ε/γ = {eps_over_gamma:.3f}")
print("\nInterpretation:")
print("The von Mises mixture captures the N=2 symmetry with moderate concentration.")
print(f"p_ID = {p_id:.3f} indicates a moderately self-organized state,")
print("comparable to intermediate systems. This supports the interpretation")
print("of Planet Nine as a projection of an information field.")
print("Plots saved: planet_nine_rose.png, planet_nine_von_mises_fit.png")
