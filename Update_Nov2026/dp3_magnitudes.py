#!/usr/bin/env python3
"""D-P3 magnitude arithmetic (all numbers computed, none hand-waved).

Constants: GM_sun, AU, year, solar-system age.
Coupling scale S = g*y0 (specific energy, m^2/s^2) of a direction-only
projected potential U = -S * Y20(theta) about axis mu0.
Secular precession rate scale: omega_dot ~ S / L, L = sqrt(GM*a*(1-e^2)).
Reference: Planet Nine tidal scale T9 = GM9/a9^3, rate ~ T9/n.
Mercury bound: supplementary perihelion advance -0.0036 +- 0.0050 "/cy
(Pitjeva, EPM2004); EPM2011: 0.002 "/cy; INPOP10a: 0.0004 "/cy.
Conservative bound used: |delta varpi_dot| < 0.01 "/cy.
"""
import numpy as np

GM = 1.32712440018e20      # m^3/s^2
AU = 1.495978707e11        # m
YR = 365.25 * 86400         # s
AGE = 4.5e9 * YR            # s
ARCSEC = 180 * 3600 / np.pi  # rad -> arcsec

def L_spec(a_au, e=0.0):
    return np.sqrt(GM * a_au * AU * (1 - e**2))

def rate_to_arcsec_cy(omega_dot):
    return omega_dot * (100 * YR) * ARCSEC

# --- 1. Minimum coupling that can matter at 300 AU within the age ---
L300 = L_spec(300)
S_min = 2 * np.pi * L300 / AGE
print(f"L(300 AU) = {L300:.3e} m^2/s")
print(f"S_min (period 4.5 Gyr at 300 AU) = {S_min:.2f} m^2/s^2")

# --- 2. Planet Nine reference ---
M9 = 6.6 * 3.986004418e14   # 6.6 Earth masses
a9 = 500 * AU
T9 = M9 / a9**3
n300 = np.sqrt(GM / (300 * AU)**3)
rate9 = T9 / n300
S9 = rate9 * L300
print(f"P9 tidal T9 = {T9:.2e} s^-2; secular rate at 300 AU = {rate9:.2e} s^-1 "
      f"(period {2*np.pi/rate9/YR/1e9:.2f} Gyr); equivalent S9 = {S9:.2f} m^2/s^2")

# --- 3. Mercury consequence of S_min and S9 (direction-only potential) ---
LM = L_spec(1.0, e=0.2056)
for S, tag in [(S_min, "S_min"), (S9, "S9 (P9-equivalent)")]:
    w = S / LM
    print(f"Mercury anomalous rate for {tag}: {rate_to_arcsec_cy(w):.3f} \"/cy "
          f"(bound 0.01 \"/cy -> excluded by x{rate_to_arcsec_cy(w)/0.01:.0f})")

# --- 4. Same check for a tidal (r^2) radial law: rate ~ T/n ---
T_min = 2 * np.pi * n300 / AGE   # tidal eigenvalue needed for 4.5-Gyr period at 300 AU
nM = np.sqrt(GM / (1.0 * AU)**3)
for T, tag in [(T_min, "T_min"), (T9, "T9")]:
    w = T / nM
    print(f"Mercury rate, tidal law, {tag}: {rate_to_arcsec_cy(w):.2e} \"/cy (negligible)")

# --- 5. Kozai critical inclination; ETNO inclinations ---
print(f"Kozai critical inclination = {np.degrees(np.arcsin(np.sqrt(2/5))):.2f} deg "
      f"(libration needs I above this; ETNO median i ~ 15 deg -> circulation only)")

# --- 6. Invariable-plane pole (Souami & Souchay 2012: i=1.58 deg, Omega=107.58 deg) ---
i_inv, Om_inv = np.radians(1.58), np.radians(107.58)
lon_pole = (np.degrees(Om_inv) - 90) % 360
lat_pole = 90 - 1.58
print(f"Invariable-plane pole: ecliptic lon {lon_pole:.2f} deg, lat {lat_pole:.2f} deg "
      f"(1.58 deg from ecliptic north)")

# --- 7. Mode-coefficient <-> Watson kappa (small-amplitude) ---
print(f"sqrt(5/16pi) = {np.sqrt(5/(16*np.pi)):.3f}; Watson kappa ~= 3*a*sqrt(5/16pi) "
      f"= {3*np.sqrt(5/(16*np.pi)):.3f}*a (first order); dipole Fisher kappa ~= |r|")
print("DONE")
