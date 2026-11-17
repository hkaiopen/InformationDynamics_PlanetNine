#!/usr/bin/env python3
"""K(x) matching numbers: tilt curve shape per case + anisotropy ratios.
Reuses the secular weight law from bodies_coupling_fit.py (Laplace weights,
borrowed from secular theory). Purpose: check whether the tilt curve is
peaked (single perturber) or flat (two perturbers / carrier step), and
compute h_perp/A = 2 sin(theta) for the measured tilts.
"""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

def b132(alpha):
    if alpha <= 0:
        return 0.0
    f = lambda psi: np.cos(psi) / (1 - 2*alpha*np.cos(psi) + alpha**2)**1.5
    val, _ = quad(f, 0, np.pi, limit=200)
    return 2*val/np.pi

agrid = np.concatenate([np.linspace(0.01, 0.95, 95), [0.97, 0.99]])
bgrid = np.array([b132(a) for a in agrid])
binterp = interp1d(agrid, bgrid, kind="cubic", fill_value="extrapolate")

INV_P = np.sin(np.radians(0.8))*np.sin(np.radians(107.58))
INV_Q = np.sin(np.radians(0.8))*np.cos(np.radians(107.58))

def incl_vec(i_deg, Om_deg):
    return np.array([np.sin(np.radians(i_deg/2))*np.sin(np.radians(Om_deg)),
                     np.sin(np.radians(i_deg/2))*np.cos(np.radians(Om_deg))])

def weight(m, a_j, a):
    if a_j < a:
        alpha = a_j/a
        return m*alpha*float(binterp(alpha))
    else:
        alpha = a/a_j
        return m*alpha**2*float(binterp(alpha))

GIANTS = [("Jupiter",317.8,5.20,1.30,100.46),("Saturn",95.2,9.58,2.49,113.67),
          ("Uranus",14.5,19.22,0.77,74.01),("Neptune",17.1,30.07,1.77,131.78)]
P9 = ("P9 BB21",6.2,380.0,16.0,100.0)
Y = ("Y light",0.1,150.0,25.0,120.0)

def forced_plane(a, extra=()):
    bodies = GIANTS + list(extra)
    num = np.zeros(2); den = 0.0
    for _, m, a_j, i, Om in bodies:
        w = weight(m, a_j, a)
        num += w*incl_vec(i, Om); den += w
    p, q = num/den
    tilt = 2*np.degrees(np.arcsin(np.clip(np.hypot(p-INV_P, q-INV_Q),0,1)))
    node = float(np.degrees(np.arctan2(p,q)) % 360)
    return tilt, node

print("tilt(a) / node(a)  [deg]  for each case (linear secular, alpha capped by grid)")
print(f"{'a':>5} {'giants':>14} {'+Y':>14} {'+P9':>14} {'+Y+P9':>14} {'MLE fit':>8}")
mle = {65:1.2,100:11.1,140:14.9,180:13.2,250:11.8,300:12.2,350:14.2}
for a in [65,100,140,180,250,300,350]:
    t0,_ = forced_plane(a)
    tY,nY = forced_plane(a, extra=[Y])
    tP,nP = forced_plane(a, extra=[P9])
    tYP,nYP = forced_plane(a, extra=[Y,P9])
    print(f"{a:5d} {t0:6.2f}        {tY:6.2f}/{nY:3.0f} {tP:6.2f}/{nP:3.0f} {tYP:6.2f}/{nYP:3.0f} {mle[a]:8.1f}")

print("\nanisotropy ratio h_perp/A = 2 sin(theta)  (vector model, h_z=0)")
for label, th in [("80-200 measured 13.3",13.3),("80-200 MLE bin 11.3",11.3),
                  ("200-400 measured 11.7 (n=29)",11.7),("200-400 Table-1 13.9",13.9),
                  ("200-400 MLE bin 13.3",13.3),("50-80 measured 3.6",3.6),
                  ("warp 13.46 (gl_axis_fit)",13.46)]:
    print(f"  {label:32s} theta={th:5.2f}  h/A = {2*np.sin(np.radians(th)):.3f}")

print("\nweight profile of Y alone (m=0.1, a=150) vs giants total, at each a:")
for a in [65,100,140,180,250,300,350]:
    wY = weight(0.1,150.0,a)
    wG = sum(weight(m,aj,a) for _,m,aj,_,_ in GIANTS)
    wP9 = weight(6.2,380.0,a)
    print(f"  a={a:3d}: wY={wY:8.3f}  wGiants={wG:8.3f}  wP9={wP9:8.3f}  wY/wG={wY/wG:6.3f}  wP9/wG={wP9/wG:6.3f}")
print("DONE")
