#!/usr/bin/env python3
"""Variants of the bodies-coupling fit: does the 200-400 AU null survive?
- V1: P9 with Siraj-like inclination i=7 deg (m=6.2, a=380).
- V2: no Planet Nine at all (Planet Y alone).
- V3: P9 omitted, light Planet Y only.
Same machinery as bodies_coupling_fit.py (imported via exec of its definitions).
"""
import numpy as np

src = open("bodies_coupling_fit.py").read()
head = src.split('print("\\n--- giants only')[0]
exec(head)

P9_LOW_I = ("Planet Nine (i=7)", 6.2, 380.0, 7.0, 100.0)


def scan(extra_fixed, label):
    print(f"\n--- {label} ---")
    for a in [65, 140, 300]:
        t, n = forced_plane(a, extra=list(extra_fixed))
        print(f"  giants+fixed at a={a}: tilt {t:.2f} deg, node {n:.0f} deg")
    n_match = 0
    for mY in [0.06, 0.1, 0.3, 1.0]:
        for aY in [100, 150, 200]:
            for iY in [10, 20, 30, 40, 60]:
                Y = ("Planet Y", mY, float(aY), float(iY), 120.0)
                t65, _ = forced_plane(65, extra=list(extra_fixed) + [Y])
                t140, n140 = forced_plane(140, extra=list(extra_fixed) + [Y])
                t300, _ = forced_plane(300, extra=list(extra_fixed) + [Y])
                if t65 < 5 and 8 <= t140 <= 20 and 90 <= n140 <= 150 and t300 < 5:
                    n_match += 1
                    print(f"  FULL MATCH (all 3 bins): m={mY} a={aY} i={iY}: "
                          f"tilt65={t65:.2f} tilt140={t140:.2f} node140={n140:.0f} tilt300={t300:.2f}")
    print(f"  full 3-bin matches: {n_match}")


scan([P9_LOW_I], "V1: P9 inclination 7 deg (Siraj-like)")
scan([], "V2: no Planet Nine (Planet Y alone)")
print("DONE")
