#!/usr/bin/env python3
"""One-off: cache the 46 Table-1 objects' state vectors and elements to
table1_vectors.npz so the fit and N-body scripts need not re-parse MPCORB."""
import numpy as np

exec(open("mean_plane_measure.py").read().split("base = df[")[0])
exec(open("mean_plane_validate.py").read().split("# Predicted bin-averaged")[0].split(
    "analyze(\"Table-1 subset 80-200 AU\"")[0])
# sub46: 46 matched objects
R, J, A, ELEM = [], [], [], []
for _, r in sub46.iterrows():
    rh, _, Jh = state_vectors(r["a"], r["e"], r["i"], r["Omega"], r["omega"], r["M"])
    R.append(rh)
    J.append(Jh)
    A.append(r["a"])
    ELEM.append([r["a"], r["e"], r["i"], r["Omega"], r["omega"], r["M"]])
np.savez("table1_vectors.npz", R=np.array(R), J=np.array(J), A=np.array(A),
         elem=np.array(ELEM))
print(f"saved table1_vectors.npz: n={len(A)}, a {min(A):.1f}-{max(A):.1f} AU")
print("DONE")
