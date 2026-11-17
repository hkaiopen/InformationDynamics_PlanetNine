# Update — November 2026: mean-plane / carrier re-analysis

This directory holds the November 2026 re-analysis of the Planet Nine
evidence. It **supersedes** the v1.0 analysis in the repository root
(`planet_nine_inversion.py`), which used the argument of perihelion (ω)
circle and a fixed-window clustering test.

## Why superseded

Three problems in v1.0, found on re-checking (November 2026):

1. The fixed-window Monte Carlo p-value (0.0107) was invalid: the window
   centre was never re-maximised in the Monte Carlo trials. The corrected
   p-value on the same data is ≈ 0.21.
2. ω is the wrong variable for a spherical-geometry signal. Perihelion
   directions on the sphere are uniform (Rayleigh R̄ = 0.096, p = 0.63).
   The real signal is a tilt of the mean plane, not a clustering of ω.
3. `MPCORB.DAT` in the repository root is a Git LFS pointer (134 bytes),
   not the data file, so a fresh clone of v1.0 cannot run.

The v1.0 files are left untouched for the record.

## What the new analysis finds

- **Mean-plane warp.** The distant Kuiper belt's mean plane is warped by
  ~15° at 80–200 AU (node ~120°), ~2.5σ (Siraj, Chyba & Tremaine 2025,
  MNRAS 543, L27). The four giant planets alone give < 0.4°. A second
  tilt of ~10–14° at 200–400 AU (1.3–2σ) matches an inclined
  Batygin–Brown-type Planet Nine.
- **Neptune crossers.** A secular quadrupole coupling alone fails the
  crosser test (ξ test p = 0.0016, −7.8σ). Adding diffusion
  (D ≈ 0.3 AU²/Myr) and a source term (full replenishment) passes it
  (ξ p = 0.51/0.56, fraction back to ~3%).
- **Detached objects.** The same diffusion lifts perihelia into the
  detached range (Δq = 35–51 AU is 0.67–0.97σ), but the rate
  overproduces by a factor of 4–11 — an open tension, stated in the
  paper's Limits section.
- **Orbit-pole width.** The pole width bounds the angular diffusion to
  D_ang ≲ D_e. The width, not the tilt, is the binding constraint.
- **Reading.** Each contradiction constrains a different term of the
  generalised Ginzburg–Landau equation. The result is a self-consistent
  reading of Planet Nine as a purely imaginary carrier: no real body to
  find, only an information carrier with its own dynamics.

Paper: Huang, Liu & Huang (2026), *Inverting the Generalised
Ginzburg-Landau Equation from the Contradictions in the Planet Nine
Evidence*, Zenodo v5, https://doi.org/10.5281/zenodo.24304862
(submitted to Planetary Research, no. 562). The Zenodo record carries
the paper, the LaTeX source, the figures, and the full analysis
outputs (including the large snapshot files not stored here).

## Contents

| Script | What it does |
|---|---|
| `rerun_nov2026.py` | Re-runs the v1.0 ω test on current MPCORB data (corrected p-value) |
| `directional_test.py`, `spherical_test_extended.py` | Spherical-geometry tests of perihelion directions and orbit poles |
| `derivation_attempt.py` | Early spherical-geometry derivation attempt |
| `dp3_magnitudes.py` | Magnitude checks for the coupling (D-P3) analysis |
| `mean_plane_signal.py`, `warp_analysis.py` | Mean-plane warp measurement |
| `gl_axis_fit.py` | Fits the coupling-kernel axis from the warp |
| `bodies_coupling_fit.py`, `bodies_coupling_variants.py` | Secular-weights (Laplace–Lagrange) coupling fits |
| `mean_plane_measure.py`, `mean_plane_validate.py` | Mean-plane measurement and validation (200–400 AU) |
| `make_table1_cache.py` | Caches the 46 Table-1 object vectors (`table1_vectors.npz`) |
| `two_body_fit.py`, `two_body_fit_80only.py` | Joint Bayesian two-body fit; 80–200 AU-only blind test |
| `nbody_check.py` | REBOUND N-body cross-check (50 Myr, three runs) |
| `kx_matching_numbers.py` | Numbers for the K(x) kernel matching |
| `dp3b_scales.py` | Scale checks for the radial-law (D-P3b) analysis |
| `crosser_simulation.py` | Secular crosser simulation (K(x)-only) |
| `crosser_distributions.py` | Crosser a / residence-time / ω′ / q distributions |
| `crosser_samepop.py`, `xi_test.py` | Same-population re-run and the ξ/ζ test |
| `complex_fit.py` | Complex (real + imaginary part) fit of the coupling kernel |
| `crosser_loss.py` | Loss-process test (a2) — overshoots, blocked |
| `joint_fit.py` | Joint fit of carrier coupling C and Planet Nine mass m₉ |
| `crosser_diffusion.py` | Diffusion + source (D+u) test — passes the ξ test |
| `crosser_fullpop.py` | Full-population run; detached-object production rate |

Data files: `etno_table_nov2026.csv` (ETNO table), `table1_vectors.npz`
(cached Table-1 vectors), `xi_objects.json` (the 17 observed crossers),
and small `*_results.json` / `*_summary.json` outputs. Large snapshot
files (tens of MB) are **not** stored here; they are in the Zenodo
record above.

## Requirements

- Python 3.8+
- `numpy`, `pandas`, `scipy`, `matplotlib`, `emcee`, `rebound`
  (only `nbody_check.py` needs REBOUND)

Install with:

```bash
pip install numpy pandas scipy matplotlib emcee rebound
```

## Usage

Most scripts read the MPCORB data file. Download `MPCORB.DAT` from the
Minor Planet Center (https://www.minorplanetcenter.net/data/) and place
it next to the script, or point the script at your copy. Example:

```bash
python rerun_nov2026.py
python two_body_fit.py
python crosser_diffusion.py
```

See the paper and the Zenodo record for the full notes, including the
known limitations (warp 2.5σ, 200–400 AU tilt 1.3–2σ, detached rate
overproduction 4–11×, D ≈ 0.3 AU²/Myr a strong fitted coefficient).

## License

CC BY-NC-SA 4.0 (same as the repository).
