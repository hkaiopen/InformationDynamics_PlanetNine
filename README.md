# Information Dynamics Inversion for Planet Nine

This repository contains a Python implementation to invert the **Information Dynamics** parameters for the hypothetical Planet Nine using the argument of perihelion (ω) distribution of extreme trans-Neptunian objects (ETNOs). The code downloads the latest MPCORB database, selects ETNOs, performs a fixed‑window clustering test, and fits a mixture of two von Mises distributions to extract the concentration parameter κ, information purity p_ID = κ/(1+κ), and the self‑organization to dissipation ratio ε/γ = p_ID/(1−p_ID).

---

## 1. Purpose

The aim is to provide a quantitative, data‑driven estimate of the dynamical state of the outer Solar System within the Information Dynamics framework. Instead of assuming a conventional massive planet, the observed N=2 symmetry in ETNO orbits is interpreted as a projection of an information field. The inverted parameters offer a direct, testable characterization of this field.

---

## 2. Method

1. **Data retrieval**  
   The script automatically downloads the Minor Planet Center’s `MPCORB.DAT` file (gzipped) if not already present.

2. **ETNO selection**  
   Objects are selected with:
   - semi‑major axis \(a > 200\) AU  
   - perihelion distance \(q = a(1-e) > 30\) AU  
   - absolute magnitude \(H > 5\)

3. **Statistical clustering test**  
   A fixed window half‑width of 60° is used to measure the fraction of ETNOs whose ω falls within two opposite windows (centered on a test longitude ωₚ₉ and ωₚ₉+180°). The best‑fitting ωₚ₉ is found by scanning 0–360°. Significance is assessed with 10,000 Monte Carlo trials against a uniform distribution.

4. **Information Dynamics parameter inversion**  
   The ω distribution is modelled as a mixture of two von Mises distributions with equal concentration κ and means separated by 180°. Maximum likelihood estimation yields κ.  
   The information purity is defined as  
   \[
   p_{\mathrm{ID}} = \frac{\kappa}{1+\kappa}
   \]  
   and the self‑organization to dissipation ratio as  
   \[
   \frac{\varepsilon}{\gamma} = \frac{p_{\mathrm{ID}}}{1-p_{\mathrm{ID}}}.
   \]  
   Bootstrap resampling (1000 trials) gives 68% confidence intervals.

5. **Visualization**  
   - A polar rose plot of the ω distribution (`planet_nine_rose.png`).  
   - A histogram with the fitted von Mises mixture density (`planet_nine_von_mises_fit.png`).

---

## 3. Requirements

- Python 3.6+
- Required packages: `numpy`, `pandas`, `matplotlib`, `scipy`, `requests`

Install them with:

```bash
pip install numpy pandas matplotlib scipy requests
```

---

## 4. Usage

Clone the repository and run the script:

```bash
git clone https://github.com/hkaiopen/InformationDynamics_PlanetNine.git
cd InformationDynamics_PlanetNine
python planet_nine_inversion.py
```

The script will:
- Download `MPCORB.DAT` if missing (about 30–40 MB, takes a few seconds).
- Process the data and output results to the console.
- Save the two plots in the current directory.

---

## 5. Output example

```
ETNO found: 43 (a > 200 AU, q > 30 AU, H > 5)

=== Fixed-window test (half-width=60°) ===
Best ω_P9 = 356.0°
Clustering fraction = 83.7% (36/43)
Monte Carlo p-value: 0.0107

=== Information Dynamics Parameter Inversion ===
Best ω_P9 = 8.2° (68% CI: [0.0, 21.8])
Concentration κ = 1.22 (68% CI: [0.80, 1.86])
Information purity p_ID = 0.550 (68% CI: [0.444, 0.651])
Self-organization/dissipation ε/γ = 1.223
```

The best ωₚ₉ from the von Mises fit (8.2°) is equivalent to 356° (the fixed‑window result) because the mixture is symmetric under adding 180°.

---

## 6. Interpretation

- **κ = 1.22**: Moderate concentration – the ETNO ω distribution is clearly non‑uniform but not extremely peaked.
- **p_ID = 0.55**: The system is moderately self‑organized (ε slightly larger than γ). This places the outer Solar System between the highly ordered interstellar object 1I/‘Oumuamua (p≈0.83) and the thermally dominated comet 2I/Borisov (p≈0.09).
- **ε/γ = 1.22**: Self‑organization is only about 22% stronger than dissipation, implying a “gentle shepherd” rather than a rigid structure.

These parameters provide a quantitative description of the hypothetical Planet Nine as an information field eigenmode, and they make specific, testable predictions for future surveys (e.g., LSST).

---

## 7. Repository contents

- `planet_nine_inversion.py` – Main Python script (as described above)
- `README.md` – This file
- (Plots are generated on the fly and not stored in the repository)

---

## 8. License

This project is open source and available under the MIT License.

---

## 9. Citation

If you use this code or the results in your own work, please cite:

> Huang, K., & Liu, H. (2026). *Information Dynamics Inversion for Planet Nine: ETNO Argument of Perihelion Clustering as a Projection of a Self-Organized Mode*. GitHub repository. https://github.com/hkaiopen/InformationDynamics_PlanetNine

---

*For any questions or suggestions, please open an issue or contact the authors.*
