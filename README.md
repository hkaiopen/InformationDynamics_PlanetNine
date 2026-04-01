# Information Dynamics Inversion for Planet Nine

This repository contains a Python implementation to invert the **Information Dynamics** parameters for the hypothetical Planet Nine using the argument of perihelion (ω) distribution of extreme trans-Neptunian objects (ETNOs). The code downloads the latest MPCORB database, selects ETNOs, performs a fixed‑window clustering test, and fits a mixture of two von Mises distributions to extract the concentration parameter κ, information purity p_ID = κ/(1+κ), and the self‑organization to dissipation ratio ε/γ = p_ID/(1−p_ID).

---

## 1. Purpose

The aim is to provide a quantitative, data‑driven estimate of the dynamical state of the outer Solar System within the Information Dynamics framework. Instead of assuming a conventional massive planet, the observed N=2 symmetry in ETNO orbits is interpreted as a projection of an information field. The inverted parameters offer a direct, testable characterization of this field.

---

## 2. Method

1. **Data source**
   The Minor Planet Center Orbit Table (MPCORB.DAT), contains orbital data of hundreds of thousands of small celestial bodies (including orbital elements, absolute brightness, latest observation data, etc.) 
   1) I have uploaded `MPCORB.DAT` file accompanying the code.
   2) You may download `MPCORB.DAT` file manually from https://www.minorplanetcenter.net/data/.
   3) The script could automatically download the Minor Planet Center’s `MPCORB.DAT` file (gzipped) if not already present.

3. **ETNO selection**  
   Objects are selected with:
   - semi‑major axis \(a > 200\) AU  
   - perihelion distance \(q = a(1-e) > 30\) AU  
   - absolute magnitude \(H > 5\)

4. **Statistical clustering test**  
   A fixed window half‑width of 60° is used to measure the fraction of ETNOs whose ω falls within two opposite windows (centered on a test longitude ωₚ₉ and ωₚ₉+180°). The best‑fitting ωₚ₉ is found by scanning 0–360°. Significance is assessed with 10,000 Monte Carlo trials against a uniform distribution.

5. **Information Dynamics parameter inversion**  
   The ω distribution is modelled as a mixture of two von Mises distributions with equal concentration κ and means separated by 180°. Maximum likelihood estimation yields κ.  
   
   The information purity is defined as:  
   > p<sub>ID</sub> = κ / (1 + κ)  
   
   and the self‑organization to dissipation ratio as:  
   > ε / γ = p<sub>ID</sub> / (1 - p<sub>ID</sub>)  
   
   Bootstrap resampling (1000 trials) gives 68% confidence intervals.

6. **Visualization**  
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

This project is licensed under the **Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License** (CC BY-NC-SA 4.0).

This license allows you to:
*   **Share** — copy and redistribute the material in any medium or format.
*   **Adapt** — remix, transform, and build upon the material.

Under the following terms:
1.  **Attribution (BY)** — You must give **appropriate credit**, provide a link to the license, and **indicate if changes were made**. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use[1](@ref).
2.  **NonCommercial (NC)** — You may **not use the material for commercial purposes**. Commercial purposes include, but are not limited to:
    *   Selling products or services that incorporate this project.
    *   Using it in paid training or courses.
    *   Integrating it into commercial software.
    *   Any use aimed at monetary compensation or private financial gain.
3.  **ShareAlike (SA)** — If you remix, transform, or build upon the material, you **must distribute your contributions under the same license** as the original (CC BY-NC-SA 4.0)[1](@ref).

**For any commercial use, you must obtain prior written permission from the author.** Please contact the author to discuss licensing options.

To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/.

---

## Support This Independent Research
 
If you feel excited on this unifying framework and would like to support further development (simulations, data analysis, outreach), any contribution — big or small — is deeply appreciated.  We definitely call for **Research Funding and Collaboration** and are deeply committed to the highest standards of **Scientific Responsibility**.

### Donation Options
- **Buy Me A Coffee**  
  [Donate via Buymeacoffee](https://buymeacoffee.com/hkaiopen)
  
- **Stripe**  
  [Donate via Stripe](https://buy.stripe.com/eVqbJ12040RT1YLgg5gnK00)
  
- **PayPal**  
  [Donate via PayPal](https://paypal.me/kevinhuangkai)  
  (Supports one-time or recurring)

- **Bank transfer**
  Please DM me on X (@KevinHuangkai) for bank details.
  
All funds go toward researching more mysteries in the universe, and expanding the Information Dynamics model.  
Thank you for supporting open, frontier science!

---

## 9. Citation

If you use this code or the results in your own work, please cite the following works:

1.  **The associated research paper (preprint):**
    > Huang, K., & Liu, H. (2026). *Information Dynamics: Gravity as a Projection of the Information Field – Inversion and Test via the Orbital Clustering of Planet Nine*. Zenodo preprint. (_Updated on Mar 31, 2026_ **https://doi.org/10.5281/zenodo.19319558 This DOI represents all versions, and will always resolve to the latest version**).

2.  **The software implementation (this repository):**
    > Huang, K. (2026). *hkaiopen/InformationDynamics_PlanetNine: v1.0* (Version v1.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.19208434

In the acknowledgements section of any publication, we also encourage you to acknowledge the use of data from the Minor Planet Center, as per their request:
> This research has made use of data and/or services provided by the International Astronomical Union's Minor Planet Center.

---

*For any questions or suggestions, please open an issue or contact the authors.*
