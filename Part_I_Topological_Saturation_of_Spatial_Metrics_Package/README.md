# The Phantom Metric (Part I — Version 8 / Unified Option A)

**Topological Saturation of Spatial Metrics: A Unified Single-Threshold Phenomenological Audit of Strong Lensing, Disk Kinematics, and 3D Dwarf Spheroidals**

* **Author:** Tomer Haimovich ([ORCID: 0009-0002-1536-1397](https://orcid.org/0009-0002-1536-1397))
* **Version & Date:** Version 8 — November 24, 2026
* **Contact:** tomerhaimovich@gmail.com
* **Part I Zenodo Archive (this package):** [https://doi.org/10.5281/zenodo.24463759](https://doi.org/10.5281/zenodo.24463759)
* **Companion Part II ($\Lambda\text{TSM}$ Cosmology & `Cobaya` MCMC Suite):** [https://doi.org/10.5281/zenodo.24463760](https://doi.org/10.5281/zenodo.24463760)
* **GitHub Repository:** [https://github.com/Tomer-Haimovich/The-Phantom-Metric](https://github.com/Tomer-Haimovich/The-Phantom-Metric)

---

## 1. Overview & Unified Cosmological Baseline

This reproducibility package accompanies **Version 8 (November 24, 2026)** of **Part I** of the Topological Saturation framework ($\Lambda\text{TSM}$). It links galactic-scale strong lensing, rotation curves, and 3D dwarf spheroidals directly to our companion covariant DBI cosmology (Part II, converged four-run `Cobaya` Planck 2018 + DESI 2024 BAO + Pantheon+ MCMC suite) under a **single shared cosmology**, a **single dimensionless structural constant**, and a **single universal galactic acceleration threshold**:

* **Shared Cosmological Baseline (Synced to Part II):**
  * $H_0 = 67.36\text{ km s}^{-1}\text{ Mpc}^{-1}$ ($h = 0.6736$; Part II `Cobaya` $\Lambda\text{TSM}$ posterior: $H_0 = 67.41 \pm 0.12\text{ km s}^{-1}\text{ Mpc}^{-1}$ on Planck+DESI+Pantheon+ and $67.32 \pm 0.12\text{ km s}^{-1}\text{ Mpc}^{-1}$ on Planck alone)
  * $\omega_b = 0.02237$
  * $\mathcal{R}_{4/3}^{\circ} \equiv \frac{16}{3}\left(1 + \frac{1}{16\pi^2}\right) \approx 5.36710706$ *(Spherically Normalized Dimensional Partition Ratio — redshift-invariant structural postulate from Part II)*
  * $\omega_{\text{phantom}} = \mathcal{R}_{4/3}^{\circ}\,\omega_b = 0.120062 \implies \Omega_{m,0} = (\omega_b + \omega_{\text{phantom}})/h^2 = 0.31391$ (Part II joint posterior: $\Omega_{m,0} = 0.3140 \pm 0.0020$, $\sigma_8 = 0.8131 \pm 0.0050$, $S_8 = 0.8318 \pm 0.0061$)
  * *Part II Model-Selection Summary (`cobaya-run --minimize`, $k = 5$ vs.\ $k = 6$ $\Lambda\text{CDM}$):* Planck 2018 alone (`Plik_lite` $TT+TE+EE$ + `SROLL2` $EE$ + Lensing, $N = 623$): $\Delta\chi^2_{\min} \approx 0.00$, $\Delta\text{AIC} \approx -2.00$, $\Delta\text{BIC} \approx -6.43$ (free $\Lambda\text{CDM}$ ratio $5.3628 \pm 0.0784$, $+0.05\sigma$); joint Planck + DESI 2024 BAO + Pantheon+ ($N = 2{,}225$): $\Delta\chi^2_{\min} = +2.10$, $\Delta\text{AIC} = +0.10$, $\Delta\text{BIC} = -5.61$ (free $\Lambda\text{CDM}$ ratio $5.2836 \pm 0.0572$, $+1.46\sigma$ driven by DESI BAO, with Pantheon+ preferring $\Lambda\text{TSM}$ by $-0.52$).
  * *Terminological note:* $\omega_{\text{phantom}}$ ($\rho_\phi$) is **not** a ghost field ($w < -1$) and **not** a dark matter particle; it denotes the effective gravitational energy density of the ghost-free, shift-symmetric DBI condensate in its dust-like limit ($w_{\text{eff}} \to 0, c_s^2 \to 0$, $M \lesssim 0.10\text{ meV}$, $r_{s,\phi}(0) \lesssim 16\text{ kpc}$).
* **Zero-Parameter Galactic Acceleration Scale & Cosmic-Time Evolution:**
  $$a_0 \equiv a_0(0) = \frac{c H_0}{\mathcal{R}_{4/3}^{\circ}} = 1.2194 \times 10^{-10}\text{ m s}^{-2}$$
  *(Part II joint posterior: $a_0 = (1.2202 \pm 0.0021) \times 10^{-10}\text{ m s}^{-2}$, a $+0.07\%$ shift leaving all galactic benchmarks invariant)*
  $$a_0(z) = \frac{c H(z)}{\mathcal{R}_{4/3}^{\circ}} = a_0(0)\sqrt{\Omega_{m,0}(1+z)^3 + \Omega_\Lambda}$$
* **Weak-Field RAR Enhancement Function:**
  $$g_{\text{eff}}(r) = g_{\text{bar}}(r)\,\nu(y), \qquad y \equiv \frac{g_{\text{bar}}(r)}{a_0(z)}, \qquad \nu(y) = \frac{1}{1 - \exp(-\sqrt{y})}$$

All four observational archives (`SLACS`, `Euclid DR1 ECLIPSE`, `SPARC`, and `LVDB`) are audited with **$0\%$ dark matter** and **zero per-galaxy tuning** in the primary locked baseline.

---

## 2. Package Contents (Part I — Version 8)

### Standalone Python Verification Suite (CLI + Google Colab Compatible)
| Script | Required Input Data (in same folder) | Outputs Generated |
| :--- | :--- | :--- |
| `SLACS_110_Code_v5.py` | `SLACS_ORIGN_DB.xlsx` | `SLACS_DB_v5.csv`, `SLACS_Phantom_Metric_Lensing_v5.csv`, `SLACS_110_Figure_v5.png` |
| `run_eclipse_golden_test_v2.py` | `eclipse_dr1_golden_1095.csv` *(or `--raw-dir` for raw Zenodo ECLIPSE)* | High-$z$ Table 2 (`Baryons_GR`, `MOND_const`, `LTSM_static`, `LTSM_Hz` with exact $\chi(\infty)$ horizon integral), Chabrier sensitivity, low-$z$ group-scale diagnostics |
| `SLACS_SPARC_Stiffness_NoiseFloor_v6.py` | `SLACS_DB_v5.csv` (or `SLACS_ORIGN_DB.xlsx`), `SPARC_Lelli2016c.txt`, `Rotmod_LTG.zip` | `SLACS_SPARC_Stiffness_NoiseFloor_v6.png` (Figure 4), `SLACS_SPARC_Stiffness_NoiseFloor_v6.csv`, dynamic $a_0(z_d)$ SLACS verification, Spearman diagnostics & HST Monte Carlo noise floor |
| `SPARC_Unified_Audit_v5.py` | `SPARC_Lelli2016c.txt`, `Rotmod_LTG.zip` | `SPARC_Unified_BTFR_v5.csv`, `SPARC_Unified_Global_v5.csv`, `SPARC_Unified_GalaxySummary_v5.csv`, `SPARC_Unified_Figure_v5.png` |
| `LVDB_Morphology_Audit_v5.py` | `comb_all.csv` | `LVDB_Unified_Audit_v5.csv`, `LVDB_Flattening_Figure_v5.png` |

### Interactive Excel Workbench (Companion Record File)
* `Phantom_Clean_Master_SLACS_SPARC-v4.xlsx`: Unified live-formula workbook covering `Control_Dashboard`, `SLACS_110` (with live $d_L^2$ mass rescaling, $E(z_d)$, $H(z_d)$, $a_0(z_d)$, and Dehnen/Hernquist/Jaffe LOS integration), `Euclid_ECLIPSE_1095` (exact $\chi(\infty)$ horizon bound, live redshift-bin summary & chart), `SPARC_175`, and `SPARC_NGC2403`. Editing cell `Control_Dashboard!C6` ($H_0$) or `C38` (`Dyn_a0_z` = `1` ON / `0` OFF) updates SLACS, SPARC, and Euclid DR1 live.

### Source Observational Catalogs Included
* `SLACS_ORIGN_DB.xlsx`: Original VizieR tables for SLACS IX ($N = 70$, Auger et al. 2009) and SLACS XIII S4TM ($N = 40$, Shu et al. 2017).
* `eclipse_dr1_golden_1095.csv`: Extracted Golden Subsample ($N = 1{,}095$ Grade A, $100\%$ consensus, non-disc, $\gamma < 0.10$, $1.5 \le n \le 8.0$, $z_d \in [0.15, 1.41]$) from the official Euclid DR1 ECLIPSE strong-lensing catalog (`v20261103`, $15{,}072$ Grade A/B lenses; DOI: [10.5281/zenodo.19654904](https://doi.org/10.5281/zenodo.19654904); Nersesian et al. 2026 `arXiv:2611.15877`, Lines et al. 2026 `arXiv:2611.15883`, Holloway et al. 2026 `arXiv:2611.15886`, Nightingale et al. 2026 `arXiv:2611.15888`).
* `SPARC_Lelli2016c.txt` & `Rotmod_LTG.zip`: Official SPARC master table and all $175$ raw `_rotmod.dat` radial rotation-curve files ($N = 3{,}391$ points; Lelli, McGaugh & Schombert 2016).
* `comb_all.csv`: Local Volume Database (LVDB; Pace 2025) containing nearby dwarf galaxies and globular/star clusters.

### Manuscript Source & Publication Figures
* `Topological_Saturation_of_Spatial_Metrics_v8.tex` (and compiled PDF): Complete Overleaf-ready two-column LaTeX source referencing Figures 1–4 (`SLACS_110_Figure_v5.png`, `SPARC_Unified_Figure_v5.png`, `LVDB_Flattening_Figure_v5.png`, and `SLACS_SPARC_Stiffness_NoiseFloor_v6.png`).

---

## 3. Quick Reproduction Instructions

Requires Python 3.9+ with `numpy`, `pandas`, `scipy`, `matplotlib`, and `openpyxl`:

```bash
pip install numpy pandas scipy matplotlib openpyxl

# 1. Strong Lensing Exact LOS Audit (SLACS N = 110)
python3 SLACS_110_Code_v5.py

# 2. High-Redshift Consistency Test on Euclid DR1 ECLIPSE (N = 1,095, exact chi(inf))
python3 run_eclipse_golden_test_v2.py

# 3. Joint Transition-Stiffness Scan, Dynamic a_0(z_d) SLACS & HST Noise-Floor MC (Figure 4)
python3 SLACS_SPARC_Stiffness_NoiseFloor_v6.py

# 4. Galactic Rotation Curves, BTFR & SPS Audit (SPARC N = 175)
python3 SPARC_Unified_Audit_v5.py

# 5. 3D Dwarf Spheroidals & Flattening Falsification Audit (LVDB N = 158)
python3 LVDB_Morphology_Audit_v5.py
```

*(Ensure all input catalogs reside in the same working directory as the scripts. Both `run_eclipse_golden_test_v2.py` and `SLACS_SPARC_Stiffness_NoiseFloor_v6.py` also execute natively inside Google Colab cells.)*

---

## 4. Summary of Verified Version 8 Empirical Results

### A. SLACS Strong Gravitational Lensing ($N = 110$, Exact LOS, $d_L^2$-Rescaled $M_*$)
* **Jaffe ($\gamma = 2.0$), SLACS IX ($N = 70$, Salpeter):**
  * Static $a_0(0)$: $R^2_{1:1} = \mathbf{0.7996}$ (vs.\ $0.6490$ pure baryonic GR), linear $R^2 = 0.8017$, mean $R_{\text{obs}}/R_{\text{calc}} = \mathbf{1.004 \pm 0.158}$, sum parity $\mathbf{99.2\%}$.
  * Dynamic $a_0(z_d)$: $R^2_{1:1} = \mathbf{0.8103}$, linear $R^2 = 0.8111$, mean $R_{\text{obs}}/R_{\text{calc}} = \mathbf{0.990}$, sum parity $\mathbf{99.1\%}$.
* **Jaffe ($\gamma = 2.0$), Combined ($N = 110$, Salp IX / Mid S4TM):**
  * Static $a_0(0)$: $R^2_{1:1} = \mathbf{0.8068}$ (vs.\ $0.7200$ GR), sum parity $\mathbf{98.4\%}$ (uniform Salpeter: $0.6572$).
  * Dynamic $a_0(z_d)$: $R^2_{1:1} = \mathbf{0.8090}$, linear $R^2 = 0.8153$, sum parity $\mathbf{96.9\%}$.
* **HST Noise Floor & Residual Diagnostics:** Monte Carlo simulations ($N_{\text{MC}} = 2{,}000$) confirm that $R^2_{1:1} = 0.7996$ sits at the HST photometric/mass resolution ceiling ($\mathbf{0.8109 \pm 0.0388}$ at $0.08\text{ dex}$; $\mathbf{0.7034 \pm 0.0659}$ at $0.10\text{ dex}$). Spearman tests show zero radial bias vs.\ $R_E/R_{\text{eff}}$ ($\rho = +0.063, p = 0.604$), mild positive trends with $z_d$ ($\rho = +0.351, p = 0.0029$ static; reduced to $\rho = +0.312, p = 0.0085$ under dynamic $a_0(z_d)$) and $\sigma_v$ ($\rho = +0.295, p = 0.013$), and median transition ratio $y_E = g_{\text{bar}}(R_E)/a_0(0) \approx 10.4$.

### B. Euclid DR1 ECLIPSE High-Redshift Consistency Test ($N = 1{,}095$, Native Salpeter)
Testing the exact asymptotic geometric bound $\beta_{\text{req}} \ge \beta_{\max}(z_d) = 1 - \chi(z_d)/\chi(\infty)$ ($I(\infty) = 3.244755$, $\Omega_{m,0} = 0.31391$) via Prugniel--Simien (1997) 3D deprojection and exact LOS geodesic integration:

| Redshift Bin ($z_d$) | $N$ | Med. $\beta_{\max}$ | Baryons GR (% Unphys) | MOND Const. (% Unphys) | $\Lambda$TSM Static (% Unphys) | $\Lambda$TSM Dynamic $a_0(z_d)$ (% Unphys) | Rel. Reduction vs. MOND |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| $0.10 \le z_d < 0.45$ | $332$ | $0.904$ | $85.8\%$ | $65.4\%$ | $64.8\%$ | **$60.2\%$** | $-7.8\%$ |
| $0.45 \le z_d < 0.70$ | $345$ | $0.849$ | $74.5\%$ | $38.3\%$ | $37.7\%$ | **$25.8\%$** | **$-32.6\%$** |
| **$0.70 \le z_d < 0.95$** | **$282$** | $0.799$ | $60.6\%$ | $18.8\%$ | $18.8\%$ | **$9.2\%$** | **$-50.9\%$** |
| **$0.95 \le z_d < 1.50$** | **$136$** | $0.756$ | $40.4\%$ | $16.2\%$ | $15.4\%$ | **$8.1\%$** | **$-50.0\%$** |
| **ALL ($0.10\text{--}1.50$)** | **$1{,}095$** | $0.840$ | $70.1\%$ | $38.7\%$ | $38.3\%$ | **$29.8\%$** | **$-23.1\%$** |

*Low-$z$ diagnostic:* Normal galaxy-scale lenses ($\theta_E < 0.9''$, $N=107$) sit near physical parity ($\text{med}\,\beta_{\text{req}} = 0.806$, $41.1\%$ unphysical), intermediate lenses ($0.9'' \le \theta_E < 1.3''$, $N=119$) reach $55.5\%$, whereas giant group-scale deflectors ($\theta_E \ge 1.3''$, $N=106$, $\text{med}\log M_* = 11.56$) drive the low-$z$ excess ($84.9\%$), marking the transition to intragroup hot gas and nonlinear DBI clustering ($\rho_\phi$, Part III).

### C. SPARC Disk Rotation Curves & BTFR ($N = 175$, $3{,}391$ Points)
* **All $175$ Galaxies (Locked SPS $\Upsilon_d=0.50, \Upsilon_b=0.70$):** Global $R^2_{1:1} = \mathbf{0.9324}$ (vs.\ $0.5540$ Newtonian); local $z \approx 0.004 \implies E(z) \approx 1.002$ ($+0.2\%$ shift, invariant under $a_0(z)$).
* **Clean ($Q<3, i\ge 30^\circ$, $N=153$):** Global $R^2_{1:1} = \mathbf{0.9331}$, per-galaxy median $\mathbf{+0.5193}$ (vs.\ $-3.4838$ Newtonian), RAR RMS $= \mathbf{0.1797\text{ dex}}$, bias $-0.0284\text{ dex}$; strict $Q=1$ ($N=93$): RMS $= \mathbf{0.1558\text{ dex}}$, bias $-0.0090\text{ dex}$.
* **Joint Transition Stiffness ($\nu_k(y) = [1-\exp(-y^{k/2})]^{-1/k}$):** Canonical $k = 1.0$ sits on the flat shared optimum ($k \approx 0.95\text{--}1.00$: $R^2_{1:1} = 0.7996$ vs.\ $0.8020$ at $k = 0.95$ on SLACS IX; exact peak $+0.5193$ on SPARC Clean), jointly falsifying $k = 2/3$ ($0.6214$ / $+0.1269$) and closely matching Simple MOND ($0.7983$ / $+0.5101$) with zero added free parameters.
* **BTFR ($V_{\text{flat}}$):** Clean ($N=123$) $R^2_{1:1} = \mathbf{0.8819}$ ($0.0671\text{ dex}$); strictly $H_0$-independent direct-$D$ ($f_D \ge 2, N=57$) $R^2_{1:1} = \mathbf{0.8832}$ ($0.0547\text{ dex}$, vs.\ $0.8433\text{--}0.8653$ for $cH_0/(2\pi)$).
* **Relaxed SPS ($\Upsilon_d \in [0.10, 1.00]$):** $R^2_{1:1} = \mathbf{0.9759}$ all / $\mathbf{0.9763}$ clean (median $\Upsilon_d = 0.428$, RAR RMS $= 0.1330\text{ dex}$).

### D. LVDB 3D Spheroids ($N = 158$) & Falsified Hypotheses
* Universal $a_0$ reduces classical dSph $\sigma_{\text{los}}$ RMS from $0.7246\text{ dex}$ (Newtonian) to $\mathbf{0.1634\text{ dex}}$ ($0.1600\text{ dex}$ on $N=24$ strictly external $0\%$ Milky Way dSphs; $0.130\text{ dex}$ on $N=14$ isolated field dwarfs; $0.235\text{ dex}$ on $N=64$ globular clusters at median $g_{\text{bar}} \approx 5.0\,a_0$).
* Morphological flattening shutoff gates ($\sqrt{1-\eta^2}$ and $\Theta(\eta_c-\eta)$), golden-ratio scaling ($a_\varphi = a_0\varphi$), and radial concentration rescaling are decisively falsified and permanently excluded.

---

## 5. Changelog (Version 8, November 24, 2026)

1. Added explicit cosmic-time evolution law $a_0(z) = cH(z)/\mathcal{R}_{4/3}^{\circ}$, verified its improvement on SLACS ($R^2_{1:1} = 0.8103$ IX / $0.8090$ Combined), and executed the first high-redshift test on Euclid DR1 ECLIPSE ($N = 1{,}095$ Golden Subsample with exact $\chi(\infty)$ horizon integration + `run_eclipse_golden_test_v2.py`).
2. Added SLACS IX residual diagnostics, HST Monte Carlo noise-floor audit, and joint SLACS+SPARC transition-stiffness scan (`SLACS_SPARC_Stiffness_NoiseFloor_v6.py` + Figure 4).
3. Locked final converged four-run `Cobaya` MCMC and deterministic `--minimize` results (`Planck 2018 Plik_lite + SROLL2 EE + Lensing`, `DESI 2024 BAO`, `Pantheon+`) into Section 1, Appendix A.1, and Appendix E.2, synced to Part II (DOI `10.5281/zenodo.24463760`).
4. Added explicit terminological clarification of $\omega_{\text{phantom}}$ / $\rho_\phi$ (ghost-free DBI condensate, not $w < -1$), unified unit/sign conventions in Appendices A–C ($M \lesssim 0.10\text{ meV}$, exact SI $T_{\mu\nu}^{(\phi)}$ and AQUAL normalization), aligned all bibliography credits and `Acknowledgments`, and linked X-COP cluster cores and Euclid group-scale lenses to Part III.
5. Released unified Excel workbook `Phantom_Clean_Master_SLACS_SPARC-v4.xlsx` with live $d_L^2$ mass rescaling, $H(z_d)$, $a_0(z_d)$, exact $\chi(\infty)$ horizon, and `Euclid_ECLIPSE_1095` sheet.

---

## 6. Citation

```bibtex
@misc{haimovich2026a,
  author       = {Haimovich, Tomer},
  title        = {Topological Saturation of Spatial Metrics: A Unified Single-Threshold Phenomenological Audit of Strong Lensing, Disk Kinematics, and 3D Dwarf Spheroidals (Version 8)},
  year         = {2026},
  month        = nov,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.24463759},
  url          = {https://doi.org/10.5281/zenodo.24463759}
}
```