# Testing a Cosmologically Anchored, Redshift-Dependent Radial Acceleration Scale (Part I — Version 9)

**Testing a Cosmologically Anchored, Redshift-Dependent Radial Acceleration Scale Across Strong Lensing, Disk Kinematics, and 3D Dwarf Spheroidals**

* **Author:** Tomer Haimovich ([ORCID: 0009-0002-1536-1397](https://orcid.org/0009-0002-1536-1397))
* **Version & Date:** Version 9 — November 25, 2026
* **Contact:** tomerhaimovich@gmail.com
* **Part I Zenodo Archive (Concept DOI — always resolves to latest version):** [https://doi.org/10.5281/zenodo.23293328](https://doi.org/10.5281/zenodo.23293328)
* **Companion Part II ($\Lambda\text{TSM}$ Cosmology & `Cobaya` MCMC Suite, Concept DOI):** [https://doi.org/10.5281/zenodo.23077417](https://doi.org/10.5281/zenodo.23077417)
* **GitHub Repository:** [https://github.com/Tomer-Haimovich/The-Phantom-Metric](https://github.com/Tomer-Haimovich/The-Phantom-Metric)

---

## 1. Overview & Cosmological Baseline

This reproducibility package accompanies **Version 9 (November 25, 2026)** of **Part I** of the Topological Saturation framework ($\Lambda\text{TSM}$). It tests a single, cosmologically anchored, redshift-dependent Radial Acceleration Relation (RAR) scale linked directly to our companion covariant Dirac--Born--Infeld (DBI) cosmology (Part II, converged four-run `Cobaya` Planck 2018 + DESI 2024 BAO + Pantheon+ MCMC suite) across **five independent observational archives** ($0 \le z \le 1.50$) with **$0\%$ dark matter** and **zero per-galaxy tuning** in the primary locked baseline:

* **Shared Flat Cosmological Baseline (Synchronized with Part II):**
  * $H_0 = 67.36\text{ km s}^{-1}\text{ Mpc}^{-1}$ ($h = 0.6736$; Part II `Cobaya` $\Lambda\text{TSM}$ posterior: $H_0 = 67.41 \pm 0.12\text{ km s}^{-1}\text{ Mpc}^{-1}$ on Planck+DESI+Pantheon+ and $67.32 \pm 0.12\text{ km s}^{-1}\text{ Mpc}^{-1}$ on Planck alone)
  * $\omega_b = 0.02237$
  * $\mathcal{R}_{4/3}^{\circ} \equiv \frac{16}{3}\left(1 + \frac{1}{16\pi^2}\right) \approx 5.36710706$ *(Spherically Normalized Dimensional Partition Ratio — redshift-invariant structural postulate from Part II)*
  * $\omega_{\text{phantom}} = \mathcal{R}_{4/3}^{\circ}\,\omega_b = 0.120062 \implies \Omega_{m,0} = (\omega_b + \omega_{\text{phantom}})/h^2 = 0.31391$ in the massless-neutrino galactic baseline ($0.31533$ including $\Sigma m_\nu = 0.06\text{ eV}$; Part II joint posterior: $\Omega_{m,0} = 0.3140 \pm 0.0020$, $\sigma_8 = 0.8131 \pm 0.0050$, $S_8 = 0.8318 \pm 0.0061$)
  * *Part II Model-Selection Summary (`cobaya-run --minimize`, $k = 5$ vs. $k = 6$ $\Lambda\text{CDM}$):* Planck 2018 alone (`Plik_lite` $TT+TE+EE$ + `SROLL2` $EE$ + lensing, $N = 623$): $\Delta\chi^2_{\text{MAP}} \approx 0.00$, $\Delta\text{AIC} \approx -2.00$, $\Delta\text{BIC} \approx -6.43$ (free $\Lambda\text{CDM}$ ratio $5.3628 \pm 0.0784$, $+0.05\sigma$); joint Planck + DESI 2024 BAO + Pantheon+ ($N = 2{,}225$): $\Delta\chi^2_{\text{MAP}} = +2.10$ (chain-minimum $+2.75$), $\Delta\text{AIC} \in [+0.10, +0.74]$, $\Delta\text{BIC} \in [-5.61, -4.96]$ (free $\Lambda\text{CDM}$ ratio $5.2836 \pm 0.0572$, $+1.46\sigma$ driven by DESI BAO, with Pantheon+ preferring $\Lambda\text{TSM}$ by $-0.52$).
  * *Terminological note:* $\omega_{\text{phantom}}$ ($\rho_\phi$) is **not** a ghost field ($w < -1$) and **not** a dark matter particle; it denotes the covariant gravitational energy density of the ghost-free, shift-symmetric DBI condensate in its dust-like limit ($w_{\text{eff}} \to 0, c_s^2 \to 0$, $M \lesssim 0.10\text{ meV}$, $r_{s,\phi}(0) \lesssim 16\text{ kpc}$).
* **Zero-Parameter Galactic Acceleration Scale & Cosmic-Time Evolution:**
  $$a_0 \equiv a_0(0) = \frac{c H_0}{\mathcal{R}_{4/3}^{\circ}} = 1.2194 \times 10^{-10}\text{ m s}^{-2}$$
  *(Part II joint posterior: $a_0(0) = (1.2202 \pm 0.0021) \times 10^{-10}\text{ m s}^{-2}$, a $+0.07\%$ shift leaving all galactic benchmarks invariant)*
  $$a_0(z) = \frac{c H(z)}{\mathcal{R}_{4/3}^{\circ}} = a_0(0)\sqrt{\Omega_{m,0}(1+z)^3 + \Omega_\Lambda}$$
* **Weak-Field RAR Enhancement Function:**
  $$g_{\text{eff}}(r) = g_{\text{bar}}(r)\,\nu(y), \qquad y \equiv \frac{g_{\text{bar}}(r)}{a_0(z)}, \qquad \nu(y) = \frac{1}{1 - \exp(-\sqrt{y})}$$

---

## 2. Package Contents (Part I — Version 9)

### Standalone Python Verification Suite (CLI + Google Colab Compatible)

| Script | Required Input Data (in same folder) | Outputs Generated |
| :--- | :--- | :--- |
| `SLACS_110_Code_v5.py` | `SLACS_ORIGN_DB.xlsx` | `SLACS_DB_v5.csv`, `SLACS_Phantom_Metric_Lensing_v5.csv`, `SLACS_110_Figure_v5.png` (Fig. 1) |
| `run_eclipse_golden_test_v2.py` | `eclipse_dr1_golden_1095.csv` *(or `--raw-dir` for raw Zenodo ECLIPSE)* | High-$z$ Table 2 (`Baryons_GR`, `MOND_const`, `LTSM_static`, `LTSM_Hz` with exact $\chi(\infty)$ horizon integral), Chabrier sensitivity, low-$z$ group-scale diagnostics |
| `run_eclipse_mc_uncertainty_audit.py` | `eclipse_dr1_golden_1095.csv` | `Euclid_ECLIPSE_MC_Audit_Results.csv` (Table 3), `Euclid_ECLIPSE_MC_Uncertainty_Figure.png` (Fig. 2; $N_{\text{MC}} = 500$ realizations/scenario across S1–S3) |
| `MUSE_DARK_III_Audit.py` | Built-in published Ciocan et al. (2026) / Logue (2026) tables | `MUSE_DARK_III_Audit_Results.csv` (Table 4), `MUSE_DARK_III_a0z_Figure.png` (Fig. 3; App. E MOND/RAR + 4-bin DC14 covariance audit) |
| `SPARC_Unified_Audit_v5.py` | `SPARC_Lelli2016c.txt`, `Rotmod_LTG.zip` | `SPARC_Unified_BTFR_v5.csv`, `SPARC_Unified_Global_v5.csv`, `SPARC_Unified_GalaxySummary_v5.csv`, `SPARC_Unified_Figure_v5.png` (Fig. 4) |
| `LVDB_Morphology_Audit_v5.py` | `comb_all.csv` | `LVDB_Unified_Audit_v5.csv`, `LVDB_Flattening_Figure_v5.png` (Fig. 5) |
| `SLACS_SPARC_Stiffness_NoiseFloor_v6.py` | `SLACS_DB_v5.csv` (or `SLACS_ORIGN_DB.xlsx`), `SPARC_Lelli2016c.txt`, `Rotmod_LTG.zip` | `SLACS_SPARC_Stiffness_NoiseFloor_v6.csv` (Table 7), `SLACS_SPARC_Stiffness_NoiseFloor_v6.png` (Fig. 6), dynamic $a_0(z_d)$ SLACS audit, Spearman diagnostics & HST MC noise floor |

### Interactive Excel Workbench
* `Phantom_Clean_Master_SLACS_SPARC-v4.xlsx`: Unified live-formula workbook covering `Control_Dashboard`, `SLACS_110` (live $d_L^2$ mass rescaling, $E(z_d)$, $H(z_d)$, $a_0(z_d)$, and Dehnen/Hernquist/Jaffe LOS integration), `Euclid_ECLIPSE_1095` (exact $\chi(\infty)$ horizon bound, live redshift-bin summary & chart), `SPARC_175`, and `SPARC_NGC2403`. Editing `Control_Dashboard!C6` ($H_0$) or `C38` (`Dyn_a0_z` = `1` ON / `0` OFF) updates SLACS, SPARC, and Euclid DR1 live.

### Source Observational Catalogs Included
* `SLACS_ORIGN_DB.xlsx`: Original VizieR tables for SLACS IX ($N = 70$, Auger et al. 2009) and SLACS XIII S4TM ($N = 40$, Shu et al. 2017).
* `eclipse_dr1_golden_1095.csv`: Extracted Golden Subsample ($N = 1{,}095$ Grade A, $100\%$ consensus, non-disc, $\gamma < 0.10$, $1.5 \le n \le 8.0$, $z_d \in [0.15, 1.41]$, $19$ with DESI spec-$z$) from the Euclid DR1 ECLIPSE strong-lensing catalog (`v20261103`, $15{,}072$ Grade A/B lenses; DOI: [10.5281/zenodo.19654904](https://doi.org/10.5281/zenodo.19654904); Nersesian et al. 2026, Lines et al. 2026, Holloway et al. 2026, Nightingale et al. 2026).
* **MUSE-DARK III Reference Data:** Intermediate-redshift ($0.33 < z < 1.44$, $N = 79$ SFGs) Appendix E self-consistent MOND/RAR measurements and 4-bin DC14 medians/uncertainties from Ciocan et al. (2026, `arXiv:2604.22613`) and Logue (2026, MNRAS Table 2), embedded directly in `MUSE_DARK_III_Audit.py`.
* `SPARC_Lelli2016c.txt` & `Rotmod_LTG.zip`: Official SPARC master table and all $175$ raw `_rotmod.dat` radial rotation-curve files ($N = 3{,}391$ points; Lelli, McGaugh & Schombert 2016).
* `comb_all.csv`: Local Volume Database (LVDB; Pace 2025) containing nearby dwarf galaxies and globular/star clusters.

### Manuscript Source & Publication Figures (Companion Record Files)
* `EPJ_Plus_Part_I_Final.tex` (and compiled PDF): Complete Overleaf-ready (`pdfLaTeX`) two-column LaTeX source referencing Figures 1–6 (`SLACS_110_Figure_v5.png`, `Euclid_ECLIPSE_MC_Uncertainty_Figure.png`, `MUSE_DARK_III_a0z_Figure.png`, `SPARC_Unified_Figure_v5.png`, `LVDB_Flattening_Figure_v5.png`, and `SLACS_SPARC_Stiffness_NoiseFloor_v6.png`).

---

## 3. Quick Reproduction Instructions

Requires Python 3.9+ with `numpy`, `pandas`, `scipy`, `matplotlib`, and `openpyxl`:

```bash
pip install numpy pandas scipy matplotlib openpyxl

# 1. Strong Lensing Exact LOS Audit (SLACS N = 110) -> Fig. 1
python3 SLACS_110_Code_v5.py

# 2. High-Redshift Consistency Test on Euclid DR1 ECLIPSE (N = 1,095, exact chi(inf)) -> Table 2
python3 run_eclipse_golden_test_v2.py

# 3. Forward Monte Carlo Uncertainty Propagation on Euclid DR1 (N_MC = 500/scenario) -> Table 3 & Fig. 2
python3 run_eclipse_mc_uncertainty_audit.py

# 4. Intermediate-Redshift Disk Kinematics Audit (MUSE-DARK III, N = 79) -> Table 4 & Fig. 3
python3 MUSE_DARK_III_Audit.py

# 5. Galactic Rotation Curves, BTFR & SPS Audit (SPARC N = 175) -> Table 5 & Fig. 4
python3 SPARC_Unified_Audit_v5.py

# 6. 3D Dwarf Spheroidals & Morphological Gating Audit (LVDB N = 158) -> Table 6 & Fig. 5
python3 LVDB_Morphology_Audit_v5.py

# 7. Joint Transition-Stiffness Scan, Dynamic a_0(z_d) SLACS & HST Noise-Floor MC -> Table 7 & Fig. 6
python3 SLACS_SPARC_Stiffness_NoiseFloor_v6.py
```

*(Ensure all input catalogs reside in the same working directory as the scripts. `run_eclipse_golden_test_v2.py`, `run_eclipse_mc_uncertainty_audit.py`, `MUSE_DARK_III_Audit.py`, and `SLACS_SPARC_Stiffness_NoiseFloor_v6.py` also execute natively inside Google Colab cells.)*

---

## 4. Summary of Verified Version 9 Empirical Results

### A. SLACS Strong Gravitational Lensing ($N = 110$, Exact LOS, $d_L^2$-Rescaled $M_*$)
* **Jaffe ($\gamma = 2.0$), SLACS IX ($N = 70$, Salpeter):**
  * Static $a_0(0)$: $R^2_{1:1} = \mathbf{0.7996}$ (vs. $0.6490$ pure baryonic GR), linear $R^2 = 0.8017$, mean $R_{\text{obs}}/R_{\text{calc}} = \mathbf{1.004 \pm 0.158}$, sum parity $\mathbf{99.2\%}$.
  * Dynamic $a_0(z_d)$: $R^2_{1:1} = \mathbf{0.8103}$, linear $R^2 = 0.8111$, mean $R_{\text{obs}}/R_{\text{calc}} = \mathbf{0.990}$, sum parity $\mathbf{99.1\%}$.
* **Jaffe ($\gamma = 2.0$), Combined ($N = 110$, Salp IX / Mid S4TM):**
  * Static $a_0(0)$: $R^2_{1:1} = \mathbf{0.8068}$ (vs. $0.7200$ GR), sum parity $\mathbf{98.4\%}$ (uniform Salpeter: $0.6572$).
  * Dynamic $a_0(z_d)$: $R^2_{1:1} = \mathbf{0.8090}$, linear $R^2 = 0.8153$, sum parity $\mathbf{96.9\%}$.
* **HST Noise Floor & Residual Diagnostics:** Monte Carlo simulations confirm that $R^2_{1:1} = 0.7996\text{--}0.8103$ sits at the HST photometric/mass resolution ceiling ($\mathbf{0.8109 \pm 0.0388}$ at $0.08\text{ dex}$; $\mathbf{0.7034 \pm 0.0659}$ at $0.10\text{ dex}$). Spearman tests show zero radial drift vs. $R_E/R_{\text{eff}}$ ($\rho = +0.063, p = 0.604$), while dynamic $a_0(z_d)$ lowers the redshift residual trend from $\rho = +0.351$ ($p = 0.0029$) to $\rho = +0.312$ ($p = 0.0085$).

### B. Euclid DR1 ECLIPSE High-Redshift Test & Forward Monte Carlo Audit ($N = 1{,}095$, $N_{\text{MC}} = 500$/Scenario)
Testing the exact asymptotic geometric bound $\beta_{\text{req}} \ge \beta_{\max}(z_d) = 1 - \chi(z_d)/\chi(\infty)$ ($I(\infty) = 3.244755$, $\Omega_{m,0} = 0.31391$) under catalog baseline and three paired Monte Carlo uncertainty scenarios (`S1`: $0.15\text{ dex}$ mass, $4\%$ photo-$z$, $3\%\,\theta_E$; `S2`: `S1` + $M_* \propto d_L^2$ coupling; `S3`: $0.20\text{ dex}$ mass, $5\%$ photo-$z$, $d_L^2$ coupling, $+0.03$ low-$z$ photo-$z$ bias; $19$ DESI spec-$z$ held fixed):

| Redshift Bin ($z_d$) | $N$ | Base MOND $\to$ Dyn | S1 MOND $\to$ Dyn (%) | S1 Rel. Red. | S2 ($d_L^2$) Rel. Red. | S3 (Heavy) Rel. Red. | $P(\text{Dyn} < \text{MOND})$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| $0.10 \le z_d < 0.45$ | $332$ | $65.4\% \to 60.2\%$ | $60.5 \pm 1.9 \to 56.3 \pm 2.0$ | $-6.9 \pm 1.8\%$ | $-6.5 \pm 1.8\%$ | $-7.2 \pm 1.9\%$ | **$100.0\%$** |
| $0.45 \le z_d < 0.70$ | $345$ | $38.3\% \to 25.8\%$ | $39.9 \pm 2.0 \to 31.2 \pm 1.9$ | $-21.9 \pm 3.3\%$ | $-19.7 \pm 3.2\%$ | $-17.0 \pm 3.0\%$ | **$100.0\%$** |
| **$0.70 \le z_d < 0.95$** | **$282$** | $18.8\% \to 9.2\%$ | $26.6 \pm 2.2 \to 15.9 \pm 1.7$ | **$-40.3 \pm 5.5\%$** | **$-36.8 \pm 5.1\%$** | **$-31.1 \pm 4.6\%$** | **$100.0\%$** |
| **$0.95 \le z_d < 1.50$** | **$136$** | $16.2\% \to 8.1\%$ | $20.6 \pm 2.4 \to 11.3 \pm 1.8$ | **$-44.8 \pm 8.4\%$** | **$-43.1 \pm 8.9\%$** | **$-39.8 \pm 8.0\%$** | **$100.0\%$** |
| **ALL ($0.10\text{--}1.50$)** | **$1{,}095$** | $38.7\% \to 29.8\%$ | $40.3 \pm 1.1 \to 32.4 \pm 1.0$ | **$-19.7 \pm 1.8\%$** | **$-18.3 \pm 1.7\%$** | **$-17.6 \pm 1.7\%$** | **$100.0\%$** |

*Low-$z$ group-scale transition:* Persists cleanly under MC scatter (`S1`: $\theta_E < 0.9''$ at $41.5 \pm 3.6\%$, $0.9''\text{--}1.3''$ at $53.3 \pm 3.4\%$, and $\theta_E \ge 1.3''$ at $74.6 \pm 3.4\%$; baseline: $41.1\%$, $55.5\%$, $84.9\%$), isolating Brightest Group Galaxies where hot X-ray intragroup gas and nonlinear DBI clustering enter.

### C. MUSE-DARK III Intermediate-Redshift Disk Kinematics ($N = 79$, $0.33 < z < 1.44$)
* **Self-Consistent MOND/RAR Comparison (Ciocan et al. 2026, App. E):** Zero-parameter $a_0(1.00) = \mathbf{2.180 \times 10^{-10}\text{ m s}^{-2}}$ matches the observed $2.19^{+0.12}_{-0.10} \times 10^{-10}\text{ m s}^{-2}$ ($95\%$ CI) within **$\mathbf{0.4\%}$** ($<0.2\sigma$). The exact linear tangent of $a_0(0)E(z)$ over $[0.333, 1.438]$ yields slope $a_1 = \mathbf{1.212}$ and intercept $\mathbf{0.993}$, reproducing Eq. (E.5) ($a_1 = 1.20 \pm 0.10$, $a_0(0) = 1.03 \pm 0.05$, $95\%$ CI) within **$1.0\%$**.
* **4-Bin DC14 Goodness-of-Fit ($0$ Free Parameters):** Shape-only $\chi^2_{\text{shape}} = \mathbf{4.12}$ ($p = 0.249$, $3$ dof); under a correlated $0.10\text{ dex}$ shared calibration prior (half the $\sim 0.2\text{ dex}$ molecular-gas budget in App. C of Ciocan et al. 2026), dynamic $\Lambda\text{TSM}$ achieves $\chi^2 = \mathbf{4.37}$ ($\mathbf{p = 0.359}$, $4$ dof; $p = 0.382$ at $0.20\text{ dex}$), whereas static MOND is rejected at $\chi^2 = 38.23$ ($p < 10^{-6}$, $\Delta\chi^2 = -33.9$ shared / $-841$ independent-bin).

### D. SPARC Disk Rotation Curves & BTFR ($N = 175$, $3{,}391$ Points)
* **Locked SPS ($\Upsilon_d = 0.50, \Upsilon_b = 0.70$):** Global $R^2_{1:1} = \mathbf{0.9324}$ all / $\mathbf{0.9331}$ clean ($Q < 3, i \ge 30^\circ, N = 153$; per-galaxy median $\mathbf{+0.5193}$ vs. $-3.4838$ Newtonian, RAR RMS $= \mathbf{0.1797\text{ dex}}$, bias $-0.0284\text{ dex}$; strict $Q = 1, N = 93$: RMS $= \mathbf{0.1558\text{ dex}}$, bias $-0.0090\text{ dex}$).
* **Joint Transition Stiffness ($\nu_k(y) = [1 - \exp(-y^{k/2})]^{-1/k}$):** Canonical $k = 1.0$ sits on the flat shared optimum ($k \approx 0.95\text{--}1.00$) of SLACS IX and SPARC Clean, jointly disfavoring $k = 2/3$ ($0.6214$ / $+0.1269$).
* **BTFR ($V_{\text{flat}}$):** Clean ($N = 123$) $R^2_{1:1} = \mathbf{0.8819}$ ($0.0671\text{ dex}$); strictly $H_0$-independent direct-distance ($f_D \ge 2, N = 57$) $R^2_{1:1} = \mathbf{0.8832}$ ($0.0547\text{ dex}$, outperforming $cH_0/(2\pi)$ at $0.8433\text{--}0.8653$).
* **Relaxed SPS ($\Upsilon_d \in [0.10, 1.00]$):** $R^2_{1:1} = \mathbf{0.9759}$ all / $\mathbf{0.9763}$ clean (median $\Upsilon_d = 0.428$, RAR RMS $= 0.1330\text{ dex}$).

### E. LVDB 3D Spheroids ($N = 158$) & Excluded Auxiliary Scalings
* Universal unscaled $a_0(0)$ lowers classical dSph $\sigma_{\text{los}}$ RMS from $0.7246\text{ dex}$ (Newtonian) to $\mathbf{0.1634\text{ dex}}$ ($\mathbf{0.1600\text{ dex}}$ on $N = 24$ strictly external $0\%$ Milky Way dSphs; $0.130\text{ dex}$ on $N = 14$ isolated field dwarfs; $0.235\text{ dex}$ on $N = 64$ compact globular clusters at median $g_{\text{bar}} \approx 5.0\,a_0$).
* Morphological flattening shutoff gates ($\sqrt{1-\eta^2}$ and $\Theta(\eta_c-\eta)$; Spearman $\rho = -0.005, p = 0.98$, optimal $\eta_c = 1.01$), golden-ratio scaling ($a_\varphi = a_0\varphi$), and radial concentration rescaling are ruled out by the data.

---

## 5. Changelog (Version 9 — November 25, 2026)

1. **Added Independent Intermediate-Redshift Disk Kinematics Audit (MUSE-DARK III):** Added `MUSE_DARK_III_Audit.py`, `MUSE_DARK_III_Audit_Results.csv`, and Figure 3 (`MUSE_DARK_III_a0z_Figure.png`), demonstrating $0.4\%$ agreement with the self-consistent MOND/RAR inference at $z \sim 1$ ($2.180$ vs. $2.19^{+0.12}_{-0.10} \times 10^{-10}\text{ m s}^{-2}$) and $\chi^2 = 4.37$ ($p = 0.359$) across the 4-bin DC14 evolution.
2. **Added Forward-Modeling Monte Carlo Uncertainty Audit for Euclid DR1 ECLIPSE:** Added `run_eclipse_mc_uncertainty_audit.py`, `Euclid_ECLIPSE_MC_Audit_Results.csv`, and Figure 2 (`Euclid_ECLIPSE_MC_Uncertainty_Figure.png`), verifying $100.0\%$ ($500/500$) paired outperformance of dynamic $a_0(z_d)$ over static MOND across all four redshift bins and three stress scenarios (up to $0.20\text{ dex}$ stellar-mass scatter, $5\%$ photo-$z$ scatter, correlated $d_L^2$ coupling, and $+0.03$ low-$z$ photo-$z$ bias).
3. **Journal-Ready Manuscript Alignment & Concept DOI Synchronization:** Updated title, tightened abstract to $\le 250$ words, removed draft-history references, added transparent generative-AI disclosure in Acknowledgements, and switched all cross-references between Parts I and II to permanent Zenodo Concept DOIs (`10.5281/zenodo.23293328` and `10.5281/zenodo.23077417`).

---

## 6. Citation

```bibtex
@misc{haimovich2026a,
  author       = {Haimovich, Tomer},
  title        = {Testing a Cosmologically Anchored, Redshift-Dependent Radial Acceleration Scale Across Strong Lensing, Disk Kinematics, and 3D Dwarf Spheroidals},
  year         = {2026},
  month        = nov,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.23293328},
  url          = {https://doi.org/10.5281/zenodo.23293328}
}
```