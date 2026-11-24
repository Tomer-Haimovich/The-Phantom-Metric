# The Topological Saturation Model ($\Lambda\text{TSM}$ — Part II, Version 8)

**The Topological Saturation Model ($\Lambda\text{TSM}$): Clean Covariant DBI Cosmology and Full `Cobaya` MCMC Audit (Planck 2018 + DESI BAO + Pantheon+)**

* **Author:** Tomer Haimovich ([ORCID: 0009-0002-1536-1397](https://orcid.org/0009-0002-1536-1397))
* **Version & Date:** Version 8 — November 24, 2026
* **Contact:** tomerhaimovich@gmail.com
* **Part II Zenodo Archive (this package):** <https://doi.org/10.5281/zenodo.24463760>
* **Companion Part I (Galactic Lensing, Kinematics & LVDB Audit):** <https://doi.org/10.5281/zenodo.24463759>
* **GitHub Repository:** <https://github.com/Tomer-Haimovich/The-Phantom-Metric>

---

## 1. Overview & Clean Covariant Cosmological Sector

This reproducibility package accompanies **Version 8 (November 24, 2026)** of **Part II** of the Topological Saturation framework ($\Lambda\text{TSM}$). It upgrades the cosmological audit from the initial synthetic Planck-like Wishart consistency check to a **fully converged four-run real-likelihood `Cobaya` MCMC suite** ($R - 1 < 0.008$) spanning **Planck 2018** (`Plik_lite` $TT+TE+EE$ + low-$\ell$ `SROLL2` $EE$ + CMB lensing, $N = 623$), **DESI 2024 BAO**, and **Pantheon+ Type Ia supernovae** ($N = 2{,}225$), tightly synchronized to Part I (`10.5281/zenodo.24463759`):

* **Covariant Action:** Minimally coupled, shift-symmetric Horndeski $G_2(\phi, X)$ Dirac--Born--Infeld (DBI) scalar condensate ($M \lesssim 0.10\text{ meV}$ in natural units $\hbar = c = 1$, $V(\phi) = 0$, $\gamma_0 \gtrsim 10^5$) in its kinetically saturated dust-like limit ($w_{\text{eff}} = \gamma^{-1} \to 0$, $c_s^2 = \gamma^{-2} \to 0$, comoving sound horizon $r_{s,\phi}(0) \lesssim 16\text{ kpc}$, gravitational Jeans wavenumber $k_J(0) \gtrsim 15.4\text{ Mpc}^{-1}$, exact luminal tensor propagation $c_T = c$, exact separate covariant conservation $\nabla^\mu T_{\mu\nu}^{(b)} = 0$ and $\nabla^\mu T_{\mu\nu}^{(\phi)} = 0$).
* **Terminological Clarification:** $\omega_{\text{phantom}}$ ($\rho_\phi$) is **not** a ghost field ($w < -1$) and **not** a dark matter particle; it denotes the covariant gravitational energy density of the ghost-free saturated DBI condensate.
* **Zero-Freedom Structural Partition Postulate ($\Delta k = -1$ vs. $\Lambda\text{CDM}$):**
  $$\mathcal{R}_{4/3}^{\circ} \equiv \frac{\omega_{\text{phantom}}}{\omega_b} = \frac{16}{3}\left(1 + \frac{1}{16\pi^2}\right) \approx 5.36710706$$
  At $\omega_b = 0.02237$, this predicts $\omega_{\text{phantom}} = 0.120062$ ($+0.05\sigma$ from Planck 2018) and $\Omega_{m,0} = 0.31391$ at $H_0 = 67.36\text{ km s}^{-1}\text{ Mpc}^{-1}$.
* **Zero-Parameter Cross-Scale Bridge to Part I:**
  $$a_0(0) = \frac{c H_0}{\mathcal{R}_{4/3}^{\circ}} = 1.2194 \times 10^{-10}\text{ m s}^{-2}, \qquad a_0(z) = \frac{c H(z)}{\mathcal{R}_{4/3}^{\circ}} = a_0(0)\sqrt{\Omega_{m,0}(1+z)^3 + \Omega_\Lambda}$$

---

## 2. Package Contents (Part II — Version 8)

### Clean Covariant `CLASS v3.4.0` Source Patch (`source/`)
* `source/background.c`: Activates $\Omega_{0,\text{phantom}} = \mathcal{R}_{4/3}^{\circ}\,\Omega_{0,b}$ consistently across all four background routines whenever $\Omega_{0,\text{cdm}} < 10^{-5}$.
* `source/perturbations.c`: Strictly identical to official `CLASS v3.4.0` (zero modified differential equations; all legacy v1 non-covariant baryon-drag terms permanently removed).

### Full Four-Run `Cobaya` MCMC Suite
* **YAML Configurations:** `ltsm_planck_only.yaml`, `lcdm_planck_only.yaml`, `ltsm_planck_desi_sn.yaml`, `lcdm_planck_desi_sn.yaml`.
* **Setup & Execution:** `00_test_classy_consistency.py`, `01_setup_cobaya_ubuntu.sh`, `02_run_all_chains.sh`.
* **Analysis & Publication Figures:**
  * `03_analyze_cobaya_results.py` & `04_regenerate_clean_figures.py`: Compute `GetDist` marginalized posteriors ($30\%$ burn-in discarded), generate `cobaya_audit_report.txt`, and render 300 DPI publication figures (`Cobaya_Planck_DESI_SNe_Corner.png` and `Cobaya_Ratio_Falsifiability.png`).
  * `chains/` (`*.minimum.txt` and `*.progress` for all four runs: `ltsm_planck_only`, `lcdm_planck_only`, `ltsm_planck_desi_sn`, `lcdm_planck_desi_sn`): Deterministic `cobaya-run --minimize` best-fit solutions and final Gelman--Rubin convergence logs.

### Independent Noiseless Baseline & BBN Verification Scripts
* `audit_phantom_full.py` + `cmb_likelihood_utils.py`: Generates `phantom_full_diagnostic.png` ($\Delta\chi^2_{\text{Wishart}} = 0.1446$, $R^2_{TT} = 1.000000$, $R^2_{TE} = R^2_{EE} = 0.999999$).
* `phase1_calibrate_and_bbn.py`: Generates `phantom_bbn_expansion_proof.png` ($0.00000\%$ BBN expansion deviation at $z = 10^9$; max $0.0217\%$ across all $z$).
* `phase2_mcmc_engine.py`: Independent 5-parameter noiseless Wishart ensemble check ($N = 4{,}480$, $\hat{R} \le 1.047$; outputs: `mcmc_audit_report.txt`, `mcmc_planck_target.npz`, `phantom_mcmc_corner.png`).

### Manuscript Source & Compiled PDF (Companion Record Files)
* `part_2_manuscript_overleaf_latex_v8.tex` (and compiled PDF, hosted alongside this archive on the same Zenodo record): Complete Overleaf-ready (`pdfLaTeX`) two-column LaTeX source referencing Figures 1–4 (`phantom_full_diagnostic.png`, `phantom_bbn_expansion_proof.png`, `Cobaya_Planck_DESI_SNe_Corner.png`, and `Cobaya_Ratio_Falsifiability.png`).

---

## 3. Summary of Verified Version 8 Results

### A. Baseline `CLASS v3.4.0` Solver Verification ($\Sigma m_\nu = 0$)
* Undiluted $2\times 2$ Wishart $\Delta\chi^2_W = \mathbf{0.1446}$ (vs. $5{,}513.03$ for the removed legacy v1 drag block); $R^2_{TT} = 1.000000$, $R^2_{TE} = R^2_{EE} = 0.999999$.
* Recombination sound horizon $r_s(z_*) = 144.5320\text{ Mpc}$ ($0.0110\%$ difference); BBN expansion deviation $\Delta H/H(z = 10^9) = \mathbf{0.00000\%}$ (max $0.0217\%$ across all $z$).

### B. Full Four-Run Real-Likelihood `Cobaya` MCMC & Minimizer Audit ($\Sigma m_\nu = 0.06\text{ eV}$)

| Parameter / Metric | Planck 2018 ($\Lambda$TSM, $k=5$) | Planck 2018 ($\Lambda$CDM, $k=6$) | Planck+DESI+SN ($\Lambda$TSM, $k=5$) | Planck+DESI+SN ($\Lambda$CDM, $k=6$) |
| :--- | :---: | :---: | :---: | :---: |
| $\omega_b \equiv \Omega_b h^2$ | $0.02236 \pm 0.00008$ | $0.02237 \pm 0.00015$ | $0.02230 \pm 0.00008$ | $0.02246 \pm 0.00013$ |
| $\omega_{\text{phantom}}$ / $\omega_{\text{cdm}}$ | $0.12001 \pm 0.00045$ | $0.11995 \pm 0.00120$ | $0.11970 \pm 0.00043$ | $0.11866 \pm 0.00086$ |
| Ratio $\omega_c / \omega_b$ | **$5.3671$ (locked)** | $5.3628 \pm 0.0784$ ($+0.05\sigma$) | **$5.3671$ (locked)** | $5.2836 \pm 0.0572$ ($+1.46\sigma$) |
| $H_0$ [$\text{km s}^{-1}\text{ Mpc}^{-1}$] | **$67.32 \pm 0.12$** | $67.35 \pm 0.54$ | **$67.41 \pm 0.12$** | $67.93 \pm 0.39$ |
| $\Omega_{m,0}$ | **$0.3156 \pm 0.0021$** | $0.3153 \pm 0.0074$ | **$0.3140 \pm 0.0020$** | $0.3072 \pm 0.0051$ |
| $\ln(10^{10} A_s)$ | $3.0515 \pm 0.0115$ | $3.0520 \pm 0.0121$ | $3.0516 \pm 0.0117$ | $3.0554 \pm 0.0126$ |
| $n_s$ | $0.9635 \pm 0.0032$ | $0.9638 \pm 0.0041$ | $0.9645 \pm 0.0032$ | $0.9670 \pm 0.0036$ |
| $\tau_{\text{reio}}$ | $0.0581 \pm 0.0054$ | $0.0585 \pm 0.0061$ | $0.0587 \pm 0.0056$ | $0.0612 \pm 0.0061$ |
| $\sigma_8$ | $0.8135 \pm 0.0049$ | $0.8135 \pm 0.0053$ | $0.8131 \pm 0.0050$ | $0.8112 \pm 0.0054$ |
| $S_8$ | $0.8343 \pm 0.0061$ | $0.8339 \pm 0.0129$ | $0.8318 \pm 0.0061$ | $0.8209 \pm 0.0098$ |
| $r_d$ [Mpc] | $147.12 \pm 0.21$ | $147.12 \pm 0.26$ | $147.26 \pm 0.20$ | $147.37 \pm 0.23$ |
| $a_0 = cH_0/\mathcal{R}_{4/3}^{\circ}$ [$10^{-10}\text{ m s}^{-2}$] | $1.2186 \pm 0.0022$ | $1.2192 \pm 0.0097$ | **$1.2202 \pm 0.0021$** | $1.2298 \pm 0.0070$ |
| Accepted Steps / $R - 1$ | $27{,}559$ / $0.0074$ | $28{,}211$ / $0.0063$ | $16{,}147$ / $0.0066$ | $28{,}360$ / $0.0079$ |
| $\chi^2_{\min}$ (`--minimize`) | $981.43$ | $981.45$ | $2405.94$ | $2403.84$ |
| $\Delta\chi^2_{\min}$ / $\Delta\text{AIC}$ / $\Delta\text{BIC}$ | **$\approx 0.00$ / $-2.00$ / $-6.43$** ($N=623$) | — | **$+2.10$ / $+0.10$ / $-5.61$** ($N=2{,}225$) | — |

*Likelihood breakdown at joint global minimum (`--minimize`):* Planck CMB $+0.26$ ($982.26$ vs. $982.00$), Pantheon+ SNe Ia $-0.52$ ($1403.90$ vs. $1404.43$, favoring $\Lambda\text{TSM}$), DESI 2024 BAO $+2.36$ ($19.78$ vs. $17.42$). Raw MCMC chain-minimum cross-check (`cobaya_audit_report.txt`): Planck alone $-0.33$ / $-2.33$ / $-6.77$; joint $+2.75$ / $+0.74$ / $-4.96$.

### C. Synced Part I Galactic Benchmarks (`10.5281/zenodo.24463759`)
* **SLACS ($N=110$):** Jaffe $R^2_{1:1} = 0.7996$ static / $0.8103$ dynamic (SLACS IX); $0.8068$ / $0.8090$ (Combined); HST noise floor $0.8109 \pm 0.0388$; canonical $k = 1.0$ sits on the flat shared optimum ($k \approx 0.95\text{--}1.00$).
* **Euclid DR1 ECLIPSE ($N=1{,}095$, exact $\chi(\infty)$):** Dynamic $a_0(z_d)$ cuts unphysical failures by $-50.9\%$ ($9.2\%$ vs. $18.8\%$) in $z_d \in [0.70, 0.95)$, $-50.0\%$ ($8.1\%$ vs. $16.2\%$) in $[0.95, 1.50)$, and $-23.1\%$ overall ($29.8\%$ vs. $38.7\%$).
* **SPARC ($N=175$):** Locked $R^2_{1:1} = 0.9324$ (clean $0.9331$, bias $-0.0284\text{ dex}$, BTFR $0.8819$ / $0.8832$ direct-$D$); relaxed $0.9759$.
* **LVDB ($N=158$):** Classical dSph RMS $0.7246 \to 0.1634\text{ dex}$ ($0.1600\text{ dex}$ external, $0.130\text{ dex}$ isolated field dwarfs); flattening and golden-ratio auxiliaries falsified.

---

## 4. Citation

```bibtex
@misc{haimovich2026tsm,
  author       = {Haimovich, Tomer},
  title        = {The Topological Saturation Model ($\Lambda$TSM): Clean Covariant DBI Cosmology and Full Cobaya MCMC Audit (Planck 2018 + DESI BAO + Pantheon+) (Version 8)},
  year         = {2026},
  month        = nov,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.24463760},
  url          = {https://doi.org/10.5281/zenodo.24463760}
}
```