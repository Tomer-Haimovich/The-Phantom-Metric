# The Topological Saturation Model (ΛTSM) & Redshift-Dependent RAR

[![Zenodo Part I](https://img.shields.io/badge/Zenodo_Part_I_(V9_Latest)-10.5281%2Fzenodo.23293328-blue.svg)](https://doi.org/10.5281/zenodo.23293328)
[![Zenodo Part II](https://img.shields.io/badge/Zenodo_Part_II_(V9_Latest)-10.5281%2Fzenodo.23077417-blue.svg)](https://doi.org/10.5281/zenodo.23077417)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0002--1536--1397-A6CE39.svg)](https://orcid.org/0009-0002-1536-1397)
[![Boltzmann Solver](https://img.shields.io/badge/CLASS-v3.4.0_(C_%2F_Cython)-success.svg)](https://github.com/lesgourg/class_public)
[![MCMC](https://img.shields.io/badge/Cobaya-4--Run_MCMC-informational.svg)](https://cobaya.readthedocs.io/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-yellow.svg)](https://www.python.org/)
[![Dark Matter Parameters](https://img.shields.io/badge/Dark_Matter_Free_Parameters-0-critical.svg)](#)

**Author:** Tomer Haimovich ([ORCID: 0009-0002-1536-1397](https://orcid.org/0009-0002-1536-1397)) | **Contact:** tomerhaimovich@gmail.com  
**Current Release:** Version 9 — November 25, 2026 *(All Zenodo links below use permanent Concept DOIs that always resolve to the latest release)*

* **Part I — Galactic Lensing, Disk Kinematics & 3D Dwarf Spheroidals (V9):**  
  *Testing a Cosmologically Anchored, Redshift-Dependent Radial Acceleration Scale Across Strong Lensing, Disk Kinematics, and 3D Dwarf Spheroidals*  
  PDF: [`Part_I_Topological_Saturation_of_Spatial_Metrics_V9.pdf`](./Part_I_Topological_Saturation_of_Spatial_Metrics_V9.pdf) | Package: [`Part_I_Topological_Saturation_of_Spatial_Metrics_Package/`](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/README.md) | Concept DOI: [10.5281/zenodo.23293328](https://doi.org/10.5281/zenodo.23293328)
* **Part II — Clean Covariant DBI Cosmology & Full `Cobaya` MCMC Audit (V9):**  
  *The Topological Saturation Model (ΛTSM): Clean Covariant DBI Cosmology and Full Cobaya MCMC Audit (Planck 2018 + DESI BAO + Pantheon+)*  
  PDF: [`Part_II_The_Topological_Saturation_Model_V9.pdf`](./Part_II_The_Topological_Saturation_Model_V9.pdf) | Package: [`Part_II_The_Topological_Saturation_Model_Package/`](./Part_II_The_Topological_Saturation_Model_Package/README.md) | Concept DOI: [10.5281/zenodo.23077417](https://doi.org/10.5281/zenodo.23077417)

---

## Executive Summary

Standard six-parameter $\Lambda\text{CDM}$ cosmology introduces Cold Dark Matter (CDM) as an unconstrained free parameter ($\omega_{\text{cdm}} = 0.12000 \pm 0.00120$ in Planck 2018) alongside hundreds of per-galaxy halo parameters across galactic dynamics and strong-lensing catalogs.

This repository hosts the Version 9 companion manuscripts, clean covariant C modifications for **`CLASS v3.4.0`**, a converged four-run real-likelihood **`Cobaya` MCMC** suite (**Planck 2018 + DESI 2024 BAO + Pantheon+**), standalone Python verification pipelines, and an interactive Excel workbench for the **Topological Saturation Model ($\Lambda\text{TSM}$)**. Formulated as a nested five-parameter ($k = 5$, $\Delta k = -1$) sub-model within the minimally coupled, shift-symmetric Horndeski $G_2(X)$ **Dirac–Born–Infeld (DBI)** class, the ghost-free, kinetically saturated condensate ($M \lesssim 0.10\text{ meV}$, $V = 0$, $\gamma(a) = \sqrt{1 + (a_*/a)^6} \gg 1$) evolves in its dust-like limit ($w_{\text{eff}} = \gamma^{-1} \to 0$, $c_s^2 = \gamma^{-2} \to 0$, $r_{s,\phi}(0) \lesssim 16\text{ kpc}$, $k_J(0) \gtrsim 15.4\text{ Mpc}^{-1}$) with exact separate covariant conservation ($\nabla^\mu T_{\mu\nu} = 0$) and luminal tensor propagation ($c_T = c$). Eliminating $\omega_{\text{cdm}}$ via a single structural phase-transition postulate, $\mathcal{R}_{4/3}^{\circ} \equiv \omega_{\text{phantom}}/\omega_b = \frac{16}{3}(1 + \frac{1}{16\pi^2}) \approx 5.36710706$, and anchoring the galactic Radial Acceleration Relation (RAR) scale to $a_0(z) = cH(z)/\mathcal{R}_{4/3}^{\circ}$ links precision CMB/BAO/SNe cosmology directly to **five independent galactic archives** (**SLACS**, **Euclid DR1 ECLIPSE**, **MUSE-DARK III**, **SPARC**, and **LVDB**) with **0% particle Cold Dark Matter and zero per-galaxy tuning**.

> **Terminological Note:** $\omega_{\text{phantom}}$ ($\rho_\phi$) denotes **neither** a phantom ghost field ($w < -1$) **nor** a non-baryonic dark matter particle; it represents the covariant gravitational energy density of the ghost-free, shift-symmetric DBI condensate in its causal dust-like limit ($0 < c_s^2 < 1$).

---

## 1. Shared Covariant Baseline & Cross-Scale Lock

In natural units ($\hbar = c = 1$), the minimally coupled, shift-symmetric Horndeski $G_2(X)$ Einstein–DBI cosmological action ($V(\phi) = 0$, timelike branch $X \equiv -\frac{1}{2}g^{\mu\nu}\nabla_\mu\phi\nabla_\nu\phi > 0$) is:

$$S = \int d^4x \sqrt{-g} \left[ \frac{R - 2\Lambda}{16\pi G} - M^4\sqrt{1 - \frac{2X}{M^4}} + M^4 + \mathcal{L}_m(\psi_m, g_{\mu\nu}) \right]$$

Every cross-scale relation between Part II (cosmology) and Part I (galaxies and strong lensing) is anchored by the single redshift-invariant structural postulate $\mathcal{R}_{4/3}^{\circ}$ (treated transparently as a unifying phenomenological ansatz rather than a perturbative QFT derivation):

| Physical Quantity | Analytical Relation | Value (Baseline / Joint Posterior) | Physical Role & Synchronization |
| :--- | :--- | :--- | :--- |
| **Dimensional Partition Ratio** | $\mathcal{R}_{4/3}^{\circ} \equiv \frac{\omega_{\text{phantom}}}{\omega_b} = \frac{16}{3}\left(1 + \frac{1}{16\pi^2}\right)$ | $5.36710706$ | Eliminates $\omega_{\text{cdm}}$ ($\Delta k = -1$ vs. $\Lambda\text{CDM}$) |
| **Condensate Density** | $\omega_{\text{phantom}} = \mathcal{R}_{4/3}^{\circ}\,\omega_b$ | $0.120062$ (at $\omega_b = 0.02237$) | Matches Planck $\omega_{\text{cdm}} = 0.12000 \pm 0.00120$ ($+0.05\sigma$) |
| **Cosmological Baseline** | $\Omega_{m,0} = (\omega_b + \omega_{\text{phantom}} + \omega_\nu)/h^2$ | $0.31391$ ($\Sigma m_\nu = 0$) / $0.31533$ ($\Sigma m_\nu = 0.06\text{ eV}$) | Shared $H_0 = 67.36\text{ km s}^{-1}\text{ Mpc}^{-1}$; Part II joint posterior $H_0 = 67.41 \pm 0.12$, $\Omega_{m,0} = 0.3140 \pm 0.0020$ |
| **Saturated DBI Limit** | $M \lesssim 0.10\text{ meV},\; \gamma(a) = \sqrt{1 + (a_*/a)^6} \gg 1$ | $w_{\text{eff}} = \gamma^{-1} \to 0,\; c_s^2 = \gamma^{-2} \to 0$ | $r_{s,\phi}(0) \lesssim 16\text{ kpc}$, $k_J(0) \gtrsim 15.4\text{ Mpc}^{-1}$, $\nabla^\mu T_{\mu\nu} = 0$, $c_T = c$ |
| **Present-Day RAR Scale** | $a_0(0) = \frac{cH_0}{\mathcal{R}_{4/3}^{\circ}}$ | $1.2194 \times 10^{-10}\text{ m s}^{-2}$ | Joint `Cobaya` posterior: $(1.2202 \pm 0.0021) \times 10^{-10}\text{ m s}^{-2}$ ($+0.07\%$) |
| **Redshift-Dependent RAR Scale** | $a_0(z) = \frac{cH(z)}{\mathcal{R}_{4/3}^{\circ}} = a_0(0)\sqrt{\Omega_{m,0}(1+z)^3 + \Omega_\Lambda}$ | Dynamic with $z$ | Tested across $0 \le z \le 1.50$ (SLACS, Euclid DR1, MUSE-DARK III, SPARC, LVDB) |
| **RAR Interpolation Function** | $g_{\text{eff}} = g_{\text{bar}}\,\nu(y),\; \nu(y) = [1 - \exp(-\sqrt{y})]^{-1}$ | $y \equiv g_{\text{bar}}/a_0(z)$, $k = 1.0$ | Shared optimum ($k \approx 0.95\text{--}1.00$); golden-ratio & flattening gates excluded |

---

## 2. Part II: Clean Covariant DBI Cosmology & Full `Cobaya` MCMC Audit (`Part_II_The_Topological_Saturation_Model_Package`)

Inside [`Part_II_The_Topological_Saturation_Model_Package/`](./Part_II_The_Topological_Saturation_Model_Package/README.md), `source/background.c` activates $\Omega_{0,\text{phantom}} = \mathcal{R}_{4/3}^{\circ}\,\Omega_{0,b}$ consistently across all four `CLASS v3.4.0` background routines whenever $\Omega_{0,\text{cdm}} < 10^{-5}$, while `source/perturbations.c` remains strictly identical to official `CLASS v3.4.0` (zero non-covariant baryon-drag operators; exact $\nabla^\mu T_{\mu\nu} = 0$).

### A. `CLASS v3.4.0` Baseline Verification ($\Sigma m_\nu = 0$, $0\%$ Particle CDM)
* **Undiluted $2\times 2$ Wishart Discrepancy ($\ell \in [50, 2500]$):** $\Delta\chi^2_W = \mathbf{0.1446}$ (vs. $5{,}513.03$ for the excluded non-covariant drag control test), $R^2_{TT} = \mathbf{1.000000}$, $R^2_{TE} = R^2_{EE} = \mathbf{0.999999}$, max $|\Delta D_\ell^{TT}| = -1.65\text{ }\mu\text{K}^2$ ($-0.034\%$ at $\ell = 277$), binned Knox RMS $= 0.0260\,\sigma/\text{bin}$.
* **Recombination & BBN Preservation:** $r_s(z_*) = 144.5320\text{ Mpc}$ ($0.0110\%$ difference), $\sigma_8 = 0.82318$ (vs. $0.82291$), and BBN expansion deviation $\Delta H/H(z = 10^9) = \mathbf{0.00000\%}$ (max $0.0217\%$ across $z \in [0, 10^9]$).

### B. Converged Four-Run Real-Likelihood `Cobaya` MCMC & MAP Audit ($\Sigma m_\nu = 0.06\text{ eV}$, $R - 1 < 0.008$)

Evaluated against official **Planck 2018** (`Plik_lite` $TT+TE+EE$ + low-$\ell$ `SROLL2` $EE$ + CMB lensing, $N = 623$), **DESI 2024 BAO**, and **Pantheon+ SNe Ia** ($N = 2{,}225$), with deterministic MAP minimization (`cobaya-run --minimize`):

| Parameter / Model-Selection Metric | Planck 2018 ($\Lambda$TSM, $k=5$) | Planck 2018 ($\Lambda$CDM, $k=6$) | Planck+DESI+SN ($\Lambda$TSM, $k=5$) | Planck+DESI+SN ($\Lambda$CDM, $k=6$) |
| :--- | :---: | :---: | :---: | :---: |
| $\omega_b \equiv \Omega_b h^2$ | $0.02236 \pm 0.00008$ | $0.02237 \pm 0.00015$ | $0.02230 \pm 0.00008$ | $0.02246 \pm 0.00013$ |
| $\omega_{\text{phantom}}$ (0-P) / $\omega_{\text{cdm}}$ (Fit) | $0.12001 \pm 0.00045$ | $0.11995 \pm 0.00120$ | $0.11970 \pm 0.00043$ | $0.11866 \pm 0.00086$ |
| Ratio $\omega_c / \omega_b$ | **$5.3671$ (locked)** | $5.3628 \pm 0.0784$ (**$+0.05\sigma$**) | **$5.3671$ (locked)** | $5.2836 \pm 0.0572$ (**$+1.46\sigma$**) |
| $H_0$ [$\text{km s}^{-1}\text{ Mpc}^{-1}$] | **$67.32 \pm 0.12$** | $67.35 \pm 0.54$ | **$67.41 \pm 0.12$** | $67.93 \pm 0.39$ |
| $\Omega_{m,0}$ (incl. $\Sigma m_\nu = 0.06\text{ eV}$) | **$0.3156 \pm 0.0021$** | $0.3153 \pm 0.0074$ | **$0.3140 \pm 0.0020$** | $0.3072 \pm 0.0051$ |
| $\sigma_8$ / $S_8$ | $0.8135 \pm 0.0049$ / $0.8343 \pm 0.0061$ | $0.8135 \pm 0.0053$ / $0.8339 \pm 0.0129$ | $0.8131 \pm 0.0050$ / $0.8318 \pm 0.0061$ | $0.8112 \pm 0.0054$ / $0.8209 \pm 0.0098$ |
| Sound Horizon $r_d$ [Mpc] | $147.12 \pm 0.21$ | $147.12 \pm 0.26$ | $147.26 \pm 0.20$ | $147.37 \pm 0.23$ |
| Derived $a_0(0)$ [$10^{-10}\text{ m s}^{-2}$] | $1.2186 \pm 0.0022$ | $1.2192 \pm 0.0097$ | **$1.2202 \pm 0.0021$** | $1.2298 \pm 0.0070$ |
| Accepted Steps / $R - 1$ | $27{,}559$ / $0.0074$ | $28{,}211$ / $0.0063$ | $16{,}147$ / $0.0066$ | $28{,}360$ / $0.0079$ |
| $\Delta\chi^2_{\text{MAP}}$ (`--minimize`) | **$\approx 0.00$** *(chain $-0.33$)* | — | **$+2.10$** *(chain $+2.75$)* | — |
| $\Delta\text{AIC}$ / $\Delta\text{BIC}$ ($\Delta k = -1$) | **$-2.00$ / $-6.43$** ($N=623$) | — | **$[+0.10, +0.74]$ / $[-5.61, -4.96]$** ($N=2{,}225$) | — |

*Joint MAP Likelihood Breakdown:* Planck CMB $+0.26$, **Pantheon+ SNe Ia $-0.52$ (favoring $\Lambda\text{TSM}$)**, and DESI 2024 BAO $+2.36$.

### Part II Publication Figures (`./Part_II_The_Topological_Saturation_Model_Package/`)

| Fig. 1: Full CMB ($TT, TE, EE$) & $P(k)$ Audit (`CLASS v3.4.0`, $\Delta\chi^2_W = 0.1446$) | Fig. 2: Primordial BBN Expansion Proof ($z = 0 \to 10^9$, $\Delta H/H = 0.00000\%$) |
| :---: | :---: |
| ![Zero-Parameter LambdaTSM Diagnostic](./Part_II_The_Topological_Saturation_Model_Package/phantom_full_diagnostic.png) | ![BBN Expansion Proof](./Part_II_The_Topological_Saturation_Model_Package/phantom_bbn_expansion_proof.png) |

| Fig. 3: Converged `Cobaya` MCMC Posterior Contours (Planck + DESI + Pantheon+) | Fig. 4: Direct Observational Test of the $\mathcal{R}_{4/3}^{\circ}$ Partition Lock (`Cobaya`) |
| :---: | :---: |
| ![Cobaya MCMC Corner Plot](./Part_II_The_Topological_Saturation_Model_Package/Cobaya_Planck_DESI_SNe_Corner.png) | ![Cobaya Ratio Falsifiability](./Part_II_The_Topological_Saturation_Model_Package/Cobaya_Ratio_Falsifiability.png) |

---

## 3. Part I: Redshift-Dependent RAR Audit Across Five Independent Galactic Archives (`Part_I_Topological_Saturation_of_Spatial_Metrics_Package`)

Inside [`Part_I_Topological_Saturation_of_Spatial_Metrics_Package/`](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/README.md), the single cosmologically anchored scale $a_0(z) = cH(z)/\mathcal{R}_{4/3}^{\circ}$ ($a_0(0) = 1.2194 \times 10^{-10}\text{ m s}^{-2}$) coupled to $\nu(y) = [1 - \exp(-\sqrt{y})]^{-1}$ is audited across five independent archives ($0 \le z \le 1.50$) with zero dark matter and zero per-galaxy tuning:

| # | Archive | Sample & Methodology | Verified Version 9 Empirical Result (Zero Free Parameters) |
| :---: | :--- | :--- | :--- |
| **1** | **SLACS** (SLACS IX + S4TM) | $N = 110$ Grade-A early-type strong lenses; exact LOS geodesic integration, Jaffe ($\gamma = 2.0$) profile, $d_L^2$-rescaled masses | Dynamic $a_0(z_d)$ achieves **$R^2_{1:1} = 0.8103$** on SLACS IX ($N = 70$, matching the $0.08\text{ dex}$ HST MC noise ceiling $0.8109 \pm 0.0388$) and **$0.8090$** across $N = 110$ (static $a_0(0)$: $0.7996$ / $0.8068$, sum parity $99.2\%$ / $98.4\%$). |
| **2** | **Euclid DR1 ECLIPSE** | $N = 1{,}095$ Golden Grade-A lenses ($z_d \in [0.10, 1.50]$), exact $\chi(\infty)$ horizon bound + paired forward MC ($N_{\text{MC}} = 500$/scenario, S1–S3 up to $0.20\text{ dex}$ mass & $5\%$ photo-$z$ scatter) | Dynamic $a_0(z_d)$ reduces geometrically impossible failures vs. static MOND by **$-23.1\%$ overall** ($29.8\%$ vs. $38.7\%$) and **$-50\%$ at $z_d \ge 0.70$** ($9.2\%$ vs. $18.8\%$; $8.1\%$ vs. $16.2\%$); outperforms static MOND in **$100.0\%$ ($500/500$)** of MC realizations in every bin while isolating a low-$z$ group-scale transition ($\theta_E \ge 1.3''$). |
| **3** | **MUSE-DARK III** | $N = 79$ intermediate-$z$ star-forming disks ($0.33 < z < 1.44$; Ciocan et al. 2026) | Zero-parameter $a_0(1.00) = \mathbf{2.180 \times 10^{-10}\text{ m s}^{-2}}$ matches self-consistent MOND/RAR modeling ($2.19^{+0.12}_{-0.10}$) within **$0.4\%$** and linearized slope ($1.212$ vs. $1.20 \pm 0.10$) within **$1.0\%$**; fits 4-bin DC14 evolution with **$\chi^2 = 4.37$, $p = 0.359$** ($0.10\text{ dex}$ shared prior; vs. $p < 10^{-6}$ for static MOND). |
| **4** | **SPARC** | $N = 175$ late-type spiral galaxies ($3{,}391$ points; locked $\Upsilon_d = 0.50, \Upsilon_b = 0.70$, signed $V_{\text{gas}}$) | Unscaled $a_0(0)$ achieves **$R^2_{1:1} = 0.9331$** clean ($0.9324$ all; strict $Q = 1$ RAR RMS $0.1558\text{ dex}$, bias $-0.0090\text{ dex}$), direct-distance BTFR **$R^2_{1:1} = 0.8832$** ($0.0547\text{ dex}$, outperforming $cH_0/(2\pi)$), and relaxed SPS $0.9763$; **golden-ratio scaling disfavored**. |
| **5** | **LVDB** | $N = 158$ 3D dwarf spheroidals & star clusters (Wolf 2010 estimator) | Lowers classical dSph RMS from $0.7246\text{ dex}$ (Newtonian) to **$0.1634\text{ dex}$** (**$0.1600\text{ dex}$** on $N = 24$ strictly external $0\%$ MW dSphs; **$0.130\text{ dex}$** isolated field dwarfs), **ruling out morphological flattening gates** ($\rho = -0.005, p = 0.98$; optimal $\eta_c = 1.01$). |

### Part I Publication Figures (`./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/`)

| Fig. 1: SLACS 110 Strong Lenses ($R^2_{1:1} = 0.7996\text{--}0.8103$) | Fig. 2: Euclid DR1 ECLIPSE Forward MC Audit ($N = 1{,}095$, $100\%$ Realizations) |
| :---: | :---: |
| ![SLACS 110 Lensing](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/SLACS_110_Figure_v5.png) | ![Euclid DR1 ECLIPSE MC Audit](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/Euclid_ECLIPSE_MC_Uncertainty_Figure.png) |

| Fig. 3: MUSE-DARK III Disk Kinematics ($0.33 < z < 1.44$, $p = 0.359$) | Fig. 4: SPARC 175 Disks, BTFR & SPS Audit ($R^2_{1:1} = 0.9331$ / $0.8832$) |
| :---: | :---: |
| ![MUSE-DARK III Evolution](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/MUSE_DARK_III_a0z_Figure.png) | ![SPARC Unified Audit](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/SPARC_Unified_Figure_v5.png) |

| Fig. 5: LVDB 158 3D Spheroids & Flattening Audit ($0.1600\text{ dex}$ RMS) | Fig. 6: Joint $\nu_k(y)$ Stiffness Scan & HST $0.08\text{ dex}$ Noise Ceiling |
| :---: | :---: |
| ![LVDB Flattening Audit](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/LVDB_Flattening_Figure_v5.png) | ![SLACS SPARC Stiffness Noise Floor](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/SLACS_SPARC_Stiffness_NoiseFloor_v6.png) |

---

## 4. Theoretical Scope, Live Falsifiability & Roadmap to Part III

Both Version 9 manuscripts explicitly delineate the current boundaries and near-term falsification criteria of the framework:
* **Cosmological Red-Line Test:** Free $\Lambda\text{CDM}$ on Planck+DESI+Pantheon+ yields $\omega_{\text{cdm}}/\omega_b = 5.2836 \pm 0.0572$ ($+1.46\sigma$ from $\mathcal{R}_{4/3}^{\circ} = 5.3671$). Upcoming **DESI DR2/DR3**, **Simons Observatory**, and **CMB-S4** precision measurements will decisively confirm or falsify the rigid partition lock.
* **Relativistic Lensing Completion:** A purely conformally coupled scalar cannot enhance null-geodesic deflection ($\tilde{\Phi} + \tilde{\Psi} = 2\Phi_N$); realizing the weak-field zero-slip lensing metric ($\Psi = \Phi = \Phi_{\text{eff}}$) assumed in Part I while keeping $c_T = c$ requires disformal or unit-timelike vector extensions (e.g., TeVeS/AeST).
* **Galaxy Groups, Cluster Cores & Scale-Separated DBI Clustering:** Pure algebraic RAR leaves a $\sim 2\text{--}3\times$ mass deficit inside X-COP cluster cores (mirroring the Euclid DR1 low-$z$ transition at $\theta_E \ge 1.3''$). Avoiding galactic double-counting motivates testing nonlinear DBI gradient stabilization and an intermediate tension window ($M \sim 0.15\text{--}0.24\text{ meV}$) across X-COP, CLASH, and the Bullet Cluster in Part III.
* **Solar System EFE & Hubble Tension:** Outer Solar System tidal quadrupole ($Q_2$) bounds under the Galactic External Field Effect require dedicated 3D PDE screening solutions, and the $H_0 = 67.36 \to 73.0\text{ km s}^{-1}\text{ Mpc}^{-1}$ tension remains open outside the present 5-parameter background.

---

## 5. Reproducibility & Step-by-Step Execution Guide

All observational catalogs, C source patches, `Cobaya` YAML/chain files, standalone Python audit scripts, and the interactive Excel workbench (`Phantom_Clean_Master_SLACS_SPARC-v4.xlsx`) are bundled directly inside the two package directories.

### Prerequisites

```bash
pip install numpy scipy pandas matplotlib openpyxl cython cobaya getdist
