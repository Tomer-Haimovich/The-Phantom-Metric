# The Phantom Metric &amp; Topological Saturation Model (&Lambda;TSM) &mdash; Version 8

[![Zenodo Part I](https://img.shields.io/badge/Zenodo_Part_I_(Galactic_Audit)-10.5281%2Fzenodo.24463759-blue.svg)](https://doi.org/10.5281/zenodo.24463759)
[![Zenodo Part II](https://img.shields.io/badge/Zenodo_Part_II_(DBI_Cosmology_%26_Cobaya)-10.5281%2Fzenodo.24463760-blue.svg)](https://doi.org/10.5281/zenodo.24463760)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0002--1536--1397-A6CE39.svg)](https://orcid.org/0009-0002-1536-1397)
[![Boltzmann Solver](https://img.shields.io/badge/CLASS-v3.4.0_(Clean_Covariant_C)-success.svg)](https://github.com/lesgourg/class_public)
[![MCMC Suite](https://img.shields.io/badge/Cobaya_MCMC-Planck%2BDESI%2BPantheon%2B-purple.svg)](https://cobaya.readthedocs.io/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-yellow.svg)](https://www.python.org/)
[![Dark Matter Free Parameters](https://img.shields.io/badge/Dark_Matter_Free_Parameters-0-critical.svg)](#)

**Author:** Tomer Haimovich ([ORCID: 0009-0002-1536-1397](https://orcid.org/0009-0002-1536-1397))  
**Version &amp; Release Date:** Version 8 &mdash; November 24, 2026  
**Preprints &amp; Permanent Archives:**
* **Part I (Galactic Strong Lensing, Disk Kinematics, Euclid DR1 ECLIPSE &amp; 3D Dwarf Spheroidals):** [`Part_I_Topological_Saturation_of_Spatial_Metrics_V8.pdf`](./Part_I_Topological_Saturation_of_Spatial_Metrics_V8.pdf) &mdash; [DOI: 10.5281/zenodo.24463759](https://doi.org/10.5281/zenodo.24463759)  
  *Reproducibility Directory:* [`./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/`](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/)
* **Part II (Clean Covariant DBI Cosmology, `CLASS v3.4.0` &amp; Full Four-Run `Cobaya` MCMC Audit):** [`Part_II_The_Topological_Saturation_Model_V8.pdf`](./Part_II_The_Topological_Saturation_Model_V8.pdf) &mdash; [DOI: 10.5281/zenodo.24463760](https://doi.org/10.5281/zenodo.24463760)  
  *Reproducibility Directory:* [`./Part_II_The_Topological_Saturation_Model_Package/`](./Part_II_The_Topological_Saturation_Model_Package/)

---

## Executive Summary

Standard six-parameter &Lambda;CDM cosmology introduces non-baryonic Cold Dark Matter (CDM) as an unconstrained free parameter (&omega;<sub>cdm</sub> &equiv; &Omega;<sub>cdm</sub><em>h</em><sup>2</sup> = 0.12000 &plusmn; 0.00120 in Planck 2018) alongside hundreds of per-galaxy dark halo parameters across galactic rotation curves and lensing systems.

This repository hosts the complete Version 8 manuscripts, clean covariant **`CLASS v3.4.0`** C source modifications, converged four-run real-likelihood **`Cobaya` MCMC** configuration and audit suite (**Planck 2018 + DESI 2024 BAO + Pantheon+**), standalone Python verification pipelines, raw observational archives, and interactive Excel workbench for the **Topological Saturation Model (&Lambda;TSM)**.

Under a **single shared cosmology** (<em>H</em><sub>0</sub> = 67.36 km s<sup>&minus;1</sup> Mpc<sup>&minus;1</sup>, &omega;<sub>b</sub> = 0.02237, &Omega;<sub>m,0</sub> = 0.31391 in the massless-neutrino baseline; 0.31533 including &Sigma;<em>m</em><sub>&nu;</sub> = 0.06 eV), a **single dimensionless structural constant** (<em>R</em><sup>&deg;</sup><sub>4/3</sub> &approx; 5.36710706), and a **single universal galactic acceleration threshold** (<em>a</em><sub>0</sub>(<em>z</em>) = <em>cH</em>(<em>z</em>)/<em>R</em><sup>&deg;</sup><sub>4/3</sub>), the framework unifies:
1. **Linear Early-Universe Cosmology (Part II):** A minimally coupled, shift-symmetric Horndeski <em>G</em><sub>2</sub>(&phi;, <em>X</em>) Dirac&ndash;Born&ndash;Infeld (DBI) scalar condensate in its kinetically saturated dust-like limit (<em>w</em><sub>eff</sub> &rarr; 0, <em>c</em><sub>s</sub><sup>2</sup> &rarr; 0, <em>c</em><sub>T</sub> = <em>c</em>, exact &nabla;<sup>&mu;</sup><em>T</em><sub>&mu;&nu;</sub> = 0, zero non-covariant baryon-drag terms) that eliminates &omega;<sub>cdm</sub> (&Delta;<em>k</em> = &minus;1), reproduces the matched CMB spectra (&Delta;&chi;<sup>2</sup><sub>Wishart</sub> = 0.1446, BBN &Delta;<em>H</em>/<em>H</em>(<em>z</em> = 10<sup>9</sup>) = 0.00000%), and achieves full observational parity with 6-parameter &Lambda;CDM across real Planck 2018, DESI 2024 BAO, and Pantheon+ likelihoods.
2. **Late-Universe Galactic Dynamics &amp; Strong Lensing (Part I):** Zero-parameter empirical audits with 0% dark matter across **SLACS** (<em>N</em> = 110), **Euclid DR1 ECLIPSE** (<em>N</em> = 1,095), **SPARC** (<em>N</em> = 175, 3,391 points), and **LVDB** (<em>N</em> = 158 3D spheroids &amp; globular clusters), while transparently documenting falsified auxiliary hypotheses and theoretical boundaries reserved for **Part III**.

> **Terminological Clarification:** Throughout this repository, &omega;<sub>phantom</sub> (&rho;<sub>&phi;</sub>) is **neither** a <em>w</em> &lt; &minus;1 ghost field **nor** a non-baryonic dark matter particle; it denotes the covariant gravitational energy density of the ghost-free, shift-symmetric saturated DBI condensate in its dust-like limit (<em>w</em><sub>eff</sub> = &gamma;<sup>&minus;1</sup> &rarr; 0, <em>c</em><sub>s</sub><sup>2</sup> = &gamma;<sup>&minus;2</sup> &rarr; 0).

---

## Repository Layout

<pre><code>.
├── Part_I_Topological_Saturation_of_Spatial_Metrics_V8.pdf       # Part I Manuscript (18 pp., v8)
├── Part_II_The_Topological_Saturation_Model_V8.pdf               # Part II Manuscript (15 pp., v8)
├── README.md                                                     # Root Overview (this file)
├── Part_I_Topological_Saturation_of_Spatial_Metrics_Package/     # Galactic &amp; Lensing Suite
│   ├── README.md
│   ├── Phantom_Clean_Master_SLACS_SPARC-v4.xlsx                  # Live-Formula Excel Workbench
│   ├── SLACS_110_Code_v5.py &amp; SLACS_ORIGN_DB.xlsx                # SLACS N=110 Exact LOS
│   ├── run_eclipse_golden_test_v2.py &amp; eclipse_dr1_golden_1095.csv
│   ├── SLACS_SPARC_Stiffness_NoiseFloor_v6.py                    # Stiffness &amp; HST Noise Floor
│   ├── SPARC_Unified_Audit_v5.py, SPARC_Lelli2016c.txt, Rotmod_LTG.zip
│   ├── LVDB_Morphology_Audit_v5.py &amp; comb_all.csv                # LVDB N=158 3D Spheroids
│   └── *_v5.csv, *_v6.csv, *_v5.png, *_v6.png                    # Verified Tables &amp; Figures
└── Part_II_The_Topological_Saturation_Model_Package/             # CLASS v3.4.0 &amp; Cobaya Suite
    ├── README.md
    ├── source/background.c &amp; source/perturbations.c              # Clean Covariant CLASS C Files
    ├── chains/*.minimum.txt &amp; chains/*.progress                  # 4-Run MAP &amp; R-1 Logs
    ├── ltsm_planck_only.yaml, lcdm_planck_only.yaml              # Planck 2018 YAML Configs
    ├── ltsm_planck_desi_sn.yaml, lcdm_planck_desi_sn.yaml        # Joint Suite YAML Configs
    ├── 00_test_classy_consistency.py ... 04_regenerate_clean_figures.py
    ├── audit_phantom_full.py, cmb_likelihood_utils.py            # Undiluted 2x2 Wishart Audit
    ├── phase1_calibrate_and_bbn.py &amp; phase2_mcmc_engine.py       # BBN Proof &amp; Ensemble Check
    └── cobaya_audit_report.txt, mcmc_audit_report.txt, *.png     # Reports &amp; 300 DPI Figures</code></pre>

---

## 1. Unified Covariant Architecture (Single Constant, Single Galactic Threshold)

The cosmological sector is governed by a minimally coupled, shift-symmetric Horndeski <em>G</em><sub>2</sub>(&phi;, <em>X</em>) DBI action (<em>M</em> = const &#8818; 0.10 meV in natural units &#8463; = <em>c</em> = 1, <em>V</em>(&phi;) = 0, signature (&minus;, +, +, +)):

```math
S = \int d^4x \sqrt{-g} \left[ \frac{R - 2\Lambda}{16\pi G} - M^4\sqrt{1 - \frac{2X}{M^4}} + M^4 + \mathcal{L}_m(\psi_m, g_{\mu\nu}) \right], \qquad X \equiv -\frac{1}{2}g^{\mu\nu}\nabla_\mu\phi\nabla_\nu\phi
```

With <em>G</em><sub>3</sub> = <em>G</em><sub>5</sub> = 0 and constant <em>G</em><sub>4</sub> = (16&pi;<em>G</em>)<sup>&minus;1</sup>, tensor propagation is strictly luminal (<em>c</em><sub>T</sub> = <em>c</em>, satisfying GW170817 identically) and minimal coupling enforces exact separate covariant conservation (&nabla;<sup>&mu;</sup><em>T</em><sup>(m)</sup><sub>&mu;&nu;</sub> = 0, &nabla;<sup>&mu;</sup><em>T</em><sup>(&phi;)</sup><sub>&mu;&nu;</sub> = 0). Exact shift-symmetric charge conservation fixes &gamma;(<em>a</em>) = &radic;[1 + (<em>a</em><sub>*</sub>/<em>a</em>)<sup>6</sup>], giving <em>w</em><sub>eff</sub> = &gamma;<sup>&minus;1</sup> &rarr; 0, <em>c</em><sub>s</sub><sup>2</sup> = &gamma;<sup>&minus;2</sup> &rarr; 0, comoving sound horizon <em>r</em><sub>s,&phi;</sub>(0) &#8818; 16 kpc, and Jeans wavenumber <em>k</em><sub>J</sub>(0) &#8819; 15.4 Mpc<sup>&minus;1</sup> for &gamma;<sub>0</sub> &#8819; 10<sup>5</sup>.

| Structural Quantity | Unified Formula | Value (Baseline) | Role &amp; Methodological Status |
| :--- | :--- | :--- | :--- |
| **Dimensional Partition Ratio** | <em>R</em><sup>&deg;</sup><sub>4/3</sub> &equiv; (16/3)(1 + 1/(16&pi;<sup>2</sup>)) | 5.36710706 | Structural phase-transition postulate (4D/3D partition + (4&pi;)<sup>&minus;2</sup> heat-kernel closure) |
| **Cosmological Condensate Density** | &omega;<sub>phantom</sub> = <em>R</em><sup>&deg;</sup><sub>4/3</sub> &middot; &omega;<sub>b</sub> | 0.120062 | Eliminates &omega;<sub>cdm</sub> (&Delta;<em>k</em> = &minus;1); +0.05&sigma; from Planck 2018 at &omega;<sub>b</sub> = 0.02237 |
| **Present-Day Galactic Threshold** | <em>a</em><sub>0</sub>(0) = <em>cH</em><sub>0</sub> / <em>R</em><sup>&deg;</sup><sub>4/3</sub> | 1.2194 &times; 10<sup>&minus;10</sup> m s<sup>&minus;2</sup> | Zero-parameter bridge at <em>H</em><sub>0</sub> = 67.36 (joint `Cobaya` posterior: 1.2202 &plusmn; 0.0021 &times; 10<sup>&minus;10</sup>) |
| **Cosmic-Time Evolution Law** | <em>a</em><sub>0</sub>(<em>z</em>) = <em>cH</em>(<em>z</em>) / <em>R</em><sup>&deg;</sup><sub>4/3</sub> = <em>a</em><sub>0</sub>(0)<em>E</em>(<em>z</em>) | <em>E</em>(<em>z</em>) = &radic;[&Omega;<sub>m,0</sub>(1 + <em>z</em>)<sup>3</sup> + &Omega;<sub>&Lambda;</sub>] | Redshift-running horizon threshold directly tested on SLACS and Euclid DR1 ECLIPSE |
| **Universal Weak-Field Transition** | &nu;(<em>y</em>) = [1 &minus; exp(&minus;&radic;<em>y</em>)]<sup>&minus;1</sup> | <em>y</em> &equiv; <em>g</em><sub>bar</sub>(<em>r</em>) / <em>a</em><sub>0</sub>(<em>z</em>) | Canonical <em>k</em> = 1.0 RAR curve (shared flat optimum <em>k</em> &approx; 0.95&ndash;1.00 of SLACS IX &amp; SPARC) |

**Empirical Falsifications &amp; Cleanups Documented in Version 8:**
* **Permanently Removed Legacy Perturbation Drag (`v1` Erratum):** Non-covariant baryon-drag terms (<em>K</em><sub>0</sub> = 2/3, &Gamma;<sub>0</sub> = 1/(8&pi;), nonlinear tanh gate in `tca_off`) violated &nabla;<sup>&mu;</sup><em>T</em><sup>(b)</sup><sub>&mu;&nu;</sub> = 0 and incurred an undiluted Wishart penalty of &Delta;&chi;<sup>2</sup><sub>W</sub> = 5,513.03. Removing them completely yields `source/perturbations.c` strictly identical to official `CLASS v3.4.0` and collapses &Delta;&chi;<sup>2</sup><sub>W</sub> to **0.1446**.
* **Falsified Galactic Auxiliary Hypotheses (Part I):** Golden-ratio disk rescaling (<em>a</em><sub>&phi;</sub> = <em>a</em><sub>0</sub>&phi;, rejected by SPARC), 2D/3D morphological flattening shutoff gates (&radic;(1 &minus; &eta;<sup>2</sup>) and &Theta;(&eta;<sub>c</sub> &minus; &eta;), rejected by LVDB classical dSphs with Spearman &rho; = &minus;0.005, <em>p</em> = 0.98, optimal &eta;<sub>c</sub> = 1.01), radial concentration rescaling, and circular baryon-fraction cluster tests are all decisively falsified and permanently excluded.

---

## 2. Part II: Clean Covariant Cosmology &amp; Full `Cobaya` MCMC Audit

### A. Real-Likelihood Four-Run `Cobaya` MCMC Suite (<em>R</em> &minus; 1 &lt; 0.008, &Sigma;<em>m</em><sub>&nu;</sub> = 0.06 eV)

All four runs use official foreground-marginalized **Planck 2018** (`Plik_lite` <em>TT</em>+<em>TE</em>+<em>EE</em> + low-&#8467; `SROLL2` <em>EE</em> + CMB lensing, <em>N</em> = 623; <em>N</em> = 650 for all 28 `SROLL2` multipoles shifts &Delta;BIC by &minus;0.04), joined with **DESI 2024 BAO** and **Pantheon+ SNe Ia** (<em>N</em> = 2,225). Primary model-selection rows report deterministic MAP optimization (`cobaya-run --minimize`, including the Gaussian <em>A</em><sub>Planck</sub> prior) alongside raw MCMC chain-minimum cross-checks:

| Parameter / Model-Selection Metric | Planck 2018 (&Lambda;TSM, <em>k</em> = 5) | Planck 2018 (&Lambda;CDM, <em>k</em> = 6) | Planck+DESI+SN (&Lambda;TSM, <em>k</em> = 5) | Planck+DESI+SN (&Lambda;CDM, <em>k</em> = 6) |
| :--- | :---: | :---: | :---: | :---: |
| &omega;<sub>b</sub> &equiv; &Omega;<sub>b</sub><em>h</em><sup>2</sup> | 0.02236 &plusmn; 0.00008 | 0.02237 &plusmn; 0.00015 | 0.02230 &plusmn; 0.00008 | 0.02246 &plusmn; 0.00013 |
| &omega;<sub>phantom</sub> (0-P) / &omega;<sub>cdm</sub> (Fit) | 0.12001 &plusmn; 0.00045 | 0.11995 &plusmn; 0.00120 | 0.11970 &plusmn; 0.00043 | 0.11866 &plusmn; 0.00086 |
| Ratio &omega;<sub>phantom</sub>/&omega;<sub>b</sub> or &omega;<sub>cdm</sub>/&omega;<sub>b</sub> | **5.3671 (locked)** | 5.3628 &plusmn; 0.0784 (+0.05&sigma;) | **5.3671 (locked)** | 5.2836 &plusmn; 0.0572 (+1.46&sigma;) |
| <em>H</em><sub>0</sub> [km s<sup>&minus;1</sup> Mpc<sup>&minus;1</sup>] | **67.32 &plusmn; 0.12** | 67.35 &plusmn; 0.54 | **67.41 &plusmn; 0.12** | 67.93 &plusmn; 0.39 |
| &Omega;<sub>m,0</sub> (incl. &Sigma;<em>m</em><sub>&nu;</sub> = 0.06 eV) | **0.3156 &plusmn; 0.0021** | 0.3153 &plusmn; 0.0074 | **0.3140 &plusmn; 0.0020** | 0.3072 &plusmn; 0.0051 |
| &sigma;<sub>8</sub> / <em>S</em><sub>8</sub> | 0.8135 &plusmn; 0.0049 / 0.8343 | 0.8135 &plusmn; 0.0053 / 0.8339 | 0.8131 &plusmn; 0.0050 / 0.8318 | 0.8112 &plusmn; 0.0054 / 0.8209 |
| Sound Horizon <em>r</em><sub>d</sub> [Mpc] | 147.12 &plusmn; 0.21 | 147.12 &plusmn; 0.26 | 147.26 &plusmn; 0.20 | 147.37 &plusmn; 0.23 |
| Derived <em>a</em><sub>0</sub> [10<sup>&minus;10</sup> m s<sup>&minus;2</sup>] | 1.2186 &plusmn; 0.0022 | 1.2192 &plusmn; 0.0097 | **1.2202 &plusmn; 0.0021** | 1.2298 &plusmn; 0.0070 |
| Accepted Steps / Gelman&ndash;Rubin <em>R</em> &minus; 1 | 27,559 / 0.0074 | 28,211 / 0.0063 | 16,147 / 0.0066 | 28,360 / 0.0079 |
| &chi;<sup>2</sup><sub>MAP</sub> (`--minimize`) / Chain Min | 981.43 / 980.97 | 981.45 / 981.30 | 2405.94 / 2405.63 | 2403.84 / 2402.88 |
| **&Delta;&chi;<sup>2</sup><sub>MAP</sub> / &Delta;AIC / &Delta;BIC** (&Delta;<em>k</em> = &minus;1) | **&approx; 0.00 / &minus;2.00 / &minus;6.43** | *(chain: &minus;0.33 / &minus;2.33 / &minus;6.77)* | **+2.10 / +0.10 / &minus;5.61** | *(chain: +2.75 / +0.74 / &minus;4.96)* |

*Joint MAP likelihood breakdown:* Planck CMB +0.26 (982.26 vs. 982.00), **Pantheon+ SNe Ia &minus;0.52 (1403.90 vs. 1404.43, favoring &Lambda;TSM)**, DESI 2024 BAO +2.36 (19.78 vs. 17.42).  
*Baseline `CLASS v3.4.0` Solver Verification (&Sigma;<em>m</em><sub>&nu;</sub> = 0):* Undiluted 2&times;2 Wishart &Delta;&chi;<sup>2</sup><sub>W</sub> = **0.1446** (<em>R</em><sup>2</sup><sub>TT</sub> = 1.000000, <em>R</em><sup>2</sup><sub>TE</sub> = <em>R</em><sup>2</sup><sub>EE</sub> = 0.999999, max &#124;&Delta;<em>D</em><sup>TT</sup><sub>&#8467;</sub>&#124; = &minus;1.65 &mu;K<sup>2</sup>), <em>r</em><sub>s</sub>(<em>z</em><sub>*</sub>) = 144.5320 Mpc (0.0110% diff), BBN &Delta;<em>H</em>/<em>H</em>(<em>z</em> = 10<sup>9</sup>) = **0.00000%** (max 0.0217%).

### B. Part II Cosmological Figures (`./Part_II_The_Topological_Saturation_Model_Package/`)

#### 1. Full Real-Likelihood `Cobaya` Marginalized Contours (Planck 2018 + DESI 2024 BAO + Pantheon+)
![Cobaya Corner Plot](./Part_II_The_Topological_Saturation_Model_Package/Cobaya_Planck_DESI_SNe_Corner.png)

#### 2. Direct Observational Test of the &Lambda;TSM Partition Lock (<em>R</em><sup>&deg;</sup><sub>4/3</sub> = 5.3671)
![Cobaya Ratio Falsifiability](./Part_II_The_Topological_Saturation_Model_Package/Cobaya_Ratio_Falsifiability.png)

#### 3. Clean `CLASS v3.4.0` 4-Panel CMB (<em>TT, TE, EE</em>) &amp; <em>P</em>(<em>k</em>) Audit + Primordial BBN Proof (<em>z</em> = 0 &rarr; 10<sup>9</sup>)
| 4-Panel CMB &amp; Linear <em>P</em>(<em>k</em>) Baseline (&Delta;&chi;<sup>2</sup><sub>W</sub> = 0.1446) | Early-Universe BBN Expansion Proof (&Delta;<em>H</em>/<em>H</em> = 0.00000%) |
| :---: | :---: |
| ![CLASS 4-Panel Diagnostic](./Part_II_The_Topological_Saturation_Model_Package/phantom_full_diagnostic.png) | ![BBN Proof](./Part_II_The_Topological_Saturation_Model_Package/phantom_bbn_expansion_proof.png) |

---

## 3. Part I: Unified Single-Threshold Galactic Audit (0% Dark Matter)

Locking the exact Part II baseline (<em>H</em><sub>0</sub> = 67.36, &Omega;<sub>m,0</sub> = 0.31391, <em>R</em><sup>&deg;</sup><sub>4/3</sub> = 5.36710706) into <em>a</em><sub>0</sub>(<em>z</em>) = <em>cH</em>(<em>z</em>)/<em>R</em><sup>&deg;</sup><sub>4/3</sub> and &nu;(<em>y</em>) = [1 &minus; exp(&minus;&radic;<em>y</em>)]<sup>&minus;1</sup> with zero per-galaxy tuning yields unified empirical parity across four independent archives:

### A. Summary of Verified Part I Empirical Benchmarks

1. **SLACS Strong Gravitational Lensing (<em>N</em> = 110, Exact LOS Geodesic Integration, <em>d</em><sub>L</sub><sup>2</sup>-Rescaled <em>M</em><sub>*</sub>):**
   * **SLACS IX (<em>N</em> = 70, Jaffe &gamma; = 2.0, Salpeter):** Static <em>a</em><sub>0</sub>(0) achieves **<em>R</em><sup>2</sup><sub>1:1</sub> = 0.7996** (vs. 0.6490 baryonic GR; linear <em>R</em><sup>2</sup> = 0.8017, mean <em>R</em><sub>obs</sub>/<em>R</em><sub>calc</sub> = 1.004 &plusmn; 0.158, sum parity **99.2%**); activating dynamic <em>a</em><sub>0</sub>(<em>z</em><sub>d</sub>) raises parity to **<em>R</em><sup>2</sup><sub>1:1</sub> = 0.8103** (99.1% parity).
   * **Combined SLACS IX + S4TM (<em>N</em> = 110, Salp IX / Mid S4TM):** Static <em>a</em><sub>0</sub>(0) achieves **<em>R</em><sup>2</sup><sub>1:1</sub> = 0.8068** (98.4% parity; uniform Salpeter 0.6572); dynamic <em>a</em><sub>0</sub>(<em>z</em><sub>d</sub>) achieves **0.8090** (linear <em>R</em><sup>2</sup> = 0.8153).
   * **HST Noise Floor &amp; Joint Stiffness:** Monte Carlo simulations confirm <em>R</em><sup>2</sup><sub>1:1</sub> = 0.7996 already sits at the HST photometric/mass resolution ceiling (**0.8109 &plusmn; 0.0388** at 0.08 dex scatter), while scanning &nu;<sub><em>k</em></sub>(<em>y</em>) = [1 &minus; exp(&minus;<em>y</em><sup><em>k</em>/2</sup>)]<sup>&minus;1/<em>k</em></sup> confirms canonical <em>k</em> = 1.0 sits on the flat shared optimum (<em>k</em> &approx; 0.95&ndash;1.00) of both SLACS IX and SPARC Clean.

2. **Euclid DR1 ECLIPSE High-Redshift Consistency Test (<em>N</em> = 1,095 Clean Grade-A Lenses, Exact &chi;(&infin;) Horizon Bound):**

| Redshift Bin (<em>z</em><sub>d</sub>) | <em>N</em> | Med. &beta;<sub>max</sub> | Baryons GR (% Unphys) | MOND Const. (% Unphys) | &Lambda;TSM Static (% Unphys) | &Lambda;TSM Dynamic <em>a</em><sub>0</sub>(<em>z</em><sub>d</sub>) (% Unphys) | Rel. Reduction vs. MOND |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0.10 &le; <em>z</em><sub>d</sub> &lt; 0.45 | 332 | 0.904 | 85.8% | 65.4% | 64.8% | **60.2%** | &minus;7.8% |
| 0.45 &le; <em>z</em><sub>d</sub> &lt; 0.70 | 345 | 0.849 | 74.5% | 38.3% | 37.7% | **25.8%** | **&minus;32.6%** |
| **0.70 &le; <em>z</em><sub>d</sub> &lt; 0.95** | **282** | 0.799 | 60.6% | 18.8% | 18.8% | **9.2%** | **&minus;50.9%** |
| **0.95 &le; <em>z</em><sub>d</sub> &lt; 1.50** | **136** | 0.756 | 40.4% | 16.2% | 15.4% | **8.1%** | **&minus;50.0%** |
| **ALL (0.10&ndash;1.50)** | **1,095** | 0.840 | 70.1% | 38.7% | 38.3% | **29.8%** | **&minus;23.1%** |

   *Group-Scale Transition Isolated at Low z:* Normal galaxy-scale lenses (&theta;<sub>E</sub> &lt; 0.9&Prime;, <em>N</em> = 107, med log <em>M</em><sub>*</sub> = 11.04) drop to 41.1% (med &beta;<sub>req</sub> = 0.806), whereas giant group-scale deflectors (&theta;<sub>E</sub> &ge; 1.3&Prime;, <em>N</em> = 106, med log <em>M</em><sub>*</sub> = 11.56) surge to 84.9%, pinpointing the transition where hot intragroup gas and nonlinear DBI clustering enter.

3. **SPARC Disk Rotation Curves &amp; BTFR (<em>N</em> = 175 Galaxies, 3,391 Points):**
   * **Locked SPS (&Upsilon;<sub>d</sub> = 0.50, &Upsilon;<sub>b</sub> = 0.70, signed <em>V</em><sub>gas</sub>):** Consolidated **<em>R</em><sup>2</sup><sub>1:1</sub> = 0.9324** across all 175 galaxies (vs. 0.5540 Newtonian); Clean (<em>Q</em> &lt; 3, <em>i</em> &ge; 30&deg;, <em>N</em> = 153): **<em>R</em><sup>2</sup><sub>1:1</sub> = 0.9331**, per-galaxy median **+0.5193** (vs. &minus;3.4838 Newtonian), RAR RMS = **0.1797 dex**, bias &minus;0.0284 dex; Strict <em>Q</em> = 1 (<em>N</em> = 93): RMS = **0.1558 dex**, bias &minus;0.0090 dex.
   * **Asymptotic BTFR:** Clean (<em>N</em> = 123) **<em>R</em><sup>2</sup><sub>1:1</sub> = 0.8819** (0.0671 dex); strictly <em>H</em><sub>0</sub>-independent direct-distance subset (<em>f</em><sub>D</sub> &ge; 2, <em>N</em> = 57) **<em>R</em><sup>2</sup><sub>1:1</sub> = 0.8832** (0.0547 dex), outperforming <em>cH</em><sub>0</sub>/(2&pi;) (0.8433&ndash;0.8653).
   * **Relaxed SPS (&Upsilon;<sub>d</sub> &isin; [0.10, 1.00], &Upsilon;<sub>b</sub> = 1.4&Upsilon;<sub>d</sub>):** **<em>R</em><sup>2</sup><sub>1:1</sub> = 0.9759** all / **0.9763** clean (median &Upsilon;<sub>d</sub> = 0.428, RAR RMS = 0.1330 dex).

4. **LVDB 3D Spheroids (<em>N</em> = 158 Systems):**
   * Universal unscaled <em>a</em><sub>0</sub> reduces gas-free classical dSph &sigma;<sub>los</sub> RMS from 0.7246 dex (Newtonian) to **0.1634 dex** (**0.1600 dex** across <em>N</em> = 24 strictly external 0% Milky Way dSphs; **0.130 dex** across <em>N</em> = 14 isolated field dwarfs), while decisively falsifying morphological flattening shutoff gates.

### B. Part I Galactic Figures (`./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/`)

| SLACS Exact LOS Strong Lensing (<em>N</em> = 110) | SPARC Rotation Curves, RAR, BTFR &amp; SPS (<em>N</em> = 175) |
| :---: | :---: |
| ![SLACS 110 Audit](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/SLACS_110_Figure_v5.png) | ![SPARC Unified Audit](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/SPARC_Unified_Figure_v5.png) |

| LVDB 3D Spheroids &amp; Flattening Falsification (<em>N</em> = 158) | Joint Stiffness Scan (<em>k</em> = 1.0) &amp; HST Noise Floor |
| :---: | :---: |
| ![LVDB Flattening Audit](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/LVDB_Flattening_Figure_v5.png) | ![Stiffness &amp; Noise Floor](./Part_I_Topological_Saturation_of_Spatial_Metrics_Package/SLACS_SPARC_Stiffness_NoiseFloor_v6.png) |

---

## 4. Falsifiable Satellite-Era Predictions &amp; Roadmap to Part III

| Frontier / Test | Target Surveys | Sharp &Lambda;TSM Benchmark / Open Requirement |
| :--- | :--- | :--- |
| **1. Cosmological Red-Line Lock** | DESI DR2/DR3 + Simons Observatory + CMB-S4 | Rigid &omega;<sub>phantom</sub>/&omega;<sub>b</sub> = 5.36710706 tested directly against the current +1.46&sigma; DESI pull (5.2836 &plusmn; 0.0572). |
| **2. Euclid DR1/DR2 Spectroscopic Blind Test** | Euclid (ECLIPSE 15,072 lenses + NISP-S) + Roman | Tests dynamic <em>a</em><sub>0</sub>(<em>z</em><sub>d</sub>) = <em>cH</em>(<em>z</em><sub>d</sub>)/<em>R</em><sup>&deg;</sup><sub>4/3</sub> and the +17.1% separation between <em>cH</em><sub>0</sub>/<em>R</em><sup>&deg;</sup><sub>4/3</sub> and <em>cH</em><sub>0</sub>/(2&pi;). |
| **3. High-<em>z</em> Disk Kinematics** | JWST NIRSpec + ALMA + KMOS3D (<em>z</em> = 1&ndash;6) | Predicts an <em>E</em>(<em>z</em>)-scaled asymptotic BTFR shift <em>V</em><sub>&infin;</sub>(<em>z</em>) &prop; [<em>E</em>(<em>z</em>)]<sup>1/4</sup> on extended outer disks (<em>g</em><sub>bar</sub> &#8818; <em>a</em><sub>0</sub>(<em>z</em>)). |
| **4. External Field Effect &amp; Solar System (<em>Q</em><sub>2</sub>)** | Euclid Weak Lensing + Gaia DR4 + Cassini/DE440 | Full 3D PDE evaluation of <em>Q</em><sub>2</sub> confronting the 3&ndash;15&sigma; Cassini/DE440 bounds (Park et al. 2026) and remote small-body constraints (Vokrouhlick&yacute; et al. 2024). |
| **5. Part III: Relativistic Lensing &amp; Cluster Collapse** | X-COP, CLASH, Euclid Groups &amp; Bullet Cluster | (i) Disformal / unit-timelike vector (TeVeS/AeST) completion beyond the Conformal Lensing No-Go limit (&Psi; = &Phi; = &Phi;<sub>eff</sub>, <em>c</em><sub>T</sub> = <em>c</em>); (ii) galactic/cluster DBI scale separation avoiding double-counting via Khoury (2015) gradient pressure (&#124;&nabla;&phi;&#124;<sup>2</sup> &rArr; <em>c</em><sub>s</sub> &rarr; <em>c</em>) and an intermediate tension window (<em>M</em> &sim; 0.15&ndash;0.24 meV, <em>M</em><sub>J</sub>(0) &sim; 2&times;10<sup>11</sup>&ndash;5.6&times;10<sup>13</sup> <em>M</em><sub>&#9737;</sub>). |

---

## 5. Step-by-Step Reproducibility Guide

### Prerequisites

<pre><code>pip install numpy scipy pandas matplotlib openpyxl cython</code></pre>

### A. Reproduce Part I Galactic &amp; High-Redshift Audits (`Part_I_Topological_Saturation_of_Spatial_Metrics_Package`)

<pre><code>cd Part_I_Topological_Saturation_of_Spatial_Metrics_Package

# 1. SLACS Exact LOS Strong Lensing Audit (N = 110)
python3 SLACS_110_Code_v5.py

# 2. Euclid DR1 ECLIPSE High-Redshift Consistency Test (N = 1,095, exact chi(inf))
python3 run_eclipse_golden_test_v2.py

# 3. Joint SLACS + SPARC Transition-Stiffness Scan &amp; HST Monte Carlo Noise Floor
python3 SLACS_SPARC_Stiffness_NoiseFloor_v6.py

# 4. Full SPARC Rotation Curves, RAR, BTFR &amp; Relaxed SPS Audit (N = 175)
python3 SPARC_Unified_Audit_v5.py

# 5. LVDB 3D Dwarf Spheroidals &amp; Morphological Flattening Falsification (N = 158)
python3 LVDB_Morphology_Audit_v5.py</code></pre>

*(You can also open `Phantom_Clean_Master_SLACS_SPARC-v4.xlsx` in Microsoft Excel to explore live-formula dashboards across SLACS, Euclid DR1 ECLIPSE, and SPARC.)*

### B. Reproduce Part II `CLASS v3.4.0` Baseline &amp; Full `Cobaya` MCMC Suite (`Part_II_The_Topological_Saturation_Model_Package`)

<pre><code># 1. Clone official CLASS v3.4.0 and apply the clean covariant Part II C files
git clone https://github.com/lesgourg/class_public.git
cd class_public
cp ../Part_II_The_Topological_Saturation_Model_Package/source/background.c source/background.c
cp ../Part_II_The_Topological_Saturation_Model_Package/source/perturbations.c source/perturbations.c

# 2. Compile CLASS and install the Python wrapper (classy)
make clean &amp;&amp; make -j4
cd python &amp;&amp; python3 setup.py install --user &amp;&amp; cd ../../Part_II_The_Topological_Saturation_Model_Package

# 3. Verify CLASS baseline CMB (Wishart dchi2 = 0.1446), BBN (0.00000%) &amp; noiseless MCMC check
python3 00_test_classy_consistency.py
python3 audit_phantom_full.py
python3 phase1_calibrate_and_bbn.py
python3 phase2_mcmc_engine.py

# 4. Install real likelihoods and reproduce the 4-run Cobaya MCMC suite (Planck 2018 + DESI + Pantheon+)
bash 01_setup_cobaya_ubuntu.sh
bash 02_run_all_chains.sh
python3 03_analyze_cobaya_results.py
python3 04_regenerate_clean_figures.py</code></pre>

---

## Citation

If you use this code, data, or manuscripts in your research, please cite Version 8 (November 24, 2026):

<pre><code>@misc{haimovich2026_part1,
  author       = {Haimovich, Tomer},
  title        = {Topological Saturation of Spatial Metrics: A Unified Single-Threshold Phenomenological Audit of Strong Lensing, Disk Kinematics, and 3D Dwarf Spheroidals (Version 8)},
  year         = {2026},
  month        = nov,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.24463759},
  url          = {https://doi.org/10.5281/zenodo.24463759}
}

@misc{haimovich2026_part2,
  author       = {Haimovich, Tomer},
  title        = {The Topological Saturation Model ($\Lambda$TSM): Clean Covariant DBI Cosmology and Full Cobaya MCMC Audit (Planck 2018 + DESI BAO + Pantheon+) (Version 8)},
  year         = {2026},
  month        = nov,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.24463760},
  url          = {https://doi.org/10.5281/zenodo.24463760}
}</code></pre>