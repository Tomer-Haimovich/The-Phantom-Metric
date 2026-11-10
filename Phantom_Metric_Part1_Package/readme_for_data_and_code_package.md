# Topological Saturation of Spatial Metrics Unifies Galactic Kinematics and Strong Gravitational Lensing — Reproducibility Package

**Author:** Tomer Haimovich (ORCID: 0009-0002-1536-1397)  
**Repository DOI:** [10.5281/zenodo.24098753](https://doi.org/10.5281/zenodo.24098753)  
**Companion Cosmological Study ($\Lambda\text{TSM}$) DOI:** [10.5281/zenodo.24098774](https://doi.org/10.5281/zenodo.24098774)

---

## Overview

This repository contains the complete, self-contained observational datasets, Python audit scripts, step-by-step physical conversion tables, and high-resolution figures accompanying the manuscript:

> **"Topological Saturation of Spatial Metrics Unifies Galactic Kinematics and Strong Gravitational Lensing"**

All analyses operate under a strict **zero-free-parameter lock** with **$0\%$ non-baryonic dark matter**, governed by the causal Hubble horizon acceleration threshold:
$$a_0 = \frac{c H_0}{2\pi} \approx 1.1288 \times 10^{-10}\text{ m s}^{-2} \quad (H_0 = 73.0\text{ km s}^{-1}\text{ Mpc}^{-1})$$
and the golden-ratio Kolmogorov–Arnold–Moser (KAM) planar resonance threshold in flattened spiral disks:
$$a_\varphi = a_0 \varphi \approx 1.8264 \times 10^{-10}\text{ m s}^{-2} \quad \left(\varphi = \frac{1+\sqrt{5}}{2}\right)$$

---

## Package Contents

### 1. Galaxy-Scale Strong Gravitational Lensing Audit ($N = 110$ Official SLACS & S4TM Lenses)
Combines all $N = 70$ Grade-A early-type lens galaxies with complete multi-band HST stellar mass and SIE lens models from **SLACS IX** (Auger et al. 2009, *ApJ* 705, 1099; Bolton et al. 2008, *ApJ* 682, 964) and all $N = 40$ confirmed Grade-A early-type lens galaxies from **SLACS XIII / S4TM** (Shu et al. 2017, *ApJ* 851, 48).

* **`SLACS_DB.csv`**: Verbatim raw observational catalog ($N = 110$) containing `lens_id`, `Survey`, spectroscopic redshifts (`z_lens`, `z_source`), HST angular Einstein radii (`theta_E_obs_arcsec`), effective radii (`R_eff_arcsec`), stellar velocity dispersions (`sigma_obs_km_s`), and HST stellar population synthesis masses (`log10_Mstar_Chabrier`, `log10_Mstar_Salpeter`).
* **`SLACS_110_Code.py`**: Standalone Python audit script that reads `SLACS_DB.csv`, computes angular diameter distances ($D_L, D_S, D_{LS}$), executes step-by-step SI physical conversions, predicts lens-plane physical Einstein ring radii $R_{E,\text{calc}} = D_L(D_{LS}/D_S)\hat{\alpha}$ with $\hat{\alpha} = 4\pi\sqrt{GM_{\text{baryon}}a_0}/c^2$ (and saturation boundaries $R_{BB} = \sqrt{GM_{\text{baryon}}/a_0}$) with zero free parameters, and outputs the processed CSV, master Excel workbook, and Figure 1.
* **`SLACS_Phantom_Metric_Lensing.csv`**: Complete processed output table with all intermediate SI conversion steps and model-to-observation ratios.
* **`SLACS_110_Figure.png`**: Publication Figure 1 mapping predicted vs. HST-observed physical Einstein ring radii ($R^2 = 0.7955, p = 5.1 \times 10^{-39}$ across $N = 110$; $R^2 = 0.8261$ across SLACS IX $N = 70$; strict 1:1 parity $R^2_{1:1} = 0.7773$; mean ratio $R_{E,\text{obs}}/R_{E,\text{calc}} = 1.016$; ensemble sum $370.3\text{ kpc}$ observed vs. $375.7\text{ kpc}$ predicted).

---

### 2. Cluster-Scale Strong Gravitational Lensing Audit ($N = 25$ CLASH Galaxy Clusters)
Evaluates the Cluster Lensing And Supernova survey with Hubble (CLASH; Zitrin et al. 2015, *ApJ* 801, 44; Postman et al. 2012, *ApJS* 199, 25) using strictly the observed $\sim 13\%$ baryonic mass fraction (intracluster X-ray gas + stellar mass) and $0\%$ dark matter.

* **`Clash_Data.csv`**: Raw observational cluster dataset from Zitrin et al. (2015).
* **`Clash_25_Code.py`**: Python audit script evaluating macroscopic saturation boundaries $R_{BB} = \sqrt{GM_{\text{baryon}}/a_0}$ against observed cluster Einstein radii and identifying unrelaxed multi-node mergers ($\Theta_{\text{relax}} = 0$).
* **`Phantom_Metric_Clusters_Predictions.csv`**: Output table of predicted vs. observed cluster Einstein radii.
* **`Clash_25_Figure.png`**: Publication Figure 2 demonstrating strict 1:1 parity ($R^2_{1:1} = 0.9710, p \approx 1.4 \times 10^{-11}$) across morphologically relaxed clusters.

---

### 3. Galactic Rotation Curve Audits ($N = 175$ SPARC Galaxies, $3,391$ Radial Data Points)
Evaluates $100\%$ of the Spitzer Photometry and Accurate Rotation Curves (SPARC) database (Lelli, McGaugh, & Schombert 2016, *AJ* 152, 157) with zero data pruning.

* **`Rotmod_LTG.zip` & `NGC2403_rotmod.dat`**: Official raw SPARC radial rotation curve files.
* **`SPARC_NGC2403_Code.py`**: Bulge-less galaxy NGC 2403 verification script ($\Upsilon_* = 0.46$ locked, $R^2 = 0.9694$).
* **`NGC2403_Phantom_Metric_ML046.csv` & `SPARC_NGC2403_Figure.png`**: Output table and Publication Figure 3 for NGC 2403.
* **`SPARC_Universal Global_175G_Code.py`**: Full 175-galaxy ($N = 3,391$ points) zero-parameter audit script with universally locked mass-to-light ratio $\Upsilon_* = 0.46$ ($R^2 = 0.9150, p < 10^{-300}$ point-wise; $p \approx 1.1 \times 10^{-95}$ galaxy-wise).
* **`SPARC_Global_Phantom_Metric_ML046.csv` & `SPARC_Data.xlsx`**: Consolidated 3,391-point kinematic audit tables (locked $\Upsilon_* = 0.46$).
* **`SPARC_Universal Global_175G_Figure.png`**: Publication Figure 4 (Locked $\Upsilon_* = 0.46$ SPARC audit).
* **`SPARC_Universal Global_175G_M-L_Code.py`**: Relaxed stellar population synthesis audit script ($\Upsilon_* \in [0.1, 1.0]$ per galaxy, $R^2 = 0.9698, p < 10^{-300}$ point-wise; $p \approx 4.8 \times 10^{-132}$ galaxy-wise).
* **`SPARC_Optimized_ML_Results.xlsx` & `SPARC_Optimized_R2.png`**: Consolidated output workbook and Publication Figure 5 for the relaxed SPARC audit.

---

## Requirements & Execution

All scripts require standard scientific Python 3 packages (`numpy`, `pandas`, `scipy`, `matplotlib`, `openpyxl`):

```bash
pip install numpy pandas scipy matplotlib openpyxl
```

To reproduce all tables and figures from scratch:

```bash
python SLACS_110_Code.py
python Clash_25_Code.py
python SPARC_NGC2403_Code.py
python "SPARC_Universal Global_175G_Code.py"
python "SPARC_Universal Global_175G_M-L_Code.py"
```