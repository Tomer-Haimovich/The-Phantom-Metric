# Topological Saturation of Spatial Metrics (Part I) — Reproducibility Suite

## Overview
This repository contains the complete Python verification scripts, observational data tables, and high-resolution figures for **Part I**: *"Topological Saturation of Spatial Metrics Unifies Galactic Kinematics and Strong Gravitational Lensing"* (Tomer Haimovich, 2026).

## Theoretical Constants (Zero Free Dark-Sector Parameters)
* **Cosmic Horizon Acceleration Threshold:** $a_0 = \frac{cH_0}{2\pi} \approx 1.1288 \times 10^{-10}\text{ m s}^{-2}$ (calibrated in the unperturbed late-universe Hubble-flow regime $D > 20\text{ Mpc}$, $H_0 = 73.0\text{ km s}^{-1}\text{ Mpc}^{-1}$).
* **Planar KAM Golden-Ratio Threshold:** $a_\varphi = a_0 \varphi \approx 1.8264 \times 10^{-10}\text{ m s}^{-2}$, where $\varphi = \frac{1+\sqrt{5}}{2} \approx 1.618034$.

## Repository Contents & Execution Guide
1. **Bulge-Less Disk Verification (NGC 2403):**
   * Script: `SPARC_NGC2403_Code.py`
   * Input Data: `NGC2403_rotmod.dat`
   * Output Table & Figure: `NGC2403_Phantom_Metric_ML046.csv`, `SPARC_NGC2403_Figure.png` ($R^2 = 0.9694$ at locked $\Upsilon_* = 0.46$).
2. **Unfiltered Universal SPARC Galactic Rotation Audit (175 Galaxies, 3,391 Measurements):**
   * Locked Baseline Script ($\Upsilon_* = 0.46$): `SPARC_Universal Global_175G_Code.py` (reads `Rotmod_LTG.zip` / `SPARC_Data.xlsx`, outputs `SPARC_Global_Phantom_Metric_ML046.csv` and `SPARC_Universal Global_175G_Figure.png`, achieving $R^2 = 0.9150$).
   * Stellar Population Synthesis Relaxed Audit ($\Upsilon_* \in [0.1, 1.0]$): `SPARC_Universal Global_175G_M-L_Code.py` (outputs `SPARC_Optimized_ML_Results.xlsx` and `SPARC_Optimized_R2.png`, achieving $R^2 = 0.9698$).
3. **SLACS Strong Gravitational Lensing Audit ($N = 100$ Early-Type Galaxies):**
   * Script: `SLACS_100_Code.py`
   * Dataset & Output: `SLACS_Phantom_Metric_Lensing.csv`, `SLACS_100_Figure.png` ($R^2 = 0.9917$, $0\%$ dark matter).
4. **CLASH Galaxy Cluster Strong-Lensing Audit ($N = 25$ Clusters: 22 Relaxed, 3 Unrelaxed Mergers):**
   * Script: `Clash_25_Code.py`
   * Input Catalog (Zitrin et al. 2015): `Clash_Data.csv`
   * Output Predictions & Figure: `Phantom_Metric_Clusters_Predictions.csv`, `Clash_25_Figure.png` (achieving strict 1:1 parity $R^2_{1:1} = 0.9710$ across all morphologically relaxed clusters using solely the $13\%$ luminous baryonic mass fraction).
