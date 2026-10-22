# Zero-Parameter Topological Saturation Model (LambdaTSM) — CLASS v3.4.0 Reproducibility Package

**Author:** Tomer Haimovich (ORCID: 0009-0002-1536-1397)  
**Zenodo DOI:** 10.5281/zenodo.23445631  
**Companion Study DOI:** 10.5281/zenodo.23436147  

## Contents of this Archive
1. `background.c` — Modified CLASS v3.4.0 background module implementing the exact zero-parameter topological density ratio $\omega_{\mathrm{phantom}}/\omega_b = \frac{16}{3}(1 + \frac{1}{16\pi^2}) \approx 5.367134$ when particle cold dark matter is absent ($\Omega_{\mathrm{cdm}} = 10^{-12}$).
2. `perturbations.c` — Modified CLASS v3.4.0 perturbations module implementing the 3D Maclaurin curvature elasticity $\mathcal{K}_0 = \frac{2}{3}$, Einstein-DBI horizon impedance $\Gamma_0 = \frac{1}{8\pi}$, causal sub-horizon gate $\Pi_{\mathrm{causal}} = k^2/(k^2+\mathcal{H}^2)$, and recombination freeze-out gate $\Theta_{\mathrm{rec}}(a) = \exp[-(a/\pi a_{\mathrm{eq}})^4]$ inside the decoupling `tca_off` regime.
3. `audit_phantom_full.py` & `phantom_full_diagnostic.png` — Automated zero-parameter baseline verification script and 4-panel CMB ($TT, TE, EE$) + Matter Power Spectrum $P(k)$ diagnostic figure ($R^2_{TT}=0.9970, R^2_{TE}=0.9951, R^2_{EE}=0.9951, \sigma_8=0.8229$).
4. `phase1_calibrate_and_bbn.py` & `phantom_bbn_expansion_proof.png` — Early-universe Big Bang Nucleosynthesis ($z = 10^9$) expansion proof ($\Delta H/H = 0.00000\%$ at $z=10^9$) and covariance calibration target (`mcmc_planck_target.npz`).
5. `phase2_mcmc_engine.py`, `mcmc_chain_checkpoint.npz`, `phantom_mcmc_corner.png`, & `mcmc_audit_report.txt` — Resumable 5-parameter affine-invariant ensemble MCMC engine ($0\%$ CDM, $N_{\mathrm{samples}}=4,800$), raw chain checkpoint, 6D posterior corner plot, and statistical audit report ($R^2_{TT}=0.9991, R^2_{TE}=0.9954, R^2_{EE}=0.9957$, Bandpower RMS $= 0.405\sigma$ per bin).
