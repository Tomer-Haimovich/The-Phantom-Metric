# Zero-Parameter Topological Saturation Model (LambdaTSM) — CLASS v3.4.0 Reproducibility Suite

## Overview
This package contains the complete C source modifications, Boltzmann solver calibration scripts, Big Bang Nucleosynthesis (BBN) expansion proofs, and Markov Chain Monte Carlo (MCMC) statistical audits for **Part II**: *"The Topological Saturation Model (LambdaTSM): Resolving the Cosmological Horizon, CMB Acoustic Power Spectra, and Late-Universe Structure via Zero-Parameter Covariant DBI Field Theory"* (Tomer Haimovich, 2026).

### Key Theoretical & Numerical Architecture
1. **4D Orthonormal Tetrad Ratio (`omega_phantom / omega_b = 5.367134`):**
   In `source/background.c`, the background topological stress-energy density is locked analytically to the baryon density via the $D^2 = 16$ independent algebraic components of the 4D orthonormal tetrad frame field $e^a_\mu$ on $T\mathcal{M}_4$ projected onto the $(D-1) = 3$ spatial hypersurface $\Sigma_3$, dressed by the 1-loop de Sitter horizon anomaly $(1 + 1/(16\pi^2))$.
2. **Covariant DBI Zero-Sound-Speed Limit ($c_s^2 = \gamma^{-2} \to 0$):**
   In the kinetically saturated DBI regime ($2|X| \to M^4, w_{\text{eff}} \to 0, c_s^2 \to 0$), the covariant scalar stress-energy conservation law in the synchronous comoving frame ($\theta_\phi = 0$) reduces identically to the comoving dust continuity equation ($\delta'_\phi = -h'/2$), requiring zero particle cold dark matter ($\omega_{\text{cdm}} = 0$).
3. **Sub-Horizon Acoustic Impedance & Smooth Covariant Expansion Gate:**
   In `source/perturbations.c`, acoustic-scale spatial metric gradients couple to the baryon velocity divergence $\theta_b$ via the parameter-free 3D Maclaurin elasticity $\mathcal{K}_0 = 2/3$, the Einstein-DBI horizon impedance $\Gamma_0 = 1/(8\pi)$, the Helmholtz-Yukawa causal projection operator $k^2/(k^2 + \mathcal{H}^2)$, and the $C^\infty$ smooth expansion gate $\frac{1}{2}[1 + \tanh(\theta_b / \theta_{\text{tol}})]$.

---

## Reproducibility & Execution Guide

### 1. Immediate Checkpoint Audit (No C Compilation Required)
To immediately verify the 5-parameter MCMC chains ($N_{\text{samples}} = 4,800$), regenerate the posterior Corner Plot (`phantom_mcmc_corner.png`), and output the full statistical audit table (`mcmc_audit_report.txt`) from the archived checkpoint (`mcmc_chain_checkpoint.npz`):
```bash
python3 phase2_mcmc_engine.py
```

### 2. Full Live Boltzmann & MCMC Execution from Scratch (`CLASS v3.4.0`)
To compile the modified `CLASS v3.4.0` Boltzmann solver and run all scripts live from scratch:
1. Clone `CLASS v3.4.0` and replace `source/background.c` and `source/perturbations.c` with the modified C files provided in this repository:
   ```bash
   git clone https://github.com/lesgourg/class_public.git
   cd class_public
   cp /path/to/LambdaTSM_ZeroParam_Package/background.c source/background.c
   cp /path/to/LambdaTSM_ZeroParam_Package/perturbations.c source/perturbations.c
   make clean && make -j4
   pip install .
   ```
2. Run the 4-panel Planck CMB ($TT, TE, EE$) and Matter Power Spectrum $P(k)$ diagnostic audit ($R^2_{TT} = 0.9970, R^2_{TE} = 0.9951, R^2_{EE} = 0.9951, \sigma_8 = 0.8229$):
   ```bash
   python3 audit_phantom_full.py
   ```
3. Run the Primordial Big Bang Nucleosynthesis ($z = 10^9$) expansion proof ($\Delta H/H = 0.00000\%$) and target calibration:
   ```bash
   python3 phase1_calibrate_and_bbn.py
   ```
4. Run the 5-parameter affine-invariant ensemble MCMC engine (supports `LAMBDATSM_MCMC_MODE=joint` [default] or `LAMBDATSM_MCMC_MODE=cmbonly`):
   ```bash
   python3 phase2_mcmc_engine.py
   ```
