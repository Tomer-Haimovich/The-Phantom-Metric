# Zero-Parameter Topological Saturation Model (LambdaTSM) — CLASS v3.4.0 Reproducibility Suite

## Overview
This package contains the complete C source modifications, Boltzmann solver calibration scripts, Big Bang Nucleosynthesis (BBN) expansion proofs, and Markov Chain Monte Carlo (MCMC) statistical audits for **The Topological Saturation Model (LambdaTSM)**.

### Key Theoretical & Numerical Architecture
1. **4D Orthonormal Tetrad Ratio (omega_phantom / omega_b = 5.367134):**
   In `source/background.c`, the background topological stress-energy density is locked analytically to the baryon density via the D^2 = 16 independent algebraic components of the 4D orthonormal tetrad frame field e^a_mu on TM_4 projected onto the (D-1) = 3 spatial hypersurface Sigma_3, dressed by the 1-loop de Sitter horizon anomaly (1 + 1/(16*pi^2)).
2. **Covariant DBI Zero-Sound-Speed Limit (c_s^2 = gamma^-2 -> 0):**
   In the kinetically saturated DBI regime (2X -> M^4, w_eff -> 0, c_s^2 -> 0), the covariant scalar stress-energy conservation law in the synchronous comoving frame (theta_phi = 0) reduces identically to the comoving dust continuity equation (delta'_phi = -h'/2), requiring zero particle cold dark matter (omega_cdm = 0).
3. **Sub-Horizon Acoustic Impedance & Smooth Covariant Expansion Gate:**
   In `source/perturbations.c`, acoustic-scale spatial metric gradients couple to the baryon velocity divergence theta_b via the parameter-free 3D Maclaurin elasticity K_0 = 2/3, the Einstein-DBI horizon impedance Gamma_0 = 1/(8*pi), the Helmholtz-Yukawa causal projection operator k^2/(k^2 + H^2), and the C^infinity smooth expansion gate 0.5*[1 + tanh(theta_b / theta_tol)].
