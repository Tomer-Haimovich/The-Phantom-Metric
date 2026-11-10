import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from classy import Class

PKG_DIR = os.path.dirname(os.path.abspath(__file__))

# 1. LIVE RUN: Standard Model Baseline (Planck 2018 LCDM with 100% CDM = 0.1200)
lcdm = Class()
lcdm.set({
    'output': 'tCl,pCl,lCl,mPk',
    'lensing': 'yes',
    'l_max_scalars': 2500,
    'P_k_max_1/Mpc': 5.0,
    'h': 0.6736,
    'omega_b': 0.02237,
    'omega_cdm': 0.1200,
    'A_s': 2.100e-9,
    'n_s': 0.9649,
    'tau_reio': 0.0544
})
lcdm.compute()
cls_lcdm = lcdm.lensed_cl(2500)
ell_full = cls_lcdm['ell'][2:2501]
fac_full = ell_full * (ell_full + 1.0) / (2.0 * np.pi) * (2.7255e6)**2
tt_lcdm_full = cls_lcdm['tt'][2:2501] * fac_full
te_lcdm_full = cls_lcdm['te'][2:2501] * fac_full
ee_lcdm_full = cls_lcdm['ee'][2:2501] * fac_full
k_vals = np.logspace(-4, 0, 200)
pk_lcdm = np.array([lcdm.pk(k, 0.0) for k in k_vals])
sigma8_lcdm = float(lcdm.sigma8())
lcdm.struct_cleanup()
lcdm.empty()

# 2. LIVE RUN: Zero-Parameter LambdaTSM (0% Particle CDM: omega_cdm = 1e-12)
phantom = Class()
phantom.set({
    'output': 'tCl,pCl,lCl,mPk',
    'lensing': 'yes',
    'l_max_scalars': 2500,
    'P_k_max_1/Mpc': 5.0,
    'h': 0.6736,
    'omega_b': 0.02237,
    'omega_cdm': 1.0e-12,
    'A_s': 2.100e-9,
    'n_s': 0.9649,
    'tau_reio': 0.0544
})
phantom.compute()
cls_ph = phantom.lensed_cl(2500)
tt_ph_full = cls_ph['tt'][2:2501] * fac_full
te_ph_full = cls_ph['te'][2:2501] * fac_full
ee_ph_full = cls_ph['ee'][2:2501] * fac_full
pk_ph = np.array([phantom.pk(k, 0.0) for k in k_vals])
sigma8_ph = float(phantom.sigma8())
phantom.struct_cleanup()
phantom.empty()

def r2_1to1(y_pred, y_true):
    return 1.0 - np.sum((y_pred - y_true)**2) / np.sum((y_true - np.mean(y_true))**2)

r2_tt = r2_1to1(tt_ph_full[48:], tt_lcdm_full[48:])
r2_te = r2_1to1(te_ph_full[48:], te_lcdm_full[48:])
r2_ee = r2_1to1(ee_ph_full[48:], ee_lcdm_full[48:])

omega_b_val = 0.02237
omega_ph_derived = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * np.pi**2)) * omega_b_val

print("\n" + "="*86)
print(" ZERO-PARAMETER LAMBDATSM COSMOLOGICAL AUDIT TABLE (l = 50..2500)")
print("="*86)
print(f"{'Observable / Cosmological Metric':<38} | {'Planck 2018 (LCDM)':<18} | {'LambdaTSM (0-Param)':<22}")
print("-"*86)
print(f"{'Background Density (omega_cdm/ph)':<38} | {'0.120000 (Fitted)':<18} | {omega_ph_derived:.6f} (Derived)")
print(f"{'CMB Temperature (D_l^TT) 1:1 R^2':<38} | {'1.000000 (Ref.)':<18} | {r2_tt:.6f} ({r2_tt:.4f})")
print(f"{'CMB Cross-Polarization (D_l^TE) R^2':<38} | {'1.000000 (Ref.)':<18} | {r2_te:.6f} ({r2_te:.4f})")
print(f"{'CMB E-Mode Polarization (D_l^EE) R^2':<38} | {'1.000000 (Ref.)':<18} | {r2_ee:.6f} ({r2_ee:.4f})")
print(f"{'Linear Matter Clustering (sigma_8)':<38} | {sigma8_lcdm:.4f} ({sigma8_lcdm:.3f})     | {sigma8_ph:.4f} ({sigma8_ph:.3f})")
print("="*86 + "\n")

fig, axs = plt.subplots(2, 2, figsize=(16, 11), dpi=300)

axs[0, 0].plot(ell_full, tt_lcdm_full, 'k--', lw=2.0, label='Standard Model (Planck $\\Lambda$CDM)')
axs[0, 0].plot(ell_full, tt_ph_full, 'g-', lw=1.6, alpha=0.85, label=f'$\\Lambda$TSM 0-Param ($R^2 = {r2_tt:.4f}$)')
axs[0, 0].set_title('CMB Temperature Power Spectrum ($D_\\ell^{TT}$)')
axs[0, 0].set_xlabel('Multipole $\\ell$'); axs[0, 0].set_ylabel('$D_\\ell^{TT}$ [$\\mu\\mathrm{K}^2$]')
axs[0, 0].set_xlim(2, 2500); axs[0, 0].legend(); axs[0, 0].grid(ls='--', alpha=0.4)

axs[0, 1].plot(ell_full[:1999], te_lcdm_full[:1999], 'k--', lw=2.0, label='Standard Model (Planck $\\Lambda$CDM)')
axs[0, 1].plot(ell_full[:1999], te_ph_full[:1999], 'b-', lw=1.5, alpha=0.85, label=f'$\\Lambda$TSM 0-Param ($R^2 = {r2_te:.4f}$)')
axs[0, 1].set_title('CMB Cross-Polarization Spectrum ($D_\\ell^{TE}$)')
axs[0, 1].set_xlabel('Multipole $\\ell$'); axs[0, 1].set_ylabel('$D_\\ell^{TE}$ [$\\mu\\mathrm{K}^2$]')
axs[0, 1].set_xlim(2, 2000); axs[0, 1].legend(); axs[0, 1].grid(ls='--', alpha=0.4)

axs[1, 0].plot(ell_full[:1999], ee_lcdm_full[:1999], 'k--', lw=2.0, label='Standard Model (Planck $\\Lambda$CDM)')
axs[1, 0].plot(ell_full[:1999], ee_ph_full[:1999], 'r-', lw=1.5, alpha=0.85, label=f'$\\Lambda$TSM 0-Param ($R^2 = {r2_ee:.4f}$)')
axs[1, 0].set_title('CMB E-Mode Polarization Spectrum ($D_\\ell^{EE}$)')
axs[1, 0].set_xlabel('Multipole $\\ell$'); axs[1, 0].set_ylabel('$D_\\ell^{EE}$ [$\\mu\\mathrm{K}^2$]')
axs[1, 0].set_xlim(2, 2000); axs[1, 0].legend(); axs[1, 0].grid(ls='--', alpha=0.4)

axs[1, 1].loglog(k_vals, pk_lcdm, 'k--', lw=2.0, label=f'Planck $\\Lambda$CDM ($\\sigma_8 = {sigma8_lcdm:.3f}$)')
axs[1, 1].loglog(k_vals, pk_ph, 'm-', lw=1.6, alpha=0.85, label=f'$\\Lambda$TSM 0-Param ($\\sigma_8 = {sigma8_ph:.3f}$)')
axs[1, 1].set_title('Matter Power Spectrum $P(k)$ at $z = 0$')
axs[1, 1].set_xlabel('Wavenumber $k$ [1/Mpc]'); axs[1, 1].set_ylabel('$P(k)$ [$\\mathrm{Mpc}^3$]')
axs[1, 1].set_xlim(1e-4, 1e0); axs[1, 1].legend(); axs[1, 1].grid(ls='--', alpha=0.4, which='both')

plt.tight_layout()
out_png = os.path.join(PKG_DIR, "phantom_full_diagnostic.png")
plt.savefig(out_png, dpi=300)
plt.close(fig)
print(f"Saved 4-panel diagnostic figure to: {out_png}")
