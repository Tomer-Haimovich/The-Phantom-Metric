import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from classy import Class
from cmb_likelihood_utils import (
    L_MIN, L_MAX_TT, L_MAX_POL, F_SKY_PLANCK,
    compute_planck_noise_dl, wishart_chi2_exact, binned_knox_chi2
)

PKG_DIR = os.path.dirname(os.path.abspath(__file__))
REF_FILE = os.path.join(PKG_DIR, "mcmc_planck_target.npz")

print("=" * 78)
print(" PHASE 1 (CLEAN v2): EXACT KNOX/WISHART TARGET CALIBRATION & BBN PROOF (z=10^9)")
print("=" * 78)

t0 = time.time()
lcdm = Class()
lcdm.set({
    'output': 'tCl,pCl,lCl,mPk',
    'lensing': 'yes',
    'l_max_scalars': L_MAX_TT,
    'P_k_max_1/Mpc': 5.0,
    'h': 0.6736,
    'omega_b': 0.02237,
    'omega_cdm': 0.1200,
    'A_s': 2.100e-9,
    'n_s': 0.9649,
    'tau_reio': 0.0544
})
lcdm.compute()
dt_single = time.time() - t0

cls_lcdm = lcdm.lensed_cl(L_MAX_TT)
ell = cls_lcdm['ell'][L_MIN:L_MAX_TT + 1].astype(float)
fac_dl = ell * (ell + 1.0) / (2.0 * np.pi)
fac = fac_dl * (2.7255e6)**2
tt_ref = cls_lcdm['tt'][L_MIN:L_MAX_TT + 1] * fac
te_ref = cls_lcdm['te'][L_MIN:L_MAX_TT + 1] * fac
ee_ref = cls_lcdm['ee'][L_MIN:L_MAX_TT + 1] * fac
sigma8_ref = float(lcdm.sigma8())
der_l = lcdm.get_current_derived_parameters(['rs_rec', 'z_rec'])

bg_lcdm = lcdm.get_background()
z_bg_l = bg_lcdm['z']
ord_l = np.argsort(z_bg_l)
z_grid = np.logspace(-2, 9, 400)
H_lcdm = np.interp(z_grid, z_bg_l[ord_l], bg_lcdm['H [1/Mpc]'][ord_l])
lcdm.struct_cleanup()
lcdm.empty()

phantom = Class()
phantom.set({
    'output': 'tCl,pCl,lCl,mPk',
    'lensing': 'yes',
    'l_max_scalars': L_MAX_TT,
    'P_k_max_1/Mpc': 5.0,
    'h': 0.6736,
    'omega_b': 0.02237,
    'omega_cdm': 1.0e-12,
    'A_s': 2.100e-9,
    'n_s': 0.9649,
    'tau_reio': 0.0544
})
phantom.compute()
cls_ph = phantom.lensed_cl(L_MAX_TT)
tt_ph = cls_ph['tt'][L_MIN:L_MAX_TT + 1] * fac
te_ph = cls_ph['te'][L_MIN:L_MAX_TT + 1] * fac
ee_ph = cls_ph['ee'][L_MIN:L_MAX_TT + 1] * fac
sigma8_ph = float(phantom.sigma8())
der_p = phantom.get_current_derived_parameters(['rs_rec', 'z_rec'])

bg_ph = phantom.get_background()
z_bg_p = bg_ph['z']
ord_p = np.argsort(z_bg_p)
H_ph = np.interp(z_grid, z_bg_p[ord_p], bg_ph['H [1/Mpc]'][ord_p])
phantom.struct_cleanup()
phantom.empty()

# Exact Knox (1995) detector noise — ZERO hidden 1.5% or foreground padding!
n_tt_dl, n_ee_dl = compute_planck_noise_dl(ell)
var_tt_knox = (2.0 / ((2.0 * ell + 1.0) * F_SKY_PLANCK)) * (tt_ref + n_tt_dl)**2
var_ee_knox = (2.0 / ((2.0 * ell + 1.0) * F_SKY_PLANCK)) * (ee_ref + n_ee_dl)**2
var_te_knox = (1.0 / ((2.0 * ell + 1.0) * F_SKY_PLANCK)) * (te_ref**2 + (tt_ref + n_tt_dl) * (ee_ref + n_ee_dl))

chi2_w0 = wishart_chi2_exact(tt_ph, te_ph, ee_ph, tt_ref, te_ref, ee_ref, ell, n_tt_dl, n_ee_dl)
knox_b0 = binned_knox_chi2(tt_ph, te_ph, ee_ph, tt_ref, te_ref, ee_ref, ell, n_tt_dl, n_ee_dl, dl=30)

np.savez(
    REF_FILE,
    ell=ell, fac=fac,
    tt_ref=tt_ref, te_ref=te_ref, ee_ref=ee_ref,
    tt_ph0=tt_ph, te_ph0=te_ph, ee_ph0=ee_ph,
    n_tt_dl=n_tt_dl, n_ee_dl=n_ee_dl,
    var_tt=var_tt_knox, var_te=var_te_knox, var_ee=var_ee_knox,
    sigma8_ref=sigma8_ref, sigma8_ph0=sigma8_ph,
    rs_rec_lcdm=der_l['rs_rec'], rs_rec_ph=der_p['rs_rec']
)

frac_H_diff_bbn = float(np.abs(H_ph[-1] - H_lcdm[-1]) / H_lcdm[-1] * 100.0)
max_H_diff_all = float(np.max(np.abs(H_ph - H_lcdm) / H_lcdm) * 100.0)
frac_rs_diff = float(np.abs(der_p['rs_rec'] - der_l['rs_rec']) / der_l['rs_rec'] * 100.0)

fig, axs = plt.subplots(1, 2, figsize=(14.5, 5.5), dpi=300)
axs[0].loglog(z_grid, H_lcdm * 2.99792458e5, 'k--', lw=2.2, label='Planck 2018 $\\Lambda$CDM')
axs[0].loglog(z_grid, H_ph * 2.99792458e5, color='#0f766e', lw=1.6, alpha=0.9, label='$\\Lambda$TSM Clean v2 (0% CDM)')
axs[0].axvline(1e9, color='#be123c', ls=':', alpha=0.8, label='BBN Epoch ($z \\sim 10^9$)')
axs[0].axvline(3400, color='#1d4ed8', ls=':', alpha=0.8, label='Equality ($z_{\\mathrm{eq}} \\approx 3400$)')
axs[0].set_title('Cosmological Expansion Rate $H(z)$ from $z=0$ to $z=10^9$', fontsize=10.8, fontweight='bold', pad=8)
axs[0].set_xlabel('Redshift $z$'); axs[0].set_ylabel('$H(z)$ [km s$^{-1}$ Mpc$^{-1}$]')
axs[0].legend(); axs[0].grid(ls=':', alpha=0.5, which='both')

rel_H_pct = (H_ph - H_lcdm) / H_lcdm * 100.0
axs[1].semilogx(z_grid, rel_H_pct, color='#7c3aed', lw=1.8, label='$\\Delta H(z)/H_{\\Lambda\\mathrm{CDM}}(z)$ [%]')
axs[1].axhline(0.0, color='k', ls='--', lw=1.1)
axs[1].set_title(
    f'Fractional Expansion Deviation (BBN $z=10^9$: {frac_H_diff_bbn:.5f}%, Max: {max_H_diff_all:.4f}%)',
    fontsize=10.8, fontweight='bold', pad=8
)
axs[1].set_xlabel('Redshift $z$'); axs[1].set_ylabel('Fractional Deviation $\\Delta H / H$ [%]')
axs[1].set_ylim(-0.05, 0.05); axs[1].legend(); axs[1].grid(ls=':', alpha=0.5, which='both')

plt.tight_layout(w_pad=2.5)
bbn_png = os.path.join(PKG_DIR, "phantom_bbn_expansion_proof.png")
plt.savefig(bbn_png, dpi=300, bbox_inches='tight', pad_inches=0.15)
plt.close(fig)

print(f"Single CLASS Evaluation Time         : {dt_single:.2f} s")
print(f"Sound Horizon r_s(z_*)               : {der_p['rs_rec']:.4f} Mpc (LCDM: {der_l['rs_rec']:.4f} Mpc, diff: {frac_rs_diff:.4f}%)")
print(f"BBN Expansion Deviation (z = 10^9)   : {frac_H_diff_bbn:.6f}% (Max across all z: {max_H_diff_all:.4f}%)")
print(f"Baseline TRUE Wishart Delta chi^2    : {chi2_w0:.4f} (UNDILUTED TT+TE+EE)")
print(f"Baseline TRUE Binned Knox chi^2      : {knox_b0['chi2_total']:.4f} across {knox_b0['n_bins']} bins (RMS={knox_b0['rms_sigma_per_bin']:.4f} sigma/bin)")
print(f"Saved Target & BBN Figure to         : {REF_FILE}, {bbn_png}")
print("=" * 78)
