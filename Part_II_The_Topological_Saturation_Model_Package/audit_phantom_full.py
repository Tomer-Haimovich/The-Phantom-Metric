import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from classy import Class
from cmb_likelihood_utils import (
    TOPO_RATIO, L_MIN, L_MAX_TT, L_MAX_POL,
    compute_planck_noise_dl, wishart_chi2_exact, binned_knox_chi2, r2_1to1
)

PKG_DIR = os.path.dirname(os.path.abspath(__file__))

BASE_PARAMS = {
    'output': 'tCl,pCl,lCl,mPk',
    'lensing': 'yes',
    'l_max_scalars': L_MAX_TT,
    'P_k_max_1/Mpc': 5.0,
    'h': 0.6736,
    'omega_b': 0.02237,
    'A_s': 2.100e-9,
    'n_s': 0.9649,
    'tau_reio': 0.0544
}


def run_class(omega_cdm_val):
    cosmo = Class()
    params = dict(BASE_PARAMS)
    params['omega_cdm'] = omega_cdm_val
    cosmo.set(params)
    cosmo.compute()
    cls = cosmo.lensed_cl(L_MAX_TT)
    ell_full = cls['ell'][2:L_MAX_TT + 1].astype(float)
    fac_full = ell_full * (ell_full + 1.0) / (2.0 * np.pi) * (2.7255e6)**2
    tt = cls['tt'][2:L_MAX_TT + 1] * fac_full
    te = cls['te'][2:L_MAX_TT + 1] * fac_full
    ee = cls['ee'][2:L_MAX_TT + 1] * fac_full
    k_vals = np.logspace(-4, 0, 200)
    pk = np.array([cosmo.pk(k, 0.0) for k in k_vals])
    s8 = float(cosmo.sigma8())
    derived = cosmo.get_current_derived_parameters(['rs_rec', 'z_rec', '100*theta_s'])
    cosmo.struct_cleanup()
    cosmo.empty()
    return ell_full, tt, te, ee, k_vals, pk, s8, derived


# 1. LIVE RUN: Planck 2018 LCDM Benchmark
ell_full, tt_l, te_l, ee_l, k_vals, pk_l, s8_l, der_l = run_class(0.1200)

# 2. LIVE RUN: Clean Covariant LambdaTSM (0% particle CDM: omega_cdm = 1e-12)
_, tt_p, te_p, ee_p, _, pk_p, s8_p, der_p = run_class(1.0e-12)

# Slice ell in [50, 2500]
i0 = L_MIN - 2  # index 48 corresponds to ell = 50
ell_eval = ell_full[i0:]
tt_l_e, te_l_e, ee_l_e = tt_l[i0:], te_l[i0:], ee_l[i0:]
tt_p_e, te_p_e, ee_p_e = tt_p[i0:], te_p[i0:], ee_p[i0:]

n_tt_dl, n_ee_dl = compute_planck_noise_dl(ell_eval)
chi2_wishart = wishart_chi2_exact(tt_p_e, te_p_e, ee_p_e, tt_l_e, te_l_e, ee_l_e, ell_eval, n_tt_dl, n_ee_dl)
knox_bin = binned_knox_chi2(tt_p_e, te_p_e, ee_p_e, tt_l_e, te_l_e, ee_l_e, ell_eval, n_tt_dl, n_ee_dl, dl=30)

r2_tt = r2_1to1(tt_p_e, tt_l_e)
pol_m = ell_eval <= L_MAX_POL
r2_te = r2_1to1(te_p_e[pol_m], te_l_e[pol_m])
r2_ee = r2_1to1(ee_p_e[pol_m], ee_l_e[pol_m])

max_idx = int(np.argmax(np.abs(tt_p_e - tt_l_e)))
max_ell = int(ell_eval[max_idx])
max_diff_uk2 = float(tt_p_e[max_idx] - tt_l_e[max_idx])
max_diff_pct = float(100.0 * max_diff_uk2 / tt_l_e[max_idx])
omega_ph_derived = TOPO_RATIO * BASE_PARAMS['omega_b']

print("\n" + "=" * 94)
print(" LAMBDATSM CLEAN COVARIANT — HONEST FULL DIAGNOSTIC AUDIT (CLASS v3.4.0, l = 50..2500)")
print("=" * 94)
print(f"{'Observable / Cosmological Metric':<44} | {'Planck 2018 (LCDM)':<18} | {'LambdaTSM (0% CDM)'}")
print("-" * 94)
print(f"{'Background Density (omega_cdm / phantom)':<44} | {'0.120000 (Fitted)':<18} | {omega_ph_derived:.6f} (Ratio Ansatz)")
print(f"{'CMB Temperature (D_l^TT) 1:1 R^2':<44} | {'1.000000':<18} | {r2_tt:.6f}")
print(f"{'CMB Cross-Polarization (D_l^TE) 1:1 R^2':<44} | {'1.000000':<18} | {r2_te:.6f}")
print(f"{'CMB E-Mode Polarization (D_l^EE) 1:1 R^2':<44} | {'1.000000':<18} | {r2_ee:.6f}")
print(f"{'Max |Delta D_l^TT| (at ell=' + str(max_ell) + ')':<44} | {'0.00 uK^2':<18} | {max_diff_uk2:+.2f} uK^2 ({max_diff_pct:+.3f}%)")
print(f"{'TRUE Undiluted 2x2 Wishart Delta chi^2':<44} | {'0.000':<18} | {chi2_wishart:.4f} (TT+TE+EE)")
print(f"{'TRUE Binned Knox chi^2 (dl=30, ' + str(knox_bin['n_bins']) + ' bins)':<44} | {'0.000':<18} | {knox_bin['chi2_total']:.4f} (RMS = {knox_bin['rms_sigma_per_bin']:.4f} sigma/bin)")
print(f"{'Linear Matter Clustering sigma_8 (z=0)':<44} | {s8_l:.5f}{'':<11} | {s8_p:.5f}")
print(f"{'Sound Horizon at Recombination r_s(z_*)':<44} | {der_l['rs_rec']:.4f} Mpc{'':<7} | {der_p['rs_rec']:.4f} Mpc (z_*={der_p['z_rec']:.2f})")
print(f"{'Acoustic Angle 100 * theta_s':<44} | {der_l['100*theta_s']:.5f}{'':<11} | {der_p['100*theta_s']:.5f}")
print("=" * 94 + "\n")

# Publication-grade 4-panel figure
fig, axs = plt.subplots(2, 2, figsize=(15, 10.5), dpi=300)

axs[0, 0].plot(ell_full, tt_l, 'k--', lw=2.0, label='Planck 2018 $\\Lambda$CDM')
axs[0, 0].plot(ell_full, tt_p, color='#0f766e', lw=1.6, alpha=0.9, label=f'$\\Lambda$TSM Clean ($R^2={r2_tt:.6f}$, $\\Delta\\chi^2_W={chi2_wishart:.2f}$)')
axs[0, 0].set_title('CMB Temperature Power Spectrum ($D_\\ell^{TT}$)', fontweight='bold', pad=8)
axs[0, 0].set_xlabel('Multipole $\\ell$'); axs[0, 0].set_ylabel('$D_\\ell^{TT}$ [$\\mu\\mathrm{K}^2$]')
axs[0, 0].set_xlim(2, 2500); axs[0, 0].legend(); axs[0, 0].grid(ls=':', alpha=0.5)

axs[0, 1].plot(ell_full[:1999], te_l[:1999], 'k--', lw=2.0, label='Planck 2018 $\\Lambda$CDM')
axs[0, 1].plot(ell_full[:1999], te_p[:1999], color='#1d4ed8', lw=1.5, alpha=0.9, label=f'$\\Lambda$TSM Clean ($R^2={r2_te:.6f}$)')
axs[0, 1].set_title('CMB Cross-Polarization Spectrum ($D_\\ell^{TE}$)', fontweight='bold', pad=8)
axs[0, 1].set_xlabel('Multipole $\\ell$'); axs[0, 1].set_ylabel('$D_\\ell^{TE}$ [$\\mu\\mathrm{K}^2$]')
axs[0, 1].set_xlim(2, 2000); axs[0, 1].legend(); axs[0, 1].grid(ls=':', alpha=0.5)

axs[1, 0].plot(ell_full[:1999], ee_l[:1999], 'k--', lw=2.0, label='Planck 2018 $\\Lambda$CDM')
axs[1, 0].plot(ell_full[:1999], ee_p[:1999], color='#be123c', lw=1.5, alpha=0.9, label=f'$\\Lambda$TSM Clean ($R^2={r2_ee:.6f}$)')
axs[1, 0].set_title('CMB E-Mode Polarization Spectrum ($D_\\ell^{EE}$)', fontweight='bold', pad=8)
axs[1, 0].set_xlabel('Multipole $\\ell$'); axs[1, 0].set_ylabel('$D_\\ell^{EE}$ [$\\mu\\mathrm{K}^2$]')
axs[1, 0].set_xlim(2, 2000); axs[1, 0].legend(); axs[1, 0].grid(ls=':', alpha=0.5)

axs[1, 1].loglog(k_vals, pk_l, 'k--', lw=2.0, label=f'Planck $\\Lambda$CDM ($\\sigma_8={s8_l:.4f}$)')
axs[1, 1].loglog(k_vals, pk_p, color='#7c3aed', lw=1.6, alpha=0.9, label=f'$\\Lambda$TSM Clean ($\\sigma_8={s8_p:.4f}$)')
axs[1, 1].set_title('Linear Matter Power Spectrum $P(k)$ at $z=0$', fontweight='bold', pad=8)
axs[1, 1].set_xlabel('Wavenumber $k$ [$\\mathrm{Mpc}^{-1}$]'); axs[1, 1].set_ylabel('$P(k)$ [$\\mathrm{Mpc}^3$]')
axs[1, 1].set_xlim(1e-4, 1e0); axs[1, 1].legend(); axs[1, 1].grid(ls=':', alpha=0.5, which='both')

plt.tight_layout()
out_png = os.path.join(PKG_DIR, "phantom_full_diagnostic.png")
plt.savefig(out_png, dpi=300, bbox_inches='tight', pad_inches=0.14)
plt.close(fig)
print(f"Saved clean 4-panel diagnostic figure to: {out_png}")
