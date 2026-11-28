import os
import sys
import math
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ==============================================================================
# 1. UNIFIED COSMOLOGY CONSTANTS (Synced with Paper 2 Baseline)
# ==============================================================================
LVDB_CSV = "comb_all.csv"
G = 6.67430e-11
M_SUN = 1.98847e30
C_LIGHT = 299792458.0
PC_TO_M = 3.08567758e16
KPC_TO_M = 3.08567758e19
MPC_TO_M = 3.08567758e22

H0_KM_S_MPC = 67.36
R_PARTITION = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * math.pi**2))   # 5.36710706...
H0_HZ = (H0_KM_S_MPC * 1e3) / MPC_TO_M
A0_UNIFIED = (C_LIGHT * H0_HZ) / R_PARTITION                      # 1.21936e-10 m/s^2
A0_LEGACY_67 = (C_LIGHT * H0_HZ) / (2.0 * math.pi)                # 1.04158e-10 m/s^2
A0_LEGACY_73 = (C_LIGHT * (73.0e3 / MPC_TO_M)) / (2.0 * math.pi)  # 1.12879e-10 m/s^2


def nu_rar(y):
    x = np.minimum(np.sqrt(np.maximum(y, 1e-30)), 50.0)
    return 1.0 / (1.0 - np.exp(-x))


# ==============================================================================
# 2. LOAD & FILTER LVDB CATALOG (comb_all.csv)
# ==============================================================================
if not os.path.exists(LVDB_CSV):
    raise FileNotFoundError(f"Required LVDB catalog not found: {LVDB_CSV}")

lvdb = pd.read_csv(LVDB_CSV)
df = lvdb[
    lvdb['vlos_sigma'].notna() & (lvdb['vlos_sigma'] > 0) &
    lvdb['mass_stellar'].notna() &
    lvdb['rhalf_sph_physical'].notna() & (lvdb['rhalf_sph_physical'] > 0)
].copy()

def classify(t):
    t = str(t)
    if t.startswith('dwarf_'):
        return 'Dwarf Galaxy'
    if t.startswith('gc_') and t != 'gc_ambiguous':
        return 'Globular/Star Cluster'
    return 'Other'

df['Category'] = df['table'].map(classify)
df = df[df['Category'].isin(['Dwarf Galaxy', 'Globular/Star Cluster'])].copy()
df['is_external'] = ~df['table'].isin(['dwarf_mw', 'gc_harris', 'gc_mw_new'])

# Exact Wolf (2010) 3D half-light radius r_1/2 = (4/3) * R_e
df['R_e_pc'] = df['rhalf_sph_physical']
df['r_half_m'] = (4.0 / 3.0) * df['R_e_pc'] * PC_TO_M

M_star = 10.0**df['mass_stellar']
M_HI = np.where(df['mass_HI'].notna(), 1.33 * (10.0**df['mass_HI']), 0.0)
df['has_HI'] = df['mass_HI'].notna()
df['M_bar_Msun'] = M_star + M_HI
df['g_bar'] = G * (0.5 * df['M_bar_Msun'] * M_SUN) / (df['r_half_m']**2)
df['g_obs'] = 3.0 * ((df['vlos_sigma'] * 1000.0)**2) / df['r_half_m']
df['eta_obs'] = 1.0 - df['ellipticity']
df['boost'] = df['g_obs'] / df['g_bar']

# Model predictions
df['g_Newt'] = df['g_bar']
df['g_UnifiedA0'] = df['g_bar'] * nu_rar(df['g_bar'] / A0_UNIFIED)
df['g_Legacy2pi67'] = df['g_bar'] * nu_rar(df['g_bar'] / A0_LEGACY_67)
df['g_Legacy2pi73'] = df['g_bar'] * nu_rar(df['g_bar'] / A0_LEGACY_73)
f_flat = np.sqrt(np.maximum(1.0 - df['eta_obs'].fillna(0.0)**2, 1e-12))
df['g_FlatSmooth'] = np.where(df['eta_obs'].notna(), df['g_bar'] * nu_rar(df['g_bar'] / (A0_UNIFIED * f_flat)), np.nan)

for col_g, col_s in [
    ('g_Newt', 'sigma_Newt'),
    ('g_UnifiedA0', 'sigma_UnifiedA0'),
    ('g_Legacy2pi67', 'sigma_Legacy2pi67'),
    ('g_Legacy2pi73', 'sigma_Legacy2pi73'),
    ('g_FlatSmooth', 'sigma_FlatSmooth'),
]:
    df[col_s] = np.sqrt(df['r_half_m'] * df[col_g] / 3.0) / 1000.0

df['log_res_UnifiedA0'] = np.log10(df['g_obs'] / df['g_UnifiedA0'])

# Exact algebraic inversion: f_obs = a_{0,eff,obs} / a_0
ratio_safe = np.where(df['g_obs'] > 1.001 * df['g_bar'], df['g_obs'] / np.maximum(df['g_obs'] - df['g_bar'], 1e-30), math.e)
df['f_obs'] = np.where(df['g_obs'] > 1.001 * df['g_bar'], (df['g_bar'] / (np.log(ratio_safe)**2)) / A0_UNIFIED, 0.0)

df.to_csv("LVDB_Unified_Audit_v5.csv", index=False)

# ==============================================================================
# 3. CONSOLE AUDIT REPORT
# ==============================================================================
dwf_all = df[df['Category'] == 'Dwarf Galaxy']
dsph_all = dwf_all[(dwf_all['M_V'] < -8.0) & (~dwf_all['has_HI'])]
dsph_ext = dsph_all[dsph_all['is_external']]
dwf_ext = dwf_all[dwf_all['is_external']]
dwf_iso = dwf_all[dwf_all['table'].isin(['dwarf_local_field', 'dwarf_local_field_distant'])]
gcs = df[df['Category'] == 'Globular/Star Cluster']

print("\n" + "=" * 106)
print(f" LVDB MORPHOLOGY & FLATTENING FALSIFICATION AUDIT v5 — UNIFIED COSMOLOGY (H0={H0_KM_S_MPC:.2f})")
print(f" a_0 = c*H0 / R_{{4/3}}^o = {A0_UNIFIED:.5e} m/s^2 | Confirmed Systems: N={len(df)} ({len(dwf_all)} Dwarfs, {len(gcs)} GCs)")
print("=" * 106)

def summarize_subset(s, title):
    s_e = s[s['eta_obs'].notna()]
    print(f"\n[{title} | N={len(s)}, w/ ellip={len(s_e)}]")
    for lbl, col_s in [
        ("Newtonian (3D OFF)", 'sigma_Newt'),
        ("Unified a_0 (cH0/R)", 'sigma_UnifiedA0'),
        ("Legacy cH0(73)/(2pi)", 'sigma_Legacy2pi73'),
        ("Legacy cH0(67)/(2pi)", 'sigma_Legacy2pi67'),
        ("Flattening sqrt(1-eta^2)", 'sigma_FlatSmooth'),
    ]:
        ss = s_e if col_s == 'sigma_FlatSmooth' else s
        if len(ss) == 0 or ss[col_s].isna().any():
            continue
        dl = np.log10(ss['vlos_sigma'].values / ss[col_s].values)
        print(f"  {lbl:26s} | sig RMS={math.sqrt(np.mean(dl**2)):.4f} dex | bias={np.mean(dl):+.4f} dex | med obs/pred={np.median(ss['vlos_sigma'].values/ss[col_s].values):.2f}x")
    if len(s_e) >= 5:
        r_res, p_res = spearmanr(s_e['eta_obs'], s_e['log_res_UnifiedA0'])
        r_f, p_f = spearmanr(s_e['eta_obs'], s_e['f_obs'])
        best_rms, best_ec = 1e9, 0.0
        for ec in np.linspace(0.10, 1.01, 92):
            gst = np.where(s_e['eta_obs'] < ec, s_e['g_UnifiedA0'], s_e['g_Newt'])
            sp_ = np.sqrt(s_e['r_half_m'] * gst / 3.0) / 1000.0
            rm = math.sqrt(np.mean(np.log10(s_e['vlos_sigma'] / sp_)**2))
            if rm < best_rms:
                best_rms, best_ec = rm, ec
        print(f"  1/0 Step Scan              | best eta_c={best_ec:.2f} (sig RMS={best_rms:.4f} dex; ec=1.01 => NEVER shuts off)")
        print(f"  Spearman vs eta (b/a)      | RAR-res: rho={r_res:+.3f} (p={p_res:.3f}) | inverted f_obs: rho={r_f:+.3f} (p={p_f:.3f})")

summarize_subset(dsph_all, "Gas-Free Classical dSphs (All, M_V < -8)")
summarize_subset(dsph_ext, "Gas-Free Classical dSphs (EXTERNAL ONLY, 0% MW)")
summarize_subset(dwf_ext,  "All External Dwarfs (M31 + Field, 0% MW)")
summarize_subset(dwf_iso,  "Isolated Field Dwarfs (0% MW, 0% M31)")
summarize_subset(gcs,      "Globular / Star Clusters")
print("=" * 106 + "\n")

# ==============================================================================
# 4. PUBLICATION-GRADE 3-PANEL FIGURE -> LVDB_Flattening_Figure_v5.png
# ==============================================================================
dsph_e = dsph_all[dsph_all['eta_obs'].notna()]
dsph_ext_e = dsph_ext[dsph_ext['eta_obs'].notna()]
dwf_all_e = dwf_all[dwf_all['eta_obs'].notna()]
gc_e = gcs[gcs['eta_obs'].notna()]

r_f_cl, p_f_cl = spearmanr(dsph_e['eta_obs'], dsph_e['f_obs'])

plt.style.use('default')
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.5, 5.5), dpi=300)

# Panel 1: Local Volume RAR (g_obs vs g_bar) colored by b/a
gb_grid = np.logspace(-14.5, -7.5, 300)
ax1.loglog(gb_grid, gb_grid, 'k--', lw=1.4, label='Newtonian ($g_{\\mathrm{obs}} = g_{\\mathrm{bar}}$)')
ax1.loglog(gb_grid, gb_grid * nu_rar(gb_grid / A0_UNIFIED), color='#0f766e', lw=2.2,
           label=r'Unified $a_0 = cH_0/\mathcal{R}_{4/3}^\circ$')
ax1.loglog(gb_grid, gb_grid * nu_rar(gb_grid / (A0_UNIFIED * math.sqrt(1.0 - 0.90**2))), color='#f59e0b', ls='-.', lw=1.6,
           label=r'Axis-Ratio $\sqrt{1-\eta^2}$ ($\eta = 0.90$)')
sc = ax1.scatter(dwf_all_e['g_bar'], dwf_all_e['g_obs'], c=dwf_all_e['eta_obs'], cmap='viridis',
                 vmin=0.2, vmax=1.0, s=42, edgecolor='black', lw=0.4, zorder=4, label=f'LVDB Dwarfs w/ $b/a$ ($N={len(dwf_all_e)}/{len(dwf_all)}$)')
ax1.scatter(gcs['g_bar'], gcs['g_obs'], marker='s', color='#dc2626', alpha=0.70, s=28, edgecolor='white', lw=0.35, zorder=3,
            label=f'Globular Clusters ($N={len(gcs)}$)')
cb = fig.colorbar(sc, ax=ax1, pad=0.02); cb.set_label(r'Axis Ratio $\eta_{\mathrm{obs}} = b/a$', fontsize=9.2)
ax1.set_xlim(1e-14, 1e-7); ax1.set_ylim(1e-13, 1e-7)
ax1.set_xlabel(r'$g_{\mathrm{bar}}(r_{1/2})$ [$\mathrm{m\,s^{-2}}$]', fontsize=10.5)
ax1.set_ylabel(r'$g_{\mathrm{obs}}(r_{1/2})$ [$\mathrm{m\,s^{-2}}$]', fontsize=10.5)
ax1.set_title('A) Local Volume RAR (LVDB, $N=158$)\nRound 3D Dwarfs Follow Universal $a_0$', fontsize=10.8, fontweight='bold', pad=8)
ax1.legend(loc='upper left', fontsize=8.0, framealpha=0.92); ax1.grid(True, which='major', ls=':', alpha=0.5)

# Panel 2: Direct Empirical Inversion f_obs(eta) vs b/a
eta_line = np.linspace(0.05, 0.999, 300)
ax2.axhline(1.0, color='#0f766e', lw=2.2, label=r'Universal $a_0$ ($f(\eta)=1$)')
ax2.plot(eta_line, np.sqrt(1.0 - eta_line**2), color='#be123c', ls='--', lw=2.2, label=r'Axis-Ratio $\sqrt{1-\eta^2}$ Curve')
other_dwf = dwf_all_e[~dwf_all_e.index.isin(dsph_e.index)]
ax2.scatter(other_dwf['eta_obs'], np.clip(other_dwf['f_obs'], 1e-2, 100), color='#94a3b8', alpha=0.60, s=32,
            label=f'Other LVDB Dwarfs ($N={len(other_dwf)}$)')
ax2.scatter(dsph_e[~dsph_e['is_external']]['eta_obs'], np.clip(dsph_e[~dsph_e['is_external']]['f_obs'], 1e-2, 100),
            color='#60a5fa', edgecolor='black', lw=0.4, s=55, zorder=5, label=f'MW Classical dSphs ($N={(~dsph_e["is_external"]).sum()}$)')
ax2.scatter(dsph_ext_e['eta_obs'], np.clip(dsph_ext_e['f_obs'], 1e-2, 100),
            color='#1d4ed8', edgecolor='white', lw=0.5, s=68, zorder=6, label=f'External dSphs, 0% MW ($N={len(dsph_ext_e)}$)')
ax2.scatter(gc_e['eta_obs'], np.clip(gc_e['f_obs'], 1e-2, 100), marker='s', color='#dc2626', edgecolor='black', s=46, zorder=7,
            label=f'GCs w/ Ellipticity ($N={len(gc_e)}$)')
for nm in ['Andromeda IX', 'Andromeda VII', 'Leo II', 'NGC 2419']:
    row = df[df['name'] == nm]
    if len(row) > 0 and pd.notna(row['eta_obs'].iloc[0]):
        ax2.annotate(nm, (row['eta_obs'].iloc[0], max(row['f_obs'].iloc[0], 0.012)),
                     xytext=(-44, 7), textcoords='offset points', fontsize=8.0, fontweight='bold',
                     arrowprops=dict(arrowstyle='->', lw=0.7, color='black'))
ax2.set_yscale('log'); ax2.set_xlim(0.20, 1.03); ax2.set_ylim(0.008, 120)
ax2.set_xlabel(r'Observed Axis Ratio $\eta_{\mathrm{obs}} = b/a$', fontsize=10.5)
ax2.set_ylabel(r'Inverted Threshold $f_{\mathrm{obs}} = a_{0,\mathrm{eff,obs}} / a_0$', fontsize=10.5)
ax2.set_title(f'B) Direct Inversion $f_{{\\mathrm{{obs}}}}(\\eta)$\nNo Drop at $\\eta\\to 1$ ($\\rho={r_f_cl:+.3f}, p={p_f_cl:.2f}$)', fontsize=10.8, fontweight='bold', pad=8)
ax2.legend(loc='lower left', fontsize=7.8, framealpha=0.92); ax2.grid(True, ls=':', alpha=0.5)

# Panel 3: 1/0 Step-Function Shutoff Cutoff Scan eta_c
ecs = np.linspace(0.20, 1.01, 82)
rms_cl, rms_ext, rms_all = [], [], []
for ec in ecs:
    for sub_s, lst in [(dsph_e, rms_cl), (dsph_ext_e, rms_ext), (dwf_all_e, rms_all)]:
        gst = np.where(sub_s['eta_obs'] < ec, sub_s['g_UnifiedA0'], sub_s['g_Newt'])
        sp_ = np.sqrt(sub_s['r_half_m'] * gst / 3.0) / 1000.0
        lst.append(math.sqrt(np.mean(np.log10(sub_s['vlos_sigma'] / sp_)**2)))

ax3.plot(ecs, rms_ext, color='#1d4ed8', lw=2.3, label=f'External dSphs, 0% MW ($N={len(dsph_ext_e)}$)')
ax3.plot(ecs, rms_cl, color='#0f766e', lw=2.0, ls='--', label=f'All Classical dSphs ($N={len(dsph_e)}$)')
ax3.plot(ecs, rms_all, color='#64748b', lw=1.8, ls='-.', label=f'All LVDB Dwarfs ($N={len(dwf_all_e)}$)')
rms_smooth_ext = math.sqrt(np.mean(np.log10(dsph_ext_e['vlos_sigma'] / dsph_ext_e['sigma_FlatSmooth'])**2))
ax3.axhline(rms_smooth_ext, color='#be123c', ls=':', lw=1.8, label=r'Smooth $\sqrt{1-\eta^2}$ (External dSphs)')
ax3.set_xlabel(r'Step Cutoff $\eta_c$ (Newtonian if $b/a \geq \eta_c$)', fontsize=10.5)
ax3.set_ylabel(r'$\sigma_{\mathrm{los}}$ RMS $\log_{10}(\sigma_{\mathrm{obs}}/\sigma_{\mathrm{pred}})$ [dex]', fontsize=10.5)
ax3.set_title('C) Step Cutoff Scan $\\Theta(\\eta_c - \\eta)$\nMonotonic Drop Until $\\eta_c > 1$ (No Cutoff)', fontsize=10.8, fontweight='bold', pad=8)
ax3.legend(loc='upper right', fontsize=8.0, framealpha=0.92); ax3.grid(True, ls=':', alpha=0.5)

plt.tight_layout()
plt.savefig("LVDB_Flattening_Figure_v5.png", dpi=300, bbox_inches='tight', pad_inches=0.12)
if 'ipykernel' in sys.modules:
    plt.show()
else:
    plt.close(fig)
print("Saved LVDB_Flattening_Figure_v5.png and LVDB_Unified_Audit_v5.csv")