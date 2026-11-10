import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.integrate import quad

# ==============================================================================
# 1. PHANTOM METRIC CONSTANTS (ZERO FREE PARAMETERS)
# ==============================================================================
G = 6.67430e-11             # m^3 kg^-1 s^-2
M_sun = 1.98847e30          # kg
kpc_to_m = 3.08567758e19    # m
c = 299792458.0             # m/s
arcsec_to_rad = np.pi / (180.0 * 3600.0)

H0_km_s_Mpc = 73.0
H0_Hz = (H0_km_s_Mpc * 1e3) / 3.08567758e22
a0 = (c * H0_Hz) / (2.0 * np.pi)  # 1.1288e-10 m/s^2

def ang_diam_dist_kpc(z, H0=73.0, Om0=0.3):
    c_km = 299792.458
    val, _ = quad(lambda x: 1.0 / np.sqrt(Om0 * (1.0 + x)**3 + (1.0 - Om0)), 0, z)
    return (c_km / H0) * val / (1.0 + z) * 1000.0

def parse_val(val):
    if pd.isna(val): return np.nan
    s = str(val).strip()
    if '⋯' in s or s == '' or s.lower() == 'nan': return np.nan
    if '±' in s: s = s.split('±')[0].strip()
    s = s.replace('$\\dagger$', '').replace('$\\ddagger$', '').replace('$ast$', '').strip()
    try: return float(s)
    except: return np.nan

# ==============================================================================
# 2. LOAD UNIFIED 110-GALAXY DATABASE FROM SLACS_DB.csv (OR SLACS_ORIGN_DB.xlsx)
# ==============================================================================
need_build = True
if os.path.exists('SLACS_DB.csv'):
    df_check = pd.read_csv('SLACS_DB.csv')
    if len(df_check) >= 110 and {'lens_id', 'z_lens', 'z_source', 'theta_E_obs_arcsec'}.issubset(df_check.columns):
        df = df_check.copy()
        need_build = False

if need_build and os.path.exists('SLACS_ORIGN_DB.xlsx'):
    t3 = pd.read_excel('SLACS_ORIGN_DB.xlsx', sheet_name='VizieR_SLACS_IX_table3').iloc[1:].reset_index(drop=True)
    t4 = pd.read_excel('SLACS_ORIGN_DB.xlsx', sheet_name='table4').iloc[1:].reset_index(drop=True)
    s1 = pd.read_excel('SLACS_ORIGN_DB.xlsx', sheet_name='VizieR_SLACS_XIII_S4TM_table1').iloc[2:].reset_index(drop=True)
    s2 = pd.read_excel('SLACS_ORIGN_DB.xlsx', sheet_name='table2').iloc[2:].reset_index(drop=True)

    t3['Name_clean'] = t3['Name'].astype(str).str.strip().str.replace('−', '-')
    t4['Name_clean'] = t4['Name'].astype(str).str.strip().str.replace('−', '-')
    slacs9 = pd.merge(t3, t4, on='Name_clean')

    s1['Target_clean'] = s1['Target'].astype(str).str.strip().str.replace('−', '-')
    s2['Target_clean'] = s2['Target'].astype(str).str.strip().str.replace('−', '-')
    s4tm = pd.merge(s1, s2, on='Target_clean')

    rows = []
    for _, r in slacs9.iterrows():
        zL = parse_val(r['zl\u200be\u200bn\u200bs'])
        zS = parse_val(r['zs\u200br\u200bc'])
        rEin_paper = parse_val(r['rEin'])
        logM_chab = parse_val(r['log[M∗Chab/M⊙]'])
        logM_salp = parse_val(r['log[M∗Salp/M⊙]'])
        re_arcsec = parse_val(r['re,I'])
        if np.isnan(re_arcsec): re_arcsec = parse_val(r['re,V'])
        sigma = parse_val(r['σa\u200bp'])
        if np.isnan(zL) or np.isnan(zS) or np.isnan(rEin_paper) or np.isnan(logM_chab):
            continue
        DL_70 = ang_diam_dist_kpc(zL, H0=70.0, Om0=0.3)
        theta_E_arcsec = (rEin_paper / DL_70) / arcsec_to_rad
        rows.append({
            'lens_id': r['Name_clean'],
            'Survey': 'SLACS_IX (Auger+2009)',
            'z_lens': zL,
            'z_source': zS,
            'theta_E_obs_arcsec': round(theta_E_arcsec, 3),
            'R_eff_arcsec': re_arcsec,
            'sigma_obs_km_s': sigma,
            'log10_Mstar_Chabrier': logM_chab,
            'log10_Mstar_Salpeter': logM_salp
        })

    for _, r in s4tm.iterrows():
        rows.append({
            'lens_id': r['Target_clean'],
            'Survey': 'SLACS_XIII_S4TM (Shu+2017)',
            'z_lens': parse_val(r['zL']),
            'z_source': parse_val(r['zS']),
            'theta_E_obs_arcsec': parse_val(r['bSIE']),
            'R_eff_arcsec': parse_val(r['Reff']),
            'sigma_obs_km_s': parse_val(r['σSDSS']),
            'log10_Mstar_Chabrier': parse_val(r['log10\u2061(M∗Chab/M⊙)']),
            'log10_Mstar_Salpeter': np.nan
        })

    df = pd.DataFrame(rows)
    df.to_csv('SLACS_DB.csv', index=False)

# ==============================================================================
# 3. ZERO-PARAMETER PHYSICAL CONVERSIONS (MEASUREMENT B: R_E IN LENS PLANE)
# ==============================================================================
DL_list, DS_list, DLS_list, DLS_over_DS_list = [], [], [], []
for _, row in df.iterrows():
    zL, zS = float(row['z_lens']), float(row['z_source'])
    DL = ang_diam_dist_kpc(zL, H0=73.0, Om0=0.3)
    DS = ang_diam_dist_kpc(zS, H0=73.0, Om0=0.3)
    Dc_L = DL * (1.0 + zL)
    Dc_S = DS * (1.0 + zS)
    DLS = (Dc_S - Dc_L) / (1.0 + zS)
    DL_list.append(round(DL, 2))
    DS_list.append(round(DS, 2))
    DLS_list.append(round(DLS, 2))
    DLS_over_DS_list.append(round((Dc_S - Dc_L) / Dc_S, 4))

df['D_L_kpc'] = DL_list
df['D_S_kpc'] = DS_list
df['D_LS_kpc'] = DLS_list
df['D_LS_over_D_S'] = DLS_over_DS_list

df['log10_Mstar_Published_Baseline'] = np.where(
    df['log10_Mstar_Salpeter'].notna(),
    0.5 * (df['log10_Mstar_Chabrier'] + df['log10_Mstar_Salpeter']),
    df['log10_Mstar_Chabrier']
).round(4)

df['m_baryonic_10_11_msun'] = (10**(df['log10_Mstar_Published_Baseline'] - 11.0)).round(4)
df['M_baryonic_absolute_sun'] = (df['m_baryonic_10_11_msun'] * 1e11).round(0)

# Replaced Hebrew column names with standard English variables
df['M_baryonic_kg'] = df['M_baryonic_absolute_sun'] * M_sun
df['GM_m3_s2'] = df['M_baryonic_kg'] * G
df['GM_over_a0_m2'] = df['GM_m3_s2'] / a0
df['R_BB_meters'] = np.sqrt(df['GM_over_a0_m2'])
df['Calculated_R_BB_kpc'] = (df['R_BB_meters'] / kpc_to_m).round(4)

theta_obs_rad = df['theta_E_obs_arcsec'] * arcsec_to_rad
df['Observed_R_BB_kpc'] = (((c**2 * theta_obs_rad) / (4.0 * np.pi * a0)) * (1.0 / df['D_LS_over_D_S']) / kpc_to_m).round(4)
df['observed_einstein_radius_kpc'] = df['Observed_R_BB_kpc']

df['R_Ein_obs_kpc'] = (df['D_L_kpc'] * df['theta_E_obs_arcsec'] * arcsec_to_rad).round(4)
alpha_rad = (4.0 * np.pi * np.sqrt(G * df['M_baryonic_kg'] * a0)) / (c**2)
df['R_Ein_calc_kpc'] = (df['D_L_kpc'] * df['D_LS_over_D_S'] * alpha_rad).round(4)
df['Ratio_Obs_to_Calc'] = (df['R_Ein_obs_kpc'] / df['R_Ein_calc_kpc']).round(4)

df = df.sort_values(by='m_baryonic_10_11_msun').reset_index(drop=True)
df.to_excel("SLACS_110_Clean_Master_Worksheet.xlsx", index=False)
df.to_csv("SLACS_Phantom_Metric_Lensing.csv", index=False)

# ==============================================================================
# 4. GENERATE SINGLE-PANEL VALIDATION FIGURE -> SLACS_110_Figure.png
# ==============================================================================
is_s4tm = df['Survey'].str.startswith('SLACS_XIII')
df_ix = df[~is_s4tm].copy()
df_s4tm = df[is_s4tm].copy()

_, _, r_re_all, p_re_all, _ = linregress(df['R_Ein_calc_kpc'], df['R_Ein_obs_kpc'])
_, _, r_re_ix, p_re_ix, _ = linregress(df_ix['R_Ein_calc_kpc'], df_ix['R_Ein_obs_kpc'])

ss_res = np.sum((df['R_Ein_obs_kpc'] - df['R_Ein_calc_kpc'])**2)
ss_tot = np.sum((df['R_Ein_obs_kpc'] - np.mean(df['R_Ein_obs_kpc']))**2)
r2_1_1 = 1.0 - (ss_res / ss_tot)

plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(10.5, 6.8))

ax.scatter(df_ix['R_Ein_calc_kpc'], df_ix['R_Ein_obs_kpc'],
           color='#00ffcc', edgecolor='white', s=75, alpha=0.9, zorder=3,
           label=f'SLACS IX Lenses ($N=70$, Auger et al. 2009; $R^2={r_re_ix**2:.4f}$)')
ax.scatter(df_s4tm['R_Ein_calc_kpc'], df_s4tm['R_Ein_obs_kpc'],
           color='#ff00ff', edgecolor='white', s=70, marker='D', alpha=0.9, zorder=3,
           label='SLACS XIII S4TM Lenses ($N=40$, Shu et al. 2017)')

max_re = max(df['R_Ein_calc_kpc'].max(), df['R_Ein_obs_kpc'].max()) * 1.06
ax.plot([0, max_re], [0, max_re], color='#ff00ff', linestyle='--', linewidth=2.2, zorder=4,
        label='Perfect Fit Line (1:1 Identity)')

ax.set_title('The Phantom Metric: Einstein Ring Predictions in SLACS Galaxies ($N = 110$ Lenses)', fontsize=14, pad=14)
ax.set_xlabel('Phantom Metric Predicted Einstein Radius $R_{E,\\mathrm{calc}}$ (kpc)', fontsize=12)
ax.set_ylabel('HST Observed Einstein Ring Radius $R_{E,\\mathrm{obs}} = D_L\\theta_E$ (kpc)', fontsize=12)
ax.set_xlim(0.4, max_re)
ax.set_ylim(0.4, max_re)
ax.grid(True, linestyle=':', alpha=0.22, color='gray')

textstr = (
    f'Combined $R^2 = {r_re_all**2:.4f}$ ($p = {p_re_all:.1e}$)\n'
    f'1:1 Parity $R^2_{{1:1}} = {r2_1_1:.4f}$ (SLACS IX $R^2 = {r_re_ix**2:.4f}$)\n'
    f'Mean $R_{{E,\\mathrm{{obs}}}}/R_{{E,\\mathrm{{calc}}}} = {df["Ratio_Obs_to_Calc"].mean():.3f}$ (Median = {df["Ratio_Obs_to_Calc"].median():.3f})\n'
    f'Ensemble Sum: {df["R_Ein_obs_kpc"].sum():.1f} kpc Obs vs. {df["R_Ein_calc_kpc"].sum():.1f} kpc Calc\n'
    f'Dark Matter = 0% (Zero Free Parameters)'
)
ax.text(0.03, 0.95, textstr, transform=ax.transAxes, verticalalignment='top',
        color='#00ffcc', fontsize=11,
        bbox=dict(boxstyle='round', facecolor='black', alpha=0.8, edgecolor='#ff00ff'))

ax.legend(loc='lower right', fontsize=10.5)
plt.tight_layout()
plt.savefig('SLACS_110_Figure.png', dpi=300)
plt.close(fig)
print("Successfully saved single-panel figure: SLACS_110_Figure.png")