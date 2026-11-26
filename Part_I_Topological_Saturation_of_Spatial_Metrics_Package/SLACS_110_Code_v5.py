import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.integrate import quad

# ==============================================================================
# 1. UNIFIED COSMOLOGY & USER CONTROLS (v5)
#    Synced identically to Paper 2 (LambdaTSM Clean Covariant Baseline):
#      H0 = 67.36 km/s/Mpc, omega_b = 0.02237,
#      R_{4/3}^o = (16/3) * (1 + 1/(16*pi^2)) = 5.36710706...,
#      omega_phantom = R_{4/3}^o * omega_b = 0.120062,
#      Omega_m0 = (omega_b + omega_phantom) / h^2 = 0.31391,
#      a_0 = c * H0 / R_{4/3}^o = 1.21936e-10 m/s^2 (0 free parameters!)
# ==============================================================================
SLACS_ORIG_XLSX = "SLACS_ORIGN_DB.xlsx"   # Primary source workbook (Auger+2009 & Shu+2017)

# Published catalog cosmology (Auger+2009 SLACS IX & Shu+2017 S4TM)
H0_CATALOG = 70.0
OM0_CATALOG = 0.30

# Active Unified Cosmology (Paper 2 Planck Baseline)
OMEGA_B = 0.02237
H0_KM_S_MPC = 67.36
H_RED = H0_KM_S_MPC / 100.0
R_PARTITION = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * math.pi**2))   # 5.36710706...
OMEGA_PHANTOM = R_PARTITION * OMEGA_B                            # 0.12006219...
OM0 = (OMEGA_B + OMEGA_PHANTOM) / (H_RED**2)                     # 0.3139097...

RESCALE_MSTAR_EXACT_DL = True            # Exact luminosity scaling: M_* *= [d_L(z; H0,Om0) / d_L(z; 70,0.30)]^2
ACTIVE_PROFILE = "jaffe"                 # Options: "hernquist" (g=1), "dehnen15" (g=1.5), "jaffe" (g=2)
W_IMF_IX = 1.0                           # 1.0 = Salpeter (SLACS IX), 0.5 = Midpoint, 0.0 = Chabrier
W_IMF_S4 = 1.0                           # 1.0 = Uniform Salpeter, 0.5 = Midpoint/Kroupa (S4TM)

# Exact CODATA / IAU Physical Constants
G = 6.67430e-11                          # m^3 kg^-1 s^-2
M_SUN = 1.98847e30                       # kg
KPC_TO_M = 3.08567758e19                 # m
MPC_TO_M = 3.08567758e22                 # m
C_LIGHT = 299792458.0                    # m/s
ARCSEC_TO_RAD = math.pi / (180.0 * 3600.0)

H0_HZ = (H0_KM_S_MPC * 1e3) / MPC_TO_M
A0_UNIFIED = (C_LIGHT * H0_HZ) / R_PARTITION          # 1.21936e-10 m/s^2
A0_LEGACY_67 = (C_LIGHT * H0_HZ) / (2.0 * math.pi)    # 1.04158e-10 m/s^2 (c*H0/2pi sensitivity check)
A0_LEGACY_73 = (C_LIGHT * (73.0e3 / MPC_TO_M)) / (2.0 * math.pi)  # 1.12879e-10 m/s^2

# Exact Dehnen (1993) 2D projected half-mass radius ratios R_eff / a_s and mass exponents m = 3 - gamma
DEHNEN_PROFILES = {
    "hernquist": {"label": "Hernquist (gamma=1.0)", "plot_label": r"Hernquist ($\gamma=1.0$)", "m": 2.0, "re_over_as": 1.815271},
    "dehnen15":  {"label": "Dehnen (gamma=1.5)",    "plot_label": r"Dehnen ($\gamma=1.5$)",    "m": 1.5, "re_over_as": 1.275663},
    "jaffe":     {"label": "Jaffe (gamma=2.0)",     "plot_label": r"Jaffe ($\gamma=2.0$)",     "m": 1.0, "re_over_as": 0.744728},
}


def trapz_compat(y, x):
    """Exact trapezoidal integration compatible with BOTH NumPy 1.x and NumPy 2.x+."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    return float(np.sum(0.5 * (y[:-1] + y[1:]) * np.diff(x)))


def ang_diam_dist_kpc(z, H0=H0_KM_S_MPC, Om0=OM0):
    c_km = C_LIGHT / 1000.0
    val, _ = quad(lambda x: 1.0 / math.sqrt(Om0 * (1.0 + x)**3 + (1.0 - Om0)), 0.0, float(z))
    return (c_km / H0) * val / (1.0 + float(z)) * 1000.0


def parse_val(val):
    if pd.isna(val):
        return np.nan
    s = str(val).strip()
    if '⋯' in s or s == '' or s.lower() == 'nan':
        return np.nan
    if '±' in s:
        s = s.split('±')[0].strip()
    s = s.replace('$\\dagger$', '').replace('$\\ddagger$', '').replace('$ast$', '').strip()
    try:
        return float(s)
    except ValueError:
        return np.nan


def solve_exact_los_RE_kpc(M_kg, R_eff_kpc, D_L_kpc, DLS_over_DS, profile_key="jaffe", use_nu=True, a_scale=A0_UNIFIED):
    """
    Exact relativistic line-of-sight geodesic lensing integral for a spherical Dehnen (1993)
    stellar profile: M(r) = M_* * (r / (r + a_s))^m, g_bar(r) = G M_* r^(m-2) / (r + a_s)^m.
    Using the exact hyperbolic substitution l = b*sinh(u) => (b/r)*dl = b*du:
        alpha(b) = (4 / c^2) * b * integral_0^inf g_eff(b * cosh(u)) du
    Solves the thin-lens equation R_E = D_d * (D_ds / D_s) * alpha(R_E) by bisection.
    """
    prof = DEHNEN_PROFILES[profile_key]
    m_exp = prof["m"]
    a_s_m = (max(float(R_eff_kpc), 0.3) / prof["re_over_as"]) * KPC_TO_M

    u_grid = np.linspace(0.0, 12.5, 320)
    cosh_u = np.cosh(u_grid)

    def alpha_rad(b_kpc):
        b_m = b_kpc * KPC_TO_M
        r_m = b_m * cosh_u
        g_bar = G * M_kg * (r_m**(m_exp - 2.0)) / ((r_m + a_s_m)**m_exp)
        if use_nu:
            x = np.sqrt(np.maximum(g_bar / a_scale, 1e-30))
            nu = np.where(x < 50.0, 1.0 / (1.0 - np.exp(-x)), 1.0)
            g_eff = g_bar * nu
        else:
            g_eff = g_bar
        return (4.0 / (C_LIGHT**2)) * b_m * trapz_compat(g_eff, u_grid)

    lo, hi = 0.05, 25.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if (D_L_kpc * DLS_over_DS * alpha_rad(mid)) > mid:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def r2_1to1(y_obs, y_pred):
    return float(1.0 - np.sum((y_obs - y_pred)**2) / np.sum((y_obs - np.mean(y_obs))**2))


# ==============================================================================
# 2. LOAD ALL 110 GRADE-A LENSES DIRECTLY FROM SLACS_ORIGN_DB.xlsx
# ==============================================================================
if not os.path.exists(SLACS_ORIG_XLSX):
    raise FileNotFoundError(f"Required original workbook not found: {SLACS_ORIG_XLSX}")

t3 = pd.read_excel(SLACS_ORIG_XLSX, sheet_name='VizieR_SLACS_IX_table3').iloc[1:].reset_index(drop=True)
t4 = pd.read_excel(SLACS_ORIG_XLSX, sheet_name='table4').iloc[1:].reset_index(drop=True)
s1 = pd.read_excel(SLACS_ORIG_XLSX, sheet_name='VizieR_SLACS_XIII_S4TM_table1').iloc[2:].reset_index(drop=True)
s2 = pd.read_excel(SLACS_ORIG_XLSX, sheet_name='table2').iloc[2:].reset_index(drop=True)

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
    if np.isnan(re_arcsec):
        re_arcsec = parse_val(r['re,V'])
    sigma = parse_val(r['σa\u200bp'])
    if np.isnan(zL) or np.isnan(zS) or np.isnan(rEin_paper) or np.isnan(logM_chab):
        continue
    DL_cat = ang_diam_dist_kpc(zL, H0=H0_CATALOG, Om0=OM0_CATALOG)
    theta_E_arcsec = (rEin_paper / DL_cat) / ARCSEC_TO_RAD
    rows.append({
        'lens_id': r['Name_clean'],
        'Survey': 'SLACS_IX (Auger+2009)',
        'z_lens': zL,
        'z_source': zS,
        'theta_E_obs_arcsec': round(theta_E_arcsec, 4),
        'R_eff_arcsec': re_arcsec,
        'sigma_obs_km_s': sigma,
        'log10_Mstar_Chabrier_H70': logM_chab,
        'log10_Mstar_Salpeter_H70': logM_salp,
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
        'log10_Mstar_Chabrier_H70': parse_val(r['log10\u2061(M∗Chab/M⊙)']),
        'log10_Mstar_Salpeter_H70': np.nan,
    })

df = pd.DataFrame(rows)
df.to_csv('SLACS_DB_v5.csv', index=False)

# ==============================================================================
# 3. EXACT UNIFIED COSMOLOGY DISTANCES, RESCALED M_* & EXACT LOS INTEGRALS
# ==============================================================================
salp_offset = round(float((df['log10_Mstar_Salpeter_H70'] - df['log10_Mstar_Chabrier_H70']).dropna().mean()), 4)

DL_list, DS_list, DLS_list, rD_list, mscale_dex_list = [], [], [], [], []
for _, row in df.iterrows():
    zL, zS = float(row['z_lens']), float(row['z_source'])
    DL_cat = ang_diam_dist_kpc(zL, H0=H0_CATALOG, Om0=OM0_CATALOG)
    DL = ang_diam_dist_kpc(zL, H0=H0_KM_S_MPC, Om0=OM0)
    DS = ang_diam_dist_kpc(zS, H0=H0_KM_S_MPC, Om0=OM0)
    DLS = (DS * (1.0 + zS) - DL * (1.0 + zL)) / (1.0 + zS)
    m_fac = (DL / DL_cat)**2 if RESCALE_MSTAR_EXACT_DL else 1.0
    DL_list.append(round(DL, 2))
    DS_list.append(round(DS, 2))
    DLS_list.append(round(DLS, 2))
    rD_list.append(round(DLS / DS, 4))
    mscale_dex_list.append(round(math.log10(m_fac), 5))

df['D_L_kpc'] = DL_list
df['D_S_kpc'] = DS_list
df['D_LS_kpc'] = DLS_list
df['D_LS_over_D_S'] = rD_list
df['Mstar_Cosmo_Shift_dex'] = mscale_dex_list
df['R_eff_kpc'] = (df['D_L_kpc'] * df['R_eff_arcsec'] * ARCSEC_TO_RAD).round(3)

# Consistent exact luminosity-distance rescaling: M_*(H0,Om0) = M_*(70,0.30) * (D_L / D_L,cat)^2
df['log10_Mstar_Active_H70'] = np.where(
    df['log10_Mstar_Salpeter_H70'].notna(),
    df['log10_Mstar_Chabrier_H70'] + W_IMF_IX * (df['log10_Mstar_Salpeter_H70'] - df['log10_Mstar_Chabrier_H70']),
    df['log10_Mstar_Chabrier_H70'] + W_IMF_S4 * salp_offset
).round(4)
df['log10_Mstar_Active'] = (df['log10_Mstar_Active_H70'] + df['Mstar_Cosmo_Shift_dex']).round(4)

# Also compute midpoint S4TM IMF (w_S4 = 0.5) column for transparent comparison
df['log10_Mstar_MidS4'] = np.where(
    df['log10_Mstar_Salpeter_H70'].notna(),
    df['log10_Mstar_Salpeter_H70'] + df['Mstar_Cosmo_Shift_dex'],
    df['log10_Mstar_Chabrier_H70'] + 0.5 * salp_offset + df['Mstar_Cosmo_Shift_dex']
).round(4)

df['M_baryonic_Msun'] = 10.0**df['log10_Mstar_Active']
df['M_baryonic_kg'] = df['M_baryonic_Msun'] * M_SUN
df['M_midS4_kg'] = (10.0**df['log10_Mstar_MidS4']) * M_SUN
df['R_Ein_obs_kpc'] = (df['D_L_kpc'] * df['theta_E_obs_arcsec'] * ARCSEC_TO_RAD).round(4)

for pkey in ["hernquist", "dehnen15", "jaffe"]:
    gr_vals, nu_vals = [], []
    for _, row in df.iterrows():
        gr_vals.append(round(solve_exact_los_RE_kpc(row['M_baryonic_kg'], row['R_eff_kpc'], row['D_L_kpc'], row['D_LS_over_D_S'], pkey, False, A0_UNIFIED), 4))
        nu_vals.append(round(solve_exact_los_RE_kpc(row['M_baryonic_kg'], row['R_eff_kpc'], row['D_L_kpc'], row['D_LS_over_D_S'], pkey, True,  A0_UNIFIED), 4))
    df[f'R_Ein_GR_{pkey}_kpc'] = gr_vals
    df[f'R_Ein_LOS_{pkey}_kpc'] = nu_vals

# Midpoint S4TM IMF (w_S4=0.5) & Legacy c*H0/(2*pi) sensitivity columns for Active Profile
gr_mid, nu_mid, nu_2pi67 = [], [], []
for _, row in df.iterrows():
    gr_mid.append(round(solve_exact_los_RE_kpc(row['M_midS4_kg'], row['R_eff_kpc'], row['D_L_kpc'], row['D_LS_over_D_S'], ACTIVE_PROFILE, False, A0_UNIFIED), 4))
    nu_mid.append(round(solve_exact_los_RE_kpc(row['M_midS4_kg'], row['R_eff_kpc'], row['D_L_kpc'], row['D_LS_over_D_S'], ACTIVE_PROFILE, True,  A0_UNIFIED), 4))
    nu_2pi67.append(round(solve_exact_los_RE_kpc(row['M_baryonic_kg'], row['R_eff_kpc'], row['D_L_kpc'], row['D_LS_over_D_S'], ACTIVE_PROFILE, True,  A0_LEGACY_67), 4))

df['R_Ein_GR_MidS4_kpc'] = gr_mid
df['R_Ein_LOS_MidS4_kpc'] = nu_mid
df['R_Ein_LOS_Legacy2pi67_kpc'] = nu_2pi67
df['R_Ein_GR_Active_kpc'] = df[f'R_Ein_GR_{ACTIVE_PROFILE}_kpc']
df['R_Ein_LOS_Active_kpc'] = df[f'R_Ein_LOS_{ACTIVE_PROFILE}_kpc']
df['Ratio_Obs_to_LOS'] = (df['R_Ein_obs_kpc'] / df['R_Ein_LOS_Active_kpc']).round(4)

df = df.sort_values(by='log10_Mstar_Active').reset_index(drop=True)
df.to_csv("SLACS_Phantom_Metric_Lensing_v5.csv", index=False)

# ==============================================================================
# 4. AUDIT SUMMARY TABLE + READY-TO-PASTE LATEX TABLE 2 BODY
# ==============================================================================
ix_mask = df['Survey'].str.startswith('SLACS_IX').values
s4_mask = ~ix_mask
all_mask = np.ones(len(df), dtype=bool)
obs_all = df['R_Ein_obs_kpc'].values

print("\n" + "=" * 112)
print(f" SLACS STRONG LENSING EXACT LOS AUDIT v5 — UNIFIED COSMOLOGY (Paper 2 Sync)")
print(f" H0 = {H0_KM_S_MPC:.2f} km/s/Mpc | Om_m0 = {OM0:.5f} | R_{{4/3}}^o = {R_PARTITION:.6f} | a_0 = c*H0/R = {A0_UNIFIED:.5e} m/s^2")
print(f" Exact M_* rescaling: mean shift = {df['Mstar_Cosmo_Shift_dex'].mean():+.4f} dex | Active Profile = {DEHNEN_PROFILES[ACTIVE_PROFILE]['label']}")
print("=" * 112)

latex_rows = []
for pkey in ["hernquist", "dehnen15", "jaffe"]:
    gr_p = df[f'R_Ein_GR_{pkey}_kpc'].values
    nu_p = df[f'R_Ein_LOS_{pkey}_kpc'].values
    print(f"[{DEHNEN_PROFILES[pkey]['label']}]")
    for lbl, mask in [("SLACS IX (N=70)", ix_mask), ("Combined (N=110)", all_mask), ("S4TM (N=40)     ", s4_mask)]:
        y = obs_all[mask]; g = gr_p[mask]; n = nu_p[mask]
        _, _, rl, pl, _ = linregress(n, y)
        rat = y / n
        par = 100.0 * (1.0 - abs(y.sum() - n.sum()) / y.sum())
        r2g, r2n = r2_1to1(y, g), r2_1to1(y, n)
        print(f"  {lbl} | GR R2_1:1={r2g:+.4f} | GR+nu R2_1:1={r2n:+.4f} | lin R2={rl**2:.4f} (p={pl:.1e}) | Ratio={rat.mean():.3f}+/-{rat.std(ddof=1):.3f} (med {np.median(rat):.3f}) | Sum {y.sum():.1f}/{n.sum():.1f} ({par:.1f}%)")
        if "S4TM" not in lbl:
            latex_rows.append(
                f"{DEHNEN_PROFILES[pkey]['label']} & {lbl.strip()} & Salpeter & ${r2g:.4f}$ & ${r2n:.4f}$ & ${rl**2:.4f}$ & ${rat.mean():.3f} \\pm {rat.std(ddof=1):.3f}$ & {par:.1f}\\% \\\\"
            )

print(f"\n[IMF Sensitivity: {DEHNEN_PROFILES[ACTIVE_PROFILE]['label']} with Midpoint S4TM IMF (w_S4 = 0.5)]")
for lbl, mask in [("S4TM (N=40, w=0.5)    ", s4_mask), ("Combined (N=110, w=0.5)", all_mask)]:
    y = obs_all[mask]; g = df['R_Ein_GR_MidS4_kpc'].values[mask]; n = df['R_Ein_LOS_MidS4_kpc'].values[mask]
    _, _, rl, pl, _ = linregress(n, y)
    rat = y / n
    par = 100.0 * (1.0 - abs(y.sum() - n.sum()) / y.sum())
    r2g, r2n = r2_1to1(y, g), r2_1to1(y, n)
    print(f"  {lbl} | GR R2_1:1={r2g:+.4f} | GR+nu R2_1:1={r2n:+.4f} | lin R2={rl**2:.4f} (p={pl:.1e}) | Ratio={rat.mean():.3f}+/-{rat.std(ddof=1):.3f} (med {np.median(rat):.3f}) | Sum {y.sum():.1f}/{n.sum():.1f} ({par:.1f}%)")
    if "Combined" in lbl:
        latex_rows.append(
            f"Jaffe ($\\gamma=2.0$) & Combined ($N=110$) & Salp (IX) / Mid (S4TM) & ${r2g:.4f}$ & ${r2n:.4f}$ & ${rl**2:.4f}$ & ${rat.mean():.3f} \\pm {rat.std(ddof=1):.3f}$ & {par:.1f}\\% \\\\"
        )

print(f"\n[Threshold Sensitivity: Legacy c*H0/(2*pi) = {A0_LEGACY_67:.4e} m/s^2 at same H0={H0_KM_S_MPC:.2f}, Jaffe]")
for lbl, mask in [("SLACS IX (N=70) ", ix_mask), ("Combined (N=110)", all_mask)]:
    y = obs_all[mask]; n = df['R_Ein_LOS_Legacy2pi67_kpc'].values[mask]
    rat = y / n
    par = 100.0 * (1.0 - abs(y.sum() - n.sum()) / y.sum())
    print(f"  {lbl} | GR+nu(2pi) R2_1:1={r2_1to1(y, n):+.4f} | Ratio={rat.mean():.3f}+/-{rat.std(ddof=1):.3f} (med {np.median(rat):.3f}) | Sum {y.sum():.1f}/{n.sum():.1f} ({par:.1f}%)")

print("\n--- READY-TO-PASTE LATEX TABLE 2 BODY ---")
for lr in latex_rows:
    print(lr)
print("=" * 112 + "\n")

# ==============================================================================
# 5. PUBLICATION-GRADE 2-PANEL FIGURE -> SLACS_110_Figure_v5.png
# ==============================================================================
df_ix = df[ix_mask]
df_s4 = df[s4_mask]

obs_ix = df_ix['R_Ein_obs_kpc'].values
pred_ix = df_ix['R_Ein_LOS_Active_kpc'].values
pred_gr_ix = df_ix['R_Ein_GR_Active_kpc'].values

pred_mid_all = df['R_Ein_LOS_MidS4_kpc'].values
pred_gr_mid_all = df['R_Ein_GR_MidS4_kpc'].values
pred_salp_all = df['R_Ein_LOS_Active_kpc'].values

_, _, r_ix, p_ix, _ = linregress(pred_ix, obs_ix)
_, _, r_mid, p_mid, _ = linregress(pred_mid_all, obs_all)
r2_11_ix = r2_1to1(obs_ix, pred_ix)
r2_11_gr_ix = r2_1to1(obs_ix, pred_gr_ix)
r2_11_mid = r2_1to1(obs_all, pred_mid_all)
r2_11_gr_mid = r2_1to1(obs_all, pred_gr_mid_all)
r2_11_salp = r2_1to1(obs_all, pred_salp_all)

plt.style.use('default')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 6.2), dpi=300)
max_re = max(obs_all.max(), pred_mid_all.max(), pred_ix.max()) * 1.07

# Panel 1: SLACS IX Homogeneous Sample (N = 70)
ax1.scatter(pred_gr_ix, obs_ix, color='#94a3b8', edgecolor='#475569', s=48, alpha=0.65, marker='s',
            label=f'Pure GR Baryons Only ($R^2_{{1:1}} = {r2_11_gr_ix:.4f}$)')
ax1.scatter(pred_ix, obs_ix, color='#0f766e', edgecolor='white', s=68, alpha=0.92, zorder=4,
            label=f'Exact LOS GR + $\\nu(a_0)$ ($R^2_{{1:1}} = {r2_11_ix:.4f}$)')
ax1.plot([0, max_re], [0, max_re], color='#be123c', linestyle='--', linewidth=1.8, zorder=5, label='1:1 Identity Line')
ax1.set_title(f'SLACS IX Homogeneous Sample ($N=70$, Unified $H_0={H0_KM_S_MPC:.2f}$)\nProfile: {DEHNEN_PROFILES[ACTIVE_PROFILE]["plot_label"]} | $a_0 = cH_0/\\mathcal{{R}}_{{4/3}}^\\circ$', fontsize=10.8, fontweight='bold', pad=10)
ax1.set_xlabel('Predicted Einstein Radius $R_{E,\\mathrm{calc}}$ (kpc)', fontsize=11)
ax1.set_ylabel('HST Observed Einstein Radius $R_{E,\\mathrm{obs}} = D_d\\theta_E$ (kpc)', fontsize=11)
ax1.set_xlim(0.4, max_re); ax1.set_ylim(0.4, max_re)
ax1.grid(True, linestyle=':', alpha=0.5)
txt1 = (
    f'Unified $a_0 = cH_0/\\mathcal{{R}}_{{4/3}}^\\circ$ (0% DM)\n'
    f'1:1 $R^2_{{1:1}} = {r2_11_ix:.4f}$ (vs GR ${r2_11_gr_ix:.4f}$)\n'
    f'Linear $R^2 = {r_ix**2:.4f}$ ($p = {p_ix:.1e}$)\n'
    f'Mean $R_{{\\mathrm{{obs}}}}/R_{{\\mathrm{{calc}}}} = {np.mean(obs_ix/pred_ix):.3f}$ (med {np.median(obs_ix/pred_ix):.3f})\n'
    f'Sum: {obs_ix.sum():.1f} Obs vs {pred_ix.sum():.1f} kpc ({100*(1-abs(obs_ix.sum()-pred_ix.sum())/obs_ix.sum()):.1f}%)'
)
ax1.text(0.04, 0.96, txt1, transform=ax1.transAxes, va='top', fontsize=9.3,
         bbox=dict(boxstyle='round,pad=0.45', facecolor='#f0fdfa', edgecolor='#0f766e', alpha=0.95))
ax1.legend(loc='lower right', fontsize=9.0, framealpha=0.95)

# Panel 2: Combined N = 110 (showing Midpoint S4TM IMF with transparent Uniform Salpeter note)
ax2.scatter(df_ix['R_Ein_LOS_MidS4_kpc'], df_ix['R_Ein_obs_kpc'], color='#0f766e', edgecolor='white', s=65, alpha=0.90, zorder=4,
            label=f'SLACS IX ($N=70$, Salpeter; $R^2_{{1:1}}={r2_11_ix:.4f}$)')
ax2.scatter(df_s4['R_Ein_LOS_MidS4_kpc'], df_s4['R_Ein_obs_kpc'], color='#7c3aed', edgecolor='white', s=58, marker='D', alpha=0.85, zorder=3,
            label=f'SLACS XIII S4TM ($N=40$, Midpoint IMF)')
ax2.plot([0, max_re], [0, max_re], color='#be123c', linestyle='--', linewidth=1.8, zorder=5, label='1:1 Identity Line')
ax2.set_title(f'Combined SLACS Sample ($N=110$, Unified $H_0={H0_KM_S_MPC:.2f}$)\nZero Free Parameters ($a_0 = {A0_UNIFIED*1e10:.3f}\\times 10^{{-10}}\\,\\mathrm{{m\\,s^{{-2}}}}$)', fontsize=10.8, fontweight='bold', pad=10)
ax2.set_xlabel('Predicted Einstein Radius $R_{E,\\mathrm{calc}}$ (kpc)', fontsize=11)
ax2.set_ylabel('HST Observed Einstein Radius $R_{E,\\mathrm{obs}} = D_d\\theta_E$ (kpc)', fontsize=11)
ax2.set_xlim(0.4, max_re); ax2.set_ylim(0.4, max_re)
ax2.grid(True, linestyle=':', alpha=0.5)
txt2 = (
    f'Combined $N=110$ Audit (Salp IX / Mid S4TM)\n'
    f'1:1 $R^2_{{1:1}} = {r2_11_mid:.4f}$ (vs GR ${r2_11_gr_mid:.4f}$)\n'
    f'Linear $R^2 = {r_mid**2:.4f}$ ($p = {p_mid:.1e}$)\n'
    f'Mean $R_{{\\mathrm{{obs}}}}/R_{{\\mathrm{{calc}}}} = {np.mean(obs_all/pred_mid_all):.3f}$ | Sum Parity {100*(1-abs(obs_all.sum()-pred_mid_all.sum())/obs_all.sum()):.1f}%\n'
    f'[Uniform Salpeter $N=110$: $R^2_{{1:1}} = {r2_11_salp:.4f}$]'
)
ax2.text(0.04, 0.96, txt2, transform=ax2.transAxes, va='top', fontsize=9.1,
         bbox=dict(boxstyle='round,pad=0.45', facecolor='#f8fafc', edgecolor='#334155', alpha=0.95))
ax2.legend(loc='lower right', fontsize=8.8, framealpha=0.95)

plt.tight_layout()
plt.savefig('SLACS_110_Figure_v5.png', dpi=300, bbox_inches='tight', pad_inches=0.14)
if 'ipykernel' in sys.modules:
    plt.show()
else:
    plt.close(fig)
print("Saved SLACS_110_Figure_v5.png, SLACS_DB_v5.csv, SLACS_Phantom_Metric_Lensing_v5.csv")