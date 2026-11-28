import os
import sys
import io
import math
import zipfile
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
import matplotlib
if 'ipykernel' not in sys.modules:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ==============================================================================
# 1. UNIFIED COSMOLOGY & CONTROLS (v5 — Synced with Paper 2 Baseline)
# ==============================================================================
ROTMOD_ZIP = "Rotmod_LTG.zip"                  # Raw SPARC radial rotation curves archive
SPARC_META_CANDIDATES = ["SPARC_Lelli2016c.txt", "SPARC_Lelli2016c.mrt"]

G = 6.67430e-11
M_SUN = 1.98847e30
KPC_TO_M = 3.08567758e19
MPC_TO_M = 3.08567758e22
C_LIGHT = 299792458.0

H0_KM_S_MPC = 67.36
H0_SPARC_HFLOW = 73.0
R_PARTITION = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * math.pi**2))   # 5.36710706...
PHI_GOLDEN = 0.5 * (1.0 + math.sqrt(5.0))                         # 1.61803399...

H0_HZ = (H0_KM_S_MPC * 1e3) / MPC_TO_M
A0_UNIFIED = (C_LIGHT * H0_HZ) / R_PARTITION                      # 1.21936e-10 m/s^2
A0_LEGACY_67 = (C_LIGHT * H0_HZ) / (2.0 * math.pi)                # 1.04158e-10 m/s^2 (c*H0/2pi)
A0_LEGACY_73 = (C_LIGHT * (73.0e3 / MPC_TO_M)) / (2.0 * math.pi)  # 1.12879e-10 m/s^2 (v1-v4 legacy)
A0_PHI_REJECT = A0_UNIFIED * PHI_GOLDEN                           # Rejected double-saturation test

YD_STD = 0.50
YB_STD = 0.70


def nu_rar(y):
    x = np.minimum(np.sqrt(np.maximum(y, 1e-30)), 50.0)
    return 1.0 / (1.0 - np.exp(-x))


def r2_1to1(y_obs, y_pred):
    return float(1.0 - np.sum((y_obs - y_pred)**2) / np.sum((y_obs - np.mean(y_obs))**2))


# ==============================================================================
# 2. PARSE SPARC MASTER TABLE (SPARC_Lelli2016c.txt / .mrt)
# ==============================================================================
meta_path = next((p for p in SPARC_META_CANDIDATES if os.path.exists(p)), None)
if meta_path is None:
    raise FileNotFoundError(f"Could not find SPARC master table ({SPARC_META_CANDIDATES}).")

with open(meta_path, "r", encoding="utf-8", errors="replace") as f:
    lines = f.readlines()
dash_idx = [i for i, ln in enumerate(lines) if ln.startswith("----")]
data_start = dash_idx[-1] + 1

meta_rows = []
for ln in lines[data_start:]:
    p = ln.strip().split()
    if len(p) >= 18:
        meta_rows.append({
            'Galaxy': p[0], 'T': int(p[1]), 'D_Mpc': float(p[2]), 'e_D': float(p[3]),
            'f_D': int(p[4]), 'Inc': float(p[5]), 'e_Inc': float(p[6]),
            'L36': float(p[7]), 'Reff_kpc': float(p[9]), 'SBeff': float(p[10]),
            'Rdisk_kpc': float(p[11]), 'SBdisk': float(p[12]),
            'MHI': float(p[13]), 'Vflat': float(p[15]), 'e_Vflat': float(p[16]), 'Q': int(p[17])
        })
meta = pd.DataFrame(meta_rows)

# ==============================================================================
# 3. BTFR AUDIT (DIRECT FROM MASTER TABLE — PUBLISHED D, CONSISTENT D, DIRECT D)
# ==============================================================================
btfr = meta[(meta['Q'] < 3) & (meta['Inc'] >= 30) & (meta['Vflat'] > 0)].copy()
btfr_dir = btfr[btfr['f_D'] >= 2].copy()
btfr_dir_q1 = btfr_dir[btfr_dir['Q'] == 1].copy()

print("\n" + "=" * 108)
print(f" SPARC UNIFIED AUDIT v5 — UNIFIED COSMOLOGY (H0={H0_KM_S_MPC:.2f}, R_{{4/3}}^o={R_PARTITION:.6f})")
print(f" a_0 = c*H0 / R_{{4/3}}^o = {A0_UNIFIED:.5e} m/s^2 | Master Table: {meta_path} (N={len(meta)} galaxies)")
print("=" * 108)
print("\n[A. BARYONIC TULLY-FISHER RELATION (BTFR) — V_flat vs (G * M_bar * a_0)^(1/4)]")

def eval_btfr(sub, label, a_sc, H0_use, rescale_hf=False):
    d_fac = np.where((sub['f_D'].values == 1) & rescale_hf, H0_SPARC_HFLOW / H0_use, 1.0)
    M_kg = (YD_STD * sub['L36'].values + 1.33 * sub['MHI'].values) * 1e9 * (d_fac**2) * M_SUN
    V_obs = sub['Vflat'].values
    V_pred = ((G * M_kg * a_sc)**0.25) / 1000.0
    lv = np.log10(V_obs / V_pred)
    print(f"  {label:50s} | R2_1:1={r2_1to1(V_obs, V_pred):.4f} | V RMS={math.sqrt(np.mean(lv**2)):.4f} dex (bias {np.mean(lv):+.4f})")
    return V_obs, V_pred

v_o_all, v_p_all = eval_btfr(btfr, f"Unified a_0 (All Clean N={len(btfr)}, consistent f_D=1)", A0_UNIFIED, H0_KM_S_MPC, True)
eval_btfr(btfr, f"Unified a_0 (All Clean N={len(btfr)}, published D)", A0_UNIFIED, H0_KM_S_MPC, False)
eval_btfr(btfr, f"Legacy c*H0(73)/(2pi) (All Clean N={len(btfr)})", A0_LEGACY_73, 73.0, False)
eval_btfr(btfr, f"Legacy c*H0(67)/(2pi) (All Clean N={len(btfr)}, rescaled)", A0_LEGACY_67, H0_KM_S_MPC, True)
v_o_dir, v_p_dir = eval_btfr(btfr_dir, f"Unified a_0 (Direct-D f_D>=2 ONLY, N={len(btfr_dir)})", A0_UNIFIED, H0_KM_S_MPC, False)
eval_btfr(btfr_dir, f"Legacy c*H0(73)/(2pi) (Direct-D N={len(btfr_dir)})", A0_LEGACY_73, 73.0, False)
eval_btfr(btfr_dir, f"Legacy c*H0(67)/(2pi) (Direct-D N={len(btfr_dir)})", A0_LEGACY_67, H0_KM_S_MPC, False)
eval_btfr(btfr_dir_q1, f"Unified a_0 (Strict Q=1 Direct-D, N={len(btfr_dir_q1)})", A0_UNIFIED, H0_KM_S_MPC, False)

btfr_out = btfr.copy()
d_fac_c = np.where(btfr_out['f_D'].values == 1, H0_SPARC_HFLOW / H0_KM_S_MPC, 1.0)
btfr_out['M_bar_Msun_Unified'] = (YD_STD * btfr_out['L36'] + 1.33 * btfr_out['MHI']) * 1e9 * (d_fac_c**2)
btfr_out['V_BTFR_Unified_km_s'] = ((G * btfr_out['M_bar_Msun_Unified'] * M_SUN * A0_UNIFIED)**0.25) / 1000.0
btfr_out.to_csv("SPARC_Unified_BTFR_v5.csv", index=False)

# ==============================================================================
# 4. FULL RADIAL ROTATION CURVES AUDIT (REQUIRES Rotmod_LTG.zip OR Rotmod_LTG/)
# ==============================================================================
def load_rotmod_files():
    records = []
    if os.path.exists(ROTMOD_ZIP):
        with zipfile.ZipFile(ROTMOD_ZIP, "r") as zf:
            names = [n for n in zf.namelist() if n.endswith("_rotmod.dat")]
            for nm in sorted(names):
                gal = os.path.basename(nm).replace("_rotmod.dat", "")
                raw = zf.read(nm).decode("utf-8", errors="replace").splitlines()
                for ln in raw:
                    if not ln.strip() or ln.strip().startswith("#"):
                        continue
                    p = ln.split()
                    if len(p) >= 6:
                        records.append((gal, float(p[0]), float(p[1]), float(p[2]), float(p[3]), float(p[4]), float(p[5])))
    elif os.path.isdir("Rotmod_LTG"):
        for fn in sorted(os.listdir("Rotmod_LTG")):
            if not fn.endswith("_rotmod.dat"):
                continue
            gal = fn.replace("_rotmod.dat", "")
            with open(os.path.join("Rotmod_LTG", fn), "r", encoding="utf-8", errors="replace") as f:
                for ln in f:
                    if not ln.strip() or ln.strip().startswith("#"):
                        continue
                    p = ln.split()
                    if len(p) >= 6:
                        records.append((gal, float(p[0]), float(p[1]), float(p[2]), float(p[3]), float(p[4]), float(p[5])))
    return pd.DataFrame(records, columns=['Galaxy', 'Rad_kpc', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul'])

sp = load_rotmod_files()

if len(sp) == 0:
    print("\n[NOTE] Rotmod_LTG.zip not found in current directory -> completed BTFR-only mode.")
    print("       Place Rotmod_LTG.zip in the same folder and re-run for the full 3,391-point radial & SPS audit.")
    # Save BTFR-only 1x2 figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.4), dpi=300)
    m_all = btfr_out['M_bar_Msun_Unified'].values
    v_all = btfr_out['Vflat'].values
    m_grid = np.logspace(7.2, 11.7, 200)
    ax1.loglog(m_grid, ((G * m_grid * M_SUN * A0_UNIFIED)**0.25)/1000.0, color='#0f766e', lw=2.3, label=r'Unified $a_0 = cH_0/\mathcal{R}_{4/3}^\circ$')
    ax1.loglog(m_grid, ((G * m_grid * M_SUN * A0_LEGACY_67)**0.25)/1000.0, color='#64748b', lw=1.6, ls='--', label=r'Legacy $cH_0(67.36)/(2\pi)$')
    ax1.scatter(m_all[btfr_out['f_D']>=2], v_all[btfr_out['f_D']>=2], color='#1d4ed8', edgecolor='white', s=52, zorder=4, label=f'Direct Distances $f_D\\geq 2$ ($N={len(btfr_dir)}$)')
    ax1.scatter(m_all[btfr_out['f_D']==1], v_all[btfr_out['f_D']==1], color='#0f766e', alpha=0.75, s=40, marker='^', zorder=3, label=f'Hubble-Flow $f_D=1$ Rescaled ($N={(btfr_out["f_D"]==1).sum()}$)')
    ax1.set_xlabel(r'Baryonic Mass $M_{\mathrm{bar}}$ [$M_\odot$]', fontsize=11)
    ax1.set_ylabel(r'Asymptotic Flat Velocity $V_{\mathrm{flat}}$ [$\mathrm{km\,s^{-1}}$]', fontsize=11)
    ax1.set_title('A) SPARC Baryonic Tully-Fisher Relation (BTFR)\nUnified $a_0$ ($H_0=67.36$, Zero Free Parameters)', fontsize=11, fontweight='bold', pad=8)
    ax1.legend(loc='upper left', fontsize=8.8); ax1.grid(True, ls=':', alpha=0.5)

    vmax_p = max(v_o_all.max(), v_p_all.max()) * 1.07
    ax2.scatter(v_p_all, v_o_all, color='#0f766e', edgecolor='white', s=52, alpha=0.88, label=f'All Clean ($N={len(btfr)}$, $R^2_{{1:1}}={r2_1to1(v_o_all, v_p_all):.4f}$)')
    ax2.scatter(v_p_dir, v_o_dir, color='#1d4ed8', edgecolor='white', s=58, marker='D', zorder=4, label=f'Direct-D Only ($N={len(btfr_dir)}$, $R^2_{{1:1}}={r2_1to1(v_o_dir, v_p_dir):.4f}$)')
    ax2.plot([0, vmax_p], [0, vmax_p], color='#be123c', ls='--', lw=1.8, label='1:1 Identity Line')
    ax2.set_xlim(15, vmax_p); ax2.set_ylim(15, vmax_p)
    ax2.set_xlabel(r'Predicted $V_{\mathrm{BTFR}} = (G M_{\mathrm{bar}} a_0)^{1/4}$ [$\mathrm{km\,s^{-1}}$]', fontsize=11)
    ax2.set_ylabel(r'Observed $V_{\mathrm{flat}}$ [$\mathrm{km\,s^{-1}}$]', fontsize=11)
    ax2.set_title('B) Predicted vs. Observed $V_{\\mathrm{flat}}$ Parity\nUnified $a_0$ Cosmology ($H_0=67.36$)', fontsize=11, fontweight='bold', pad=8)
    ax2.legend(loc='lower right', fontsize=8.8); ax2.grid(True, ls=':', alpha=0.5)
    plt.tight_layout()
    plt.savefig("SPARC_Unified_BTFR_Figure_v5.png", dpi=300, bbox_inches='tight', pad_inches=0.12)
    plt.close(fig)
    print("Saved SPARC_Unified_BTFR_v5.csv and SPARC_Unified_BTFR_Figure_v5.png")
    sys.exit(0)

# Merge radial curves with SPARC master metadata
sp = sp.merge(meta, on='Galaxy', how='left')
assert sp['Q'].notna().all(), "Some Rotmod galaxies failed to match SPARC_Lelli2016c!"

# Exact signed baryonic velocity squared: V_bar^2 = Vgas*|Vgas| + Yd*Vdisk*|Vdisk| + Yb*Vbul*|Vbul|
vbar_sq = (
    sp['Vgas'].values * np.abs(sp['Vgas'].values)
    + YD_STD * sp['Vdisk'].values * np.abs(sp['Vdisk'].values)
    + YB_STD * sp['Vbul'].values * np.abs(sp['Vbul'].values)
)
sp['V_baryonic'] = np.sqrt(np.maximum(vbar_sq, 0.0))
sp['a_baryonic_SI'] = np.where(sp['Rad_kpc'].values > 0, ((sp['V_baryonic'].values * 1000.0)**2) / (sp['Rad_kpc'].values * KPC_TO_M), 0.0)

# Predictions at published SPARC distances
sp['V_calc_Unified'] = np.where(sp['a_baryonic_SI'] > 0, sp['V_baryonic'] * np.sqrt(nu_rar(sp['a_baryonic_SI'] / A0_UNIFIED)), 0.0)
sp['V_calc_Legacy67'] = np.where(sp['a_baryonic_SI'] > 0, sp['V_baryonic'] * np.sqrt(nu_rar(sp['a_baryonic_SI'] / A0_LEGACY_67)), 0.0)
sp['V_calc_Legacy73'] = np.where(sp['a_baryonic_SI'] > 0, sp['V_baryonic'] * np.sqrt(nu_rar(sp['a_baryonic_SI'] / A0_LEGACY_73)), 0.0)
sp['V_calc_PhiRej'] = np.where(sp['a_baryonic_SI'] > 0, sp['V_baryonic'] * np.sqrt(nu_rar(sp['a_baryonic_SI'] / A0_PHI_REJECT)), 0.0)

# Predictions with consistent f_D=1 Hubble-flow distance rescaling (D_new = D * 73/67.36 -> V_bar_new = V_bar * sqrt(s), g_bar invariant!)
d_scale_pt = np.where(sp['f_D'].values == 1, H0_SPARC_HFLOW / H0_KM_S_MPC, 1.0)
sp['Rad_kpc_Rescaled'] = sp['Rad_kpc'] * d_scale_pt
sp['V_baryonic_Rescaled'] = sp['V_baryonic'] * np.sqrt(d_scale_pt)
sp['V_calc_Unified_Rescaled'] = sp['V_calc_Unified'] * np.sqrt(d_scale_pt)

def audit_radial_subset(sub, v_col, r_col, label):
    y = sub['Vobs'].values
    p = sub[v_col].values
    r_m = sub[r_col].values * KPC_TO_M
    ab = sub['a_baryonic_SI'].values
    r2_g = r2_1to1(y, p)
    m = (y > 0) & (p > 0) & (ab > 0)
    g_obs = ((y[m] * 1000.0)**2) / r_m[m]
    g_calc = ((p[m] * 1000.0)**2) / r_m[m]
    ld = np.log10(g_obs / g_calc)
    rms = math.sqrt(np.mean(ld**2))
    bias = float(np.mean(ld))
    over = float(np.mean(p[m] > y[m]) * 100.0)
    chi2s, gr2s = [], []
    for _, grp in sub.groupby('Galaxy'):
        yg = grp['Vobs'].values; pg = grp[v_col].values; eg = np.maximum(grp['errV'].values, 0.1)
        chi2s.append(np.mean(((yg - pg) / eg)**2))
        if len(grp) >= 3 and np.var(yg) > 0:
            gr2s.append(r2_1to1(yg, pg))
    print(f"  {label:48s} | Gals={sub['Galaxy'].nunique():3d} | Pts={len(sub):4d} | R2_1:1={r2_g:.4f} | MedGalR2={np.median(gr2s):+.3f} | RMS={rms:.4f} | Bias={bias:+.4f} | Over={over:4.1f}% | MedChi2={np.median(chi2s):5.2f}")
    return r2_g, rms, bias, float(np.median(gr2s))

print("\n[B. FULL RADIAL SPARC AUDIT — LOCKED SPS (Yd=0.50, Yb=0.70)]")
m_all = np.ones(len(sp), dtype=bool)
m_cln = (sp['Q'] < 3) & (sp['Inc'] >= 30)
m_str = (sp['Q'] == 1) & (sp['Inc'] >= 30)

r2_all_u, rms_all_u, _, _ = audit_radial_subset(sp[m_all], 'V_calc_Unified', 'Rad_kpc', "Unified a_0 (All 175, Published D)")
r2_cln_u, rms_cln_u, _, _ = audit_radial_subset(sp[m_cln], 'V_calc_Unified', 'Rad_kpc', "Unified a_0 (Clean 153: Q<3,i>=30, Published D)")
r2_str_u, rms_str_u, _, _ = audit_radial_subset(sp[m_str], 'V_calc_Unified', 'Rad_kpc', "Unified a_0 (Strict 93: Q=1,i>=30, Published D)")
audit_radial_subset(sp[m_all], 'V_calc_Unified_Rescaled', 'Rad_kpc_Rescaled', "Unified a_0 (All 175, Rescaled f_D=1 -> 67.36)")
audit_radial_subset(sp[m_cln], 'V_calc_Unified_Rescaled', 'Rad_kpc_Rescaled', "Unified a_0 (Clean 153, Rescaled f_D=1 -> 67.36)")
audit_radial_subset(sp[m_cln], 'V_baryonic', 'Rad_kpc', "Newtonian Baryons Only (Clean 153)")
audit_radial_subset(sp[m_cln], 'V_calc_Legacy73', 'Rad_kpc', "Legacy c*H0(73)/(2pi) (Clean 153)")
audit_radial_subset(sp[m_cln], 'V_calc_Legacy67', 'Rad_kpc', "Legacy c*H0(67.36)/(2pi) (Clean 153)")
audit_radial_subset(sp[m_cln], 'V_calc_PhiRej', 'Rad_kpc', "Rejected Golden Ratio a0*phi (Clean 153)")

# Relaxed SPS Audit (Yd in [0.10, 1.00], Yb = 1.4 * Yd, unweighted SSE)
gal_summary = []
v_relaxed_all = np.zeros(len(sp))
for gal, idx in sp.groupby('Galaxy').groups.items():
    sub = sp.loc[idx]
    vg = sub['Vgas'].values; vd = sub['Vdisk'].values; vb = sub['Vbul'].values
    r_m = sub['Rad_kpc'].values * KPC_TO_M
    yo = sub['Vobs'].values

    def sse_yd(yd):
        yb = 1.4 * yd
        vsq = vg * np.abs(vg) + yd * vd * np.abs(vd) + yb * vb * np.abs(vb)
        vbar = np.sqrt(np.maximum(vsq, 0.0))
        ab = np.where(r_m > 0, ((vbar * 1000.0)**2) / r_m, 0.0)
        vc = np.where(ab > 0, vbar * np.sqrt(nu_rar(ab / A0_UNIFIED)), 0.0)
        return float(np.sum((yo - vc)**2))

    res = minimize_scalar(sse_yd, bounds=(0.10, 1.00), method='bounded')
    yd_opt = float(res.x)
    yb_opt = 1.4 * yd_opt
    vsq = vg * np.abs(vg) + yd_opt * vd * np.abs(vd) + yb_opt * vb * np.abs(vb)
    vbar = np.sqrt(np.maximum(vsq, 0.0))
    ab = np.where(r_m > 0, ((vbar * 1000.0)**2) / r_m, 0.0)
    vc = np.where(ab > 0, vbar * np.sqrt(nu_rar(ab / A0_UNIFIED)), 0.0)
    v_relaxed_all[idx] = vc
    gal_summary.append({
        'Galaxy': gal, 'T': int(sub['T'].iloc[0]), 'Q': int(sub['Q'].iloc[0]),
        'Inc': float(sub['Inc'].iloc[0]), 'f_D': int(sub['f_D'].iloc[0]),
        'N_pts': len(sub), 'Yd_opt': round(yd_opt, 4),
        'at_lower_bound': yd_opt <= 0.105, 'at_upper_bound': yd_opt >= 0.995,
        'R2_locked': round(r2_1to1(yo, sub['V_calc_Unified'].values), 4) if len(sub)>=3 and np.var(yo)>0 else np.nan,
        'R2_relaxed': round(r2_1to1(yo, vc), 4) if len(sub)>=3 and np.var(yo)>0 else np.nan,
    })

sp['V_calc_Unified_RelaxedSPS'] = v_relaxed_all
gs_df = pd.DataFrame(gal_summary)
r2_rel_all = r2_1to1(sp['Vobs'].values, sp['V_calc_Unified_RelaxedSPS'].values)
r2_rel_cln = r2_1to1(sp.loc[m_cln, 'Vobs'].values, sp.loc[m_cln, 'V_calc_Unified_RelaxedSPS'].values)
print(f"\n[C. RELAXED SPS AUDIT (Yd in [0.10, 1.00], Yb=1.4*Yd, Unified a_0)]")
print(f"  All 175: R2_1:1={r2_rel_all:.4f} | Median Yd={gs_df['Yd_opt'].median():.3f} | Lower bound (<=0.105): {gs_df['at_lower_bound'].sum()}/175 | Upper bound (>=0.995): {gs_df['at_upper_bound'].sum()}/175")
gs_cln = gs_df[(gs_df['Q'] < 3) & (gs_df['Inc'] >= 30)]
print(f"  Clean 153: R2_1:1={r2_rel_cln:.4f} | Median Yd={gs_cln['Yd_opt'].median():.3f} | Lower bound: {gs_cln['at_lower_bound'].sum()}/153 | Upper bound: {gs_cln['at_upper_bound'].sum()}/153")

sp.to_csv("SPARC_Unified_Global_v5.csv", index=False)
gs_df.to_csv("SPARC_Unified_GalaxySummary_v5.csv", index=False)

# ==============================================================================
# 5. PUBLICATION-GRADE 2x2 FIGURE -> SPARC_Unified_Figure_v5.png
# ==============================================================================
plt.style.use('default')
fig, axes = plt.subplots(2, 2, figsize=(13.6, 10.6), dpi=300)

sp_cln_df = sp[m_cln & (sp['Vobs'] > 0) & (sp['a_baryonic_SI'] > 0)]
g_bar_p = sp_cln_df['a_baryonic_SI'].values
g_obs_p = ((sp_cln_df['Vobs'].values * 1000.0)**2) / (sp_cln_df['Rad_kpc'].values * KPC_TO_M)

# Panel 1: Radial Acceleration Relation (Clean N=153)
ax = axes[0, 0]
gb_grid = np.logspace(-13.0, -8.3, 300)
ax.loglog(gb_grid, gb_grid, 'k--', lw=1.4, label='Newtonian ($g_{\\mathrm{obs}}=g_{\\mathrm{bar}}$)')
ax.loglog(gb_grid, gb_grid * nu_rar(gb_grid / A0_UNIFIED), color='#0f766e', lw=2.3, label=r'Unified $a_0$ ($a_0 = cH_0/\mathcal{R}_{4/3}^\circ$)')
ax.loglog(gb_grid, gb_grid * nu_rar(gb_grid / A0_LEGACY_67), color='#64748b', ls=':', lw=1.6, label=r'Legacy $cH_0(67.36)/(2\pi)$')
ax.scatter(g_bar_p, g_obs_p, s=6, alpha=0.28, color='#0284c7', rasterized=True, label=f'Clean SPARC ($N=153$, {len(sp_cln_df)} pts)')
ax.set_xlim(1e-13, 5e-9); ax.set_ylim(1e-12, 5e-9)
ax.set_xlabel(r'Baryonic Acceleration $g_{\mathrm{bar}}$ [$\mathrm{m\,s^{-2}}$]', fontsize=10.5)
ax.set_ylabel(r'Observed Acceleration $g_{\mathrm{obs}}$ [$\mathrm{m\,s^{-2}}$]', fontsize=10.5)
ax.set_title(f'A) Clean SPARC RAR (RMS = {rms_cln_u:.4f} dex, Strict Q=1: {rms_str_u:.4f} dex)', fontsize=10.8, fontweight='bold')
ax.legend(loc='upper left', fontsize=8.4, framealpha=0.92); ax.grid(True, ls=':', alpha=0.5)

# Panel 2: V_obs vs V_calc (All 175 Locked vs Relaxed SPS)
ax = axes[0, 1]
vmax_v = max(sp['Vobs'].max(), sp['V_calc_Unified'].max()) * 1.05
ax.scatter(sp['V_calc_Unified'], sp['Vobs'], s=8, alpha=0.25, color='#0f766e', rasterized=True,
           label=f'Locked SPS ($R^2_{{1:1}}={r2_all_u:.4f}$, Clean ${r2_cln_u:.4f}$)')
ax.scatter(sp['V_calc_Unified_RelaxedSPS'], sp['Vobs'], s=6, alpha=0.18, color='#7c3aed', rasterized=True,
           label=f'Relaxed SPS ($R^2_{{1:1}}={r2_rel_all:.4f}$)')
ax.plot([0, vmax_v], [0, vmax_v], color='#be123c', ls='--', lw=1.8, label='1:1 Identity Line')
ax.set_xlim(0, vmax_v); ax.set_ylim(0, vmax_v)
ax.set_xlabel(r'Predicted Velocity $V_{\mathrm{calc}}$ [$\mathrm{km\,s^{-1}}$]', fontsize=10.5)
ax.set_ylabel(r'Observed Velocity $V_{\mathrm{obs}}$ [$\mathrm{km\,s^{-1}}$]', fontsize=10.5)
ax.set_title('B) Global 1:1 Velocity Parity (All 175 SPARC, 3,391 pts)', fontsize=10.8, fontweight='bold')
ax.legend(loc='lower right', fontsize=8.5, framealpha=0.92); ax.grid(True, ls=':', alpha=0.5)

# Panel 3: BTFR (Direct-D vs Rescaled Hubble-Flow)
ax = axes[1, 0]
m_all_b = btfr_out['M_bar_Msun_Unified'].values
v_all_b = btfr_out['Vflat'].values
m_grid = np.logspace(7.2, 11.7, 200)
ax.loglog(m_grid, ((G * m_grid * M_SUN * A0_UNIFIED)**0.25)/1000.0, color='#0f766e', lw=2.3, label=r'Unified $a_0 = cH_0/\mathcal{R}_{4/3}^\circ$')
ax.scatter(m_all_b[btfr_out['f_D']>=2], v_all_b[btfr_out['f_D']>=2], color='#1d4ed8', edgecolor='white', s=48, zorder=4, label=f'Direct Distances $f_D\\geq 2$ ($N={len(btfr_dir)}$)')
ax.scatter(m_all_b[btfr_out['f_D']==1], v_all_b[btfr_out['f_D']==1], color='#0f766e', alpha=0.72, s=36, marker='^', zorder=3, label=f'Hubble-Flow $f_D=1$ Rescaled')
ax.set_xlabel(r'Baryonic Mass $M_{\mathrm{bar}}$ [$M_\odot$]', fontsize=10.5)
ax.set_ylabel(r'Asymptotic $V_{\mathrm{flat}}$ [$\mathrm{km\,s^{-1}}$]', fontsize=10.5)
ax.set_title(f'C) SPARC BTFR ($R^2_{{1:1}}={r2_1to1(v_o_all, v_p_all):.4f}$; Direct-D: ${r2_1to1(v_o_dir, v_p_dir):.4f}$)', fontsize=10.8, fontweight='bold')
ax.legend(loc='upper left', fontsize=8.4, framealpha=0.92); ax.grid(True, ls=':', alpha=0.5)

# Panel 4: Relaxed SPS Disk M/L Distribution
ax = axes[1, 1]
ax.hist(gs_df['Yd_opt'], bins=18, color='#0f766e', alpha=0.80, edgecolor='black', lw=0.7)
ax.axvline(0.50, color='#be123c', ls='--', lw=2.0, label='Standard SPS ($\\Upsilon_d = 0.50$)')
ax.axvline(gs_df['Yd_opt'].median(), color='#1d4ed8', ls='-', lw=2.0, label=f'Recovered Median ($\\Upsilon_d = {gs_df["Yd_opt"].median():.3f}$)')
ax.set_xlabel(r'Optimal Disk Mass-to-Light Ratio $\Upsilon_d$ [$M_\odot/L_\odot$]', fontsize=10.5)
ax.set_ylabel('Number of Galaxies', fontsize=10.5)
ax.set_title('D) Relaxed SPS Distribution (175 SPARC Galaxies)', fontsize=10.8, fontweight='bold')
ax.legend(loc='upper right', fontsize=8.6, framealpha=0.92); ax.grid(True, axis='y', ls=':', alpha=0.5)

plt.tight_layout()
plt.savefig("SPARC_Unified_Figure_v5.png", dpi=300, bbox_inches='tight', pad_inches=0.12)
if 'ipykernel' in sys.modules:
    plt.show()
else:
    plt.close(fig)
print("Saved SPARC_Unified_Figure_v5.png, SPARC_Unified_Global_v5.csv, SPARC_Unified_GalaxySummary_v5.csv")