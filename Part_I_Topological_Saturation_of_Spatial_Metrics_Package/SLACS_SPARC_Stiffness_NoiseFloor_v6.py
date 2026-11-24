#!/usr/bin/env python3
"""
SLACS IX & SPARC Joint Transition-Stiffness, Dynamic a_0(z) & HST Noise-Floor Audit (v6)
========================================================================================
Standalone, reviewer-ready script (CLI + Google Colab compatible) that pulls
directly from the Part I Zenodo package files in the current directory and
generates Figure 4 ('SLACS_SPARC_Stiffness_NoiseFloor_v6.png') and its
accompanying CSV summary ('SLACS_SPARC_Stiffness_NoiseFloor_v6.csv').

Inputs auto-detected in the script/working directory:
  - SLACS (N = 70 SLACS IX + N = 40 S4TM):
      Priority: SLACS_DB_v5.csv -> Phantom_Clean_Master_SLACS_SPARC-v4/v3/v2.xlsx -> SLACS_ORIGN_DB.xlsx
  - SPARC Clean (Q < 3, i >= 30 deg, N = 153):
      Priority: SPARC_Lelli2016c.txt + Rotmod_LTG.zip -> Phantom_Clean_Master_SLACS_SPARC-v4/v3/v2.xlsx (all 175 fallback)

Usage:
  python3 SLACS_SPARC_Stiffness_NoiseFloor_v6.py
  python3 SLACS_SPARC_Stiffness_NoiseFloor_v6.py --H0 67.36 --Om 0.31391 --mc-trials 2000
"""

import argparse
import io
import math
import os
import sys
import zipfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, linregress
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Exact CODATA / IAU Physical Constants ---
C_M_S = 299792458.0
G_SI = 6.67430e-11
M_SUN_KG = 1.98847e30
MPC_M = 3.085677581491367e22
KPC_M = MPC_M / 1.0e3
ARCSEC_RAD = math.pi / (180.0 * 3600.0)

# --- Lambda-TSM Structural Constant ---
R_LTSM = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * math.pi**2))  # ~= 5.36710706
JAFFE_RE_OVER_AS = 0.744728  # Exact 2D half-mass projection ratio for Jaffe (gamma=2.0)
SALP_CHAB_OFFSET = 0.2464    # Mean Auger+2009 Salpeter - Chabrier offset [dex]


def build_cosmology(h0_kms, omega_m):
    omega_l = 1.0 - omega_m
    h0_si = (h0_kms * 1000.0) / MPC_M
    z_grid = np.linspace(0.0, 5.0, 50001)
    ez_grid = np.sqrt(omega_m * (1.0 + z_grid)**3 + omega_l)
    dz = z_grid[1] - z_grid[0]
    chi_grid = np.concatenate(
        [[0.0], np.cumsum(0.5 * (1.0 / ez_grid[:-1] + 1.0 / ez_grid[1:]) * dz)]
    ) * (C_M_S / h0_si)  # meters

    def chi(z):
        return np.interp(z, z_grid, chi_grid)

    def ez(z):
        return np.sqrt(omega_m * (1.0 + np.asarray(z, dtype=float))**3 + omega_l)

    return chi, ez, h0_si


def _pick_col(df, candidates):
    lower_map = {str(c).strip().lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    for cand in candidates:
        cl = cand.lower()
        for k, orig in lower_map.items():
            if cl in k:
                return orig
    return None


def load_slacs_all(base_dir):
    """Load SLACS IX (N=70) and, when available, S4TM (N=40) from package files."""
    # 1. Try SLACS_DB_v5.csv or SLACS_Phantom_Metric_Lensing_v5.csv
    for fname in ("SLACS_DB_v5.csv", "SLACS_Phantom_Metric_Lensing_v5.csv"):
        p = os.path.join(base_dir, fname)
        if os.path.exists(p):
            df = pd.read_csv(p)
            c_surv = _pick_col(df, ["Survey", "survey", "sample"])
            c_zl = _pick_col(df, ["z_lens", "zl", "z_d"])
            c_zs = _pick_col(df, ["z_source", "zs", "z_s"])
            c_th = _pick_col(df, ["theta_E_obs_arcsec", "theta_e_arcsec", "theta_e", "thetae"])
            c_re = _pick_col(df, ["R_eff_arcsec", "r_eff_arcsec", "re_i", "r_eff"])
            c_ms = _pick_col(df, ["log10_Mstar_Salpeter", "logm_salp", "logm_salpeter", "mstar_salp", "salpeter"])
            c_mc = _pick_col(df, ["log10_Mstar_Chabrier", "logm_chab", "logm_chabrier", "chabrier"])
            c_sg = _pick_col(df, ["sigma_obs_km_s", "sigma_v", "sigma", "vdisp"])
            if all(c is not None for c in (c_zl, c_zs, c_th, c_re)) and (c_ms is not None or c_mc is not None):
                already_rescaled = ("v5" in fname.lower() and "lensing" in fname.lower())
                if c_surv is not None:
                    is_ix = df[c_surv].astype(str).str.contains("SLACS_IX|Auger", case=False, na=False).values
                else:
                    is_ix = np.ones(len(df), dtype=bool)
                ix_df = df[is_ix].dropna(subset=[c_zl, c_zs, c_th, c_re, c_ms]).copy()
                if len(ix_df) >= 65:
                    res = {
                        "source": fname,
                        "ix": {
                            "z_lens": ix_df[c_zl].astype(float).values,
                            "z_source": ix_df[c_zs].astype(float).values,
                            "theta_E": ix_df[c_th].astype(float).values,
                            "R_eff": ix_df[c_re].astype(float).values,
                            "logM_salp": ix_df[c_ms].astype(float).values,
                            "sigma": ix_df[c_sg].astype(float).values if c_sg is not None else np.full(len(ix_df), np.nan),
                        },
                        "s4": None,
                        "masses_at_cat70": not already_rescaled,
                    }
                    if c_surv is not None and c_mc is not None:
                        s4_df = df[~is_ix].dropna(subset=[c_zl, c_zs, c_th, c_re, c_mc]).copy()
                        if len(s4_df) >= 35:
                            res["s4"] = {
                                "z_lens": s4_df[c_zl].astype(float).values,
                                "z_source": s4_df[c_zs].astype(float).values,
                                "theta_E": s4_df[c_th].astype(float).values,
                                "R_eff": s4_df[c_re].astype(float).values,
                                "logM_chab": s4_df[c_mc].astype(float).values,
                            }
                    return res

    # 2. Try unified Master Excel workbooks (v4 -> v3 -> v2)
    for pat in (
        "Phantom_Clean_Master_SLACS_SPARC-v4.xlsx",
        "Phantom_Clean_Master_SLACS_SPARC-v3.xlsx",
        "Phantom_Clean_Master_SLACS_SPARC-v2.xlsx",
    ):
        p = os.path.join(base_dir, pat)
        if os.path.exists(p):
            df = pd.read_excel(p, sheet_name="SLACS_110")
            is_ix = df["Survey"].astype(str).str.contains("SLACS_IX", na=False)
            ix = df[is_ix].copy()
            s4 = df[~is_ix].copy()
            return {
                "source": f"{pat} [SLACS_110]",
                "ix": {
                    "z_lens": ix["z_lens"].astype(float).values,
                    "z_source": ix["z_source"].astype(float).values,
                    "theta_E": ix["theta_E_obs_arcsec"].astype(float).values,
                    "R_eff": ix["R_eff_arcsec"].astype(float).values,
                    "logM_salp": ix["log10_Mstar_Salpeter"].astype(float).values,
                    "sigma": ix["sigma_obs_km_s"].astype(float).values,
                },
                "s4": {
                    "z_lens": s4["z_lens"].astype(float).values,
                    "z_source": s4["z_source"].astype(float).values,
                    "theta_E": s4["theta_E_obs_arcsec"].astype(float).values,
                    "R_eff": s4["R_eff_arcsec"].astype(float).values,
                    "logM_chab": s4["log10_Mstar_Chabrier"].astype(float).values,
                } if len(s4) >= 35 else None,
                "masses_at_cat70": True,
            }

    # 3. Try raw SLACS_ORIGN_DB.xlsx (Auger+2009 Table 3/4 sheet)
    p_orig = os.path.join(base_dir, "SLACS_ORIGN_DB.xlsx")
    if os.path.exists(p_orig):
        xls = pd.ExcelFile(p_orig)
        for sh in xls.sheet_names:
            df = pd.read_excel(p_orig, sheet_name=sh)
            c_zl = _pick_col(df, ["z_lens", "zl", "z_l"])
            c_zs = _pick_col(df, ["z_source", "zs", "z_s"])
            c_ms = _pick_col(df, ["log10_Mstar_Salpeter", "Salp", "logM_Salp"])
            c_re = _pick_col(df, ["r_e_I", "R_eff_arcsec", "re_I", "r_eff"])
            c_th = _pick_col(df, ["theta_E", "r_Ein", "thetaE"])
            c_sg = _pick_col(df, ["sigma", "vdisp"])
            if all(c is not None for c in (c_zl, c_zs, c_ms, c_re, c_th)) and len(df) >= 65:
                sub = df.dropna(subset=[c_zl, c_zs, c_ms, c_re, c_th]).copy()
                zl = sub[c_zl].astype(float).values
                th_raw = sub[c_th].astype(float).values
                if "ein" in str(c_th).lower() and np.median(th_raw) > 2.0:
                    chi70, _, _ = build_cosmology(70.0, 0.30)
                    dl70_kpc = (chi70(zl) / (1.0 + zl)) / KPC_M
                    th_as = (th_raw / dl70_kpc) / ARCSEC_RAD
                else:
                    th_as = th_raw
                return {
                    "source": f"SLACS_ORIGN_DB.xlsx [{sh}]",
                    "ix": {
                        "z_lens": zl,
                        "z_source": sub[c_zs].astype(float).values,
                        "theta_E": th_as,
                        "R_eff": sub[c_re].astype(float).values,
                        "logM_salp": sub[c_ms].astype(float).values,
                        "sigma": sub[c_sg].astype(float).values if c_sg is not None else np.full(len(sub), np.nan),
                    },
                    "s4": None,
                    "masses_at_cat70": True,
                }

    sys.exit(
        "ERROR: Could not find SLACS data (SLACS_DB_v5.csv, Phantom_Clean_Master_SLACS_SPARC-v4.xlsx, "
        "or SLACS_ORIGN_DB.xlsx) in the working directory."
    )


def load_sparc_clean(base_dir):
    """Load SPARC rotation curves (Clean Q < 3, i >= 30 deg, N = 153 from raw archive; fallback to Excel)."""
    txt_path = os.path.join(base_dir, "SPARC_Lelli2016c.txt")
    zip_path = os.path.join(base_dir, "Rotmod_LTG.zip")

    if os.path.exists(txt_path) and os.path.exists(zip_path):
        meta = {}
        with open(txt_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        start_idx = 0
        for idx, ln in enumerate(lines):
            if ln.strip().startswith("---"):
                start_idx = idx + 1
        for ln in lines[start_idx:]:
            parts = ln.split()
            if len(parts) < 18:
                continue
            try:
                gname = parts[0].strip()
                inc = float(parts[5])
                qflag = int(parts[17])
                meta[gname] = (qflag, inc)
            except ValueError:
                continue

        rows = []
        with zipfile.ZipFile(zip_path, "r") as zf:
            for member in zf.namelist():
                if not member.endswith("_rotmod.dat"):
                    continue
                gname = os.path.basename(member).replace("_rotmod.dat", "")
                if gname not in meta:
                    continue
                qflag, inc = meta[gname]
                if not (qflag < 3 and inc >= 30.0):
                    continue
                with io.TextIOWrapper(zf.open(member), encoding="utf-8", errors="ignore") as mf:
                    for ln in mf:
                        s = ln.strip()
                        if not s or s.startswith("#"):
                            continue
                        p = s.split()
                        if len(p) < 6:
                            continue
                        r_kpc = float(p[0])
                        if r_kpc <= 0:
                            continue
                        rows.append((gname, r_kpc, float(p[1]), float(p[2]), float(p[3]), float(p[4]), float(p[5])))
        df_sp = pd.DataFrame(rows, columns=["Galaxy", "Rad_kpc", "Vobs", "errV", "Vgas", "Vdisk", "Vbulge"])
        return df_sp, f"SPARC Clean (N={df_sp['Galaxy'].nunique()})"

    for pat in (
        "Phantom_Clean_Master_SLACS_SPARC-v4.xlsx",
        "Phantom_Clean_Master_SLACS_SPARC-v3.xlsx",
        "Phantom_Clean_Master_SLACS_SPARC-v2.xlsx",
    ):
        p = os.path.join(base_dir, pat)
        if os.path.exists(p):
            df_sp = pd.read_excel(p, sheet_name="SPARC_175")
            df_sp = df_sp[["Galaxy", "Rad_kpc", "Vobs", "errV", "Vgas", "Vdisk", "Vbulge"]].dropna()
            return df_sp, f"SPARC All (N={df_sp['Galaxy'].nunique()})"

    sys.exit("ERROR: Could not find SPARC_Lelli2016c.txt + Rotmod_LTG.zip or Master Excel workbook.")


def _solve_jaffe_batch(zd, zs, th_e, r_eff_as, logm, chi_new, ez_new, a0_0, mode="nu_k", k_val=1.0, dynamic_z=False):
    d_l = chi_new(zd) / (1.0 + zd)
    d_s = chi_new(zs) / (1.0 + zs)
    d_ls = (chi_new(zs) - chi_new(zd)) / (1.0 + zs)
    beta = d_ls / d_s
    re_obs_kpc = (d_l * th_e * ARCSEC_RAD) / KPC_M
    reff_kpc = (d_l * r_eff_as * ARCSEC_RAD) / KPC_M
    a_s_m = (reff_kpc / JAFFE_RE_OVER_AS) * KPC_M
    m_kg = (10.0**logm) * M_SUN_KG

    u_nodes = np.linspace(0.0, 12.5, 321)
    du = u_nodes[1] - u_nodes[0]
    w_u = np.ones_like(u_nodes) * du
    w_u[0] *= 0.5
    w_u[-1] *= 0.5
    cosh_u = np.cosh(u_nodes)

    re_calc = np.zeros(len(zd))
    for i in range(len(zd)):
        a0_i = a0_0 * (float(ez_new(zd[i])) if dynamic_z else 1.0)
        pref = (4.0 * d_l[i] * beta[i]) / (C_M_S**2)

        def f_val(r_m):
            r = r_m * cosh_u
            g_n = G_SI * m_kg[i] / (r * (r + a_s_m[i]))
            y = np.maximum(g_n / a0_i, 1e-12)
            if mode == "gr":
                nu = 1.0
            elif mode == "simple":
                nu = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 / y))
            else:
                nu = (1.0 - np.exp(-(y**(0.5 * k_val))))**(-1.0 / k_val)
            return pref * np.sum(g_n * nu * w_u) - 1.0

        lo, hi = -2.0, 2.5
        if f_val((10.0**lo) * KPC_M) < 0.0:
            re_calc[i] = 0.0
            continue
        for _ in range(42):
            mid = 0.5 * (lo + hi)
            if f_val((10.0**mid) * KPC_M) > 0.0:
                lo = mid
            else:
                hi = mid
        re_calc[i] = 10.0**(0.5 * (lo + hi))

    r2_1to1 = 1.0 - np.sum((re_obs_kpc - re_calc)**2) / np.sum((re_obs_kpc - np.mean(re_obs_kpc))**2)
    lr = linregress(re_calc, re_obs_kpc)
    parity = 1.0 - abs(np.sum(re_obs_kpc) - np.sum(re_calc)) / np.sum(re_obs_kpc)
    return re_obs_kpc, re_calc, reff_kpc, a_s_m, m_kg, float(r2_1to1), float(lr.rvalue**2), float(parity)


def eval_slacs_ix(slacs_pkg, h0_kms, omega_m, k_grid, n_mc=2000, seed=42):
    chi_new, ez_new, h0_si = build_cosmology(h0_kms, omega_m)
    chi_cat, _, _ = build_cosmology(70.0, 0.30)
    a0_0 = (C_M_S * h0_si) / R_LTSM

    ix = slacs_pkg["ix"]
    zd = ix["z_lens"]
    zs = ix["z_source"]
    th_e = ix["theta_E"]
    r_eff_as = ix["R_eff"]
    sigma = ix["sigma"]
    if slacs_pkg["masses_at_cat70"]:
        dlogm_ix = 2.0 * np.log10((chi_new(zd) * (1.0 + zd)) / (chi_cat(zd) * (1.0 + zd)))
        logm = ix["logM_salp"] + dlogm_ix
    else:
        logm = ix["logM_salp"].copy()

    re_obs, re_gr, reff_kpc, a_s_m, m_kg, r2_gr, _, par_gr = _solve_jaffe_batch(
        zd, zs, th_e, r_eff_as, logm, chi_new, ez_new, a0_0, mode="gr"
    )
    _, re_simp, _, _, _, r2_simp, _, par_simp = _solve_jaffe_batch(
        zd, zs, th_e, r_eff_as, logm, chi_new, ez_new, a0_0, mode="simple"
    )

    slacs_scan = []
    re_k1 = None
    for kv in k_grid:
        _, rc, _, _, _, r2, _, par = _solve_jaffe_batch(
            zd, zs, th_e, r_eff_as, logm, chi_new, ez_new, a0_0, mode="nu_k", k_val=kv, dynamic_z=False
        )
        slacs_scan.append((kv, r2, par))
        if abs(kv - 1.0) < 1e-9:
            re_k1 = rc

    # Dynamic a_0(z_d) comparison on SLACS IX (N=70)
    _, re_dyn, _, _, _, r2_dyn, lr2_dyn, par_dyn = _solve_jaffe_batch(
        zd, zs, th_e, r_eff_as, logm, chi_new, ez_new, a0_0, mode="nu_k", k_val=1.0, dynamic_z=True
    )
    rho_z_dyn, p_z_dyn = spearmanr(zd, re_obs / re_dyn)

    # Combined N=110 dynamic a_0(z_d) check (if S4TM present)
    comb_dyn = None
    if slacs_pkg["s4"] is not None:
        s4 = slacs_pkg["s4"]
        zd_4, zs_4 = s4["z_lens"], s4["z_source"]
        dlogm_4 = 2.0 * np.log10((chi_new(zd_4) * (1.0 + zd_4)) / (chi_cat(zd_4) * (1.0 + zd_4))) if slacs_pkg["masses_at_cat70"] else 0.0
        logm_4_mid = s4["logM_chab"] + 0.5 * SALP_CHAB_OFFSET + dlogm_4
        re_obs_c = np.concatenate([re_obs, _solve_jaffe_batch(zd_4, zs_4, s4["theta_E"], s4["R_eff"], logm_4_mid, chi_new, ez_new, a0_0, mode="nu_k", k_val=1.0, dynamic_z=True)[0]])
        re_dyn_4 = _solve_jaffe_batch(zd_4, zs_4, s4["theta_E"], s4["R_eff"], logm_4_mid, chi_new, ez_new, a0_0, mode="nu_k", k_val=1.0, dynamic_z=True)[1]
        re_dyn_c = np.concatenate([re_dyn, re_dyn_4])
        r2_c = 1.0 - np.sum((re_obs_c - re_dyn_c)**2) / np.sum((re_obs_c - np.mean(re_obs_c))**2)
        lr_c = linregress(re_dyn_c, re_obs_c)
        par_c = 1.0 - abs(np.sum(re_obs_c) - np.sum(re_dyn_c)) / np.sum(re_obs_c)
        comb_dyn = {
            "N": len(re_obs_c),
            "r2": float(r2_c),
            "lin_r2": float(lr_c.rvalue**2),
            "mean_ratio": float(np.mean(re_obs_c / re_dyn_c)),
            "parity": float(par_c),
        }

    # Residual diagnostics at static k = 1.0
    ratio = re_obs / re_k1
    rho_rad, p_rad = spearmanr(re_obs / reff_kpc, ratio)
    rho_z, p_z = spearmanr(zd, ratio)
    valid_sig = np.isfinite(sigma) & (sigma > 0)
    rho_sig, p_sig = spearmanr(sigma[valid_sig], ratio[valid_sig])
    g_n_at_re = G_SI * m_kg / ((re_obs * KPC_M) * (re_obs * KPC_M + a_s_m))
    med_y_e = float(np.median(g_n_at_re / a0_0))

    # Monte Carlo HST measurement noise-floor realizations (0.08 dex & 0.10 dex mass scatter)
    _, re_k1_hi, _, _, _, _, _, _ = _solve_jaffe_batch(
        zd, zs, th_e, r_eff_as, logm + 0.05, chi_new, ez_new, a0_0, mode="nu_k", k_val=1.0, dynamic_z=False
    )
    dln_re_dlogm = np.log(re_k1_hi / re_k1) / 0.05
    rng = np.random.default_rng(seed)
    mc_results = {}
    for sig_dex in (0.08, 0.10):
        draws = np.empty(n_mc)
        denom = np.sum((re_k1 - np.mean(re_k1))**2)
        for j in range(n_mc):
            dm = rng.normal(0.0, sig_dex, size=len(zd))
            re_sim = re_k1 * np.exp(dln_re_dlogm * dm)
            draws[j] = 1.0 - np.sum((re_k1 - re_sim)**2) / denom
        mc_results[sig_dex] = draws

    return {
        "N": len(zd),
        "r2_gr": r2_gr,
        "par_gr": par_gr,
        "r2_simple": r2_simp,
        "par_simple": par_simp,
        "scan": slacs_scan,
        "dyn_ix": {
            "r2": r2_dyn, "lin_r2": lr2_dyn, "parity": par_dyn,
            "mean_ratio": float(np.mean(re_obs / re_dyn)),
            "rho_z": rho_z_dyn, "p_z": p_z_dyn,
        },
        "dyn_comb": comb_dyn,
        "diag": {
            "rho_rad": rho_rad, "p_rad": p_rad,
            "rho_z": rho_z, "p_z": p_z,
            "rho_sig": rho_sig, "p_sig": p_sig,
            "med_y_e": med_y_e,
        },
        "mc": mc_results,
    }


def eval_sparc_stiffness(df_sp, h0_kms, k_grid):
    h0_si = (h0_kms * 1000.0) / MPC_M
    a0 = (C_M_S * h0_si) / R_LTSM

    v_obs = df_sp["Vobs"].values.astype(float)
    v_gas = df_sp["Vgas"].values.astype(float)
    v_disk = df_sp["Vdisk"].values.astype(float)
    v_bulge = df_sp["Vbulge"].values.astype(float)
    rad_m = df_sp["Rad_kpc"].values.astype(float) * KPC_M
    gals = df_sp["Galaxy"].values
    uniq_gals = df_sp["Galaxy"].unique()

    vbar_sq = np.maximum(0.0, np.sign(v_gas) * v_gas**2 + 0.50 * v_disk**2 + 0.70 * v_bulge**2)
    g_bar = (vbar_sq * 1.0e6) / rad_m
    g_obs = (v_obs**2 * 1.0e6) / rad_m
    valid_rar = (vbar_sq > 0) & (g_obs > 0)

    def eval_nu(mode, k_val=1.0):
        y = np.maximum(g_bar / a0, 1e-12)
        if mode == "newt":
            nu = np.ones_like(y)
        elif mode == "simple":
            nu = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 / y))
        else:
            nu = (1.0 - np.exp(-(y**(0.5 * k_val))))**(-1.0 / k_val)
        v_calc = np.sqrt(vbar_sq * nu)
        g_calc = g_bar * nu
        r2_gals = []
        for g in uniq_gals:
            m = (gals == g)
            if np.sum(m) >= 3:
                vo, vc = v_obs[m], v_calc[m]
                tot = np.sum((vo - np.mean(vo))**2)
                if tot > 0:
                    r2_gals.append(1.0 - np.sum((vo - vc)**2) / tot)
        log_diff = np.log10(g_obs[valid_rar] / g_calc[valid_rar])
        rms_dex = float(np.sqrt(np.mean(log_diff**2)))
        return float(np.median(r2_gals)), rms_dex

    r2_newt, _ = eval_nu("newt")
    r2_simp, rms_simp = eval_nu("simple")
    scan = []
    for kv in k_grid:
        mr2, rms = eval_nu("nu_k", kv)
        scan.append((kv, mr2, rms))

    return {
        "N_gal": len(uniq_gals),
        "r2_newt": r2_newt,
        "r2_simple": r2_simp,
        "rms_simple": rms_simp,
        "scan": scan,
    }


def make_figure(slacs_res, sparc_res, sparc_label, out_png):
    k_vals = np.array([x[0] for x in slacs_res["scan"]])
    slacs_r2 = np.array([x[1] for x in slacs_res["scan"]])
    sparc_r2 = np.array([x[1] for x in sparc_res["scan"]])

    idx_k1 = int(np.argmin(np.abs(k_vals - 1.0)))
    slacs_k1 = slacs_r2[idx_k1]
    sparc_k1 = sparc_r2[idx_k1]

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.linewidth": 1.0,
    })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.9), dpi=300)

    # --- Panel A: Joint Stiffness Scan ---
    c_slacs = "#0F5C6E"
    c_sparc = "#C45B28"

    mask_plot = k_vals >= 0.60
    ax1.plot(k_vals[mask_plot], slacs_r2[mask_plot], "-o", color=c_slacs, lw=2.2, ms=5.5,
             label=f"SLACS IX Strong Lensing (N={slacs_res['N']}, $R^2_{{1:1}}$)")
    ax1.plot(k_vals[mask_plot], sparc_r2[mask_plot], "-s", color=c_sparc, lw=2.2, ms=5.0,
             label=f"{sparc_label} (Med. Gal. $R^2_{{1:1}}$)")

    ax1.axvline(1.0, color="#222222", ls="--", lw=1.4, label=r"Canonical RAR ($k = 1.0$, 0 free params)")
    ax1.axvline(2.0 / 3.0, color="#888888", ls=":", lw=1.3, label=r"Softened ($k = 2/3$, jointly rejected)")

    ax1.scatter([1.0], [slacs_k1], s=95, color=c_slacs, edgecolor="white", zorder=6)
    ax1.scatter([1.0], [sparc_k1], s=95, color=c_sparc, edgecolor="white", zorder=6)

    ax1.axhline(slacs_res["r2_simple"], color=c_slacs, ls="-.", lw=1.0, alpha=0.55,
                label=f"Simple MOND: SLACS {slacs_res['r2_simple']:.4f} / SPARC {sparc_res['r2_simple']:+.4f}")
    ax1.axhline(sparc_res["r2_simple"], color=c_sparc, ls="-.", lw=1.0, alpha=0.55)

    ax1.annotate(f"k=1.0\n{slacs_k1:.4f}", xy=(1.0, slacs_k1), xytext=(1.13, slacs_k1 + 0.03),
                 fontsize=8.5, fontweight="bold", color=c_slacs,
                 arrowprops=dict(arrowstyle="->", color=c_slacs, lw=1.0))
    ax1.annotate(f"k=1.0\n{sparc_k1:+.4f}", xy=(1.0, sparc_k1), xytext=(1.08, sparc_k1 + 0.08),
                 fontsize=8.5, fontweight="bold", color=c_sparc,
                 arrowprops=dict(arrowstyle="->", color=c_sparc, lw=1.0))

    ax1.set_title("A) Joint Stiffness Scan: $\\nu_k(y)=[1-\\exp(-y^{k/2})]^{-1/k}$", fontsize=10.5, fontweight="bold", pad=8)
    ax1.set_xlabel("Transition Stiffness Exponent $k$")
    ax1.set_ylabel("Strict 1:1 Parity Metric ($R^2_{1:1}$)")
    ax1.set_xlim(0.60, 2.05)
    ymin = min(np.min(sparc_r2[mask_plot]) - 0.08, 0.0)
    ax1.set_ylim(ymin - 0.04, 0.92)
    ax1.grid(True, ls="--", alpha=0.30)
    ax1.legend(loc="lower right", fontsize=7.8, framealpha=0.92)

    # --- Panel B: Monte Carlo HST Noise-Floor Distribution ---
    mc08 = slacs_res["mc"][0.08]
    mc10 = slacs_res["mc"][0.10]

    bins = np.linspace(0.45, 0.94, 38)
    ax2.hist(mc08, bins=bins, density=True, alpha=0.58, color="#1F77B4", edgecolor="#124A73", lw=0.6,
             label=f"HST 0.08 dex scatter: ${np.mean(mc08):.4f} \\pm {np.std(mc08):.4f}$")
    ax2.hist(mc10, bins=bins, density=True, alpha=0.45, color="#D97736", edgecolor="#8F491B", lw=0.6,
             label=f"HST 0.10 dex scatter: ${np.mean(mc10):.4f} \\pm {np.std(mc10):.4f}$")

    ax2.axvline(slacs_k1, color="#0B3C49", ls="--", lw=2.3,
                label=f"Observed SLACS IX ($k=1.0$): $R^2_{{1:1}} = {slacs_k1:.4f}$")
    ax2.axvline(slacs_res["r2_gr"], color="#7F7F7F", ls=":", lw=1.8,
                label=f"Pure Baryonic GR ($\\nu=1$): $R^2_{{1:1}} = {slacs_res['r2_gr']:.4f}$")

    ax2.set_title("B) SLACS IX HST Noise-Floor Audit (Monte Carlo)", fontsize=10.5, fontweight="bold", pad=8)
    ax2.set_xlabel("Simulated Strict 1:1 $R^2_{1:1}$ Under Zero-Scatter Gravity")
    ax2.set_ylabel("Probability Density")
    ax2.set_xlim(0.45, 0.94)
    ax2.grid(True, ls="--", alpha=0.30)
    ax2.legend(loc="upper left", fontsize=8.0, framealpha=0.92)

    plt.tight_layout(w_pad=2.4)
    plt.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    parser = argparse.ArgumentParser(description="SLACS IX & SPARC Transition-Stiffness, Dynamic a_0(z) & HST Noise-Floor Audit (v6).")
    parser.add_argument("--data-dir", default=base_dir, help="Directory containing the Part I package data files.")
    parser.add_argument("--H0", type=float, default=67.36, help="Hubble constant H0 [km/s/Mpc] (default: 67.36).")
    parser.add_argument("--Om", type=float, default=0.31391, help="Matter density Omega_m,0 (default: 0.31391).")
    parser.add_argument("--mc-trials", type=int, default=2000, help="Monte Carlo realizations (default: 2000).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    args, _ = parser.parse_known_args()

    k_grid = np.array([0.50, 2.0 / 3.0, 0.75, 0.85, 0.90, 0.95, 1.00, 1.05, 1.15, 1.25, 1.50, 1.75, 2.00])

    slacs_pkg = load_slacs_all(args.data_dir)
    sparc_df, sparc_label = load_sparc_clean(args.data_dir)

    print("================================================================================")
    print("  SLACS IX & SPARC JOINT TRANSITION-STIFFNESS & HST NOISE-FLOOR AUDIT (v6)")
    print("================================================================================")
    print(f"Cosmology  : H0 = {args.H0:.2f} km/s/Mpc, Omega_m = {args.Om:.5f}, R_4/3 = {R_LTSM:.8f}")
    print(f"SLACS src  : {slacs_pkg['source']} (IX N = {len(slacs_pkg['ix']['z_lens'])})")
    print(f"SPARC src  : {sparc_label}\n")

    slacs_res = eval_slacs_ix(slacs_pkg, args.H0, args.Om, k_grid, n_mc=args.mc_trials, seed=args.seed)
    sparc_res = eval_sparc_stiffness(sparc_df, args.H0, k_grid)

    print(f"{'k':>7s} | {'SLACS R^2_1:1':>13s} | {'SLACS Parity':>12s} | {'SPARC Med R^2':>13s} | {'SPARC RMS (dex)':>15s}")
    print("-" * 72)
    csv_rows = []
    for (kv, sr2, spar), (_, pr2, prms) in zip(slacs_res["scan"], sparc_res["scan"]):
        tag = " <-- PEAK" if abs(kv - 1.0) < 1e-9 else ""
        print(f"{kv:7.3f} | {sr2:13.4f} | {spar*100:11.2f}% | {pr2:+13.4f} | {prms:15.4f}{tag}")
        csv_rows.append({
            "Model": f"nu_k (k={kv:.4f})", "k": round(float(kv), 4),
            "SLACS_IX_R2_1to1": round(sr2, 4), "SLACS_IX_SumParity_pct": round(spar * 100.0, 2),
            "SPARC_MedGal_R2_1to1": round(pr2, 4), "SPARC_RAR_RMS_dex": round(prms, 4),
        })

    print("-" * 72)
    print(f"{'GR/Newt':>7s} | {slacs_res['r2_gr']:13.4f} | {slacs_res['par_gr']*100:11.2f}% | {sparc_res['r2_newt']:+13.4f} | {'---':>15s}")
    print(f"{'Simple':>7s} | {slacs_res['r2_simple']:13.4f} | {slacs_res['par_simple']*100:11.2f}% | {sparc_res['r2_simple']:+13.4f} | {sparc_res['rms_simple']:15.4f}")

    csv_rows.append({"Model": "Baryons_GR_Newtonian", "k": np.nan, "SLACS_IX_R2_1to1": round(slacs_res["r2_gr"], 4), "SLACS_IX_SumParity_pct": round(slacs_res["par_gr"]*100, 2), "SPARC_MedGal_R2_1to1": round(sparc_res["r2_newt"], 4), "SPARC_RAR_RMS_dex": np.nan})
    csv_rows.append({"Model": "Simple_MOND", "k": np.nan, "SLACS_IX_R2_1to1": round(slacs_res["r2_simple"], 4), "SLACS_IX_SumParity_pct": round(slacs_res["par_simple"]*100, 2), "SPARC_MedGal_R2_1to1": round(sparc_res["r2_simple"], 4), "SPARC_RAR_RMS_dex": round(sparc_res["rms_simple"], 4)})

    d = slacs_res["diag"]
    print("\n=== SLACS IX Residual Diagnostics (k = 1.0, Jaffe, Static a_0(0)) ===")
    print(f"  Radial bias (R_E/R_eff) : Spearman rho = {d['rho_rad']:+.3f}, p = {d['p_rad']:.4f}")
    print(f"  Redshift trend (z_L)    : Spearman rho = {d['rho_z']:+.3f}, p = {d['p_z']:.4f}")
    print(f"  Dispersion trend (sigma): Spearman rho = {d['rho_sig']:+.3f}, p = {d['p_sig']:.4f}")
    print(f"  Median transition ratio : y_E = g_bar(R_E)/a_0 = {d['med_y_e']:.2f}")

    dx = slacs_res["dyn_ix"]
    print("\n=== SLACS Dynamic a_0(z_d) = c*H(z_d)/R_4/3 Audit (Jaffe, D_L^2-rescaled) ===")
    print(f"  SLACS IX (N=70, Salp)   : R^2_1:1 = {dx['r2']:.4f} | Lin R^2 = {dx['lin_r2']:.4f} | "
          f"Mean obs/calc = {dx['mean_ratio']:.3f} | Sum Parity = {dx['parity']*100:.1f}% | "
          f"Spearman(z_L) rho = {dx['rho_z']:+.3f} (p = {dx['p_z']:.4f})")
    if slacs_res["dyn_comb"] is not None:
        dc = slacs_res["dyn_comb"]
        print(f"  Combined (N={dc['N']}, Salp/Mid): R^2_1:1 = {dc['r2']:.4f} | Lin R^2 = {dc['lin_r2']:.4f} | "
              f"Mean obs/calc = {dc['mean_ratio']:.3f} | Sum Parity = {dc['parity']*100:.1f}%")

    print("\n=== HST Measurement Noise-Floor Monte Carlo (N_MC = {}) ===".format(args.mc_trials))
    for sig, arr in slacs_res["mc"].items():
        print(f"  Sigma = {sig:.2f} dex -> Expected R^2_1:1 = {np.mean(arr):.4f} +/- {np.std(arr):.4f}")

    out_png = os.path.join(args.data_dir, "SLACS_SPARC_Stiffness_NoiseFloor_v6.png")
    out_csv = os.path.join(args.data_dir, "SLACS_SPARC_Stiffness_NoiseFloor_v6.csv")
    pd.DataFrame(csv_rows).to_csv(out_csv, index=False)
    make_figure(slacs_res, sparc_res, sparc_label, out_png)
    print(f"\nSaved figure : {out_png}")
    print(f"Saved table  : {out_csv}")


if __name__ == "__main__":
    main()
