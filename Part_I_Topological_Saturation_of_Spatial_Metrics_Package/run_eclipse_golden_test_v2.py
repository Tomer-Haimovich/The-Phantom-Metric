#!/usr/bin/env python3
"""
Euclid DR1 ECLIPSE High-Redshift Strong-Lensing Consistency Test (Lambda-TSM vs. MOND & GR)
===========================================================================================
Standalone, reviewer-ready reproduction script (CLI + Google Colab compatible) for the
Euclid DR1 ECLIPSE Golden Subsample (N = 1,095 Grade A lenses with 100% visual-inspection
consensus, non-disc deflectors, external shear gamma < 0.10, and 1.5 <= n_sersic <= 8.0).

Physics & Method:
  For each lens, the 3D stellar mass profile is deprojected from the observed 2D Sersic
  parameters (M_*, R_eff, n_sersic) via the Prugniel-Simien (1997) approximation, and the
  line-of-sight (LOS) gravitational deflection integral is evaluated numerically at R_E
  under four gravity prescriptions:
    1. Baryons_GR  : Pure General Relativity from surviving stellar mass only
    2. MOND_const  : Static RAR/MOND with constant a_0 = 1.20e-10 m s^-2
    3. LTSM_static : Static Lambda-TSM floor a_0(0) = c H_0 / R_{4/3}
    4. LTSM_Hz     : Dynamic Lambda-TSM floor a_0(z_d) = c H(z_d) / R_{4/3}
  where R_{4/3} = (16/3) * (1 + 1 / (16 * pi^2)) ~= 5.36710706.

  Because full spectroscopic source redshifts (z_s) are not yet released for the DR1
  catalog, each model is tested against the exact asymptotic geometric bound:
      beta_req = (D_ds / D_s)_req = c^2 D_d theta_E^2 / (4 G M_eff(<R_E))
      beta_max(z_d) = lim_{z_s -> inf} (D_ds / D_s) = 1 - chi(z_d) / chi(inf)
  where chi(inf) = (c / H_0) * int_0^inf dz / E(z) is evaluated exactly via the
  singularity-free substitution s = (1 + z)^(-1/2) in [0, 1] (< 1e-9 accuracy).
  Any lens with beta_req >= beta_max(z_d) is geometrically impossible under that model.

Usage (external reviewer — zero machine-specific paths):
  1. Place this script and 'eclipse_dr1_golden_1095.csv' in the same folder and run:
       python3 run_eclipse_golden_test_v2_colab.py
  2. Optional flags:
       python3 run_eclipse_golden_test_v2_colab.py --H0 67.36 --Om 0.31391
       python3 run_eclipse_golden_test_v2_colab.py --csv /path/to/eclipse_dr1_golden_1095.csv
       python3 run_eclipse_golden_test_v2_colab.py --raw-dir /path/to/raw_eclipse_dr1_zenodo_dir
       python3 run_eclipse_golden_test_v2_colab.py --export-csv eclipse_dr1_golden_1095.csv

Data Citation:
  Euclid Collaboration: Nersesian et al. (ECLIPSE D, DR1 v20261103);
  Lines et al. (2026, arXiv:2611.15883); Holloway et al. (2026);
  Nightingale et al. (2026, arXiv:2611.15888).
"""

import argparse
import csv
import io
import math
import os
import sys
import zipfile
import numpy as np
from scipy.special import gammainc

# --- Exact SI Physical Constants (IAU / CODATA) ---
C_M_S = 299792458.0
G_SI = 6.67430e-11
M_SUN_KG = 1.98847e30
MPC_M = 3.085677581491367e22
ARCSEC_RAD = math.pi / (180.0 * 3600.0)

# --- Theory Constants ---
R_LTSM = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * math.pi**2))  # ~= 5.36710706
A0_MOND_CONST = 1.20e-10  # m s^-2


def build_cosmology(h0_kms, omega_m):
    """Build exact flat-FLRW distance interpolators on s = (1+z)^(-1/2) in [0, 1].

    The substitution s = (1+z)^(-1/2) maps z in [0, inf) cleanly to s in (0, 1],
    giving int_0^z dz'/E(z') = int_s^1 2 ds' / sqrt(Omega_m + Omega_L * s'^6),
    which is completely smooth and exact all the way to z = inf (s = 0).
    """
    omega_l = 1.0 - omega_m
    h0_si = (h0_kms * 1000.0) / MPC_M

    s_grid = np.linspace(0.0, 1.0, 20001)  # s=0 -> z=inf, s=1 -> z=0
    integrand = 2.0 / np.sqrt(omega_m + omega_l * (s_grid**6))
    ds = s_grid[1] - s_grid[0]
    cum_from_0 = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[:-1] + integrand[1:]) * ds)])
    I_inf = float(cum_from_0[-1])
    I_of_s = I_inf - cum_from_0  # monotonically decreasing from I_inf (s=0) to 0 (s=1)

    # Monotonically increasing arrays for fast inversion I -> s
    I_asc = I_of_s[::-1]
    s_asc = s_grid[::-1]

    def I_of_z(z):
        s = (1.0 + np.asarray(z, dtype=float))**(-0.5)
        return np.interp(s, s_grid, I_of_s)

    def chi_of_z(z):
        return (C_M_S / h0_si) * I_of_z(z)

    def ez_of_z(z):
        return np.sqrt(omega_m * (1.0 + np.asarray(z, dtype=float))**3 + omega_l)

    def dd_of_z(z):
        return chi_of_z(z) / (1.0 + np.asarray(z, dtype=float))

    def beta_max_of_z(zd):
        # Exact asymptotic horizon limit lim_{z_s -> inf} (D_ds / D_s) = 1 - chi(z_d) / chi(inf)
        return 1.0 - I_of_z(zd) / I_inf

    def invert_zs(zd, beta):
        # Exact inversion of beta = (chi_s - chi_d) / chi_s => I_s = I_d / (1 - beta)
        bmax = float(beta_max_of_z(zd))
        if beta <= 0.0 or beta >= bmax:
            return float("nan")
        I_s = float(I_of_z(zd)) / (1.0 - beta)
        if I_s >= I_inf:
            return float("nan")
        s_val = float(np.interp(I_s, I_asc, s_asc))
        if s_val <= 1e-6:
            return float("nan")
        return (1.0 / (s_val**2)) - 1.0

    return h0_si, I_inf, ez_of_z, dd_of_z, beta_max_of_z, invert_zs


def load_from_golden_csv(csv_path):
    """Load pre-selected N=1,095 Golden Subsample CSV."""
    rows = []
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append({
                "uid": r["unique_identifier"],
                "z_d": float(r["z_d"]),
                "has_desi_specz": int(float(r["has_desi_specz"])),
                "logM": float(r["log10_Mstar_Salp"]),
                "theta_E": float(r["theta_E_arcsec"]),
                "R_eff": float(r["R_eff_arcsec"]),
                "n_sersic": float(r["n_sersic"]),
                "gamma_ext": float(r["gamma_ext"]),
            })
    return rows


def load_from_raw_eclipse_dir(raw_dir):
    """Rebuild the Golden Subsample directly from the raw Zenodo ECLIPSE DR1 catalogs."""
    key_path = os.path.join(raw_dir, "eclipse_dr1_catalogue_key_v20261103.csv")
    obj_path = os.path.join(raw_dir, "ECLIPSE_object_catalogue", "ECLIPSE_object_catalogue.csv")
    zip_path = os.path.join(raw_dir, "ECLIPSE_lens_models", "catalogue_csvs.zip")

    with open(key_path, "r", encoding="utf-8") as f:
        key = {r["unique_identifier"]: r for r in csv.DictReader(f)}
    with open(obj_path, "r", encoding="utf-8") as f:
        obj = {r["Unique_Identifier"]: r for r in csv.DictReader(f)}

    shear_map = {}
    with zipfile.ZipFile(zip_path) as zf:
        with io.TextIOWrapper(zf.open("dr1_catalogue_v5_full/csv/lens_mass.csv")) as f:
            for r in csv.DictReader(f):
                try:
                    shear_map[r["unique_identifier"]] = math.hypot(
                        float(r["shear_gamma_1"]), float(r["shear_gamma_2"])
                    )
                except ValueError:
                    pass

    rows = []
    for uid, kr in key.items():
        if (
            kr["final_sample"] != "True"
            or kr["lens_grade"] != "A"
            or kr["lens_model_consensus_100"] != "Success"
        ):
            continue
        orow = obj.get(uid)
        if not orow or orow["disc_deflector"] != "False" or orow["large_system_warning"] != "False":
            continue
        gamma = shear_map.get(uid, 99.0)
        if gamma >= 0.10:
            continue
        try:
            raw_spec = orow["desi_deflector_specz"]
            z_spec = float(raw_spec) if raw_spec not in ("", "nan", "None") else -1.0
            zd = z_spec if z_spec > 0 else float(kr["eclipse_median_phz"])
            logm = float(kr["eclipse_log_survived_mass_med"])
            theta_e = float(kr["eclipse_effective_einstein_radius"])
            reff = float(kr["eclipse_effective_radius"])
            n_ser = float(kr["eclipse_sersic_index"])
        except ValueError:
            continue
        if not (
            0.05 < zd < 1.5
            and 9.5 < logm < 13.0
            and 0.2 < theta_e < 4.0
            and 0.05 < reff < 10.0
            and 1.5 <= n_ser <= 8.0
        ):
            continue
        rows.append({
            "uid": uid,
            "z_d": zd,
            "has_desi_specz": int(z_spec > 0),
            "logM": logm,
            "theta_E": theta_e,
            "R_eff": reff,
            "n_sersic": n_ser,
            "gamma_ext": gamma,
        })
    return rows


def nu_rar(g_n, a0):
    """McGaugh / Lambda-TSM Radial Acceleration Relation interpolation function."""
    y = np.maximum(g_n / a0, 1e-12)
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))


def evaluate_models(rows, h0_kms=67.36, omega_m=0.31391, export_csv=None):
    """Compute exact LOS geodesic lensing integrals and print summary tables."""
    h0_si, I_inf, ez_of_z, dd_of_z, beta_max_of_z, invert_zs = build_cosmology(h0_kms, omega_m)

    zd = np.array([r["z_d"] for r in rows])
    logm = np.array([r["logM"] for r in rows])
    the = np.array([r["theta_E"] for r in rows])
    reff = np.array([r["R_eff"] for r in rows])
    n_ser = np.array([r["n_sersic"] for r in rows])
    n_spec = sum(r["has_desi_specz"] for r in rows)

    a0_ltsm_0 = (C_M_S * h0_si) / R_LTSM
    print("====================================================================================")
    print("  EUCLID DR1 ECLIPSE GOLDEN SUBSAMPLE — LAMBDA-TSM vs. MOND & GR CONSISTENCY TEST")
    print("====================================================================================")
    print(f"Cosmology : H0 = {h0_kms:.2f} km/s/Mpc, Omega_m = {omega_m:.5f}, Omega_L = {1.0 - omega_m:.5f}, I(inf) = {I_inf:.6f}")
    print(f"Constants : R_4/3 = {R_LTSM:.8f}, a0_LTSM(0) = {a0_ltsm_0:.4e} m/s^2, a0_MOND = {A0_MOND_CONST:.2e} m/s^2")
    print(f"Sample    : N = {len(rows)} (Grade A, final_sample, 100% consensus, non-disc, gamma < 0.10)")
    print(f"            DESI spec-z = {n_spec} | z_d range = [{np.min(zd):.2f}, {np.max(zd):.2f}] | med(z_d) = {np.median(zd):.3f}")
    print(f"            med(theta_E) = {np.median(the):.3f}\" | med(R_eff) = {np.median(reff):.3f}\" | med(n) = {np.median(n_ser):.2f}\n")

    # Vectorized Prugniel-Simien 3D deprojection + LOS integral over u in [0, 9]
    u_grid = np.linspace(0.0, 9.0, 91)
    w_u = np.ones_like(u_grid) * (u_grid[1] - u_grid[0])
    w_u[0] *= 0.5
    w_u[-1] *= 0.5
    cosh_u = np.cosh(u_grid)[None, :]

    dd_m = dd_of_z(zd)
    re_m = dd_m * the * ARCSEC_RAD
    reff_m = dd_m * reff * ARCSEC_RAD
    bmax = beta_max_of_z(zd)

    bn = 2.0 * n_ser - 1.0 / 3.0 + 4.0 / (405.0 * n_ser) + 46.0 / (25515.0 * n_ser**2)
    pn = 1.0 - 0.6097 / n_ser + 0.05563 / (n_ser**2)
    a_gam = n_ser * (3.0 - pn)

    r_grid = re_m[:, None] * cosh_u
    x_arg = bn[:, None] * (r_grid / reff_m[:, None])**(1.0 / n_ser[:, None])
    f_enc = gammainc(a_gam[:, None], x_arg)

    def run_imf_suite(imf_label, dlogm):
        mstar_kg = (10.0**(logm + dlogm)) * M_SUN_KG
        g_n = G_SI * (mstar_kg[:, None] * f_enc) / (r_grid**2)
        a0_ltsm_z = a0_ltsm_0 * ez_of_z(zd)

        models = {
            "Baryons_GR":  g_n,
            "MOND_const":  g_n * nu_rar(g_n, A0_MOND_CONST),
            "LTSM_static": g_n * nu_rar(g_n, a0_ltsm_0),
            "LTSM_Hz":     g_n * nu_rar(g_n, a0_ltsm_z[:, None]),
        }

        betas = {}
        for mname, g_tot in models.items():
            m_eff = (re_m**2 / G_SI) * np.sum(g_tot * w_u[None, :], axis=1)
            betas[mname] = (C_M_S**2 * dd_m * (the * ARCSEC_RAD)**2) / (4.0 * G_SI * m_eff)

        print(f"================ IMF: {imf_label} (dlogM = {dlogm:+.2f} dex) ================")
        bins = [
            ("z_d [0.10, 0.45)", 0.10, 0.45),
            ("z_d [0.45, 0.70)", 0.45, 0.70),
            ("z_d [0.70, 0.95)", 0.70, 0.95),
            ("z_d [0.95, 1.50)", 0.95, 1.50),
            ("ALL (0.10-1.50)",  0.10, 1.50),
        ]
        for tag, b_lo, b_hi in bins:
            mask = (zd >= b_lo) & (zd < b_hi)
            nb = int(np.sum(mask))
            print(f"\n--- {tag} | N = {nb} | med(beta_max) = {np.median(bmax[mask]):.3f} ---")
            print(f'{"Model":12s} | {"med(beta_req)":>13s} | {"% Unphys (>=bmax)":>17s} | {"med(z_s,imp)":>12s} | {"IQR(z_s)":>14s}')
            unphys_map = {}
            for mname in ("Baryons_GR", "MOND_const", "LTSM_static", "LTSM_Hz"):
                bm = betas[mname][mask]
                zm = zd[mask]
                unphys = 100.0 * float(np.mean(bm >= bmax[mask]))
                unphys_map[mname] = unphys
                zs_valid = [invert_zs(z_i, b_i) for z_i, b_i in zip(zm, bm)]
                zs_valid = [z for z in zs_valid if not math.isnan(z)]
                med_zs = float(np.median(zs_valid)) if zs_valid else float("nan")
                q25, q75 = (np.percentile(zs_valid, [25, 75]) if len(zs_valid) >= 4 else (float("nan"), float("nan")))
                print(f"{mname:12s} | {np.median(bm):13.3f} | {unphys:16.1f}% | {med_zs:12.2f} | [{q25:.2f}, {q75:.2f}]")
            if unphys_map["MOND_const"] > 0:
                rel_red = 100.0 * (unphys_map["MOND_const"] - unphys_map["LTSM_Hz"]) / unphys_map["MOND_const"]
                print(f"  -> LTSM_Hz vs MOND_const failure reduction: -{unphys_map['MOND_const'] - unphys_map['LTSM_Hz']:.1f} pp ({rel_red:.1f}% relative)")
        return betas

    betas_salp = run_imf_suite("Salpeter (ECLIPSE native)", 0.0)
    print()
    run_imf_suite("Chabrier (-0.22 dex)", -0.22)

    # Low-z diagnostic split (galaxy-scale vs group-scale massive deflectors)
    print("\n================ LOW-Z DIAGNOSTIC (z_d in [0.10, 0.45), LTSM_Hz, Salpeter) ================")
    low_m = (zd >= 0.10) & (zd < 0.45)
    b_l = betas_salp["LTSM_Hz"]
    splits = [
        ('theta_E < 0.9" (galaxy-scale)',    low_m & (the < 0.9)),
        ('0.9" <= theta_E < 1.3"',           low_m & (the >= 0.9) & (the < 1.3)),
        ('theta_E >= 1.3" (group/massive)',  low_m & (the >= 1.3)),
        ('R_eff < 0.8"',                     low_m & (reff < 0.8)),
        ('0.8" <= R_eff < 1.5"',             low_m & (reff >= 0.8) & (reff < 1.5)),
        ('R_eff >= 1.5"',                    low_m & (reff >= 1.5)),
    ]
    for label, sm in splits:
        print(f"{label:32s} | N={int(np.sum(sm)):3d} | med(logM)={np.median(logm[sm]):.2f} | "
              f"med(beta_req)={np.median(b_l[sm]):.3f} | Unphys={100.0 * np.mean(b_l[sm] >= bmax[sm]):5.1f}%")

    if export_csv:
        def fmt_zs(z_val):
            return "" if math.isnan(z_val) else round(z_val, 4)

        with open(export_csv, "w", newline="", encoding="utf-8") as f:
            wr = csv.writer(f)
            wr.writerow([
                "unique_identifier", "z_d", "has_desi_specz", "log10_Mstar_Salp",
                "theta_E_arcsec", "R_eff_arcsec", "n_sersic", "gamma_ext",
                "Dd_Mpc", "beta_max", "beta_req_GR", "beta_req_MOND", "beta_req_LTSM",
                "zs_imp_GR", "zs_imp_MOND", "zs_imp_LTSM",
            ])
            for i, r in enumerate(rows):
                wr.writerow([
                    r["uid"], round(r["z_d"], 4), r["has_desi_specz"], round(r["logM"], 4),
                    round(r["theta_E"], 4), round(r["R_eff"], 4), round(r["n_sersic"], 3), round(r["gamma_ext"], 4),
                    round(float(dd_m[i] / MPC_M), 2), round(float(bmax[i]), 4),
                    round(float(betas_salp["Baryons_GR"][i]), 4),
                    round(float(betas_salp["MOND_const"][i]), 4),
                    round(float(betas_salp["LTSM_Hz"][i]), 4),
                    fmt_zs(invert_zs(r["z_d"], betas_salp["Baryons_GR"][i])),
                    fmt_zs(invert_zs(r["z_d"], betas_salp["MOND_const"][i])),
                    fmt_zs(invert_zs(r["z_d"], betas_salp["LTSM_Hz"][i])),
                ])
        print(f"\nExported updated Golden Subsample CSV -> {export_csv}")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    default_csv = os.path.join(script_dir, "eclipse_dr1_golden_1095.csv")
    if not os.path.exists(default_csv):
        default_csv = "eclipse_dr1_golden_1095.csv"

    parser = argparse.ArgumentParser(
        description="Reproduce the Euclid DR1 ECLIPSE Golden Subsample (N=1,095) test for Lambda-TSM vs. MOND & GR."
    )
    parser.add_argument("--csv", default=default_csv, help="Path to eclipse_dr1_golden_1095.csv")
    parser.add_argument("--raw-dir", default=None, help="Optional: rebuild from raw Zenodo ECLIPSE DR1 directory")
    parser.add_argument("--H0", type=float, default=67.36, help="Hubble constant H0 [km/s/Mpc] (default: 67.36)")
    parser.add_argument("--Om", type=float, default=0.31391, help="Flat matter density Omega_m,0 (default: 0.31391)")
    parser.add_argument("--export-csv", default=None, help="Optional: export updated Golden Subsample CSV")
    args, _ = parser.parse_known_args()

    if args.raw_dir:
        rows = load_from_raw_eclipse_dir(args.raw_dir)
    else:
        if not os.path.exists(args.csv):
            sys.exit(
                f"ERROR: Could not find '{args.csv}'. Place 'eclipse_dr1_golden_1095.csv' next to this "
                "script, pass --csv <path>, or pass --raw-dir <raw_eclipse_dr1_folder>."
            )
        rows = load_from_golden_csv(args.csv)

    evaluate_models(rows, h0_kms=args.H0, omega_m=args.Om, export_csv=args.export_csv)


if __name__ == "__main__":
    main()
