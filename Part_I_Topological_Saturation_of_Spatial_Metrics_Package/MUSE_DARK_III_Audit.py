#!/usr/bin/env python3
"""
MUSE_DARK_III_Audit.py
Zero-parameter empirical audit of the LambdaTSM redshift-running acceleration law
    a0(z) = [c * H(z)] / R_{4/3}^circ = a0(0) * sqrt(Omega_m0 * (1 + z)^3 + Omega_Lambda)
against the intermediate-redshift MUSE-DARK III RAR measurements (0.33 < z < 1.44; N = 79 SFGs).

Sources:
  - Ciocan et al. (2026), A&A 709, L16 / arXiv:2604.22613v1 (MUSE-DARK III)
  - B.-I. Ciocan tabulated 4-bin DC14 values published in Logue (2026), MNRAS Table 2
  - McGaugh, Lelli & Schombert (2016), PRL 117, 201101 (SPARC z = 0 anchor)

Outputs:
  - MUSE_DARK_III_Audit_Results.csv
  - MUSE_DARK_III_a0z_Figure.png (300 DPI publication-ready 2-panel figure)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_dist


# ------------------------------------------------------------------------------
# 1. Unified LambdaTSM Cosmological Baseline (Part I & Part II v8)
# ------------------------------------------------------------------------------
C_SI     = 299792458.0
MPC_M    = 3.085677581491367e22
H0_KMS   = 67.36
H0_SI    = H0_KMS * 1e3 / MPC_M
R43      = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * np.pi**2))  # 5.36710706...
OM0      = 0.31391
OL0      = 1.0 - OM0

A0_LTSM_0 = (C_SI * H0_SI / R43) * 1e10           # 1.21936 (in 10^-10 m s^-2)
A0_2PI_0  = (C_SI * H0_SI / (2.0 * np.pi)) * 1e10 # 1.04158 (in 10^-10 m s^-2)
A0_MOND_0 = 1.2000                                # Empirical SPARC


def E_z(z):
    return np.sqrt(OM0 * (1.0 + z)**3 + OL0)


# ------------------------------------------------------------------------------
# 2. MUSE-DARK III 4-Bin Data (Ciocan et al. 2026 Fig. 3 / Logue 2026 Table 2)
# ------------------------------------------------------------------------------
Z_LO   = np.array([0.333, 0.671, 0.980, 1.112])
Z_HI   = np.array([0.671, 0.980, 1.112, 1.438])
Z_MEAN = np.array([0.578, 0.846, 1.058, 1.319])
Z_MID  = 0.5 * (Z_LO + Z_HI)
A0_OBS = np.array([1.993, 2.201, 2.567, 2.715])   # x10^-10 m s^-2 (posterior medians)
A0_ERR = np.array([0.081, 0.105, 0.115, 0.136])   # 1-sigma posterior std

LOG_OBS = np.log10(A0_OBS)
LOG_ERR = A0_ERR / (A0_OBS * np.log(10.0))


def evaluate_model(name, pred_vals, z_eval, sigma_sys=0.10):
    log_pred = np.log10(pred_vals)
    r = LOG_OBS - log_pred
    w = 1.0 / (LOG_ERR**2)

    # (a) Independent bins at model's own normalization (4 dof)
    chi2_own = float(np.sum(w * r**2))
    p_own = float(1.0 - chi2_dist.cdf(chi2_own, 4))

    # (b) Free amplitude / shape-only test (3 dof)
    delta_free = float(np.sum(w * r) / np.sum(w))
    chi2_shape = float(np.sum(w * (r - delta_free)**2))
    p_shape = float(1.0 - chi2_dist.cdf(chi2_shape, 3))

    # (c) Shared calibration offset prior (0.10 dex, fully correlated across bins, 4 dof)
    C10 = np.diag(LOG_ERR**2) + (sigma_sys**2) * np.ones((4, 4))
    chi2_sh10 = float(r @ np.linalg.inv(C10) @ r)
    p_sh10 = float(1.0 - chi2_dist.cdf(chi2_sh10, 4))
    delta_sh10 = float(np.sum(w * r) / (np.sum(w) + 1.0 / (sigma_sys**2)))

    # (d) Appendix C molecular-gas systematic budget (0.20 dex shared prior, 4 dof)
    C20 = np.diag(LOG_ERR**2) + (0.20**2) * np.ones((4, 4))
    chi2_sh20 = float(r @ np.linalg.inv(C20) @ r)
    p_sh20 = float(1.0 - chi2_dist.cdf(chi2_sh20, 4))

    is_dyn = ("Dynamic" in name)
    return {
        "Model": name,
        "a0_z0_1e10": round(float(pred_vals[0] / E_z(z_eval[0]) if is_dyn else pred_vals[0]), 4),
        "chi2_own_4dof": round(chi2_own, 3),
        "p_own": p_own,
        "chi2_shape_3dof": round(chi2_shape, 3),
        "p_shape": round(p_shape, 4),
        "delta_free_dex": round(delta_free, 4),
        "chi2_shared_010_4dof": round(chi2_sh10, 3),
        "p_shared_010": round(p_sh10, 4),
        "delta_shared_010_dex": round(delta_sh10, 4),
        "chi2_shared_020_4dof": round(chi2_sh20, 3),
        "p_shared_020": round(p_sh20, 4),
    }


def main():
    models = [
        ("LambdaTSM Dynamic: cH(z)/R43",    A0_LTSM_0 * E_z(Z_MEAN)),
        ("MOND Dynamic: 1.20 * E(z)",       A0_MOND_0 * E_z(Z_MEAN)),
        ("Legacy 2pi Dynamic: cH(z)/(2pi)", A0_2PI_0  * E_z(Z_MEAN)),
        ("LambdaTSM Static: cH0/R43",       np.full(4, A0_LTSM_0)),
        ("MOND Static: 1.20",               np.full(4, A0_MOND_0)),
    ]

    df = pd.DataFrame([evaluate_model(n, v, Z_MEAN) for n, v in models])
    df.to_csv("MUSE_DARK_III_Audit_Results.csv", index=False)

    print("=== MUSE-DARK III 4-BIN ZERO-PARAMETER AUDIT (at <z>) ===")
    print(df.to_string(index=False))

    # Linearized tangent over 0.333 < z < 1.438 vs. Appendix E Eq. (E.5)
    z_rng = np.linspace(0.333, 1.438, 1000)
    a1_ltsm, a0_lin_ltsm = np.polyfit(z_rng, A0_LTSM_0 * E_z(z_rng), 1)
    print("\n=== SELF-CONSISTENT MOND/RAR COMPARISON (MUSE-DARK III Appendix E) ===")
    print(f"  LambdaTSM a0(z=1.00) = {A0_LTSM_0 * E_z(1.0):.3f} vs. App. E Eq.(E.4) 2.19 (+0.12/-0.10, 95% CI)")
    print(f"  LambdaTSM tangent on [0.333, 1.438]: a0_lin(0) = {a0_lin_ltsm:.3f}, a1 = {a1_ltsm:.3f}")
    print(f"  App. E Eq.(E.5) MOND RAR linear fit: a0(0)     = 1.03 +/- 0.05, a1 = 1.20 +/- 0.10 (95% CI)")

    # --------------------------------------------------------------------------
    # Publication Figure (300 DPI)
    # --------------------------------------------------------------------------
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.linewidth": 1.0,
    })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 5.0), dpi=300)

    # Panel A: Redshift evolution a0(z)
    z_plot = np.linspace(0.0, 1.50, 400)
    y_ltsm = A0_LTSM_0 * E_z(z_plot)
    y_2pi  = A0_2PI_0  * E_z(z_plot)
    y_stat = np.full_like(z_plot, A0_MOND_0)

    z_muse = np.linspace(0.333, 1.438, 200)
    y_mond_e5    = 1.03 + 1.20 * z_muse
    y_mond_e5_hi = (1.03 + 0.025) + (1.20 + 0.05) * z_muse
    y_mond_e5_lo = (1.03 - 0.025) + (1.20 - 0.05) * z_muse
    y_dc14_e4    = 1.00 + 1.59 * z_muse

    ax1.fill_between(Z_MEAN, A0_OBS * 10**(-0.10), A0_OBS * 10**(0.10),
                     color="#95a5a6", alpha=0.18,
                     label=r"Shared $\pm 0.10$ dex cal./gas band (App. C)")

    ax1.fill_between(z_muse, y_mond_e5_lo, y_mond_e5_hi, color="#2980b9", alpha=0.22)
    ax1.plot(z_muse, y_mond_e5, color="#2980b9", ls="-", lw=2.0,
             label=r"MUSE-DARK III App. E (MOND RAR, $a_1=1.20$)")
    ax1.plot(z_muse, y_dc14_e4, color="#8e44ad", ls=":", lw=1.7,
             label=r"MUSE-DARK III Main (DC14, $a_1=1.59$)")

    ax1.plot(z_plot, y_ltsm, color="#0f766e", lw=2.8, zorder=6,
             label=r"$\Lambda$TSM Dynamic: $cH(z)/\mathcal{R}_{4/3}^{\circ}$ (0 free params)")
    ax1.plot(z_plot, y_2pi, color="#d97706", ls="-.", lw=1.6, zorder=4,
             label=r"Legacy $2\pi$ Dynamic: $cH(z)/(2\pi)$")
    ax1.plot(z_plot, y_stat, color="#b91c1c", ls="--", lw=1.8, zorder=4,
             label=r"Static MOND ($a_0 = 1.20 \times 10^{-10}\ \mathrm{m\,s^{-2}}$)")

    ax1.errorbar([0.0], [1.20], yerr=[0.24], fmt="*", ms=11, color="#dc2626",
                 ecolor="#dc2626", elinewidth=1.4, capsize=3.5, zorder=8,
                 label=r"SPARC $z=0$ ($1.20 \pm 0.24$ syst.)")
    ax1.errorbar([0.98], [2.19], yerr=[[0.05], [0.06]], fmt="D", ms=6.5, color="#1d4ed8",
                 ecolor="#1d4ed8", elinewidth=1.5, capsize=3.5, zorder=9,
                 label=r"App. E MOND RAR ($z\sim 1$: $2.19^{+0.12}_{-0.10}$ 95% CI)")
    ax1.errorbar(Z_MEAN, A0_OBS, yerr=A0_ERR, xerr=[Z_MEAN - Z_LO, Z_HI - Z_MEAN],
                 fmt="o", ms=6.0, color="#111827", ecolor="#374151", elinewidth=1.3,
                 capsize=3.0, zorder=7, label="MUSE-DARK III 4-Bin DC14 (Ciocan+26)")

    ax1.set_xlim(-0.06, 1.50)
    ax1.set_ylim(0.85, 3.45)
    ax1.set_xlabel("Redshift $z$", fontsize=10.5)
    ax1.set_ylabel(r"Characteristic Acceleration $a_0(z)\ [10^{-10}\ \mathrm{m\,s^{-2}}]$", fontsize=10.5)
    ax1.set_title(r"A) Redshift Evolution: SPARC ($z=0$) to MUSE-DARK III ($z\leq 1.44$)",
                  fontsize=10.5, fontweight="bold")
    ax1.grid(True, ls=":", alpha=0.45)
    ax1.legend(loc="upper left", fontsize=7.4, frameon=True, framealpha=0.93)

    # Panel B: 4-bin goodness-of-fit comparison (clean, zero text overlap)
    labels_b = [
        r"$\Lambda$TSM Dyn." + "\n" + r"$cH(z)/\mathcal{R}_{4/3}^{\circ}$",
        r"MOND$\times E(z)$" + "\n" + r"$1.20\,E(z)$",
        r"Legacy $2\pi$ Dyn." + "\n" + r"$cH(z)/(2\pi)$",
        r"$\Lambda$TSM Static" + "\n" + r"$1.219$",
        r"MOND Static" + "\n" + r"$1.200$",
    ]
    chi2_vals = df["chi2_shared_010_4dof"].values
    p_vals    = df["p_shared_010"].values
    colors_b  = ["#0f766e", "#2563eb", "#d97706", "#9ca3af", "#b91c1c"]

    x_pos = np.arange(len(labels_b))
    bars = ax2.bar(x_pos, chi2_vals, width=0.56, color=colors_b,
                   edgecolor="#1f2937", linewidth=0.9, zorder=3)

    ax2.axhline(9.49, color="#374151", ls="--", lw=1.25, zorder=2,
                label=r"$p = 0.05$ rejection threshold ($\chi^2 = 9.49$, 4 dof)")

    for i, (bar, cval, pval) in enumerate(zip(bars, chi2_vals, p_vals)):
        ptxt = f"$p={pval:.3f}$" if pval >= 0.001 else r"$p<10^{-6}$"
        # Place short-bar labels neatly above the bar (for i=2, sit cleanly just above 9.49 or inside-white halo)
        y_txt = cval + 0.9 if cval > 10 else cval + 0.65
        ax2.text(bar.get_x() + bar.get_width() / 2, y_txt,
                 f"$\\chi^2={cval:.2f}$\n({ptxt})",
                 ha="center", va="bottom", fontsize=7.9, fontweight="bold",
                 color="#111827",
                 bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.82),
                 zorder=5)

    box_txt = (
        "Self-Consistent MOND RAR (App. E):\n"
        r"• Obs $a_0(z\sim 1) = 2.19^{+0.12}_{-0.10}$  vs.  $\Lambda$TSM $2.180$ ($0.4\%$)" + "\n"
        r"• Obs slope $a_1 = 1.20\pm 0.10$  vs.  $\Lambda$TSM $1.212$ ($1.0\%$)" + "\n"
        r"• $\Delta\chi^2$ vs. Static MOND: $-33.9$ (shared) / $-841$ (own)"
    )
    ax2.text(0.03, 0.96, box_txt, transform=ax2.transAxes, va="top", ha="left",
             fontsize=8.0,
             bbox=dict(boxstyle="round,pad=0.42", facecolor="#f0fdfa", edgecolor="#0f766e", alpha=0.96),
             zorder=6)

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(labels_b, fontsize=8.3)
    ax2.set_ylim(0, 47.0)
    ax2.set_ylabel(r"$\chi^2$ (4-Bin Shared $0.10$ dex Calibration Prior, 4 dof)", fontsize=10.0)
    ax2.set_title(r"B) MUSE-DARK III 4-Bin Goodness-of-Fit (0 Free Parameters)",
                  fontsize=10.5, fontweight="bold")
    ax2.grid(True, axis="y", ls=":", alpha=0.45, zorder=0)
    ax2.legend(loc="upper left", bbox_to_anchor=(0.02, 0.72), fontsize=8.0, frameon=True, framealpha=0.95)

    plt.tight_layout(pad=1.2)
    plt.savefig("MUSE_DARK_III_a0z_Figure.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("\nSaved MUSE_DARK_III_Audit_Results.csv and MUSE_DARK_III_a0z_Figure.png.")


if __name__ == "__main__":
    main()
