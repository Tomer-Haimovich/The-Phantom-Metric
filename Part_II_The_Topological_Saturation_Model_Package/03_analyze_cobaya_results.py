import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from getdist.mcsamples import loadMCSamples
from getdist import plots

os.makedirs("results", exist_ok=True)

R_PARTITION = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * math.pi**2))  # 5.36710715...

RUNS = {
    "ltsm_p18": ("chains/ltsm_planck_only",     "LambdaTSM (k=5) | Planck 2018", 5),
    "lcdm_p18": ("chains/lcdm_planck_only",     "LambdaCDM (k=6) | Planck 2018", 6),
    "ltsm_all": ("chains/ltsm_planck_desi_sn",  "LambdaTSM (k=5) | Planck+DESI+SNe", 5),
    "lcdm_all": ("chains/lcdm_planck_desi_sn",  "LambdaCDM (k=6) | Planck+DESI+SNe", 6),
}

# Effective data-point counts for BIC reference:
#   Planck 2018 lite + sroll2 + lensing ~ 623 bandpowers (N_p18)
#   + DESI 2024 BAO (12) + Pantheon+ (1590) ~ 2225 data points (N_all)
N_EFF = {"p18": 623, "all": 2225}

samples = {}
stats = {}

def load_chain(prefix):
    if not os.path.exists(prefix + ".1.txt") and not os.path.exists(prefix + ".txt"):
        return None
    return loadMCSamples(prefix, settings={"ignore_rows": 0.3})

for key, (prefix, label, k_cosmo) in RUNS.items():
    s = load_chain(prefix)
    if s is None:
        print(f"[SKIP] Chain not found yet: {prefix}")
        continue
    samples[key] = s
    mstats = s.getMargeStats()
    pnames = [p.name for p in s.paramNames.names]
    if "chi2" in pnames:
        chi2_min = float(np.min(s.samples[:, pnames.index("chi2")]))
    else:
        chi2_min = 2.0 * float(np.min(s.loglikes))

    row = {"label": label, "k": k_cosmo, "chi2_min": chi2_min, "N_samples": len(s.samples)}
    for par in ["omega_b", "omega_cdm", "H0", "logA", "n_s", "tau_reio", "Omega_m", "sigma8", "S8", "rs_drag", "dm_baryon_ratio", "a0_unified_1e10"]:
        if par in pnames:
            p = mstats.parWithName(par)
            row[par] = (p.mean, p.err)
    stats[key] = row

def fmt(row, par, dec=4):
    if par not in row:
        return "---"
    m, e = row[par]
    return f"{m:.{dec}f} +/- {e:.{dec}f}"

lines_out = []
lines_out.append("=" * 108)
lines_out.append(" REAL-SKY COBAYA MCMC AUDIT REPORT — LambdaTSM (k=5) vs. LambdaCDM (k=6)")
lines_out.append(f" Structural Partition Lock: R_{{4/3}}^o = (16/3)*(1 + 1/(16*pi^2)) = {R_PARTITION:.6f}")
lines_out.append("=" * 108)

for tag, k_ltsm, k_lcdm, n_key, title in [
    ("Planck 2018 Only", "ltsm_p18", "lcdm_p18", "p18", "A) PLANCK 2018 (Plik-lite + SROLL2 EE + Lensing)"),
    ("Planck + DESI + Pantheon+", "ltsm_all", "lcdm_all", "all", "B) PLANCK 2018 + DESI 2024 BAO + PANTHEON+ SNe Ia"),
]:
    if k_ltsm not in stats or k_lcdm not in stats:
        continue
    rl = stats[k_ltsm]
    rc = stats[k_lcdm]
    dchi2 = rl["chi2_min"] - rc["chi2_min"]
    daic = dchi2 + 2.0 * (5 - 6)
    dbic = dchi2 + (5 - 6) * math.log(N_EFF[n_key])
    rat_m, rat_e = rc["dm_baryon_ratio"]
    z_ratio = (R_PARTITION - rat_m) / rat_e

    lines_out.append(f"\n--- {title} ---")
    lines_out.append(f"  Best-fit chi2 : LambdaTSM = {rl['chi2_min']:.2f} | LambdaCDM = {rc['chi2_min']:.2f} | Delta_chi2 = {dchi2:+.3f}")
    lines_out.append(f"  Model Select  : Delta_AIC = {daic:+.2f} | Delta_BIC (N={N_EFF[n_key]}) = {dbic:+.2f}  [negative favors LambdaTSM]")
    lines_out.append(f"  Ratio Test    : LambdaCDM omega_c/omega_b = {rat_m:.4f} +/- {rat_e:.4f} vs R_{{4/3}}^o = {R_PARTITION:.4f} ({z_ratio:+.2f} sigma)")
    lines_out.append(f"  {'Parameter':18s} | {'LambdaTSM (k=5)':24s} | {'LambdaCDM (k=6)':24s}")
    lines_out.append("  " + "-" * 72)
    for par, dec in [
        ("omega_b", 5), ("omega_cdm", 5), ("H0", 2), ("Omega_m", 4),
        ("logA", 4), ("n_s", 4), ("tau_reio", 4), ("sigma8", 4),
        ("S8", 4), ("rs_drag", 2), ("dm_baryon_ratio", 4), ("a0_unified_1e10", 4)
    ]:
        lines_out.append(f"  {par:18s} | {fmt(rl, par, dec):24s} | {fmt(rc, par, dec):24s}")

lines_out.append("\n--- READY-TO-PASTE LATEX TABLE BODY FOR PAPER 2 ---")
if all(k in stats for k in ["ltsm_p18", "lcdm_p18", "ltsm_all", "lcdm_all"]):
    for par, tex, dec in [
        ("omega_b", r"$\omega_b$", 5),
        ("omega_cdm", r"$\omega_{\text{cdm}}$ / $\omega_{\text{phantom}}$", 5),
        ("H0", r"$H_0\ [\text{km s}^{-1}\text{ Mpc}^{-1}]$", 2),
        ("Omega_m", r"$\Omega_m$", 4),
        ("n_s", r"$n_s$", 4),
        ("sigma8", r"$\sigma_8$ ($\Sigma m_\nu=0.06\text{ eV}$)", 4),
        ("S8", r"$S_8$", 4),
        ("rs_drag", r"$r_d\ [\text{Mpc}]$", 2),
        ("a0_unified_1e10", r"$a_0 = cH_0/\mathcal{R}_{4/3}^{\circ}\ [10^{-10}\text{ m s}^{-2}]$", 4),
    ]:
        cells = []
        for k in ["lcdm_p18", "ltsm_p18", "lcdm_all", "ltsm_all"]:
            m, e = stats[k][par]
            cells.append(f"${m:.{dec}f} \\pm {e:.{dec}f}$")
        lines_out.append(f"{tex} & " + " & ".join(cells) + r" \\")

report_txt = "\n".join(lines_out)
print(report_txt)
with open("results/cobaya_audit_report.txt", "w", encoding="utf-8") as f:
    f.write(report_txt + "\n")

if "ltsm_all" in samples and "lcdm_all" in samples:
    g = plots.get_subplot_plotter(width_inch=9.5)
    g.settings.axes_fontsize = 10
    g.settings.lab_fontsize = 11
    g.settings.legend_fontsize = 10
    g.triangle_plot(
        [samples["lcdm_all"], samples["ltsm_all"]],
        ["omega_b", "H0", "Omega_m", "n_s", "sigma8", "a0_unified_1e10"],
        filled=True,
        legend_labels=[
            r"$\Lambda\mathrm{CDM}$ ($k=6$, Planck+DESI+Pantheon+)",
            r"$\Lambda\mathrm{TSM}$ ($k=5$, $\omega_{\mathrm{ph}}=\mathcal{R}_{4/3}^\circ\omega_b$)",
        ],
        contour_colors=["#64748b", "#0f766e"],
        line_args=[{"lw": 1.4, "color": "#475569"}, {"lw": 1.8, "color": "#0f766e"}],
    )
    g.export("results/Cobaya_Planck_DESI_SNe_Corner.png", dpi=300)
    print("\nSaved results/Cobaya_Planck_DESI_SNe_Corner.png")

if "lcdm_p18" in samples or "lcdm_all" in samples:
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=300)
    for key, col, lbl in [
        ("lcdm_p18", "#0284c7", r"$\Lambda\mathrm{CDM}$: Planck 2018 Only"),
        ("lcdm_all", "#7c3aed", r"$\Lambda\mathrm{CDM}$: Planck 2018 + DESI BAO + Pantheon+"),
    ]:
        if key not in samples:
            continue
        dens = samples[key].get1DDensity("dm_baryon_ratio")
        xs = np.linspace(dens.x.min(), dens.x.max(), 400)
        ys = dens.Prob(xs)
        m, e = stats[key]["dm_baryon_ratio"]
        z = (R_PARTITION - m) / e
        ax.plot(xs, ys, color=col, lw=2.3, label=f"{lbl}\n(${m:.4f} \\pm {e:.4f}$, diff ${z:+.2f}\\sigma$)")
        ax.fill_between(xs, 0, ys, color=col, alpha=0.16)

    ax.axvline(R_PARTITION, color="#be123c", ls="--", lw=2.2,
               label=r"$\Lambda\mathrm{TSM}$ Rigid Lock: $\mathcal{R}_{4/3}^\circ = 5.3671$")
    ax.set_xlabel(r"Dark-to-Baryon Density Ratio $\omega_{\mathrm{cdm}} / \omega_b$", fontsize=11)
    ax.set_ylabel("Normalized Posterior Probability", fontsize=11)
    ax.set_title(r"Direct Observational Test of the $\Lambda\mathrm{TSM}$ Partition Lock", fontsize=11.5, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left", fontsize=8.8, framealpha=0.95)
    ax.grid(True, ls=":", alpha=0.5)
    plt.tight_layout()
    plt.savefig("results/Cobaya_Ratio_Falsifiability.png", dpi=300)
    plt.close(fig)
    print("Saved results/Cobaya_Ratio_Falsifiability.png")
