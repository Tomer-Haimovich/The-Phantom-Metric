import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from getdist import loadMCSamples, plots

os.makedirs("results", exist_ok=True)
R_LOCK = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * math.pi**2))

s_cdm_pl  = loadMCSamples("chains/lcdm_planck_only",    settings={"ignore_rows": 0.3})
s_tsm_all = loadMCSamples("chains/ltsm_planck_desi_sn", settings={"ignore_rows": 0.3})
s_cdm_all = loadMCSamples("chains/lcdm_planck_desi_sn", settings={"ignore_rows": 0.3})

latex_map = {
    "omega_b": r"\Omega_b h^2",
    "H0":      r"H_0\ [\mathrm{km\,s^{-1}\,Mpc^{-1}}]",
    "Omega_m": r"\Omega_{m,0}",
    "n_s":     r"n_s",
    "sigma8":  r"\sigma_8",
    "S8":      r"S_8",
}
for s in (s_tsm_all, s_cdm_all):
    for name, lbl in latex_map.items():
        p = s.paramNames.parWithName(name)
        if p is not None:
            p.label = lbl

# 1. Clean Triangle Plot (no collinear a0 line, clean ticks)
g = plots.get_subplot_plotter(subplot_size=2.2)
g.settings.axes_fontsize = 10
g.settings.lab_fontsize = 11.5
g.settings.legend_fontsize = 11.5
g.settings.alpha_filled_add = 0.78
g.settings.linewidth_contour = 1.4
g.settings.linewidth = 1.8
g.settings.axis_tick_max_labels = 4
g.settings.tight_layout = True

g.triangle_plot(
    [s_cdm_all, s_tsm_all],
    params=["omega_b", "H0", "Omega_m", "n_s", "sigma8", "S8"],
    filled=True,
    contour_colors=["#475569", "#0f766e"],
    legend_labels=[
        r"$\Lambda\mathrm{CDM}$ ($k=6$, Planck + DESI + Pantheon+)",
        r"$\Lambda\mathrm{TSM}$ ($k=5$, $\omega_{\mathrm{ph}} = \mathcal{R}_{4/3}^{\circ}\,\omega_b$)",
    ],
    legend_loc="upper right",
)
out_corner = os.path.abspath("results/Cobaya_Planck_DESI_SNe_Corner.png")
g.export(out_corner, dpi=300)
plt.close("all")
print("Saved clean corner plot:", out_corner)

# 2. Clean Partition-Lock Falsifiability Plot (legend in empty upper-right space)
def extract_1d_kde(samples_obj, param_name):
    dens = samples_obj.get1DDensity(param_name)
    xs = np.linspace(4.98, 5.95, 650)
    ys = np.array([dens(x) for x in xs], dtype=float)
    ys /= np.max(ys)
    return xs, ys, float(samples_obj.mean(param_name)), float(samples_obj.std(param_name))

x_pl, y_pl, mu_pl, sd_pl = extract_1d_kde(s_cdm_pl,  "dm_baryon_ratio")
x_al, y_al, mu_al, sd_al = extract_1d_kde(s_cdm_all, "dm_baryon_ratio")
z_pl = (R_LOCK - mu_pl) / sd_pl
z_al = (R_LOCK - mu_al) / sd_al

plt.style.use("default")
fig, ax = plt.subplots(figsize=(8.8, 5.3), dpi=300)

ax.plot(x_pl, y_pl, color="#0284c7", lw=2.3,
        label=f"$\\Lambda$CDM: Planck 2018 Only\n({mu_pl:.4f} $\\pm$ {sd_pl:.4f}, diff {z_pl:+.2f}$\\sigma$)")
ax.fill_between(x_pl, 0.0, y_pl, color="#0284c7", alpha=0.14)

ax.plot(x_al, y_al, color="#7c3aed", lw=2.3,
        label=f"$\\Lambda$CDM: Planck + DESI BAO + Pantheon+\n({mu_al:.4f} $\\pm$ {sd_al:.4f}, diff {z_al:+.2f}$\\sigma$)")
ax.fill_between(x_al, 0.0, y_al, color="#7c3aed", alpha=0.16)

ax.axvline(R_LOCK, color="#be123c", ls="--", lw=2.4, zorder=5,
           label=f"$\\Lambda$TSM Rigid Lock: $\\mathcal{{R}}_{{4/3}}^{{\\circ}} = {R_LOCK:.4f}$")

ax.set_xlim(5.00, 5.95)
ax.set_ylim(0.0, 1.15)
ax.set_xlabel(r"Dark-to-Baryon Density Ratio $\omega_{\mathrm{cdm}} / \omega_b$", fontsize=12)
ax.set_ylabel("Normalized Posterior Probability", fontsize=12)
ax.set_title(r"Direct Observational Test of the $\Lambda\mathrm{TSM}$ Partition Lock", fontsize=13, fontweight="bold", pad=10)
ax.grid(True, ls=":", alpha=0.45)
ax.legend(loc="upper right", fontsize=9.4, frameon=True, framealpha=0.96, edgecolor="#cbd5e1")

plt.tight_layout()
out_ratio = os.path.abspath("results/Cobaya_Ratio_Falsifiability.png")
plt.savefig(out_ratio, dpi=300, bbox_inches="tight", pad_inches=0.14)
plt.close(fig)
print("Saved clean ratio plot :", out_ratio)
