"""
Exact CMB Likelihood & Convergence Diagnostics (LambdaTSM Clean Rebuild v2)
- Zero hidden variance floors (no artificial 1.5% or foreground padding)
- Zero illegal /30 or /300 chi^2 dilution
- Exact correlated 2x2 TT+TE+EE Wishart log-likelihood (Hamimeche & Lewis 2008 / Knox 1995)
- Properly binned (Delta ell = 30) diagonal Knox chi^2 diagnostic
- Gelman-Rubin R-hat convergence diagnostic across split chains
"""
import numpy as np

# Exact Planck 2018 high-l noise specifications (Knox 1995)
F_SKY_PLANCK = 0.65
THETA_FWHM_ARCMIN = 5.5
NOISE_T_UK_ARCMIN = 27.0
NOISE_P_UK_ARCMIN = 55.0
L_MIN = 50
L_MAX_TT = 2500
L_MAX_POL = 2000

# Exact 4D/3D phenomenological saturation ratio ansatz
TOPO_RATIO = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * np.pi**2))


def compute_planck_noise_dl(ell_arr):
    """Exact Knox (1995) beam-deconvolved white detector noise in D_l [uK^2]."""
    ell = np.asarray(ell_arr, dtype=float)
    fac_dl = ell * (ell + 1.0) / (2.0 * np.pi)
    theta_rad = (THETA_FWHM_ARCMIN / 60.0) * (np.pi / 180.0)
    beam_window = np.exp(ell * (ell + 1.0) * (theta_rad**2) / (8.0 * np.log(2.0)))
    n_tt_cl = ((NOISE_T_UK_ARCMIN / 60.0) * (np.pi / 180.0))**2 * beam_window
    n_ee_cl = ((NOISE_P_UK_ARCMIN / 60.0) * (np.pi / 180.0))**2 * beam_window
    return n_tt_cl * fac_dl, n_ee_cl * fac_dl


def wishart_chi2_exact(tt_mod, te_mod, ee_mod, tt_obs, te_obs, ee_obs,
                       ell_arr, n_tt_dl, n_ee_dl, f_sky=F_SKY_PLANCK, lmax_pol=L_MAX_POL):
    """
    Exact undiluted correlated 2x2 Wishart CMB log-likelihood (-2 ln L / L_max):
      - For ell <= lmax_pol: full 2x2 covariance [[TT+N_TT, TE], [TE, EE+N_EE]]
      - For ell >  lmax_pol: 1x1 TT-only Wishart
    Identically 0.0 when model == target; strictly positive otherwise.
    """
    ell = np.asarray(ell_arr, dtype=float)
    pol_mask = ell <= lmax_pol
    tt_only_mask = ~pol_mask

    # 2x2 correlated TT+TE+EE block
    a_m = tt_mod[pol_mask] + n_tt_dl[pol_mask]
    d_m = ee_mod[pol_mask] + n_ee_dl[pol_mask]
    b_m = te_mod[pol_mask]
    det_m = a_m * d_m - b_m * b_m

    a_o = tt_obs[pol_mask] + n_tt_dl[pol_mask]
    d_o = ee_obs[pol_mask] + n_ee_dl[pol_mask]
    b_o = te_obs[pol_mask]
    det_o = a_o * d_o - b_o * b_o

    if np.any(det_m <= 0.0) or np.any(det_o <= 0.0):
        return np.inf

    tr_term = (d_m * a_o + a_m * d_o - 2.0 * b_m * b_o) / det_m
    chi2_2x2 = np.sum((2.0 * ell[pol_mask] + 1.0) * f_sky * (tr_term + np.log(det_m / det_o) - 2.0))

    # 1x1 TT-only block (ell > lmax_pol)
    c_m = tt_mod[tt_only_mask] + n_tt_dl[tt_only_mask]
    c_o = tt_obs[tt_only_mask] + n_tt_dl[tt_only_mask]
    if np.any(c_m <= 0.0) or np.any(c_o <= 0.0):
        return np.inf
    chi2_1x1 = np.sum((2.0 * ell[tt_only_mask] + 1.0) * f_sky * (c_o / c_m + np.log(c_m / c_o) - 1.0))

    return float(max(0.0, chi2_2x2 + chi2_1x1))


def binned_knox_chi2(tt_mod, te_mod, ee_mod, tt_obs, te_obs, ee_obs,
                     ell_arr, n_tt_dl, n_ee_dl, f_sky=F_SKY_PLANCK, dl=30, lmax_pol=L_MAX_POL):
    """
    Properly binned (Delta ell = dl) diagonal Knox chi^2.
    For a bin of width dl, the variance of the bin mean is (1/dl^2) * sum(var_ell),
    so chi^2 is NOT diluted by /dl!
    """
    ell = np.asarray(ell_arr, dtype=float)
    var_tt = (2.0 / ((2.0 * ell + 1.0) * f_sky)) * (tt_obs + n_tt_dl)**2
    var_ee = (2.0 / ((2.0 * ell + 1.0) * f_sky)) * (ee_obs + n_ee_dl)**2
    var_te = (1.0 / ((2.0 * ell + 1.0) * f_sky)) * (te_obs**2 + (tt_obs + n_tt_dl) * (ee_obs + n_ee_dl))

    def _bin_chi2(m, o, v, mask):
        diff = m[mask] - o[mask]
        vv = v[mask]
        nb = len(diff) // dl
        tot = 0.0
        for i in range(nb):
            s = slice(i * dl, (i + 1) * dl)
            mean_diff = np.mean(diff[s])
            var_mean = np.sum(vv[s]) / float(dl * dl)
            tot += (mean_diff * mean_diff) / var_mean
        return float(tot), int(nb)

    pol_mask = ell <= lmax_pol
    all_mask = np.ones_like(ell, dtype=bool)
    c_tt, n_tt = _bin_chi2(tt_mod, tt_obs, var_tt, all_mask)
    c_te, n_te = _bin_chi2(te_mod, te_obs, var_te, pol_mask)
    c_ee, n_ee = _bin_chi2(ee_mod, ee_obs, var_ee, pol_mask)
    n_tot = n_tt + n_te + n_ee
    c_tot = c_tt + c_te + c_ee
    rms_per_bin = float(np.sqrt(c_tot / n_tot)) if n_tot > 0 else np.nan
    return {
        "chi2_tt": c_tt, "chi2_te": c_te, "chi2_ee": c_ee,
        "chi2_total": c_tot, "n_bins": n_tot, "rms_sigma_per_bin": rms_per_bin
    }


def r2_1to1(y_pred, y_true):
    y_p = np.asarray(y_pred, dtype=float)
    y_t = np.asarray(y_true, dtype=float)
    return float(1.0 - np.sum((y_p - y_t)**2) / np.sum((y_t - np.mean(y_t))**2))


def gelman_rubin_rhat(chain_post_burn):
    """
    Compute split-walker Gelman-Rubin R-hat for an ensemble chain of shape
    (n_steps, n_walkers, n_dim). Splits walkers into 2 independent sub-ensembles
    and splits each in half along time -> 4 sub-chains.
    R-hat < 1.05 indicates solid convergence.
    """
    n_steps, n_walkers, n_dim = chain_post_burn.shape
    if n_steps < 10 or n_walkers < 4:
        return np.full(n_dim, np.nan)
    half_t = n_steps // 2
    half_w = n_walkers // 2
    subchains = [
        chain_post_burn[:half_t, :half_w, :].reshape(-1, n_dim),
        chain_post_burn[half_t:2 * half_t, :half_w, :].reshape(-1, n_dim),
        chain_post_burn[:half_t, half_w:2 * half_w, :].reshape(-1, n_dim),
        chain_post_burn[half_t:2 * half_t, half_w:2 * half_w, :].reshape(-1, n_dim),
    ]
    m = len(subchains)
    n = subchains[0].shape[0]
    means = np.array([sc.mean(axis=0) for sc in subchains])
    vars_ = np.array([sc.var(axis=0, ddof=1) for sc in subchains])
    W = np.mean(vars_, axis=0)
    B = n * np.var(means, axis=0, ddof=1)
    var_hat = ((n - 1.0) / n) * W + (1.0 / n) * B
    with np.errstate(divide='ignore', invalid='ignore'):
        rhat = np.sqrt(var_hat / W)
    return rhat
