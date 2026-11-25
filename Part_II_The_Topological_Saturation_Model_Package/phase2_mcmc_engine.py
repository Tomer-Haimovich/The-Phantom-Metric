import os
import sys
import time
import argparse
import multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cmb_likelihood_utils import (
    TOPO_RATIO, L_MIN, L_MAX_TT, L_MAX_POL,
    wishart_chi2_exact, binned_knox_chi2, r2_1to1, gelman_rubin_rhat
)

PKG_DIR = os.path.dirname(os.path.abspath(__file__))
REF_FILE = os.path.join(PKG_DIR, "mcmc_planck_target.npz")
CHK_FILE = os.path.join(PKG_DIR, "mcmc_chain_checkpoint.npz")
TMP_CHK_FILE = os.path.join(PKG_DIR, "mcmc_chain_checkpoint.tmp.npz")
REPORT_FILE = os.path.join(PKG_DIR, "mcmc_audit_report.txt")
CORNER_PNG = os.path.join(PKG_DIR, "phantom_mcmc_corner.png")

PARAM_NAMES = [r"$\omega_b$", r"$h$", r"$\ln(10^{10}A_s)$", r"$n_s$", r"$\tau_{\mathrm{reio}}$", r"$\sigma_8$ (Derived)"]
PLANCK_TRUTH = np.array([0.02237, 0.6736, np.log(21.00), 0.9649, 0.0544])

# Wide flat physical priors (zero informative priors on omega_b, h, A_s, n_s in cmbonly mode!)
BOUNDS_MIN = np.array([0.01900, 0.6000, 2.8500, 0.9000, 0.0100])
BOUNDS_MAX = np.array([0.02600, 0.7600, 3.2500, 1.0400, 0.1200])

# True independent external priors:
# 1) Low-l Planck SimAll optical depth (l < 30, independent of high-l >= 50)
TAU_PRIOR_MEAN, TAU_PRIOR_SIGMA = 0.0544, 0.0073
# 2) Independent primordial BBN Deuterium prior (Cooke+2018 / Fields+2020; ONLY used in --mode cmbbbn)
BBN_WB_MEAN, BBN_WB_SIGMA = 0.02233, 0.00036

# Worker-local globals (populated safely via spawn initializer to avoid OpenMP fork deadlocks)
_W_REF = None
_W_MODE = "cmbonly"
_W_CLASS = None


def _init_worker(ref_dict, mode, omp_threads):
    global _W_REF, _W_MODE, _W_CLASS
    os.environ["OMP_NUM_THREADS"] = str(omp_threads)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    from classy import Class
    _W_REF = ref_dict
    _W_MODE = mode
    _W_CLASS = Class


def _eval_one(theta, ref, mode, cls_type):
    omega_b, h, lnAs, n_s, tau_reio = theta
    if np.any(theta < BOUNDS_MIN) or np.any(theta > BOUNDS_MAX):
        return -np.inf, np.full(6, np.nan)

    A_s = 1.0e-10 * np.exp(lnAs)
    cosmo = cls_type()
    try:
        cosmo.set({
            'output': 'tCl,pCl,lCl,mPk',
            'lensing': 'yes',
            'l_max_scalars': L_MAX_TT,
            'P_k_max_1/Mpc': 5.0,
            'h': float(h),
            'omega_b': float(omega_b),
            'omega_cdm': 1.0e-12,
            'A_s': float(A_s),
            'n_s': float(n_s),
            'tau_reio': float(tau_reio)
        })
        cosmo.compute()
        cls = cosmo.lensed_cl(L_MAX_TT)
        fac = ref['fac']
        tt = cls['tt'][L_MIN:L_MAX_TT + 1] * fac
        te = cls['te'][L_MIN:L_MAX_TT + 1] * fac
        ee = cls['ee'][L_MIN:L_MAX_TT + 1] * fac
        s8 = float(cosmo.sigma8())
    except Exception:
        return -np.inf, np.full(6, np.nan)
    finally:
        try:
            cosmo.struct_cleanup()
            cosmo.empty()
        except Exception:
            pass

    ell = ref['ell']
    chi2_wishart = wishart_chi2_exact(
        tt, te, ee, ref['tt_ref'], ref['te_ref'], ref['ee_ref'],
        ell, ref['n_tt_dl'], ref['n_ee_dl']
    )
    if not np.isfinite(chi2_wishart):
        return -np.inf, np.full(6, np.nan)

    knox_b = binned_knox_chi2(
        tt, te, ee, ref['tt_ref'], ref['te_ref'], ref['ee_ref'],
        ell, ref['n_tt_dl'], ref['n_ee_dl'], dl=30
    )

    chi2_lowE = ((tau_reio - TAU_PRIOR_MEAN) / TAU_PRIOR_SIGMA)**2
    if mode == "cmbbbn":
        chi2_bbn = ((omega_b - BBN_WB_MEAN) / BBN_WB_SIGMA)**2
        chi2_total = chi2_wishart + chi2_lowE + chi2_bbn
    else:
        chi2_total = chi2_wishart + chi2_lowE

    log_prob = -0.5 * chi2_total
    pol_m = ell <= L_MAX_POL
    r2_tt = r2_1to1(tt, ref['tt_ref'])
    r2_te = r2_1to1(te[pol_m], ref['te_ref'][pol_m])
    r2_ee = r2_1to1(ee[pol_m], ref['ee_ref'][pol_m])

    return float(log_prob), np.array([s8, r2_tt, r2_te, r2_ee, chi2_wishart, knox_b['rms_sigma_per_bin']], dtype=float)


def _worker_eval(theta):
    return _eval_one(theta, _W_REF, _W_MODE, _W_CLASS)


def load_reference():
    if not os.path.exists(REF_FILE):
        print("ERROR: Please run phase1_calibrate_and_bbn.py first!")
        sys.exit(1)
    d = np.load(REF_FILE)
    return {k: d[k] for k in d.files}


def smooth_1d(arr, sigma=1.1):
    radius = int(3.0 * sigma + 0.5)
    x = np.arange(-radius, radius + 1, dtype=float)
    k = np.exp(-0.5 * (x / sigma)**2)
    k /= np.sum(k)
    return np.convolve(np.pad(arr, radius, mode='edge'), k, mode='valid')


def smooth_2d(mat, sigma=1.0):
    out = np.apply_along_axis(smooth_1d, 0, mat, sigma)
    return np.apply_along_axis(smooth_1d, 1, out, sigma)


def generate_corner_and_report(chain, log_probs, metrics, step_done, acc_rate, ref, mode="cmbonly"):
    burn = max(1, int(0.30 * step_done)) if step_done >= 10 else 0
    chain_pb = chain[burn:]
    rhat_5 = gelman_rubin_rhat(chain_pb)

    flat_theta = chain_pb.reshape(-1, 5)
    flat_logp = log_probs[burn:].reshape(-1)
    flat_met = metrics[burn:].reshape(-1, 6)

    valid = np.isfinite(flat_logp) & np.isfinite(flat_met[:, 0])
    flat_theta = flat_theta[valid]
    flat_logp = flat_logp[valid]
    flat_met = flat_met[valid]
    if len(flat_logp) < 10:
        return

    s8_ref = float(ref['sigma8_ref'])
    data_6d = np.column_stack([flat_theta, flat_met[:, 0]])
    truth_6d = np.append(PLANCK_TRUTH, s8_ref)
    best_idx = int(np.argmax(flat_logp))
    best_6d = data_6d[best_idx]
    best_m = flat_met[best_idx]

    wb_m, wb_s = float(np.mean(flat_theta[:, 0])), float(np.std(flat_theta[:, 0]))
    wph_m, wph_s = TOPO_RATIO * wb_m, TOPO_RATIO * wb_s

    ell = ref['ell']
    pol_m = ell <= L_MAX_POL
    r2_tt0 = r2_1to1(ref['tt_ph0'], ref['tt_ref'])
    r2_te0 = r2_1to1(ref['te_ph0'][pol_m], ref['te_ref'][pol_m])
    r2_ee0 = r2_1to1(ref['ee_ph0'][pol_m], ref['ee_ref'][pol_m])
    chi2_w0 = wishart_chi2_exact(
        ref['tt_ph0'], ref['te_ph0'], ref['ee_ph0'],
        ref['tt_ref'], ref['te_ref'], ref['ee_ref'],
        ell, ref['n_tt_dl'], ref['n_ee_dl']
    )

    mode_str = "CMB-Only Undiluted Wishart (k=5)" if mode == "cmbonly" else "CMB + BBN D/H Undiluted Wishart (k=5)"
    lines = [
        "=" * 112,
        f" LAMBDATSM CLEAN COVARIANT v2 — HONEST 5-PARAMETER MCMC AUDIT (0% CDM | Steps: {step_done} | Samples: {len(flat_logp)})",
        f" Likelihood: {mode_str} | Acceptance Rate: {acc_rate*100:.1f}%",
        "=" * 112,
        f"{'Parameter / Metric':<28} | {'Planck 2018 (k=6)':<22} | {'LambdaTSM Posterior (k=5)':<28} | {'Best-Fit':<10} | {'R-hat'}",
        "-" * 112,
        f"{'omega_b (Baryon Density)':<28} | {'0.02237 +/- 0.00015':<22} | {wb_m:.5f} +/- {wb_s:.5f}{'':<10} | {best_6d[0]:<10.5f} | {rhat_5[0]:.3f}",
        f"{'omega_cdm / omega_phantom':<28} | {'0.12000 +/- 0.00120':<22} | {wph_m:.5f} +/- {wph_s:.5f} (0-P){'':<4} | {TOPO_RATIO*best_6d[0]:<10.5f} | {rhat_5[0]:.3f}",
        f"{'h (Reduced Hubble)':<28} | {'0.6736  +/- 0.0054':<22} | {np.mean(flat_theta[:,1]):.4f}  +/- {np.std(flat_theta[:,1]):.4f}{'':<11} | {best_6d[1]:<10.4f} | {rhat_5[1]:.3f}",
        f"{'ln(10^10 A_s)':<28} | {'3.0445  +/- 0.0140':<22} | {np.mean(flat_theta[:,2]):.4f}  +/- {np.std(flat_theta[:,2]):.4f}{'':<11} | {best_6d[2]:<10.4f} | {rhat_5[2]:.3f}",
        f"{'n_s (Spectral Index)':<28} | {'0.9649  +/- 0.0042':<22} | {np.mean(flat_theta[:,3]):.4f}  +/- {np.std(flat_theta[:,3]):.4f}{'':<11} | {best_6d[3]:<10.4f} | {rhat_5[3]:.3f}",
        f"{'tau_reio (Optical Depth)':<28} | {'0.0544  +/- 0.0073':<22} | {np.mean(flat_theta[:,4]):.4f}  +/- {np.std(flat_theta[:,4]):.4f}{'':<11} | {best_6d[4]:<10.4f} | {rhat_5[4]:.3f}",
        f"{'sigma_8 (Derived z=0)':<28} | {f'{s8_ref:.4f}  +/- 0.0060':<22} | {np.mean(flat_met[:,0]):.4f}  +/- {np.std(flat_met[:,0]):.4f}{'':<11} | {best_6d[5]:<10.4f} |   --",
        "-" * 112,
        f"Fixed Baseline (Planck params) : R2_TT={r2_tt0:.6f} | R2_TE={r2_te0:.6f} | R2_EE={r2_ee0:.6f} | TRUE Wishart Delta chi^2 = {chi2_w0:.4f}",
        f"Chain Best-Fit (5-param)       : R2_TT={best_m[1]:.6f} | R2_TE={best_m[2]:.6f} | R2_EE={best_m[3]:.6f} | TRUE Wishart Delta chi^2 = {best_m[4]:.4f}",
        f"TRUE Properly-Binned Knox RMS  : {best_m[5]:.4f} sigma/bin (dl=30, NO /30 dilution, NO hidden 1.5% noise floor)",
        "=" * 112,
    ]
    report_text = "\n".join(lines)
    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write(report_text + "\n")
    print("\n" + report_text)

    ndim_plot = 6
    fig, axes = plt.subplots(ndim_plot, ndim_plot, figsize=(13.5, 13.5), dpi=240)
    for i in range(ndim_plot):
        for j in range(ndim_plot):
            ax = axes[i, j]
            if j > i:
                ax.axis('off')
                continue
            x = data_6d[:, j]
            q1, q16, q50, q84, q99 = np.percentile(x, [0.5, 16, 50, 84, 99.5])
            span_x = max(q99 - q1, 1e-5)
            xlim_lo = min(q1, truth_6d[j]) - 0.22 * span_x
            xlim_hi = max(q99, truth_6d[j]) + 0.22 * span_x

            if i == j:
                counts, bins, _ = ax.hist(x, bins=26, range=(xlim_lo, xlim_hi), density=True, color='#0f766e', alpha=0.40, edgecolor='none')
                centers = 0.5 * (bins[:-1] + bins[1:])
                ax.plot(centers, smooth_1d(counts, sigma=1.1), color='#134e4a', lw=1.8)
                ax.axvline(truth_6d[j], color='#be123c', ls='--', lw=1.5)
                ax.axvline(best_6d[j], color='#1d4ed8', ls='-', lw=1.5)
                ax.set_xlim(xlim_lo, xlim_hi)
                fmt_t = f"{q50:.5f}_{{-{q50-q16:.5f}}}^{{+{q84-q50:.5f}}}" if j == 0 else f"{q50:.4f}_{{-{q50-q16:.4f}}}^{{+{q84-q50:.4f}}}"
                ax.set_title(f"{PARAM_NAMES[j]}\n${fmt_t}$", fontsize=9.2)
                ax.set_yticks([])
            else:
                y = data_6d[:, i]
                yq1, yq99 = np.percentile(y, [0.5, 99.5])
                span_y = max(yq99 - yq1, 1e-5)
                ylim_lo = min(yq1, truth_6d[i]) - 0.22 * span_y
                ylim_hi = max(yq99, truth_6d[i]) + 0.22 * span_y
                H, xedges, yedges = np.histogram2d(x, y, bins=26, range=[[xlim_lo, xlim_hi], [ylim_lo, ylim_hi]], density=True)
                H_s = smooth_2d(H.T, sigma=1.0)
                hmax = float(H_s.max())
                if hmax > 0:
                    H_flat = np.sort(H_s.ravel())[::-1]
                    cdf = np.cumsum(H_flat)
                    cdf /= cdf[-1] if cdf[-1] > 0 else 1.0
                    lvl68 = min(float(H_flat[np.searchsorted(cdf, 0.68)]), hmax * 0.90)
                    lvl95 = min(float(H_flat[np.searchsorted(cdf, 0.95)]), lvl68 * 0.50)
                    X, Y = 0.5 * (xedges[:-1] + xedges[1:]), 0.5 * (yedges[:-1] + yedges[1:])
                    if 0.0 <= lvl95 < lvl68 < hmax:
                        ax.contourf(X, Y, H_s, levels=[lvl95, lvl68, hmax * 1.01], colors=['#99f6e4', '#0f766e'], alpha=0.80)
                        ax.contour(X, Y, H_s, levels=[lvl95, lvl68], colors=['#14b8a6', '#134e4a'], linewidths=1.0)
                ax.axvline(truth_6d[j], color='#be123c', ls='--', lw=1.1, alpha=0.75)
                ax.axhline(truth_6d[i], color='#be123c', ls='--', lw=1.1, alpha=0.75)
                ax.plot(best_6d[j], best_6d[i], 'o', color='#1d4ed8', ms=4.5)
                ax.set_xlim(xlim_lo, xlim_hi); ax.set_ylim(ylim_lo, ylim_hi)

            if i == ndim_plot - 1:
                ax.set_xlabel(PARAM_NAMES[j], fontsize=9.5)
            else:
                ax.set_xticklabels([])
            if j == 0 and i > 0:
                ax.set_ylabel(PARAM_NAMES[i], fontsize=9.5)
            elif j > 0:
                ax.set_yticklabels([])

    fig.suptitle(f"$\\Lambda$TSM Clean v2 — 5-Parameter Undiluted Wishart MCMC ({mode_str}, $N={len(flat_logp)}$)", fontsize=13.5, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.985])
    plt.savefig(CORNER_PNG, dpi=240)
    plt.close(fig)
    print(f"--> Saved Corner Plot: {CORNER_PNG}")


def run_mcmc(total_steps=300, nwalkers=16, force_fresh=False, mode="cmbonly", n_procs=None, omp_threads=1):
    ref = load_reference()
    ndim = 5
    rng = np.random.default_rng(42)

    # Guardrail: warn immediately if target was generated with the old buggy perturbations.c!
    chi2_chk = wishart_chi2_exact(
        ref['tt_ph0'], ref['te_ph0'], ref['ee_ph0'],
        ref['tt_ref'], ref['te_ref'], ref['ee_ref'],
        ref['ell'], ref['n_tt_dl'], ref['n_ee_dl']
    )
    if chi2_chk > 100.0:
        print("=" * 90)
        print(f" WARNING: mcmc_planck_target.npz has Baseline Wishart Delta chi^2 = {chi2_chk:.1f}!")
        print(" This means CLASS was NOT yet recompiled with the clean v2 perturbations.c!")
        print(" Recompile CLASS (copy C files -> rm -rf build -> make clean && make -j4 -> pip install)")
        print(" and re-run phase1_calibrate_and_bbn.py BEFORE running MCMC.")
        print("=" * 90)

    if n_procs is None:
        cpu_n = os.cpu_count() or 4
        n_procs = max(1, min(nwalkers // 2, cpu_n))

    init_center = np.array([0.02255, 0.6780, 3.0550, 0.9720, 0.0580])
    init_spread = np.array([0.00035, 0.0090, 0.0200, 0.0120, 0.0100])

    total_props = 0
    total_acc = 0

    ctx = mp.get_context("spawn")
    pool = ctx.Pool(processes=n_procs, initializer=_init_worker, initargs=(ref, mode, omp_threads)) if n_procs > 1 else None
    serial_cls = None
    if pool is None:
        from classy import Class
        serial_cls = Class

    def eval_batch(theta_batch):
        if pool is not None:
            return pool.map(_worker_eval, list(theta_batch))
        return [_eval_one(th, ref, mode, serial_cls) for th in theta_batch]

    try:
        if os.path.exists(CHK_FILE) and not force_fresh:
            chk = np.load(CHK_FILE)
            chain = list(chk['chain'])
            log_probs = list(chk['log_probs'])
            metrics = list(chk['metrics'])
            total_props = int(chk['total_props']) if 'total_props' in chk else 1
            total_acc = int(chk['total_acc']) if 'total_acc' in chk else 0
            start_step = len(chain)
            pos = chain[-1].copy()
            cur_logp = log_probs[-1].copy()
            cur_met = metrics[-1].copy()
            nwalkers = pos.shape[0]
            print(f"[CHECKPOINT] Resuming from step {start_step}/{total_steps} using {n_procs} worker procs (OMP={omp_threads}).")
            if start_step >= total_steps:
                generate_corner_and_report(np.array(chain), np.array(log_probs), np.array(metrics), start_step, total_acc / max(1, total_props), ref, mode=mode)
                return
        else:
            print(f"[FRESH PARALLEL MCMC] {nwalkers} walkers x {total_steps} steps | Mode: {mode} | Procs: {n_procs} (OMP={omp_threads})")
            pos = init_center + rng.normal(0.0, init_spread, size=(nwalkers, ndim))
            pos = np.clip(pos, BOUNDS_MIN + 1e-5, BOUNDS_MAX - 1e-5)
            results = eval_batch(pos)
            cur_logp = np.array([r[0] for r in results])
            cur_met = np.array([r[1] for r in results])
            chain, log_probs, metrics = [pos.copy()], [cur_logp.copy()], [cur_met.copy()]
            start_step = 1
            np.savez(TMP_CHK_FILE, chain=np.array(chain), log_probs=np.array(log_probs), metrics=np.array(metrics), total_props=0, total_acc=0)
            os.replace(TMP_CHK_FILE, CHK_FILE)

        a_scale = 2.0
        half = nwalkers // 2
        t_start = time.time()

        for step in range(start_step + 1, total_steps + 1):
            for first, second in [(slice(0, half), slice(half, nwalkers)), (slice(half, nwalkers), slice(0, half))]:
                active = pos[first]
                comp = pos[second]
                n_act = len(active)
                u = rng.uniform(0.0, 1.0, size=n_act)
                z = ((a_scale - 1.0) * u + 1.0)**2 / a_scale
                partner_idx = rng.integers(0, len(comp), size=n_act)
                proposals = comp[partner_idx] + z[:, None] * (active - comp[partner_idx])
                prop_res = eval_batch(proposals)
                prop_logp = np.array([r[0] for r in prop_res])
                prop_met = np.array([r[1] for r in prop_res])

                log_accept = (ndim - 1.0) * np.log(z) + (prop_logp - cur_logp[first])
                accept = np.isfinite(prop_logp) & (np.log(rng.uniform(0.0, 1.0, size=n_act)) < log_accept)
                active[accept] = proposals[accept]
                cur_logp[first][accept] = prop_logp[accept]
                cur_met[first][accept] = prop_met[accept]
                pos[first] = active
                total_props += n_act
                total_acc += int(np.sum(accept))

            chain.append(pos.copy())
            log_probs.append(cur_logp.copy())
            metrics.append(cur_met.copy())

            if step % 5 == 0 or step == total_steps:
                np.savez(TMP_CHK_FILE, chain=np.array(chain), log_probs=np.array(log_probs), metrics=np.array(metrics),
                         total_props=total_props, total_acc=total_acc)
                os.replace(TMP_CHK_FILE, CHK_FILE)
                elapsed = time.time() - t_start
                steps_done_now = step - start_step
                eta_s = (elapsed / max(1, steps_done_now)) * (total_steps - step)
                print(f"[Step {step:4d}/{total_steps}] Acc: {100*total_acc/max(1,total_props):4.1f}% | Best logP: {np.max(cur_logp):8.3f} | Elapsed: {elapsed:6.1f}s | ETA: {eta_s/60.0:4.1f} min")

            if step % 25 == 0 or step == total_steps:
                generate_corner_and_report(np.array(chain), np.array(log_probs), np.array(metrics), step, total_acc / max(1, total_props), ref, mode=mode)
    finally:
        if pool is not None:
            pool.close()
            pool.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LambdaTSM Clean v2 — Parallel Undiluted 5-Parameter Wishart MCMC")
    parser.add_argument("--fresh", action="store_true", help="Ignore existing checkpoint and run fresh from scratch")
    parser.add_argument("--test", action="store_true", help="Fast sanity check (20 steps x 12 walkers)")
    parser.add_argument("--steps", type=int, default=300, help="Total MCMC steps per walker (default: 300)")
    parser.add_argument("--walkers", type=int, default=16, help="Number of affine-invariant walkers (default: 16)")
    parser.add_argument("--procs", type=int, default=None, help="Parallel worker processes (default: auto up to walkers/2)")
    parser.add_argument("--omp", type=int, default=1, help="OMP_NUM_THREADS per worker process (default: 1)")
    parser.add_argument("--mode", type=str, default="cmbonly", choices=["cmbonly", "cmbbbn"],
                        help="cmbonly = Pure CMB + SimAll tau (default); cmbbbn = CMB + SimAll tau + Cooke+2018 BBN D/H")
    args = parser.parse_args()
    steps = 20 if args.test else args.steps
    walkers = 12 if args.test else args.walkers
    run_mcmc(total_steps=steps, nwalkers=walkers, force_fresh=args.fresh or args.test,
             mode=args.mode, n_procs=args.procs, omp_threads=args.omp)
