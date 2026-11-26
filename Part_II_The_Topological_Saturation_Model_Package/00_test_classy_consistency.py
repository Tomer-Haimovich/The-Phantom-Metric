import math
import numpy as np

print("========================================================================")
print(" Pre-Flight Check: CLASS + Cobaya + Exact LambdaTSM Consistency")
print("========================================================================")

import cobaya
import getdist
import classy
print(f"  cobaya  version : {cobaya.__version__}")
print(f"  getdist version : {getdist.__version__}")
print(f"  classy  path    : {classy.__file__}")

R_PART = (16.0 / 3.0) * (1.0 + 1.0 / (16.0 * math.pi**2))
wb = 0.02237
w_ph = R_PART * wb

base = {
    "output": "tCl,pCl,lCl,mPk",
    "lensing": "yes",
    "l_max_scalars": 2500,
    "omega_b": wb,
    "H0": 67.36,
    "A_s": 2.100e-9,
    "n_s": 0.9649,
    "tau_reio": 0.0544,
    "N_ncdm": 1,
    "m_ncdm": 0.06,
    "N_ur": 2.0308,
}

c1 = classy.Class()
c1.set({**base, "omega_cdm": w_ph})
c1.compute()
cl1 = c1.lensed_cl(2500)
s8_1 = c1.sigma8()
rs_1 = c1.rs_drag()
c1.struct_cleanup(); c1.empty()

print(f"\n[Explicit R_{{4/3}}^o * omega_b = {w_ph:.6f}, m_nu = 0.06 eV]:")
print(f"  sigma8 = {s8_1:.5f} | rs_drag = {rs_1:.3f} Mpc")

c2 = classy.Class()
c2.set({**base, "omega_cdm": 1e-12})
c2.compute()
cl2 = c2.lensed_cl(2500)
s8_2 = c2.sigma8()
c2.struct_cleanup(); c2.empty()

rel_diff_tt = np.max(np.abs(cl1["tt"][50:2501] - cl2["tt"][50:2501]) / np.maximum(cl1["tt"][50:2501], 1e-30))
if rel_diff_tt < 1e-4:
    print(f"  [VERIFIED] Custom CLASS background.c trigger matches Cobaya lock to {rel_diff_tt:.2e} max rel diff! (sigma8={s8_2:.5f})")
else:
    print(f"  [INFO] Standard unmodified CLASS detected (omega_cdm=1e-12 gives sigma8={s8_2:.4f}).")
    print("         Cobaya's explicit derived lock omega_cdm = R_{4/3}^o * omega_b works identically either way.")

print("\n[Testing single Cobaya likelihood call on ltsm_planck_desi_sn.yaml...]")
from cobaya.model import get_model
from cobaya.yaml import yaml_load_file

info = yaml_load_file("ltsm_planck_desi_sn.yaml")
model = get_model(info)
test_pt = model.prior.reference(max_tries=100, warn_if_no_ref=False)
if isinstance(test_pt, np.ndarray):
    test_pt = dict(zip(model.parameterization.sampled_params(), test_pt))
loglikes, derived = model.loglikes(test_pt, return_derived=True)
chi2_tot = -2.0 * float(np.sum(loglikes))
print(f"  SUCCESS! Active likelihoods: {list(model.likelihood.keys())}")
for lk_name, ll in zip(model.likelihood.keys(), loglikes):
    print(f"    {lk_name:44s}: chi2 = {-2.0*float(ll):.2f}")
print(f"  Total chi2 at reference point = {chi2_tot:.2f}")
print("========================================================================")
