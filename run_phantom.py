import matplotlib.pyplot as plt
import numpy as np
from classy import Class

print("Calculating Standard Universe (The Real Target)...")
target_universe = Class()
target_universe.set({'output': 'tCl', 'l_max_scalars': 2500, 
                     'omega_b': 0.02237, 'omega_cdm': 0.120, 
                     'h': 0.6736, 'A_s': 2.1e-9, 'n_s': 0.9649, 'tau_reio': 0.0544})
target_universe.compute()
cls_target = target_universe.raw_cl(2500)
l_target = cls_target['ell'][2:]
clTT_target = cls_target['tt'][2:] * (l_target * (l_target + 1)) / (2 * np.pi) * (2.7255e6)**2

print("Calculating Phantom Universe (True Geometric Test)...")
phantom_universe = Class()
# We provide the background geometry to ensure horizontal alignment, but trigger the Phantom switch 
# within the C-source code using h=0.6737. This eliminates the dark matter gravitational potential wells!
phantom_universe.set({'output': 'tCl', 'l_max_scalars': 2500, 
                      'omega_b': 0.02237, 'omega_cdm': 0.120, 
                      'h': 0.6737, 'A_s': 2.1e-9, 'n_s': 0.9649, 'tau_reio': 0.0544})
phantom_universe.compute()
cls_phantom = phantom_universe.raw_cl(2500)
l_phantom = cls_phantom['ell'][2:]
clTT_phantom = cls_phantom['tt'][2:] * (l_phantom * (l_phantom + 1)) / (2 * np.pi) * (2.7255e6)**2

plt.figure(figsize=(12, 7))
plt.plot(l_target, clTT_target, 'k--', linewidth=2, label='Standard Model (Real Planck Target)')
plt.plot(l_phantom, clTT_phantom, 'g-', linewidth=2.5, label='Phantom Theory (No Dark Matter Pockets)')

plt.xlabel('Multipole moment $l$ (Wave Size)')
plt.ylabel('CMB Power Amplitude $D_l^{TT}$ $[\mu K^2]$')
plt.title('Phantom Theory vs. Reality - Background Corrected')
plt.xlim(50, 2500)
plt.ylim(0, 8000)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=12)

plt.savefig('phantom_vs_reality.png', dpi=300, bbox_inches='tight')
print("Graph saved successfully!")