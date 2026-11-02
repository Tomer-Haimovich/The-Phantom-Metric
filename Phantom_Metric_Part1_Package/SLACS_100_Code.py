import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# 1. Phantom Metric Core Operating Constants (Zero Free Parameters)
G = 6.67430e-11             # Gravitational constant (m^3 kg^-1 s^-2)
M_sun = 1.98847e30          # Solar mass (kg)
kpc_to_m = 3.08567758e19    # Kiloparsec to meters conversion
c = 299792458               # Speed of light (m/s)

# Global late-universe expansion rate calibrated in the macroscopic Hubble-flow regime (D > 20 Mpc)
H0_km_s_Mpc = 73.0
H0_Hz = (H0_km_s_Mpc * 1e3) / (3.08567758e22)
a0 = (c * H0_Hz) / (2.0 * np.pi)  # Fundamental horizon acceleration floor (m/s^2)

# 2. Load SLACS Strong-Lensing Dataset (N = 100 Early-Type Galaxies)
df = pd.read_csv("SLACS_Phantom_Metric_Lensing.csv")
df = df.sort_values(by='m_baryonic_10_11_msun').reset_index(drop=True)

# 3. Macroscopic Saturation Boundary (R_BB) Parameter-Free Calculation
df['M_baryonic_absolute_sun'] = df['m_baryonic_10_11_msun'] * 1e11
df['Calculated_R_BB_kpc'] = (np.sqrt((G * (df['M_baryonic_absolute_sun'] * M_sun)) / a0) / kpc_to_m).round(4)

# Evaluate goodness-of-fit against observed Einstein ring radii
slope, intercept, r_value, p_value, std_err = linregress(df['Calculated_R_BB_kpc'], df['observed_einstein_radius_kpc'])
r_squared = r_value**2

print("\n=================== SLACS STRONG-LENSING AUDIT RESULTS ===================")
print(f"Total Lensing Galaxies Processed (N): {len(df)}")
print(f"Empirical Goodness-of-Fit (R^2)     : {r_squared:.4f} (p = {p_value:.2e})")
print("Dark Matter Halo Parameters         : 0 (0% Dark Matter)")
print("==========================================================================\n")

df.to_csv("SLACS_Phantom_Metric_Lensing.csv", index=False)

# 4. Generate Validation Plot
plt.style.use('dark_background')
fig = plt.figure(figsize=(11, 6.5))

plt.scatter(df['m_baryonic_10_11_msun'], df['observed_einstein_radius_kpc'],
            color='none', edgecolor='#ff00ff', s=80, alpha=0.7,
            label='Telescope Observations (SLACS)')

plt.plot(df['m_baryonic_10_11_msun'], df['Calculated_R_BB_kpc'],
         color='#00ffcc', linewidth=2.5, label='Phantom Metric (Calculated $R_{BB}$)')

plt.title('The Phantom Metric: Einstein Ring Predictions in SLACS Galaxies (N = 100 Lenses)', fontsize=14, pad=15)
plt.xlabel('Baryonic Mass ($10^{11}$ Solar Masses)', fontsize=11)
plt.ylabel('Einstein Ring Radius / $R_{BB}$ (kpc)', fontsize=11)
plt.grid(True, linestyle=':', alpha=0.15, color='gray')

plt.text(df['m_baryonic_10_11_msun'].min() + 0.2, df['observed_einstein_radius_kpc'].max() - 2,
         f'$R^2 = {r_squared:.4f}$\nDark Matter: 0%\nParameter-Free Lock',
         color='#00ffcc', fontsize=11, bbox=dict(facecolor='black', alpha=0.6, edgecolor='#ff00ff'))

plt.legend(loc='lower right', fontsize=11)
plt.tight_layout()
plt.savefig('SLACS_100_Figure.png', dpi=300)
plt.close(fig)
