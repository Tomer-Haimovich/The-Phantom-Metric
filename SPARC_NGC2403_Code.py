import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. Theoretical Constants (STRICTLY LOCKED)
A_0 = 1.12879e-10           # Base minimal acceleration floor (m/s^2)
PHI = 1.618033988749895     # Exact Topological Dispersion Ratio (phi)
TOPO_FACTOR = PHI           # Topological multiplier
ML_LOCKED = 0.46            # Rigidly locked M/L ratio (Zero local curve-fitting)
KPC_TO_METERS = 3.086e19    # Kiloparsec to meters conversion

def topological_velocity_model_phi(radius, v_gas, v_disk, v_bulge):
    # Calculate visible baryonic mass using the universally locked M/L ratio
    # Applying standard geometric squaring identically to Excel logic
    v_baryonic_sq = (v_gas**2 + ML_LOCKED * v_disk**2 + ML_LOCKED * v_bulge**2)
    v_baryonic_sq = np.where(v_baryonic_sq < 0, 0, v_baryonic_sq)
    v_baryonic = np.sqrt(v_baryonic_sq)

    # Calculate classical Newtonian acceleration
    radius_meters = radius * KPC_TO_METERS
    acceleration_baryonic = ((v_baryonic * 1000) ** 2) / radius_meters
    acceleration_baryonic = np.maximum(acceleration_baryonic, 1e-20)

    # Apply Topological Metric Deformation Function (nu) based on The Phantom Metric
    exponent = np.sqrt(acceleration_baryonic / (A_0 * TOPO_FACTOR))
    exponent = np.minimum(exponent, 50)  # Prevent overflow
    nu = 1.0 / (1.0 - np.exp(-exponent))

    # Calculate final observed velocity integrating the metric deformation
    v_final_sq = (v_baryonic ** 2) * nu
    return np.sqrt(np.maximum(v_final_sq, 0)), v_baryonic

def load_sparc_file_fixed(filename):
    # Parser for raw SPARC observational .dat files
    with open(filename, 'r') as f:
        lines = f.readlines()
    start_idx = 0
    for i, line in enumerate(lines):
        if line.strip() and not line.startswith('#') and line.split()[0].replace('.', '', 1).isdigit():
            start_idx = i
            break
    data = pd.read_csv(filename, sep=r'\s+', skiprows=start_idx, header=None)
    
    if len(data.columns) == 5:
        data.columns = ['Radius_kpc', 'V_observed', 'V_err', 'V_gas', 'V_disk']
        data['V_bulge'] = 0.0
    elif len(data.columns) >= 6:
        data = data.iloc[:, :6]
        data.columns = ['Radius_kpc', 'V_observed', 'V_err', 'V_gas', 'V_disk', 'V_bulge']
    return data

# 2. Execution and Calculation
df = load_sparc_file_fixed('NGC2403_rotmod.dat')
df['V_model_predicted'], df['V_baryonic_final'] = topological_velocity_model_phi(
    df['Radius_kpc'], df['V_gas'], df['V_disk'], df['V_bulge']
)

# 3. Statistical Validation (R-squared)
ss_res = np.sum((df['V_observed'] - df['V_model_predicted']) ** 2)
ss_tot = np.sum((df['V_observed'] - np.mean(df['V_observed'])) ** 2)
r_squared = 1 - (ss_res / ss_tot)

print(f"--- EXACT PLANAR DISPERSION TOPOLOGY (NGC 2403) ---")
print(f"Base Acceleration (a_0): {A_0} m/s^2 (LOCKED)")
print(f"Mass-to-Light Ratio: {ML_LOCKED} (LOCKED)")
print(f"Empirical Goodness-of-Fit (R^2): {r_squared:.4f}")

# Export the raw data to CSV for transparency
df.to_csv("NGC2403_Phantom_Metric_ML046.csv", index=False)
print("Saved data to NGC2403_Phantom_Metric_ML046.csv")

# 4. Plotting
plt.style.use('default')
plt.figure(figsize=(10, 6))
plt.plot(df['Radius_kpc'], df['V_observed'], 'ro', label='SPARC Telescope Observations')
plt.plot(df['Radius_kpc'], df['V_baryonic_final'], 'g--', label=f'Newtonian Baseline (M/L={ML_LOCKED})')
plt.plot(df['Radius_kpc'], df['V_model_predicted'], 'b', linewidth=2, label=f'Phantom Metric Framework (R^2 = {r_squared:.4f})')
plt.title(f'Galactic Rotation Curve: NGC 2403\nParameter-Free Verification (Locked M/L={ML_LOCKED})')
plt.xlabel('Radius (kpc)')
plt.ylabel('Rotation Velocity (km/s)')
plt.legend()
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig('figure1_ngc2403_locked.png', dpi=300)
plt.show()