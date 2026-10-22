import pandas as pd
import numpy as np
from scipy.optimize import minimize_scalar
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

# 1. Load Data
# Ensure your dataset is named 'SPARC_Data.xlsx' and data is on the first sheet.
# Required columns: 'Galaxy', 'Rad_kpc', 'Vobs', 'Vgas', 'Vdisk', 'Vbulge'
file_path = "SPARC_Data.xlsx"
df = pd.read_excel(file_path, sheet_name=0)

# 2. Physics Constants
a0 = 1.128790e-10
phi = 1.618034
a_phi = a0 * phi

def calc_v_model(upsilon, row):
    V_gas = row['Vgas']
    V_disk = row['Vdisk']
    V_bulge = row['Vbulge']
    r_meters = row['Rad_kpc'] * 3.086e19
    
    # Calculate total baryonic velocity squared
    V_baryonic_sq = V_gas**2 + upsilon * V_disk**2 + upsilon * V_bulge**2
    V_baryonic_sq = max(0, V_baryonic_sq)
    V_baryonic = np.sqrt(V_baryonic_sq)
    
    # Calculate Newtonian baryonic acceleration
    a_baryonic = (V_baryonic * 1000)**2 / r_meters
    
    if a_baryonic == 0:
        return 0
    
    # Apply Topological Metric Deformation (nu)
    exponent = np.sqrt(a_baryonic / a_phi)
    nu = 1 / (1 - np.exp(-exponent))
    
    # Calculate predicted velocity
    V_mod = np.sqrt(V_baryonic_sq * nu)
    return V_mod

def error_for_galaxy(upsilon, galaxy_data):
    errors = []
    for _, row in galaxy_data.iterrows():
        v_mod = calc_v_model(upsilon, row)
        v_obs = row['Vobs']
        errors.append((v_mod - v_obs)**2)
    return np.sum(errors)

# 3. Optimize Mass-to-Light (M/L) ratio for each galaxy
optimal_upsilons = {}
optimized_preds = []
obs_all = []
df_optimized = df.copy()

print("Optimizing M/L for all galaxies... Please wait.")
for gal in df['Galaxy'].unique():
    gal_data = df[df['Galaxy'] == gal]
    
    # Restrict M/L between physical bounds (0.1 and 1.0)
    res = minimize_scalar(error_for_galaxy, bounds=(0.1, 1.0), args=(gal_data,), method='bounded')
    optimal_upsilons[gal] = res.x
    
    # Update DataFrame and Predictions
    for idx, row in gal_data.iterrows():
        opt_v = calc_v_model(res.x, row)
        optimized_preds.append(opt_v)
        obs_all.append(row['Vobs'])
        df_optimized.at[idx, 'Optimal_ML'] = res.x
        df_optimized.at[idx, 'V_Model_Opt'] = opt_v

# 4. Calculate Final R^2
r2_optimized = r2_score(obs_all, optimized_preds)
print(f"Optimization Complete! New R^2 is: {r2_optimized:.4f}")

# 5. Export new Excel with results
df_optimized.to_excel("SPARC_Optimized_ML_Results.xlsx", index=False)
print("Saved optimized data to 'SPARC_Optimized_ML_Results.xlsx'")

# 6. Generate Plot
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(12, 8))

# Scatter Plot
ax.scatter(obs_all, optimized_preds, color='mediumspringgreen', alpha=0.5, 
           label='SPARC Empirical Observations', edgecolors='teal', s=40)

# 1:1 Line
ax.plot([0, 400], [0, 400], color='magenta', linestyle='--', linewidth=2, label='Perfect Fit Line (1:1)')

# Formatting
ax.set_title(f'Universal Spacetime Manifold Verification\nConsolidated $R^2$ = {r2_optimized:.4f} (Variable Optimal M/L)', fontsize=14, pad=15)
ax.set_xlabel('Observed Velocity ($V_{obs}$ km/s)', fontsize=12)
ax.set_ylabel('Model Predicted Velocity ($V_{calc}$ km/s)', fontsize=12)
ax.grid(True, alpha=0.2)
ax.legend(loc='upper left', fontsize=10, framealpha=1, edgecolor='white')

# Save Plot
plt.savefig("SPARC_Optimized_R2.png", format='png', bbox_inches='tight', dpi=300)
print("Saved plot to 'SPARC_Optimized_R2.png'")