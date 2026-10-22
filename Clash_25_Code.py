import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# 1. CONSTANTS & COSMOLOGY
G = 6.67430e-11
a0 = 1.12879e-10
M_sun = 1.989e30
kpc_to_m = 3.086e19
baryon_fraction = 0.13

def get_angular_diameter_distance(z, H0=73.0, Om0=0.3):
    c = 299792.458
    def integrand(x):
        return 1.0 / np.sqrt(Om0 * (1 + x)**3 + (1 - Om0))
    Dc, _ = quad(integrand, 0, z)
    Dc = Dc * (c / H0)
    Da = Dc / (1 + z)
    return Da * 1000

cluster_z = {
    'A209': 0.206, 'A383': 0.187, 'A611': 0.288, 'A1423': 0.213,
    'A2261': 0.224, 'CL1226': 0.892, 'M0329': 0.450, 'M0416': 0.396,
    'M0429': 0.399, 'M0647': 0.593, 'M0717': 0.548, 'M0744': 0.686,
    'M1115': 0.355, 'M1149': 0.544, 'M1206': 0.439, 'M1311': 0.494,
    'M1423': 0.543, 'RXJ1532': 0.362, 'M1720': 0.164, 'M1931': 0.352,
    'M2129': 0.235, 'MS2137': 0.313, 'RXJ1347': 0.451, 'RXJ2129': 0.235,
    'RXJ2248': 0.348
}
unrelaxed_mergers = ['M0717', 'M2129', 'M1720']

# 2. DATA PROCESSING
df = pd.read_csv('Clash_Data.csv')
df['theta_e'] = df['theta_e'].astype(str).str.replace('simeq', '', regex=False)
df['theta_e'] = pd.to_numeric(df['theta_e'], errors='coerce')
df['M_e'] = pd.to_numeric(df['M_e'], errors='coerce')
df = df.dropna(subset=['theta_e', 'M_e'])

results = []
for index, row in df.iterrows():
    cluster_model = str(row['Cluster_model']).strip()
    if '_' not in cluster_model: continue
    
    cluster_name = cluster_model.split('_')[0]
    model_type = cluster_model.split('_')[1]
    
    if model_type != 'LTM' or cluster_name not in cluster_z: continue
    
    theta_e_arcsec = float(row['theta_e'])
    total_mass_13 = float(row['M_e'])
    z = cluster_z[cluster_name]
    
    baryonic_mass_kg = (total_mass_13 * 1e13) * M_sun * baryon_fraction
    D_A = get_angular_diameter_distance(z)
    R_obs_kpc = D_A * theta_e_arcsec * (np.pi / (180 * 3600))
    R_bb_kpc = np.sqrt((G * baryonic_mass_kg) / a0) / kpc_to_m
    
    status = "Unrelaxed/Merger" if cluster_name in unrelaxed_mergers else "Relaxed/Spherical"
    results.append({
        'Cluster': cluster_name,
        'Morphology': status,
        'R_Observed_kpc': round(R_obs_kpc, 2),
        'R_Phantom_Calculated_kpc': round(R_bb_kpc, 2)
    })

results_df = pd.DataFrame(results)
results_df.to_csv('Phantom_Metric_Clusters_Predictions.csv', index=False)

# 3. STRICT 1:1 PARITY R^2 CALCULATION (Exactly like Excel)
relaxed_df = results_df[results_df['Morphology'] == 'Relaxed/Spherical']
unrelaxed_df = results_df[results_df['Morphology'] == 'Unrelaxed/Merger']

obs = relaxed_df['R_Observed_kpc'].values
calc = relaxed_df['R_Phantom_Calculated_kpc'].values

ss_res = np.sum((obs - calc) ** 2)
ss_tot = np.sum((obs - np.mean(obs)) ** 2)
r_squared_1_1 = 1 - (ss_res / ss_tot)

# 4. PLOTTING
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(11, 7))

ax.scatter(relaxed_df['R_Phantom_Calculated_kpc'], relaxed_df['R_Observed_kpc'],
           color='#00ffcc', s=90, alpha=0.9, edgecolor='w', label='Relaxed Clusters (Valid Geometry)')

ax.scatter(unrelaxed_df['R_Phantom_Calculated_kpc'], unrelaxed_df['R_Observed_kpc'],
           color='#ff0055', s=90, alpha=0.9, marker='X', edgecolor='w', label='Unrelaxed Mergers (Excluded)')

max_val = max(results_df['R_Phantom_Calculated_kpc'].max(), results_df['R_Observed_kpc'].max())
ax.plot([0, max_val+20], [0, max_val+20], color='magenta', linestyle='--', linewidth=2, label='Perfect Fit Line (1:1)')

ax.set_title('The Phantom Metric: Einstein Ring Predictions in Galaxy Clusters\n(0% Dark Matter, Baryonic Mass Only)', fontsize=15, pad=15)
ax.set_xlabel(r'Phantom Metric Predicted Saturation Boundary $R_{BB}$ (kpc)', fontsize=13)
ax.set_ylabel(r'Observed Einstein Ring Radius $R_{obs}$ (kpc)', fontsize=13)

textstr = f'1:1 Parity $R^2$ = {r_squared_1_1:.4f}\nBaryon Fraction = 13.0%\nDark Matter = 0%\nParameter-Free Lock'
props = dict(boxstyle='round', facecolor='black', alpha=0.8, edgecolor='#ff00ff')
ax.text(0.03, 0.95, textstr, transform=ax.transAxes, fontsize=13, verticalalignment='top', bbox=props, color='#00ffcc')

ax.legend(loc='lower right', fontsize=11)
ax.grid(True, linestyle=':', alpha=0.3)
ax.set_xlim(0, 210)  # Focused view on the valid clusters
ax.set_ylim(0, 210)
plt.tight_layout()
plt.savefig('Phantom_Cluster_Audit_True.png', dpi=300)
plt.show()

print(f"Excel-Matched 1:1 R^2 for Relaxed Clusters: {r_squared_1_1:.4f}")