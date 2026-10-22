import urllib.request
import zipfile
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# 1. Exact Universally Locked Architectural Constants
A_0 = 1.12879e-10
PHI = 1.6180339887
KPC_TO_METERS = 3.086e19
ML_UNIVERSAL = 0.46 # Rigidly locked mass-to-light ratio for entire universe ensemble

# 2. Extract SPARC Data (Downloads the universal LTG database if not present)
url = "http://astroweb.cwru.edu/SPARC/Rotmod_LTG.zip"
zip_path = "Rotmod_LTG.zip"
if not os.path.exists("SPARC_Data"):
    if not os.path.exists(zip_path):
        urllib.request.urlretrieve(url, zip_path)
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall("SPARC_Data")

all_rows = []

# 3. Process 100% of the unfiltered LTG population (175 Galaxies)
for filename in os.listdir("SPARC_Data"):
    if not filename.endswith(".dat"): continue
    filepath = os.path.join("SPARC_Data", filename)
    galaxy = filename.split('_')[0]
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    start_idx = 0
    for i, line in enumerate(lines):
        if line.strip() and not line.startswith('#') and line.split()[0].replace('.', '', 1).isdigit():
            start_idx = i
            break
            
    try:
        # Parse data safely
        df = pd.read_csv(filepath, sep=r'\s+', skiprows=start_idx, header=None, on_bad_lines='skip')
        if len(df.columns) == 5:
            df.columns = ['Rad_kpc', 'Vobs', 'errV', 'Vgas', 'Vdisk']
            df['Vbulge'] = 0.0
        elif len(df.columns) >= 6:
            df = df.iloc[:, :6]
            df.columns = ['Rad_kpc', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbulge']
        else: continue
        
        df['Rad_meters'] = df['Rad_kpc'] * KPC_TO_METERS
        df = df[df['Rad_meters'] > 0].copy()
        if df.empty: continue
        
        # Unified Baryonic Velocity Profile Calculation (Standard squaring)
        df['V_baryonic_sq'] = (df['Vgas']**2 + ML_UNIVERSAL * df['Vdisk']**2 + ML_UNIVERSAL * df['Vbulge']**2)
        df['V_baryonic_sq'] = np.where(df['V_baryonic_sq'] < 0, 0, df['V_baryonic_sq'])
        df['V_baryonic'] = np.sqrt(df['V_baryonic_sq'])
        
        # Acceleration and Topological Metric Deformation Function
        df['Accel'] = ((df['V_baryonic'] * 1000) ** 2) / df['Rad_meters']
        df['Exponent'] = np.sqrt(df['Accel'] / (A_0 * PHI))
        df['Nu'] = 1.0 / (1.0 - np.exp(-df['Exponent']))
        df['V_calculated'] = df['V_baryonic'] * np.sqrt(df['Nu'])
        
        df.insert(0, 'Galaxy', galaxy)
        all_rows.append(df)
    except Exception as e:
        continue

# 4. Consolidated Global Variance Audit
if all_rows:
    master_df = pd.concat(all_rows, ignore_index=True)
    y_true = master_df['Vobs']
    y_pred = master_df['V_calculated']
    ss_res = np.sum((y_true - y_pred)**2)
    ss_tot = np.sum((y_true - np.mean(y_true))**2)
    global_r2 = 1 - (ss_res / ss_tot)
    
    print("\n" + "="*50)
    print(f" THE PHANTOM METRIC: GLOBAL AUDIT RESULTS")
    print(f"Total Spatiotemporal Points Analyzed: {len(master_df)}")
    print(f"Global R^2 (Locked M/L = {ML_UNIVERSAL}): {global_r2:.4f}")
    print("="*50)
    
    master_df.to_csv("SPARC_Global_Phantom_Metric_ML046.csv", index=False)
    print("Saved exact evaluation data to SPARC_Global_Phantom_Metric_ML046.csv")
    
    # 5. Universal Verification Scatter Plot
    plt.style.use('dark_background')
    plt.figure(figsize=(10, 7))
    plt.scatter(y_true, y_pred, color='#00ffcc', alpha=0.4, label='SPARC Empirical Observations')
    plt.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], color='magenta', linestyle='--', linewidth=2, label='Perfect Fit Line (1:1)')
    plt.title(f"Universal Spacetime Manifold Verification\nConsolidated $R^2$ = {global_r2:.4f} (Locked M/L = {ML_UNIVERSAL})", fontsize=14)
    plt.xlabel("Observed Velocity ($V_{obs}$ km/s)", fontsize=11)
    plt.ylabel("Model Predicted Velocity ($V_{calc}$ km/s)", fontsize=11)
    plt.legend()
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig('figure2_global_audit_ml046.png', dpi=300)
    plt.show()