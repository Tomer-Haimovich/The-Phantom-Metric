#!/usr/bin/env bash
set -e

echo "========================================================================"
echo " LambdaTSM Real-Sky Cobaya Setup (Ubuntu / WSL)"
echo " Installs Cobaya + GetDist + Official Planck 2018 Native + DESI + SNe"
echo "========================================================================"

python3 --version
if ! python3 -c "import classy; print('Found classy:', classy.__file__)" 2>/dev/null; then
  echo "[WARNING] 'classy' is not currently importable in this Python environment."
  echo "          Activate the venv/conda environment where CLASS v3.4.0 is compiled,"
  echo "          or let cobaya-install install standard classy automatically."
fi

echo ""
echo "[Step 1/3] Installing/upgrading Cobaya and GetDist..."
pip install --upgrade cobaya getdist pyyaml scipy matplotlib

echo ""
echo "[Step 2/3] Downloading official Planck 2018 (native), DESI 2024 BAO, and Pantheon+ data..."
mkdir -p ./cobaya_packages ./chains ./results
cobaya-install ltsm_planck_desi_sn.yaml -p ./cobaya_packages --no-set-global

echo ""
echo "[Step 3/3] Running pre-flight consistency check..."
python3 00_test_classy_consistency.py

echo ""
echo "========================================================================"
echo " Setup Complete! Next step: run  bash 02_run_all_chains.sh"
echo "========================================================================"
