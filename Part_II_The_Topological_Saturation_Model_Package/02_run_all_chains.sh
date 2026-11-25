#!/usr/bin/env bash
set -e

NPROCS=${1:-4}
export OMP_NUM_THREADS=${2:-2}

echo "========================================================================"
echo " Running Matched Real-Sky Cobaya MCMC Audits (LambdaTSM k=5 vs LambdaCDM k=6)"
echo " MPI processes: ${NPROCS} | OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
echo "========================================================================"

run_yaml() {
  local yfile="$1"
  echo ""
  echo ">>> Starting: ${yfile} ($(date)) <<<"
  if command -v mpirun >/dev/null 2>&1 && python3 -c "import mpi4py" 2>/dev/null; then
    mpirun -n "${NPROCS}" cobaya-run "${yfile}"
  else
    echo "[INFO] mpirun or mpi4py not found; running single-process cobaya-run."
    echo "       For 4x faster parallel Gelman-Rubin chains, install: sudo apt install openmpi-bin && pip install mpi4py"
    cobaya-run "${yfile}"
  fi
  echo ">>> Finished: ${yfile} ($(date)) <<<"
}

# 1. Pure Planck 2018 (TTTEEE_lite_native + lowl.EE_sroll2 + lensing.native)
run_yaml ltsm_planck_only.yaml
run_yaml lcdm_planck_only.yaml

# 2. Combined Planck 2018 + DESI 2024 BAO + Pantheon+ SNe Ia
run_yaml ltsm_planck_desi_sn.yaml
run_yaml lcdm_planck_desi_sn.yaml

echo ""
echo "========================================================================"
echo " All 4 Cobaya MCMC runs finished! Running analysis & plot generator..."
echo "========================================================================"
python3 03_analyze_cobaya_results.py
