#!/bin/bash
# SBATCH script for strong-scaling test of STEP_BOUSSI_ETDAB_CN.
#
# The job allocates enough resources for the MAXIMUM proc count, then
# runs the same bsnsq_scaling_exec binary with decreasing -n values.
# Each srun call is an independent MPI job (its own MPI_INIT/FINALIZE).
# OMP_NUM_THREADS is forced to 1 for all runs to isolate MPI scaling.

#SBATCH -A phy220056
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -J bsnsq_scaling
#SBATCH -o bsnsq_scaling.o%j
#SBATCH -e bsnsq_scaling.e%j
#SBATCH -p shared
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=end

module load intel
module load impi
module list

# ---- Configuration -----------------------------------------------
MAX_PROCS=128            # must equal --nodes * --ntasks-per-node

# Proc counts to test (descending).  Remove entries larger than the
# total cores allocated, or add smaller ones as desired.
PROC_LIST="128 64 32 16 8 4 2"

# Force OMP to 1 thread for all runs, overriding any environment default
export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo ""

# ---- Clean output directory so each run starts fresh ---------------
rm -f ./output/*.dat ./output/*.info ./output/*.output

# ---- Scaling loop --------------------------------------------------
for NPROCS in $PROC_LIST; do

    # Skip if requested procs exceed what is allocated
    if [ "$NPROCS" -gt "$MAX_PROCS" ]; then
        echo ">>> Skipping NPROCS=$NPROCS (exceeds allocation of $MAX_PROCS)"
        continue
    fi

    echo "======================================================================="
    echo ">>> SCALING TEST: NPROCS = $NPROCS  (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
    echo "======================================================================="

    # Each srun is an independent MPI invocation; bsnsq_scaling does
    # its own MPI_INIT inside setup_environment and MPI_FINALIZE at end.
    srun --mpi=pmi2 \
         -n $NPROCS \
         ./bin/bsnsq_scaling_exec

    echo ">>> Done NPROCS = $NPROCS"
    echo ""

    # Brief pause between runs to let the system settle
    sleep 2

done

echo "======================================================================="
echo ">>> ALL SCALING RUNS COMPLETE"
echo "======================================================================="
