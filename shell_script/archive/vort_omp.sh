#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Stampede2 KNL nodes
#
#   *** OpenMP Job on Normal Queue ***
# 
# Last revised: 20 Oct 2017
#
# Notes:
#
#   -- Launch this script by executing
#   -- Copy/edit this script as desired.  Launch by executing
#      "sbatch knl.openmp.slurm" on a Stampede2 login node.
#
#   -- OpenMP codes run on a single node (upper case N = 1).
#        OpenMP ignores the value of lower case n,
#        but slurm needs a plausible value to schedule the job.
#
#   -- Default value of OMP_NUM_THREADS is 1; be sure to change it!
#
#   -- Increase thread count gradually while looking for optimal setting.
#        If there is sufficient memory available, the optimal setting
#        is often 68 (1 thread per core) or 136 (2 threads per core).

#----------------------------------------------------

#SBATCH -J knl_strong           # Job name
#SBATCH -o knl_strong.o%j       # Name of stdout output file
#SBATCH -e knl_strong.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for OpenMP)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for OpenMP)
#SBATCH -t 00:05:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=jinge@berkeley.edu
#SBATCH --mail-type=all    # Send email at begin and end of job


# Other commands must follow all #SBATCH directives...

module list
pwd
date

# Set thread count (default value is 1)...

export OMP_NUM_THREADS=1

# Launch OpenMP code...

# Start Time
start_time="$(date -u +%s.%N)"

# vort2:
#echo "%read.input" | ./bin/copy_exec
#time ./bin/vort2_exec   # Do not use ibrun or any other MPI launcher

# vort:
make swipe
date
echo "%read.input" | ./bin/init_exec
echo "%read.input" | ./bin/copy_exec
echo "%read.input" | ./bin/vort_exec
date

end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Total of $elapsed seconds elapsed for process"

# ---------------------------------------------------

