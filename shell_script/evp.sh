#!/bin/bash
# Set thread count (default value is 1)...
export OMP_NUM_THREADS=4

# Launch OpenMP code...
make swipe
echo "%read.input" | ./bin/init_exec
#./bin/evp_linear2_exec   # Do not use ibrun or any other MPI launcher
./bin/evp_linear_exec
#./bin/check_u2o_exec

# ---------------------------------------------------

