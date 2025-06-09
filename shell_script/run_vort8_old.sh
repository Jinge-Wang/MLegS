make swipe
cp ./converg/new_perturb.input ./
cp ./converg/cor_perturb.input ./

srun --mpi=pmi2 -n 8 ./bin/init_exec
{   echo "new_perturb" 
    echo "T"
} | srun --mpi=pmi2 -n 8 ./bin/addperturb_non_exec
printf '%s\n' "new_perturb" "T" | srun --mpi=pmi2 -n 8 ./bin/addperturb_non_exec
printf '%s\n' "cor_perturb" "F" "T" | srun --mpi=pmi2 -n 8 ./bin/addperturb_non_exec
echo "cor_perturb" | srun --mpi=pmi2 -n 8 ./bin/vort8_exec

# printf '%s\n' "new_perturb" "T" | srun --mpi=pmi2 -n 8 ./bin/addperturb_non_exec
# echo "new_perturb" | srun --mpi=pmi2 -n 8 ./bin/vort8_exec

# ---------------------------------------------------
