#! /bin/bash -x
#SBATCH --job-name="NSBsolver"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
###Multiply nodes and ntasks-per-node and use that for ntasks
#SBATCH --ntasks=48
### Karins Nodes specs: comp 25 has 28 cores, comp 178-179 ea. have 25 cores
##SBATCH --export=ALL
#SBATCH --time=002:59:00
##SBATCH --time=143:59:59
#SBATCH --exclusive
#SBATCH -o out_NSsolver.out
#SBATCH -e err_NSsolver.err
##SBATCH -p kleiderman
# Go to the directory from which our job was launched
cd $SLURM_SUBMIT_DIR

# Create a short JOBID basd on the one provided by the scheduler 
JOBID='echo $SLURM_JOBID'

# Save a copy of our environment and script
cat $0 > script.$JOBID
printenv > env.$JOBID

# Load the necessary modules
module load PrgEnv/devtoolset-9
module load MPI/openmpi/4.1.1/gcc
module load Lib/petsc/3.13-gnu-ompi-dbg

# Run the job
# The echo will go into the standard output for this job
# The standard output file will end up in the directory
# from which the job was launched
echo "running job"
srun run_main -u_sys_ksp_type bcgs -u_sys_ksp_rtol 1e-6 -u_sys_ksp_initial_guess_nonzero 1 -v_sys_ksp_type bcgs -v_sys_ksp_rtol 1e-6 -v_sys_ksp_initial_guess_nonzero 1 -p_sys_ksp_type preonly -p_sys_pc_type lu -p_sys_pc_factor_mat_solver_type mumps
echo "job has finished"
