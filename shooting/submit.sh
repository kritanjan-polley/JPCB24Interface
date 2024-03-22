#!/bin/bash
#SBATCH -A lr_statmech
#SBATCH --partition=lr5
#SBATCH --output=/global/scratch/users/kpolley/shooting/NaICl/slurm.%a.%j.out
#SBATCH -J ozone_NaICl_shooting
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --qos=lr_lowprio
#SBATCH -t 0:20:00
#SBATCH --cpus-per-task=2
#SBATCH --array=1-1000

echo Running on host $SLURM_SUBMIT_HOST at `date` from directory $SLURM_SUBMIT_DIR
echo Slurm job id: $SLURM_JOB_ID, task id : $SLURM_ARRAY_TASK_ID, processor id : $SLURM_PROCID
echo Nodes assigned to this job : $SLURM_JOB_NODELIST
echo Number of CPUs for this job : $SLURM_CPUS_PER_TASK with $SLURM_NPROCS processors

module load python/3.6

numNodes=$SLURM_NPROCS
ncpu=$SLURM_CPUS_PER_TASK
numNodes=10
tps=30
px=2
py=5
pz=1

vel=`python gauss_velocity.py`
velZ=$vel
location=`python random.choice.py`
locationZ=$location
export omp_num_threads=${ncpu}
mpirun -np ${numNodes} lmp -sf omp -pk omp ${ncpu} -in prod.onlyWater.lmp -v myseed $RANDOM -v tps ${tps} -v px ${px} -v py ${py} -v pz ${pz} -v velZ ${velZ} -v loc ${locationZ}
