#!/bin/bash


# for i in $(seq -25.0 1.0 25.0)
# for i in $(seq -28.0 1.0 -14.0) $(seq 14.0 1.0 28.0)
for i in -16.0 -25.0 16.0 24.0 4.0 9.0
do

if (( $(echo "$i > -40.0" | bc -l) ))
then
myPart=lr5
numNodes=14
px=2
py=7
pz=1
else
myPart=csd_lr6_96
numNodes=20
px=4
py=5
pz=1
fi

if [ "${myPart}" == "csd_lr6_96" ]
then
myqos=condo_statmech
else
myqos=lr_lowprio
fi


cat >plumedJob_NH_${i}.sh << REALEND
#!/bin/bash
#SBATCH -A lr_statmech
#SBATCH --partition=${myPart}
#SBATCH --output=slurm_${i}.out
#SBATCH -J Popeye_prod               # Job name
#SBATCH --nodes=1                 # Total number of mpi tasks requested
#SBATCH --ntasks=${numNodes}
#SBATCH --qos=${myqos}
#SBATCH -t 28:00:00           # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --exclude=n0110.lr4,n0070.lr4,n0008.lr4,n0008.lr4,n0029.lr4,n0019.lr4,n0031.lr4,n0030.lr4,n0028.lr4,n0018.lr4,n0007.lr4,n0135.lr5
#SBATCH --cpus-per-task=2

echo "with ${numNodes} nodes and partitions px(${px}), py(${py}) and pz(${pz})"

REALEND


echo '
echo Running on host `hostname` at `date`
echo "slurm job id: $SLURM_JOB_ID"
' >> plumedJob_NH_${i}.sh


cat >> plumedJob_NH_${i}.sh << eof
env OMP_NUM_THREADS=2
mpirun -np ${numNodes} lmp -sf omp -pk omp 2 -in prod.pol.slab.lmp -v location ${i} -v myseed $RANDOM -v px ${px} -v py ${py} -v pz ${pz}
eof

sbatch --requeue plumedJob_NH_${i}.sh 

sleep 0.1
done
