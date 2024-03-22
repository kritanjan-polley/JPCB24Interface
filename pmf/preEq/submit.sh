#!/bin/bash

for i in $(seq -35.0 1.0 35.0)
do

if (( $(echo "$i > 0.0" | bc -l) ))
then
numNodes=20
px=4
py=5
pz=1
else
numNodes=25
px=5
py=5
pz=1
fi

cat >job_eq_${i}.sh << REALEND
#!/bin/bash
#SBATCH -A m1864
#SBATCH -J Obelix_dang_preEq    # Job name
#SBATCH -o slurm_${i}.out       # output
#SBATCH -C cpu
#SBATCH --qos=shared
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${numNodes}
#SBATCH --cpus-per-task=2
#SBATCH --exclude=nid004151,nid004170,nid004158,nid004123,nid004159,nid004132,nid004162,nid004122,nid004157,nid004171,nid004109,nid004153,nid004160,nid004152

# module load python

echo "with ${numNodes} nodes and partitions px(${px}), py(${py}) and pz(${pz})"
REALEND


echo '
echo Running on host `hostname` at `date`
echo "slurm job id: $SLURM_JOB_ID"
' >> job_eq_${i}.sh


cat >> job_eq_${i}.sh << eof
env OMP_NUM_THREADS=2
mpirun -np ${numNodes} lmp_pm -sf omp -pk omp 2 -in min.slab.lmp -v myseed $RANDOM -v px ${px} -v py ${py} -v pz ${pz} -v location ${i}
bash prep_min_data.sed $i
python3 polarizer.py -q -f drude.dff data.slabPlumed.tip4p.${i}.min data.pol.slabPlumed.${i}.min
python3 ForceField.py
python3 get-pair-potentials.py
mpirun -np ${numNodes} lmp_pm -sf omp -pk omp 2 -in min.pol.slab.lmp -v myseed $RANDOM -v px ${px} -v py ${py} -v pz ${pz} -v location ${i}
rm data.slabPlumed.tip4p.${i}.min data.pol.slabPlumed.${i}.min
eof

sbatch --requeue  job_eq_${i}.sh

sleep 0.1
done
