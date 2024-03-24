#!/bin/bash
#$ -cwd
#$ -N st.obelix
#$ -pe mpi 64
#$ -S /bin/bash
## #$ -q ibnet ##
#$ -l h_rt=24:00:00
#$ -l h=!(compute-0-37.local)
#$ -o slurm.st.out
#$ -e sge.err

module load python-3.6

echo Running on host `hostname` on queue $QUEUE
echo "SGE job id: $JOB_ID" 
echo Time is `date`
echo Directory is `pwd`
echo This job has allocated $NSLOTS processors

python FPE_St.py
