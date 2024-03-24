#!/bin/bash
#$ -cwd
#$ -N slurm.asterix_desorption
#$ -pe smp 1
#$ -t 1-499
#$ -S /bin/bash
# #$ -q ibnet ##
#$ -l h_rt=12:00:00
#$ -l h=!(compute-0-37.local)
#$ -o slurm.$TASK_ID.out
#$ -e sge.err

echo Running on host `hostname` on queue $QUEUE
echo "SGE job id: $JOB_ID" 
echo Time is `date`
echo Directory is `pwd`
echo This job has allocated $NSLOTS processors

~/julia-1.8.5/bin/julia -O3 checkJuliaSerial.jl
