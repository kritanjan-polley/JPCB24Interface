#!/bin/bash

for i in $(seq 1 1 990)
do

cat > job_${i}.sh << REALEND
#!/bin/bash
#$ -cwd
#$ -N asterix_desorption
#$ -pe smp 1
#$ -S /bin/bash
##$ -q ${myPart} ## ibnet ##
#$ -o slurm_${i}.out
#$ -e sge.err
#$ -l h_rt=5:00:00
#$ -l h=!(compute-0-37.local)
REALEND


echo '

echo Running on host `hostname` on queue $QUEUE
echo "SGE job id: $JOB_ID" 
echo Time is `date`
echo Directory is `pwd`
echo This job has allocated $NSLOTS processors

' >> job_${i}.sh


cat >> job_${i}.sh << eof
~/julia-1.8.5/bin/julia -O3 checkJuliaSerial.jl
eof

qsub  job_${i}.sh

sleep 0.1
done
