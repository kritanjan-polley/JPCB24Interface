#!/bin/bash
start_time="$(date -u +%s.%N)"

# files=semi.3.txt
endtime=30

for i in $(seq 0 1 ${endtime})
do
files=semi.$i.txt
if test -f "$files"; then
   echo "file already exist for "$i""
else
   touch semi.${i}.txt
   echo 'Time Econserve TotEng KinEng PotEng Temp c_slabCOM[3] c_o3COM[3]' >> semi.${i}.txt
fi
done

echo ' Are you sure?[y/n]'
read decision
if [ "$decision" == "n" ]; then
   echo 'ok then...'
else
   for i in $(seq 0 1 ${endtime})
   do 
     picosecond=$((i*1000))
     echo "${picosecond}"
     for j in intermediate.shooting.*.out
     do
       # awk '$1 ~ '"/^${picosecond}/" ${j} >> semi.${i}.txt
       grep  '^[[:blank:]]'${picosecond}'[[:blank:]]\s' ${j} >> semi.${i}.txt
     done
   awk '!seen[$0]++' semi.${i}.txt > foo && mv foo semi.${i}.txt
   done
fi 
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "It took $elapsed seconds"
