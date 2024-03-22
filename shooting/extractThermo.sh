#!/bin/bash
start_time="$(date -u +%s.%N)"

j=0
for i in $(seq 0 1 50)
do 
  awk '{print $1,$8,$9}' semi.$i.txt > final.$i.txt  
done


end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "It took $elapsed seconds"
