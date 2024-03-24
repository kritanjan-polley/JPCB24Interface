#!/bin/bash
start_time=$(gdate +%s.%N)

for i in $(seq 1 1 30)
do
echo $i
awk  'NF<3' probability_${i}.txt  > temp
sed '/^[[:blank:]]*$/ d' temp > temp2
awk '{print $2}' temp2 > probability_${i}.txt
done

rm temp temp2

end_time=$(gdate +%s.%N)
elapsed=$(echo "scale=9; $end_time - $start_time" | bc)
echo 'It took' $elapsed 'seconds'