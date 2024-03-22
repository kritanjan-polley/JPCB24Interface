#!/bin/bash
start_time="$(date -u +%s.%N)"

j=0
for i in intermediate.shooting.*
do 
j=$(($j + 1))
awk '/[[:blank:]]*Time [[:blank:]]*Econserve/{flag=1; next} /[[:blank:]]*Loop [[:blank:]]*time/{flag=0} flag' $i >  foo
# awk '{print $1,$8}' foo >  bar ##ca.$j.out
grep '[0-9][0-9]*' foo > $i
rm foo
done

end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "It took $elapsed seconds"
