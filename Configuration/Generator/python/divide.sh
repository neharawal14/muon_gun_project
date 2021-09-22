#!/bin/bash
min_pt=(4 9 19 39 49 )
max_pt=(5 10 20 40 50 )
len=${#min_pt[@]}
echo "here is the length $len"
for i in {0,1,2,3,4}
do
  echo "running i $i"
  for j in {2nd,3rd,4th}
  do
    max=${max_pt[$i]}
    min=${min_pt[$i]}
    echo "running j $j"
    echo "PT20_${i} ex_${j}"
    echo "${max_pt[$i]}"
    cp SingleMuPt20_${j}smallbin.py SingleMuPt${max_pt[$i]}_${j}smallbin.py
    sed -i "s/20/$max/g" SingleMuPt${max_pt[$i]}_${j}smallbin.py
    sed -i "s/19/$min/g" SingleMuPt${max_pt[$i]}_${j}smallbin.py 
  done
done
