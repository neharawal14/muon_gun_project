#!/bin/bash
min_pt=(4 9 29 74 99 199 399)
max_pt=(5 10 30 75 100 200 400)
len=${#min_pt[@]}
echo "here is the length $len"
for i in {0,1,2,3,4,5,6}
do
  echo "running i "$i 
  echo "ie. pt " ${max_pt[$i]}
  for j in {1stfull,2ndfull_1st,2ndfull_2nd,3rdfull_1st,3rdfull_2nd}
  do
    max=${max_pt[$i]}
    min=${min_pt[$i]}
    echo "running j $j"
    echo "PT50_${i} ex_${j}"
    echo "${max_pt[$i]}"
    cp SingleMuPt50_${j}.py SingleMuPt${max_pt[$i]}_${j}.py
    sed -i "s/50/$max/g" SingleMuPt${max_pt[$i]}_${j}.py
    sed -i "s/49/$min/g" SingleMuPt${max_pt[$i]}_${j}.py 
  done
done
