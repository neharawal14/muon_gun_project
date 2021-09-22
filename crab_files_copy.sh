bin_number=(1stfull 2ndfull_1st 2ndfull_2nd 3rdfull_1st 3rdfull_2nd)
pt_value=(10 30 75 100 200 400)
for i in {0,1,2,3,4,5} 
do
  for j in {1stfull,2ndfull_1st,2ndfull_2nd,3rdfull_1st,3rdfull_2nd}
  do
    cp CRAB_SUBMISSION/Pt_5/crab_digi_5_${j}.py CRAB_SUBMISSION/Pt_${pt_value[$i]}/crab_digi_${pt_value[$i]}_${j}.py
    sed -i "s/5/${pt_value[$i]}/g" CRAB_SUBMISSION/Pt_${pt_value[$i]}/crab_digi_${pt_value[$i]}_${j}.py
  done
done
