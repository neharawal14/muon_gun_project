#!/bin/bash
pt_bin=(5 10 20 40 50 75 100)
for i in ${!pt_bin[*]}
do
  cp SingleMuPt${pt_bin[$i]}_2ndfull.py SingleMuPt${pt_bin[$i]}_3rdfull.py
  sed -i '14 s/2nd full/3rd full/' SingleMuPt${pt_bin[$i]}_3rdfull.py
  sed -i '8 s/1.4/2.4/' SingleMuPt${pt_bin[$i]}_3rdfull.py
  sed -i '10 s/0.9/1.4/' SingleMuPt${pt_bin[$i]}_3rdfull.py
done
