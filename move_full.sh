#!/bin/bash
list=(5 10 20 40 50 75 100)
folders=(5_local 10_10kevents_FLO_BSset_final 20 40 50_10kevents_FLO 100)
for i in "${!folders[@]}"
do
 cp SingleMuPt${list[$i]}_2ndfull_GEN_* /afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt${folders[$i]}/fullbin_2nd/
 cp SingleMuPt${list[$i]}_3rdfull_GEN_* /afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt${folders[$i]}/fullbin_3rd/
done
