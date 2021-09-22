#!/bin/bash
pt_list=(75 100) 
folder_list=(75 100)
mv Muon_gun_pt75_1stbin_full.root /afs/cern.ch/user/n/nrawal/eos/analyze_programs/SingleMuPt75/fullbin/
for i in ${!pt_list[*]} 
do
#mv Muon_gun_pt${pt_list[$i]}_1stbin_full.root /afs/cern.ch/user/n/nrawal/eos/analyze_programs/SingleMuPt${folder_list[$i]}/fullbin/
mv Muon_gun_pt${pt_list[$i]}_2ndbin_full.root /afs/cern.ch/user/n/nrawal/eos/analyze_programs/SingleMuPt${folder_list[$i]}/fullbin_2nd/
mv Muon_gun_pt${pt_list[$i]}_3rdbin_full.root /afs/cern.ch/user/n/nrawal/eos/analyze_programs/SingleMuPt${folder_list[$i]}/fullbin_3rd/
#mv Muon_gun_pt${pt_list[$i]}_2ndbin_small.root /afs/cern.ch/user/n/nrawal/eos/analyze_programs/SingleMuPt${folder_list[$i]}/2ndbin_small/
#mv Muon_gun_pt${pt_list[$i]}_3rdbin_small.root /afs/cern.ch/user/n/nrawal/eos/analyze_programs/SingleMuPt${folder_list[$i]}/3rdbin_small/
#mv Muon_gun_pt${pt_list[$i]}_4thbin_small.root /afs/cern.ch/user/n/nrawal/eos/analyze_programs/SingleMuPt${folder_list[$i]}/4thbin_small/
done
