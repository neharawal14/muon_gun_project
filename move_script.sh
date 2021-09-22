#!/bin/bash
#LIST={SingleMuPt40 SingleMuPt50 SingleMuPt75 SingleMuPt100}
LIST=5
#for i in {0}
#do
  for j in {1st,2nd,3rd,4th,5th}
  do
     mv SingleMuPt5_${j}bin_GEN_* /afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/${j}bin_new/
  done
#done
