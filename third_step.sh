#!/bin/sh
cmsDriver.py SingleMuPt${1}_${2}bin_GEN_DIGI_L1 --runUnscheduled --conditions auto:phase1_2018_realistic -s RAW2DIGI,L1Reco,RECO -n 10000 --era Run2_2018 --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-RECO
