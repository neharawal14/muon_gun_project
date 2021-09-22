#!/bin/sh
#cmsenv
#cmsDriver.py Configuration/Generator/python/SingleMuPt${1}_${2}bin --conditions auto:phase1_2018_realistic -n 10000 --era Run2_2018 --eventcontent FEVTDEBUG --relval 9000,100 -s GEN,SIM --datatier GEN-SIM --beamspot Realistic25ns13TeVEarly2018Collision
cmsDriver.py SingleMuPt${1}_${2}bin_GEN --conditions auto:phase1_2018_realistic -s DIGI:pdigi_valid,L1,DIGI2RAW -n 10000 --era Run2_2018 --eventcontent FEVTDEBUGHLT
#cmsDriver.py SingleMuPt${1}_${2}bin_GEN_DIGI_L1 --runUnscheduled --conditions auto:phase1_2018_realistic -s RAW2DIGI,L1Reco,RECO -n 10000 --era Run2_2018 --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-RECO
