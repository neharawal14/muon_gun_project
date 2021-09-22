from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'SingleMuPt100_3rdfull_1st_large_step2'
config.General.workArea = 'crab_projects_pt100_3rdfull_1st_withoutPU'
config.General.transferOutputs = True
config.General.failureLimit=1
config.General.transferLogs=True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'SingleMuPt100_3rdfull_1st_GEN_SIM_DIGI_L1_DIGI2RAW.py'

#config.Data.outputPrimaryDataset = 'MinBias'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#NJOBS = 2000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True
#config.Data.outputDatasetTag = 'ZMuMu_step2'
#config.Data.outputPrimaryDataset = 'ZMuMu_step2_ferrico'
config.Data.inputDataset = '/SingleMuPt100_3rdfull_1st_large_GEN-SIM_step1_neha_withoutPU/nrawal-SingleMuPt100_3rdfull_1st_large_GEN-SIM_step1_withoutPU-c0efa3a78dca96311b96a607e84a72ab/USER'
#config.Data.inputDataset = '/ZMuMu_GEN-SIM_step1_ferrico/ferrico-ZMuMu_GEN-SIM_step1-9eadee95878022f078e16d6b70fe376c/USER'
config.Data.inputDBS = 'phys03'

config.Site.storageSite = 'T2_US_Florida'
