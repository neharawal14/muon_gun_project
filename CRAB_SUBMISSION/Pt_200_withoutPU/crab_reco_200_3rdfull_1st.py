from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'SingleMuPt200_3rdfull_1st_large_step3'
config.General.workArea = 'crab_projects_pt200_3rdfull_1st_withoutPU'
config.General.transferOutputs = True
config.General.failureLimit=1
config.General.transferLogs=True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'SingleMuPt200_3rdfull_1st_GEN_DIGI_L1_RAW2DIGI_RECO.py'

#config.Data.outputPrimaryDataset = 'MinBias'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.publication = True
config.Data.inputDataset = '/SingleMuPt200_3rddfull_1st_GEN-SIM_step1_neha/nrawal-crab_SingleMuPt200_3rdfull_1st_large_step2-b92895268740895c325a8071b664b5b5/USER'
#config.Data.inputDataset = '/ZMuMu_GEN-SIM_step1_ferrico/ferrico-ZMuMu_GEN-SIM_step1-9eadee95878022f078e16d6b70fe376c/USER'
config.Data.inputDBS = 'phys03'

config.Site.storageSite = 'T2_US_Florida'
