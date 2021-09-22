from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config

config = config()

config.General.requestName = 'SingleMuPt50_1stfull_GEN-SIM_step1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.failureLimit=1
config.General.transferLogs=True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'SingleMuPt50_full_GEN_SIM.py'

#config.Data.outputPrimaryDataset = 'MinBias'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 50
NJOBS = 200  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True
config.Data.outputDatasetTag = 'SingleMuPt50_1stfull_GEN-SIM_step1'
config.Data.outputPrimaryDataset = 'SingleMuPt50_1stfull_GEN-SIM_step1_neha'
config.Site.storageSite = 'T2_US_Florida'
