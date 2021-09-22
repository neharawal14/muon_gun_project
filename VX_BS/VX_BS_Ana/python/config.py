import FWCore.ParameterSet.Config as cms

from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.GlobalTag.globaltag='102X_upgrade2018_realistic_v18'

process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet(

	SkipEvent = cms.untracked.vstring('ProductNotFound')

	)

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root',' with the source file you want to use
				fileNames = cms.untracked.vstring(
#			     'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5/Endcap/SingleMuPt5_endcap_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
           'file:/afs/cern.ch/user/n/nrawal/work/Muon_gun_Vx_BS_latest/CMSSW_11_0_3/src/CRAB_SUBMISSION/Pt_5/SingleMuPt5_3rdfull_2nd_GEN_DIGI_L1_PU_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/1stbin/SingleMuPt5_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/2ndbin/SingleMuPt5_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/3rdbin/SingleMuPt5_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/4thbin/SingleMuPt5_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/5thbin/SingleMuPt5_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/1stbin_full/SingleMuPt5_1stbin_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/2ndbin/SingleMuPt10_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/1stbin/SingleMuPt10_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/3rdbin/SingleMuPt10_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/4thbin/SingleMuPt10_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/5thbin/SingleMuPt10_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/1stbin_full/SingleMuPt10_1stbin_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#			    'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/Endcap/SingleMuPt20_endcap_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',


#          'file:/afs/cern.ch/user/n/nrawal/work/Muon_gun_Vx_BS_latest/CMSSW_11_0_3/src/SingleMuPt50_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#        'file:/afs/cern.ch/user/n/nrawal/work/Muon_gun_Vx_BS_latest/CMSSW_11_0_3/src/output.root',

#         'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/SingleMuPt5_new_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root' 


#          'file:/afs/cern.ch/user/n/nrawal/work/Muon_gun_Vx_BS_latest/CMSSW_11_0_3/src/SingleMuPt50_2ndfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
          #    Pt : 5 
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/1stbin/SingleMuPt5_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/2ndbin/SingleMuPt5_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/3rdbin/SingleMuPt5_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/4thbin/SingleMuPt5_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/5thbin/SingleMuPt5_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',


#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/1stbin_new/SingleMuPt5_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/2ndbin_new/SingleMuPt5_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/3rdbin_new/SingleMuPt5_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/4thbin_new/SingleMuPt5_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/5thbin_new/SingleMuPt5_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

# For bigger full bins like Filippo

#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/fullbin/SingleMuPt5_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/fullbin_2nd/SingleMuPt5_2ndfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/fullbin_3rd/SingleMuPt5_3rdfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

# small bins for eta<0.2
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/2nd_small/SingleMuPt5_2ndsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/3rd_small/SingleMuPt5_3rdsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/4th_small/SingleMuPt5_4thsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

#           Pt : 10

#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/1stbin_new/SingleMuPt10_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/2ndbin_new/SingleMuPt10_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/3rdbin_new/SingleMuPt10_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/4thbin_new/SingleMuPt10_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/5thbin_new/SingleMuPt10_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

# For bigger full bins like Filippo

#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/fullbin/SingleMuPt10_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/fullbin_2nd/SingleMuPt10_2ndfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/fullbin_3rd/SingleMuPt10_3rdfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

# small bins for eta<0.2
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/2nd_small/SingleMuPt10_2ndsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/3rd_small/SingleMuPt10_3rdsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/4th_small/SingleMuPt10_4thsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',






#			    'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5/Endcap/SingleMuPt5_endcap_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
# 			    'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/Endcap/SingleMuPt50_endcap_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',


#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/1stbin/SingleMuPt20_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/2ndbin/SingleMuPt20_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/3rdbin/SingleMuPt20_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/4thbin/SingleMuPt20_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/5thbin/SingleMuPt20_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

# for bigger full bins like Filippo

#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/fullbin/SingleMuPt20_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/fullbin_2nd/SingleMuPt20_2ndfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/fullbin_3rd/SingleMuPt20_3rdfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

# small bins for eta<0.2
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/2nd_small/SingleMuPt20_2ndsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/3rd_small/SingleMuPt20_3rdsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/4th_small/SingleMuPt20_4thsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',


#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt30/1stbin/SingleMuPt30_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt30/2ndbin/SingleMuPt30_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt30/3rdbin/SingleMuPt30_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt30/4thbin/SingleMuPt30_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt30/5thbin/SingleMuPt30_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',


#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/1stbin/SingleMuPt40_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/2ndbin/SingleMuPt40_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#			       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/3rdbin/SingleMuPt40_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/4thbin/SingleMuPt40_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#	       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/5thbin/SingleMuPt40_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',



# for bigger full bins like Filippo

#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/fullbin/SingleMuPt40_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/fullbin_2nd/SingleMuPt40_2ndfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/fullbin_3rd/SingleMuPt40_3rdfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

# small bins for eta<0.2
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/2nd_small/SingleMuPt40_2ndsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/3rd_small/SingleMuPt40_3rdsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/4th_small/SingleMuPt40_4thsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',







#	        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/1stbin/SingleMuPt50_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

#	        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/2ndbin/SingleMuPt50_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

#        	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/3rdbin/SingleMuPt50_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
        
#        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/4thbin/SingleMuPt50_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
      
#        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/5thbin/SingleMuPt50_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

# for bigger full bins like Filippo

#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/fullbin/SingleMuPt50_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/fullbin_2nd/SingleMuPt50_2ndfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/fullbin_3rd/SingleMuPt50_3rdfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

# small bins for eta<0.2
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/2nd_small/SingleMuPt50_2ndsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/3rd_small/SingleMuPt50_3rdsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/4th_small/SingleMuPt50_4thsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',





#	        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt75/1stbin/SingleMuPt75_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

#	        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt75/2ndbin/SingleMuPt75_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

#        	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt75/3rdbin/SingleMuPt75_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
        
#        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt75/4thbin/SingleMuPt75_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

      
#        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt75/5thbin/SingleMuPt75_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',


# for bigger full bins like Filippo

#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt75/fullbin/SingleMuPt75_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt75/fullbin_2nd/SingleMuPt75_2ndfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt75/fullbin_3rd/SingleMuPt75_3rdfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',




#	        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/1stbin/SingleMuPt100_1stbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

#	        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/2ndbin/SingleMuPt100_2ndbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

#        	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/3rdbin/SingleMuPt100_3rdbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
        
#        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/4thbin/SingleMuPt100_4thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
      
#        'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/5thbin/SingleMuPt100_5thbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',


# for bigger full bins like Filippo

#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/fullbin/SingleMuPt100_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/fullbin_2nd/SingleMuPt100_2ndfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/fullbin_3rd/SingleMuPt100_3rdfull_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',




#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/fullbin/SingleMuPt5_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/fullbin/SingleMuPt10_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/fullbin/SingleMuPt20_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/fullbin/SingleMuPt40_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/fullbin/SingleMuPt50_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/fullbin/SingleMuPt100_full_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',



#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/2nd_small/SingleMuPt5_2ndsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/3rd_small/SingleMuPt5_3rdsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#          'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/4th_small/SingleMuPt5_4thsmallbin_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',


#			    'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt5_local/SingleMuPt5_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

#         'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_10kevents_FLO_BSset_final/SingleMuPt10_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',      
#			    'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_20kevents_FLO_BSet_final/SingleMuPt10_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',      

#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt20/SingleMuPt20_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root'
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt30/SingleMuPt30_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root'
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt40/SingleMuPt40_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root'
#               'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt50_10kevents_FLO/SingleMuPt50_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',

#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt75/SingleMuPt75_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root'
# 			    'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt100/SingleMuPt100_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
# 	   		    'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt200_local/SingleMuPt200_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt300_correct/SingleMuPt300_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root'
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt400/SingleMuPt400_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root'
#		       	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt600/SingleMuPt600_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root'
#	        	'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt1000/SingleMuPt1000_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root'





#				'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt200/again/SingleMuPt200_new_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#               'file:/afs/cern.ch/user/n/nrawal/eos/ntuple_muongun/SingleMuPt10_again_crab/SingleMuPt10_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',	                                          		



#						    'file:/afs/cern.ch/work/n/nrawal/Muon_gun_Vx_BS_latest/CMSSW_11_0_3/src/SingleMuPt200_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
#                                                    'file:/afs/cern.ch/work/n/nrawal/Muon_gun_Vx_BS_latest/CMSSW_11_0_3/src/ZMM_13TeV_TuneCUETP8M1_cfi_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root',
# 						    'file:/afs/cern.ch/work/f/ferrico/private/VxBS_Constraint_11/CMSSW_11_0_3/src/MuonGun/step3_RAW2DIGI_L1Reco_RECO_3000.root',',
                                                    #'file:/afs/cern.ch/work/n/nrawal/Muon_gun_Vx_BS_latest/CMSSW_11_0_3/src/SingleMuPt10_pythia8_cfi_GEN_DIGI_L1_RAW2DIGI_L1Reco_RECO.root'
# 						    'file:/afs/cern.ch/work/f/ferrico/private/VxBS_Constraint_11/CMSSW_11_0_3/src/ZMuMu/step3_RAW2DIGI_L1Reco_RECO.root',
############# 1cm
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_99.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_98.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_97.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_96.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_95.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_94.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_93.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_92.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_91.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_90.root',',

# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_89.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_88.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_87.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_86.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_85.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_84.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_83.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_82.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_81.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_80.root',',

# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_79.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_78.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_77.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_76.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_75.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_74.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_73.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_72.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_71.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_70.root',',

# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_69.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_68.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_67.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_66.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_65.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_64.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_63.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_62.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_61.root',',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200916_143215/0000/step3_RAW2DIGI_L1Reco_RECO_60.root',',
############# 1cm

# 						    'file:/afs/cern.ch/work/f/ferrico/private/VxBS_Constraint_11/CMSSW_11_0_3/src/MuonGun/step3_RAW2DIGI_L1Reco_RECO.root',',


############# 10um
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_86.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_8.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_79.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_78.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_76.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_74.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_65.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_61.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_59.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_5.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_45.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_40.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_4.root',
## 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_37.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_21.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_200.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_2.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_197.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_196.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_186.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_185.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_184.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_182.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_181.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_180.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_18.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_179.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_163.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_162.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_160.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_158.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_157.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_156.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_155.root',
## 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_153.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_151.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_150.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_15.root',
##  						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_147.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_144.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_142.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_140.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_139.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_137.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_133.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_132.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_125.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_124.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_119.root',
# 						    '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200925_132943/0000/step3_RAW2DIGI_L1Reco_RECO_5u_1cm_117.root',
############## 10um




############# 2.5x10-5
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_84.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_95.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_92.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_91.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_90.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_89.root',
	# # '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_88.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_87.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_86.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_85.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_191.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_190.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_189.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_188.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_187.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_186.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_185.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_184.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_183.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_182.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_181.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_180.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_179.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_178.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_177.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_176.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_175.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_174.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_173.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_172.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_170.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_169.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_168.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_167.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_165.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_164.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_163.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_161.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_160.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_158.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_157.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_156.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_155.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_154.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_153.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_152.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_151.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_150.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_149.root',
# '/store/user/ferrico/ZMuMu_GEN-SIM_step1_ferrico/crab_ZMuMu_step3/200927_090244/0000/step3_RAW2DIGI_L1Reco_RECO_xyzError_148.root',
############# 2.5x10-5

                )
                            )

process.TFileService = cms.Service("TFileService",
#			            fileName = cms.string("Muon_gun_pt5.root")
#    				    	fileName = cms.string("Muon_gun_pt10.root")
#			            fileName = cms.string("Muon_gun_pt20.root")
#			            fileName = cms.string("Muon_gun_pt30.root")
#			            fileName = cms.string("Muon_gun_pt40.root")
#                 fileName = cms.string("Muon_gun_pt50.root")
#	                fileName = cms.string("Muon_gun_pt75.root")
#      				    fileName = cms.string("Muon_gun_pt100.root")
#	       			    fileName = cms.string("Muon_gun_pt200.root")
#	       			    fileName = cms.string("Muon_gun_pt300.root")
#	       			    fileName = cms.string("Muon_gun_pt400.root")
#	       			    fileName = cms.string("Muon_gun_pt600.root")
#                 fileName = cms.string("Muon_gun_pt1000.root")

#                  fileName = cms.string("Muon_gun_pt5_another.root")

#                 fileName = cms.string("Muon_gun_pt5_1stbin.root")
#                  fileName = cms.string("Muon_gun_pt5_2ndbin.root")
#                 fileName = cms.string("Muon_gun_pt5_3rdbin.root")
#                 fileName = cms.string("Muon_gun_pt5_4thbin.root")
#                 fileName = cms.string("Muon_gun_pt5_5thbin.root")

#                  fileName = cms.string("Muon_gun_pt5_1stbin_full.root")
#                  fileName = cms.string("Muon_gun_pt5_2ndbin_full.root")
#                  fileName = cms.string("Muon_gun_pt5_3rdbin_full.root")

#                  fileName = cms.string("Muon_gun_pt5_2ndbin_small.root")
#                  fileName = cms.string("Muon_gun_pt5_3rdbin_small.root")
#                  fileName = cms.string("Muon_gun_pt5_4thbin_small.root")


#                 fileName = cms.string("Muon_gun_pt10_1stbin.root")
#                 fileName = cms.string("Muon_gun_pt10_2ndbin.root")
#                 fileName = cms.string("Muon_gun_pt10_3rdbin.root")
#                 fileName = cms.string("Muon_gun_pt10_4thbin.root")
#                 fileName = cms.string("Muon_gun_pt10_5thbin.root")
#                 fileName = cms.string("Muon_gun_pt20_1stbin.root")
#                 fileName = cms.string("Muon_gun_pt10_3rdbin.root")

#                  fileName = cms.string("Muon_gun_pt10_1stbin_full.root")
#                  fileName = cms.string("Muon_gun_pt10_2ndbin_full.root")
#                  fileName = cms.string("Muon_gun_pt10_3rdbin_full.root")

#                  fileName = cms.string("Muon_gun_pt10_2ndbin_small.root")
#                  fileName = cms.string("Muon_gun_pt10_3rdbin_small.root")
#                  fileName = cms.string("Muon_gun_pt10_4thbin_small.root")



#                 fileName = cms.string("Muon_gun_pt5_endcap.root")
#                 fileName = cms.string("Muon_gun_pt20_endcap.root")
#                 fileName = cms.string("Muon_gun_pt50_endcap.root")
#			       	    fileName = cms.string("Muon_gun_pt5_2000events_tmp.root")


#                  fileName = cms.string("Muon_gun_pt20_1stbin.root")
#                  fileName = cms.string("Muon_gun_pt20_2ndbin.root")
#                  fileName = cms.string("Muon_gun_pt20_3rdbin.root")
#                  fileName = cms.string("Muon_gun_pt20_4thbin.root")
#                  fileName = cms.string("Muon_gun_pt20_5thbin.root")

#                  fileName = cms.string("Muon_gun_pt20_1stbin_full.root")
#                  fileName = cms.string("Muon_gun_pt20_2ndbin_full.root")
#                  fileName = cms.string("Muon_gun_pt20_3rdbin_full.root")

#                  fileName = cms.string("Muon_gun_pt20_2ndbin_small.root")
#                  fileName = cms.string("Muon_gun_pt20_3rdbin_small.root")
#                  fileName = cms.string("Muon_gun_pt20_4thbin_small.root")




#                 fileName = cms.string("Muon_gun_pt30_1stbin.root")
#                 fileName = cms.string("Muon_gun_pt30_2ndbin.root")
#                 fileName = cms.string("Muon_gun_pt30_3rdbin.root")
#                 fileName = cms.string("Muon_gun_pt30_4thbin.root")
#                fileName = cms.string("Muon_gun_pt30_5thbin.root")

#                 fileName = cms.string("Muon_gun_pt40_1stbin.root")
#                 fileName = cms.string("Muon_gun_pt40_2ndbin.root")
#                 fileName = cms.string("Muon_gun_pt40_3rdbin.root")
#                 fileName = cms.string("Muon_gun_pt40_4thbin.root")
#                 fileName = cms.string("Muon_gun_pt40_5thbin.root")

#                  fileName = cms.string("Muon_gun_pt40_1stbin_full.root")
#                  fileName = cms.string("Muon_gun_pt40_2ndbin_full.root")
#                  fileName = cms.string("Muon_gun_pt40_3rdbin_full.root")

#                  fileName = cms.string("Muon_gun_pt40_2ndbin_small.root")
#                  fileName = cms.string("Muon_gun_pt40_3rdbin_small.root")
#                  fileName = cms.string("Muon_gun_pt40_4thbin_small.root")



#                 fileName = cms.string("Muon_gun_pt50_1stbin.root")
#                 fileName = cms.string("Muon_gun_pt50_2ndbin.root")
#                 fileName = cms.string("Muon_gun_pt50_3rdbin.root")
#                 fileName = cms.string("Muon_gun_pt50_4thbin.root")
#                 fileName = cms.string("Muon_gun_pt50_5thbin.root")

#                  fileName = cms.string("Muon_gun_pt50_1stbin_full.root")
#                  fileName = cms.string("Muon_gun_pt50_2ndbin_full.root")
#                  fileName = cms.string("Muon_gun_pt50_3rdbin_full.root")

#                  fileName = cms.string("Muon_gun_pt50_2ndbin_small.root")
#                  fileName = cms.string("Muon_gun_pt50_3rdbin_small.root")
#                  fileName = cms.string("Muon_gun_pt50_4thbin_small.root")



#                 fileName = cms.string("Muon_gun_pt75_1stbin.root")
#                 fileName = cms.string("Muon_gun_pt75_2ndbin.root")
#                 fileName = cms.string("Muon_gun_pt75_3rdbin.root")
#                 fileName = cms.string("Muon_gun_pt75_4thbin.root")
#                 fileName = cms.string("Muon_gun_pt75_5thbin.root")

#                  fileName = cms.string("Muon_gun_pt75_1stbin_full.root")
#                  fileName = cms.string("Muon_gun_pt75_2ndbin_full.root")
#                  fileName = cms.string("Muon_gun_pt75_3rdbin_full.root")


#                 fileName = cms.string("Muon_gun_pt100_1stbin.root")
#                 fileName = cms.string("Muon_gun_pt100_2ndbin.root")
#                  fileName = cms.string("Muon_gun_pt100_3rdbin.root")
#                 fileName = cms.string("Muon_gun_pt100_4thbin.root")
#                 fileName = cms.string("Muon_gun_pt100_5thbin.root")

#                  fileName = cms.string("Muon_gun_pt100_1stbin_full.root")
#                  fileName = cms.string("Muon_gun_pt100_2ndbin_full.root")
#                  fileName = cms.string("Muon_gun_pt100_3rdbin_full.root")

#                  fileName = cms.string("Muon_gun_pt100_2ndbin_small.root")
#                  fileName = cms.string("Muon_gun_pt100_3rdbin_small.root")
#                  fileName = cms.string("Muon_gun_pt100_4thbin_small.root")




#                 fileName = cms.string("Muon_gun_pt5_fullbin.root")
#                 fileName = cms.string("Muon_gun_pt10_fullbin.root")
#                 fileName = cms.string("Muon_gun_pt20_fullbin.root")
#                 fileName = cms.string("Muon_gun_pt40_fullbin.root")
#                 fileName = cms.string("Muon_gun_pt50_fullbin.root")
#                 fileName = cms.string("Muon_gun_pt100_fullbin.root")

#                 fileName = cms.string("Muon_gun_pt50_full2_PU.root")
#                 fileName = cms.string("Muon_gun_pt50_full2_newPU_100.root")
                 fileName = cms.string("output_5_3rdfull_2nd_tmp.root")
#                fileName = cms.string("Muon_gun_pt50_test_crab.root")
)


process.demo = cms.EDAnalyzer('VX_BS_Ana',
							
							prunedgenParticlesSrc = cms.untracked.InputTag("genParticles"),
							beamSpotSrc  = cms.untracked.InputTag("offlineBeamSpot"),
							vertexSrc    = cms.untracked.InputTag("offlinePrimaryVertices"),
							muonSrc = cms.untracked.InputTag("muons"),
							FLORIDASrc = cms.untracked.InputTag("globalMuonsFLORIDA","UpdatedAtVtx"),
  							globalSrc = cms.untracked.InputTag("globalMuonsFLORIDA"),
#  globalSrc : globalMuons used when analyzing only track information : globalMuonsFLORIDA is used for BS studies 
#							globalSrc = cms.untracked.InputTag("globalMuons"),
							stdAloneVtxSrc = cms.untracked.InputTag("standAloneMuons","UpdatedAtVtx"),
							stdAloneSrc = cms.untracked.InputTag("standAloneMuons"),

# 							globalSrc = cms.untracked.InputTag("globalMuons"),


                              )

process.p = cms.Path(process.demo)
