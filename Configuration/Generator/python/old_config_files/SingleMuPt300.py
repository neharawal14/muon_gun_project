import FWCore.ParameterSet.Config as cms
generator = cms.EDFilter("Pythia8PtGun",
                         PGunParameters = cms.PSet(
        MaxPt = cms.double(300.01),
        MinPt = cms.double(299.99),
        ParticleID = cms.vint32(-13),
        AddAntiParticle = cms.bool(False),
        MaxEta = cms.double(0.02),
        MaxPhi = cms.double(0.01),
        MinEta = cms.double(0.01),
        MinPhi = cms.double(-0.01) ## in radians
        ),
                         Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
                         psethack = cms.string('single mu pt 300'),
                         firstRun = cms.untracked.uint32(1),
                         PythiaParameters = cms.PSet(parameterSets = cms.vstring())
                         )
