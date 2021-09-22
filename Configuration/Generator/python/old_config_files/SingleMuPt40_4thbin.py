import FWCore.ParameterSet.Config as cms
generator = cms.EDFilter("Pythia8PtGun",
                         PGunParameters = cms.PSet(
        MaxPt = cms.double(40.01),
        MinPt = cms.double(39.99),
        ParticleID = cms.vint32(-13),
        AddAntiParticle = cms.bool(False),
        MaxEta = cms.double(0.8),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(0.6),
        MinPhi = cms.double(-3.14159265359) ## in radians
        ),
                         Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
                         psethack = cms.string('single mu pt 40 4th bin'),
                         firstRun = cms.untracked.uint32(1),
                         PythiaParameters = cms.PSet(parameterSets = cms.vstring())
                         )
