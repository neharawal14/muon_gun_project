import FWCore.ParameterSet.Config as cms
generator = cms.EDFilter("Pythia8PtGun",
                         PGunParameters = cms.PSet(
        MaxPt = cms.double(75.01),
        MinPt = cms.double(74.99),
        ParticleID = cms.vint32(-13),
        AddAntiParticle = cms.bool(False),
        MaxEta = cms.double(0.9),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-0.9),
        MinPhi = cms.double(-3.14159265359) ## in radians
        ),
                         Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
                         psethack = cms.string('single mu pt 75 1stbin'),
                         firstRun = cms.untracked.uint32(1),
                         PythiaParameters = cms.PSet(parameterSets = cms.vstring())
                         )
