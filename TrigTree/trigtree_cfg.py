import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("START3X_V26A::All") # the one for jan29 MC rereco
#process.GlobalTag.globaltag = cms.string("GR10_P_V5::All") # data above 133332

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#     'file:run135149_1evt.root'
     '/store/relval/CMSSW_3_5_5/RelValZEE/GEN-SIM-RECO/MC_3XY_V25-v1/0006/84C40920-D937-DF11-86C0-001A9281172E.root'
    )
)
# require physics declared
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )


# configure HLT
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')


process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15), 
                                           maxd0 = cms.double(2) 
                                           )


process.demo = cms.EDAnalyzer('TrigTree',
fillElec=cms.untracked.bool(False)
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("TrigTree.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )


process.p = cms.Path(
   # process.hltLevel1GTSeed* #MC->block
    process.scrapingVeto*
   # process.hltPhysicsDeclared* #MC->block
    process.primaryVertexFilter*
    process.demo
    )
#process.p = cms.Path(process.demo)
