import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'rfio:/castor/cern.ch/user/l/lovedeep/Fall10_DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia.root'
#        'file:store-mc-Summer10-Wenu-GEN-SIM-RECO-START36_V9_S09-v1-0045-C8491FBA-0D7B-DF11-8B6E-E41F13180DC8.root'
    )
)

process.demo = cms.EDAnalyzer('Te'
)


process.p = cms.Path(process.demo)
