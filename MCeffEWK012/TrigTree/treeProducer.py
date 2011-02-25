import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("MC_38Y_V14::All")  
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/relval/CMSSW_3_8_7/RelValWE/GEN-SIM-RECO/START38_V13-v1/0016/2079B40C-7CFC-DF11-9C66-0018F3D096E6.root',
'/store/relval/CMSSW_3_8_7/RelValWE/GEN-SIM-RECO/START38_V13-v1/0016/2E02EB01-82FC-DF11-AAD3-001731EF61B4.root',
'/store/relval/CMSSW_3_8_7/RelValWE/GEN-SIM-RECO/START38_V13-v1/0016/4A0E1782-80FC-DF11-9DAB-0018F3D0961E.root',
'/store/relval/CMSSW_3_8_7/RelValWE/GEN-SIM-RECO/START38_V13-v1/0016/8228D48C-7AFC-DF11-9A09-0030486791BA.root',
'/store/relval/CMSSW_3_8_7/RelValWE/GEN-SIM-RECO/START38_V13-v1/0017/32912E8D-93FC-DF11-8ACF-002618943915.root'
#        'rfio:/castor/cern.ch/user/l/lovedeep/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-START38_V12-v2_GEN-SIM-RECO_DE3C12ED-01E4-DF11-95EA-0024E8768439.root'
#'/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/802/F6FB8B26-8563-DF11-8488-000423D98950.root',
    )
)



process.demo = cms.EDAnalyzer('TrigTree',   genPartLabel=cms.InputTag("genParticles"),  IsData = cms.untracked.bool(False) )

process.TFileService = cms.Service("TFileService", 
                                  fileName = cms.string("MCtree.root"),
                                 closeFileFast = cms.untracked.bool(True)
                                )


process.p = cms.Path(
    process.demo
    )
