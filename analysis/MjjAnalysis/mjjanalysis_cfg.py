import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        $inputFileNames
#   'file:/hdfs/store/mc/Fall10/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/START38_V12-v1/0029/08829E4D-C3F8-DF11-BDCB-00238BDFD154.root'
'rfio:/castor/cern.ch/user/l/lovedeep/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-START38_V12-v2_GEN-SIM-RECO_DE3C12ED-01E4-DF11-95EA-0024E8768439.root'
    )
)

process.TFileService = cms.Service("TFileService", 
#        fileName = cms.untracked.string('$outputFileName')
      fileName = cms.string("histo.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.demo = cms.EDAnalyzer('MjjAnalysis'
)


process.p = cms.Path(process.demo)#+process.out)
