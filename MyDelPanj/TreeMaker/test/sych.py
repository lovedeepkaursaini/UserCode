import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     $inputFileNames
#'file:/scratch/anil79/patTuple_PF2PAT_442.root'
#'file:/hdfs/store/user/anil79/PF2PAT_DYJets_Mar4-ColinPF2PAT/1/ColinPF2PAT-F00C5B98-713D-E111-9243-90E6BA0D0990.root'
#'file:/hdfs/store/user/anil79/newElePF2PAT_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1_Feb25-PF2PAT/1/PF2PAT-F00C5B98-713D-E111-9243-90E6BA0D0990.root'
    )
)
from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *
process.myPatElectrons = selectedPatElectrons.clone()
process.myPatElectrons.src = cms.InputTag("selectedPatElectronsPFlow")

process.demo = cms.EDAnalyzer('SychZJet'
)
process.TFileService = cms.Service("TFileService",
                             fileName = cms.string('$outputFileName'),
#                                   fileName = cms.string("SychzeeTest.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )
from DelPanj.TreeMaker.trigFilter_cff import *
process.trigfil = trigFilterData.clone()
process.p = cms.Path(
            process.trigfil*
	    process.demo
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
 
