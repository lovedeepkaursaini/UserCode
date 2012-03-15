import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     $inputFileNames
#'file:/hdfs/store/user/anil79/newElePF2PAT_DoubleElectron_Run2011A-08Nov2011-v1_Feb25-PF2PAT/1/PF2PAT-FE18939A-B71B-E111-939B-001A928116CE.root'
#'file:/hdfs/store/user/anil79/PF2PAT_DYJets_Mar4-ColinPF2PAT/1/ColinPF2PAT-F00C5B98-713D-E111-9243-90E6BA0D0990.root'
#'file:/hdfs/store/user/anil79/newElePF2PAT_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1_Feb25-PF2PAT/1/PF2PAT-C2CCC0F0-8F3D-E111-BC90-E0CB4E553663.root'
)
)
baseJetSel = cms.PSet(
  Jets=cms.InputTag("selectedPatJetsPFlow")
)

from DelPanj.TreeMaker.eSelLvdp2011_cff import *
process.tree = cms.EDAnalyzer(
            'TreeMaker',
                    fillPUweightInfo_ = cms.bool(True),
#                    DontDoPUReweight_ = cms.bool(False),#MC
                    DontDoPUReweight_ = cms.bool(True),#Data
                    fillEventInfo_ = cms.bool(True),
                    fillGenInfo_   = cms.bool(False),
                    fillMuonInfo_  = cms.bool(False),
                    fillElecInfo_  = cms.bool(True),
                    fillElecIsoInfo_ = cms.bool(False),
                    fillJetInfo_   = cms.bool(True),
                    fillMetInfo_   = cms.bool(False),
                    fillTrigInfo_  = cms.bool(True),
                    fillPhotInfo_  = cms.bool(False),
                    fillZJetPlant_ = cms.bool(False),
                    genPartLabel=cms.InputTag("genParticles"),
                    patMuons=cms.InputTag("selectedPatMuonsPFlow"),
                    patElectrons = cms.InputTag("selectedPatElectronsPFlow"),
                    leadElecPset_ = eSelLvdp2011,
                    subLeadElecPset_ = eSelLvdp2011,
                    patJetLabel =cms.InputTag("selectedPatJetsPFlow"),
                    patMet=cms.InputTag("patMETs"),
                    beamSpotLabel=cms.InputTag("offlineBeamSpot"),
                    patJetPfAk05 = baseJetSel,
                    outFileName=cms.string('outputFileName.root')
                    )


from DelPanj.TreeMaker.zeeFilter_cff import *
from DelPanj.TreeMaker.trigFilter_cff import *

process.zee = zeeFilterLvdp2011NoIdNoIsoNoMWin
process.zee.scoreHistos = cms.double(1)
process.zeeIdIsoMw = zeeFilterLvdp2011NoMWin
process.trigfil = trigFilterData.clone()

process.TFileService = cms.Service("TFileService",
      fileName = cms.string('$outputFileName'),
#      fileName = cms.string("zee.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(
        #process.zee
	process.trigfil*process.zee*
        process.zeeIdIsoMw*process.tree##Trigger Applied.
	)
  
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
 



