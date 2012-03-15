import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     $inputFileNames
#'file:/hdfs/store/user/anil79/newElePF2PAT_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1_Feb25-PF2PAT/1/PF2PAT-C2CCC0F0-8F3D-E111-BC90-E0CB4E553663.root'
)
)
baseJetSel = cms.PSet(
  Jets=cms.InputTag("selectedPatJetsPFlow")
)

from DelPanj.TreeMaker.eSelLvdp2011_cff import *
process.tree = cms.EDAnalyzer(
            'TreeMaker',
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
from DelPanj.TreeMaker.evtCounter_cfi import *

process.diele = zeeFilterLvdp2011NoPtEtaNoIdNoIsoNoMWin
process.zee = zeeFilterLvdp2011NoIdNoIsoNoMWin
process.zee.scoreHistos = cms.double(1)
process.zeeId = zeeFilterLvdp2011NoIsoNoMWin.clone()  
process.zeeIdIso = zeeFilterLvdp2011NoMWin
process.zeeIdIsoMw = zeeFilterLvdp2011
process.zeeIdIsoMw.scoreHistos = cms.double(0)
process.trigfil = trigFilterData.clone()

#print process.zeeId.leadElecPset_

process.countAll = evtCounter.clone()
process.countAll.instance=0
process.countKin = evtCounter.clone()
process.countKin.instance=1
process.countId = evtCounter.clone()
process.countId.instance=2
process.countIdIso = evtCounter.clone()
process.countIdIso.instance=3
process.countIdIsoMw = evtCounter.clone()
process.countIdIsoMw.instance=4
process.countIdIsoMwTrig = evtCounter.clone()
process.countIdIsoMwTrig.instance=5

verbose =0 ;
if verbose:
	print "leadElecPSet_=", process.zeeIdIsoMw.leadElecPset_
	print "subLeadElecPSet_=", process.zeeIdIsoMw.subLeadElecPset_
	print "zMassLowerLimit_=",  process.zeeIdIsoMw.zMassLowerLimit
	print "zMassUpperLimit_=",  process.zeeIdIsoMw.zMassUpperLimit

process.TFileService = cms.Service("TFileService",
      fileName = cms.string('$outputFileName'),
#      fileName = cms.string("zee.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(
        #process.countAll
        #+(process.diele*process.countKin)##for debugging
	process.trigfil*process.zee#*process.countKin) ## no cuts on ele (except kinematic)
        #+(process.zeeId*process.countId)#+##id cuts applied to ele
        #+(process.zeeIdIso*process.countIdIso)#+##id, iso cuts applied.
        #+(process.zeeIdIsoMw*process.countIdIsoMw)##all the cuts applied.
	#*process.trigfil*process.countIdIsoMwTrig*process.tree##Trigger Applied.
	)
  
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
 



