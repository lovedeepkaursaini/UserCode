import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso

process = cms.Process("p")

MC_flag = False
GLOBAL_TAG = 'START44_V12::All'

HLTProcessName = "HLT"
OUTPUT_FILE_NAME = "TnPTree_DYtoEE.root"


ELECTRON_ET_CUT_MIN = 10.0
ELECTRON_COLL = "gsfElectrons"
ELECTRON_CUTS = "ecalDrivenSeed==1 && (abs(superCluster.eta)<2.5) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
####

PHOTON_COLL = "photons"
PHOTON_CUTS = "hadronicOverEm<0.15 && (abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.566) && ((isEB && sigmaIetaIeta<0.01) || (isEE && sigmaIetaIeta<0.03)) && (superCluster.energy*sin(superCluster.position.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
####
SUPERCLUSTER_COLL_EB = "correctedHybridSuperClusters"
SUPERCLUSTER_COLL_EE = "correctedMulti5x5SuperClustersWithPreshower"
SUPERCLUSTER_CUTS = "abs(eta)<2.5 && et>" + str(ELECTRON_ET_CUT_MIN)


JET_COLL = "ak5PFJets"
JET_CUTS = "abs(eta)<2.4 && chargedHadronEnergyFraction>0 && electronEnergyFraction<0.1 && nConstituents>1 && neutralHadronEnergyFraction<0.99 && neutralEmEnergyFraction<0.99"

#########################################################
#TRIGGER TO BE TESTED                                                 #
#######################################################################
HLTEle17Path1 = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1"
HLTEle17Path2 = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2"
HLTEle17Path3 = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3"
HLTEle17Path4 = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4"
HLTEle17Path5 = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5"
HLTEle17Path6 = "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6"
HLTEle17Path7 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6"
HLTEle17Path8 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7"
HLTEle17Path9 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8"
HLTEle17Path10 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9"
HLTEle17Path11 = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10"
HLTEle8Path1 = HLTEle17Path1
HLTEle8Path2 = HLTEle17Path2
HLTEle8Path3 = HLTEle17Path3
HLTEle8Path4 = HLTEle17Path4
HLTEle8Path5 = HLTEle17Path5
HLTEle8Path6 = HLTEle17Path6
HLTEle8Path7 = HLTEle17Path7
HLTEle8Path8 = HLTEle17Path8
HLTEle8Path9 = HLTEle17Path9
HLTEle8Path10 = HLTEle17Path10
HLTEle8Path11 = HLTEle17Path11
hltEle17Filter  = "hltEle17CaloIdLCaloIsoVLPixelMatchFilter"
hltEle8Filter   = "hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"
hltEle17Filter7 = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter"
hltEle8Filter7  = "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter"

hltTagsPassingEle17HLT= cms.VInputTag(
    cms.InputTag(HLTEle17Path1,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path2,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path3,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path4,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path5,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path6,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path7,hltEle17Filter7, HLTProcessName),
    cms.InputTag(HLTEle17Path8,hltEle17Filter7, HLTProcessName),
    cms.InputTag(HLTEle17Path9,hltEle17Filter7, HLTProcessName),
    cms.InputTag(HLTEle17Path10,hltEle17Filter7, HLTProcessName),
    cms.InputTag(HLTEle17Path11,hltEle17Filter7, HLTProcessName),
     )

hltTagsPassingEle8HLT =cms.VInputTag(
    cms.InputTag(HLTEle8Path1, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path2, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path3, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path4, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path5, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path6, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path7, hltEle8Filter7 , HLTProcessName),
    cms.InputTag(HLTEle8Path8, hltEle8Filter7 , HLTProcessName),
    cms.InputTag(HLTEle8Path9, hltEle8Filter7 , HLTProcessName),
    cms.InputTag(HLTEle8Path10, hltEle8Filter7 , HLTProcessName),
    cms.InputTag(HLTEle8Path11, hltEle8Filter7 , HLTProcessName),
     )
hltTagsPassingNotEle17HLT =cms.VInputTag(
    cms.InputTag(HLTEle17Path1,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path2,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path3,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path4,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path5,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path6,hltEle17Filter,  HLTProcessName),
    cms.InputTag(HLTEle17Path7,hltEle17Filter7, HLTProcessName),
    cms.InputTag(HLTEle17Path8,hltEle17Filter7, HLTProcessName),
    cms.InputTag(HLTEle17Path9,hltEle17Filter7, HLTProcessName),
    cms.InputTag(HLTEle17Path10,hltEle17Filter7, HLTProcessName),
    cms.InputTag(HLTEle17Path11,hltEle17Filter7, HLTProcessName),
     )
hltTagsPassingEle8NotEle17HLT =cms.VInputTag(
    cms.InputTag(HLTEle8Path1, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path2, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path3, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path4, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path5, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path6, hltEle8Filter ,  HLTProcessName),
    cms.InputTag(HLTEle8Path7, hltEle8Filter7 , HLTProcessName),
    cms.InputTag(HLTEle8Path8, hltEle8Filter7 , HLTProcessName),
    cms.InputTag(HLTEle8Path9, hltEle8Filter7 , HLTProcessName),
    cms.InputTag(HLTEle8Path10, hltEle8Filter7 , HLTProcessName),
    cms.InputTag(HLTEle8Path11, hltEle8Filter7 , HLTProcessName),
     )

######################################
#########################################################
#NEUTRAL TRIGGER
#########################################################

HLTEle17Mass30Path1 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1"
HLTEle17Mass30Path2 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2"
HLTEle17Mass30Path3 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3"
HLTEle17Mass30Path4 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v2"
HLTEle17Mass30Path5 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v3"
HLTEle17Mass30Path6 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v4"
HLTEle17Mass30Path7 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v5"
HLTEle17Mass30Path8 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v6"

HLTEle17Mass30Path9 = "HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v8"
HLTEle17Mass30Path10 ="HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v9"
HLTEle17Mass30Path11 ="HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v10"

HLTEle8Mass30Path1 = HLTEle17Mass30Path1
HLTEle8Mass30Path2 = HLTEle17Mass30Path2
HLTEle8Mass30Path3 = HLTEle17Mass30Path3
HLTEle8Mass30Path4 = HLTEle17Mass30Path4
HLTEle8Mass30Path5 = HLTEle17Mass30Path5
HLTEle8Mass30Path6 = HLTEle17Mass30Path6
HLTEle8Mass30Path7 = HLTEle17Mass30Path7
HLTEle8Mass30Path8 = HLTEle17Mass30Path8
HLTEle8Mass30Path9 = HLTEle17Mass30Path9
HLTEle8Mass30Path10 = HLTEle17Mass30Path10
HLTEle8Mass30Path11 = HLTEle17Mass30Path11

nuetHLTEle17Mass30FilterSC8 = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter"
nuetHLTEle8Mass30FilterSC8 =  "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter"

nuetHLTEle17Mass30FilterEle8 = "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter"
nuetHLTEle8Mass30FilterEle8 =  "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter"
hltTagsPassingEle8Mass30HLTSC = cms.VInputTag(
     cms.InputTag(HLTEle8Mass30Path1, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path2, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path3, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path4, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path5, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path6, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path7, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path8, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path9, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path10, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path11, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     )

hltTagsPassingEle8Mass30HLTGsf = cms.VInputTag(
     cms.InputTag(HLTEle8Mass30Path1, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path2, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path3, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path4, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path5, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path6, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path7, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path8, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path9, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path10, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path11, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     )
hltTagsPassingEle17Mass30HLT =cms.VInputTag(
     cms.InputTag(HLTEle17Mass30Path1, nuetHLTEle17Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path2, nuetHLTEle17Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path3, nuetHLTEle17Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path4, nuetHLTEle17Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path5, nuetHLTEle17Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path6, nuetHLTEle17Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path7, nuetHLTEle17Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path8, nuetHLTEle17Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path9, nuetHLTEle17Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path10, nuetHLTEle17Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle17Mass30Path11, nuetHLTEle17Mass30FilterSC8, HLTProcessName),
     )

hltTagsPassingEle8Mass30HLT =cms.VInputTag(
     cms.InputTag(HLTEle8Mass30Path1, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path2, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path3, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path4, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path5, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path6, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path7, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path8, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path9, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path10, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path11, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     )

hltTagsPassingEle8Mass30HLTWPMedium =cms.VInputTag(
     cms.InputTag(HLTEle8Mass30Path1, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path2, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path3, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path4, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path5, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path6, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path7, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path8, nuetHLTEle8Mass30FilterEle8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path9, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path10, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     cms.InputTag(HLTEle8Mass30Path11, nuetHLTEle8Mass30FilterSC8, HLTProcessName),
     )

########################

##    ___            _           _      
##   |_ _|_ __   ___| |_   _  __| | ___ 
##    | || '_ \ / __| | | | |/ _` |/ _ \
##    | || | | | (__| | |_| | (_| |  __/
##   |___|_| |_|\___|_|\__,_|\__,_|\___|
##
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi") 
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = GLOBAL_TAG
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

##  |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
##  | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
##  |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
##  |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
##  
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
       $inputFileNames
#'file:/hdfs/store/mc/Fall11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/AODSIM/PU_S6_START44_V9B-v1/0000/FE941796-9929-E111-BF06-00A0D1EC3950.root'
#'file:/hdfs/store/data/Run2011A/DoubleElectron/AOD/03Oct2011-v1/0000/00108EDF-28EF-E011-8414-0026189438F9.root'
#"file:/hdfs/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/193/336/04FCF820-9D97-E111-BD56-BCAEC5364C42.root",
#"file:/hdfs/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/193/336/062C556B-9C97-E111-A297-001D09F2B2CF.root",
#"file:/hdfs/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/193/336/F02F4FE0-8097-E111-B906-002481E0DEC6.root",
                                                  )
)

#if not MC_flag:
 #       import FWCore.PythonUtilities.LumiList as LumiList
  #      import FWCore.ParameterSet.Types as CfgTypes
   #     myLumis = LumiList.LumiList(filename = '/afs/hep.wisc.edu/home/anil79/TnP/CMSSW_5_2_4_patch4/src/PhysicsTools/TagAndProbe/test/excbadCert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt').getCMSSWString().split(',')
     #   process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    #    process.source.lumisToProcess.extend(myLumis)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    
process.source.inputCommands = cms.untracked.vstring("keep *","drop *_MEtoEDMConverter_*_*")


process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)
process.JetsForRho = cms.Sequence( process.kt6PFJets )

##   ____                         ____ _           _            
##  / ___| _   _ _ __   ___ _ __ / ___| |_   _ ___| |_ ___ _ __ 
##  \___ \| | | | '_ \ / _ \ '__| |   | | | | / __| __/ _ \ '__|
##   ___) | |_| | |_) |  __/ |  | |___| | |_| \__ \ ||  __/ |   
##  |____/ \__,_| .__/ \___|_|   \____|_|\__,_|___/\__\___|_|   
##  

#  SuperClusters  ################
process.superClusters = cms.EDProducer("SuperClusterMerger",
   src = cms.VInputTag(cms.InputTag( SUPERCLUSTER_COLL_EB ,"", "RECO"),
                       cms.InputTag( SUPERCLUSTER_COLL_EE ,"", "RECO") )  
)

process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
   src = cms.InputTag("superClusters"),
   particleType = cms.int32(11),
)

#   Get the above SC's Candidates and place a cut on their Et and eta
process.goodSuperClusters = cms.EDFilter("CandViewSelector",
      src = cms.InputTag("superClusterCands"),
      cut = cms.string( SUPERCLUSTER_CUTS ),
      filter = cms.bool(True)
)                                         
                                         

#### remove real jets (with high hadronic energy fraction) from SC collection
##### this improves the purity of the probe sample without affecting efficiency

process.JetsToRemoveFromSuperCluster = cms.EDFilter("CaloJetSelector",   
    src = cms.InputTag("ak5CaloJets"),
    cut = cms.string('pt>5 && energyFractionHadronic > 0.15')
)
process.goodSuperClustersClean = cms.EDProducer("CandViewCleaner",
    srcObject = cms.InputTag("goodSuperClusters"),
    module_label = cms.string(''),
    srcObjectsToRemove = cms.VInputTag(cms.InputTag("JetsToRemoveFromSuperCluster")),
    deltaRMin = cms.double(0.1)
)

#  Photons!!! ################ 
process.goodPhotons = cms.EDFilter(
    "PhotonSelector",
    src = cms.InputTag( PHOTON_COLL ),
    cut = cms.string(PHOTON_CUTS)
    )


process.PassingEle8Mass30HLTSC = cms.EDProducer("trgMatchedCandidateProducer",    
    InputProducer = cms.InputTag("goodSuperClustersClean" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),                                   
    hltTags = hltTagsPassingEle8Mass30HLTSC,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.sc_sequence = cms.Sequence(
    process.superClusters +
    process.superClusterCands +
    process.goodSuperClusters +
    process.JetsToRemoveFromSuperCluster +
    process.goodSuperClustersClean +
    process.PassingEle8Mass30HLTSC +
    process.goodPhotons 
    )


##    ____      __ _____ _           _                   
##   / ___|___ / _| ____| | ___  ___| |_ _ __ ___  _ __  
##  | |  _/ __| |_|  _| | |/ _ \/ __| __| '__/ _ \| '_ \ 
##  | |_| \__ \  _| |___| |  __/ (__| |_| | | (_) | | | |
##   \____|___/_| |_____|_|\___|\___|\__|_|  \___/|_| |_|
##  
#  GsfElectron ################ 

process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')


process.goodElectrons = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( ELECTRON_COLL ),
    cut = cms.string( ELECTRON_CUTS )    
)

process.PassingEle8Mass30HLTGsf = cms.EDProducer("trgMatchedPatElectronProducer",    
    InputProducer = cms.InputTag("goodElectrons" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle8Mass30HLTGsf,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.GsfMatchedSuperClusterCands = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("PassingEle8Mass30HLTSC"),
   ReferenceElectronCollection = cms.untracked.InputTag("PassingEle8Mass30HLTGsf"),
   deltaR =  cms.untracked.double(0.3)
)


process.GsfMatchedSuperClusterCandsNoTrig = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
   deltaR =  cms.untracked.double(0.3)
)

process.GsfMatchedPhotonCands = process.GsfMatchedSuperClusterCands.clone()
process.GsfMatchedPhotonCands.src = cms.InputTag("goodPhotons")

            

##    _____ _           _                     ___    _ 
##   | ____| | ___  ___| |_ _ __ ___  _ __   |_ _|__| |
##   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \   | |/ _` |
##   | |___| |  __/ (__| |_| | | (_) | | | |  | | (_| |
##   |_____|_|\___|\___|\__|_|  \___/|_| |_| |___\__,_|
##   
# Electron ID  ######

process.PassingWPMedium = cms.EDProducer("ElectronCutBasedCandidateProducer",
   ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
   IsoDepElectron = cms.VInputTag(cms.InputTag('elPFIsoDepositChargedPFIso'),
                                  cms.InputTag('elPFIsoDepositGammaPFIso'),
                                  cms.InputTag('elPFIsoDepositNeutralPFIso')),
   IsoValElectronPF = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
   conversionsInputTag     = cms.InputTag("allConversions"),
   beamSpotInputTag        = cms.InputTag("offlineBeamSpot"),
   rhoIsoInputTag          = cms.InputTag("kt6PFJets", "rho"),
   primaryVertexInputTag   = cms.InputTag("offlinePrimaryVertices"),
   nameIDBOOL = cms.string('medium'),
)
#///meri bari
process.PassingWPTight = cms.EDProducer("ElectronCutBasedCandidateProducer",
   ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
   IsoDepElectron = cms.VInputTag(cms.InputTag('elPFIsoDepositChargedPFIso'),
                                  cms.InputTag('elPFIsoDepositGammaPFIso'),
                                  cms.InputTag('elPFIsoDepositNeutralPFIso')),
   IsoValElectronPF = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
   conversionsInputTag     = cms.InputTag("allConversions"),
   beamSpotInputTag        = cms.InputTag("offlineBeamSpot"),
   rhoIsoInputTag          = cms.InputTag("kt6PFJets", "rho"),
   primaryVertexInputTag   = cms.InputTag("offlinePrimaryVertices"),
   nameIDBOOL = cms.string('tight'),
)

process.PassingWPLoose = cms.EDProducer("ElectronCutBasedCandidateProducer",
   ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
   IsoDepElectron = cms.VInputTag(cms.InputTag('elPFIsoDepositChargedPFIso'),
                                  cms.InputTag('elPFIsoDepositGammaPFIso'),
                                  cms.InputTag('elPFIsoDepositNeutralPFIso')),
   IsoValElectronPF = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
   conversionsInputTag     = cms.InputTag("allConversions"),
   beamSpotInputTag        = cms.InputTag("offlineBeamSpot"),
   rhoIsoInputTag          = cms.InputTag("kt6PFJets", "rho"),
   primaryVertexInputTag   = cms.InputTag("offlinePrimaryVertices"),
   nameIDBOOL = cms.string('loose'),
)
process.PassingWPVeto = cms.EDProducer("ElectronCutBasedCandidateProducer",
   ReferenceElectronCollection = cms.untracked.InputTag("goodElectrons"),
   IsoDepElectron = cms.VInputTag(cms.InputTag('elPFIsoDepositChargedPFIso'),
                                  cms.InputTag('elPFIsoDepositGammaPFIso'),
                                  cms.InputTag('elPFIsoDepositNeutralPFIso')),
   IsoValElectronPF = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                    cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
   conversionsInputTag     = cms.InputTag("allConversions"),
   beamSpotInputTag        = cms.InputTag("offlineBeamSpot"),
   rhoIsoInputTag          = cms.InputTag("kt6PFJets", "rho"),
   primaryVertexInputTag   = cms.InputTag("offlinePrimaryVertices"),
   nameIDBOOL = cms.string('veto'),
)
#///ithe tak 
##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
##   
# Trigger  ##################

process.PassingEle17Mass30HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle17Mass30HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)
#meri bari
process.PassingEle17Mass30HLTTight = cms.EDProducer("trgMatchedGsfElectronProducer",
    InputProducer = cms.InputTag("PassingWPTight" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle17Mass30HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)
)
process.PassingEle17Mass30HLTLoose = cms.EDProducer("trgMatchedGsfElectronProducer",
    InputProducer = cms.InputTag("PassingWPLoose" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle17Mass30HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)
)
process.PassingEle17Mass30HLTVeto = cms.EDProducer("trgMatchedGsfElectronProducer",
    InputProducer = cms.InputTag("PassingWPVeto" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle17Mass30HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)
)
#ithe tak



process.PassingEle17HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    noHltFiring = cms.untracked.bool(True),
    hltTags =hltTagsPassingEle17HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingEle8HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    noHltFiring = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle8HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingNotEle17HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    noHltFiring = cms.untracked.bool(True),
    antiSelect = cms.untracked.bool(True),
    hltTags = hltTagsPassingNotEle17HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingEle8NotEle17HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingNotEle17HLT" ),
    isTriggerFilter = cms.untracked.bool(True),
    noHltFiring = cms.untracked.bool(True),
    hltTags = hltTagsPassingEle8NotEle17HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingEle8Mass30HLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle8Mass30HLT,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

process.PassingEle8Mass30HLTWPMedium = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag("PassingWPMedium" ),
    isTriggerFilter = cms.untracked.bool(True),
    matchUnprescaledTriggerOnly = cms.untracked.bool(False),
    hltTags = hltTagsPassingEle8Mass30HLTWPMedium,
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)



##    _____      _                        _  __     __             
##   | ____|_  _| |_ ___ _ __ _ __   __ _| | \ \   / /_ _ _ __ ___ 
##   |  _| \ \/ / __/ _ \ '__| '_ \ / _` | |  \ \ / / _` | '__/ __|
##   | |___ >  <| ||  __/ |  | | | | (_| | |   \ V / (_| | |  \__ \
##   |_____/_/\_\\__\___|_|  |_| |_|\__,_|_|    \_/ \__,_|_|  |___/
##   
## Here we show how to use a module to compute an external variable
## process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
## ak5PFResidual.useCondDB = False

process.superClusterDRToNearestJet = cms.EDProducer("DeltaRNearestJetComputer",
    probes = cms.InputTag("goodSuperClusters"),
       # ^^--- NOTA BENE: if probes are defined by ref, as in this case, 
       #       this must be the full collection, not the subset by refs.
    objects = cms.InputTag(JET_COLL),
    objectSelection = cms.string(JET_CUTS + " && pt > 20.0"),
)
process.JetMultiplicityInSCEvents = cms.EDProducer("CandMultiplicityCounter",
    probes = cms.InputTag("goodSuperClusters"),
    objects = cms.InputTag(JET_COLL),
    objectSelection = cms.string(JET_CUTS + " && pt > 20.0"),
)

process.PhotonDRToNearestJet = process.superClusterDRToNearestJet.clone()
process.PhotonDRToNearestJet.probes =cms.InputTag("goodPhotons")
process.JetMultiplicityInPhotonEvents = process.JetMultiplicityInSCEvents.clone()
process.JetMultiplicityInPhotonEvents.probes = cms.InputTag("goodPhotons")

process.GsfDRToNearestJet = process.superClusterDRToNearestJet.clone()
process.GsfDRToNearestJet.probes = cms.InputTag( ELECTRON_COLL )
process.JetMultiplicityInGsfEvents = process.JetMultiplicityInSCEvents.clone()
process.JetMultiplicityInGsfEvents.probes = cms.InputTag( ELECTRON_COLL )

process.ext_ToNearestJet_sequence = cms.Sequence(
    #process.ak5PFResidual + 
    process.superClusterDRToNearestJet +
    process.JetMultiplicityInSCEvents +
    process.PhotonDRToNearestJet +
    process.JetMultiplicityInPhotonEvents +    
    process.GsfDRToNearestJet +
    process.JetMultiplicityInGsfEvents
    )


##    _____             ____        __ _       _ _   _             
##   |_   _|_ _  __ _  |  _ \  ___ / _(_)_ __ (_) |_(_) ___  _ __  
##     | |/ _` |/ _` | | | | |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \ 
##     | | (_| | (_| | | |_| |  __/  _| | | | | | |_| | (_) | | | |
##     |_|\__,_|\__, | |____/ \___|_| |_|_| |_|_|\__|_|\___/|_| |_|
##              |___/
## 


process.TagMedium = process.PassingEle17Mass30HLT.clone()
process.TagTight = process.PassingEle17Mass30HLTTight.clone()
process.TagLoose = process.PassingEle17Mass30HLTLoose.clone()
process.TagVeto = process.PassingEle17Mass30HLTVeto.clone()
process.TagMedium.InputProducer = cms.InputTag( "PassingWPMedium" )
process.TagTight.InputProducer = cms.InputTag( "PassingWPTight" )
process.TagLoose.InputProducer = cms.InputTag( "PassingWPLoose" )
process.TagVeto.InputProducer = cms.InputTag( "PassingWPVeto" )

process.TagMatchedSuperClusterCandsClean = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("TagMedium"),
   deltaR =  cms.untracked.double(0.3)
)
process.TagMatchedSuperClusterCandsCleanTight = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("TagTight"),
   deltaR =  cms.untracked.double(0.3)
)
process.TagMatchedSuperClusterCandsCleanLoose = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("TagLoose"),
   deltaR =  cms.untracked.double(0.3)
)
process.TagMatchedSuperClusterCandsCleanVeto = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("TagVeto"),
   deltaR =  cms.untracked.double(0.3)
)

process.TagMatchedPhotonCands = process.TagMatchedSuperClusterCandsClean.clone()
process.TagMatchedPhotonCands.src     = cms.InputTag("goodPhotons")

process.TagMatchedPhotonCandsTight = process.TagMatchedSuperClusterCandsCleanTight.clone()
process.TagMatchedPhotonCandsTight.src     = cms.InputTag("goodPhotons")
process.TagMatchedPhotonCandsLoose = process.TagMatchedSuperClusterCandsCleanLoose.clone()
process.TagMatchedPhotonCandsLoose.src     = cms.InputTag("goodPhotons")
process.TagMatchedPhotonCandsVeto = process.TagMatchedSuperClusterCandsCleanVeto.clone()
process.TagMatchedPhotonCandsVeto.src     = cms.InputTag("goodPhotons")

process.ele_sequence = cms.Sequence(
    process.pfParticleSelectionSequence + 
    process.eleIsoSequence + 
    process.goodElectrons +
    process.PassingEle8Mass30HLTGsf +
    process.GsfMatchedSuperClusterCands +
    process.GsfMatchedSuperClusterCandsNoTrig +
    process.GsfMatchedPhotonCands +
    process.kt6PFJets *
    process.PassingWPMedium +
    process.PassingWPTight +
    process.PassingWPLoose +
    process.PassingWPVeto +
    process.PassingEle17Mass30HLT +
    process.PassingEle17Mass30HLTTight +
    process.PassingEle17Mass30HLTLoose +
    process.PassingEle17Mass30HLTVeto +
    process.PassingEle17HLT +
    process.PassingEle8HLT +
    process.PassingNotEle17HLT +
    process.PassingEle8NotEle17HLT +
    process.PassingEle8Mass30HLT +
    process.PassingEle8Mass30HLTWPMedium +
    process.TagMedium +
    process.TagMatchedSuperClusterCandsClean +
    process.TagMatchedPhotonCands+
    process.TagTight +
    process.TagMatchedSuperClusterCandsCleanTight +
    process.TagMatchedPhotonCandsTight+
    process.TagLoose +
    process.TagMatchedSuperClusterCandsCleanLoose +
    process.TagMatchedPhotonCandsLoose+
    process.TagVeto +
    process.TagMatchedSuperClusterCandsCleanVeto +
    process.TagMatchedPhotonCandsVeto

    )


##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
#  Tag & probe selection ######
process.tagSC = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingEle17Mass30HLT goodSuperClustersClean"), # charge conjugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("40 < mass < 10000"),
)

process.tagPhoton = process.tagSC.clone()
process.tagPhoton.decay = cms.string("PassingEle17Mass30HLT goodPhotons")
process.GsfGsf = process.tagSC.clone()
process.GsfGsf.decay = cms.string("goodElectrons goodElectrons")
process.tagGsf = process.tagSC.clone()
process.tagGsf.decay = cms.string("PassingEle17Mass30HLT goodElectrons")
process.tagWPMedium = process.tagSC.clone()
process.tagWPMedium.decay = cms.string("PassingEle17Mass30HLT PassingWPMedium")
process.Zsignal = process.tagSC.clone()
process.Zsignal.decay = cms.string("PassingEle17HLT PassingEle8HLT")
process.elecMet = process.tagSC.clone()
process.elecMet.decay = cms.string("pfMet PassingWPMedium")
process.elecMet.cut = cms.string("mt > 0")

process.CSVarsTagGsf = cms.EDProducer("ColinsSoperVariablesComputer",
    parentBoson = cms.InputTag("tagGsf")
)
process.CSVarsGsfGsf = process.CSVarsTagGsf.clone()
process.CSVarsGsfGsf.parentBoson = cms.InputTag("GsfGsf")

process.tagTightSC = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingEle17Mass30HLTTight goodSuperClustersClean"), # charge conjugate states are implied
    checkCharge = cms.bool(False),
    cut   = cms.string("40 < mass < 10000"),
)
process.tagTightGsf = process.tagTightSC.clone()
process.tagTightGsf.decay = cms.string("PassingEle17Mass30HLTTight goodElectrons")
process.tagWPTight = process.tagTightSC.clone()
process.tagWPTight.decay = cms.string("PassingEle17Mass30HLTTight PassingWPTight")

process.tagLooseSC = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingEle17Mass30HLTLoose goodSuperClustersClean"), # charge conjugate states are implied
    checkCharge = cms.bool(False),
    cut   = cms.string("40 < mass < 10000"),
)
process.tagLooseGsf = process.tagLooseSC.clone()
process.tagLooseGsf.decay = cms.string("PassingEle17Mass30HLTLoose goodElectrons")
process.tagWPLoose = process.tagLooseSC.clone()
process.tagWPLoose.decay = cms.string("PassingEle17Mass30HLTLoose PassingWPLoose")

process.tagVetoSC = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingEle17Mass30HLTVeto goodSuperClustersClean"), # charge conjugate states are implied
    checkCharge = cms.bool(False),
    cut   = cms.string("40 < mass < 10000"),
)
process.tagVetoGsf = process.tagVetoSC.clone()
process.tagVetoGsf.decay = cms.string("PassingEle17Mass30HLTVeto goodElectrons")
process.tagWPVeto = process.tagVetoSC.clone()
process.tagWPVeto.decay = cms.string("PassingEle17Mass30HLTVeto PassingWPVeto")


process.allTagsAndProbes = cms.Sequence(
    process.Zsignal +
    process.tagSC +
    process.tagPhoton +
    process.tagGsf +
    process.GsfGsf +
    process.tagWPMedium +
    process.elecMet + 
    process.CSVarsTagGsf +
    process.CSVarsGsfGsf+
    process.tagTightSC +
    process.tagTightGsf +
    process.tagWPTight+
    process.tagLooseSC +
    process.tagLooseGsf +
    process.tagWPLoose+
    process.tagVetoSC +
    process.tagVetoGsf +
    process.tagWPVeto

)

##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        
process.McMatchTag = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("TagTight"),
    distMin = cms.double(0.2),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)
process.McMatchSC = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("goodSuperClustersClean"),
    distMin = cms.double(0.2),
    matched = cms.InputTag("genParticles")
)
process.McMatchPhoton = process.McMatchSC.clone()
process.McMatchPhoton.src = cms.InputTag("goodPhotons")
process.McMatchGsf = process.McMatchTag.clone()
process.McMatchGsf.src = cms.InputTag("goodElectrons")
process.McMatchWPMedium = process.McMatchTag.clone()
process.McMatchWPMedium.src = cms.InputTag("PassingWPMedium")
    
process.mc_sequence = cms.Sequence(
   process.McMatchTag +
   process.McMatchSC +
   process.McMatchPhoton +
   process.McMatchGsf + 
   process.McMatchWPMedium 
)

############################################################################
##    _____           _       _ ____            _            _   _  ____  ##
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |/ ___| ##
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | |  _  ##
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | |_| | ##
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_|\____| ##
##              |___/                                                     ##
##                                                                        ##
############################################################################
##    ____                      _     _           
##   |  _ \ ___ _   _ ___  __ _| |__ | | ___  ___ 
##   | |_) / _ \ | | / __|/ _` | '_ \| |/ _ \/ __|
##   |  _ <  __/ |_| \__ \ (_| | |_) | |  __/\__ \
##   |_| \_\___|\__,_|___/\__,_|_.__/|_|\___||___/
##
## I define some common variables for re-use later.
## This will save us repeating the same code for each efficiency category
ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt  = cms.string("pt"),
    mass  = cms.string("mass"),
)   

ProbeVariablesToStore = cms.PSet(
    probe_gsfEle_eta = cms.string("eta"),
    probe_gsfEle_abseta = cms.string("abs(eta)"),
    probe_gsfEle_pt  = cms.string("pt"),
    probe_gsfEle_et  = cms.string("et"),
    probe_gsfEle_e  = cms.string("energy"),
  ## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    probe_sc_eta    = cms.string("superCluster.eta"),
    probe_sc_abseta    = cms.string("abs(superCluster.eta)"),
)


TagVariablesToStore = cms.PSet(
    gsfEle_eta = cms.string("eta"),
    gsfEle_abseta = cms.string("abs(eta)"),
    gsfEle_pt  = cms.string("pt"),
    gsfEle_et  = cms.string("et"),
    gsfEle_e  = cms.string("energy"),
    ## super cluster quantities
    sc_energy = cms.string("superCluster.energy"),
    sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    sc_eta    = cms.string("superCluster.eta"),
    sc_abseta    = cms.string("abs(superCluster.eta)"),
)

CommonStuffForGsfElectronProbe = cms.PSet(
    variables = cms.PSet(ProbeVariablesToStore),
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool (True),
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
          mass60to120 = cms.string("60 < mass < 120")
    ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags     =  cms.PSet(
          passingGsf = cms.InputTag("goodElectrons"),
          isWPMedium = cms.InputTag("PassingWPMedium"),
          passingHLT = cms.InputTag("PassingEle17Mass30HLT"),     
          isWPTight = cms.InputTag("PassingWPTight"),
          passingHLTTight = cms.InputTag("PassingEle17Mass30HLTTight"),
          isWPLoose = cms.InputTag("PassingWPLoose"),
          passingHLTLoose = cms.InputTag("PassingEle17Mass30HLTLoose"),
          isWPVeto = cms.InputTag("PassingWPVeto"),
          passingHLTVeto = cms.InputTag("PassingEle17Mass30HLTVeto"),
    ),    
)

CommonStuffForSuperClusterProbe = CommonStuffForGsfElectronProbe.clone()
CommonStuffForSuperClusterProbe.variables = cms.PSet(
    probe_eta = cms.string("eta"),
    probe_abseta = cms.string("abs(eta)"),
    probe_pt  = cms.string("pt"),
    probe_et  = cms.string("et"),
    probe_e  = cms.string("energy"),
    )


if MC_flag:
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(MC_flag),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_flag),
        checkMotherInUnbiasEff = cms.bool(MC_flag),
        mcVariables = cms.PSet(
        probe_eta = cms.string("eta"),
        probe_abseta = cms.string("abs(eta)"),
        probe_pt  = cms.string("pt"),
        probe_et  = cms.string("et"),
        probe_e  = cms.string("energy"),
        probe_mass  = cms.string("mass"),
        ),
        mcFlags     =  cms.PSet(
        probe_flag = cms.string("pt>0")
        ),      
        )
else:
     mcTruthCommonStuff = cms.PSet(
         isMC = cms.bool(False)
         )


##    ____   ____       __     ____      __ 
##   / ___| / ___|      \ \   / ___|___ / _|
##   \___ \| |      _____\ \ | |  _/ __| |_ 
##    ___) | |___  |_____/ / | |_| \__ \  _|
##   |____/ \____|      /_/   \____|___/_|  
##
## super cluster --> gsf electron
process.SuperClusterToGsfElectron = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration                 
    tagProbePairs = cms.InputTag("tagSC"),
    arbitration   = cms.string("Random2"),                      
    flags = cms.PSet(
        probe_passingGsf = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig"),        
        probe_passingHLT = cms.InputTag("TagMatchedSuperClusterCandsClean")
    ),
    probeMatches  = cms.InputTag("McMatchSC"),
    allProbes     = cms.InputTag("goodSuperClustersClean")
)
process.SuperClusterToGsfElectron.variables.probe_dRjet = cms.InputTag("superClusterDRToNearestJet")
process.SuperClusterToGsfElectron.variables.probe_nJets = cms.InputTag("JetMultiplicityInSCEvents")
#process.SuperClusterToGsfElectron.tagVariables.dRjet = cms.InputTag("GsfDRToNearestJet")
#process.SuperClusterToGsfElectron.tagVariables.nJets = cms.InputTag("JetMultiplicityInGsfEvents")


process.NonIsoTagSuperClusterToGsfElectron = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration                 
    tagProbePairs = cms.InputTag("nonisotagSC"),
    arbitration   = cms.string("Random2"),                      
    flags = cms.PSet(
        probe_passingGsf = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig"),        
        probe_passingHLT = cms.InputTag("TagMatchedSuperClusterCandsClean")
    ),
    probeMatches  = cms.InputTag("McMatchSC"),
    allProbes     = cms.InputTag("goodSuperClustersClean")
)
process.SuperClusterToGsfElectron.variables.probe_dRjet = cms.InputTag("superClusterDRToNearestJet")
process.SuperClusterToGsfElectron.variables.probe_nJets = cms.InputTag("JetMultiplicityInSCEvents")
#process.SuperClusterToGsfElectron.tagVariables.dRjet = cms.InputTag("GsfDRToNearestJet")
#process.SuperClusterToGsfElectron.tagVariables.nJets = cms.InputTag("JetMultiplicityInGsfEvents")



## good photon --> gsf electron
process.PhotonToGsfElectron = process.SuperClusterToGsfElectron.clone()
process.PhotonToGsfElectron.tagProbePairs = cms.InputTag("tagPhoton")
process.PhotonToGsfElectron.flags = cms.PSet(
    probe_passingGsf = cms.InputTag("GsfMatchedPhotonCands"),
    probe_passingHLT = cms.InputTag("TagMatchedPhotonCands"),
    )
process.PhotonToGsfElectron.probeMatches  = cms.InputTag("McMatchPhoton")
process.PhotonToGsfElectron.allProbes     = cms.InputTag("goodPhotons")
process.PhotonToGsfElectron.variables.probe_dRjet = cms.InputTag("PhotonDRToNearestJet")
process.PhotonToGsfElectron.variables.probe_nJets = cms.InputTag("JetMultiplicityInPhotonEvents")
process.PhotonToGsfElectron.variables.probe_trackiso = cms.string("trkSumPtHollowConeDR03")
process.PhotonToGsfElectron.variables.probe_ecaliso = cms.string("ecalRecHitSumEtConeDR03")
process.PhotonToGsfElectron.variables.probe_hcaliso = cms.string("hcalTowerSumEtConeDR03")
process.PhotonToGsfElectron.variables.probe_HoverE  = cms.string("hadronicOverEm")
process.PhotonToGsfElectron.variables.probe_sigmaIetaIeta = cms.string("sigmaIetaIeta")
process.PhotonToGsfElectron.tagVariables.dRjet = cms.InputTag("GsfDRToNearestJet")
process.PhotonToGsfElectron.tagVariables.nJets = cms.InputTag("JetMultiplicityInGsfEvents")


##   ____      __       __    ___                 ___    _ 
##  / ___|___ / _|      \ \  |_ _|___  ___       |_ _|__| |
## | |  _/ __| |_   _____\ \  | |/ __|/ _ \       | |/ _` |
## | |_| \__ \  _| |_____/ /  | |\__ \ (_) |  _   | | (_| |
##  \____|___/_|        /_/  |___|___/\___/  ( ) |___\__,_|
##                                           |/            
##  gsf electron --> isolation, electron id  etc.
############# Needed for pileup re-weighting ##########
process.pileupReweightingProducer = cms.EDProducer("PileupWeightProducer",
                                                   FirstTime = cms.untracked.bool(True)
)

process.GsfElectronToId = cms.EDAnalyzer("TagProbeFitTreeProducer",
#    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
    tagProbePairs = cms.InputTag("tagTightGsf"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
#        probe_isWP95 = cms.InputTag("WP95MatchedSuperClusterCandsClean"),
        probe_isWPMedium = cms.InputTag("PassingWPMedium"),
    ),
    probeMatches  = cms.InputTag("McMatchGsf"),
#    allProbes     = cms.InputTag("PassingEle8Mass30HLTGsf")
    allProbes     = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig")
)
process.GsfElectronToId.PUWeightSrc = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.GsfElectronToIdTight = cms.EDAnalyzer("TagProbeFitTreeProducer",
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
    tagProbePairs = cms.InputTag("tagTightGsf"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
        probe_isWPTight = cms.InputTag("PassingWPTight"),
    ),
    probeMatches  = cms.InputTag("McMatchGsf"),
    allProbes     = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig")
)
process.GsfElectronToIdTight.PUWeightSrc = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.GsfElectronToIdLoose = cms.EDAnalyzer("TagProbeFitTreeProducer",
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
    tagProbePairs = cms.InputTag("tagTightGsf"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
        probe_isWPLoose = cms.InputTag("PassingWPLoose"),
    ),
    probeMatches  = cms.InputTag("McMatchGsf"),
    allProbes     = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig")
)
process.GsfElectronToIdLoose.PUWeightSrc = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.GsfElectronToIdVeto = cms.EDAnalyzer("TagProbeFitTreeProducer",
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
    tagProbePairs = cms.InputTag("tagTightGsf"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
        probe_isWPVeto = cms.InputTag("PassingWPVeto"),
    ),
    probeMatches  = cms.InputTag("McMatchGsf"),
    allProbes     = cms.InputTag("GsfMatchedSuperClusterCandsNoTrig")
)
process.GsfElectronToIdVeto.PUWeightSrc = cms.InputTag("pileupReweightingProducer","pileupWeights")

process.Zsignalcounter = cms.EDAnalyzer("TagProbeFitTreeProducer",
#    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
    tagProbePairs = cms.InputTag("Zsignal"),
    arbitration   = cms.string("OnePair"),
    flags = cms.PSet(
#        probe_isWP95 = cms.InputTag("WP95MatchedSuperClusterCandsClean"),
        probe_isWPMedium = cms.InputTag("PassingWPMedium"),
    ),
    probeMatches  = cms.InputTag("McMatchGsf"),
#    allProbes     = cms.InputTag("PassingEle8Mass30HLTGsf")
    allProbes     = cms.InputTag("PassingWPMedium")
)


#process.GsfElectronToId.variables.probe_dRjet = cms.InputTag("GsfDRToNearestJet")
#process.GsfElectronToId.variables.probe_nJets = cms.InputTag("JetMultiplicityInGsfEvents")
#process.GsfElectronToId.tagVariables.dRjet = cms.InputTag("GsfDRToNearestJet")
#process.GsfElectronToId.tagVariables.nJets = cms.InputTag("JetMultiplicityInGsfEvents")
#process.GsfElectronToId.pairVariables.costheta = cms.InputTag("CSVarsTagGsf","costheta")
#process.GsfElectronToId.pairVariables.sin2theta = cms.InputTag("CSVarsTagGsf","sin2theta")
#process.GsfElectronToId.pairVariables.tanphi = cms.InputTag("CSVarsTagGsf","tanphi")


process.GsfElectronPlusGsfElectron = process.GsfElectronToId.clone()
process.GsfElectronPlusGsfElectron.tagProbePairs = cms.InputTag("GsfGsf")
process.GsfElectronPlusGsfElectron.tagMatches = cms.InputTag("McMatchGsf")
process.GsfElectronPlusGsfElectron.pairVariables.costheta = cms.InputTag("CSVarsGsfGsf","costheta")
process.GsfElectronPlusGsfElectron.pairVariables.sin2theta = cms.InputTag("CSVarsGsfGsf","sin2theta")
process.GsfElectronPlusGsfElectron.pairVariables.tanphi = cms.InputTag("CSVarsGsfGsf","tanphi")


process.GsfElectronPlusMet = process.GsfElectronToId.clone()
process.GsfElectronPlusMet.tagProbePairs = cms.InputTag("elecMet")
process.GsfElectronPlusMet.tagVariables = cms.PSet()
process.GsfElectronPlusMet.pairVariables =  cms.PSet(ZVariablesToStore)
process.GsfElectronPlusMet.pairFlags =  cms.PSet( isMTabove40 = cms.string("mt > 40") )
process.GsfElectronPlusMet.isMC = cms.bool(False)


##    ___    _       __    _   _ _   _____ 
##   |_ _|__| |      \ \  | | | | | |_   _|
##    | |/ _` |  _____\ \ | |_| | |   | |  
##    | | (_| | |_____/ / |  _  | |___| |  
##   |___\__,_|      /_/  |_| |_|_____|_|
##
##  offline selection --> HLT. First specify which quantities to store in the TP tree. 
if MC_flag:
    HLTmcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(MC_flag),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_flag),
        checkMotherInUnbiasEff = cms.bool(MC_flag),
        mcVariables = cms.PSet(
          probe_eta = cms.string("eta"),
          probe_abseta = cms.string("abs(eta)"),
          probe_phi  = cms.string("phi"),
          probe_et  = cms.string("et"),
          probe_charge = cms.string("charge"),
        ),
        mcFlags     =  cms.PSet(
          probe_flag = cms.string("pt>0")
        ),      
        )
else:
     HLTmcTruthCommonStuff = cms.PSet(
         isMC = cms.bool(False)
         )

##  WPMedium --> HLTEle17
process.WPMediumToHLTEle17 = cms.EDAnalyzer("TagProbeFitTreeProducer",
    HLTmcTruthCommonStuff,                                
    variables = cms.PSet(
      probe_gsfEle_eta = cms.string("eta"),
      probe_gsfEle_abseta = cms.string("abs(eta)"),
      probe_gsfEle_et  = cms.string("et"),
      probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
      probe_sc_eta    = cms.string("superCluster.eta"), 
      probe_sc_abseta    = cms.string("abs(superCluster.eta)"), 
    ),
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool (True),                                                        
    tagProbePairs = cms.InputTag("tagWPMedium"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet( 
        probe_passingHLT = cms.InputTag("PassingEle17HLT")        
    ),
    probeMatches  = cms.InputTag("McMatchWPMedium"),
    allProbes     = cms.InputTag("PassingWPMedium")
)

##  WPMedium --> HLTEle8NotEle17
process.WPMediumToHLTEle8NotEle17 = process.WPMediumToHLTEle17.clone()
process.WPMediumToHLTEle8NotEle17.tagProbePairs = cms.InputTag("tagWPMedium")
process.WPMediumToHLTEle8NotEle17.probeMatches  = cms.InputTag("McMatchWPMedium")
process.WPMediumToHLTEle8NotEle17.allProbes     = cms.InputTag("PassingWPMedium")
process.WPMediumToHLTEle8NotEle17.flags =cms.PSet( 
        probe_passingHLT = cms.InputTag("PassingEle8NotEle17HLT")        
    ) 

process.tree_sequence = cms.Sequence(
    process.SuperClusterToGsfElectron +
    process.GsfElectronToId +
    process.GsfElectronToIdTight +
    process.GsfElectronToIdLoose +
    process.GsfElectronToIdVeto +
    process.Zsignalcounter +
    process.WPMediumToHLTEle17 +
    process.WPMediumToHLTEle8NotEle17 
)    

##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##
process.out = cms.OutputModule("PoolOutputModule", 
       fileName = cms.untracked.string("PFEE.root"),
       SelectEvents = cms.untracked.PSet( 
       SelectEvents = cms.vstring("p")
       )
    )
process.outpath = cms.EndPath(process.out)
process.outpath.remove(process.out)

if MC_flag:
    process.p = cms.Path(
        process.JetsForRho + process.sc_sequence + process.eIDSequence + process.ele_sequence + 
        process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.mc_sequence + 
        process.tree_sequence
        )
else:
    process.p = cms.Path(
        process.sc_sequence +
        process.ele_sequence + 
        process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.pileupReweightingProducer +
        process.mc_sequence +
        process.tree_sequence
        )
    
process.TFileService = cms.Service(
#    "TFileService", fileName = cms.string( OUTPUT_FILE_NAME )
    "TFileService", fileName = cms.string( '$outputFileName') 
)
