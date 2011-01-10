import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

##    ___            _           _      
##   |_ _|_ __   ___| |_   _  __| | ___ 
##    | || '_ \ / __| | | | |/ _` |/ _ \
##    | || | | | (__| | |_| | (_| |  __/
##   |___|_| |_|\___|_|\__,_|\__,_|\___|

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START38_V14::All')#FT_R_38X_V14A::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load("RecoEgamma.EgammaTools.correctedElectronsProducer_cfi")

########################
isData = False
MC_flag = False
##runs: 138564-140401
#HLTPath = "HLT_Photon15_Cleaned_L1R"

##runs: 141956-144114
#HLTPath = "HLT_Ele15_SW_CaloEleId_L1R"

##runs: 146428-147116
HLTPath = "HLT_Ele17_SW_CaloEleId_L1R"

##runs: 147196-148058
##HLTPath = "HLT_Ele17_SW_TightEleId_L1R"

##runs: 148819-149064
#HLTPath = "HLT_Ele17_SW_TighterEleIdIsol_L1R_v2"

##runs: 149181-149442
#HLTPath = "HLT_Ele17_SW_TighterEleIdIsol_L1R_v3"
########################

OUTPUT_FILE_NAME = "testNewWrite.root"



JET_COLL = "ak5PFJets"
JET_CUTS = "pt>30 && abs(eta)<2.4 && chargedHadronEnergyFraction>0 && chargedEmEnergyFraction<0.99 && nConstituents>1 && neutralHadronEnergyFraction<0.99 && neutralEmEnergyFraction<0.99"

##   ____             _ ____                           
##  |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
##  | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
##  |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
##  |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
##  

readFiles = cms.untracked.vstring()

process.source = cms.Source("PoolSource", fileNames=readFiles)
readFiles.extend( [
       'rfio:/castor/cern.ch/user/l/lovedeep/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-START38_V12-v2_GEN-SIM-RECO_DE3C12ED-01E4-DF11-95EA-0024E8768439.root'
#        '/store/relval/CMSSW_3_8_2/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0019/726DE59F-27B0-DF11-AA89-001BFCDBD1BA.root',
 #       '/store/relval/CMSSW_3_8_2/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0019/06C3C534-BFAF-DF11-BE7A-0030486790BA.root',
  #      '/store/relval/CMSSW_3_8_2/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0018/F4080025-93AF-DF11-BE81-00261894390C.root',
   #     '/store/relval/CMSSW_3_8_2/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0018/DE22B8A8-94AF-DF11-B1AA-003048678A80.root',
    #    '/store/relval/CMSSW_3_8_2/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0018/80F4CF0A-A2AF-DF11-A5F7-001A92810A92.root',
     #   '/store/relval/CMSSW_3_8_2/RelValZEE/GEN-SIM-RECO/START38_V9-v1/0018/3820A6AB-92AF-DF11-834F-0018F3D09710.root'
        #'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/715/24043F7E-19CA-DF11-A0D5-003048F1BF66.root',

]);


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source.inputCommands = cms.untracked.vstring("keep *","drop *_MEtoEDMConverter_*_*")



#  Photons!!! ################ 
process.goodPhotons = cms.EDFilter("PhotonSelector",
                                   src = cms.InputTag("photons"),
                                   cut = cms.string("hadronicOverEm<0.15"
                                                    " && (superCluster.energy*sin(superCluster.position.theta)>20.0)"
                                                    " && (abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.566)")
                                   )


process.FilteredPhotons = cms.EDFilter("PhotonRefSelector",
         src = cms.InputTag("goodPhotons"),
         cut = cms.string(process.goodPhotons.cut.value() +
                          " && ( (isEB && (sigmaIetaIeta<0.01) && (hadronicOverEm<0.15))"
                          " || (isEE && (sigmaIetaIeta<0.03) && (hadronicOverEm<0.15)))"
                          )
)




##   ____                         ____ _           _            
##  / ___| _   _ _ __   ___ _ __ / ___| |_   _ ___| |_ ___ _ __ 
##  \___ \| | | | '_ \ / _ \ '__| |   | | | | / __| __/ _ \ '__|
##   ___) | |_| | |_) |  __/ |  | |___| | |_| \__ \ ||  __/ |   
##  |____/ \__,_| .__/ \___|_|   \____|_|\__,_|___/\__\___|_|   
##  

#  SuperClusters  ################
process.superClusters = cms.EDProducer("SuperClusterMerger",
   src = cms.VInputTag(cms.InputTag("hybridSuperClusters","", "RECO"),
                       cms.InputTag("multi5x5SuperClustersWithPreshower","", "RECO"))  
)

process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
   src = cms.InputTag("superClusters"),
   particleType = cms.int32(11),
)

#   Get the above SC's Candidates and place a cut on their Et and eta
process.goodSuperClusters = cms.EDFilter("CandViewSelector",
      src = cms.InputTag("superClusterCands"),
      cut = cms.string("et>20.0 && abs(eta)<2.5 && !(1.4442< abs(eta) <1.566)"),
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



## process.superClusters = cms.EDFilter("EgammaHLTRecoEcalCandidateProducers",
##    scHybridBarrelProducer =  cms.InputTag("hybridSuperClusters","", "RECO"),
##    scIslandEndcapProducer =  cms.InputTag("multi5x5SuperClustersWithPreshower","", "RECO"),    
##    recoEcalCandidateCollection = cms.string("")
## )


process.sc_sequence = cms.Sequence( process.goodPhotons *
                                    process.FilteredPhotons *
                                    process.superClusters *
                                    process.superClusterCands *
                                    process.goodSuperClusters *
                                    process.JetsToRemoveFromSuperCluster *
                                    process.goodSuperClustersClean
                                    )


##    ____      __ _____ _           _                   
##   / ___|___ / _| ____| | ___  ___| |_ _ __ ___  _ __  
##  | |  _/ __| |_|  _| | |/ _ \/ __| __| '__/ _ \| '_ \ 
##  | |_| \__ \  _| |___| |  __/ (__| |_| | | (_) | | | |
##   \____|___/_| |_____|_|\___|\___|\__|_|  \___/|_| |_|
##  

#  GsfElectron ################



process.GsfConvRej= cms.EDProducer("ConRej",
   probes = cms.InputTag("gsfElectrons"),
   IsData = cms.untracked.bool(isData)
)


process.PassingGsf = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("GsfConvRej"),
    cut = cms.string("(abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.566)"
                     " && ecalDrivenSeed==1 && (gsfTrack.trackerExpectedHitsInner.numberOfHits <= 0)"# && abs(dist)>=0.02 && abs(dcot)>=0.02"
                     " && (ecalEnergy*sin(superClusterPosition.theta)>20.0)")
# && (hadronicOverEm<0.15)")    
)


process.GsfMatchedSuperClusterCands = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("PassingGsf"),
   deltaR =  cms.untracked.double(0.3)
)


process.GsfMatchedPhotonCands = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("FilteredPhotons"),
   ReferenceElectronCollection = cms.untracked.InputTag("PassingGsf"),
   deltaR =  cms.untracked.double(0.3)
)

      


##     ___           _       _   _             
##    |_ _|___  ___ | | __ _| |_(_) ___  _ __  
##     | |/ __|/ _ \| |/ _` | __| |/ _ \| '_ \ 
##     | |\__ \ (_) | | (_| | |_| | (_) | | | |
##    |___|___/\___/|_|\__,_|\__|_|\___/|_| |_|

                                         
#  Isolation ################ 
process.PassingIsolation = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("GsfConvRej"),
    cut = cms.string(process.PassingGsf.cut.value() +
         " && (( isEB && ( dr03TkSumPt/p4.Pt < 0.09 && dr03EcalRecHitSumEt/p4.Pt < 0.07 && dr03HcalTowerSumEt/p4.Pt < 0.10 ))"
         " || (isEE && (dr03TkSumPt/p4.Pt < 0.04 && dr03EcalRecHitSumEt/p4.Pt < 0.05  && dr03HcalTowerSumEt/p4.Pt < 0.025 )))")
#         " && (( isEB && ( (dr03TkSumPt + max(0., dr03EcalRecHitSumEt - 1.) + dr03HcalTowerSumEt)/(p4.Pt) < 0.15 ))"
#         " || (isEE && ((dr03TkSumPt + dr03EcalRecHitSumEt + dr03HcalTowerSumEt)/(p4.Pt) < 0.1 )))")
)

##    _____ _           _                     ___    _ 
##   | ____| | ___  ___| |_ _ __ ___  _ __   |_ _|__| |
##   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \   | |/ _` |
##   | |___| |  __/ (__| |_| | | (_) | | | |  | | (_| |
##   |_____|_|\___|\___|\__|_|  \___/|_| |_| |___\__,_|
##   

# Electron ID  ######

process.PassingId = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("GsfConvRej"),
    cut = cms.string(process.PassingIsolation.cut.value() +
#                     " && (gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1)"
                     " && ((isEB"
                                   " && (sigmaIetaIeta<0.01)"
                                   " && ( -0.06<deltaPhiSuperClusterTrackAtVtx<0.06 )"
                                   " && ( -0.004<deltaEtaSuperClusterTrackAtVtx<0.004 )"
                                   " && (hadronicOverEm<0.04)"
                                   ")"
                     " || (isEE"
                                   " && (sigmaIetaIeta<0.03)"
                                   " && ( -0.03<deltaPhiSuperClusterTrackAtVtx<0.03 )"
                                   " && ( -0.007<deltaEtaSuperClusterTrackAtVtx<0.007 )"
                                   " && (hadronicOverEm<0.025) "
                                   "))"
                     )
)


process.PassingId80 = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag("GsfConvRej"),
    cut = cms.string(process.PassingGsf.cut.value() +
                     " && (gsfTrack.trackerExpectedHitsInner.numberOfHits <= 0)"
                     " && ((isEB"
                     " && ( dr03TkSumPt/p4.Pt <0.09 && dr03EcalRecHitSumEt/p4.Pt < 0.07 && dr03HcalTowerSumEt/p4.Pt  < 0.1 )"
#                                   " && ( (dr03TkSumPt + max(0., dr03EcalRecHitSumEt - 1.) + dr03HcalTowerSumEt)/(p4.Pt) < 0.07 )"
                                   " && (sigmaIetaIeta<0.01)"
                                   " && ( -0.06<deltaPhiSuperClusterTrackAtVtx<0.06 )"
                                   " && ( -0.004<deltaEtaSuperClusterTrackAtVtx<0.004 )"
                                   " && (hadronicOverEm<0.04)"
                                   ")"
                     " || (isEE"
                      " && ( dr03TkSumPt/p4.Pt <0.04 && dr03EcalRecHitSumEt/p4.Pt < 0.05 && dr03HcalTowerSumEt/p4.Pt  < 0.025 )"
#                                   " && ( (dr03TkSumPt + dr03EcalRecHitSumEt + dr03HcalTowerSumEt)/(p4.Pt) < 0.06 )"
                                   " && (sigmaIetaIeta<0.03)"
                                   " && ( -0.03<deltaPhiSuperClusterTrackAtVtx<0.03 )"
                                   " && ( -0.007<deltaEtaSuperClusterTrackAtVtx<0.007 )"
                                   " && (hadronicOverEm<0.025) "
                                   "))"
                     )
)



                         
##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
##   

# Trigger  ##################
process.PassingHLT = cms.EDProducer("trgMatchedGsfElectronProducer",                     
    InputProducer = cms.InputTag("PassingId"),                          
    hltTag = cms.untracked.InputTag(HLTPath,"","HLT"),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT")
)

process.badSuperClustersClean = cms.EDProducer("CandViewCleaner",
    srcObject = cms.InputTag("goodSuperClustersClean"),
    module_label = cms.string(''),
    srcObjectsToRemove = cms.VInputTag(cms.InputTag("PassingHLT")),
    deltaRMin = cms.double(0.1)
)

##    _____      _                        _  __     __             
##   | ____|_  _| |_ ___ _ __ _ __   __ _| | \ \   / /_ _ _ __ ___ 
##   |  _| \ \/ / __/ _ \ '__| '_ \ / _` | |  \ \ / / _` | '__/ __|
##   | |___ >  <| ||  __/ |  | | | | (_| | |   \ V / (_| | |  \__ \
##   |_____/_/\_\\__\___|_|  |_| |_|\__,_|_|    \_/ \__,_|_|  |___/
##   

## Here we show how to use a module to compute an external variable
#process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.load('RecoJets.Configuration.RecoPFJets_cff')
## process.ak5PFJets.doAreaFastjet = True
## process.kt6PFJets.doRhoFastjet = True
process.load('JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff')
## process.L1Fastjet.algorithm = cms.string('AK5Calo') #DUMMY THESE are DUMMY,
## process.L1Fastjet.era = 'Spring10' #DUMMY
## process.L1Fastjet.level = cms.string('L2Relative') #DUMMY
## process.L1Fastjet.useCondDB = cms.untracked.bool(False)
## process.offsetCorrection = cms.Sequence(process.ak5PFJets *
##                                         process.kt6PFJets * process.ak5PFJetsL1)

#from JetMETCorrections.Configuration.DefaultJEC_cff import *
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

# to run the offset corrections on PFJets and L2L3 on top of them
process.PUcorrAk5PFJets=process.ak5PFJets.clone()
process.PUcorrAk5PFJets.doAreaFastjet = True
process.PUcorrKt6PFJets=process.kt6PFJets.clone()
process.PUcorrKt6PFJets.doRhoFastjet = True
process.PUcorrL1Fastjet=process.L1Fastjet.clone();
process.PUcorrL1Fastjet.algorithm = cms.string('AK5Calo') #DUMMY THESE are DUMMY,
process.PUcorrL1Fastjet.era = 'Spring10'                  #DUMMY
process.PUcorrL1Fastjet.level = cms.string('L2Relative')  #DUMMY
process.PUcorrL1Fastjet.useCondDB = cms.untracked.bool(False)
process.PUcorrL1Fastjet.srcMedianPt = 'PUcorrKt6PFJets'
process.PUcorrAk5PFJetsL1 = process.ak5PFJetsL1.clone()
process.PUcorrAk5PFJetsL1.src = 'PUcorrAk5PFJets'
process.PUcorrAk5PFJetsL1.correctors = ['PUcorrL1Fastjet']

process.PUcorrAk5PFJetsL1L2L3Residual   = process.ak5PFJetsL2L3.clone(src = 'PUcorrAk5PFJetsL1', correctors = ['ak5PFL2L3Residual'])
process.PUcorrAk5PFJetsL1L2L3  = process.ak5PFJetsL2L3.clone(src = 'PUcorrAk5PFJetsL1', correctors = ['ak5PFL2L3'])

process.offsetCorrection = cms.Sequence(process.PUcorrAk5PFJets *  process.PUcorrKt6PFJets * process.PUcorrAk5PFJetsL1)

# data sequences use residual corrections
process.CaloJetSequenceData = cms.Sequence( process.ak5CaloJetsL2L3Residual )
process.PFJetAK5SequenceData = cms.Sequence( process.ak5PFJetsL2L3Residual * process.offsetCorrection * process.PUcorrAk5PFJetsL1L2L3Residual )
#process.ourJetSequenceData = cms.Sequence( process.PFJetAK5SequenceData )
process.ourJetSequenceData = cms.Sequence( process.CaloJetSequenceData * process.PFJetAK5SequenceData )



process.cleanJets = cms.EDProducer(
    "JetViewCleaner",
    srcObject = cms.InputTag("PUcorrAk5PFJetsL1L2L3Residual"),
    srcObjectsToRemove = cms.VInputTag( cms.InputTag("gsfElectrons")),
    deltaRMin = cms.double(0.5)
    )

process.superClusterDRToNearestJet = cms.EDProducer("DeltaRNearestJetComputer",
    probes = cms.InputTag("goodSuperClusters"),
       # ^^--- NOTA BENE: if probes are defined by ref, as in this case, 
       #       this must be the full collection, not the subset by refs.
    objects = cms.InputTag("cleanJets"),
    objectSelection = cms.string(JET_CUTS),
)


process.JetMultiplicityInSCEvents = cms.EDProducer("CandMultiplicityCounter",
    probes = cms.InputTag("goodSuperClusters"),
    objects = cms.InputTag("cleanJets"),
    objectSelection = cms.string(JET_CUTS),
)


process.JetHiPtInSCEvents = cms.EDProducer("CandPtReader",
    probes = cms.InputTag("goodSuperClusters"),
    objects = cms.InputTag("cleanJets"),
    objectSelection = cms.string(JET_CUTS),
)


process.GsfDRToNearestJet = cms.EDProducer("DeltaRNearestJetComputer",
    probes = cms.InputTag("GsfConvRej"),
    objects = cms.InputTag("cleanJets"),
    objectSelection = cms.string(JET_CUTS),
)


process.JetMultiplicityInGsfEvents = cms.EDProducer("CandMultiplicityCounter",
    probes = cms.InputTag("GsfConvRej"),
    objects = cms.InputTag("cleanJets"),
    objectSelection = cms.string(JET_CUTS),
)


process.JetHiPtInGsfEvents = cms.EDProducer("CandPtReader",
    probes = cms.InputTag("GsfConvRej"),
    objects = cms.InputTag("cleanJets"),
    objectSelection = cms.string(JET_CUTS),
)


process.ext_ToNearestJet_sequence = cms.Sequence(
    process.ourJetSequenceData +    #process.ak5CaloL2L3 + 
    process.cleanJets+
    process.superClusterDRToNearestJet +
    process.JetMultiplicityInSCEvents +
    process.JetHiPtInSCEvents+
    process.GsfDRToNearestJet +
    process.JetMultiplicityInGsfEvents+
    process.JetHiPtInGsfEvents
    )


##    _____             ____        __ _       _ _   _             
##   |_   _|_ _  __ _  |  _ \  ___ / _(_)_ __ (_) |_(_) ___  _ __  
##     | |/ _` |/ _` | | | | |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \ 
##     | | (_| | (_| | | |_| |  __/  _| | | | | | |_| | (_) | | | |
##     |_|\__,_|\__, | |____/ \___|_| |_|_| |_|_|\__|_|\___/|_| |_|
##              |___/                                              

process.Tag = process.PassingHLT.clone()

process.TagMatchedSuperClusterCandsClean = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("goodSuperClustersClean"),
   ReferenceElectronCollection = cms.untracked.InputTag("Tag"),
   deltaR =  cms.untracked.double(0.3)
)


process.TagMatchedPhotonCands = cms.EDProducer("ElectronMatchedCandidateProducer",
   src     = cms.InputTag("FilteredPhotons"),
   ReferenceElectronCollection = cms.untracked.InputTag("Tag"),
   deltaR =  cms.untracked.double(0.3)
)

process.IsoMatchedSuperClusterCandsClean = process.TagMatchedSuperClusterCandsClean.clone()
process.IsoMatchedSuperClusterCandsClean.ReferenceElectronCollection = cms.untracked.InputTag("PassingIsolation")

process.IdMatchedSuperClusterCandsClean = process.TagMatchedSuperClusterCandsClean.clone()
process.IdMatchedSuperClusterCandsClean.ReferenceElectronCollection = cms.untracked.InputTag("PassingId")

process.Id80MatchedSuperClusterCandsClean = process.TagMatchedSuperClusterCandsClean.clone()
process.Id80MatchedSuperClusterCandsClean.ReferenceElectronCollection = cms.untracked.InputTag("PassingId80")

process.IsoMatchedPhotonCands = process.GsfMatchedPhotonCands.clone()
process.IsoMatchedPhotonCands.ReferenceElectronCollection = cms.untracked.InputTag("PassingIsolation")

process.IdMatchedPhotonCands = process.GsfMatchedPhotonCands.clone()
process.IdMatchedPhotonCands.ReferenceElectronCollection = cms.untracked.InputTag("PassingId")

process.Id80MatchedPhotonCands = process.GsfMatchedPhotonCands.clone()
process.Id80MatchedPhotonCands.ReferenceElectronCollection = cms.untracked.InputTag("PassingId80")


process.ele_sequence = cms.Sequence(
    process.GsfConvRej *(process.PassingGsf * process.GsfMatchedSuperClusterCands +
    process.GsfMatchedPhotonCands +# process.GsfConvRej+
    process.PassingIsolation + process.PassingId +
    process.PassingId80 +
    process.PassingHLT + process.Tag*
    process.TagMatchedSuperClusterCandsClean *
    process.badSuperClustersClean *
    process.TagMatchedPhotonCands *
    process.IsoMatchedSuperClusterCandsClean *
    process.IdMatchedSuperClusterCandsClean *
    process.Id80MatchedSuperClusterCandsClean *
    process.IsoMatchedPhotonCands *
    process.IdMatchedPhotonCands *
    process.Id80MatchedPhotonCands
    )
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
    decay = cms.string("Tag goodSuperClustersClean"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("60 < mass < 120"),
)


process.tagPhoton = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("Tag FilteredPhotons"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("60 < mass < 120"),
)

process.SCSC = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("badSuperClustersClean badSuperClustersClean"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("60 < mass < 120"),
)

process.GsfGsf = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("PassingGsf PassingGsf"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                                   
    cut   = cms.string("60 < mass < 120"),
)

process.tagConR = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("Tag GsfConvRej"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                                   
    cut   = cms.string("60 < mass < 120"),
)

process.tagGsf = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("Tag PassingGsf"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                                   
    cut   = cms.string("60 < mass < 120"),
)


process.tagIso = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("Tag PassingIsolation"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                                   
    cut   = cms.string("60 < mass < 120"),
)

process.tagId = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("Tag PassingId"), # charge coniugate states are implied
    checkCharge = cms.bool(False),
    cut   = cms.string("60 < mass < 120"),
)

process.tagId80 = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("Tag PassingId80"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                                  
    cut   = cms.string("60 < mass < 120"),
)


process.tagHLT = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("Tag PassingHLT"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                                   
    cut   = cms.string("60 < mass < 120"),
)



process.allTagsAndProbes = cms.Sequence(
    process.tagSC + process.SCSC + process.tagPhoton +
    process.tagConR +process.tagGsf + process.GsfGsf +
    process.tagIso + process.tagId +  process.tagId80 +process.tagHLT
)


##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        

process.McMatchTag = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("Tag"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)


process.McMatchSC = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("goodSuperClustersClean"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)


process.McMatchPhoton = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("FilteredPhotons"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)


process.McMatchSCbad = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("badSuperClustersClean"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)


process.McMatchGsf = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("PassingGsf"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)

process.McMatchIso = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("PassingIsolation"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)

process.McMatchId = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("PassingId"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)

process.McMatchHLT = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("PassingHLT"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)




process.mc_sequence = cms.Sequence(
   process.McMatchTag + process.McMatchSC + process.McMatchPhoton +
   process.McMatchGsf + process.McMatchIso +
   process.McMatchId  + process.McMatchHLT
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


## I define some common variables for re-use later.
## This will save us repeating the same code for each efficiency category

ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    pt  = cms.string("pt"),
    phi  = cms.string("phi"),
    et  = cms.string("et"),
    e  = cms.string("energy"),
    rapidity  = cms.string("rapidity"),
    mass  = cms.string("mass"),
)   



ProbeVariablesToStore = cms.PSet(
    probe_gsfEle_abseta = cms.string("abs(eta)"),
    probe_gsfEle_eta = cms.string("eta"),
    probe_gsfEle_pt  = cms.string("pt"),
    probe_gsfEle_phi  = cms.string("phi"),
    probe_gsfEle_et  = cms.string("et"),
    probe_gsfEle_e  = cms.string("energy"),
    probe_gsfEle_charge = cms.string("charge"),
    probe_gsfEle_rapidity  = cms.string("rapidity"),
    ## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    probe_sc_ecalEnergyToEt    = cms.string("ecalEnergy*sin(superClusterPosition.theta)"),    
    probe_sc_eta    = cms.string("superCluster.eta"),
    probe_sc_abseta    = cms.string("abs(superCluster.eta)"),
    probe_sc_phi    = cms.string("superCluster.phi"),

)


TagVariablesToStore = cms.PSet(
    gsfEle_eta = cms.string("eta"),
    gsfEle_pt  = cms.string("pt"),
    gsfEle_phi  = cms.string("phi"),
    gsfEle_et  = cms.string("et"),
    gsfEle_e  = cms.string("energy"),
    gsfEle_charge = cms.string("charge"),
    gsfEle_rapidity  = cms.string("rapidity"),
    ## super cluster quantities
    sc_energy = cms.string("superCluster.energy"),
    sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    sc_eta    = cms.string("superCluster.eta"),
    sc_phi    = cms.string("superCluster.phi"),
    sc_rawEnergy = cms.string("superCluster.rawEnergy"), 
)


ProbePhotonVariablesToStore = cms.PSet(
        probe_abseta = cms.string("abs(eta)"),
        probe_eta = cms.string("eta"),
        probe_phi  = cms.string("phi"),
        probe_et  = cms.string("et"),
)


ProbeSuperClusterVariablesToStore = cms.PSet(
    probe_sc_eta = cms.string("eta"),
    probe_sc_abseta = cms.string("abs(eta)"),
    probe_sc_pt  = cms.string("pt"),
    probe_sc_phi  = cms.string("phi"),
    probe_sc_et  = cms.string("et"),
    probe_sc_e  = cms.string("energy"),
)


TagSuperClusterVariablesToStore = cms.PSet(
    sc_eta = cms.string("eta"),
    sc_pt  = cms.string("pt"),
    sc_phi  = cms.string("phi"),
    sc_et  = cms.string("et"),
    sc_e  = cms.string("energy"),
)





CommonStuffForSuperClusterProbe = cms.PSet(
   variables = cms.PSet(ProbeSuperClusterVariablesToStore),
   ignoreExceptions =  cms.bool (False),
   #fillTagTree      =  cms.bool (True),
   addRunLumiInfo   =  cms.bool (True),
   addEventVariablesInfo   =  cms.bool (True),
   pairVariables =  cms.PSet(ZVariablesToStore),
   pairFlags     =  cms.PSet(
          mass60to120 = cms.string("60 < mass < 120")
    ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags     =  cms.PSet(
          flag = cms.string("pt>0")
    ),    
)






CommonStuffForGsfElectronProbe = cms.PSet(
    variables = cms.PSet(ProbeVariablesToStore),
    ignoreExceptions =  cms.bool (False),
    #fillTagTree      =  cms.bool (True),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool (True),
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
          mass60to120 = cms.string("60 < mass < 120")
    ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags     =  cms.PSet(
          flag = cms.string("pt>0")
    ),    
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
        probe_pt  = cms.string("pt"),
        probe_phi  = cms.string("phi"),
        probe_et  = cms.string("et"),
        probe_e  = cms.string("energy"),
        probe_p  = cms.string("p"),
        probe_px  = cms.string("px"),
        probe_py  = cms.string("py"),
        probe_pz  = cms.string("pz"),
        probe_theta  = cms.string("theta"),    
        probe_vx     = cms.string("vx"),
        probe_vy     = cms.string("vy"),
        probe_vz     = cms.string("vz"),   
        probe_charge = cms.string("charge"),
        probe_rapidity  = cms.string("rapidity"),    
        probe_mass  = cms.string("mass"),
        probe_mt  = cms.string("mt"),    
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

## super cluster --> gsf electron
process.SCToGsf = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    CommonStuffForSuperClusterProbe, mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration                 
    tagProbePairs = cms.InputTag("tagSC"),
    arbitration   = cms.string("Random2"),                      
    flags = cms.PSet(
        probe_passing = cms.InputTag("GsfMatchedSuperClusterCands"),
        probe_passingGsf = cms.InputTag("GsfMatchedSuperClusterCands"),        
        probe_passingIso = cms.InputTag("IsoMatchedSuperClusterCandsClean"),
        probe_passingId = cms.InputTag("IdMatchedSuperClusterCandsClean"), 
        probe_passingId80 = cms.InputTag("Id80MatchedSuperClusterCandsClean"),       
        probe_passingALL = cms.InputTag("TagMatchedSuperClusterCandsClean")
    ),
    probeMatches  = cms.InputTag("McMatchSC"),
    allProbes     = cms.InputTag("goodSuperClustersClean")
)
process.SCToGsf.variables.probe_dRjet = cms.InputTag("superClusterDRToNearestJet")
process.SCToGsf.variables.probe_nJets = cms.InputTag("JetMultiplicityInSCEvents")
process.SCToGsf.variables.probe_highPtJet = cms.InputTag("JetHiPtInSCEvents")





process.SCSCtoTagSC = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
   variables = cms.PSet(ProbeSuperClusterVariablesToStore),
   ignoreExceptions =  cms.bool (False),
   addRunLumiInfo   =  cms.bool (True),
   addEventVariablesInfo   =  cms.bool (True),
   pairVariables =  cms.PSet(ZVariablesToStore),
   pairFlags     =  cms.PSet(
          mass60to120 = cms.string("60 < mass < 120")
    ),
    tagVariables   =  cms.PSet(TagSuperClusterVariablesToStore),
    tagFlags     =  cms.PSet(
          flag = cms.string("pt>0")
    ),                                         
    isMC = cms.bool(False),
    #mcTruthCommonStuff,
    # choice of tag and probe pairs, and arbitration                      
    tagProbePairs = cms.InputTag("SCSC"),
    arbitration   = cms.string("Random2"),
    massForArbitration = cms.double(91.1876),
    flags = cms.PSet(
          probe_passing = cms.InputTag("TagMatchedSuperClusterCandsClean")
    ),
    probeMatches  = cms.InputTag("McMatchSCbad"),         
    allProbes     = cms.InputTag("badSuperClustersClean")
)


## good photon --> gsf electron
process.PhotonToGsf = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    mcTruthCommonStuff,
    CommonStuffForSuperClusterProbe,
    # choice of tag and probe pairs, and arbitration                 
    tagProbePairs = cms.InputTag("tagPhoton"),
    arbitration   = cms.string("Random2"),                      
    flags = cms.PSet(
        probe_passing = cms.InputTag("GsfMatchedPhotonCands"),
        probe_passingALL = cms.InputTag("TagMatchedPhotonCands"),
        probe_passingIso = cms.InputTag("IsoMatchedPhotonCands"),
        probe_passingId = cms.InputTag("IdMatchedPhotonCands"),
        probe_passingId80 = cms.InputTag("Id80MatchedPhotonCands")
    ),
    probeMatches  = cms.InputTag("McMatchPhoton"),
    allProbes     = cms.InputTag("FilteredPhotons")
)
process.PhotonToGsf.variables=ProbePhotonVariablesToStore


process.SCSCbad = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
   #######mcTruthCommonStuff,
   variables = cms.PSet(ProbeSuperClusterVariablesToStore),
   ignoreExceptions =  cms.bool (False),
   addRunLumiInfo   =  cms.bool (True),
   addEventVariablesInfo   =  cms.bool (True),
   pairVariables =  cms.PSet(ZVariablesToStore),
   pairFlags     =  cms.PSet(
          mass60to120 = cms.string("60 < mass < 120")
          ),
   tagVariables   =  cms.PSet(TagSuperClusterVariablesToStore),
   tagFlags     =  cms.PSet(
          flag = cms.string("pt>0")
   ),                                         
   isMC = cms.bool(False),
   # choice of tag and probe pairs, and arbitration                      
   tagProbePairs = cms.InputTag("SCSC"),
   arbitration   = cms.string("Random2"),
   massForArbitration = cms.double(91.1876),
   flags = cms.PSet(
          probe_passing = cms.InputTag("TagMatchedSuperClusterCandsClean")
   ),
   #probeMatches  = cms.InputTag("McMatchSCbad"),         
   allProbes     = cms.InputTag("badSuperClustersClean")
)

process.GsfGsfToIso = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ########mcTruthCommonStuff,
    CommonStuffForGsfElectronProbe,
    isMC = cms.bool(False), 
    tagProbePairs = cms.InputTag("GsfGsf"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
        probe_passing = cms.InputTag("PassingIsolation")
    ),
    #probeMatches  = cms.InputTag("McMatchGsf"),
    allProbes     = cms.InputTag("PassingGsf")
)


##     ____      __       __    ___           
##    / ___|___ / _|      \ \  |_ _|___  ___  
##   | |  _/ __| |_   _____\ \  | |/ __|/ _ \ 
##   | |_| \__ \  _| |_____/ /  | |\__ \ (_) |
##    \____|___/_|        /_/  |___|___/\___/ 
##   
##  gsf electron --> isolation

process.GsfToIso = cms.EDAnalyzer("TagProbeFitTreeProducer",
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,                        
    tagProbePairs = cms.InputTag("tagGsf"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
        probe_passing = cms.InputTag("PassingIsolation"),
        probe_passingIso = cms.InputTag("PassingIsolation"),
        probe_passingId = cms.InputTag("PassingId"),
        probe_passingId80 = cms.InputTag("PassingId80"),        
        probe_passingALL = cms.InputTag("PassingHLT")        
    ),
    probeMatches  = cms.InputTag("McMatchGsf"),
    allProbes     = cms.InputTag("PassingGsf")
)
process.GsfToIso.variables.probe_dRjet = cms.InputTag("GsfDRToNearestJet")
process.GsfToIso.variables.probe_nJets = cms.InputTag("JetMultiplicityInGsfEvents")
process.GsfToIso.variables.probe_highPtJet = cms.InputTag("JetHiPtInGsfEvents")

##    ___                 __    ___    _ 
##   |_ _|___  ___        \ \  |_ _|__| |
##    | |/ __|/ _ \   _____\ \  | |/ _` |
##    | |\__ \ (_) | |_____/ /  | | (_| |
##   |___|___/\___/       /_/  |___\__,_|
##   
##  isolation --> Id

process.IsoToId = cms.EDAnalyzer("TagProbeFitTreeProducer",
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,                              
    tagProbePairs = cms.InputTag("tagIso"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
        probe_passing = cms.InputTag("PassingId"),
        probe_passingId = cms.InputTag("PassingId"),
        probe_passingId80 = cms.InputTag("PassingId80"),        
        probe_passingALL = cms.InputTag("PassingHLT")         
    ),
    probeMatches  = cms.InputTag("McMatchIso"),
    allProbes     = cms.InputTag("PassingIsolation")
)
process.IsoToId.variables.probe_dRjet = cms.InputTag("GsfDRToNearestJet")
process.IsoToId.variables.probe_nJets = cms.InputTag("JetMultiplicityInGsfEvents")
process.IsoToId.variables.probe_highPtJet = cms.InputTag("JetHiPtInGsfEvents")

##    ___    _       __    _   _ _   _____ 
##   |_ _|__| |      \ \  | | | | | |_   _|
##    | |/ _` |  _____\ \ | |_| | |   | |  
##    | | (_| | |_____/ / |  _  | |___| |  
##   |___\__,_|      /_/  |_| |_|_____|_|  

##  Id --> HLT
process.IdToHLT = cms.EDAnalyzer("TagProbeFitTreeProducer",
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,                             
    tagProbePairs = cms.InputTag("tagId"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
        probe_passing = cms.InputTag("PassingHLT"),
        probe_passingId80 = cms.InputTag("PassingId80")        
    ),
    probeMatches  = cms.InputTag("McMatchId"),
    allProbes     = cms.InputTag("PassingId")
)



process.Id80ToHLT = cms.EDAnalyzer("TagProbeFitTreeProducer",
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,
    tagProbePairs = cms.InputTag("tagId80"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
        probe_passing = cms.InputTag("PassingHLT"),
        probe_passingId80 = cms.InputTag("PassingId80")
    ),
    probeMatches  = cms.InputTag("McMatchId"),
    allProbes     = cms.InputTag("PassingId80")
)


process.IdToHLT.variables.probe_dRjet = cms.InputTag("GsfDRToNearestJet")
process.IdToHLT.variables.probe_nJets = cms.InputTag("JetMultiplicityInGsfEvents")
process.IdToHLT.variables.probe_highPtJet = cms.InputTag("JetHiPtInGsfEvents")

process.Id80ToHLT.variables.probe_dRjet = cms.InputTag("GsfDRToNearestJet")
process.Id80ToHLT.variables.probe_nJets = cms.InputTag("JetMultiplicityInGsfEvents")
process.Id80ToHLT.variables.probe_highPtJet = cms.InputTag("JetHiPtInGsfEvents")


 
process.tree_sequence = cms.Sequence(
    process.SCToGsf +
    process.SCSCbad +
    process.PhotonToGsf +
    process.GsfToIso +
    process.GsfGsfToIso +
    process.IsoToId + 
    process.IdToHLT+ process.Id80ToHLT
)    

##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##
process.tagAndProbe = cms.Path(
#    process.gsfElectrons +    
    process.sc_sequence + process.ele_sequence +
    process.ext_ToNearestJet_sequence + 
    process.allTagsAndProbes +
#    process.mc_sequence + 
    process.tree_sequence
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(OUTPUT_FILE_NAME)
                                   )
