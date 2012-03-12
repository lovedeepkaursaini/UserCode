import FWCore.ParameterSet.Config as cms

# DQM monitor module for EWK-WMuNu
ewkZeeDQM = cms.EDAnalyzer(
    "EwkElecDQM",
    TrigTag = cms.untracked.InputTag("TriggerResults::HLT"),
    ElecTag = cms.untracked.InputTag("gsfElectrons"),
    METTag = cms.untracked.InputTag("met"),
    JetTag = cms.untracked.InputTag("ak5CaloJets"),
    TrigTags = cms.untracked.vstring("Ele","DoubleEle"),
    TagMulti = cms.untracked.vint32(2,1),
    VetoTheseTags=cms.untracked.vstring("Mu","Jet"),
    
    NJetMax = cms.untracked.int32(99999),
    NJetMin = cms.untracked.int32(0),
    NEleMax = cms.untracked.int32(9999),
    NEleMin = cms.untracked.int32(2),
    MetMax  = cms.untracked.double(9999),
    MetMin  = cms.untracked.double(0),
    MassMax = cms.untracked.double(9999),
    MassMin = cms.untracked.double(0),
    ElePtCut = cms.untracked.double(15),
    EleEtaCut = cms.untracked.double(2.5),
    SieieBarrel = cms.untracked.double(0.01),
    SieieEndcap = cms.untracked.double(0.028),
    DetainBarrel = cms.untracked.double(0.0071),
    DetainEndcap = cms.untracked.double(0.0066),
    EcalIsoCutBarrel = cms.untracked.double(5.7),
    EcalIsoCutEndcap = cms.untracked.double(5.0),
    HcalIsoCutBarrel = cms.untracked.double(8.1),
    HcalIsoCutEndcap = cms.untracked.double(3.4),
    TrkIsoCutBarrel = cms.untracked.double(7.2),
    TrkIsoCutEndcap = cms.untracked.double(5.1),
    
    JetPtCut = cms.untracked.double(30),
    JetEtaCut = cms.untracked.double(3)
    )


ewkWeNuDQM = ewkZeeDQM.clone()
ewkWeNuDQM.TrigTags = cms.untracked.vstring("Ele")
ewkWeNuDQM.TagMulti = cms.untracked.vint32(1)
ewkWeNuDQM.VetoTheseTags=cms.untracked.vstring("Mu","Jet", "DoubleEle")
ewkWeNuDQM.NEleMin = cms.untracked.int32(1)
ewkWeNuDQM.NEleMax = cms.untracked.int32(1)
ewkWeNuDQM.MetMin = cms.untracked.double(0)



