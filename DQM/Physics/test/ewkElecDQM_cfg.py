import FWCore.ParameterSet.Config as cms

process = cms.Process("EwkDQM")
process.load("DQM.Physics.ewkElecDQM_cfi")
process.ewkZeeDQM.TrigTag = cms.untracked.InputTag("TriggerResults::HLT")
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.DQM.collectorHost = ''

process.dqmSaver.workflow = cms.untracked.string('/Physics/EWK/Elec')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    #input = cms.untracked.int32(1000)
    )

## Source
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import pickRelValInputFiles
#process.source = cms.Source(
#    "PoolSource",
#    fileNames = cms.untracked.vstring(
#    pickRelValInputFiles(

  #  cmsswVersion  = 'CMSSW_4_1_3'
  #  , relVal        = 'RelValWE'
  #  , globalTag     = 'START311_V2'
  #  , numberOfFiles = 1
  #  )
  #  )
  #  )

process.source = cms.Source(
"PoolSource",
                               fileNames = cms.untracked.vstring(
#'file:/tmp/lovedeep/929A72DD-DD0E-E011-8BD6-00215E2216AA.root'
#'rfio:/castor/cern.ch/user/l/lovedeep/Spring11_DYToEE_M-20_TuneZ2_7TeV-pythia6_AODSIM_PU_S2_START311_V2-v2_0000_008DBF23-BE5D-E011-B028-E41F13181CF8.root'
  'rfio:/castor/cern.ch/user/l/lovedeep/store_mc_Summer11_DYToEE_M-20_TuneZ2_7TeV-pythia6_AODSIM_PU_S3_START42_V11-v2_0000_FEB5D990-917C-E011-A37F-003048D3C8D6.root'
                # 'rfio:/castor/cern.ch/user/l/lovedeep/Spring11_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_AODSIM_PU_S1_START311_V1G1-v1_0004_007E0957-964E-E011-BA72-485B39800BBB.root'
    #pickRelValInputFiles( cmsswVersion  = 'CMSSW_4_1_3'
    #                         , relVal        = 'RelValZEE'
    #                         , globalTag     = 'START311_V2'
     #                        , numberOfFiles = 1
     #                        )
        )
                            )

process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring('detailedInfo'),
    detailedInfo = cms.untracked.PSet(
    default = cms.untracked.PSet( limit = cms.untracked.int32(10)),
    threshold = cms.untracked.string('INFO')
    )
    )
process.p = cms.Path(process.ewkZeeDQM+process.dqmSaver) # 
#process.p = cms.Path(process.ewkWeNuDQM+process.dqmSaver)


                                
