// -*- C++ -*-
//
// Package:    ElecTrig
// Class:      ElecTrig
// 
/**\class ElecTrig ElecTrig.cc UserCode/ElecTrig/src/ElecTrig.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lovedeep Kaur Saini,40 4-B28,+41227671589,
//         Created:  Tue May 11 00:25:55 CEST 2010
// $Id$
//
//

#include <boost/algorithm/string.hpp>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Common/interface/TriggerNames.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "UserCode/ElecTrig/interface/ElecSelector.h"
#include "UserCode/TrigTree/interface/ElecTree.h"
#include "UserCode/TrigTree/interface/HltTree.h"
#include "UserCode/TrigTree/interface/SCTree.h"
#include "UserCode/TrigTree/interface/L1Tree.h"
#include "UserCode/TrigTree/interface/JetUtil.h"
#include "UserCode/TrigTree/interface/EvtInfoTree.h"

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "UserCode/TrigTree/interface/ElecTree.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//
// class declaration
//

class TrigTree : public edm::EDAnalyzer {
public:
  explicit TrigTree(const edm::ParameterSet&);
  ~TrigTree();
  
  // private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------

  TTree* tree;
  EvtInfoTree* evtInfoTree_; 
  HltTree* hltTree_;
  ElecTree* elecTree_;
  SCTree* scTree_;
  L1Tree* l1Tree_;
  
  edm::InputTag elecTag_;
  const std::vector<std::string> elecTrig_;
  bool fillElectrons_;
  bool fillSuperClusters_;
  bool fillHltInfo_;
  bool fillL1Info_;
  bool fillEvtInfo_;

};

//
// constructors and destructor
//
TrigTree::TrigTree(const edm::ParameterSet& iConfig):
  elecTag_(iConfig.getUntrackedParameter<edm::InputTag> 
	   ("ElecTag", edm::InputTag("gsfElectrons"))),
  fillElectrons_(iConfig.getUntrackedParameter<bool> ("fillElec",1)),
  fillSuperClusters_(iConfig.getUntrackedParameter<bool> ("fillSC",1)),
  fillHltInfo_(iConfig.getUntrackedParameter<bool> ("fillHlt",1)),
  fillL1Info_(iConfig.getUntrackedParameter<bool> ("fillL1",1)),
  fillEvtInfo_(iConfig.getUntrackedParameter<bool> ("fillEvt",1)){
  
  edm::Service<TFileService> fs;
  TFileDirectory treeDir = fs->mkdir("Summary");
  tree=new TTree("TrigTree","TrigTree");  
  if(fillEvtInfo_)evtInfoTree_=new EvtInfoTree("TrigAnls",tree);
  if(fillElectrons_)  elecTree_=new ElecTree("TrigAnls",tree,iConfig);
  if(fillHltInfo_)hltTree_=new HltTree("TrigAnls",tree);
  if(fillSuperClusters_)scTree_=new SCTree("TrigAnls",tree,iConfig);
  if(fillL1Info_)l1Tree_=new L1Tree("TrigAnls",tree);
}



TrigTree::~TrigTree(){
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TrigTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  if(fillEvtInfo_)evtInfoTree_->Clear();
  if(fillElectrons_) elecTree_->Clear();
  if(fillHltInfo_)hltTree_->Clear();
  if(fillL1Info_)l1Tree_->Clear();
  if(fillSuperClusters_)scTree_->Clear();

  //evt info
  if(fillEvtInfo_)evtInfoTree_->EvtInfoCollector(iEvent);
  
  //Trigger Info
  if(fillHltInfo_)hltTree_->Fill(iEvent);
  
  //L1 Info
  if(fillL1Info_)l1Tree_->Fill(iEvent,iSetup);
  
  //Fill SuperClusters
  if(fillSuperClusters_)scTree_->Fill(iEvent,iSetup);
 

  // Beam spot
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  if (!iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), beamSpotHandle)) {
    std::cout<<">>> No beam spot found !!!"<<std::endl;
    return;
  }
 
  //gsf electron info
  edm::Handle<reco::GsfElectronCollection> eleCollH;
  if (!iEvent.getByLabel(elecTag_, eleCollH)) 
    {std::cout<<"No ELectron Collection\n"; return;}
  
  //Sort the Electron Collection
  reco::GsfElectronCollection eleColl(*(eleCollH.product()));
  std::sort(eleColl.begin(),eleColl.end(),eleScEtGreater());
  reco::GsfElectronCollection::const_iterator gIter;
  if(fillElectrons_)
    for(gIter=eleColl.begin(); gIter!=eleColl.end();gIter++)
      elecTree_->ElecInfoCollector(iEvent,iSetup,*gIter,beamSpotHandle);
     
  
  tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
TrigTree::beginJob(){

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrigTree::endJob() {

}

//define this as a plug-in
DEFINE_FWK_MODULE(TrigTree);
