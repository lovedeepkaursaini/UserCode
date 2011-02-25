#include "UserCode/TrigTree/interface/HltTree.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h" 

HltTree::HltTree(std::string name,TTree* tree)
{
  name_=name;
  tree_=tree;

  tree_->Branch("trigResults",&trigResult_);
  tree_->Branch("trigName",&trigName_);
  
  tree_->Branch("hltEle10LWPt",      &hltEle10LWPt_);
  tree_->Branch("hltEle10LWEta",     &hltEle10LWEta_);
  tree_->Branch("hltEle10LWPhi",     &hltEle10LWPhi_);
  tree_->Branch("hltEle15LWPt",      &hltEle15LWPt_);
  tree_->Branch("hltEle15LWEta",     &hltEle15LWEta_);
  tree_->Branch("hltEle15LWPhi",     &hltEle15LWPhi_);
  tree_->Branch("hltPhoton10Pt",     &hltPhoton10Pt_);
  tree_->Branch("hltPhoton10Eta",    &hltPhoton10Eta_);
  tree_->Branch("hltPhoton10Phi",    &hltPhoton10Phi_);
  //Create New Path ;)-------------------------------
  tree_->Branch("hltPhoton10now15Pt",     &hltPhoton10now15Pt_);
  tree_->Branch("hltPhoton10now15Eta",    &hltPhoton10now15Eta_);
  tree_->Branch("hltPhoton10now15Phi",    &hltPhoton10now15Phi_);
  //-------------------------------------------------
  tree_->Branch("hltPhoton15Pt",     &hltPhoton15Pt_);
  tree_->Branch("hltPhoton15Eta",    &hltPhoton15Eta_);
  tree_->Branch("hltPhoton15Phi",    &hltPhoton15Phi_);
}

void
HltTree::Fill(const edm::Event& iEvent)
{
  using namespace edm;
  
  edm::Handle<edm::TriggerResults> trigResults;

  edm::InputTag trigTag("TriggerResults::HLT");
  if (not iEvent.getByLabel(trigTag, trigResults)) {
    std::cout << ">>> TRIGGER collection does not exist !!!\n";
    return;
  }

  const edm::TriggerNames & trigNames = iEvent.triggerNames(*trigResults);

  for (unsigned int i=0; i<trigResults->size(); i++)
  {
    std::string trigName = trigNames.triggerName(i);
    trigName_.push_back(trigName);
    int trigResult = trigResults->accept(i); //bool not to use
    trigResult_.push_back(trigResult);
  }
  
  // get HLT candiates
  edm::Handle<trigger::TriggerEvent> trgEvent;
  iEvent.getByLabel(InputTag("hltTriggerSummaryAOD","","HLT"), trgEvent);   
  const trigger::TriggerObjectCollection& TOC(trgEvent->getObjects());   
  
  // HLT_Ele10_LW_L1R
  edm::InputTag myLastFilter("hltL1NonIsoHLTNonIsoSingleElectronLWEt10PixelMatchFilter","","HLT");
  
  // filterIndex must be less than the size of trgEvent or you get a CMSException: _M_range_check
  if ( trgEvent->filterIndex(myLastFilter) < trgEvent->sizeFilters() ) 
  {
    const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myLastFilter) ) );

    for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) 
    {
      trigger::size_type hltf = keys[hlto];
      const trigger::TriggerObject& L3obj(TOC[hltf]);
      
      hltEle10LWPt_.push_back( L3obj.pt() );
      hltEle10LWEta_.push_back( L3obj.eta() );
      hltEle10LWPhi_.push_back( L3obj.phi() );
    }
  }

  // HLT_Ele15_LW_L1R
  myLastFilter = edm::InputTag("hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter","","HLT");
  
  // filterIndex must be less than the size of trgEvent or you get a CMSException: _M_range_check
  if ( trgEvent->filterIndex(myLastFilter) < trgEvent->sizeFilters() ) 
  {
    const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myLastFilter) ) );

    for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) 
    {
      trigger::size_type hltf = keys[hlto];
      const trigger::TriggerObject& L3obj(TOC[hltf]);
      
      hltEle15LWPt_.push_back( L3obj.pt() );
      hltEle15LWEta_.push_back( L3obj.eta() );
      hltEle15LWPhi_.push_back( L3obj.phi() );
    }
  }


  // HLT_Photon10_L1R
  myLastFilter = edm::InputTag("hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter","","HLT");
  if ( trgEvent->filterIndex(myLastFilter) < trgEvent->sizeFilters() ) 
  {
    const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myLastFilter) ) );

    for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) 
    {
      trigger::size_type hltf = keys[hlto];
      const trigger::TriggerObject& L3obj(TOC[hltf]);
      
      hltPhoton10Pt_.push_back( L3obj.pt() );
      hltPhoton10Eta_.push_back( L3obj.eta() );
      hltPhoton10Phi_.push_back( L3obj.phi() );
//create new trig...
      if(L3obj.pt()>15.)      
      {
      hltPhoton10now15Pt_.push_back( L3obj.pt() );
      hltPhoton10now15Eta_.push_back( L3obj.eta() );
      hltPhoton10now15Phi_.push_back( L3obj.phi() );
      }
    }
  }


  // HLT_Photon15_L1R
  myLastFilter = edm::InputTag("hltL1NonIsoHLTNonIsoSinglePhotonEt15HcalIsolFilter","","HLT");
  if ( trgEvent->filterIndex(myLastFilter) < trgEvent->sizeFilters() ) 
  {
    const trigger::Keys& keys( trgEvent->filterKeys( trgEvent->filterIndex(myLastFilter) ) );

    for ( unsigned int hlto = 0; hlto < keys.size(); hlto++ ) 
    {
      trigger::size_type hltf = keys[hlto];
      const trigger::TriggerObject& L3obj(TOC[hltf]);
      
      hltPhoton15Pt_.push_back( L3obj.pt() );
      hltPhoton15Eta_.push_back( L3obj.eta() );
      hltPhoton15Phi_.push_back( L3obj.phi() );
    }
  }

  
}


void
HltTree::Clear(){
  trigResult_.clear();
  trigName_.clear();
  
  hltEle10LWPt_.clear();
  hltEle10LWEta_.clear();
  hltEle10LWPhi_.clear();
  hltEle15LWPt_.clear();
  hltEle15LWEta_.clear();
  hltEle15LWPhi_.clear();
  hltPhoton10Pt_.clear();
  hltPhoton10Eta_.clear();
  hltPhoton10Phi_.clear();

  //New Trigger Drama
  hltPhoton10now15Pt_.clear();
  hltPhoton10now15Eta_.clear();
  hltPhoton10now15Phi_.clear();



  hltPhoton15Pt_.clear();
  hltPhoton15Eta_.clear();
  hltPhoton15Phi_.clear();
}
