#include "UserCode/TrigTree/interface/HltTree.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h" 
HltTree::HltTree(std::string name,TTree* tree){
  name_=name;
  tree_=tree;

  tree_->Branch("trigResults",&trigResult_);
  tree_->Branch("trigName",&trigName_);

}

void
HltTree::Fill(const edm::Event& iEvent){
  
  edm::Handle<edm::TriggerResults> trigResults;
//  edm::TriggerNames trigNames;
  edm::InputTag trigTag("TriggerResults::HLT");
  if (not iEvent.getByLabel(trigTag, trigResults)) {
    std::cout << ">>> TRIGGER collection does not exist !!!\n";
    return;
  }

  const edm::TriggerNames & trigNames = iEvent.triggerNames(*trigResults);
//old trigNames.init(*trigResults); 
  for (unsigned int i=0; i<trigResults->size(); i++){
    std::string trigName = trigNames.triggerName(i);
    trigName_.push_back(trigName);
    bool trigResult = trigResults->accept(i);
    trigResult_.push_back(trigResult);
  }
  
}


void
HltTree::Clear(){
  trigResult_.clear();
  trigName_.clear();
}
