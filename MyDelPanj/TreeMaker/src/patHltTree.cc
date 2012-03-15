#include "DelPanj/TreeMaker/interface/patHltTree.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h" 

patHltTree::patHltTree(std::string name,TTree* tree)
{
  name_=name;
  tree_=tree;

  tree_->Branch("trigResults",&trigResult_);
  tree_->Branch("trigName",&trigName_);
  
}

void
patHltTree::Fill(const edm::Event& iEvent)
{
  Clear();
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
   size_t found=trigName.find("Ele");
   if (found==std::string::npos)continue;
    
    trigName_.push_back(trigName);
    int trigResult = trigResults->accept(i); //bool not to use
    trigResult_.push_back(trigResult);
  }
}


void
patHltTree::Clear(){
  trigResult_.clear();
  trigName_.clear();
}


