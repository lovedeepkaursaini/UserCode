#include "UserCode/TrigTree/interface/L1Tree.h"

L1Tree::L1Tree(std::string name,TTree* tree){
  name_=name;
  tree_=tree;

  tree_->Branch("L1trigResults",&L1trigResult_);
  tree_->Branch("L1trigErrCode",&L1trigErrCode_);
  tree_->Branch("L1trigName",&L1trigName_);

}

void
L1Tree::Fill(const edm::Event& iEvent,const edm::EventSetup& iSetup){
  
  m_l1GtUtils.retrieveL1EventSetup(iSetup);
      // L1 Trigger
      m_l1GtUtils.retrieveL1EventSetup(iSetup);
      int iErrorCode = -1;

      std::vector< std::string > l1TrigNames;
      l1TrigNames.push_back("L1_SingleIsoEG5");
      l1TrigNames.push_back("L1_SingleIsoEG8");
      l1TrigNames.push_back("L1_SingleIsoEG10");
      l1TrigNames.push_back("L1_SingleIsoEG12");
      l1TrigNames.push_back("L1_SingleIsoEG15");
      l1TrigNames.push_back("L1_SingleEG2");
      l1TrigNames.push_back("L1_SingleEG5");
      l1TrigNames.push_back("L1_SingleEG8");
      l1TrigNames.push_back("L1_SingleEG10");
      l1TrigNames.push_back("L1_SingleEG12");
      l1TrigNames.push_back("L1_SingleEG15");
      l1TrigNames.push_back("L1_SingleEG20"); 

      std::map<std::string,bool> l1TriggersOfInterest_;

      for (unsigned int n = 0; n < l1TrigNames.size(); n++)
	{ 
	  bool decision = 
	    m_l1GtUtils.decisionBeforeMask
	    (iEvent, l1TrigNames.at(n), iErrorCode);
	  
	  L1trigName_.push_back(l1TrigNames.at(n));
	  L1trigResult_.push_back(decision);
	  L1trigErrCode_.push_back(iErrorCode);

	  /*if iErrCode=0 (fine), =1 
	  (Trig. doesn't exist in this L1 menu), 
	  =else, an error with trigger.*/
	}
      
}


void
L1Tree::Clear(){
  L1trigResult_.clear();
  L1trigErrCode_.clear();
  L1trigName_.clear();
}
