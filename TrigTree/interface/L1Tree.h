#ifndef __CRAZY_PUNJABI7__
#define __CRAZY_PUNJABI7__


#include<iostream>
#include<vector>
#include<string>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"


class L1Tree{

 public:
  L1Tree(std::string name,TTree* tree);
  void Fill(const edm::Event& iEvent,const edm::EventSetup& iSetup);
  void Clear();
 private:
  L1Tree(){};
  L1GtUtils m_l1GtUtils;
  
  std::string name_;
  TTree* tree_;
  std::vector<bool> L1trigResult_;
  std::vector<bool> L1trigErrCode_;
  std::vector<std::string> L1trigName_;
};

























#endif
