#ifndef hlttree_h
#define hlttree_h

#include<iostream>
#include<string>
#include<vector>

#include "TTree.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"

class HltTree{

 public:
  HltTree(std::string name,TTree* tree);
  void Fill(const edm::Event& iEvent);
  void Clear();
 private:
  HltTree(){};
 
  std::string name_;
  TTree* tree_;
  std::vector<int> trigResult_;
  std::vector<std::string> trigName_;
  
  std::vector<double> hltEle10LWPt_;
  std::vector<double> hltEle10LWEta_;
  std::vector<double> hltEle10LWPhi_;
  std::vector<double> hltEle15LWPt_;
  std::vector<double> hltEle15LWEta_;
  std::vector<double> hltEle15LWPhi_;
  std::vector<double> hltPhoton10Pt_ ;
  std::vector<double> hltPhoton10Eta_ ;
  std::vector<double> hltPhoton10Phi_ ;
  //----------------------------------------
  std::vector<double> hltPhoton10now15Pt_ ;
  std::vector<double> hltPhoton10now15Eta_ ;
  std::vector<double> hltPhoton10now15Phi_ ;
  //-----------------------------------------
  std::vector<double> hltPhoton15Pt_ ;
  std::vector<double> hltPhoton15Eta_ ;
  std::vector<double> hltPhoton15Phi_ ;
};

#endif
