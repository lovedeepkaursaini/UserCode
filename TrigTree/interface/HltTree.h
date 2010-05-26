#ifndef __CRAZY_PUNJABI5__
#define __CRAZY_PUNJABI5__


#include<iostream>
#include<vector>
#include<string>
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
  std::vector<bool> trigResult_;
  std::vector<std::string> trigName_;
};

























#endif
