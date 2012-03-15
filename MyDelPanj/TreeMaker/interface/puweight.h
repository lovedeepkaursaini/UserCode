#ifndef __puweight__
#define __puweight__

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


class puweight{

 public:
  puweight(std::string name, TTree* tree);
  ~puweight();
  void Fill(const edm::Event& iEvent,bool data); 
  void Clear();
 private:
  puweight(){};//Don't allow user
  void SetBranches();
  void AddBranch(double* x, std::string name);
    
  TTree *tree_;
  std::string name_;

  edm::LumiReWeighting LumiWeights_;

std::vector< float > Data2011_;
std::vector<float> Fall2011_;

double PUweight_;
};

#endif

