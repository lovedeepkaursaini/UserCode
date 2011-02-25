#ifndef _GEN_INFO_TREE_H_
#define _GEN_INFO_TREE_H_

#include <memory>
#include <vector>
#include <string>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

//
// class declaration
//

class genInfoTree{

 public:
  genInfoTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~genInfoTree();
  void Fill(const edm::Event& iEvent);
  void SetBranches();
  void Clear();
  
  
  // ----------member data ---------------------------
  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<int>*, std::string name);
  
  edm::InputTag genPartLabel_;
  
  
  TTree* tree_;

  std::vector<double> genMuPx_;
  std::vector<double> genMuPy_;
  std::vector<double> genMuPz_;
  std::vector<double> genMuE_;
  std::vector<double> genMuP_;
  std::vector<double> genMuTheta_;
  std::vector<double> genMuPt_;
  std::vector<double> genMuEta_;
  std::vector<double> genMuPhi_;
  std::vector<double> genMuEt_;
  std::vector<double> genMuQ_;
  std::vector<double> genElePx_;
  std::vector<double> genElePy_;
  std::vector<double> genElePz_;
  std::vector<double> genEleE_;
  std::vector<double> genEleP_;
  std::vector<double> genEleTheta_;
  std::vector<double> genElePt_;
  std::vector<double> genEleEta_;
  std::vector<double> genElePhi_;
  std::vector<double> genEleEt_;
  std::vector<double> genEleQ_;
  std::vector<double> genJetPx_;
  std::vector<double> genJetPy_;
  std::vector<double> genJetPz_;
  std::vector<double> genJetE_;
  std::vector<double> genJetP_;
  std::vector<double> genJetTheta_;
  std::vector<double> genJetPt_;
  std::vector<double> genJetEta_; 
  std::vector<double> genJetPhi_; 
  std::vector<double> genJetEt_;
  std::vector<double> genJetQ_;
};

#endif
