#ifndef __ZJET_PLANT_H_
#define __ZJET_PLANT_H_
/*
  Author: Anil P Singh
  NCU (TW)
*/
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DelPanj/TreeMaker/interface/baseTree.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DelPanj/TreeMaker/interface/eSelector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/eHist.h"
using namespace std;
using namespace edm;



class patZJetPlant  : public baseTree{

 public:
  patZJetPlant(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~patZJetPlant();
  void Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup) ; 
  void SetBranches();
  void Clear();
  
 private:
  /*  bool DescendingOrder(const TLorentzVector& l1, const TLorentzVector& l2){
    return l1.Pt()>l2.Pt();
  }
  
  bool AscendingOrder(const TLorentzVector& l1, const TLorentzVector& l2){
    return l1.Pt()<l2.Pt();
    }*/
  patZJetPlant();
  edm::InputTag patJetLabel_;
  edm::InputTag patElecLabel_;
  eSelector ewp1_;//for lead ele.
  eSelector ewp2_;//for sublead ele
  
  //Following Variables Go Inside the tree
  double pxGenZ  ;
  double pyGenZ  ;
  double pzGenZ  ;
  double enGenZ  ;
  double pxD1GenZ  ;
  double pyD1GenZ  ;
  double pzD1GenZ  ;
  double enD1GenZ  ;
  double pxD2GenZ  ;
  double pyD2GenZ  ;
  double pzD2GenZ  ;
  double enD2GenZ  ;
  
  double pxPatZ  ;
  double pyPatZ  ;
  double pzPatZ  ;
  double enPatZ  ;
  double pxD1PatZ  ;
  double pyD1PatZ  ;
  double pzD1PatZ  ;
  double enD1PatZ  ;
  int    passD1PatZ  ;
  double pxD2PatZ  ;
  double pyD2PatZ  ;
  double pzD2PatZ  ;
  double enD2PatZ  ;
  int    passD2PatZ  ;
  
  std::vector<double> genJetPx;
  std::vector<double> genJetPy;
  std::vector<double> genJetPz;
  std::vector<double> genJetEn;
  std::vector<double> patJetPx;
  std::vector<double> patJetPy;
  std::vector<double> patJetPz;
  std::vector<double> patJetEn;
};

#endif
