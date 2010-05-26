#ifndef __CRAZY_PUNJABI4__
#define __CRAZY_PUNJABI4__

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h" 
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

class ElecTree{

 public:
  ElecTree(std::string name, TTree* tree
	   /*,const edm::ParameterSet iConfig*/);
  
  ~ElecTree();
  void ElecInfoCollector(const reco::GsfElectron& elec); 
  void Clear();
 private:
  ElecTree(){};//Don't allow user
  void SetBranches();
  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  //void AddBranch(bool* x, std::string name);
  //void AddBranch(double* x, std::string name);
  
  
  TTree *tree_;
  std::string name_;
  double nElec_;
  std::vector<double> elept_;
  std::vector<double> eleeta_;
  std::vector<double> elephi_;
  std::vector<double> eleet_;

  std::vector<double> dhElecClsTrkCalo_;
  std::vector<double> dhSeedClsTrkCalo_;
  std::vector<double> dhSuperClsTrkVtx_;
  
  std::vector<double> dPhiElecClsTrkCalo_;
  std::vector<double> dPhiSeedClsTrkCalo_;
  std::vector<double> dPhiSuperClsTrkVtx_;
  
  std::vector<double> trkPosXVtx_;
  std::vector<double> trkPosYVtx_;
  std::vector<double> trkPosZVtx_;
  
  std::vector<double> trkMomXVtx_;
  std::vector<double> trkMomYVtx_;
  std::vector<double> trkMomZVtx_;
  
  std::vector<double> trkPosXCalo_;
  std::vector<double> trkPosYCalo_;
  std::vector<double> trkPosZCalo_;
  
  std::vector<double> trkMomXCalo_;
  std::vector<double> trkMomYCalo_;
  std::vector<double> trkMomZCalo_;

  std::vector<double> eEleClsOverPout_;
  std::vector<double> eSeedClsOverP_;
  std::vector<double> eSeedClsOverPout_;
  std::vector<double> eSuperClsOverP_;
  
  std::vector<double> eleMomErr_;

  // isolation dR 0.3     
  std::vector<double> eledr03EcalRecHitSumEt_;
  std::vector<double> eledr03HcalDepth1TowerSumEt_;
  std::vector<double> eledr03HcalDepth2TowerSumEt_;
  std::vector<double> eledr03HcalTowerSumEt_;
  std::vector<double> eledr03TkSumPt_;
  
  // isolation dR 0.4
  std::vector<double> eledr04EcalRecHitSumEt_;
  std::vector<double> eledr04HcalDepth1TowerSumEt_;
  std::vector<double> eledr04HcalDepth2TowerSumEt_;
  std::vector<double> eledr04HcalTowerSumEt_;
  std::vector<double> eledr04TkSumPt_;
  
  std::vector<double> e1x5_;
  std::vector<double> e2x5Max_;
  std::vector<double> e5x5_;
  std::vector<double> eler1x5_;
  std::vector<double> eler2x5max_;
  
  std::vector<double> scpreshowerenergy_;
  std::vector<double> scetawidth_;
  std::vector<double> scphiwidth_;
  std::vector<double> eleenergy_;
  
  std::vector<double> elehcalDepth1OverEcal_;
  std::vector<double> elehcalDepth2OverEcal_;
  std::vector<double> elehcalOverEcal_;

  std::vector<double> elesigmaEtaEta_;
  std::vector<double> elesigmaIetaIeta_;
  
  std::vector<double> elebasicClustersSize_;
  std::vector<double> elenumberOfBrems_;
  std::vector<double> elefbrem_;
  std::vector<double> elescPixCharge_;
  std::vector<double> electfcharge_;

  std::vector<double> elecharge_;
  std::vector<double> eleisGsfScPixChargeConsistent_;
  std::vector<double> eleisGsfCtfChargeConsistent_;
  std::vector<double> eleisGsfCtfScPixChargeConsistent_;

  std::vector<double> elevertexChi2_;
  std::vector<double> elevertexNdof_;
  std::vector<double> elevertexNormalizedChi2_;
  std::vector<double> elevx_;
  std::vector<double> elevy_;
  std::vector<double> elevz_;
  
  std::vector<double> elep_;
  std::vector<double> elepx_;
  std::vector<double> elepy_;
  std::vector<double> elepz_;
};

#endif
