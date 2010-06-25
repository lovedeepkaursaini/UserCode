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
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"


class ElecTree{

 public:
  ElecTree(std::string name, TTree* tree,const edm::ParameterSet& cfg
	   /*,const edm::ParameterSet iConfig*/);
  
  ~ElecTree();
  void ElecInfoCollector(const edm::Event& iEvent,const edm::EventSetup& iSetup,const reco::GsfElectron& elec,edm::Handle<reco::BeamSpot> BSH); 
  void Clear();
 private:
  ElecTree(){};//Don't allow user
  void SetBranches();
  void ConversionRejection(const edm::Event&, const edm::EventSetup&, const reco::GsfElectron& elec);
  void AddBranch(double* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<int>*, std::string name);
  //void AddBranch(bool* x, std::string name);
  //void AddBranch(double* x, std::string name);
  
  edm::InputTag recHitsEBTag_, recHitsEETag_;
  edm::InputTag dcsTag_;
  bool isData_;
  
  TTree *tree_;
  std::string name_;
  double nElec_;
  const reco::BasicCluster *bc;
  std::vector<double> elept_;
  std::vector<double> eleeta_;
  std::vector<double> elephi_;
  std::vector<double> eleet_;
  std::vector<double> scet_;
  std::vector<double> sceta_;
  std::vector<double> scphi_;
  std::vector<double> scenergy_;

  std::vector<double> elermax3x3_;
  std::vector<double> eleisBarrel_;
  std::vector<double> eleisEndcap_;
  std::vector<double> eleisEcalDriven_;
  std::vector<double> eleisTrackerDriven_;
  std::vector<double> eleecalenergy_;
  std::vector<double> eleecalenergyerror_;
  std::vector<double> elescr9_;
  std::vector<double> elesceseedoveresupercluster_;
  std::vector<double> eledxy_;
  std::vector<double> elegsfcharge_;

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
  
  // relative isolations
  std::vector<double> eleRelIsoEcal_ ;
  std::vector<double> eleRelIsoHcal_ ;
  std::vector<double> eleRelIsoTrk_ ;
  std::vector<double> eleRelIsoComb_ ;
  
  // conversion rejection
  std::vector<double> eleMissingHits_ ;
  std::vector<double> eleDist_ ;
  std::vector<double> eleDeltaCotTheta_ ;
  std::vector<double> eleConvRadius_ ;

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
  std::vector<double> elevertexX_;
  std::vector<double> elevertexY_;
  std::vector<double> elevertexZ_;
  std::vector<double> elevertexTIP_;
  
  std::vector<double> elep_;
  std::vector<double> elepx_;
  std::vector<double> elepy_;
  std::vector<double> elepz_;

  std::vector<double> eleambiguousTracks_;
  std::vector<double> elefoundHits_;
  std::vector<double> elelostHits_;
  std::vector<double> elechi2_;
};

#endif
