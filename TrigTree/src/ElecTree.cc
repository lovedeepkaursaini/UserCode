#include <iostream>
#include "UserCode/TrigTree/interface/ElecTree.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/Common/interface/Handle.h"

// for conversion finder
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"


ElecTree::ElecTree(std::string name, TTree* tree,const edm::ParameterSet & iConfig)
{
  nElec_=0;
  name_=name;
  tree_=tree;
  SetBranches();
  
  recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEcalRecHitsEB"));
  recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEcalRecHitsEE"));
  
  dcsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("DcsTag",edm::InputTag("scalersRawToDigi"));
  
  isData_ = iConfig.getUntrackedParameter<bool>("IsData");
  
  std::cout << "Data: " << isData_ << std::endl;
}


ElecTree::~ElecTree(){
  delete tree_;
}


void
ElecTree::ElecInfoCollector(const edm::Event& iEvent,const edm::EventSetup& iSetup,const reco::GsfElectron& ele,edm::Handle<reco::BeamSpot> beamSpotHandle){

  nElec_=nElec_+1;
  elept_.push_back(ele.pt());
  eleeta_.push_back(ele.eta());
  elephi_.push_back(ele.phi());
  eleet_.push_back(ele.et());
  elep_.push_back(ele.p()) ;
  elepx_.push_back(ele.px()) ;
  elepy_.push_back(ele.py()) ;
  elepz_.push_back(ele.pz()) ;
  eleenergy_.push_back(ele.energy());
  
  EcalClusterLazyTools lazyTools(iEvent, iSetup, recHitsEBTag_, recHitsEETag_);

  reco::GsfTrackRef track = ele.gsfTrack();
  reco::SuperClusterRef supercluster = ele.superCluster();

  scenergy_.push_back(supercluster->energy());
  scet_.push_back(supercluster->energy()/cosh(supercluster->eta()));
  scpreshowerenergy_.push_back(supercluster->preshowerEnergy());
  scetawidth_.push_back(supercluster->etaWidth());
  scphiwidth_.push_back(supercluster->phiWidth());
  sceta_.push_back(supercluster->eta());
  scphi_.push_back(supercluster->phi());
  
  //std::cout<<"(ele.et()):  "<<(ele.et())<<"(supercluster->energy()/cosh(supercluster->eta())):  "<<(supercluster->energy()/cosh(supercluster->eta()))<<std::endl;
  const reco::CaloClusterPtr bcref = supercluster->seed();
  const reco::BasicCluster *bc = bcref.get();
  
  float ee3x3 = lazyTools.e3x3(*bc);
  float eeMax = lazyTools.eMax(*bc);
  float rmax3x3=eeMax/ee3x3; 
  elermax3x3_.push_back(rmax3x3);

  //divide det. region
  if (fabs(ele.eta()) < 1.4442 )
    {
      eleisBarrel_.push_back(1);
      eleisEndcap_.push_back(0);
    }
  else if (fabs(ele.eta())> 1.566 && fabs(ele.eta())< 2.5)
    {
      eleisBarrel_.push_back(0);
      eleisEndcap_.push_back(1);
    }
  else 
    {
      eleisBarrel_.push_back(0);
      eleisEndcap_.push_back(0);
    } 
  
  //see if e is ecal/tracker driven
  if(ele.ecalDrivenSeed())
    {
      eleisEcalDriven_.push_back(1);
      eleisTrackerDriven_.push_back(0);
    }
  else if(ele.trackerDrivenSeed())
    {
      eleisEcalDriven_.push_back(0);
      eleisTrackerDriven_.push_back(1);
    }
  else
    {
      eleisEcalDriven_.push_back(0);
      eleisTrackerDriven_.push_back(0);
    }
  eleecalenergy_.push_back(ele.ecalEnergy());
  eleecalenergyerror_.push_back(ele.ecalEnergyError());



  dhElecClsTrkCalo_.push_back(ele.deltaEtaEleClusterTrackAtCalo());
  dhSeedClsTrkCalo_.push_back(ele.deltaEtaSeedClusterTrackAtCalo());
  dhSuperClsTrkVtx_.push_back(ele.deltaEtaSuperClusterTrackAtVtx());
 
  dPhiElecClsTrkCalo_.push_back(ele.deltaPhiEleClusterTrackAtCalo());
  dPhiSeedClsTrkCalo_.push_back(ele.deltaPhiSeedClusterTrackAtCalo());
  dPhiSuperClsTrkVtx_.push_back(ele.deltaPhiSuperClusterTrackAtVtx());
 
  trkPosXVtx_.push_back(ele.trackPositionAtVtx().X());
  trkPosYVtx_.push_back(ele.trackPositionAtVtx().Y());
  trkPosZVtx_.push_back(ele.trackPositionAtVtx().Z());

  trkMomXVtx_.push_back(ele.trackMomentumAtVtx().X());
  trkMomYVtx_.push_back(ele.trackMomentumAtVtx().Y());
  trkMomZVtx_.push_back(ele.trackMomentumAtVtx().Z());

  trkPosXCalo_.push_back(ele.trackPositionAtCalo().X());
  trkPosYCalo_.push_back(ele.trackPositionAtCalo().Y());
  trkPosZCalo_.push_back(ele.trackPositionAtCalo().Z());

  trkMomXCalo_.push_back(ele.trackMomentumAtCalo().X());
  trkMomYCalo_.push_back(ele.trackMomentumAtCalo().Y());
  trkMomZCalo_.push_back(ele.trackMomentumAtCalo().Z());

  eEleClsOverPout_.push_back(ele.eEleClusterOverPout());
  eSeedClsOverP_.push_back(ele.eSeedClusterOverP());
  eSeedClsOverPout_.push_back(ele.eSeedClusterOverPout());
  eSuperClsOverP_.push_back(ele.eSuperClusterOverP());
  
  eleMomErr_.push_back(ele.electronMomentumError());

  // isolation dR 0.3     
  eledr03EcalRecHitSumEt_.push_back(ele.dr03EcalRecHitSumEt());
  eledr03HcalDepth1TowerSumEt_.push_back(ele.dr03HcalDepth1TowerSumEt());
  eledr03HcalDepth2TowerSumEt_.push_back(ele.dr03HcalDepth2TowerSumEt());
  eledr03HcalTowerSumEt_.push_back(ele.dr03HcalTowerSumEt());
  eledr03TkSumPt_.push_back(ele.dr03TkSumPt());
  
  // isolation dR 0.4
  eledr04EcalRecHitSumEt_.push_back(ele.dr04EcalRecHitSumEt());
  eledr04HcalDepth1TowerSumEt_.push_back(ele.dr04HcalDepth1TowerSumEt());
  eledr04HcalDepth2TowerSumEt_.push_back(ele.dr04HcalDepth2TowerSumEt());
  eledr04HcalTowerSumEt_.push_back(ele.dr04HcalTowerSumEt());
  eledr04TkSumPt_.push_back(ele.dr04TkSumPt());

  // relative isolations
  eleRelIsoEcal_.push_back(ele.dr03EcalRecHitSumEt()/ele.p4().Pt());
  eleRelIsoHcal_.push_back(ele.dr03HcalTowerSumEt()/ele.p4().Pt());
  eleRelIsoTrk_.push_back(ele.dr03TkSumPt()/ele.p4().Pt());
  if ( fabs(ele.eta()) < 1.479 ) // barrel
    eleRelIsoComb_.push_back( ( ele.dr03TkSumPt() + std::max(0., ele.dr03EcalRecHitSumEt() - 1.) + ele.dr03HcalTowerSumEt() ) / ele.p4().Pt() ) ;
  else
    eleRelIsoComb_.push_back( ( ele.dr03TkSumPt() + ele.dr03EcalRecHitSumEt() + ele.dr03HcalTowerSumEt() ) / ele.p4().Pt() ) ;
  
  ConversionRejection(iEvent, iSetup, ele);


  e1x5_.push_back(ele.e1x5());
  e2x5Max_.push_back(ele.e2x5Max()) ;
  e5x5_.push_back(ele.e5x5()) ;
  eler1x5_.push_back(ele.e1x5()/ele.e5x5());
  eler2x5max_.push_back(ele.e2x5Max()/ele.e5x5());

  elescr9_.push_back(ee3x3/ele.energy()); // energy in 3x3 around seed cluster / SC energy
  elesceseedoveresupercluster_.push_back(eeMax/ele.energy()); // energy of max crystal / SC energy
  
  elehcalDepth1OverEcal_.push_back(ele.hcalDepth1OverEcal());
  elehcalDepth2OverEcal_.push_back(ele.hcalDepth2OverEcal());
  elehcalOverEcal_.push_back(ele.hcalOverEcal());
  
  // sigma eta eta, sigma ieta ieta
  elesigmaEtaEta_.push_back(ele.sigmaEtaEta());
  elesigmaIetaIeta_.push_back(ele.sigmaIetaIeta());
  
  // misc + charge
  elebasicClustersSize_.push_back(ele.basicClustersSize()) ;
  elenumberOfBrems_.push_back(ele.numberOfBrems());
  elefbrem_.push_back(ele.fbrem()) ;
  elescPixCharge_.push_back(ele.scPixCharge()) ;
  if (ele.closestCtfTrackRef().isNonnull())
    {
      electfcharge_.push_back(ele.closestCtfTrackRef()->charge());
    }
  else electfcharge_.push_back(0.);
  elecharge_.push_back(ele.charge());
  eleisGsfScPixChargeConsistent_.push_back(ele.isGsfScPixChargeConsistent());
  eleisGsfCtfChargeConsistent_.push_back(ele.isGsfCtfChargeConsistent());
  eleisGsfCtfScPixChargeConsistent_.push_back(ele.isGsfCtfScPixChargeConsistent());
  
  // vertex
  elevertexChi2_.push_back(ele.vertexChi2());
  elevertexNdof_.push_back(ele.vertexNdof());
  elevertexNormalizedChi2_.push_back(ele.vertexNormalizedChi2());
  elevx_.push_back(ele.vx());
  elevy_.push_back(ele.vy());
  elevz_.push_back(ele.vz());
  elevertexX_.push_back(ele.vertex().x());
  elevertexY_.push_back(ele.vertex().y());
  elevertexZ_.push_back(ele.vertex().z());
 // vertex TIP
  double vertexTIP =
    (ele.vertex().x()-beamSpotHandle->position().x()) * (ele.vertex().x()-beamSpotHandle->position().x()) +
    (ele.vertex().y()-beamSpotHandle->position().y()) * (ele.vertex().y()-beamSpotHandle->position().y()) ;
  vertexTIP = sqrt(vertexTIP) ;
  elevertexTIP_.push_back(vertexTIP);
  
  // track related distributions
  eledxy_.push_back(track->dxy(beamSpotHandle->position()));
  elegsfcharge_.push_back(track->charge());
  eleambiguousTracks_.push_back(ele.ambiguousGsfTracksSize());
  elefoundHits_.push_back(track->numberOfValidHits() );
  elelostHits_.push_back(track->numberOfLostHits() );
  elechi2_.push_back(track->normalizedChi2() );
}



void
ElecTree::ConversionRejection(const edm::Event& iEvent, const edm::EventSetup& iSetup, const reco::GsfElectron& ele)
{
  double evt_bField;
  
  edm::Handle<reco::TrackCollection> tracks_h;
  iEvent.getByLabel("generalTracks", tracks_h);
  
  if (isData_) 
  {
    edm::Handle<DcsStatusCollection> dcsHandle;
    iEvent.getByLabel(dcsTag_, dcsHandle);
    
    // scale factor = 3.801/18166.0 which are average values taken over a stable two week period
    float currentToBFieldScaleFactor = 2.09237036221512717e-04;
    float current = (*dcsHandle)[0].magnetCurrent();
    evt_bField = current*currentToBFieldScaleFactor;
  }
  else 
  {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
        
    evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
  }
  
  ConversionFinder convFinder;
  ConversionInfo convInfo = convFinder.getConversionInfo(ele, tracks_h, evt_bField);

  eleMissingHits_.push_back( ele.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits() );
  
  eleDist_.push_back(convInfo.dist() );
  eleDeltaCotTheta_.push_back(convInfo.dcot() ) ;
  eleConvRadius_.push_back(convInfo.radiusOfConversion() ) ;
}


void
ElecTree::SetBranches(){

  AddBranch(&nElec_,"Num");
  AddBranch(&elept_,  "Pt");
  AddBranch(&eleeta_, "Eta");
  AddBranch(&elephi_, "Phi");
  AddBranch(&eleet_,"Et");
  AddBranch(&scet_,"scEt");
  AddBranch(&sceta_,"scEta");
  AddBranch(&scphi_,"scphi");
  AddBranch(&scenergy_,"scenergy");

  AddBranch(&elermax3x3_,"elermax3x3");
  AddBranch(&eleisBarrel_,"eleisBarrel");
  AddBranch(&eleisEndcap_,"eleisEndcap");
  AddBranch(&eleisEcalDriven_,"eleisEcalDriven");
  AddBranch(&eleisTrackerDriven_,"eleisTrackerDriven");
  AddBranch(&eleecalenergy_,"eleecalenergy");
  AddBranch(&eleecalenergyerror_,"eleecalenergyerror");
  AddBranch(&elescr9_,"elescr9");
  AddBranch(&elesceseedoveresupercluster_,"elesceseedoveresupercluster");
  
  AddBranch(&dhElecClsTrkCalo_,"DhElecClsTrkAtCalo");
  AddBranch(&dhSeedClsTrkCalo_,"DhSeedClsTrkAtCalo");
  AddBranch(&dhSuperClsTrkVtx_,"DhSuperClsTrkAtVtx");
 
  AddBranch(&dPhiElecClsTrkCalo_,"DphiElecClsTrkAtCalo");
  AddBranch(&dPhiSeedClsTrkCalo_,"DphiSeedClsTrkAtCalo");
  AddBranch(&dPhiSuperClsTrkVtx_,"DphiSuperClsTrkAtVtx");
 
  AddBranch(&trkPosXVtx_,"PositionXTrkAtVtx");
  AddBranch(&trkPosYVtx_,"PositionYTrkAtVtx");
  AddBranch(&trkPosZVtx_,"PositionZTrkAtVtx");

  AddBranch(&trkMomXVtx_,"MomentumXTrkAtVtx");
  AddBranch(&trkMomYVtx_,"MomentumYTrkAtVtx");
  AddBranch(&trkMomZVtx_,"MomentumZTrkAtVtx");

  AddBranch(&trkPosXCalo_,"PositionXTrkAtCalo");
  AddBranch(&trkPosYCalo_,"PositionYTrkAtCalo");
  AddBranch(&trkPosZCalo_,"PositionZTrkAtCalo");

  AddBranch(&trkMomXCalo_,"MomentumXTrkAtCalo");
  AddBranch(&trkMomYCalo_,"MomentumYTrkAtCalo");
  AddBranch(&trkMomZCalo_,"MomentumZTrkAtCalo");

  AddBranch(&eEleClsOverPout_,"eEleClsOverPout");
  AddBranch(&eSeedClsOverP_,"eSeedClsOverP");
  AddBranch(&eSeedClsOverPout_,"eSeedClsOverPout");
  AddBranch(&eSuperClsOverP_,"eSuperClsOverP");
  
  AddBranch(&eleMomErr_,"eleMomErr");

  // isolation dR 0.3     
  AddBranch(&eledr03EcalRecHitSumEt_,"eledr03EcalRecHitSumEt");
  AddBranch(&eledr03HcalDepth1TowerSumEt_,"eledr03HcalDepth1TowerSumEt");
  AddBranch(&eledr03HcalDepth2TowerSumEt_,"eledr03HcalDepth2TowerSumEt");
  AddBranch(&eledr03HcalTowerSumEt_,"eledr03HcalTowerSumEt");
  AddBranch(&eledr03TkSumPt_,"eledr03TkSumPt");
  
  // isolation dR 0.4
  AddBranch(&eledr04EcalRecHitSumEt_,"eledr04EcalRecHitSumEt");
  AddBranch(&eledr04HcalDepth1TowerSumEt_,"eledr04HcalDepth1TowerSumEt");
  AddBranch(&eledr04HcalDepth2TowerSumEt_,"eledr04HcalDepth2TowerSumEt");
  AddBranch(&eledr04HcalTowerSumEt_,"eledr04HcalTowerSumEt");
  AddBranch(&eledr04TkSumPt_,"eledr04TkSumPt");

  // relative isolations
  AddBranch(&eleRelIsoEcal_, "eleRelIsoEcal");
  AddBranch(&eleRelIsoHcal_, "eleRelIsoHcal");
  AddBranch(&eleRelIsoTrk_,  "eleRelIsoTrk");
  AddBranch(&eleRelIsoComb_,  "eleRelIsoComb");
  
  // conversion rejection
  AddBranch(&eleMissingHits_, "eleMissingHits");
  AddBranch(&eleDist_, "eleDist" );
  AddBranch(&eleDeltaCotTheta_, "eleDeltaCotTheta");
  AddBranch(&eleConvRadius_, "eleConvRadius");


  AddBranch(&e1x5_,"e1x5");
  AddBranch(&e2x5Max_,"e2x5Max");
  AddBranch(&e5x5_,"e5x5");
  AddBranch(&eler1x5_,"eler1x5");
  AddBranch(&eler2x5max_,"eler2x5max");
  
  AddBranch(&scpreshowerenergy_,"scpreshowerenergy");
  AddBranch(&scetawidth_,"scetawidth");
  AddBranch(&scphiwidth_,"scphiwidth");
  AddBranch(&eleenergy_,"eleenergy");
  
  AddBranch(&elehcalDepth1OverEcal_,"elehcalDepth1OverEcal");
  AddBranch(&elehcalDepth2OverEcal_,"elehcalDepth2OverEcal");
  AddBranch(&elehcalOverEcal_,"elehcalOverEcal");

  AddBranch(&elesigmaEtaEta_,"elesigmaEtaEta");
  AddBranch(&elesigmaIetaIeta_,"elesigmaIetaIeta");
  
  AddBranch(&elebasicClustersSize_,"elebasicClustersSize");
  AddBranch(&elenumberOfBrems_,"elenumberOfBrems");
  AddBranch(&elefbrem_,"elefbrem");
  AddBranch(&elescPixCharge_,"elescPixCharge");
  AddBranch(&electfcharge_,"electfcharge");

  AddBranch(&elecharge_,"elecharge");
  AddBranch(&eleisGsfScPixChargeConsistent_,"eleisGsfScPixChargeConsistent");
  AddBranch(&eleisGsfCtfChargeConsistent_,"eleisGsfCtfChargeConsistent");
  AddBranch(&eleisGsfCtfScPixChargeConsistent_,"eleisGsfCtfScPixChargeConsistent");

  AddBranch(&elevertexChi2_,"elevertexChi2");
  AddBranch(&elevertexNdof_,"elevertexNdof");
  AddBranch(&elevertexNormalizedChi2_,"elevertexNormalizedChi2");
  AddBranch(&elevx_,"elevx");
  AddBranch(&elevy_,"elevy");
  AddBranch(&elevz_,"elevz");
  AddBranch(&elevertexX_,"elevertexX");
  AddBranch(&elevertexY_,"elevertexY");
  AddBranch(&elevertexZ_,"elevertexZ");
  AddBranch(&elevertexTIP_,"elevertexTIP");
  
  AddBranch(&elep_,"elep");
  AddBranch(&elepx_,"elepx");
  AddBranch(&elepy_,"elepy");
  AddBranch(&elepz_,"elepz");

  AddBranch(&eledxy_,"eledxy");
  AddBranch(&elegsfcharge_,"elegsfcharge");
  AddBranch(&eleambiguousTracks_,"eleambiguousTracks");
  AddBranch(&elefoundHits_,"elefoundHits");
  AddBranch(&elelostHits_,"elelostHits");
  AddBranch(&elechi2_,"elechi2");
}

void
ElecTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName="Electron"+name;
  tree_->Branch(brName.c_str(),vec);
}

void 
ElecTree::AddBranch(double* x, std::string name){
  std::string brName="Electron"+name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}

void 
ElecTree::Clear(){
  nElec_=0;

  elept_.clear();
  eleeta_.clear();
  elephi_.clear();
  eleet_.clear();
  scet_.clear();
  sceta_.clear();
  scphi_.clear();
  scenergy_.clear();

  elermax3x3_.clear();
  eleisBarrel_.clear();
  eleisEndcap_.clear();
  eleisEcalDriven_.clear();
  eleisTrackerDriven_.clear();
  eleecalenergy_.clear();
  eleecalenergyerror_.clear();
  elescr9_.clear();
  elesceseedoveresupercluster_.clear();
  
  dhElecClsTrkCalo_.clear();
  dhSeedClsTrkCalo_.clear();
  dhSuperClsTrkVtx_.clear();
  dPhiElecClsTrkCalo_.clear();
  dPhiSeedClsTrkCalo_.clear();
  dPhiSuperClsTrkVtx_.clear();
  trkPosXVtx_.clear();
  trkPosYVtx_.clear();
  trkPosZVtx_.clear();
  trkMomXVtx_.clear();
  trkMomYVtx_.clear();
  trkMomZVtx_.clear();
  trkPosXCalo_.clear();
  trkPosYCalo_.clear();
  trkPosZCalo_.clear();
  trkMomXCalo_.clear();
  trkMomYCalo_.clear();
  trkMomZCalo_.clear();

  eEleClsOverPout_.clear();
  eSeedClsOverP_.clear();
  eSeedClsOverPout_.clear();
  eSuperClsOverP_.clear();
  
  eleMomErr_.clear();

  // isolation dR 0.3     
  eledr03EcalRecHitSumEt_.clear();
  eledr03HcalDepth1TowerSumEt_.clear();
  eledr03HcalDepth2TowerSumEt_.clear();
  eledr03HcalTowerSumEt_.clear();
  eledr03TkSumPt_.clear();
  
  // isolation dR 0.4
  eledr04EcalRecHitSumEt_.clear();
  eledr04HcalDepth1TowerSumEt_.clear();
  eledr04HcalDepth2TowerSumEt_.clear();
  eledr04HcalTowerSumEt_.clear();
  eledr04TkSumPt_.clear();
 
  // relative isolations
  eleRelIsoComb_.clear();
  eleRelIsoEcal_.clear();
  eleRelIsoHcal_.clear();
  eleRelIsoTrk_.clear();
  
  // conversion rejection
  eleMissingHits_.clear();
  eleDist_.clear();
  eleDeltaCotTheta_.clear();
  eleConvRadius_.clear();
 
  e1x5_.clear();
  e2x5Max_.clear();
  e5x5_.clear();
  eler1x5_.clear();
  eler2x5max_.clear();
  
  scpreshowerenergy_.clear();
  scetawidth_.clear();
  scphiwidth_.clear();
  eleenergy_.clear();
  
  elehcalDepth1OverEcal_.clear();
  elehcalDepth2OverEcal_.clear();
  elehcalOverEcal_.clear();
  
  elesigmaEtaEta_.clear();
  elesigmaIetaIeta_.clear();
  
  elebasicClustersSize_.clear();
  elenumberOfBrems_.clear();
  elefbrem_.clear();
  elescPixCharge_.clear();
  electfcharge_.clear();
  
  elecharge_.clear();
  eleisGsfScPixChargeConsistent_.clear();
  eleisGsfCtfChargeConsistent_.clear();
  eleisGsfCtfScPixChargeConsistent_.clear();
  
  elevertexChi2_.clear();
  elevertexNdof_.clear();
  elevertexNormalizedChi2_.clear();
  elevx_.clear();
  elevy_.clear();
  elevz_.clear();
  elevertexX_.clear();
  elevertexY_.clear();
  elevertexZ_.clear();
  elevertexTIP_.clear();
 
  elep_.clear();
  elepx_.clear();
  elepy_.clear();
  elepz_.clear();
  
  eledxy_.clear();
  elegsfcharge_.clear();
  eleambiguousTracks_.clear();
  elefoundHits_.clear();
  elelostHits_.clear();
  elechi2_.clear();
}
