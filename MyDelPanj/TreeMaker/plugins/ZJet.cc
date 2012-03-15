// -*- C++ -*-
//
// Package:    ZJet
// Class:      ZJet
// 
/**\class ZJet ZJet.cc DelPanj/ZJet/src/ZJet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lovedeep Kaur (Panjab U)
//         Created:  Fri Dec 30 09:58:44 CST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h" 
#include "DelPanj/TreeMaker/interface/eSelector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/eHist.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


//
// class declaration
//

namespace zJet{
  TH1D* prettyHistogram(edm::Service<TFileService> fs, std::string name, std::string title, std::string xtitle, double nbins, double rlow, double rhigh ){
    TH1D* hist = fs->make<TH1D>(name.c_str(),title.c_str(),nbins,rlow,rhigh);
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    return hist;
  }


  TH1D* prettyHistogram(edm::Service<TFileService> fs, std::string name, std::string title, std::string xtitle, double nbins, float* xbins){
    TH1D* hist = fs->make<TH1D>(name.c_str(),title.c_str(),nbins,xbins);
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    return hist;
  }

  double Phi_0_2pi(double x) {
    while (x >= 2*M_PI) x -= 2*M_PI;
    while (x <     0.)  x += 2*M_PI;
    return x;
  }

  double radius(double eta1, double phi1,double eta2, double phi2){

    const double TWOPI= 2.0*M_PI;

    phi1=Phi_0_2pi(phi1);
    phi2=Phi_0_2pi(phi2);

    double dphi=Phi_0_2pi(phi1-phi2);
    dphi = TMath::Min(dphi,TWOPI-dphi);
    double deta = eta1-eta2;

    return sqrt(deta*deta+dphi*dphi);
  }

  double max(double one, double two){ if (one<two)return two;else return one;}
}

class ZJet : public edm::EDAnalyzer {
   public:
      explicit ZJet(const edm::ParameterSet&);
      ~ZJet();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);


  TH1D* z0Hist_;
  TH1D* z1Hist_;
  TH1D* z2Hist_;
  TH1D* z3Hist_;
  TH1D* z4Hist_;
  TH1D* z5Hist_;
  TH1D* z6Hist_;
  TH1D* z7Hist_;
  TH1D* z8Hist_;
  TH1D* z9Hist_;
  TH1D* z10Hist_;
  TH1D* z11Hist_;

  TH1D* hNjetGoodZX; 
  TH1D* hEta1 ; 
  TH1D* FirstJetPt_Z1jet;
  TH1D* SecondJetPt_Z2jet;
  TH1D* ThirdJetPt_Z3jet;
  TH1D* FourthJetPt_Z4jet;
  TH1D* FifthJetPt_Z5jet;
  TH1D* SixthJetPt_Z6jet;
  
  TH1D* hEffi_ ; 
  TH1D* hNjetGoodIncZX ; 
  TH1D* ZMass_Zinc0jet;
  TH1D* ZMass_Zinc1jet;
  TH1D* ZMass_Zinc2jet;
  TH1D* ZMass_Zinc3jet;
  TH1D* ZMass_Zinc4jet;
  TH1D* ZMass_Zinc5jet;
  TH1D* ZMass_Zinc6jet;
  
  TH1D* ZMass_Zexc0jet;
  TH1D* ZMass_Zexc1jet;
  TH1D* ZMass_Zexc2jet;
  TH1D* ZMass_Zexc3jet;
  TH1D* ZMass_Zexc4jet;
  TH1D* ZMass_Zexc5jet;
  TH1D* ZMass_Zexc6jet;
  
  TH1D* ZPt_Zinc0jet;
  TH1D* ZPt_Zinc1jet;
  TH1D* ZPt_Zinc2jet;
  TH1D* ZPt_Zinc3jet;
  TH1D* ZPt_Zinc4jet;
  TH1D* ZPt_Zinc5jet;
  TH1D* ZPt_Zinc6jet;
  
  TH1D* ZPt_Zexc0jet;
  TH1D* ZPt_Zexc1jet;
  TH1D* ZPt_Zexc2jet;
  TH1D* ZPt_Zexc3jet;
  TH1D* ZPt_Zexc4jet;
  TH1D* ZPt_Zexc5jet;
  TH1D* ZPt_Zexc6jet;
  
  TH1D* ZRapidity_Zinc0jet;
  TH1D* ZRapidity_Zinc1jet;
  TH1D* ZRapidity_Zinc2jet;
  TH1D* ZRapidity_Zinc3jet;
  TH1D* ZRapidity_Zinc4jet;
  TH1D* ZRapidity_Zinc5jet;
  TH1D* ZRapidity_Zinc6jet;
  
  
  TH1D* ZRapidity_Zexc0jet;
  TH1D* ZRapidity_Zexc1jet;
  TH1D* ZRapidity_Zexc2jet;
  TH1D* ZRapidity_Zexc3jet;
  TH1D* ZRapidity_Zexc4jet;
  TH1D* ZRapidity_Zexc5jet;
  TH1D* ZRapidity_Zexc6jet;
  
  
  TH1D* FirstJetY_Z1jet;
  TH1D* SecondJetY_Z2jet;
  TH1D* ThirdJetY_Z3jet;
  TH1D* FourthJetY_Z4jet;
  TH1D* FifthJetY_Z5jet;
  TH1D* SixthJetY_Z6jet;
  
  TH1D* hAccept ; 
  TH1D* hCounter;
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZJet::ZJet(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  z0Hist_ = fs->make<TH1D>("hZ0Mass","hZ0Mass",50,0,200);
  z1Hist_ = fs->make<TH1D>("hZ1Mass","hZ1Mass",50,0,200);
  z2Hist_ = fs->make<TH1D>("hZ2Mass","hZ2Mass",50,0,200);
  z3Hist_ = fs->make<TH1D>("hZ3Mass","hZ3Mass",50,0,200);
  z4Hist_ = fs->make<TH1D>("hZ4Mass","hZ4Mass",50,0,200);
  z5Hist_ = fs->make<TH1D>("hZ5Mass","hZ5Mass",50,0,200);

  hNjetGoodZX = zJet::prettyHistogram(fs,"h_Njet_GoodZX","N_{Jets}(Exc) [Z+Jets]","NJets_{Exc}",10,0,10);
  
  float jpt[10] = {30,50,70,90,120,150,180,210,300,400};
  float jpt2[9] = {30,50,70,90,120,150,200,270,400};
  float jpt3[6] = {30,50,70,90,120,400};
  float jpt4[6] = {30,50,70,90,120,400};
  
  hEta1 = new TH1D("hEta1","hEta2",30,-3,3);

  FirstJetPt_Z1jet=zJet::prettyHistogram(fs,"FirstJetPt_Z1jet","First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
  SecondJetPt_Z2jet=zJet::prettyHistogram(fs,"SecondJetPt_Z2jet","Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
  ThirdJetPt_Z3jet=zJet::prettyHistogram(fs,"ThirdJetPt_Z3jet","Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  FourthJetPt_Z4jet=zJet::prettyHistogram(fs,"FourthJetPt_Z4jet","Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt4);
  FifthJetPt_Z5jet=zJet::prettyHistogram(fs,"FifthJetPt_Z5jet","Fifth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",9,jpt);
  SixthJetPt_Z6jet=zJet::prettyHistogram(fs,"SixthJetPt_Z6jet","#geq 6th Jets (Z+Jets)","p_{T}(Jet) [GeV]",9,jpt);
  
  hEffi_ = zJet::prettyHistogram(fs,"hEffi_","Weight","Weight",50,0,1);

  hNjetGoodIncZX = zJet::prettyHistogram(fs,"h_Njet_GoodIncZX","N_{Jets}(Exc) [Z+Jets]","NJets_{Inc}",10,0,10);
  ZMass_Zinc0jet=zJet::prettyHistogram(fs,"ZMass_Zinc0jet","Z Invariant Mass (N_{jets} #geq 0)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zinc1jet=zJet::prettyHistogram(fs,"ZMass_Zinc1jet","Z Invariant Mass (N_{jets} #geq 1)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zinc2jet=zJet::prettyHistogram(fs,"ZMass_Zinc2jet","Z Invariant Mass (N_{jets} #geq 2)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zinc3jet=zJet::prettyHistogram(fs,"ZMass_Zinc3jet","Z Invariant Mass (N_{jets} #geq 3)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zinc4jet=zJet::prettyHistogram(fs,"ZMass_Zinc4jet","Z Invariant Mass (N_{jets} #geq 4)","M_{ee} [GeV]",40,60.,120.);
  ZMass_Zinc5jet=zJet::prettyHistogram(fs,"ZMass_Zinc5jet","Z Invariant Mass (N_{jets} #geq 5)","M_{ee} [GeV]",40,60.,120.);
  ZMass_Zinc6jet=zJet::prettyHistogram(fs,"ZMass_Zinc6jet","Z Invariant Mass (N_{jets} #geq 6)","M_{ee} [GeV]",40,60.,120.);
  
  ZMass_Zexc0jet=zJet::prettyHistogram(fs,"ZMass_Zexc0jet","Z Invariant Mass (N_{jets} = 0)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zexc1jet=zJet::prettyHistogram(fs,"ZMass_Zexc1jet","Z Invariant Mass (N_{jets} = 1)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zexc2jet=zJet::prettyHistogram(fs,"ZMass_Zexc2jet","Z Invariant Mass (N_{jets} = 2)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zexc3jet=zJet::prettyHistogram(fs,"ZMass_Zexc3jet","Z Invariant Mass (N_{jets} = 3)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zexc4jet=zJet::prettyHistogram(fs,"ZMass_Zexc4jet","Z Invariant Mass (N_{jets} = 4)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zexc5jet=zJet::prettyHistogram(fs,"ZMass_Zexc5jet","Z Invariant Mass (N_{jets} = 5)","M_{ee} [GeV]",60,60.,120.);
  ZMass_Zexc6jet=zJet::prettyHistogram(fs,"ZMass_Zexc6jet","Z Invariant Mass (N_{jets} = 6)","M_{ee} [GeV]",20,60.,120.);
  
  ZPt_Zinc0jet=zJet::prettyHistogram(fs,"ZPt_Zinc0jet","Z p_{T} (N_{jets} #geq 0)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zinc1jet=zJet::prettyHistogram(fs,"ZPt_Zinc1jet","Z p_{T} (N_{jets} #geq 1)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zinc2jet=zJet::prettyHistogram(fs,"ZPt_Zinc2jet","Z p_{T} (N_{jets} #geq 2)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zinc3jet=zJet::prettyHistogram(fs,"ZPt_Zinc3jet","Z p_{T} (N_{jets} #geq 3)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zinc4jet=zJet::prettyHistogram(fs,"ZPt_Zinc4jet","Z p_{T} (N_{jets} #geq 4)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zinc5jet=zJet::prettyHistogram(fs,"ZPt_Zinc5jet","Z p_{T} (N_{jets} #geq 5)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zinc6jet=zJet::prettyHistogram(fs,"ZPt_Zinc6jet","Z p_{T} (N_{jets} #geq 6)","p_{T}(Z) [GeV]",25,0.,500.);
  
  
  ZPt_Zexc0jet=zJet::prettyHistogram(fs,"ZPt_Zexc0jet","Z p_{T} (N_{jets} = 0)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zexc1jet=zJet::prettyHistogram(fs,"ZPt_Zexc1jet","Z p_{T} (N_{jets} = 1)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zexc2jet=zJet::prettyHistogram(fs,"ZPt_Zexc2jet","Z p_{T} (N_{jets} = 2)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zexc3jet=zJet::prettyHistogram(fs,"ZPt_Zexc3jet","Z p_{T} (N_{jets} = 3)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zexc4jet=zJet::prettyHistogram(fs,"ZPt_Zexc4jet","Z p_{T} (N_{jets} = 4)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zexc5jet=zJet::prettyHistogram(fs,"ZPt_Zexc5jet","Z p_{T} (N_{jets} = 5)","p_{T}(Z) [GeV]",25,0.,500.);
  ZPt_Zexc6jet=zJet::prettyHistogram(fs,"ZPt_Zexc6jet","Z p_{T} (N_{jets} = 6)","p_{T}(Z) [GeV]",25,0.,500.);
  
  
  ZRapidity_Zinc0jet=zJet::prettyHistogram(fs,"ZRapidity_Zinc0jet","ZRapidity (N_{jets} #geq 0)","Y(Z)",30,-3.,3.);
  ZRapidity_Zinc1jet=zJet::prettyHistogram(fs,"ZRapidity_Zinc1jet","ZRapidity (N_{jets} #geq 1)","Y(Z)",30,-3.,3.);
  ZRapidity_Zinc2jet=zJet::prettyHistogram(fs,"ZRapidity_Zinc2jet","ZRapidity (N_{jets} #geq 2)","Y(Z)",30,-3.,3.);
  ZRapidity_Zinc3jet=zJet::prettyHistogram(fs,"ZRapidity_Zinc3jet","ZRapidity (N_{jets} #geq 3)","Y(Z)",30,-3.,3.);
  ZRapidity_Zinc4jet=zJet::prettyHistogram(fs,"ZRapidity_Zinc4jet","ZRapidity (N_{jets} #geq 4)","Y(Z)",30,-3.,3.);
  ZRapidity_Zinc5jet=zJet::prettyHistogram(fs,"ZRapidity_Zinc5jet","ZRapidity (N_{jets} #geq 5)","Y(Z)",30,-3.,3.);
  ZRapidity_Zinc6jet=zJet::prettyHistogram(fs,"ZRapidity_Zinc6jet","ZRapidity (N_{jets} #geq 6)","Y(Z)",30,-3.,3.);
   
   
  ZRapidity_Zexc0jet=zJet::prettyHistogram(fs,"ZRapidity_Zexc0jet","ZRapidity (N_{jets} = 0)","Y(Z)",10,0.,3.);
  ZRapidity_Zexc1jet=zJet::prettyHistogram(fs,"ZRapidity_Zexc1jet","ZRapidity (N_{jets} = 1)","Y(Z)",10,0.,3.);
  ZRapidity_Zexc2jet=zJet::prettyHistogram(fs,"ZRapidity_Zexc2jet","ZRapidity (N_{jets} = 2)","Y(Z)",10,0.,3.);
  ZRapidity_Zexc3jet=zJet::prettyHistogram(fs,"ZRapidity_Zexc3jet","ZRapidity (N_{jets} = 3)","Y(Z)",10,0.,3.);
  ZRapidity_Zexc4jet=zJet::prettyHistogram(fs,"ZRapidity_Zexc4jet","ZRapidity (N_{jets} = 4)","Y(Z)",10,0.,3.);
  ZRapidity_Zexc5jet=zJet::prettyHistogram(fs,"ZRapidity_Zexc5jet","ZRapidity (N_{jets} = 5)","Y(Z)",10,0.,3.);
  ZRapidity_Zexc6jet=zJet::prettyHistogram(fs,"ZRapidity_Zexc6jet","ZRapidity (N_{jets} = 6)","Y(Z)",10,0.,3.);
  
  
  FirstJetY_Z1jet=zJet::prettyHistogram(fs,"FirstJetY_Z1jet","First Leading Jet (Z+Jets)","Y (Jet)",10,0.,3.);
  SecondJetY_Z2jet=zJet::prettyHistogram(fs,"SecondJetY_Z2jet","Second Leading Jet (Z+Jets)","Y(Jet) ",10,0.,3.);
  ThirdJetY_Z3jet=zJet::prettyHistogram(fs,"ThirdJetY_Z3jet","Third Leading (Z+Jets)","Y(Jets) ",10,0.,3.);
  FourthJetY_Z4jet=zJet::prettyHistogram(fs,"FourthJetY_Z4jet","Fourth Leading Jet (Z+Jets)","Y(Jet) ",10,0.,3.);
  FifthJetY_Z5jet=zJet::prettyHistogram(fs,"FifthJetY_Z5jet","Fifth Leading Jet (Z+Jets)","Y(Jet) ",10,0.,3.);
  SixthJetY_Z6jet=zJet::prettyHistogram(fs,"SixthJetY_Z6jet","#geq 6th Jets (Z+Jets)","Y(Jet) ",10,0.,3.);
  
  hAccept = zJet::prettyHistogram(fs,"hAccept","Accept","A",50,0.,1.);
  hCounter = zJet::prettyHistogram(fs,"hCounter","Counter","A",10,0.,10.);
}


ZJet::~ZJet()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  hCounter->Fill(1);
  using namespace edm;
  bool usePatElec = 0;
  bool useGsfElec = 1;

  edm::Handle<pat::ElectronCollection> patElecHandle;
  edm::Handle<reco::GsfElectronCollection> gsfElecHandle;
  //std::cout<<"0"<<std::endl;
  if(usePatElec){
    if(not iEvent.getByLabel("selectedPatElectronsPFlow",patElecHandle)){
      std::cout<<"FATAL EXCEPTION: "<<"patElectrons Not Found: "<<std::endl;
      return;
    }
  }
  //std::cout<<"1"<<std::endl;
  if(useGsfElec){
    if(not iEvent.getByLabel("gsfElectrons", gsfElecHandle)){
      std::cout<<"FATAL EXCEPTION: "<<"gsfElectrons Not Found: "<<std::endl;
      return;
    }
  }
  
  pat::Electron* e1Ptr=0; 
  pat::Electron* e2Ptr=0; 
  
  //std::cout<<"4"<<std::endl;
  
  if(usePatElec){
    pat::ElectronCollection patEleColl(*(patElecHandle.product()));
    if(patEleColl.size()<2)return;

    e1Ptr = new pat::Electron(patEleColl[0]);
    e2Ptr = new pat::Electron(patEleColl[1]);
  }
  
  if(useGsfElec){
    reco::GsfElectronCollection gsfEleColl(*(gsfElecHandle.product())); 
    //std::cout<<"Found: "<<gsfEleColl.size()<<std::endl;
    if(gsfEleColl.size()<2)return;
    e1Ptr = new pat::Electron(gsfEleColl[0]);
    e2Ptr = new pat::Electron(gsfEleColl[1]);
  }

  //std::cout<<"5"<<std::endl;

  pat::Electron e1(*e1Ptr);
  pat::Electron e2(*e2Ptr);

  double pt1     = e1.pt();
  double eta1     = e1.superCluster()->eta();
  double phi1     = e1.phi();
  double tiso1    = e1.dr03TkSumPt();
  double eiso1    = e1.dr03EcalRecHitSumEt();
  double hiso1    = e1.dr03HcalTowerSumEt();
  double sieie1   = e1.sigmaIetaIeta();
  double delphi1  = e1.deltaPhiSuperClusterTrackAtVtx();
  double detain1  = e1.deltaEtaSuperClusterTrackAtVtx();
  double dcot1    = e1.convDcot();
  double dist1    = e1.convDist();
  double hoe1     = e1.hadronicOverEm();
  double nmshits1 = e1.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits();
  
  double pt2     = e2.pt();
  double eta2     = e2.superCluster()->eta();
  double phi2     = e2.phi();
  double tiso2    = e2.dr03TkSumPt();
  double eiso2    = e2.dr03EcalRecHitSumEt();
  double hiso2    = e2.dr03HcalTowerSumEt();
  double sieie2   = e2.sigmaIetaIeta();
  double delphi2  = e2.deltaPhiSuperClusterTrackAtVtx();
  double detain2  = e2.deltaEtaSuperClusterTrackAtVtx();
  double dcot2    = e2.convDcot();
  double dist2    = e2.convDist();
  double hoe2     = e2.hadronicOverEm();
  double nmshits2 = e2.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits();
    
  TLorentzVector l1 = Part2LorVec(e1);
  TLorentzVector l2 = Part2LorVec(e2);
  
  TLorentzVector p= l1+l2;
  z0Hist_ ->Fill(p.M());//Fill before any cut.
  if((p.M()<60)||(p.M()>120))return;
  
  
  //Evaluate Kinematic Cuts
  bool bothPassPtEta = pt1 > 20 
    && fabs(eta1)<2.4 && !(fabs(eta1)>1.446 && fabs(eta1)<1.566)  
    && pt2 > 20 
    && fabs(eta2)<2.4 && !(fabs(eta2)>1.446 && fabs(eta2)<1.566);
  
  // !(-0.02<convDist<0.02 && -0.02<convDcot<0.02) 
  
  //Evaluate Electron Identification cuts
  bool e1PassId = 0;  
  if(e1.isEB()){    
    e1PassId = sieie1 <0.01&& fabs(delphi1)<0.06&& fabs(detain1)<0.004&& 
      !(fabs(dcot1)<0.02&& fabs(dist1)<0.02)&& hoe1 <0.04&& nmshits1<=0;
  } else if(e1.isEE()){
    e1PassId = sieie1 <0.03&& fabs(delphi1)<0.03&& fabs(detain1)<0.007&& 
      !(fabs(dcot1)<0.02&& fabs(dist1)<0.02)&& hoe1 <0.15&& nmshits1<=0;
  }
  
  bool e2PassId = 0;
  if(e2.isEB()){    
    e2PassId = sieie2 <0.01&& fabs(delphi2)<0.06&& fabs(detain2)<0.004&& 
      !(fabs(dcot2)<0.02&& fabs(dist2)<0.02)&& hoe2 <0.04&& nmshits2<=0;
  } 
  else if(e2.isEE()){
    e2PassId = sieie2<0.03&& fabs(delphi2)<0.03&& fabs(detain2)<0.007&& 
      !(fabs(dcot2)<0.02&& fabs(dist2)<0.02)&& hoe2 <0.15&& nmshits2<=0;
  }
  
  bool bothPassId = e1PassId&&e2PassId;
  
  //Evaluate Isolation Cuts:
  double lepIsoRho = -9999999999999;
  //edm::Handle<double> rhoLepIso= -999999999;
  //const edm::InputTag eventrhoLepIso("kt6PFJets", "rho");
  //iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
  //if( *rhoLepIso == *rhoLepIso) lepIsoRho = *rhoLepIso;
  //else  lepIsoRho =  -999999.9;
  
  Double_t Atracker[2] = {0., 0.}; //   barrel/endcap
  Double_t Aecal[2]    = {0.101, 0.046}; //   barrel/endcap
  Double_t Ahcal[2]    = {0.021 , 0.040 }; //   barrel/endcap
  enum detID{barrel=0, endcap=1};
  
  double combRelRhoSubtractedIso1 = 0;
  double combRelRhoSubtractedIso2 = 0;
  
  
  if(e1.isEB())combRelRhoSubtractedIso1 = (tiso1 + zJet::max(0. ,eiso1 - Aecal[barrel]*(lepIsoRho)) + zJet::max(0.,hiso1 - Ahcal[barrel]*(lepIsoRho)) )/pt1;
  
  else if(e1.isEE())combRelRhoSubtractedIso1 = (tiso1 + zJet::max(0. ,eiso1 - Aecal[endcap]*(lepIsoRho)) + zJet::max(0.,hiso1 - Ahcal[endcap]*(lepIsoRho)) )/pt1;
  
  if(e2.isEB())combRelRhoSubtractedIso2 = (tiso2 + zJet::max(0. ,eiso2 - Aecal[barrel]*(lepIsoRho)) + zJet::max(0.,hiso2 - Ahcal[barrel]*(lepIsoRho)) )/pt2;
  
  else if(e2.isEE())combRelRhoSubtractedIso2 = (tiso2 + zJet::max(0. ,eiso2 - Aecal[endcap]*(lepIsoRho)) + zJet::max(0.,hiso2 - Ahcal[endcap]*(lepIsoRho)) )/pt2;

  
  bool bothPassIso = combRelRhoSubtractedIso1<0.15 && combRelRhoSubtractedIso2<0.15;
  
  //Conversion Rejection
  bool conv1=0;
  bool conv2=0;
  
  
  if(e1.isEB() && !(fabs(dcot1)<0.02&& fabs(dist1)<=0.02)&& nmshits1<=0) conv1=1;
  if(e1.isEE() && !(fabs(dcot1)<0.02&& fabs(dist1)<=0.02)&& nmshits1<=0) conv1=1;
  if(e2.isEB() && !(fabs(dcot2)<0.02&& fabs(dist2)<=0.02)&& nmshits2<=0) conv2=1;
  if(e2.isEE() && !(fabs(dcot2)<0.02&& fabs(dist2)<=0.02)&& nmshits2<=0) conv2=1;
  bool bothPassConv = conv1 && conv2;
  
  
  
  
  
  //Apply Cuts and Count the Events.
  z1Hist_ ->Fill(p.M());
  if(bothPassPtEta)
    z2Hist_ ->Fill(p.M());
  if(bothPassPtEta &&  bothPassId)
    z3Hist_ ->Fill(p.M());  
  if(bothPassPtEta && bothPassId && bothPassIso)
    z4Hist_ ->Fill(p.M()); 
  if(bothPassPtEta && bothPassId && bothPassIso && bothPassConv)
    z5Hist_ ->Fill(p.M());
 

 if(!(bothPassPtEta && bothPassId && bothPassIso && bothPassConv))return;
 
 std::vector<TLorentzVector> jets;
 bool usePatJets=0;
 bool useRecJets=1;

 if(usePatJets){
   edm::Handle<std::vector<pat::Jet> > patJetHandle;
   if(not iEvent.getByLabel("selectedPatJetsPFlow",patJetHandle)){
     std::cout<<"PAT Jets Not Found" <<std::endl;
     return;
   }

   const std::vector<pat::Jet>*  patJets = patJetHandle.product();
   std::vector<pat::Jet>::const_iterator jet =patJets->begin();
   for(;jet!=patJets->end();jet++){
     double jetPt = jet->pt();
     if(jetPt<30)continue;
     double jetEta = jet->eta();
     double jetPhi = jet->phi();
     double delR_e1 = zJet::radius(jetEta,jetPhi,eta1,phi1);
     double delR_e2 = zJet::radius(jetEta,jetPhi,eta2,phi2);
     if((delR_e1<0.3)||(delR_e2<0.3))continue;
     TLorentzVector* theJet = new TLorentzVector();
     theJet->SetPxPyPzE(jet->px(),jet->py(),jet->pz(),jet->energy());
     jets.push_back(*theJet);
   }
 }
  
 if(useRecJets){
   edm::Handle<std::vector<reco::PFJet> > recJetHandle;
   if(not iEvent.getByLabel("ak5PFJets",recJetHandle)){
     std::cout<<"RECO Jets Not Found" <<std::endl;
     return; 
   }
   const std::vector<reco::PFJet>* recoJets = recJetHandle.product();
   std::vector<reco::PFJet>::const_iterator jet =recoJets->begin();
   for(;jet!=recoJets->end();jet++){
     double jetPt = jet->pt();
     if(jetPt<30)continue;
     double jetEta = jet->eta();
     double jetPhi = jet->phi();
     double delR_e1 = zJet::radius(jetEta,jetPhi,eta1,phi1);
     double delR_e2 = zJet::radius(jetEta,jetPhi,eta2,phi2);
     if((delR_e1<0.3)||(delR_e2<0.3))continue;
     TLorentzVector* theJet = new TLorentzVector();
     theJet->SetPxPyPzE(jet->px(),jet->py(),jet->pz(),jet->energy());
     jets.push_back(*theJet);
    }
 }
  
  
  double weight = 1;
  

 hNjetGoodZX->Fill(jets.size(),weight);
  
   // if(jets.size()>=0){
      hNjetGoodIncZX->Fill(0.,weight);
      ZMass_Zinc0jet->Fill(p.M(),weight);
      ZPt_Zinc0jet->Fill(p.Pt(),weight);
      ZRapidity_Zinc0jet->Fill(fabs(p.Rapidity()),weight);
      if(jets.size()==0){
    	ZMass_Zexc0jet->Fill(p.M(),weight);
    	ZPt_Zexc0jet->Fill(p.Pt(),weight);
    	ZRapidity_Zexc0jet->Fill(fabs(p.Rapidity()),weight);
 	hCounter->Fill(0.,weight);
      }
    //}





   if(jets.size()>=1){
     
     hNjetGoodIncZX->Fill(1., weight);
     ZMass_Zinc1jet->Fill(p.M(),weight);
     ZPt_Zinc1jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc1jet->Fill(fabs(p.Rapidity()),weight);
     if(jets.size()==1){
       ZMass_Zexc1jet->Fill(p.M(),weight);
       ZPt_Zexc1jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc1jet->Fill(fabs(p.Rapidity()),weight);
       FirstJetPt_Z1jet->Fill(jets[0].Pt(),weight);
       FirstJetY_Z1jet->Fill(fabs(jets[0].Eta()),weight);
       hCounter->Fill(1.,weight);
     }
   }
   
     //std::cout<<"3: "<<std::endl;
   if(jets.size()>=2){
     
     hNjetGoodIncZX->Fill(2., weight);
     ZMass_Zinc2jet->Fill(p.M(),weight);
     ZPt_Zinc2jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc2jet->Fill(fabs(p.Rapidity()),weight);
     if(jets.size()==2){
       ZMass_Zexc2jet->Fill(p.M(),weight);
       ZPt_Zexc2jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc2jet->Fill(fabs(p.Rapidity()),weight);
       SecondJetPt_Z2jet->Fill(jets[1].Pt(),weight);
       SecondJetY_Z2jet->Fill(fabs(jets[1].Eta()),weight);
       hCounter->Fill(2.,weight);
     }
   }
     //std::cout<<"4: "<<std::endl;

   if(jets.size()>=3){
     
     hNjetGoodIncZX->Fill(3., weight);
     ZMass_Zinc3jet->Fill(p.M(),weight);
     ZPt_Zinc3jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc3jet->Fill(fabs(p.Rapidity()),weight);
     if(jets.size()==3){
       ZMass_Zexc3jet->Fill(p.M(),weight);
       ZPt_Zexc3jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc3jet->Fill(fabs(p.Rapidity()),weight);
       ThirdJetPt_Z3jet->Fill(jets[2].Pt(),weight);
       ThirdJetY_Z3jet->Fill(fabs(jets[2].Eta()),weight);
       hCounter->Fill(3.,weight);
     }
   }
   
  //std::cout<<"5: "<<std::endl;

   if(jets.size()>=4){
     
     hNjetGoodIncZX->Fill(4., weight);
     ZMass_Zinc4jet->Fill(p.M(),weight);
     ZPt_Zinc4jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc4jet->Fill(fabs(p.Rapidity()),weight);
     if(jets.size()==4){
       ZMass_Zexc4jet->Fill(p.M(),weight);
       ZPt_Zexc4jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc4jet->Fill(fabs(p.Rapidity()),weight);
       FourthJetPt_Z4jet->Fill(jets[3].Pt(),weight);
       FourthJetY_Z4jet->Fill(fabs(jets[3].Eta()),weight);
       hCounter->Fill(4.,weight);
     }
   }
  //std::cout<<"6: "<<std::endl;

   if(jets.size()>=5){
     
     hNjetGoodIncZX->Fill(5., weight);
     ZMass_Zinc5jet->Fill(p.M(),weight);
     ZPt_Zinc5jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc5jet->Fill(p.Rapidity(),weight);
     if(jets.size()==5){
       ZMass_Zexc5jet->Fill(p.M(),weight);
       ZPt_Zexc5jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc5jet->Fill(p.Rapidity(),weight);
       FifthJetPt_Z5jet->Fill(jets[4].Pt(),weight);
       FifthJetY_Z5jet->Fill(jets[4].Eta(),weight);
       hCounter->Fill(5.,weight);
     }
   }
   
   if(jets.size()>=6){
     
     hNjetGoodIncZX->Fill(6., weight);
     ZMass_Zinc6jet->Fill(p.M(),weight);
     ZPt_Zinc6jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc6jet->Fill(p.Rapidity(),weight);
     SixthJetPt_Z6jet->Fill(jets[5].Pt(),weight);
     SixthJetY_Z6jet->Fill(jets[5].Eta(),weight);
     if(jets.size()==6){
       ZMass_Zexc6jet->Fill(p.M(),weight);
       ZPt_Zexc6jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc6jet->Fill(p.Rapidity(),weight);
       hCounter->Fill(6.,weight);
     }
     
   }
 }


// ------------ method called once each job just before starting event loop  ------------
void 
ZJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZJet::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ZJet::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZJet::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZJet::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZJet::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZJet);


/*
  
    bool id1=0;
    bool id2=0;
    
    if(e1.isEB() && sieie1 <0.01&& fabs(delphi1)<0.06&& fabs(detain1)<0.004 && hoe1 <0.04) id1=1;
    if(e1.isEE() && sieie1 <0.03&& fabs(delphi1)<0.03&& fabs(detain1)<0.007&&hoe1 <0.15) id1=1;
    if(e2.isEB() && sieie2 <0.01&& fabs(delphi2)<0.06&& fabs(detain2)<0.004&&hoe2 <0.04) id2=1;
    if(e2.isEE() && sieie1 <0.03&& fabs(delphi2)<0.03&& fabs(detain2)<0.007&&hoe2 <0.15) id2=1;
  
//Our Current Isolation Cuts: Record

  bool e1PassIso = 0; 
  if(e1.isEB()){  e1PassIso = tiso1/pt1<0.09 && eiso1/pt1<0.07 && hiso1/pt1<0.10;}
  else if(e1.isEE()){ e1PassIso=tiso1/pt1<0.04 && eiso1/pt1<0.05 && hiso1/pt1<0.025;}
  
  bool e2PassIso = 0;
  if(e2.isEB()){e2PassIso = tiso2/pt2<0.09 && eiso2/pt2<0.07 && hiso2/pt2<0.10;}
  else if(e2.isEE()){e2PassIso=tiso2/pt2<0.04 && eiso2/pt2<0.05 && hiso2/pt2<0.025;}
   



 */
