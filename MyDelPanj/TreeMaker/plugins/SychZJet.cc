// system include files
#include <memory>

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

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
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"





#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1I.h"
#include "TH1F.h"

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"


//
// class declaration
//
using namespace edm;
using namespace std;
using namespace reco;

//
// class declaration
//

namespace SychzJet{
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

}

class SychZJet : public edm::EDAnalyzer {
   public:
      explicit SychZJet(const edm::ParameterSet&);
      ~SychZJet();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  edm::LumiReWeighting LumiWeights_;
  edm::InputTag vertexSrc_;
std::vector< float > TrueDist2011_;
std::vector<float> WLumi_;


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
  TH1D* hCounter ; 
  TH1D* hWeights ;

  TH1D* hNjetGoodZX; 
  TH1D* hNjetGoodIncZX ; 
  TH1D* ZMass_Zinc0jet;
  TH1D* ZMass_Zinc1jet;
  TH1D* ZMass_Zinc2jet;
  TH1D* ZMass_Zinc3jet;
  TH1D* ZMass_Zinc4jet;
  TH1D* ZMass_Zinc5jet;
  TH1D* ZMass_Zinc6jet;
  
  TH1D* FirstJetPt_Z1jet;
  TH1D* SecondJetPt_Z2jet;
  TH1D* ThirdJetPt_Z3jet;
  TH1D* FourthJetPt_Z4jet;
  TH1D* FifthJetPt_Z5jet;
  TH1D* SixthJetPt_Z6jet;

  TH1D* FirstJetY_Z1jet;
  TH1D* SecondJetY_Z2jet;
  TH1D* ThirdJetY_Z3jet;
  TH1D* FourthJetY_Z4jet;
  TH1D* FifthJetY_Z5jet;
  TH1D* SixthJetY_Z6jet;
  
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
SychZJet::SychZJet(const edm::ParameterSet& iConfig)

{


Double_t Fall2011[50] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08
  };


double Data[50]={2.44414e+06,1.25579e+06,3.42957e+06,5.18618e+07,2.54054e+08,4.36586e+08,4.86031e+08,4.63551e+08,4.18993e+08,3.84891e+08,3.65724e+08,3.41505e+08,3.22711e+08,3.0886e+08,2.87693e+08,2.5129e+08,1.99438e+08,1.40551e+08,8.66577e+07,4.6234e+07,2.12058e+07,8.37396e+06,2.88178e+06,882886,246537,63900.8,15571.6,3628.24,840.61,211.248,66.7507,32.1445,26.92,33.3738,47.1181,66.9628,92.522,123.311,158.278,195.606,232.736,266.603,294.026,312.195,319.142,314.095,297.617,271.501,238.455,1.93299e+07};

 for( int i=0; i<50; ++i) {
      TrueDist2011_.push_back(Data[i]);
      WLumi_.push_back(Fall2011[i]);
   }
 LumiWeights_ = edm::LumiReWeighting(WLumi_,TrueDist2011_);


   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  z0Hist_ = fs->make<TH1D>("hZ0Mass","hZ0Mass",50,0.,200.);
  z1Hist_ = fs->make<TH1D>("hZ1Mass","hZ1Mass",50,0.,200.);
  z2Hist_ = fs->make<TH1D>("hZ2Mass","hZ2Mass",50,0.,200.);
  z3Hist_ = fs->make<TH1D>("hZ3Mass","hZ3Mass",50,0.,200.);
  z4Hist_ = fs->make<TH1D>("hZ4Mass","hZ4Mass",50,0.,200.);
  z5Hist_ = fs->make<TH1D>("hZ5Mass","hZ5Mass",42,70.,112.);
  z6Hist_ = fs->make<TH1D>("hZ6Mass","hZ6Mass",42,70.,112.);
  z7Hist_ = fs->make<TH1D>("hZ7Mass","hZ7Mass",42,70.,112.);
  z8Hist_ = fs->make<TH1D>("hZ8Mass","hZ8Mass",42,70.,112.);
  z9Hist_ = fs->make<TH1D>("hZ9Mass","hZ9Mass",42,70.,112.);
  
  hNjetGoodZX = SychzJet::prettyHistogram(fs,"h_Njet_GoodZX","N_{Jets}(Exc) [Z+Jets]","NJets_{Exc}",10,0,10);

  hNjetGoodIncZX = SychzJet::prettyHistogram(fs,"h_Njet_GoodIncZX","N_{Jets}(Inc) [Z+Jets]","NJets_{Inc}",10,0,10);
  ZMass_Zinc0jet=SychzJet::prettyHistogram(fs,"ZMass_Zinc0jet","Z Invariant Mass (N_{jets} #geq 0)","M_{ee} [GeV]",42,70.,112.);
  ZMass_Zinc1jet=SychzJet::prettyHistogram(fs,"ZMass_Zinc1jet","Z Invariant Mass (N_{jets} #geq 1)","M_{ee} [GeV]",42,70.,112.);
  ZMass_Zinc2jet=SychzJet::prettyHistogram(fs,"ZMass_Zinc2jet","Z Invariant Mass (N_{jets} #geq 2)","M_{ee} [GeV]",42,70.,112.);
  ZMass_Zinc3jet=SychzJet::prettyHistogram(fs,"ZMass_Zinc3jet","Z Invariant Mass (N_{jets} #geq 3)","M_{ee} [GeV]",42,70.,112.);
  ZMass_Zinc4jet=SychzJet::prettyHistogram(fs,"ZMass_Zinc4jet","Z Invariant Mass (N_{jets} #geq 4)","M_{ee} [GeV]",42,70.,112.);
  ZMass_Zinc5jet=SychzJet::prettyHistogram(fs,"ZMass_Zinc5jet","Z Invariant Mass (N_{jets} #geq 5)","M_{ee} [GeV]",42,70.,112.);
  ZMass_Zinc6jet=SychzJet::prettyHistogram(fs,"ZMass_Zinc6jet","Z Invariant Mass (N_{jets} #geq 6)","M_{ee} [GeV]",42,70.,112.);

  hCounter = SychzJet::prettyHistogram(fs,"hCounter","Counter","EventCounter",16,0.5,16.5);
  hWeights = SychzJet::prettyHistogram(fs,"hWeights","Weight","Weight",50,0.,10.);
  //binning jetpT





  float jpt[10] = {30,50,70,90,120,150,180,240,300,400};
  float jpt2[9] = {30,50,70,90,120,150,200,270,400};
  float jpt3[6] = {30,50,70,90,120,400};
 // float jpt4[6] = {30,50,70,90,120,400};
 FirstJetPt_Z1jet=SychzJet::prettyHistogram(fs,"FirstJetPt_Z1jet","First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
 SecondJetPt_Z2jet=SychzJet::prettyHistogram(fs,"SecondJetPt_Z2jet","Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
 ThirdJetPt_Z3jet=SychzJet::prettyHistogram(fs,"ThirdJetPt_Z3jet","Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
 FourthJetPt_Z4jet=SychzJet::prettyHistogram(fs,"FourthJetPt_Z4jet","Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
 FifthJetPt_Z5jet=SychzJet::prettyHistogram(fs,"FifthJetPt_Z5jet","Fifth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
 SixthJetPt_Z6jet=SychzJet::prettyHistogram(fs,"SixthJetPt_Z6jet","#geq 6th Jets (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);

 FirstJetY_Z1jet=SychzJet::prettyHistogram(fs,"FirstJetY_Z1jet","First Leading Jet (Z+Jets)","#eta (Jet)",9,0.,3.);
 SecondJetY_Z2jet=SychzJet::prettyHistogram(fs,"SecondJetY_Z2jet","Second Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
 ThirdJetY_Z3jet=SychzJet::prettyHistogram(fs,"ThirdJetY_Z3jet","Third Leading (Z+Jets)","#eta (Jet) ",9,0.,3.);
 FourthJetY_Z4jet=SychzJet::prettyHistogram(fs,"FourthJetY_Z4jet","Fourth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
 FifthJetY_Z5jet=SychzJet::prettyHistogram(fs,"FifthJetY_Z5jet","Fifth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
 SixthJetY_Z6jet=SychzJet::prettyHistogram(fs,"SixthJetY_Z6jet","#geq 6th Jets (Z+Jets)","#eta (Jet) ",9,0.,3.);




}


SychZJet::~SychZJet()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SychZJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

  std::vector<PileupSummaryInfo>::const_iterator PVI;


  float npT=-1.;
  float npIT=-1.;


  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

    int BX = PVI->getBunchCrossing();

    if(BX == 0) {
      npT = PVI->getTrueNumInteractions();
      npIT = PVI->getPU_NumInteractions();
    }
  }

  double MyWeight = LumiWeights_.weight( npT );
hWeights->Fill(MyWeight);
//cout<<MyWeight<<endl;
//  double weight=1.;
  double weight=MyWeight;
  hCounter->Fill(1,weight);
  using namespace edm;
  bool usePatElec = 1;
  bool useGsfElec = 0;

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
  
//  if(usePatElec){
    pat::ElectronCollection patEleColl(*(patElecHandle.product()));
  if(usePatElec){
    std::sort(patEleColl.begin(),patEleColl.end(),PtGreater());
    //std::cout<<"Found: "<<patEleColl.size()<<std::endl;
    if(patEleColl.size()<2)return;

    e1Ptr = new pat::Electron(patEleColl[0]);
    e2Ptr = new pat::Electron(patEleColl[1]);
  }
    if(patEleColl.size()>=2)  hCounter->Fill(2,weight);
  
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

double q1 = e1.threeCharge();
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
  
double q2 = e2.threeCharge();
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
  //cout<<q1*q2<<endl;  

  TLorentzVector l1 = Part2LorVec(e1);
  TLorentzVector l2 = Part2LorVec(e2);
  
  TLorentzVector p= l1+l2;
hCounter->Fill(3,weight);

  z3Hist_ ->Fill(p.M(),weight);//Fill before any cut.

if(!(q1*q2<0))return;
hCounter->Fill(4,weight);
  z4Hist_ ->Fill(p.M(),weight);//Fill before any cut.
//  if((p.M()<60)||(p.M()>120))return;
  if(!(fabs(p.M()-91)<20) ) return;  
hCounter->Fill(5,weight);
  z5Hist_ ->Fill(p.M(),weight);
  //Evaluate Kinematic Cuts
  bool bothPassPtEta = pt1 > 20 
    && fabs(eta1)<2.4 && !(fabs(eta1)>1.446 && fabs(eta1)<1.566)  
    && pt2 > 20 
    && fabs(eta2)<2.4 && !(fabs(eta2)>1.446 && fabs(eta2)<1.566);
  
  //Evaluate Electron Identification cuts
  bool e1PassId = 0;  
  if(e1.isEB()){    
    e1PassId = sieie1 <0.01&& fabs(delphi1)<0.06&& fabs(detain1)<0.004&& hoe1 <0.04 ; 
  } else if(e1.isEE()){
    e1PassId = sieie1 <0.03&& fabs(delphi1)<0.03&& fabs(detain1)<0.007&& hoe1 <0.15 ; 
  }
  
  bool e2PassId = 0;
  if(e2.isEB()){    
    e2PassId = sieie2 <0.01&& fabs(delphi2)<0.06&& fabs(detain2)<0.004&& hoe2 <0.04 ;
  } 
  else if(e2.isEE()){
    e2PassId = sieie2<0.03&& fabs(delphi2)<0.03&& fabs(detain2)<0.007&& hoe2 <0.15 ;
  }
  
  bool bothPassId = e1PassId&&e2PassId;
  
  //Evaluate Isolation Cuts:
  double lepIsoRho;// = -9999999999999;
  edm::Handle<double> rhoLepIso;//= -999999999;
  const edm::InputTag eventrhoLepIso("kt6PFJetsForIsolation", "rho");
  iEvent.getByLabel(eventrhoLepIso, rhoLepIso);
  if( *rhoLepIso == *rhoLepIso) lepIsoRho = *rhoLepIso;
  else  lepIsoRho =  999999.9;
  
//  Double_t Atracker[2] = {0., 0.}; //   barrel/endcap
  Double_t Aecal[2]    = {0.101, 0.046}; //   barrel/endcap
  Double_t Ahcal[2]    = {0.021 , 0.040 }; //   barrel/endcap
  enum detID{barrel=0, endcap=1};
  
  double combRelRhoSubtractedIso1 = 0;
  double combRelRhoSubtractedIso2 = 0;
  
  
  if(e1.isEB())combRelRhoSubtractedIso1 = (tiso1 + max(0. ,eiso1 - Aecal[barrel]*(lepIsoRho)) + max(0.,hiso1 - Ahcal[barrel]*(lepIsoRho)) )/pt1;
  
  else if(e1.isEE())combRelRhoSubtractedIso1 = (tiso1 + max(0. ,eiso1 - Aecal[endcap]*(lepIsoRho)) + max(0.,hiso1 - Ahcal[endcap]*(lepIsoRho)) )/pt1;
  
  if(e2.isEB())combRelRhoSubtractedIso2 = (tiso2 + max(0. ,eiso2 - Aecal[barrel]*(lepIsoRho)) + max(0.,hiso2 - Ahcal[barrel]*(lepIsoRho)) )/pt2;
  
  else if(e2.isEE())combRelRhoSubtractedIso2 = (tiso2 + max(0. ,eiso2 - Aecal[endcap]*(lepIsoRho)) + max(0.,hiso2 - Ahcal[endcap]*(lepIsoRho)) )/pt2;

  
  bool bothPassCombDetIso = combRelRhoSubtractedIso1<0.15 && combRelRhoSubtractedIso2<0.15;


//Matteo Det-Iso
bool IsoPass1 = 0;
bool IsoPass2 = 0 ;
  if(e1.isEB()){ 
	IsoPass1 = (tiso1/pt1)< 0.09
	&& (max(0. ,eiso1 - Aecal[barrel]*(lepIsoRho))/pt1)<  0.07
	&&  (max(0.,hiso1 - Ahcal[barrel]*(lepIsoRho))/pt1)< 0.1; 
	}
  else if(e1.isEE()){ 
        IsoPass1 = (tiso1/pt1)< 0.04
        && (max(0. ,eiso1 - Aecal[barrel]*(lepIsoRho))/pt1)<  0.05
        &&  (max(0.,hiso1 - Ahcal[barrel]*(lepIsoRho))/pt1)< 0.025; 
        }
  if(e2.isEB()){ 
        IsoPass2 = (tiso2/pt2)< 0.09
        && (max(0. ,eiso2 - Aecal[barrel]*(lepIsoRho))/pt2)<  0.07
        &&  (max(0.,hiso2 - Ahcal[barrel]*(lepIsoRho))/pt2)< 0.1;
        }
  else if(e2.isEE()){ 
        IsoPass2 = (tiso2/pt2)< 0.04
        && (max(0. ,eiso2 - Aecal[barrel]*(lepIsoRho))/pt2)<  0.05
        &&  (max(0.,hiso2 - Ahcal[barrel]*(lepIsoRho))/pt2)< 0.025;
        }

  bool bothPassDetIso = IsoPass1 && IsoPass2 ;
	

//////////////////////////
  
  //Conversion Rejection
  bool conv1=0;
  bool conv2=0;
  
double CombRelpfIso1=  (  e1.chargedHadronIso() +  e1.neutralHadronIso() + e1.photonIso()) /pt1;
double CombRelpfIso2=  (  e2.chargedHadronIso() +  e2.neutralHadronIso() + e2.photonIso()) /pt2;
//cout<<CombRelpfIso1<<endl; 
  bool bothPassPfIso = CombRelpfIso1  < 0.2 && CombRelpfIso2 < 0.2;


  if(e1.isEB() && !(fabs(dcot1)<0.02&& fabs(dist1)<0.02)&& nmshits1<=0) conv1=1;
else  if(e1.isEE() && !(fabs(dcot1)<0.02&& fabs(dist1)<0.02)&& nmshits1<=0) conv1=1;
  if(e2.isEB() && !(fabs(dcot2)<0.02&& fabs(dist2)<0.02)&& nmshits2<=0) conv2=1;
else  if(e2.isEE() && !(fabs(dcot2)<0.02&& fabs(dist2)<0.02)&& nmshits2<=0) conv2=1;
  bool bothPassConv = conv1 && conv2;
  
  
  
  //Apply Cuts and Count the Events.
 
  if(bothPassPtEta)
{hCounter->Fill(6,weight);
    z6Hist_ ->Fill(p.M(),weight);}
  if(bothPassPtEta &&  bothPassId)
{hCounter->Fill(7,weight);
    z7Hist_ ->Fill(p.M(),weight);  }
  if(bothPassPtEta && bothPassId && bothPassPfIso)
{hCounter->Fill(8,weight);
    z8Hist_ ->Fill(p.M(),weight); }
  if(bothPassPtEta && bothPassId && bothPassPfIso)
{
    z1Hist_ ->Fill(p.M(),weight); }

  if(bothPassPtEta && bothPassId && bothPassPfIso && bothPassConv)
{hCounter->Fill(9,weight);
    z9Hist_ ->Fill(p.M(),weight);}
 

 if(!(bothPassPtEta && bothPassId && bothPassPfIso && bothPassConv))return;
hCounter->Fill(10,weight);
 
 std::vector<TLorentzVector> jets;
 bool usePatJets=1;
 bool useRecJets=0;

 if(usePatJets){
   edm::Handle<std::vector<pat::Jet> > patJetHandle;
   if(not iEvent.getByLabel("selectedPatJetsPFlow",patJetHandle)){
     std::cout<<"PAT Jets Not Found" <<std::endl;
     return;
   }

   const std::vector<pat::Jet>*  patJets = patJetHandle.product();
   std::vector<pat::Jet>::const_iterator jet =patJets->begin();
   for(;jet!=patJets->end();jet++){
     bool jetID = jet->chargedHadronEnergyFraction() > 0 && 
	jet->chargedEmEnergyFraction()<0.99 && 
	jet->nConstituents()>1 && 
	jet->neutralHadronEnergyFraction()<0.99 && 
	jet->neutralEmEnergyFraction()<0.99 &&
	jet->chargedMultiplicity()>0.0;
     if(!(jetID)) continue;
     double jetPt = jet->pt();
     if(jetPt<=30)continue;
     double jetEta = jet->eta();
     if(abs(jetEta)>2.4)continue;
     double jetPhi = jet->phi();
     double delR_e1 = SychzJet::radius(jetEta,jetPhi,eta1,phi1);
     double delR_e2 = SychzJet::radius(jetEta,jetPhi,eta2,phi2);
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
     if(abs(jetEta)>2.5)continue;
     double jetPhi = jet->phi();
     double delR_e1 = SychzJet::radius(jetEta,jetPhi,eta1,phi1);
     double delR_e2 = SychzJet::radius(jetEta,jetPhi,eta2,phi2);
     if((delR_e1<0.3)||(delR_e2<0.3))continue;
     TLorentzVector* theJet = new TLorentzVector();
     theJet->SetPxPyPzE(jet->px(),jet->py(),jet->pz(),jet->energy());
     jets.push_back(*theJet);
    }
 }
  
  



 hNjetGoodZX->Fill(jets.size(),weight);
  
   // if(jets.size()>=0){
      hNjetGoodIncZX->Fill(0.,weight);
      ZMass_Zinc0jet->Fill(p.M(),weight);
    //}


   if(jets.size()>=1){
     hCounter->Fill(11,weight);

     hNjetGoodIncZX->Fill(1., weight);
     ZMass_Zinc1jet->Fill(p.M(),weight);
     FirstJetPt_Z1jet->Fill(jets.at(0).Pt(),weight);
     FirstJetY_Z1jet->Fill(jets.at(0).Eta(),weight);
   }
   
     //std::cout<<"3: "<<std::endl;
   if(jets.size()>=2){
     hCounter->Fill(12,weight);
     hNjetGoodIncZX->Fill(2., weight);
     ZMass_Zinc2jet->Fill(p.M(),weight);
     SecondJetPt_Z2jet->Fill(jets.at(1).Pt(),weight);
     SecondJetY_Z2jet->Fill(jets.at(1).Eta(),weight);
   }
     //std::cout<<"4: "<<std::endl;

   if(jets.size()>=3){
     hCounter->Fill(13,weight);
     hNjetGoodIncZX->Fill(3., weight);
     ZMass_Zinc3jet->Fill(p.M(),weight);
     ThirdJetPt_Z3jet->Fill(jets.at(2).Pt(),weight);
     ThirdJetY_Z3jet->Fill(jets.at(2).Eta(),weight);
   }
   
  //std::cout<<"5: "<<std::endl;

   if(jets.size()>=4){
     hCounter->Fill(14,weight);
     hNjetGoodIncZX->Fill(4., weight);
     ZMass_Zinc4jet->Fill(p.M(),weight);
     FourthJetPt_Z4jet->Fill(jets.at(3).Pt(),weight);
     FourthJetY_Z4jet->Fill(jets.at(3).Eta(),weight);
   }
  //std::cout<<"6: "<<std::endl;

   if(jets.size()>=5){
     hCounter->Fill(15,weight);
     hNjetGoodIncZX->Fill(5., weight);
     ZMass_Zinc5jet->Fill(p.M(),weight);
     FifthJetPt_Z5jet->Fill(jets.at(4).Pt(),weight);
     FifthJetY_Z5jet->Fill(jets.at(4).Eta(),weight);
   }
   
   if(jets.size()>=6){
     hCounter->Fill(16,weight);
     hNjetGoodIncZX->Fill(6., weight);
     ZMass_Zinc6jet->Fill(p.M(),weight);
     SixthJetPt_Z6jet->Fill(jets.at(5).Pt(),weight);
     SixthJetY_Z6jet->Fill(jets.at(5).Eta(),weight);
     }
}


// ------------ method called once each job just before starting event loop  ------------
void 
SychZJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SychZJet::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
SychZJet::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SychZJet::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SychZJet::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SychZJet::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SychZJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SychZJet);


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
