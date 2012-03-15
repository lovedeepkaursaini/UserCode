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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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

namespace ResultsOnlyzJet{
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

class ResultsOnlyZJet : public edm::EDAnalyzer {
   public:
      explicit ResultsOnlyZJet(const edm::ParameterSet&);
      ~ResultsOnlyZJet();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
bool isData;
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


TH1D* hNjetGoodZX ;


TH1D*   hCounter ;
TH1D*   hWeights ;

TH1D*   FirstJetPt_Z1jet;
TH1D*   SecondJetPt_Z2jet;
TH1D*   ThirdJetPt_Z3jet;
TH1D*   FourthJetPt_Z4jet;
TH1D*   FifthJetPt_Z5jet;
TH1D*   SixthJetPt_Z6jet;
  
TH1D*   FirstJetY_Z1jet;
TH1D*   SecondJetY_Z2jet;
TH1D*   ThirdJetY_Z3jet;
TH1D*   FourthJetY_Z4jet;
TH1D*  FifthJetY_Z5jet;//prettyHistogram(fs,"FifthJetY_Z5jet","Fifth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
 TH1D*  SixthJetY_Z6jet;//prettyHistogram(fs,"SixthJetY_Z6jet","#geq 6th Jets (Z+Jets)","#eta (Jet) ",9,0.,3.);
  
TH1D* hFillnPV;//= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPV","nPV","nPV",50,0.,50.);
TH1D* hFillnPVexc1j;//= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc1j","nPV, N_{jets} == 1)","nPV",50,0.,50.);
TH1D* hFillnPVexc2j;//= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc2j","nPV, N_{jets} == 2)","nPV",50,0.,50.);
TH1D* hFillnPVexc3j;//= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc3j","nPV, N_{jets} == 3)","nPV",50,0.,50.);
TH1D* hFillnPVexc4j;//= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc4j","nPV, N_{jets} == 4)","nPV",50,0.,50.);
TH1D* hFillnPVexc5j;//= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc5j","nPV, N_{jets} == 5)","nPV",50,0.,50.);
TH1D* hFillnPVexc6j;//= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc6j","nPV, N_{jets} == 6)","nPV",50,0.,50.);


  TH1D* ZMassOSMWexc0j;//prettyHistogram(fs,"ZMassOSMWexc0j","Z Invariant Mass (N_{jets} == 0)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc0j_EB;//prettyHistogram(fs,"ZMassOSMWexc0j_EB","Z Invariant Mass (EB, N_{jets} == 0)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc0j_EE;//prettyHistogram(fs,"ZMassOSMWexc0j_EE","Z Invariant Mass (EE, N_{jets} == 0)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc0j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc0j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 0)","M_{ee} [GeV]",20,71.,111.);

  TH1D* ZMassOSMWexc1j;//prettyHistogram(fs,"ZMassOSMWexc1j","Z Invariant Mass (N_{jets} == 1)","M_{ee} [GeV]",20,71.,111.);
TH1D*   ZMassOSMWexc1j_EB;//prettyHistogram(fs,"ZMassOSMWexc1j_EB","Z Invariant Mass (EB, N_{jets} == 1)","M_{ee} [GeV]",20,71.,111.);
TH1D*   ZMassOSMWexc1j_EE;//prettyHistogram(fs,"ZMassOSMWexc1j_EE","Z Invariant Mass (EE, N_{jets} == 1)","M_{ee} [GeV]",20,71.,111.);
TH1D*   ZMassOSMWexc1j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc1j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 1)","M_{ee} [GeV]",20,71.,111.);

  TH1D* ZMassOSMWexc2j;//prettyHistogram(fs,"ZMassOSMWexc2j","Z Invariant Mass (N_{jets} == 2)","M_{ee} [GeV]",20,71.,111.);
    TH1D* ZMassOSMWexc2j_EB;//prettyHistogram(fs,"ZMassOSMWexc2j_EB","Z Invariant Mass (EB, N_{jets} == 2)","M_{ee} [GeV]",20,71.,111.);
    TH1D* ZMassOSMWexc2j_EE;//prettyHistogram(fs,"ZMassOSMWexc2j_EE","Z Invariant Mass (EE, N_{jets} == 2)","M_{ee} [GeV]",20,71.,111.);
    TH1D* ZMassOSMWexc2j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc2j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 2)","M_{ee} [GeV]",20,71.,111.);

    TH1D* ZMassOSMWexc3j;//prettyHistogram(fs,"ZMassOSMWexc3j","Z Invariant Mass (N_{jets} == 3)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc3j_EB;//prettyHistogram(fs,"ZMassOSMWexc3j_EB","Z Invariant Mass (EB, N_{jets} == 3)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc3j_EE;//prettyHistogram(fs,"ZMassOSMWexc3j_EE","Z Invariant Mass (EE, N_{jets} == 3)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc3j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc3j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 3)","M_{ee} [GeV]",20,71.,111.);

 TH1D*  ZMassOSMWexc4j;//prettyHistogram(fs,"ZMassOSMWexc4j","Z Invariant Mass (N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc4j_EB;//prettyHistogram(fs,"ZMassOSMWexc4j_EB","Z Invariant Mass (EB, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc4j_EE;//prettyHistogram(fs,"ZMassOSMWexc4j_EE","Z Invariant Mass (EE, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc4j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc4j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);

 TH1D*  ZMassOSMWexc5j;//prettyHistogram(fs,"ZMassOSMWexc4j","Z Invariant Mass (N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc5j_EB;//prettyHistogram(fs,"ZMassOSMWexc4j_EB","Z Invariant Mass (EB, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc5j_EE;//prettyHistogram(fs,"ZMassOSMWexc4j_EE","Z Invariant Mass (EE, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc5j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc4j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);

 TH1D*  ZMassOSMWexc6j;//prettyHistogram(fs,"ZMassOSMWexc4j","Z Invariant Mass (N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc6j_EB;//prettyHistogram(fs,"ZMassOSMWexc4j_EB","Z Invariant Mass (EB, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc6j_EE;//prettyHistogram(fs,"ZMassOSMWexc4j_EE","Z Invariant Mass (EE, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc6j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc4j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);


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
ResultsOnlyZJet::ResultsOnlyZJet(const edm::ParameterSet& iConfig):
     isData(iConfig.getParameter<bool>("isData"))
{

Double_t Fall2011[50] = {
    0.003388501,    0.010357558,    0.024724258,    0.042348605,    0.058279812,    0.068851751,    0.072914824,    0.071579609,    0.066811668,    0.060672356,
    0.054528356,    0.04919354,    0.044886042,    0.041341896,    0.0384679,    0.035871463,    0.03341952,    0.030915649,    0.028395374,    0.025798107,
    0.023237445,    0.020602754,    0.0180688,    0.015559693,    0.013211063,    0.010964293,    0.008920993,    0.007080504,    0.005499239,    0.004187022,    0.003096474,    0.002237361,    0.001566428,    0.001074149,    0.000721755,    0.000470838,    0.00030268,    0.000184665,    0.000112883,    6.74043E-05,    3.82178E-05,    2.22847E-05,    1.20933E-05,    6.96173E-06,    3.4689E-06,    1.96172E-06,    8.49283E-07,    5.02393E-07,    2.15311E-07,    9.56938E-08  };


double Data[50]={2.44414e+06,1.25579e+06,3.42957e+06,5.18618e+07,2.54054e+08,4.36586e+08,4.86031e+08,4.63551e+08,4.18993e+08,3.84891e+08,3.65724e+08,3.41505e+08,3.22711e+08,3.0886e+08,2.87693e+08,2.5129e+08,1.99438e+08,1.40551e+08,8.66577e+07,4.6234e+07,2.12058e+07,8.37396e+06,2.88178e+06,882886,246537,63900.8,15571.6,3628.24,840.61,211.248,66.7507,32.1445,26.92,33.3738,47.1181,66.9628,92.522,123.311,158.278,195.606,232.736,266.603,294.026,312.195,319.142,314.095,297.617,271.501,238.455,1.93299e+07};

 for( int i=0; i<50; ++i) {
      TrueDist2011_.push_back(Data[i]);
      WLumi_.push_back(Fall2011[i]);
   }
 LumiWeights_ = edm::LumiReWeighting(WLumi_,TrueDist2011_);


   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  
 hNjetGoodZX = ResultsOnlyzJet::prettyHistogram(fs,"h_Njet_GoodZX","N_{Jets}(Exc) [Z+Jets]","NJets_{Exc}",10,0,10);



  hCounter = ResultsOnlyzJet::prettyHistogram(fs,"hCounter","Counter","EventCounter",16,0.5,16.5);
hCounter->GetXaxis()->SetBinLabel(1,"Trig");
hCounter->GetXaxis()->SetBinLabel(2,"EleSel");
hCounter->GetXaxis()->SetBinLabel(3,"MW");
hCounter->GetXaxis()->SetBinLabel(4,"nJ==0");
hCounter->GetXaxis()->SetBinLabel(5,"nJ==1");
hCounter->GetXaxis()->SetBinLabel(6,"nJ==2");
hCounter->GetXaxis()->SetBinLabel(7,"nJ==3");
hCounter->GetXaxis()->SetBinLabel(8,"nJ==4");
hCounter->GetXaxis()->SetBinLabel(9,"nJ==5");
hCounter->GetXaxis()->SetBinLabel(10,"nJ==6");
hCounter->GetXaxis()->SetBinLabel(11,"-");



  hWeights = ResultsOnlyzJet::prettyHistogram(fs,"hWeights","Weight","Weight",50,0.,10.);
  //binning jetpT
  float jpt[10] = {30,50,70,90,120,150,180,240,300,400};
  float jpt2[9] = {30,50,70,90,120,150,200,270,400};
  float jpt3[6] = {30,50,70,90,120,400};
  // float jpt4[6] = {30,50,70,90,120,400};
  FirstJetPt_Z1jet=ResultsOnlyzJet::prettyHistogram(fs,"FirstJetPt_Z1jet","First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
  SecondJetPt_Z2jet=ResultsOnlyzJet::prettyHistogram(fs,"SecondJetPt_Z2jet","Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
  ThirdJetPt_Z3jet=ResultsOnlyzJet::prettyHistogram(fs,"ThirdJetPt_Z3jet","Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  FourthJetPt_Z4jet=ResultsOnlyzJet::prettyHistogram(fs,"FourthJetPt_Z4jet","Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  FifthJetPt_Z5jet=ResultsOnlyzJet::prettyHistogram(fs,"FifthJetPt_Z5jet","Fifth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  SixthJetPt_Z6jet=ResultsOnlyzJet::prettyHistogram(fs,"SixthJetPt_Z6jet","#geq 6th Jets (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  
  FirstJetY_Z1jet=ResultsOnlyzJet::prettyHistogram(fs,"FirstJetY_Z1jet","First Leading Jet (Z+Jets)","#eta (Jet)",9,0.,3.);
  SecondJetY_Z2jet=ResultsOnlyzJet::prettyHistogram(fs,"SecondJetY_Z2jet","Second Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  ThirdJetY_Z3jet=ResultsOnlyzJet::prettyHistogram(fs,"ThirdJetY_Z3jet","Third Leading (Z+Jets)","#eta (Jet) ",9,0.,3.);
  FourthJetY_Z4jet=ResultsOnlyzJet::prettyHistogram(fs,"FourthJetY_Z4jet","Fourth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  FifthJetY_Z5jet=ResultsOnlyzJet::prettyHistogram(fs,"FifthJetY_Z5jet","Fifth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  SixthJetY_Z6jet=ResultsOnlyzJet::prettyHistogram(fs,"SixthJetY_Z6jet","#geq 6th Jets (Z+Jets)","#eta (Jet) ",9,0.,3.);
  
hFillnPV= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPV","nPV","nPV",50,0.,50.);
hFillnPVexc1j= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc1j","nPV, N_{jets} == 1)","nPV",50,0.,50.);
hFillnPVexc2j= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc2j","nPV, N_{jets} == 2)","nPV",50,0.,50.);
hFillnPVexc3j= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc3j","nPV, N_{jets} == 3)","nPV",50,0.,50.);
hFillnPVexc4j= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc4j","nPV, N_{jets} == 4)","nPV",50,0.,50.);
hFillnPVexc5j= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc5j","nPV, N_{jets} == 5)","nPV",50,0.,50.);
hFillnPVexc6j= ResultsOnlyzJet::prettyHistogram(fs,"hFillnPVexc6j","nPV, N_{jets} == 6)","nPV",50,0.,50.);


  ZMassOSMWexc0j=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc0j","Z Invariant Mass (N_{jets} == 0)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc0j_EB=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc0j_EB","Z Invariant Mass (EB, N_{jets} == 0)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc0j_EE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc0j_EE","Z Invariant Mass (EE, N_{jets} == 0)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc0j_EBEE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc0j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 0)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc1j=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc1j","Z Invariant Mass (N_{jets} == 1)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc1j_EB=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc1j_EB","Z Invariant Mass (EB, N_{jets} == 1)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc1j_EE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc1j_EE","Z Invariant Mass (EE, N_{jets} == 1)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc1j_EBEE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc1j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 1)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc2j=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc2j","Z Invariant Mass (N_{jets} == 2)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc2j_EB=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc2j_EB","Z Invariant Mass (EB, N_{jets} == 2)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc2j_EE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc2j_EE","Z Invariant Mass (EE, N_{jets} == 2)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc2j_EBEE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc2j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 2)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc3j=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc3j","Z Invariant Mass (N_{jets} == 3)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc3j_EB=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc3j_EB","Z Invariant Mass (EB, N_{jets} == 3)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc3j_EE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc3j_EE","Z Invariant Mass (EE, N_{jets} == 3)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc3j_EBEE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc3j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 3)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc4j=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc4j","Z Invariant Mass (N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc4j_EB=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc4j_EB","Z Invariant Mass (EB, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc4j_EE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc4j_EE","Z Invariant Mass (EE, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc4j_EBEE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc4j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 4)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc5j=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc5j","Z Invariant Mass (N_{jets} == 5)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc5j_EB=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc5j_EB","Z Invariant Mass (EB, N_{jets} == 5)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc5j_EE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc5j_EE","Z Invariant Mass (EE, N_{jets} == 5)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc5j_EBEE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc5j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 5)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc6j=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc6j","Z Invariant Mass (N_{jets} == 6)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc6j_EB=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc6j_EB","Z Invariant Mass (EB, N_{jets} == 6)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc6j_EE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc6j_EE","Z Invariant Mass (EE, N_{jets} == 6)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc6j_EBEE=ResultsOnlyzJet::prettyHistogram(fs,"ZMassOSMWexc6j_EBEE","Z Invariant Mass (EBEE, N_{jets} == 6)","M_{ee} [GeV]",20,71.,111.);



}


ResultsOnlyZJet::~ResultsOnlyZJet()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ResultsOnlyZJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
double MyWeight;
if(isData) MyWeight = 1.;
else 
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

 MyWeight = LumiWeights_.weight( npT );
  hWeights->Fill(MyWeight);
  //cout<<MyWeight<<endl;
  //  double weight=1.;
}
  double weight=MyWeight;
  hCounter->Fill(1,weight);
  using namespace edm;
  bool usePatElec = 1;

  int nPV=0;

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);
  for(unsigned int ind=0;ind<recVtxs->size();ind++) {
    if (!((*recVtxs)[ind].isFake()) && ((*recVtxs)[ind].ndof()>4)
        && (fabs((*recVtxs)[ind].z())<=15.0) &&
        ((*recVtxs)[ind].position().Rho()<=2.0) ) {
      nPV++;
     }
  }

  hFillnPV->Fill(nPV,weight);
  edm::Handle<pat::ElectronCollection> patElecHandle;
  if(usePatElec){
    if(not iEvent.getByLabel("selectedPatElectronsPFlow",patElecHandle)){
      std::cout<<"FATAL EXCEPTION: "<<"patElectrons Not Found: "<<std::endl;
      return;
    }
  }

  pat::Electron* e1Ptr=0; 
  pat::Electron* e2Ptr=0; 
  
  pat::ElectronCollection patEleColl(*(patElecHandle.product()));
  if(usePatElec){
    std::sort(patEleColl.begin(),patEleColl.end(),PtGreater());
    //std::cout<<"Found: "<<patEleColl.size()<<std::endl;
    if(patEleColl.size()<2)return;
    
    e1Ptr = new pat::Electron(patEleColl[0]);
    e2Ptr = new pat::Electron(patEleColl[1]);
  }
  
  pat::Electron e1(*e1Ptr);
  pat::Electron e2(*e2Ptr);
  
  double q1 = e1.threeCharge();
  double pt1     = e1.pt();
  double eta1     = e1.superCluster()->eta();
  double phi1     = e1.phi();
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
  } else if(e2.isEE()){
    e2PassId = sieie2<0.03&& fabs(delphi2)<0.03&& fabs(detain2)<0.007&& hoe2 <0.15 ;
  }
  
  bool bothPassId = e1PassId&&e2PassId;
  
  double CombRelpfIso1=  (  e1.chargedHadronIso() +  e1.neutralHadronIso() + e1.photonIso()) /pt1;
  double CombRelpfIso2=  (  e2.chargedHadronIso() +  e2.neutralHadronIso() + e2.photonIso()) /pt2;
  //cout<<CombRelpfIso1<<endl; 
  bool bothPassPfIso = CombRelpfIso1  < 0.2 && CombRelpfIso2 < 0.2;

  bool conv1=0;
  bool conv2=0;
  if(e1.isEB() && !(fabs(dcot1)<0.02&& fabs(dist1)<0.02)&& nmshits1<=0) conv1=1;
  else  if(e1.isEE() && !(fabs(dcot1)<0.02&& fabs(dist1)<0.02)&& nmshits1<=0) conv1=1;
  if(e2.isEB() && !(fabs(dcot2)<0.02&& fabs(dist2)<0.02)&& nmshits2<=0) conv2=1;
  else  if(e2.isEE() && !(fabs(dcot2)<0.02&& fabs(dist2)<0.02)&& nmshits2<=0) conv2=1;
  bool bothPassConv = conv1 && conv2;
  
  //Apply Cuts and Count the Events.
  if(!(bothPassPtEta && bothPassId && bothPassPfIso && bothPassConv))return;
  hCounter->Fill(2,weight);
 
  bool OSign = (q1*q2)<0;
  bool MWanalysis = (fabs(p.M()-91)<20);

  if(!(OSign && MWanalysis))return;
  hCounter->Fill(3,weight);

  edm::Handle<std::vector<pat::Jet> > patJetHandle;
  if(not iEvent.getByLabel("selectedPatJetsPFlow",patJetHandle)){
    std::cout<<"PAT Jets Not Found" <<std::endl;
    return;
  }
  
  const std::vector<pat::Jet>*  patJets = patJetHandle.product();
  std::vector<pat::Jet>::const_iterator jet =patJets->begin();
  
  std::vector<TLorentzVector> jets;
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
      double delR_e1 = ResultsOnlyzJet::radius(jetEta,jetPhi,eta1,phi1);
      double delR_e2 = ResultsOnlyzJet::radius(jetEta,jetPhi,eta2,phi2);
      if((delR_e1<0.3)||(delR_e2<0.3))continue;
      TLorentzVector* theJet = new TLorentzVector();
      theJet->SetPxPyPzE(jet->px(),jet->py(),jet->pz(),jet->energy());
      jets.push_back(*theJet);
  }
  
  hNjetGoodZX->Fill(jets.size(),weight);


   if(jets.size()==0)	
      {
        hCounter->Fill(4,weight);
        ZMassOSMWexc0j->Fill(p.M(),weight);
        if(e1.isEB() && e2.isEB())ZMassOSMWexc0j_EB->Fill(p.M(),weight);
        else if(e1.isEE() && e2.isEE())ZMassOSMWexc0j_EE->Fill(p.M(),weight);
        else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc0j_EBEE->Fill(p.M(),weight);
      }

   else   if(jets.size()==1){
     hCounter->Fill(5,weight);
     hFillnPVexc1j->Fill(nPV,weight);
     FirstJetPt_Z1jet->Fill(jets.at(0).Pt(),weight);
     FirstJetY_Z1jet->Fill(jets.at(0).Eta(),weight);
        ZMassOSMWexc1j->Fill(p.M(),weight);
        if(e1.isEB() && e2.isEB())ZMassOSMWexc1j_EB->Fill(p.M(),weight);
        else if(e1.isEE() && e2.isEE())ZMassOSMWexc1j_EE->Fill(p.M(),weight);
        else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc1j_EBEE->Fill(p.M(),weight);
   }
  
   else   if(jets.size()==2){
     hCounter->Fill(6,weight);
     hFillnPVexc2j->Fill(nPV,weight);
     SecondJetPt_Z2jet->Fill(jets.at(1).Pt(),weight);
     SecondJetY_Z2jet->Fill(jets.at(1).Eta(),weight);
        ZMassOSMWexc2j->Fill(p.M(),weight);
        if(e1.isEB() && e2.isEB())ZMassOSMWexc2j_EB->Fill(p.M(),weight);
        else if(e1.isEE() && e2.isEE())ZMassOSMWexc2j_EE->Fill(p.M(),weight);
        else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc2j_EBEE->Fill(p.M(),weight);
   }

else   if(jets.size()==3){
     hCounter->Fill(7,weight);
     hFillnPVexc3j->Fill(nPV,weight);
     ThirdJetPt_Z3jet->Fill(jets.at(2).Pt(),weight);
     ThirdJetY_Z3jet->Fill(jets.at(2).Eta(),weight);
        ZMassOSMWexc3j->Fill(p.M(),weight);
        if(e1.isEB() && e2.isEB())ZMassOSMWexc3j_EB->Fill(p.M(),weight);
        else if(e1.isEE() && e2.isEE())ZMassOSMWexc3j_EE->Fill(p.M(),weight);
        else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc3j_EBEE->Fill(p.M(),weight);
   }
   
else   if(jets.size()==4){
     hCounter->Fill(8,weight);
     hFillnPVexc4j->Fill(nPV,weight);
     FourthJetPt_Z4jet->Fill(jets.at(3).Pt(),weight);
     FourthJetY_Z4jet->Fill(jets.at(3).Eta(),weight);
        ZMassOSMWexc4j->Fill(p.M(),weight);
        if(e1.isEB() && e2.isEB())ZMassOSMWexc4j_EB->Fill(p.M(),weight);
        else if(e1.isEE() && e2.isEE())ZMassOSMWexc4j_EE->Fill(p.M(),weight);
        else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc4j_EBEE->Fill(p.M(),weight);
   }

else   if(jets.size()==5){
     hCounter->Fill(9,weight);
     hFillnPVexc5j->Fill(nPV,weight);
     FifthJetPt_Z5jet->Fill(jets.at(4).Pt(),weight);
     FifthJetY_Z5jet->Fill(jets.at(4).Eta(),weight);
        ZMassOSMWexc5j->Fill(p.M(),weight);
        if(e1.isEB() && e2.isEB())ZMassOSMWexc5j_EB->Fill(p.M(),weight);
        else if(e1.isEE() && e2.isEE())ZMassOSMWexc5j_EE->Fill(p.M(),weight);
        else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc5j_EBEE->Fill(p.M(),weight);
   }
   
else   if(jets.size()==6){
     hCounter->Fill(10,weight);
     hFillnPVexc6j->Fill(nPV,weight);
     SixthJetPt_Z6jet->Fill(jets.at(5).Pt(),weight);
     SixthJetY_Z6jet->Fill(jets.at(5).Eta(),weight);
        ZMassOSMWexc6j->Fill(p.M(),weight);
        if(e1.isEB() && e2.isEB())ZMassOSMWexc6j_EB->Fill(p.M(),weight);
        else if(e1.isEE() && e2.isEE())ZMassOSMWexc6j_EE->Fill(p.M(),weight);
        else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc6j_EBEE->Fill(p.M(),weight);
     }

}


// ------------ method called once each job just before starting event loop  ------------
void 
ResultsOnlyZJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ResultsOnlyZJet::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ResultsOnlyZJet::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ResultsOnlyZJet::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ResultsOnlyZJet::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ResultsOnlyZJet::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ResultsOnlyZJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ResultsOnlyZJet);


