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
#include "DelPanj/TreeMaker/interface/myLib.h"
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

namespace DeepzJet{
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

class DeepZJet : public edm::EDAnalyzer {
   public:
      explicit DeepZJet(const edm::ParameterSet&);
      ~DeepZJet();

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

TH1D*   ZMass_Zexc0jet;
TH1D*   ZMass_Zexc1jet;
TH1D*   ZMass_Zexc2jet;
TH1D*   ZMass_Zexc3jet;
TH1D*   ZMass_Zexc4jet;
TH1D*   ZMass_Zexc5jet;
TH1D*   ZMass_Zexc6jet;

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
  
TH1D* hFillnPV;//= DeepzJet::prettyHistogram(fs,"hFillnPV","nPV","nPV",50,0.,50.);
TH1D* hFillnPVexc1j;//= DeepzJet::prettyHistogram(fs,"hFillnPVexc1j","nPV, N_{jets} #eq 1)","nPV",50,0.,50.);
TH1D* hFillnPVexc2j;//= DeepzJet::prettyHistogram(fs,"hFillnPVexc2j","nPV, N_{jets} #eq 2)","nPV",50,0.,50.);
TH1D* hFillnPVexc3j;//= DeepzJet::prettyHistogram(fs,"hFillnPVexc3j","nPV, N_{jets} #eq 3)","nPV",50,0.,50.);
TH1D* hFillnPVexc4j;//= DeepzJet::prettyHistogram(fs,"hFillnPVexc4j","nPV, N_{jets} #eq 4)","nPV",50,0.,50.);
TH1D* hFillnPVexc5j;//= DeepzJet::prettyHistogram(fs,"hFillnPVexc5j","nPV, N_{jets} #eq 5)","nPV",50,0.,50.);
TH1D* hFillnPVexc6j;//= DeepzJet::prettyHistogram(fs,"hFillnPVexc6j","nPV, N_{jets} #eq 6)","nPV",50,0.,50.);

  TH1D* ZMassOSMW;//prettyHistogram(fs,"ZMassOSMW","Z Invariant Mass (OS, MW:71-111, inc)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWdeep;//prettyHistogram(fs,"ZMassOSMWdeep","Z Invariant Mass (OS, MW:50-200, inc)","M_{ee} [GeV]",150,50.,200.);
  TH1D* ZMassSSMW;//prettyHistogram(fs,"ZMassSSMW","Z Invariant Mass (SS, MW:71-111, inc)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassSSMWdeep;//prettyHistogram(fs,"ZMassSSMWdeep","Z Invariant Mass (SS, MW:50-200, inc)","M_{ee} [GeV]",150,50.,200.);

  TH1D* ZMassOSMWexc0j;//prettyHistogram(fs,"ZMassOSMWexc0j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWdeepexc0j;//prettyHistogram(fs,"ZMassOSMWdeepexc0j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 0)","M_{ee} [GeV]",150,50.,200.);
  TH1D* ZMassSSMWexc0j;//prettyHistogram(fs,"ZMassSSMWexc0j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassSSMWdeepexc0j;//prettyHistogram(fs,"ZMassSSMWdeepexc0j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 0)","M_{ee} [GeV]",150,50.,200.);
  TH1D* ZMassOSMWexc0j_EB;//prettyHistogram(fs,"ZMassOSMWexc0j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc0j_EE;//prettyHistogram(fs,"ZMassOSMWexc0j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc0j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc0j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);

  TH1D* ZMassOSMWexc1j;//prettyHistogram(fs,"ZMassOSMWexc1j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWdeepexc1j;//prettyHistogram(fs,"ZMassOSMWdeepexc1j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 1)","M_{ee} [GeV]",150,50.,200.);
TH1D*   ZMassSSMWexc1j;//prettyHistogram(fs,"ZMassSSMWexc1j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
TH1D*   ZMassSSMWdeepexc1j;//prettyHistogram(fs,"ZMassSSMWdeepexc1j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 1)","M_{ee} [GeV]",150,50.,200.);
TH1D*   ZMassOSMWexc1j_EB;//prettyHistogram(fs,"ZMassOSMWexc1j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
TH1D*   ZMassOSMWexc1j_EE;//prettyHistogram(fs,"ZMassOSMWexc1j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
TH1D*   ZMassOSMWexc1j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc1j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);

  TH1D* ZMassOSMWexc2j;//prettyHistogram(fs,"ZMassOSMWexc2j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  TH1D*   ZMassOSMWdeepexc2j;//prettyHistogram(fs,"ZMassOSMWdeepexc2j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 2)","M_{ee} [GeV]",150,50.,200.);
    TH1D* ZMassSSMWexc2j;//prettyHistogram(fs,"ZMassSSMWexc2j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
    TH1D* ZMassSSMWdeepexc2j;//prettyHistogram(fs,"ZMassSSMWdeepexc2j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 2)","M_{ee} [GeV]",150,50.,200.);
    TH1D* ZMassOSMWexc2j_EB;//prettyHistogram(fs,"ZMassOSMWexc2j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
    TH1D* ZMassOSMWexc2j_EE;//prettyHistogram(fs,"ZMassOSMWexc2j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
    TH1D* ZMassOSMWexc2j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc2j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);

    TH1D* ZMassOSMWexc3j;//prettyHistogram(fs,"ZMassOSMWexc3j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
    TH1D* ZMassOSMWdeepexc3j;//prettyHistogram(fs,"ZMassOSMWdeepexc3j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 3)","M_{ee} [GeV]",150,50.,200.);
  TH1D* ZMassSSMWexc3j;//prettyHistogram(fs,"ZMassSSMWexc3j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassSSMWdeepexc3j;//prettyHistogram(fs,"ZMassSSMWdeepexc3j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 3)","M_{ee} [GeV]",150,50.,200.);
  TH1D* ZMassOSMWexc3j_EB;//prettyHistogram(fs,"ZMassOSMWexc3j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc3j_EE;//prettyHistogram(fs,"ZMassOSMWexc3j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc3j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc3j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);

 TH1D*  ZMassOSMWexc4j;//prettyHistogram(fs,"ZMassOSMWexc4j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
 TH1D*  ZMassOSMWdeepexc4j;//prettyHistogram(fs,"ZMassOSMWdeepexc4j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 4)","M_{ee} [GeV]",150,50.,200.);
  TH1D* ZMassSSMWexc4j;//prettyHistogram(fs,"ZMassSSMWexc4j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassSSMWdeepexc4j;//prettyHistogram(fs,"ZMassSSMWdeepexc4j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 4)","M_{ee} [GeV]",150,50.,200.);
  TH1D* ZMassOSMWexc4j_EB;//prettyHistogram(fs,"ZMassOSMWexc4j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc4j_EE;//prettyHistogram(fs,"ZMassOSMWexc4j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
  TH1D* ZMassOSMWexc4j_EBEE;//prettyHistogram(fs,"ZMassOSMWexc4j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);

TH1D*   hPUFirstJetPt_Z1jet;
TH1D*   hPUSecondJetPt_Z2jet;//prettyHistogram(fs,"hPUSecondJetPt_Z2jet","hPU, Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
TH1D*   hPUThirdJetPt_Z3jet;//prettyHistogram(fs,"hPUThirdJetPt_Z3jet","hPU, Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
TH1D*   hPUFourthJetPt_Z4jet;//prettyHistogram(fs,"hPUFourthJetPt_Z4jet","hPU, Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  
TH1D*   hPUFirstJetY_Z1jet;//prettyHistogram(fs,"hPUFirstJetY_Z1jet","hPU, First Leading Jet (Z+Jets)","#eta (Jet)",9,0.,3.);
TH1D*   hPUSecondJetY_Z2jet;//prettyHistogram(fs,"hPUSecondJetY_Z2jet","hPU, Second Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
TH1D*   hPUThirdJetY_Z3jet;//prettyHistogram(fs,"hPUThirdJetY_Z3jet","hPU, Third Leading (Z+Jets)","#eta (Jet) ",9,0.,3.);
TH1D*  hPUFourthJetY_Z4jet;//prettyHistogram(fs,"hPUFourthJetY_Z4jet","hPU, Fourth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
   
TH1D*   hPUZMassOSMWexc1j;//prettyHistogram(fs,"hPUZMassOSMWexc1j","hPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
TH1D*   hPUZMassOSMWexc2j;//prettyHistogram(fs,"hPUZMassOSMWexc2j","hPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
TH1D*   hPUZMassOSMWexc3j;//prettyHistogram(fs,"hPUZMassOSMWexc3j","hPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
TH1D*   hPUZMassOSMWexc4j;//prettyHistogram(fs,"hPUZMassOSMWexc4j","hPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);

TH1D*  lPUFirstJetPt_Z1jet;//prettyHistogram(fs,"lPUFirstJetPt_Z1jet","lPU, First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
 TH1D*  lPUSecondJetPt_Z2jet;//prettyHistogram(fs,"lPUSecondJetPt_Z2jet","lPU, Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
 TH1D*  lPUThirdJetPt_Z3jet;//prettyHistogram(fs,"lPUThirdJetPt_Z3jet","lPU, Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
 TH1D*  lPUFourthJetPt_Z4jet;//prettyHistogram(fs,"lPUFourthJetPt_Z4jet","lPU, Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  
 TH1D*  lPUFirstJetY_Z1jet;//prettyHistogram(fs,"lPUFirstJetY_Z1jet","lPU, First Leading Jet (Z+Jets)","#eta (Jet)",9,0.,3.);
 TH1D*  lPUSecondJetY_Z2jet;//prettyHistogram(fs,"lPUSecondJetY_Z2jet","lPU, Second Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
 TH1D*  lPUThirdJetY_Z3jet;//prettyHistogram(fs,"lPUThirdJetY_Z3jet","lPU, Third Leading (Z+Jets)","#eta (Jet) ",9,0.,3.);
 TH1D* lPUFourthJetY_Z4jet;//prettyHistogram(fs,"lPUFourthJetY_Z4jet","lPU, Fourth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
   
  TH1D* lPUZMassOSMWexc1j;//prettyHistogram(fs,"lPUZMassOSMWexc1j","lPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  TH1D* lPUZMassOSMWexc2j;//prettyHistogram(fs,"lPUZMassOSMWexc2j","lPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  TH1D* lPUZMassOSMWexc3j;//prettyHistogram(fs,"lPUZMassOSMWexc3j","lPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  TH1D* lPUZMassOSMWexc4j;//prettyHistogram(fs,"lPUZMassOSMWexc4j","lPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);

  TH1D* vhPUFirstJetPt_Z1jet;//prettyHistogram(fs,"vhPUFirstJetPt_Z1jet","vhPU, First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
  TH1D* vhPUSecondJetPt_Z2jet;//prettyHistogram(fs,"vhPUSecondJetPt_Z2jet","vhPU, Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
 TH1D*  vhPUThirdJetPt_Z3jet;//prettyHistogram(fs,"vhPUThirdJetPt_Z3jet","vhPU, Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
 TH1D*  vhPUFourthJetPt_Z4jet;//prettyHistogram(fs,"vhPUFourthJetPt_Z4jet","vhPU, Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  
  TH1D* vhPUFirstJetY_Z1jet;//prettyHistogram(fs,"vhPUFirstJetY_Z1jet","vhPU, First Leading Jet (Z+Jets)","#eta (Jet)",9,0.,3.);
  TH1D* vhPUSecondJetY_Z2jet;//prettyHistogram(fs,"vhPUSecondJetY_Z2jet","vhPU, Second Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  TH1D* vhPUThirdJetY_Z3jet;//prettyHistogram(fs,"vhPUThirdJetY_Z3jet","vhPU, Third Leading (Z+Jets)","#eta (Jet) ",9,0.,3.);
 TH1D* vhPUFourthJetY_Z4jet;//prettyHistogram(fs,"vhPUFourthJetY_Z4jet","vhPU, Fourth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
   
  TH1D* vhPUZMassOSMWexc1j;//prettyHistogram(fs,"vhPUZMassOSMWexc1j","vhPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  TH1D* vhPUZMassOSMWexc2j;//prettyHistogram(fs,"vhPUZMassOSMWexc2j","vhPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  TH1D* vhPUZMassOSMWexc3j;//prettyHistogram(fs,"vhPUZMassOSMWexc3j","vhPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  TH1D* vhPUZMassOSMWexc4j;//prettyHistogram(fs,"vhPUZMassOSMWexc4j","vhPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
 


TH1D* cosThetaStarBoostToCM_exc1j ;//= new cosThetaStarBoostToCM_exc1j();
TH1D* cosThetaStarZBoostToCM_exc1j ; 
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
DeepZJet::DeepZJet(const edm::ParameterSet& iConfig):
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
  
 hNjetGoodZX = DeepzJet::prettyHistogram(fs,"h_Njet_GoodZX","N_{Jets}(Exc) [Z+Jets]","NJets_{Exc}",10,0,10);

  ZMass_Zexc0jet=DeepzJet::prettyHistogram(fs,"ZMass_Zexc0jet","Z Invariant Mass (N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);
  ZMass_Zexc1jet=DeepzJet::prettyHistogram(fs,"ZMass_Zexc1jet","Z Invariant Mass (N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  ZMass_Zexc2jet=DeepzJet::prettyHistogram(fs,"ZMass_Zexc2jet","Z Invariant Mass (N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  ZMass_Zexc3jet=DeepzJet::prettyHistogram(fs,"ZMass_Zexc3jet","Z Invariant Mass (N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  ZMass_Zexc4jet=DeepzJet::prettyHistogram(fs,"ZMass_Zexc4jet","Z Invariant Mass (N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
  ZMass_Zexc5jet=DeepzJet::prettyHistogram(fs,"ZMass_Zexc5jet","Z Invariant Mass (N_{jets} #eq 5)","M_{ee} [GeV]",20,71.,111.);
  ZMass_Zexc6jet=DeepzJet::prettyHistogram(fs,"ZMass_Zexc6jet","Z Invariant Mass (N_{jets} #eq 6)","M_{ee} [GeV]",20,71.,111.);


  hCounter = DeepzJet::prettyHistogram(fs,"hCounter","Counter","EventCounter",16,0.5,16.5);
hCounter->GetXaxis()->SetBinLabel(1,"Trig");
hCounter->GetXaxis()->SetBinLabel(2,"EleSel");
hCounter->GetXaxis()->SetBinLabel(3,"OS,MW");
hCounter->GetXaxis()->SetBinLabel(4,"nJ==0");
hCounter->GetXaxis()->SetBinLabel(5,"nJ==1");
hCounter->GetXaxis()->SetBinLabel(6,"nJ==2");
hCounter->GetXaxis()->SetBinLabel(7,"nJ==3");
hCounter->GetXaxis()->SetBinLabel(8,"nJ==4");
hCounter->GetXaxis()->SetBinLabel(9,"nJ==5");
hCounter->GetXaxis()->SetBinLabel(10,"nJ==6");
hCounter->GetXaxis()->SetBinLabel(11,"-");



  hWeights = DeepzJet::prettyHistogram(fs,"hWeights","Weight","Weight",50,0.,10.);
  //binning jetpT
  float jpt[10] = {30,50,70,90,120,150,180,240,300,400};
  float jpt2[9] = {30,50,70,90,120,150,200,270,400};
  float jpt3[6] = {30,50,70,90,120,400};
  // float jpt4[6] = {30,50,70,90,120,400};
  FirstJetPt_Z1jet=DeepzJet::prettyHistogram(fs,"FirstJetPt_Z1jet","First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
  SecondJetPt_Z2jet=DeepzJet::prettyHistogram(fs,"SecondJetPt_Z2jet","Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
  ThirdJetPt_Z3jet=DeepzJet::prettyHistogram(fs,"ThirdJetPt_Z3jet","Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  FourthJetPt_Z4jet=DeepzJet::prettyHistogram(fs,"FourthJetPt_Z4jet","Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  FifthJetPt_Z5jet=DeepzJet::prettyHistogram(fs,"FifthJetPt_Z5jet","Fifth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  SixthJetPt_Z6jet=DeepzJet::prettyHistogram(fs,"SixthJetPt_Z6jet","#geq 6th Jets (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  
  FirstJetY_Z1jet=DeepzJet::prettyHistogram(fs,"FirstJetY_Z1jet","First Leading Jet (Z+Jets)","#eta (Jet)",9,0.,3.);
  SecondJetY_Z2jet=DeepzJet::prettyHistogram(fs,"SecondJetY_Z2jet","Second Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  ThirdJetY_Z3jet=DeepzJet::prettyHistogram(fs,"ThirdJetY_Z3jet","Third Leading (Z+Jets)","#eta (Jet) ",9,0.,3.);
  FourthJetY_Z4jet=DeepzJet::prettyHistogram(fs,"FourthJetY_Z4jet","Fourth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  FifthJetY_Z5jet=DeepzJet::prettyHistogram(fs,"FifthJetY_Z5jet","Fifth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  SixthJetY_Z6jet=DeepzJet::prettyHistogram(fs,"SixthJetY_Z6jet","#geq 6th Jets (Z+Jets)","#eta (Jet) ",9,0.,3.);
  
hFillnPV= DeepzJet::prettyHistogram(fs,"hFillnPV","nPV","nPV",50,0.,50.);
hFillnPVexc1j= DeepzJet::prettyHistogram(fs,"hFillnPVexc1j","nPV, N_{jets} #eq 1)","nPV",50,0.,50.);
hFillnPVexc2j= DeepzJet::prettyHistogram(fs,"hFillnPVexc2j","nPV, N_{jets} #eq 2)","nPV",50,0.,50.);
hFillnPVexc3j= DeepzJet::prettyHistogram(fs,"hFillnPVexc3j","nPV, N_{jets} #eq 3)","nPV",50,0.,50.);
hFillnPVexc4j= DeepzJet::prettyHistogram(fs,"hFillnPVexc4j","nPV, N_{jets} #eq 4)","nPV",50,0.,50.);
hFillnPVexc5j= DeepzJet::prettyHistogram(fs,"hFillnPVexc5j","nPV, N_{jets} #eq 5)","nPV",50,0.,50.);
hFillnPVexc6j= DeepzJet::prettyHistogram(fs,"hFillnPVexc6j","nPV, N_{jets} #eq 6)","nPV",50,0.,50.);

  ZMassOSMW=DeepzJet::prettyHistogram(fs,"ZMassOSMW","Z Invariant Mass (OS, MW:71-111, inc)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWdeep=DeepzJet::prettyHistogram(fs,"ZMassOSMWdeep","Z Invariant Mass (OS, MW:50-200, inc)","M_{ee} [GeV]",150,50.,200.);
  ZMassSSMW=DeepzJet::prettyHistogram(fs,"ZMassSSMW","Z Invariant Mass (SS, MW:71-111, inc)","M_{ee} [GeV]",20,71.,111.);
  ZMassSSMWdeep=DeepzJet::prettyHistogram(fs,"ZMassSSMWdeep","Z Invariant Mass (SS, MW:50-200, inc)","M_{ee} [GeV]",150,50.,200.);

  ZMassOSMWexc0j=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc0j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWdeepexc0j=DeepzJet::prettyHistogram(fs,"ZMassOSMWdeepexc0j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 0)","M_{ee} [GeV]",150,50.,200.);
  ZMassSSMWexc0j=DeepzJet::prettyHistogram(fs,"ZMassSSMWexc0j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);
  ZMassSSMWdeepexc0j=DeepzJet::prettyHistogram(fs,"ZMassSSMWdeepexc0j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 0)","M_{ee} [GeV]",150,50.,200.);
  ZMassOSMWexc0j_EB=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc0j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc0j_EE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc0j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc0j_EBEE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc0j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 0)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc1j=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc1j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWdeepexc1j=DeepzJet::prettyHistogram(fs,"ZMassOSMWdeepexc1j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 1)","M_{ee} [GeV]",150,50.,200.);
  ZMassSSMWexc1j=DeepzJet::prettyHistogram(fs,"ZMassSSMWexc1j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  ZMassSSMWdeepexc1j=DeepzJet::prettyHistogram(fs,"ZMassSSMWdeepexc1j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 1)","M_{ee} [GeV]",150,50.,200.);
  ZMassOSMWexc1j_EB=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc1j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc1j_EE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc1j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc1j_EBEE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc1j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc2j=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc2j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWdeepexc2j=DeepzJet::prettyHistogram(fs,"ZMassOSMWdeepexc2j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 2)","M_{ee} [GeV]",150,50.,200.);
  ZMassSSMWexc2j=DeepzJet::prettyHistogram(fs,"ZMassSSMWexc2j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  ZMassSSMWdeepexc2j=DeepzJet::prettyHistogram(fs,"ZMassSSMWdeepexc2j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 2)","M_{ee} [GeV]",150,50.,200.);
  ZMassOSMWexc2j_EB=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc2j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc2j_EE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc2j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc2j_EBEE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc2j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc3j=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc3j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWdeepexc3j=DeepzJet::prettyHistogram(fs,"ZMassOSMWdeepexc3j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 3)","M_{ee} [GeV]",150,50.,200.);
  ZMassSSMWexc3j=DeepzJet::prettyHistogram(fs,"ZMassSSMWexc3j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  ZMassSSMWdeepexc3j=DeepzJet::prettyHistogram(fs,"ZMassSSMWdeepexc3j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 3)","M_{ee} [GeV]",150,50.,200.);
  ZMassOSMWexc3j_EB=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc3j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc3j_EE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc3j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc3j_EBEE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc3j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);

  ZMassOSMWexc4j=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc4j","Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWdeepexc4j=DeepzJet::prettyHistogram(fs,"ZMassOSMWdeepexc4j","Z Invariant Mass (OS, MW:50-200, N_{jets} #eq 4)","M_{ee} [GeV]",150,50.,200.);
  ZMassSSMWexc4j=DeepzJet::prettyHistogram(fs,"ZMassSSMWexc4j","Z Invariant Mass (SS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
  ZMassSSMWdeepexc4j=DeepzJet::prettyHistogram(fs,"ZMassSSMWdeepexc4j","Z Invariant Mass (SS, MW:50-200, N_{jets} #eq 4)","M_{ee} [GeV]",150,50.,200.);
  ZMassOSMWexc4j_EB=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc4j_EB","Z Invariant Mass (EB, OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc4j_EE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc4j_EE","Z Invariant Mass (EE, OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);
  ZMassOSMWexc4j_EBEE=DeepzJet::prettyHistogram(fs,"ZMassOSMWexc4j_EBEE","Z Invariant Mass (EBEE, OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);



  hPUFirstJetPt_Z1jet=DeepzJet::prettyHistogram(fs,"hPUFirstJetPt_Z1jet","hPU, First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
  hPUSecondJetPt_Z2jet=DeepzJet::prettyHistogram(fs,"hPUSecondJetPt_Z2jet","hPU, Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
  hPUThirdJetPt_Z3jet=DeepzJet::prettyHistogram(fs,"hPUThirdJetPt_Z3jet","hPU, Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  hPUFourthJetPt_Z4jet=DeepzJet::prettyHistogram(fs,"hPUFourthJetPt_Z4jet","hPU, Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  
  hPUFirstJetY_Z1jet=DeepzJet::prettyHistogram(fs,"hPUFirstJetY_Z1jet","hPU, First Leading Jet (Z+Jets)","#eta (Jet)",9,0.,3.);
  hPUSecondJetY_Z2jet=DeepzJet::prettyHistogram(fs,"hPUSecondJetY_Z2jet","hPU, Second Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  hPUThirdJetY_Z3jet=DeepzJet::prettyHistogram(fs,"hPUThirdJetY_Z3jet","hPU, Third Leading (Z+Jets)","#eta (Jet) ",9,0.,3.);
 hPUFourthJetY_Z4jet=DeepzJet::prettyHistogram(fs,"hPUFourthJetY_Z4jet","hPU, Fourth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
   
  hPUZMassOSMWexc1j=DeepzJet::prettyHistogram(fs,"hPUZMassOSMWexc1j","hPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  hPUZMassOSMWexc2j=DeepzJet::prettyHistogram(fs,"hPUZMassOSMWexc2j","hPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  hPUZMassOSMWexc3j=DeepzJet::prettyHistogram(fs,"hPUZMassOSMWexc3j","hPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  hPUZMassOSMWexc4j=DeepzJet::prettyHistogram(fs,"hPUZMassOSMWexc4j","hPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);

 lPUFirstJetPt_Z1jet=DeepzJet::prettyHistogram(fs,"lPUFirstJetPt_Z1jet","lPU, First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
  lPUSecondJetPt_Z2jet=DeepzJet::prettyHistogram(fs,"lPUSecondJetPt_Z2jet","lPU, Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
  lPUThirdJetPt_Z3jet=DeepzJet::prettyHistogram(fs,"lPUThirdJetPt_Z3jet","lPU, Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  lPUFourthJetPt_Z4jet=DeepzJet::prettyHistogram(fs,"lPUFourthJetPt_Z4jet","lPU, Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  
  lPUFirstJetY_Z1jet=DeepzJet::prettyHistogram(fs,"lPUFirstJetY_Z1jet","lPU, First Leading Jet (Z+Jets)","#eta (Jet)",9,0.,3.);
  lPUSecondJetY_Z2jet=DeepzJet::prettyHistogram(fs,"lPUSecondJetY_Z2jet","lPU, Second Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  lPUThirdJetY_Z3jet=DeepzJet::prettyHistogram(fs,"lPUThirdJetY_Z3jet","lPU, Third Leading (Z+Jets)","#eta (Jet) ",9,0.,3.);
 lPUFourthJetY_Z4jet=DeepzJet::prettyHistogram(fs,"lPUFourthJetY_Z4jet","lPU, Fourth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
   
  lPUZMassOSMWexc1j=DeepzJet::prettyHistogram(fs,"lPUZMassOSMWexc1j","lPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  lPUZMassOSMWexc2j=DeepzJet::prettyHistogram(fs,"lPUZMassOSMWexc2j","lPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  lPUZMassOSMWexc3j=DeepzJet::prettyHistogram(fs,"lPUZMassOSMWexc3j","lPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  lPUZMassOSMWexc4j=DeepzJet::prettyHistogram(fs,"lPUZMassOSMWexc4j","lPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);

  vhPUFirstJetPt_Z1jet=DeepzJet::prettyHistogram(fs,"vhPUFirstJetPt_Z1jet","vhPU, First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
  vhPUSecondJetPt_Z2jet=DeepzJet::prettyHistogram(fs,"vhPUSecondJetPt_Z2jet","vhPU, Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
  vhPUThirdJetPt_Z3jet=DeepzJet::prettyHistogram(fs,"vhPUThirdJetPt_Z3jet","vhPU, Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  vhPUFourthJetPt_Z4jet=DeepzJet::prettyHistogram(fs,"vhPUFourthJetPt_Z4jet","vhPU, Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
  
  vhPUFirstJetY_Z1jet=DeepzJet::prettyHistogram(fs,"vhPUFirstJetY_Z1jet","vhPU, First Leading Jet (Z+Jets)","#eta (Jet)",9,0.,3.);
  vhPUSecondJetY_Z2jet=DeepzJet::prettyHistogram(fs,"vhPUSecondJetY_Z2jet","vhPU, Second Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
  vhPUThirdJetY_Z3jet=DeepzJet::prettyHistogram(fs,"vhPUThirdJetY_Z3jet","vhPU, Third Leading (Z+Jets)","#eta (Jet) ",9,0.,3.);
 vhPUFourthJetY_Z4jet=DeepzJet::prettyHistogram(fs,"vhPUFourthJetY_Z4jet","vhPU, Fourth Leading Jet (Z+Jets)","#eta (Jet) ",9,0.,3.);
   
  vhPUZMassOSMWexc1j=DeepzJet::prettyHistogram(fs,"vhPUZMassOSMWexc1j","vhPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 1)","M_{ee} [GeV]",20,71.,111.);
  vhPUZMassOSMWexc2j=DeepzJet::prettyHistogram(fs,"vhPUZMassOSMWexc2j","vhPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 2)","M_{ee} [GeV]",20,71.,111.);
  vhPUZMassOSMWexc3j=DeepzJet::prettyHistogram(fs,"vhPUZMassOSMWexc3j","vhPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 3)","M_{ee} [GeV]",20,71.,111.);
  vhPUZMassOSMWexc4j=DeepzJet::prettyHistogram(fs,"vhPUZMassOSMWexc4j","vhPU, Z Invariant Mass (OS, MW:71-111, N_{jets} #eq 4)","M_{ee} [GeV]",20,71.,111.);

 cosThetaStarBoostToCM_exc1j = DeepzJet::prettyHistogram(fs,"cosThetaStarBoostToCM_exc1j","cosThetaStarBoostToCM_exc1j","|cosThetaStar|",20,0.,1.);

 cosThetaStarZBoostToCM_exc1j = DeepzJet::prettyHistogram(fs,"cosThetaStarZBoostToCM_exc1j","cosThetaStarZBoostToCM_exc1j","|cosThetaStar (ZBoostToCM)|",20,0.,1.);

}


DeepZJet::~DeepZJet()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DeepZJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  bool SSign = (q1*q2)>0;
  bool MWanalysis = (fabs(p.M()-91)<20);
  bool MWfordeep=(p.M()>50) && (p.M()<200);

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
      //for z daughters
      double delR_e1 = 9999999999.; double delR_e2 = 9999999999.;
      if(OSign && MWanalysis) {
	delR_e1 = DeepzJet::radius(jetEta,jetPhi,eta1,phi1);
	delR_e2 = DeepzJet::radius(jetEta,jetPhi,eta2,phi2);
      }
      if((delR_e1<0.3)||(delR_e2<0.3))continue;
      TLorentzVector* theJet = new TLorentzVector();
      theJet->SetPxPyPzE(jet->px(),jet->py(),jet->pz(),jet->energy());
      jets.push_back(*theJet);
  }
  
  hNjetGoodZX->Fill(jets.size(),weight);


  bool exc0j = jets.size()==0;
  bool exc1j = jets.size()==1;
  bool exc2j = jets.size()==2;
  bool exc3j = jets.size()==3;
  bool exc4j = jets.size()==4;
  bool exc5j = jets.size()==5;

  bool lowPU= nPV<=8;
  bool highPU= nPV>8 && nPV<=15;
  bool vhighPU= nPV>16;


  if(OSign && MWanalysis) ZMassOSMW->Fill(p.M(),weight);
  if(OSign && MWfordeep) ZMassOSMWdeep->Fill(p.M(),weight);
  if(SSign && MWanalysis) ZMassSSMW->Fill(p.M(),weight);
  if(SSign && MWfordeep) ZMassSSMWdeep->Fill(p.M(),weight);

  if(exc0j){
    if(OSign && MWanalysis) 
      { 
	ZMassOSMWexc0j->Fill(p.M(),weight);
	if(e1.isEB() && e2.isEB())ZMassOSMWexc0j_EB->Fill(p.M(),weight);
	else if(e1.isEE() && e2.isEE())ZMassOSMWexc0j_EE->Fill(p.M(),weight);
	else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc0j_EBEE->Fill(p.M(),weight);
      }
    if(OSign && MWfordeep) ZMassOSMWdeepexc0j->Fill(p.M(),weight);
    if(SSign && MWanalysis) ZMassSSMWexc0j->Fill(p.M(),weight);
    if(SSign && MWfordeep) ZMassSSMWdeepexc0j->Fill(p.M(),weight);
  }

  if(exc1j){
    if(OSign && MWanalysis) 
      { 
	ZMassOSMWexc1j->Fill(p.M(),weight);
	if(e1.isEB() && e2.isEB())ZMassOSMWexc1j_EB->Fill(p.M(),weight);
	else if(e1.isEE() && e2.isEE())ZMassOSMWexc1j_EE->Fill(p.M(),weight);
	else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc1j_EBEE->Fill(p.M(),weight);
      }
    if(OSign && MWfordeep) ZMassOSMWdeepexc1j->Fill(p.M(),weight);
    if(SSign && MWanalysis) ZMassSSMWexc1j->Fill(p.M(),weight);
    if(SSign && MWfordeep) ZMassSSMWdeepexc1j->Fill(p.M(),weight);
  }

  if(exc2j){
    if(OSign && MWanalysis) 
      { 
	ZMassOSMWexc2j->Fill(p.M(),weight);
	if(e1.isEB() && e2.isEB())ZMassOSMWexc2j_EB->Fill(p.M(),weight);
	else if(e1.isEE() && e2.isEE())ZMassOSMWexc2j_EE->Fill(p.M(),weight);
	else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc2j_EBEE->Fill(p.M(),weight);
      }
    if(OSign && MWfordeep) ZMassOSMWdeepexc2j->Fill(p.M(),weight);
    if(SSign && MWanalysis) ZMassSSMWexc2j->Fill(p.M(),weight);
    if(SSign && MWfordeep) ZMassSSMWdeepexc2j->Fill(p.M(),weight);
  }

  if(exc3j){
    if(OSign && MWanalysis) 
      { 
	ZMassOSMWexc3j->Fill(p.M(),weight);
	if(e1.isEB() && e2.isEB())ZMassOSMWexc3j_EB->Fill(p.M(),weight);
	else if(e1.isEE() && e2.isEE())ZMassOSMWexc3j_EE->Fill(p.M(),weight);
	else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc3j_EBEE->Fill(p.M(),weight);
      }
    if(OSign && MWfordeep) ZMassOSMWdeepexc3j->Fill(p.M(),weight);
    if(SSign && MWanalysis) ZMassSSMWexc3j->Fill(p.M(),weight);
    if(SSign && MWfordeep) ZMassSSMWdeepexc3j->Fill(p.M(),weight);
  }

  if(exc4j){
    if(OSign && MWanalysis) 
      { 
	ZMassOSMWexc4j->Fill(p.M(),weight);
	if(e1.isEB() && e2.isEB())ZMassOSMWexc4j_EB->Fill(p.M(),weight);
	else if(e1.isEE() && e2.isEE())ZMassOSMWexc4j_EE->Fill(p.M(),weight);
	else if((e1.isEB() && e2.isEE()) || (e1.isEE() && e2.isEB()))ZMassOSMWexc4j_EBEE->Fill(p.M(),weight);
      }
    if(OSign && MWfordeep) ZMassOSMWdeepexc4j->Fill(p.M(),weight);
    if(SSign && MWanalysis) ZMassSSMWexc4j->Fill(p.M(),weight);
    if(SSign && MWfordeep) ZMassSSMWdeepexc4j->Fill(p.M(),weight);
  }


  if(!(OSign && MWanalysis))return;
  hCounter->Fill(3,weight);

  if(lowPU){
   if(jets.size()==1){
     lPUZMassOSMWexc1j->Fill(p.M(),weight);
     lPUFirstJetPt_Z1jet->Fill(jets.at(0).Pt(),weight);
     lPUFirstJetY_Z1jet->Fill(jets.at(0).Eta(),weight);
   }
   else if(jets.size()==2){
     lPUZMassOSMWexc2j->Fill(p.M(),weight);
     lPUSecondJetPt_Z2jet->Fill(jets.at(1).Pt(),weight);
     lPUSecondJetY_Z2jet->Fill(jets.at(1).Eta(),weight);
   }
   else if(jets.size()==3){
     lPUZMassOSMWexc3j->Fill(p.M(),weight);
     lPUThirdJetPt_Z3jet->Fill(jets.at(2).Pt(),weight);
     lPUThirdJetY_Z3jet->Fill(jets.at(2).Eta(),weight);
   }
   else if(jets.size()==4){
     lPUZMassOSMWexc4j->Fill(p.M(),weight);
     lPUFourthJetPt_Z4jet->Fill(jets.at(3).Pt(),weight);
     lPUFourthJetY_Z4jet->Fill(jets.at(3).Eta(),weight);
   }
  }
  else if(highPU){
   if(jets.size()==1){
     hPUZMassOSMWexc1j->Fill(p.M(),weight);
     hPUFirstJetPt_Z1jet->Fill(jets.at(0).Pt(),weight);
     hPUFirstJetY_Z1jet->Fill(jets.at(0).Eta(),weight);
   }
   else if(jets.size()==2){
     hPUZMassOSMWexc2j->Fill(p.M(),weight);
     hPUSecondJetPt_Z2jet->Fill(jets.at(1).Pt(),weight);
     hPUSecondJetY_Z2jet->Fill(jets.at(1).Eta(),weight);
   }
   else if(jets.size()==3){
     hPUZMassOSMWexc3j->Fill(p.M(),weight);
     hPUThirdJetPt_Z3jet->Fill(jets.at(2).Pt(),weight);
     hPUThirdJetY_Z3jet->Fill(jets.at(2).Eta(),weight);
   }
   else if(jets.size()==4){
     hPUZMassOSMWexc4j->Fill(p.M(),weight);
     hPUFourthJetPt_Z4jet->Fill(jets.at(3).Pt(),weight);
     hPUFourthJetY_Z4jet->Fill(jets.at(3).Eta(),weight);
   }
  }
  else if(vhighPU){
   if(jets.size()==1){
     vhPUZMassOSMWexc1j->Fill(p.M(),weight);
     vhPUFirstJetPt_Z1jet->Fill(jets.at(0).Pt(),weight);
     vhPUFirstJetY_Z1jet->Fill(jets.at(0).Eta(),weight);
   }
   else if(jets.size()==2){
     vhPUZMassOSMWexc2j->Fill(p.M(),weight);
     vhPUSecondJetPt_Z2jet->Fill(jets.at(1).Pt(),weight);
     vhPUSecondJetY_Z2jet->Fill(jets.at(1).Eta(),weight);
   }
   else if(jets.size()==3){
     vhPUZMassOSMWexc3j->Fill(p.M(),weight);
     vhPUThirdJetPt_Z3jet->Fill(jets.at(2).Pt(),weight);
     vhPUThirdJetY_Z3jet->Fill(jets.at(2).Eta(),weight);
   }
   else if(jets.size()==4){
     vhPUZMassOSMWexc4j->Fill(p.M(),weight);
     vhPUFourthJetPt_Z4jet->Fill(jets.at(3).Pt(),weight);
     vhPUFourthJetY_Z4jet->Fill(jets.at(3).Eta(),weight);
   }
  }


   if(jets.size()==0)	hCounter->Fill(4,weight);

   else   if(jets.size()==1){
     hCounter->Fill(5,weight);
     hFillnPVexc1j->Fill(nPV,weight);
     FirstJetPt_Z1jet->Fill(jets.at(0).Pt(),weight);
     FirstJetY_Z1jet->Fill(jets.at(0).Eta(),weight);
double costhetaZ = eiko::cosThetaStar_ZBoostToCM(p,jets[0]);
cosThetaStarZBoostToCM_exc1j->Fill(costhetaZ,weight);

double costheta = eiko::cosThetaStar_BoostToCM(p,jets[0]);
cosThetaStarBoostToCM_exc1j->Fill(costheta,weight);


   }
  
   else   if(jets.size()==2){
     hCounter->Fill(6,weight);
     hFillnPVexc2j->Fill(nPV,weight);
     SecondJetPt_Z2jet->Fill(jets.at(1).Pt(),weight);
     SecondJetY_Z2jet->Fill(jets.at(1).Eta(),weight);
   }

else   if(jets.size()==3){
     hCounter->Fill(7,weight);
     hFillnPVexc3j->Fill(nPV,weight);
     ThirdJetPt_Z3jet->Fill(jets.at(2).Pt(),weight);
     ThirdJetY_Z3jet->Fill(jets.at(2).Eta(),weight);
   }
   
else   if(jets.size()==4){
     hCounter->Fill(8,weight);
     hFillnPVexc4j->Fill(nPV,weight);
     FourthJetPt_Z4jet->Fill(jets.at(3).Pt(),weight);
     FourthJetY_Z4jet->Fill(jets.at(3).Eta(),weight);
   }

else   if(jets.size()==5){
     hCounter->Fill(9,weight);
     hFillnPVexc5j->Fill(nPV,weight);
     FifthJetPt_Z5jet->Fill(jets.at(4).Pt(),weight);
     FifthJetY_Z5jet->Fill(jets.at(4).Eta(),weight);
   }
   
else   if(jets.size()==6){
     hCounter->Fill(10,weight);
     hFillnPVexc6j->Fill(nPV,weight);
     SixthJetPt_Z6jet->Fill(jets.at(5).Pt(),weight);
     SixthJetY_Z6jet->Fill(jets.at(5).Eta(),weight);
     }
}


// ------------ method called once each job just before starting event loop  ------------
void 
DeepZJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DeepZJet::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
DeepZJet::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DeepZJet::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DeepZJet::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DeepZJet::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DeepZJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepZJet);


