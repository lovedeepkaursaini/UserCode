// -*- C++ -*-
//
// Package:    WJet
// Class:      WJet
// 
/**\class WJet WJet.cc DelPanj/WJet/src/WJet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lovedeep Kaur (Panjab U)
//         Created:  Sat Dec 31 05:10:54 CST 2011
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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Math/interface/LorentzVector.h" 
#include "DelPanj/TreeMaker/interface/eSelector.h"
#include "DelPanj/TreeMaker/interface/utils.h"
#include "DelPanj/TreeMaker/interface/eHist.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class WJet : public edm::EDAnalyzer {
   public:
      explicit WJet(const edm::ParameterSet&);
      ~WJet();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

  TH1D* e0HistB_;
  TH1D* e0HistE_;

  TH1D* e1HistB_;
  TH1D* e1HistE_;

  TH1D* e2HistB_;
  TH1D* e2HistE_;
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
WJet::WJet(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  e0HistB_ = fs->make<TH1D>("hE0PtB","hE0PtB",50,0,200);
  e0HistE_ = fs->make<TH1D>("hE0PtE","hE0PtE",50,0,200);

  e1HistB_ = fs->make<TH1D>("hE1PtB","hE1PtB",50,0,200);
  e1HistE_ = fs->make<TH1D>("hE1PtE","hE1PtE",50,0,200);

  e2HistB_ = fs->make<TH1D>("hE2PtB","hE2PtB",50,0,200);
  e2HistE_ = fs->make<TH1D>("hE2PtE","hE2PtE",50,0,200);


}


WJet::~WJet()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
WJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
 using namespace edm;
  edm::Handle<pat::ElectronCollection> patElecHandle;
  if(not iEvent.getByLabel("selectedPatElectrons",patElecHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "<<std::endl;}
//	     <<patElecLabel_<<std::endl; return false;}
  pat::ElectronCollection eleColl(*(patElecHandle.product()));
  if(eleColl.size()<1)return;
  
  pat::Electron e1 = eleColl[0];
  double pt1       = e1.pt();
  double eta1      = e1.eta();
  double tiso1     = e1.dr03TkSumPt();
  double eiso1     = e1.dr03EcalRecHitSumEt();
  double hiso1     = e1.dr03HcalTowerSumEt();
  double sieie1    = e1.sigmaIetaIeta();
  double delphi1   = e1.deltaPhiSuperClusterTrackAtVtx();
  double detain1   = e1.deltaEtaSuperClusterTrackAtVtx();
  double dcot1     = e1.convDcot();
  double dist1     = e1.convDist();
  double hoe1      = e1.hadronicOverEm();
  double nmshits1  = e1.gsfTrack().get()->trackerExpectedHitsInner().numberOfHits();
  
  bool passPtEta  = pt1 > 20 && fabs(eta1)<2.5 && !(fabs(eta1)>1.4442 && fabs(eta1)<1.566)  ;
  
  bool e1PassId = 0;  
  if(e1.isEB() && sieie1 <0.01&& fabs(delphi1)<0.06&& fabs(detain1)<0.004 && hoe1 <0.04)e1PassId=1;
  else if(e1.isEE() && sieie1 <0.03&& fabs(delphi1)<0.03&& fabs(detain1)<0.007&&hoe1 <0.15) e1PassId=1;
    
  bool e1PassIso = 0; 
  if(e1.isEB()){  e1PassIso = tiso1/pt1<0.09 && eiso1/pt1<0.07 && hiso1/pt1<0.10;}
  else if(e1.isEE()){ e1PassIso=tiso1/pt1<0.04 && eiso1/pt1<0.05 && hiso1/pt1<0.025;}
  
  bool e1PassConv = 0;
  if(e1.isEB() && !(fabs(dcot1)<=0.02&& fabs(dist1)<=0.02)&& nmshits1<=0)e1PassConv=1;
  else if(e1.isEE() && !(fabs(dcot1)<=0.02&& fabs(dist1)<=0.02)&& nmshits1<=0)e1PassConv=1;


  if(!(passPtEta))return;
  if(e1.isEB())e0HistB_->Fill(pt1);
  else if(e1.isEE())e0HistE_->Fill(pt1);

  if(!( e1PassId&&e1PassIso))return;
  if(e1.isEB())e1HistB_->Fill(pt1);
  else if(e1.isEE())e1HistE_->Fill(pt1);

  if(!(e1PassConv))return;
  if(e1.isEB())e2HistB_->Fill(pt1);
  else if(e1.isEE())e2HistE_->Fill(pt1);
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
WJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
WJet::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
WJet::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
WJet::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
WJet::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
WJet::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
WJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WJet);
