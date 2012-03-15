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

class check : public edm::EDAnalyzer {
   public:
      explicit check(const edm::ParameterSet&);
      ~check();

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
check::check(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  z0Hist_ = fs->make<TH1D>("hZ0Mass","hZ0Mass",10,0,10.);
}


check::~check()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
check::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
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
  
  
  //std::cout<<"4"<<std::endl;
  
  if(usePatElec){
    pat::ElectronCollection patEleColl(*(patElecHandle.product()));
    std::sort(patEleColl.begin(),patEleColl.end(),PtGreater());
  z0Hist_->Fill(patEleColl.size());
  }
  

 }


// ------------ method called once each job just before starting event loop  ------------
void 
check::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
check::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
check::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
check::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
check::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
check::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
check::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(check);


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
