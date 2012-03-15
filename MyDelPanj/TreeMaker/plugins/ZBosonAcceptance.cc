// -*- C++ -*-
//
// Package:    ZBosonAcceptance
// Class:      ZBosonAcceptance
// 
/**\class ZBosonAcceptance ZBosonAcceptance.cc DelPanj/ZBosonAcceptance/src/ZBosonAcceptance.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lovedeep Kaur (Panjab U)
//         Created:  Tue Dec 27 10:37:50 CST 2011
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

#include <memory>
#include <vector>
#include <string>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "TLorentzVector.h"

//
// class declaration
//

class ZBosonAcceptance : public edm::EDAnalyzer {
   public:
      explicit ZBosonAcceptance(const edm::ParameterSet&);
      ~ZBosonAcceptance();

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
  TH1D* hZPtBefore_;
  TH1D* hZPtAfter_;

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
ZBosonAcceptance::ZBosonAcceptance(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
edm::Service<TFileService> fs;
hZPtBefore_ = fs->make<TH1D>("hZPtBefore_","hZPtBefore_",20,0,200);
hZPtAfter_  = fs->make<TH1D>("hZPtAfter_","hZPtAfter_",20,0,200);
}


ZBosonAcceptance::~ZBosonAcceptance()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZBosonAcceptance::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  std::cout<<"1\n";
   using namespace edm;
   edm::Handle<reco::GenParticleCollection> genParticleHandle;
   if(not iEvent.getByLabel("genParticles", genParticleHandle))
     {
       std::cout<<
        "GenAnalyzer: Generator Level Information not found\n"
		<<std::endl;
     }
    std::vector<TLorentzVector*> ele;

   const reco::GenParticleCollection* genColl= &(*genParticleHandle);
   reco::GenParticleCollection::const_iterator geni = genColl->begin();
   for(; geni!=genColl->end();geni++){
     reco::GenParticle gen = *geni;
     //Look out for the GenMuons/Electrons                                                                                                                    
     if(gen.status()==1){
       double pid = fabs(gen.pdgId());
       if(!(pid==11))continue;
       TLorentzVector* e1 = new TLorentzVector();
       e1->SetPtEtaPhiM(gen.pt(),gen.eta(),gen.phi(),gen.mass()) ;
       ele.push_back(e1);      
     }
   }
   //   std::cout<<"2\n";

     if(ele.size()<2)return;

     TLorentzVector zBefore(*(ele[0])+*(ele[1]));
     //std::cout<<"213\n";

     if (( zBefore.M()<60)||(zBefore.M()>120) )return;
     hZPtBefore_->Fill(zBefore.Pt());
     //std::cout<<"23\n";
//     if (( zBefore.M()<75)||(zBefore.M()>105) )return;
     if(!(((ele[0])->Pt()>20)&&((ele[1])->Pt()>20))) return;
     //std::cout<<"3\n";

     double eta1_ = ele[0]->Eta();
     bool isBarrel1 = fabs(eta1_)<1.4442;
     bool isEndCap1 = fabs(eta1_)>1.566 && fabs(eta1_)<2.5;

     double eta2_ = ele[1]->Eta();
     bool isBarrel2 = fabs(eta2_)<1.4442;
     bool isEndCap2 = fabs(eta2_)>1.566 && fabs(eta2_)<2.5;
   
     if(!((isBarrel1||isEndCap1) && (isBarrel2||isEndCap2)))return;
     //std::cout<<"4\n";

     hZPtAfter_->Fill(zBefore.Pt());
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZBosonAcceptance::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZBosonAcceptance::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ZBosonAcceptance::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZBosonAcceptance::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZBosonAcceptance::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZBosonAcceptance::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZBosonAcceptance::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZBosonAcceptance);
