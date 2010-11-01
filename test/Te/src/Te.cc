// -*- C++ -*-
//
// Package:    Te
// Class:      Te
// 
/**\class Te Te.cc test/Te/src/Te.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  AnilSingh
//         Created:  Sun Oct 24 14:19:17 IST 2010
// $Id$
//
//

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h" 
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class Te : public edm::EDAnalyzer {
   public:
      explicit Te(const edm::ParameterSet&);
      ~Te();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
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
Te::Te(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


Te::~Te()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
Te::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
 
  edm::Handle<edm::TriggerResults> trigResults;

  edm::InputTag trigTag("TriggerResults::HLT");
  if (not iEvent.getByLabel(trigTag, trigResults)) {
    std::cout << ">>> TRIGGER collection does not exist !!!\n";
    return;
  }

  const edm::TriggerNames & trigNames = iEvent.triggerNames(*trigResults);

  for (unsigned int i=0; i<trigResults->size(); i++)
  {
    std::string trigName = trigNames.triggerName(i);
    std::cout<<trigName<<std::endl;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
Te::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Te::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(Te);
