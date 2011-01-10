// -*- C++ -*-
//
// Package:    ConRej
// Class:      ConRej
// 
/**\class ConRej ConRej.cc test/ConRej/src/ConRej.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lovedeep Kaur Saini,,,
//         Created:  Wed Oct 20 16:26:05 CEST 2010
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"


//
// class declaration
//

class ConRej : public edm::EDProducer {
   public:
      explicit ConRej(const edm::ParameterSet&);
      ~ConRej();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      edm::InputTag probes_;                 
      bool isData_;
      
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
ConRej::ConRej(const edm::ParameterSet& iConfig):
  probes_(iConfig.getParameter<edm::InputTag>("probes")),
  isData_(iConfig.getUntrackedParameter<bool>("IsData"))
{
  produces<std::vector<reco::GsfElectron> >();
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
*/
   //now do what ever other initialization is needed
  
}


ConRej::~ConRej()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ConRej::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
using namespace std;
  std::auto_ptr<std::vector<reco::GsfElectron> > 
    pos(new std::vector<reco::GsfElectron>());

  double evt_bField;
  edm::Handle<reco::TrackCollection> tracks_h;
  if(!(iEvent.getByLabel("generalTracks", tracks_h))) {std::cout<<"No track Collection\n"; return;} ;

    edm::Handle<reco::GsfElectronCollection> elHandle;
    if(!(iEvent.getByLabel(probes_, elHandle)))    {std::cout<<"No ELectron Collection\n"; return;} ;
    if(isData_){
    edm::Handle<DcsStatusCollection> dcsHandle;
    iEvent.getByLabel("scalersRawToDigi", dcsHandle);

    float currentToBFieldScaleFactor = 2.09237036221512717e-04; // this is for data only!!! must do some other magic for monte-carlo   
    float current = (*dcsHandle)[0].magnetCurrent();
    evt_bField = current*currentToBFieldScaleFactor;
    }
    else
    {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
    }
  double dist = 0.0;
  double dcot = 0.0;
  double convradius = 0.0;
  double passConvRej = 0.0;

    const reco::GsfElectronCollection* electronCollection = elHandle.product();
    reco::GsfElectronCollection::const_iterator eleIt = electronCollection->begin();

      for (eleIt=electronCollection->begin(); eleIt!=electronCollection->end(); eleIt++) {
//    const reco::GsfElectronCollection* electronCollection = elHandle.product();
         ConversionFinder convFinder;
	  ConversionInfo convInfo = convFinder.getConversionInfo(*eleIt, tracks_h, evt_bField);
	  dist = convInfo.dist();
	  dcot = convInfo.dcot();
	  convradius = convInfo.radiusOfConversion();
	  if( fabs(dist)>=0.02 && fabs(dcot)>=0.02)  passConvRej = 1.0;
          if(passConvRej)pos->push_back(*eleIt);
     }      
         iEvent.put(pos);

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
ConRej::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ConRej::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ConRej);
