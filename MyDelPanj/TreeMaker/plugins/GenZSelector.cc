
// system include files
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

class GenZSelector : public edm::EDFilter {
   public:
      explicit GenZSelector(const edm::ParameterSet&);
      ~GenZSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  //  TH1D* hZPtAll_;
  // TH1D* hZPtFid_;
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
GenZSelector::GenZSelector(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


GenZSelector::~GenZSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
GenZSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  if(not iEvent.getByLabel("genParticles", genParticleHandle))
    {
      std::cout<<
	"GenAnalyzer: Generator Level Information not found\n"
              <<std::endl;
    }

  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;
  std::vector<double> en;
 
  const reco::GenParticleCollection* genColl= &(*genParticleHandle);
  reco::GenParticleCollection::const_iterator geni = genColl->begin();
  for(; geni!=genColl->end();geni++){
    reco::GenParticle gen = *geni;
    //Look out for the GenMuons/Electrons
    if(gen.status()==1){
      double pid = fabs(gen.pdgId());
      
      if(!(pid==11))continue;
      if(fabs(gen.eta())>2.5)continue;
      if(fabs(gen.px())<10)continue;
      px.push_back(gen.px());
      py.push_back(gen.py());
      pz.push_back(gen.pz());
      en.push_back(gen.energy());
    }
}
    if(px.size()<2)return false;
    
      return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenZSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenZSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
GenZSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
GenZSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
GenZSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
GenZSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenZSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenZSelector);
