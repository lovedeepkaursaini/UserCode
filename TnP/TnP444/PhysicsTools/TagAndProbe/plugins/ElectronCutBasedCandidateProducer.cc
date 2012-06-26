#include "PhysicsTools/TagAndProbe/interface/ElectronCutBasedCandidateProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

ElectronCutBasedCandidateProducer::ElectronCutBasedCandidateProducer(const edm::ParameterSet 
&params)
{

  const edm::InputTag allelectrons("gsfElectrons");
  electronCollection_ = 
    params.getUntrackedParameter<edm::InputTag>("ReferenceElectronCollection", 
						allelectrons);
  inputTagIsoDepElectrons_ = params.getParameter< std::vector<edm::InputTag> >("IsoDepElectron");    
  inputTagIsoValElectronsPFId_   = params.getParameter< std::vector<edm::InputTag> >("IsoValElectronPF");
  conversionsInputTag_    = params.getParameter<edm::InputTag>("conversionsInputTag");
  beamSpotInputTag_       = params.getParameter<edm::InputTag>("beamSpotInputTag");
  rhoIsoInputTag_          = params.getParameter<edm::InputTag>("rhoIsoInputTag");
  primaryVertexInputTag_  = params.getParameter<edm::InputTag>("primaryVertexInputTag");  

  nameIDBOOL_ = params.getParameter<std::string>( "nameIDBOOL" );

  produces< edm::PtrVector<reco::GsfElectron> >();
  produces< edm::RefToBaseVector<reco::GsfElectron> >();
}




ElectronCutBasedCandidateProducer::~ElectronCutBasedCandidateProducer()
{

}


//
// member functions
//


// ------------ method called to produce the data  ------------

void ElectronCutBasedCandidateProducer::produce(edm::Event &event, 
			      const edm::EventSetup &eventSetup)
{
   // Create the output collection
  std::auto_ptr< edm::RefToBaseVector<reco::GsfElectron> > 
    outColRef( new edm::RefToBaseVector<reco::GsfElectron> );
  std::auto_ptr< edm::PtrVector<reco::GsfElectron> > 
    outColPtr( new edm::PtrVector<reco::GsfElectron> );

  // Read electrons
  edm::Handle<edm::View<reco::GsfElectron> > velectrons;
  event.getByLabel(electronCollection_, velectrons);
  edm::Handle<reco::GsfElectronCollection> electrons;
  event.getByLabel("gsfElectrons", electrons);

  const edm::PtrVector<reco::GsfElectron>& ptrVect = velectrons->ptrVector();
  const edm::RefToBaseVector<reco::GsfElectron>& refs = velectrons->refVector();
  unsigned int counter=0;


/////// Pileup density "rho" for lepton isolation subtraction /////
//  double RhoForLeptonIsolation = -999999.9;
//  edm::Handle<double> rhoLepIso;
//  const edm::InputTag eventrhoLepIso("kt6PFJets", "rho");
//  event.getByLabel(eventrhoLepIso, rhoLepIso);
//  if( *rhoLepIso == *rhoLepIso) RhoForLeptonIsolation = *rhoLepIso;


  
  // get the iso deposits. 3 (charged hadrons, photons, neutral hadrons)
  unsigned nTypes=3;
  IsoDepositMaps electronIsoDep(nTypes);

  for (size_t j = 0; j<inputTagIsoDepElectrons_.size(); ++j) {
    event.getByLabel(inputTagIsoDepElectrons_[j], electronIsoDep[j]);
  }
  IsoDepositVals electronIsoValPFId(nTypes);
  const IsoDepositVals * electronIsoVals = &electronIsoValPFId;

  for (size_t j = 0; j<inputTagIsoValElectronsPFId_.size(); ++j) {
    event.getByLabel(inputTagIsoValElectronsPFId_[j], electronIsoValPFId[j]);
  }

    // conversions
    edm::Handle<reco::ConversionCollection> conversions_h;
    event.getByLabel(conversionsInputTag_, conversions_h);
    // beam spot
    edm::Handle<reco::BeamSpot> beamspot_h;
    event.getByLabel(beamSpotInputTag_, beamspot_h);
    const reco::BeamSpot &beamSpot = *(beamspot_h.product());

    // vertices
    edm::Handle<reco::VertexCollection> vtx_h;
    event.getByLabel(primaryVertexInputTag_, vtx_h);
    // rho for isolation
    edm::Handle<double> rhoIso_h;
    event.getByLabel(rhoIsoInputTag_, rhoIso_h);
    double rhoIso = *(rhoIso_h.product());


  for(edm::View<reco::GsfElectron>::const_iterator gsfIt = velectrons->begin();
      gsfIt != velectrons->end(); ++gsfIt, ++counter){


  // Loop over electrons
    unsigned nele=electrons->size();
   for(unsigned iele=0; iele<nele;++iele) {
//  unsigned iele = 0;    

      reco::GsfElectronRef myElectronRef(electrons,iele);
   if (fabs(gsfIt->p4().Pt() - myElectronRef->pt()) < 0.0001 ){
//      double isoval = 0.0;
//      isoval = gsfIt->dr03TkSumPt() + gsfIt->dr03HcalTowerSumEt();
//      isoval += 
//gsfIt->hadronicOverEm()*gsfIt->superCluster()->energy()/fabs(cosh(gsfIt->superCluster()->eta()));
//      if (gsfIt->isEB()) isoval += std::max(gsfIt->dr03EcalRecHitSumEt() - 1.0,0.0);
//      if (gsfIt->isEE()) isoval += gsfIt->dr03EcalRecHitSumEt();
//            isoval = std::max(isoval-0.3*0.3*pi*RhoForLeptonIsolation ,0.0);
    
      
 
      double iso_ch =  (*(*electronIsoVals)[0])[myElectronRef];
      double iso_em = (*(*electronIsoVals)[1])[myElectronRef];
      double iso_nh = (*(*electronIsoVals)[2])[myElectronRef];
//std::cout<<iso_ch<<"\t"<<iso_em<<"\t"<<iso_nh<<std::endl;
      bool veto       = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO, myElectronRef, 
                      conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
      bool loose      = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, myElectronRef, 
                      conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
      bool medium     = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, myElectronRef, 
                      conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);
      bool tight      = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, myElectronRef, 
                      conversions_h, beamSpot, vtx_h, iso_ch, iso_em, iso_nh, rhoIso);

//      std::cout<<"medium: "<<medium<<" run: "<<event.run()<<" event: "<<event.id().event()<<std::endl;
        // eop/fbrem cuts for extra tight ID
      bool fbremeopin = EgammaCutBasedEleId::PassEoverPCuts(myElectronRef);

        // cuts to match tight trigger requirements
      bool trigtight = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, myElectronRef);

        // for 2011 WP70 trigger
        bool trigwp70 = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERWP70, myElectronRef);
 
//      isoval = charged+photon+neutral;
//      isoval /= gsfIt->p4().Pt();
            
//      if( (gsfIt->isEB() && isoval < combrelisoCutEB_) || (gsfIt->isEE() && isoval < 
//combrelisoCutEE_) ){
	std::string nam = nameIDBOOL_;
	bool pass = false;
	if(nam=="veto") pass = veto;
	else if (nam=="loose")pass= loose;
	else if(nam=="medium")pass= medium;
	else if(nam=="tight")pass= tight;
	if (pass) {
	outColRef->push_back( refs[counter] );
	outColPtr->push_back( ptrVect[counter]  );
	} // end if loop
    }
  } // end candidate loop
}
  event.put(outColRef);
  event.put(outColPtr);
}




// ------ method called once each job just before starting event loop  ---



void ElectronCutBasedCandidateProducer::beginJob() {}



void ElectronCutBasedCandidateProducer::endJob() {}



//define this as a plug-in
DEFINE_FWK_MODULE( ElectronCutBasedCandidateProducer );



