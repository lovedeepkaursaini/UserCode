// -*- C++ -*-
//
// Package:    PileupWeightProducer
// Class:      PileupWeightProducer
// 
/**\class PileupWeightProducer PileupWeightProducer.cc 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ricardo Vasquez Sierra,6 R-025,+41227672274,
//         Created:  Mon Nov 21 15:05:26 CET 2011
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
#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include <vector>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//
// class declaration
//

class PileupWeightProducer : public edm::EDProducer {
   public:
      explicit PileupWeightProducer(const edm::ParameterSet&);
      ~PileupWeightProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  bool firsttime_;
  std::string pileupMC_;
  std::string pileupData_;
  edm::LumiReWeighting LumiWeights_;
  edm::Lumi3DReWeighting LumiWeightsNominal_;
  edm::Lumi3DReWeighting LumiWeightsUp_;
  edm::Lumi3DReWeighting LumiWeightsDown_;
std::vector< float > Data2011_;
std::vector<float> Fall2011_;
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
PileupWeightProducer::PileupWeightProducer(const edm::ParameterSet& iConfig)
{
//   firsttime_= iConfig.existsAs<bool>("FirstTime") ? iConfig.getParameter<bool>("FirstTime") : true ;
//   pileupMC_ = iConfig.existsAs<std::string>("PileupMCFile") ? iConfig.getParameter<std::string>("PileupMCFile") : "PUMC_dist.root" ;
//   pileupData_ = iConfig.existsAs<std::string>("PileupDataFile") ? iConfig.getParameter<std::string>("PileupDataFile") : "PUData_dist.root" ;

  firsttime_ =  iConfig.getUntrackedParameter<bool>("FirstTime");

  //register your products
  
  produces<std::vector<float> >( "pileupWeights" ).setBranchAlias( "pileupWeights" );
  
  

  if ( firsttime_ )
    {
//      std::cout<< " Initializing with the following files MC: " << pileupMC_ << " data: " << pileupData_ << std::endl;
Double_t Fall2011[50] = {
   0.003388501,    0.010357558,    0.024724258,    0.042348605,    0.058279812,    0.068851751,    0.072914824,    0.071579609,    0.066811668,    0.060672356,    0.054528356,    0.04919354,    0.044886042,    0.041341896,    0.0384679,    0.035871463,    0.03341952,    0.030915649,    0.028395374,   0.025798107,    0.023237445,    0.020602754,    0.0180688,    0.015559693,    0.013211063,    0.010964293,    0.008920993,    0.007080504,   0.005499239,    0.004187022,    0.003096474,    0.002237361,    0.001566428,    0.001074149,    0.000721755,   0.000470838,    0.00030268,    0.000184665,   0.000112883,    6.74043E-05,    3.82178E-05,    2.22847E-05,    1.20933E-05,    6.96173E-06,    3.4689E-06,    1.96172E-06,    8.49283E-07,   5.02393E-07,    2.15311E-07,    9.56938E-08  };

double Data[50]={2.90319e+06,1.09975e+06,6.52617e+06,1.11487e+08,3.87573e+08,5.46905e+08,5.4739e+08,4.93312e+08,4.37187e+08,4.03341e+08,3.72317e+08,3.46568e+08,3.27998e+08,3.0256e+08,2.59765e+08,1.9859e+08,1.31263e+08,7.38331e+07,3.50021e+07,1.39562e+07,4.73169e+06,1.40154e+06,377317,96352,24194.2,6086.59,1528.98,382.112,101.523,39.0797,32.526,43.5174,63.6907,91.3667,125.885,166.022,209.478,252.849,291.962,322.505,340.794,344.501,333.146,308.193,272.744,230.904,187.004,144.883,107.381,2.50586e+07};



 for( int i=0; i<50; ++i) {
      Data2011_.push_back(Data[i]);
      Fall2011_.push_back(Fall2011[i]);
   }
 LumiWeights_ = edm::LumiReWeighting(Fall2011_,Data2011_);
}
//      LumiWeights_ = edm::LumiReWeighting(pileupMC_, pileupData_, "pileup", "pileup");}
/*      LumiWeightsNominal_.weight3D_set( pileupMC_, pileupData_, "pileup", "pileup");
      LumiWeightsUp_.weight3D_set( pileupMC_, pileupData_, "pileup", "pileup");
      LumiWeightsDown_.weight3D_set( pileupMC_, pileupData_, "pileup", "pileup");

      LumiWeightsNominal_.weight3D_init(1.0);
      LumiWeightsUp_.weight3D_init(1.08);
      LumiWeightsDown_.weight3D_init(0.92);
    }
  else 
    {
      std::cout<< " Initializing with Weight3D.root " << std::endl; 
      LumiWeightsNominal_.weight3D_init("Weight3D.root");
      LumiWeightsUp_.weight3D_init("Weight3DscaleUp.root");
      LumiWeightsDown_.weight3D_init("Weight3DscaleDown.root");
    }
*/


  
}


PileupWeightProducer::~PileupWeightProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called to produce the data  ------------
void
PileupWeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   std::auto_ptr<std::vector<float> > pileupWeights( new std::vector<float> );

/*   edm::EventBase* iEventB = dynamic_cast<edm::EventBase*>(&iEvent);
  double MyWeight = LumiWeights_.weight( (*iEventB) );

   double nominalWeight3D = LumiWeightsNominal_.weight3D( (*iEventB) );
   double weight3DUp = LumiWeightsUp_.weight3D( (*iEventB) );
   double weight3DDown = LumiWeightsDown_.weight3D( (*iEventB) );

   pileupWeights->push_back( MyWeight );

   pileupWeights->push_back( nominalWeight3D );
   pileupWeights->push_back( weight3DUp );
   pileupWeights->push_back( weight3DDown );

   iEvent.put(pileupWeights, "pileupWeights");
*/   
///yar kuch kerna paina...mera apna style he theek e
Handle<std::vector< PileupSummaryInfo > >  PupInfo;
iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

std::vector<PileupSummaryInfo>::const_iterator PVI;

float Tnpv = -1;
for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

   int BX = PVI->getBunchCrossing();

   if(BX == 0) { 
     Tnpv = PVI->getTrueNumInteractions();
     continue;
   }

}
double MyWeight = LumiWeights_.weight( Tnpv );
   pileupWeights->push_back( MyWeight );
   iEvent.put(pileupWeights, "pileupWeights");

}

// ------------ method called once each job just before starting event loop  ------------
void 
PileupWeightProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PileupWeightProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
PileupWeightProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PileupWeightProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PileupWeightProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PileupWeightProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PileupWeightProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PileupWeightProducer);
