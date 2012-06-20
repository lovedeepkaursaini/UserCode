#ifndef PhysicsTools_TagAndProbe_ElectronCutBasedCandidateProducer_h
#define PhysicsTools_TagAndProbe_ElectronCutBasedCandidateProducer_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include <iostream>
#include <string>
#include <map>

// forward declarations

class ElectronCutBasedCandidateProducer : public edm::EDProducer 
{
 public:
  explicit ElectronCutBasedCandidateProducer(const edm::ParameterSet&);
  ~ElectronCutBasedCandidateProducer();

  typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

 private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------

  edm::InputTag electronCollection_;
  std::vector<edm::InputTag> inputTagIsoDepElectrons_;
  std::vector<edm::InputTag> inputTagIsoValElectronsPFId_;   
  edm::InputTag               conversionsInputTag_;
  edm::InputTag               beamSpotInputTag_;
  edm::InputTag               rhoIsoInputTag_;
  edm::InputTag               primaryVertexInputTag_;
std::string nameIDBOOL_;
};

#endif
