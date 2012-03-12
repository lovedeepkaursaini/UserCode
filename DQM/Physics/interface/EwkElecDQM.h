#ifndef EwkElecDQM_H
#define EwkElecDQM_H

/** \class EwkElecDQM
 *
 *  DQM offline fwor workEWK Electrons
 *
 */

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DQM/Physics/interface/dqmGlobal.h"
#include "DQM/Physics/interface/dqmElectron.h"
#include "DQM/Physics/interface/dqmJet.h"
#include "DQM/Physics/interface/dqmBoson.h"
#include "DQM/Physics/interface/dqmWorkpoint.h"


class DQMStore;
class MonitorElement;
class EwkElecDQM : public edm::EDAnalyzer {
public:
  EwkElecDQM (const edm::ParameterSet &);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);
  bool triggerDecision(const edm::Event&);
  double calcDeltaPhi(double phi1, double phi2);


private:
  edm::InputTag trigTag_; 
  edm::InputTag elecTag_;
  edm::InputTag metTag_;
  edm::InputTag jetTag_;

  const std::vector<std::string> tags_;
  const std::vector<int> occurences_;
  const std::vector<std::string> vetoes_;

  double eJetMinX_;
  double eMinX_; 
  unsigned int nGoodElectrons;

  ewkdqm::gWorkpoint gwp_;
  ewkdqm::eWorkpoint ewp_;
  ewkdqm::jWorkpoint jwp_;

  DQMStore* theDbe;

  ewkdqm::GlobalHist* globalHist_;
  
  // ewkdqm::GlobalHist* globalHistBefore_;
  // ewkdqm::GlobalHist* globalHistAfter_;

  ewkdqm::ElectronHist* elecHistBefore_;
  ewkdqm::ElectronHist* elecHistAfter_;

  ewkdqm::BosonHist* zBosonBeforeCuts_;
  ewkdqm::BosonHist* zBosonAfterCuts_;
  
  ewkdqm::BosonHist* wBosonBeforeCuts_;
  ewkdqm::BosonHist* wBosonAfterCuts_;
  
  ewkdqm::BosonHist* zBosonEE_;
  ewkdqm::BosonHist* zBosonEB_;
  ewkdqm::BosonHist* zBosonBB_;
  
  ewkdqm::BosonHist* zBoson1jet_;
  ewkdqm::BosonHist* zBoson2jet_;
  ewkdqm::BosonHist* zBoson3jet_;

  ewkdqm::JetHist* jetHist_;
  ewkdqm::JetHist* jetHistLead_;
  ewkdqm::JetHist* jetHist2ndLead_;
  ewkdqm::JetHist* jetHist3rdLead_;
  ewkdqm::JetHist* jetHist4thLead_;



  bool isValidHltConfig_;
  HLTConfigProvider  hltConfigProvider_;

  double rec_;
  double recpteta_;
  double recptetaid_;
  double recptetaidiso_;
  double recptetaidisonjet_;
  double recptetaidisonjetmet_;
  double recptetaidisonjetmettrig_;


  };
#endif
