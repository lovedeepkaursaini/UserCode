#include "DQM/Physics/interface/EwkElecDQM.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/View.h"
#include <algorithm>  
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include <TLorentzVector.h>
using namespace edm;
using namespace std;
using namespace reco;


struct HighestPt{
   bool operator()( TLorentzVector j1, TLorentzVector j2 ) const{
    return j1.Pt() > j2.Pt() ;
   }
};




int 
NumberOfOccurences(std::string tag, std::string parent){
  if(tag==" " || parent==" ")return 0;
  int tagSize = tag.size();
  
  size_t found = parent.find(tag);
  int occurences = 0;
  
  while(found!=std::string::npos){
    unsigned int parSize = parent.size();
    size_t numFinder = found+tagSize;
    if(parent[found-1]=='_')
      if((parent[numFinder]<='9')&&(parent[numFinder]>='0'))
    	occurences++;
    if(found+tagSize>=parSize) 
      break;
    parent = parent.substr(found+tagSize,parSize);
    found = parent.find(tag);
  }
  return occurences;
}


//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



EwkElecDQM::EwkElecDQM( const ParameterSet & cfg ) :
  //Object Tags.
  trigTag_(cfg.getUntrackedParameter<edm::InputTag> ("TrigTag", edm::InputTag("TriggerResults::HLT"))),
  elecTag_(cfg.getUntrackedParameter<edm::InputTag> ("ElecTag", edm::InputTag("gsfElectrons"))), 
  metTag_(cfg.getUntrackedParameter<edm::InputTag> ("METTag", edm::InputTag("met"))),
  jetTag_(cfg.getUntrackedParameter<edm::InputTag> ("JetTag", edm::InputTag("sisCone5CaloJets"))),
  //Tags to create required pathnames.
  tags_(cfg.getUntrackedParameter<std::vector< std::string > >("TrigTags")),
  occurences_(cfg.getUntrackedParameter<std::vector<int> >("TagMulti")),
  vetoes_(cfg.getUntrackedParameter<std::vector< std::string > >("VetoTheseTags")),
  //working points: selection criteria.
  gwp_(cfg),
  ewp_(cfg),
  jwp_(cfg)
{
  isValidHltConfig_ = false;
}


//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



void EwkElecDQM::beginRun(const Run& r, const EventSetup& iSetup) {
   bool isConfigChanged = false;
   isValidHltConfig_ = hltConfigProvider_.init( r, iSetup, trigTag_.process(), isConfigChanged );
}



//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////




void EwkElecDQM::beginJob() {
  theDbe = Service<DQMStore>().operator->();
  theDbe->setCurrentFolder("Physics/EwkElecDQM");
  
  globalHist_  = new ewkdqm::GlobalHist("NOCUTS", theDbe);
  elecHistBefore_ =  new ewkdqm::ElectronHist("BEFORECUTS",theDbe);
  elecHistAfter_  =  new ewkdqm::ElectronHist("AFTERCUTS", theDbe);

  zBosonBeforeCuts_ = new ewkdqm::BosonHist("Z_BEFORECUTS", theDbe);
  zBosonAfterCuts_ = new ewkdqm::BosonHist("Z_AFTERCUTS", theDbe);
  
  wBosonBeforeCuts_ = new ewkdqm::BosonHist("W_BEFORECUTS", theDbe);
  wBosonAfterCuts_ = new ewkdqm::BosonHist("W_AFTERCUTS", theDbe);
  
  zBosonEE_ = new ewkdqm::BosonHist("Z_AFTERCUTS_EE", theDbe);
  zBosonEB_ = new ewkdqm::BosonHist("Z_AFTERCUTS_EB", theDbe);
  zBosonBB_ = new ewkdqm::BosonHist("Z_AFTERCUTS_BB", theDbe);
  
  zBoson1jet_ = new ewkdqm::BosonHist("Z_AFTERCUTS_1Jet", theDbe);
  zBoson2jet_ = new ewkdqm::BosonHist("Z_AFTERCUTS_2Jet", theDbe);
  zBoson3jet_ = new ewkdqm::BosonHist("Z_AFTERCUTS_3Jet", theDbe);
  
  jetHistLead_=      new ewkdqm::JetHist("LEAD_JETS",  theDbe);
  jetHist2ndLead_ =  new ewkdqm::JetHist("2ndLEAD_JETS",  theDbe);
  jetHist3rdLead_ =  new ewkdqm::JetHist("3rdLEAD_JETS",  theDbe);
  jetHist4thLead_ =  new ewkdqm::JetHist("4thLEAD_JETS",  theDbe);
  



   rec_=0;
   recpteta_=0;
   recptetaid_=0;
   recptetaidiso_=0;
   recptetaidisonjet_=0;
   recptetaidisonjetmet_=0;
   recptetaidisonjetmettrig_=0;
}


//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////




void EwkElecDQM::endJob() {
  double emin = gwp_.neleMinX_; 
  LogVerbatim("")<<"Evt with "<<emin<<" Elec(s) (no cut) = "<<rec_<<"  [events]";
  LogVerbatim("")<<"Evt with "<<emin<<" Elec(s) (passing pt, eta) = "<<recpteta_<<"  [events]; AbsEff(RelEff) = " <<recpteta_/rec_<<"("<<recpteta_/rec_<<")";
  LogVerbatim("")<<"Evt with "<<emin<<" Elec(s) (passing pt,eta,Id) = "<<recptetaid_<<" [events]; AbsEff(RelEff) =  " <<recptetaid_/rec_<<"("<<recptetaid_/recpteta_<<")";
  LogVerbatim("")<<"Evt with "<<emin<<" Elec(s) (passing pt,eta,Id,Iso) = "<<recptetaidiso_<<" [events] AbsEff(RelEff) =  " <<recptetaidiso_/rec_<<"("<<recptetaidiso_/recptetaid_<<")";
  LogVerbatim("")<<"Evt with "<<emin<<" Elec(s) (passing pt,eta,Id,Iso,njet) = "<<recptetaidisonjet_<<" [events] AbsEff(RelEff) =  "<<recptetaidisonjet_/rec_<<"("<<recptetaidisonjet_/recptetaidiso_<<")";
  LogVerbatim("")<<"Evt with "<<emin<<" Elec(s) (passing pt,eta,Id,Iso,njet,met) = "<<recptetaidisonjetmet_<<" [events] AbsEff(RelEff) =  "<<recptetaidisonjetmet_/rec_<<"("<<recptetaidisonjetmet_/recptetaidisonjet_<<")";
  LogVerbatim("")<<"Evt with "<<emin<<" Elec(s) (passing pt,eta,Id,Iso,njet,met,trig) = "<<recptetaidisonjetmettrig_<<" [events] AbsEff(RelEff) =  "      <<recptetaidisonjetmettrig_/rec_<<"("<<recptetaidisonjetmettrig_/recptetaidisonjetmet_<<")";

  }



//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



void EwkElecDQM::endRun(const Run& r, const EventSetup&) {
}


//////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



void EwkElecDQM::analyze (const Event & ev, const EventSetup &) {
  
  
  //Trigger
  bool trigger_fired= triggerDecision(ev);

  // Beam spot
  Handle<reco::BeamSpot> beamSpotHandle;
  if (!ev.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle))return;
  
  // Electron collection
  Handle<View<GsfElectron> > electronCollection;
  if (!ev.getByLabel(elecTag_, electronCollection))	return;
  unsigned int electronCollectionSize = electronCollection->size();
  
  // Jet collection
  Handle<View<Jet> > jetCollection;
  if (!ev.getByLabel(jetTag_, jetCollection)) return;
  unsigned int jetCollectionSize = jetCollection->size();

  // MET
  Handle<View<MET> > metCollection;
  if (!ev.getByLabel(metTag_, metCollection)) return;
  const MET& met = metCollection->at(0);
  
  double met_px =met.px();
  double met_py =met.py();
  double met_et = sqrt(met_px*met_px+met_py*met_py);
  
  
  
  //Step1: Loop Over Electrons, and collect the ones passing selection criteria.
  
  double nelec =0,npteta =0, nid=0, niso=0; //counters for electrons at each cut level
 
  std::vector<ewkdqm::Electron> electrons;
  std::vector<ewkdqm::Electron> goodElectrons;
  for (unsigned int i=0; i<electronCollectionSize; i++) {
    
    const GsfElectron& elec = electronCollection->at(i);
    ewkdqm::Electron el(elec);
    electrons.push_back(el);//
    nelec++;
    
    if(!(el.PassPtX(ewp_)&& el.PassEtaX(ewp_)))continue;
    //electrons.push_back(el);//store the crude electrons for n-1 plots    
    npteta++;
    
    if(!(el.PassSIeIeX(ewp_)&& el.PassDetaInX(ewp_)))continue;
    nid++;
    
    if(!(el.PassEisoX(ewp_)&& el.PassHisoX(ewp_) && el.PassTisoX(ewp_)))continue;
    niso++; 
    goodElectrons.push_back(el);//guys passing all the cuts
  }


  
  //Step2: Loop over Jets, and remove the electron candidates from jet list.

  std::vector<TLorentzVector> jets;
  for (unsigned int i=0; i<jetCollectionSize; i++) {
    const Jet& jet = jetCollection->at(i);
    float jet_current_et = jet.et();
    float jet_current_eta = jet.eta();
    bool tooClose2Electron = false;
    
    //loop over electrons to see if any of them is too close 
    //to the current jet.
    for(unsigned int j=0; j<goodElectrons.size();j++){
      double et =  (goodElectrons[j]).pt();
      double dphi = fabs(jet.phi()-(goodElectrons[j]).phi()); 
      if (dphi < 0) dphi = -dphi;
      if (dphi > 3.1415926)dphi = 2 * 3.1415926 - dphi;
      double deta = fabs(jet.eta()-(goodElectrons[j]).eta());
      double delR = std::sqrt(std::pow(deta,2)+std::pow(dphi,2));
      if(delR>0.6)continue;
      if(et>0.0)tooClose2Electron = true;
    }
    
    //jet selection criteria
    if(tooClose2Electron) continue;
    if (jet_current_et < jwp_.ptX_) continue;
    if (fabs(jet_current_eta) > jwp_.etaX_)continue;

    // TLorentzVector jetvec(jet.px(),jet.py(),jet.pz(),jet.energy());
    jets.push_back(TLorentzVector(jet.px(),jet.py(),jet.pz(),jet.energy()));
  }    

  
  //Step3: Sort All the selected collections according to pT
  std::sort(jets.begin(),jets.end(), HighestPt());
  std::reverse(jets.begin(),jets.end());
  
  std::sort(electrons.begin(),electrons.end());
  std::reverse(electrons.begin(),electrons.end());
  
  std::sort(goodElectrons.begin(),goodElectrons.end());
  std::reverse(goodElectrons.begin(),goodElectrons.end());
  
  int nele = electrons.size();
  int neleGood = goodElectrons.size();
  int njets = jets.size();

  ewkdqm::Global g(nele,neleGood, njets,met_et,trigger_fired);
  globalHist_->Fill(g);
  
  for(unsigned int i=0; i<electrons.size(); i++){
    elecHistBefore_->Fill(electrons[i]);
    if(g.PassAllX(gwp_))
      elecHistAfter_->FillAfterCuts((electrons[i]),ewp_);
  }
   
  //Step4: The Boson bussiness goes here.

  ewkdqm::Boson bos(electrons, met_px, met_py); 
  ewkdqm::Boson goodBos(goodElectrons, met_px, met_py); 
    
  ////Some Boson Histograms.
  if(bos.isZ()) zBosonBeforeCuts_->Fill(bos);
  else if(bos.isW())wBosonBeforeCuts_->Fill(bos);
  
  if(g.PassAllX(gwp_)){
    
    if(goodBos.isW())wBosonAfterCuts_->Fill(goodBos);  
    else if(goodBos.isZ()){
      zBosonAfterCuts_->Fill(goodBos);
      if(goodBos.isEE())zBosonEE_->Fill(goodBos);
      else if(goodBos.isEB())zBosonEB_->Fill(goodBos);
      else if(goodBos.isBB())zBosonBB_->Fill(goodBos);    
      
      if(jets.size()){
	jetHistLead_->Fill(jets[0]);
	zBoson1jet_->Fill(goodBos);
      }

      if(jets.size()>1){
	jetHist2ndLead_->Fill(jets[1]);
	zBoson2jet_->Fill(goodBos);
      }
      
      if(jets.size()>2){
	jetHist3rdLead_->Fill(jets[2]);
	zBoson3jet_->Fill(goodBos);
      }
    }
  }
  //Step8: Event Counting.
  int numEleMin = 1;
  if(nele >= numEleMin) rec_++;
  if(npteta >= numEleMin)recpteta_++;
  if(nid >= numEleMin) recptetaid_++;
  if(niso >= numEleMin)recptetaidiso_++;
  //((niso >= numEleMin)&& g.PassNjetsX(gwp_))recptetaidisonjet_++;
  //((niso >= numEleMin) && g.PassNjetsX(gwp_) && g.PassMetX(gwp_))recptetaidisonjetmet_++;
  //((niso >= numEleMin)&& g.PassAllX(gwp_))recptetaidisonjetmettrig_++;
  //if((niso >= numEleMin)&& njetx && metx && trig && zrej) recptetaidisonjetmettrigzrej_++;
  
 return;
  }



// This always returns only a positive deltaPhi
double EwkElecDQM::calcDeltaPhi(double phi1, double phi2) {

  double deltaPhi = phi1 - phi2;

  if (deltaPhi < 0) deltaPhi = -deltaPhi;

  if (deltaPhi > 3.1415926) {
    deltaPhi = 2 * 3.1415926 - deltaPhi;
  }

  return deltaPhi;
}
bool
EwkElecDQM::triggerDecision(const edm::Event& ev){

   edm::Handle<TriggerResults> triggerResults;
   if (!ev.getByLabel(trigTag_, triggerResults)) {
     return 0;
   }
   const std::vector<std::string>& triggerNames = hltConfigProvider_.triggerNames();
   std::vector<std::string> unprescaledPaths;
   
   for(size_t ts=0; ts<triggerNames.size(); ts++){
     
     bool prescaled=false;
     std::string trig = triggerNames[ts];
     const unsigned int prescaleSize=
       hltConfigProvider_.prescaleSize();
     for(unsigned int ps=0; ps<prescaleSize;ps++){
       const unsigned int prescaleValue=hltConfigProvider_.prescaleValue(ps,trig);
       if(prescaleValue !=1 )prescaled = true;
     }
     
     //check that tags_ & occurences_ have same size                                                          
     if(!(tags_.size()==occurences_.size())){
       LogVerbatim("") << ">>> Check Your Cfg";
       exit(0);
     }
     
     bool interestingPath =0;
     for(unsigned int i=0; i!=tags_.size(); i++){
       std::string tag = tags_[i];
       int targetOccurences=occurences_[i];
       int actualOccurences = NumberOfOccurences(tag, trig);
       interestingPath = interestingPath||(targetOccurences==actualOccurences);
     }
     
     
     bool keepThisPath = 1;
     for(unsigned int i=0; i!=vetoes_.size(); i++){
       keepThisPath = keepThisPath&&(NumberOfOccurences(vetoes_[i], trig)==0);
     }
     
     
     if(interestingPath&&(!prescaled)&& keepThisPath)
       unprescaledPaths.push_back(trig);
   } 

   bool trigger_fired = false;
   for(unsigned int i=0; i!=unprescaledPaths.size(); i++){
     std::string hltPath_ = unprescaledPaths[i];
     unsigned int  triggerIndex = hltConfigProvider_.triggerIndex(hltPath_);
     if (triggerIndex < triggerResults->size()) trigger_fired
       = trigger_fired || triggerResults->accept(triggerIndex);
   }
   /*   
   LogVerbatim("") << ">>> Trigger bit: " << trigger_fired << " for one of ( " ;
   for (unsigned int k = 0; k < unprescaledPaths.size(); k++){
     LogVerbatim("") << unprescaledPaths[k] << " ";
   }
   LogVerbatim("") << ")";
   */
   return trigger_fired;
 }
 
