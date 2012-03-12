#ifndef __EWKDQMELECTRON__H_
#define __EWKDQMELECTRON__H_

/*

Anil Singh
Panjab University

*/

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DQM/Physics/interface/dqmGlobal.h"
#include "DQM/Physics/interface/dqmWorkpoint.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"


//
namespace ewkdqm{
  class Electron;
  class ElectronHist;
}

class ewkdqm::Electron : public reco::GsfElectron{

 public:
  Electron (const reco::GsfElectron& e);
  ~Electron ();
  
  bool InBarrel();
  bool InEndcap();

  //basic functions
  bool PassPtX(ewkdqm::eWorkpoint& wp);
  bool PassEtaX(ewkdqm::eWorkpoint& wp);
  bool PassSIeIeX(ewkdqm::eWorkpoint& wp);
  bool PassDetaInX(ewkdqm::eWorkpoint& wp);
  bool PassEisoX(ewkdqm::eWorkpoint& wp);
  bool PassHisoX(ewkdqm::eWorkpoint& wp);
  bool PassTisoX(ewkdqm::eWorkpoint& wp);
 
  bool PassPtEtaX(eWorkpoint& wp);
  bool PassIdX(eWorkpoint& wp);
  bool PassIsoX(eWorkpoint& wp);
  
  
  //- - - - O v e r l o a d e d   O p e r a t o r s - - - - - - - 
  
  bool operator < (const ewkdqm::Electron& e) const; //for sorting
  
 
   
 private:
  Electron();
};
//----i would like this function to return a lorentz vector :)

//- - - - - - - - - - - - - - - - c o n s t r u c t o r - - - - - - - - -
ewkdqm::Electron::Electron(const reco::GsfElectron& e): reco::GsfElectron(e){
  
}

//- - - - - - - - - - - - - - - - d e s t r u c t o r - - - - - - - - - -
ewkdqm::Electron::~Electron(){}






//- - - - - - - - - - - - - - - - m e m b e r  f u n c t i o n s - - - - - - -

bool ewkdqm::Electron::InBarrel(){
  double eta = this->eta();
  return fabs(eta) < 1.446;
}


bool ewkdqm::Electron::InEndcap(){
  double eta =  fabs(this->eta());
  return (eta>=1.566 && eta<2.5);
}


bool ewkdqm::Electron::PassPtX(eWorkpoint& wp){
  bool ptx = this->pt() >= wp.ptX_;
  return ptx;
}


bool ewkdqm::Electron::PassEtaX(ewkdqm::eWorkpoint& wp){
  bool etax = fabs(this->eta())<= wp.etaX_;
  return etax;
}


bool ewkdqm::Electron::PassSIeIeX(ewkdqm::eWorkpoint& wp){
  bool var=0; 
  double sieie = this->scSigmaIEtaIEta();
  if(this->InBarrel()) var = (sieie < wp.sIeIeX_barrel_);
  else if(this->InEndcap())var = (sieie < wp.sIeIeX_endcap_);
  return var;
}


bool ewkdqm::Electron::PassDetaInX(ewkdqm::eWorkpoint& wp){
  
  double deta = 99999;
  if(this->superCluster().isNonnull())
    deta = fabs(this->deltaEtaSuperClusterTrackAtVtx());
  if(this->InBarrel())
    return deta < wp.dEtaInX_barrel_;
  else if(this->InEndcap())
    return deta < wp.dEtaInX_endcap_;
  else return false;
}


bool  ewkdqm::Electron::PassEisoX(eWorkpoint& wp){
  bool var=0;  
  double ecaliso_ =  this->dr03EcalRecHitSumEt();
  if(this->InBarrel())var = (ecaliso_< wp.eIsoX_barrel_);
  else if(this->InEndcap())var = (ecaliso_< wp.eIsoX_endcap_);
  return var;
}


bool ewkdqm::Electron::PassHisoX(eWorkpoint& wp){
  bool var=0;  
  double hcaliso_ =  this->dr03HcalTowerSumEt();
  if(this->InBarrel())var = (hcaliso_< wp.hIsoX_barrel_);
  else if(this->InEndcap())var = (hcaliso_< wp.hIsoX_endcap_);
  return var;
}


bool ewkdqm::Electron::PassTisoX(eWorkpoint& wp){
  bool var=0;  
  double trkiso_  =  this->dr04TkSumPt();  
  if(this->InBarrel())var = (trkiso_< wp.tIsoX_barrel_);
  else if(this->InEndcap())var = (trkiso_< wp.tIsoX_endcap_);
  return var;
}


bool ewkdqm::Electron::operator < (const ewkdqm::Electron& e) const{
  double pt_   =  this->pt();
  return (pt_< e.pt());
}

//composite functions.
bool ewkdqm::Electron::PassIsoX(eWorkpoint& wp){return (PassEisoX(wp)&&PassTisoX(wp)&&PassHisoX(wp));}
bool ewkdqm::Electron::PassIdX(eWorkpoint& wp){return (PassDetaInX(wp) && PassSIeIeX(wp));}
bool ewkdqm::Electron::PassPtEtaX(eWorkpoint& wp){return (PassPtX(wp) && PassEtaX(wp));}



class ewkdqm::ElectronHist{
 public:
  ElectronHist(std::string cutLevel, DQMStore* theDbe);
  ~ElectronHist();
  void Fill(ewkdqm::Electron& e);
  void FillAfterCuts(ewkdqm::Electron& e, ewkdqm::eWorkpoint& ewp);
 private:
  ElectronHist();//dont allow the default constructor.
  MonitorElement* InitHist( std::string qty,  std::string cutLevel, 
			    std::string title,int nbins, int rlow,
			    int rhigh, DQMStore* theDbe);
  
  MonitorElement* pt_;
  MonitorElement* eta_;
  MonitorElement* sieiebarrel_;
  MonitorElement* sieieendcap_;
  MonitorElement* detainbarrel_;
  MonitorElement* detainendcap_;
  MonitorElement* ecalisobarrel_;
  MonitorElement* ecalisoendcap_;
  MonitorElement* hcalisobarrel_;
  MonitorElement* hcalisoendcap_;
  MonitorElement* trkisobarrel_;
  MonitorElement* trkisoendcap_;
};


//default constructor
ewkdqm::ElectronHist::ElectronHist(){}

//overloaded constructor
ewkdqm::ElectronHist::ElectronHist(std::string cutLevel, DQMStore* theDbe){
  
  pt_  =  InitHist( "PT_ELE",cutLevel,"Electron transverse momentum [GeV]",100,0.,100.,theDbe);
  eta_ =  InitHist("ETA_ELE",cutLevel,"Electron pseudo-rapidity",100,-3,3, theDbe);
  sieiebarrel_  =InitHist("SIEIEBARREL",cutLevel,"Electron #sigma_{i#etai#eta} (barrel)", 70,  0., 0.07, theDbe);
  detainbarrel_ =InitHist("DETAINBARREL",cutLevel,"Electron #Delta#eta_{in} (barrel)", 40, -0.02, 0.02, theDbe);
  ecalisobarrel_=InitHist("ECALISOBARREL",cutLevel,"Absolute electron ECAL isolation variable (barrel)" ,50,0.,50.,theDbe);
  hcalisobarrel_=InitHist("HCALISOBARREL",cutLevel,"Absolute electron HCAL isolation variable (barrel)",50,0.,50., theDbe);
  trkisobarrel_ =InitHist("TRKISOBARREL",cutLevel,"Absolute electron track isolation variable (barrel)",50,0.,50., theDbe);
  
//  pt_  =  InitHist( "PT", cutLevel,  "Electron transverse momentum [GeV]", 100, 0., 100., theDbe);
//  eta_ =  InitHist( "ETA",cutLevel,  "Electron pseudo-rapidity", 50,-2.5,2.5, theDbe);
  sieieendcap_  =InitHist("SIEIEENDCAP",  cutLevel,"Electron #sigma_{i#etai#eta} (endcap)", 70,  0.,   0.07, theDbe);
  detainendcap_ =InitHist("DETAINENDCAP", cutLevel,"Electron #Delta#eta_{in} (endcap)",     40, -0.02, 0.02, theDbe);
  ecalisoendcap_=InitHist("ECALISOENDCAP",cutLevel,"Absolute electron ECAL isolation variable (endcap)",50,0.,50.,theDbe);
  hcalisoendcap_=InitHist("HCALISOENDCAP",cutLevel,"Absolute electron HCAL isolation variable (endcap)",50,0.,50., theDbe);
  trkisoendcap_ =InitHist("TRKISOENDCAP", cutLevel,"Absolute electron track isolation variable (endcap)",50,0.,50.,theDbe);
  
}

//destructor
ewkdqm::ElectronHist::~ElectronHist(){
  delete pt_;
  delete eta_;
  delete sieiebarrel_;
  delete sieieendcap_;
  delete detainbarrel_;
  delete detainendcap_;
  delete ecalisobarrel_;
  delete ecalisoendcap_;
  delete hcalisobarrel_;
  delete hcalisoendcap_;
  delete trkisobarrel_;
  delete trkisoendcap_;
}


//Set Histrogram Name, range and title
MonitorElement* 
ewkdqm::ElectronHist::InitHist(std::string qty,  std::string cutLevel,
			       std::string title,int nbins, int rlow,
			       int rhigh, DQMStore* theDbe){
  MonitorElement* m;
  std::string name =qty+"_"+cutLevel;
  m = theDbe->book1D(name.c_str(),title.c_str(),nbins,rlow,rhigh);
  return m;
  
}



//N-1 plots
void 
ewkdqm::ElectronHist::FillAfterCuts(ewkdqm::Electron& e, ewkdqm::eWorkpoint& ewp){
  
  //This function fills N-1 plots after full selection
  //criteria.
  bool ptx=e.PassPtX(ewp);
  bool etax=e.PassEtaX(ewp);
  bool sieiex=e.PassSIeIeX(ewp);
  bool detax=e.PassDetaInX(ewp);
  bool eisox=e.PassEisoX(ewp);
  bool hisox=e.PassHisoX(ewp);
  bool tisox=e.PassTisoX(ewp);
  bool globalx = 1;//relic from old code
  bool fillPt    = etax && sieiex && detax && eisox && hisox && tisox && globalx;  
  bool fillEta   = ptx && sieiex && detax && eisox && hisox && tisox  && globalx;  
  bool fillSieie = etax && ptx && detax && eisox && hisox && tisox    && globalx;  
  bool fillDeta  = etax && ptx && sieiex && eisox && hisox && tisox   && globalx;  
  bool fillEiso  = etax && ptx && sieiex && detax && hisox && tisox   && globalx;  
  bool fillHiso  = etax && ptx && sieiex && detax && eisox && tisox   && globalx;  
  bool fillTiso  = etax && ptx && sieiex && detax && hisox && eisox   && globalx;  


  double deta = 999999;
  if(e.superCluster().isNonnull())
         deta  = e.deltaEtaSuperClusterTrackAtVtx();
  double sieie = e.scSigmaIEtaIEta();
  double eiso  = e.dr03EcalRecHitSumEt();
  double hiso  = e.dr03HcalTowerSumEt();
  double tiso  = e.dr04TkSumPt();


  if(fillPt)
    pt_->Fill(e.pt());
  
  if(fillEta)  
    eta_->Fill(e.eta());
  
  if(fillDeta){
    if(e.InBarrel())detainbarrel_ -> Fill(deta);
    else if(e.InEndcap())detainendcap_ -> Fill(deta);
  }
  
  if(fillSieie){
    if(e.InBarrel())sieiebarrel_ -> Fill(sieie);
    else if(e.InEndcap())sieieendcap_ -> Fill(sieie);
  }
  

  if(fillEiso){
    if(e.InBarrel()) ecalisobarrel_ -> Fill(eiso);
    else if(e.InEndcap())ecalisoendcap_ -> Fill(eiso);
  }
  

  if(fillHiso){
    if(e.InBarrel()) hcalisobarrel_ -> Fill(hiso);
    else if(e.InEndcap())hcalisoendcap_ -> Fill(hiso);
  }
  
  if(fillTiso){
    if(e.InBarrel()) trkisobarrel_ -> Fill(tiso);
    else if(e.InEndcap())trkisoendcap_ -> Fill(tiso);
    }
}

//simple histogram fill.
void ewkdqm::ElectronHist::Fill(ewkdqm::Electron& e){
 
  double deta = 999999;
  if(e.superCluster().isNonnull())
         deta  = e.deltaEtaSuperClusterTrackAtVtx();
  double sieie = e.scSigmaIEtaIEta();
  double eiso  = e.dr03EcalRecHitSumEt();
  double hiso  = e.dr03HcalTowerSumEt();
  double tiso  = e.dr04TkSumPt();
 
  pt_->Fill(e.pt());
  eta_->Fill(e.eta());
  
  if(e.InBarrel()){
    sieiebarrel_  -> Fill(sieie);
    detainbarrel_ -> Fill(deta);
    ecalisobarrel_-> Fill(eiso);
    hcalisobarrel_-> Fill(hiso);
    trkisobarrel_ -> Fill(tiso);
  }
  
  if(e.InEndcap()){
    sieieendcap_  -> Fill(sieie);
    detainendcap_ -> Fill(deta);
    ecalisoendcap_-> Fill(eiso);
    hcalisoendcap_-> Fill(hiso);
    trkisoendcap_ -> Fill(tiso);  
  }
}
#endif
