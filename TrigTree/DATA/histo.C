#ifndef PANJAB_HISTO_C
#define PANJAB_HISTO_C

#include "TH1D.h"
#include "TFile.h"
#include <string>
#include <map>
#include <sstream>
#include "JetUtilMC.h"
#include <iostream>

namespace util{

std::vector<std::string> cutParser(std::string cutString){
  std::vector<std::string> substrings;
  std::string::const_iterator itr=cutString.begin();
  std::stringstream s1;
  for(; itr!=cutString.end(); itr++){
    if((*itr)!=';')    s1<<*itr;
    else{
      substrings.push_back(s1.str());
      s1.str("");
    }
    
    if(itr==(cutString.end()-1)){
      substrings.push_back(s1.str());
      s1.str("");
    }
  } 
  return substrings;
}

}

class histo{
public:
  histo(std::string trigNam, std::vector<std::string> selVec,
	TFile* rootfile,std::string name);
  std::map<std::string, TH1D*> hisMap;
private:
  TFile* file;
  
  histo(){}
};
histo::histo(std::string trigNam, std::vector<std::string> selecVec, 
	     TFile* rootfile, std::string name){
  // file=rootfile;
  // file->cd();
  //  TDirectory* dir=gDirectory->mkdir(trigNam.c_str());
  //  TDirectory* dir=gDirectory->mkdir((name+trigNam).c_str());
  // dir->cd();
  for(unsigned int i=0; i!=selecVec.size(); i++){
    std::vector<std::string> cutNams=util::cutParser(selecVec[i]);
    std::string key=trigNam+"_"+ cutNams[cutNams.size()-1];
    std::string keyEt=key+"_Et";
    std::string keyEta=key+"_Eta";
    std::string keyPhi=key+"_Phi";
    std::string keyHoE=key+"_HoE";
    std::string keySihih=key+"_Sihih";
    std::string keydEta=key+"_dEta";
    std::string keydPhi=key+"_dPhi";
    std::string keyTrkIso=key+"_TrkIso";
    std::string keyEcalIso=key+"_EcalIso";
    std::string keyHcalIso=key+"_HcalIso";
    std::string keySpikes=key+"_Spikes";
    std::string keyMisHits=key+"_MisHits";
    std::string keyDist=key+"_Dist";
    std::string keydCotTheta=key+"_dCotTheta";
    std::string keyHltMatch=key+"_HltMatch";
    std::string keyEcalDriven=key+"_EcalDriven";
    std::string keyTrkDriven=key+"_TrkDriven";
    TH1D* hEt  = new TH1D(keyEt.c_str(),keyEt.c_str(),100,0,60);
    TH1D* hEta = new TH1D(keyEta.c_str(),keyEta.c_str(),100,-2.5,2.5);
    TH1D* hPhi = new TH1D(keyPhi.c_str(),keyPhi.c_str(),100,-3.5,3.5);
    TH1D* hHoE = new TH1D(keyHoE.c_str(),keyHoE.c_str(),100,0,0.1);
    TH1D* hSihih = new TH1D(keySihih.c_str(),keySihih.c_str(),100,0,0.1);
    TH1D* hdEta = new TH1D(keydEta.c_str(),keydEta.c_str(),100,-0.05,0.05);
    TH1D* hdPhi = new TH1D(keydPhi.c_str(),keydPhi.c_str(),100,-1.2,1.2);
    TH1D* hTrkIso = new TH1D(keyTrkIso.c_str(),keyTrkIso.c_str(),100,0.,5.);
    TH1D* hEcalIso = new TH1D(keyEcalIso.c_str(),keyEcalIso.c_str(),100,0.,5.);
    TH1D* hHcalIso = new TH1D(keyHcalIso.c_str(),keyHcalIso.c_str(),100,0.,5.);
    TH1D* hSpikes = new TH1D(keySpikes.c_str(),keySpikes.c_str(),100,0.,1.);
    TH1D* hMisHits = new TH1D(keyMisHits.c_str(),keyMisHits.c_str(),10,-1.,9.);
    TH1D* hDist = new TH1D(keyDist.c_str(),keyDist.c_str(),100,-0.1,0.1);
    TH1D* hdCotTheta = new TH1D(keydCotTheta.c_str(),keydCotTheta.c_str(),100,-0.1,0.1);
    TH1D* hHltMatch = new TH1D(keyHltMatch.c_str(),keyHltMatch.c_str(),100,0.,4.);
    TH1D* hEcalDriven = new TH1D(keyEcalDriven.c_str(),keyEcalDriven.c_str(),2,0,2);
    TH1D* hTrkDriven = new TH1D(keyTrkDriven.c_str(),keyTrkDriven.c_str(),2,0,2);
    hisMap[keyEt]=hEt;
    hisMap[keyEta]=hEta;
    hisMap[keyPhi]=hPhi;
    hisMap[keyHoE]=hHoE;
    hisMap[keySihih]=hSihih;
    hisMap[keydEta]=hdEta;
    hisMap[keydPhi]=hdPhi;
    hisMap[keyTrkIso]=hTrkIso;
    hisMap[keyEcalIso]=hEcalIso;
    hisMap[keyHcalIso]=hHcalIso;
    hisMap[keySpikes]=hSpikes;
    hisMap[keyMisHits]=hMisHits;
    hisMap[keyDist]=hDist;
    hisMap[keydCotTheta]=hdCotTheta;
    hisMap[keyHltMatch]=hHltMatch;
    hisMap[keyEcalDriven]=hEcalDriven;
    hisMap[keyTrkDriven]=hTrkDriven;
  }
}



class eleWorkPoint{

public:
  eleWorkPoint(std::string name,double spikex, double etx ,
	       unsigned int mishitsx, double distx, double dcotthetax,
	       double trkisox,double ecalisox, double hcalisox,
	       double sihsihx ,double dphix,double detax,double hoex,
	       double hltmatchx)
    :name_(name),spikex_(spikex),etx_(etx),
     mishitsx_(mishitsx),distx_(distx),dcotthetax_(dcotthetax),
     trkisox_(trkisox),ecalisox_(ecalisox),hcalisox_(hcalisox),
     sihsihx_(sihsihx),dphix_(dphix),detax_(detax),hoex_(hoex),
     hltmatchx_(hltmatchx){};
  double cutval(std::string);
  std::string name(){return name_;}
private:
  eleWorkPoint(){}
  std::string name_;  double spikex_;  double etx_;
  unsigned int mishitsx_; double distx_; double dcotthetax_;
  double trkisox_;  double ecalisox_;  double hcalisox_;
  double sihsihx_;  double dphix_;  double detax_; double hoex_;
  double hltmatchx_;
};

double 
eleWorkPoint::cutval(std::string cut)
{
  if(cut=="SpikeRemovalX") return spikex_;
  if(cut=="EtX") return etx_;
  if(cut=="MisHitsX") return mishitsx_;
  if(cut=="DistX") return distx_;
  if(cut=="dCotThetaX") return dcotthetax_;
  if(cut=="TrkIsoX") return trkisox_;
  if(cut=="EcalIsoX") return ecalisox_;
  if(cut=="HcalIsoX") return hcalisox_;
  if(cut=="SihSihX") return sihsihx_;
  if(cut=="dPhiX") return dphix_;
  if(cut=="dEtaX") return detax_;
  if(cut=="HoverEX") return hoex_;
  if(cut=="HltMatchX") return hltmatchx_;
  return -99999;
}



#endif
