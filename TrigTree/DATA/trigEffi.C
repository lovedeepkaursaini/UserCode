#define trigEffi_cxx
#include "trigEffi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2.h> 
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <vector>
#include "TDirectory.h"
using namespace std;

void fillSelHisto(scColl& Collection, histo* h){
  
  scColl::iterator iter=Collection.begin();
  for(;iter!=Collection.end();iter++){
    std::string name =iter->first;
    std::vector<xyz> coll=iter->second;
    if(coll.size())(h->hisMap[name+"_Et"])->Fill(coll[0].Pt());
    if(coll.size())(h->hisMap[name+"_Eta"])->Fill(coll[0].Eta());
    if(coll.size())(h->hisMap[name+"_HoE"])->Fill(coll[0].HoverE());
    if(coll.size())(h->hisMap[name+"_Sihih"])->Fill(coll[0].sigmaIetaIeta());
    if(coll.size())(h->hisMap[name+"_nBC"])->Fill(coll[0].nBasicClusters());
  }
}



void 
trigEffi::Loop()
{
 
   TFile * file = new TFile("EffiHisto_data.root","recreate");
        
   std::vector<std::string> selVec; 
   selVec.push_back("NoCut");
   selVec.push_back("NoCut;EBEtX");
   selVec.push_back("NoCut;EBEtX;EBHoverEX");
   selVec.push_back("NoCut;EBEtX;EBHoverEX;EBSihSihX");
   selVec.push_back("NoCut;EBEtX;EBHoverEX;EBSihSihX;EBnBClstX");

   scWorkPoint* scWP95=new scWorkPoint("scWP95",10,0.05,0.01,1);
//   scWorkPoint* scWP2=new scWorkPoint("scWP2",15,0.03,0.0155,1);
//   scWorkPoint* scWP3=new scWorkPoint("scWP3",20,0.02,0.0055,1);

   Loop(selVec,scWP95,file);
//   Loop(selVec,scWP2,file);
//   Loop(selVec,scWP3,file);

   file->cd();
   file->Write();
   file->Close();
   
}

void 
trigEffi::Loop(std::vector<std::string> selVec, scWorkPoint* scWP, TFile* file){
  file->cd();
  std::string name="";
  std::string dirNam=scWP->name();
  TDirectory* dir=gDirectory->mkdir(dirNam.c_str());
  dir->cd();

  histo* test=new histo("NoTrig",selVec,file,name);
  histo* testL1_2=new histo("L1_SingleEG2",selVec,file,name);
  histo* testL1_5=new histo("L1_SingleEG5",selVec,file,name);
  histo* testL1_8=new histo("L1_SingleEG8",selVec,file,name);
  histo* testL1_10=new histo("L1_SingleEG10",selVec,file,name);
  histo* testL1_12=new histo("L1_SingleEG12",selVec,file,name);
  histo* testL1_15=new histo("L1_SingleEG15",selVec,file,name);
  histo* testL1_20=new histo("L1_SingleEG20",selVec,file,name);
  histo* testL1Iso_5=new histo("L1_SingleIsoEG5",selVec,file,name);
  histo* testL1Iso_8=new histo("L1_SingleIsoEG8",selVec,file,name);
  histo* testL1Iso_10=new histo("L1_SingleIsoEG10",selVec,file,name);
  histo* testL1Iso_12=new histo("L1_SingleIsoEG12",selVec,file,name);
  histo* testL1Iso_15=new histo("L1_SingleIsoEG15",selVec,file,name);
  histo* testHltEle10=new histo("HLT_Ele10_LW_L1R",selVec,file,name);
  histo* testHltEle15=new histo("HLT_Ele15_LW_L1R",selVec,file,name);
  histo* testHltEle20=new histo("HLT_Ele20_LW_L1R",selVec,file,name);

  histo* testHltEle10EleId=new histo("HLT_Ele10_LW_EleId_L1R",selVec,file,name);
  histo* testHltEle15SiStrip=new histo("HLT_Ele15_SiStrip_L1R",selVec,file,name);

  histo* testHltPhot10=new histo("HLT_Photon10_L1R",selVec,file,name);
  histo* testHltPhot15=new histo("HLT_Photon15_L1R",selVec,file,name);
  histo* testHltPhot20=new histo("HLT_Photon20_L1R",selVec,file,name);
  histo* testHltPhot30=new histo("HLT_Photon30_L1R",selVec,file,name);
  
  histo* testHltPho15TrackIso=new histo("HLT_Photon15_TrackIso_L1R",selVec,file,name);
  histo* testHltPho15LooseEcalIso=new histo("HLT_Photon15_LooseEcalIso_L1R",selVec,file,name);

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=1; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if(!(EvtInfo_BunchXing==1 || EvtInfo_BunchXing==1786))continue;
    
    if(jentry%10000==0)  std::cout<<"Jentry: "<<jentry<<std::endl;
    
    trigStuff.clear();
    
    for(unsigned int i=0; i!=trigName->size(); i++)
      trigStuff[trigName->at(i)]=trigResults->at(i);
    
    for(unsigned int i=0; i!=L1trigName->size(); i++)
      trigStuff[L1trigName->at(i)]=L1trigResults->at(i);
    
    trigStuff["NoTrig"]=1;
    
    scColl noTrigColl=MkCollection("NoTrig",selVec,scWP);
    scColl L1EG2TrigColl=MkCollection("L1_SingleEG2",selVec,scWP);
    scColl L1EG5TrigColl=MkCollection("L1_SingleEG5",selVec,scWP);
    scColl L1EG8TrigColl=MkCollection("L1_SingleEG8",selVec,scWP);
    scColl L1EG10TrigColl=MkCollection("L1_SingleEG10",selVec,scWP);
    scColl L1EG12TrigColl=MkCollection("L1_SingleEG12",selVec,scWP);
    scColl L1EG15TrigColl=MkCollection("L1_SingleEG15",selVec,scWP);
    scColl L1EG20TrigColl=MkCollection("L1_SingleEG20",selVec,scWP);
    scColl L1EGIso5TrigColl=MkCollection("L1_SingleIsoEG5",selVec,scWP);
    scColl L1EGIso8TrigColl=MkCollection("L1_SingleIsoEG8",selVec,scWP);
    scColl L1EGIso10TrigColl=MkCollection("L1_SingleIsoEG10",selVec,scWP);
    scColl L1EGIso12TrigColl=MkCollection("L1_SingleIsoEG12",selVec,scWP);
    scColl L1EGIso15TrigColl=MkCollection("L1_SingleIsoEG15",selVec,scWP);
    
    scColl HltEle10TrigColl=MkCollection("HLT_Ele10_LW_L1R",selVec,scWP);
    scColl HltEle15TrigColl=MkCollection("HLT_Ele15_LW_L1R",selVec,scWP);
    scColl HltEle20TrigColl=MkCollection("HLT_Ele20_LW_L1R",selVec,scWP);
    
    scColl HltEle10EleIdTrigColl=MkCollection("HLT_Ele10_LW_EleId_L1R",selVec,scWP);
    scColl HltEle15SiStripTrigColl=MkCollection("HLT_Ele15_SiStrip_L1R",selVec,scWP);
    
    scColl HltPhot10TrigColl=MkCollection("HLT_Photon10_L1R",selVec,scWP);
    scColl HltPhot15TrigColl=MkCollection("HLT_Photon15_L1R",selVec,scWP);
    scColl HltPhot20TrigColl=MkCollection("HLT_Photon20_L1R",selVec,scWP);
    scColl HltPhot30TrigColl=MkCollection("HLT_Photon30_L1R",selVec,scWP);
    
    scColl HLTPhot15TrackIsoColl=MkCollection("HLT_Photon15_TrackIso_L1R",selVec,scWP);
    scColl HLTPhot15LooseEcalIsoColl=MkCollection("HLT_Photon15_LooseEcalIso_L1R",selVec,scWP);
    
    fillSelHisto(noTrigColl,test);
    fillSelHisto(L1EG2TrigColl,testL1_2);
    fillSelHisto(L1EG5TrigColl,testL1_5);
    fillSelHisto(L1EG8TrigColl,testL1_8);
    fillSelHisto(L1EG10TrigColl,testL1_10);
    fillSelHisto(L1EG12TrigColl,testL1_12);
    fillSelHisto(L1EG15TrigColl,testL1_15);
    fillSelHisto(L1EG20TrigColl,testL1_20);
    fillSelHisto(L1EGIso5TrigColl,testL1Iso_5);
    fillSelHisto(L1EGIso8TrigColl,testL1Iso_8);
    fillSelHisto(L1EGIso10TrigColl,testL1Iso_10);
    fillSelHisto(L1EGIso12TrigColl,testL1Iso_12);
    fillSelHisto(L1EGIso15TrigColl,testL1Iso_15);

    fillSelHisto(HltEle10TrigColl,testHltEle10);
    fillSelHisto(HltEle15TrigColl,testHltEle15);
    fillSelHisto(HltEle20TrigColl,testHltEle20);
    
    fillSelHisto(HltEle10EleIdTrigColl,testHltEle10EleId);
    fillSelHisto(HltEle15SiStripTrigColl,testHltEle15SiStrip);
    
    fillSelHisto(HltPhot10TrigColl,testHltPhot10);
    fillSelHisto(HltPhot15TrigColl,testHltPhot15);
    fillSelHisto(HltPhot20TrigColl,testHltPhot20);
    fillSelHisto(HltPhot30TrigColl,testHltPhot30);
    
    fillSelHisto(HLTPhot15TrackIsoColl,testHltPho15TrackIso);
    fillSelHisto(HLTPhot15LooseEcalIsoColl,testHltPho15LooseEcalIso);
  }
}

scColl 
trigEffi::MkCollection(std::string trigName, std::vector<std::string> selecVec,scWorkPoint* scWP){
  scColl collection;
  for(unsigned int i=0; i!=selecVec.size(); i++){
    std::vector<std::string> cutNams=cutParser(selecVec[i]);
    std::string key=trigName+"_"+ cutNams[cutNams.size()-1];
    std::vector<xyz> SCcoll = SelectElectrons(trigName,selecVec[i],scWP);
    collection[key]=SCcoll;
  }
  return collection;
}

std::vector<xyz>
trigEffi::SelectElectrons(std::string trigName, std::string cutString, scWorkPoint* scWP){
  std::vector<xyz> ScVector; 

  for(int i=0; i !=scEBNum; i++){
    double EBEt=scEBEt->at(i);
    double EBEta=scEBEta->at(i);
    double EBPhi=scEBPhi->at(i);
    double EBEnergy1=scEBenergy1->at(i);
    double EBEnergy3x3=scEBenergy3x3->at(i);
    double EBHoverE=scEBHoverE->at(i);
    double EBsigmaIetaIeta=scEBsigmaIetaIeta->at(i);
    double EBnBasicClusters=scEBnBasicClusters->at(i);
    
    if((EBEnergy1/EBEnergy3x3)>0.9)continue;
    std::map<std::string, bool> cutRcrd;
    cutRcrd["EBEtX"]=(EBEt > scWP->cutval("EBEtX"));
    cutRcrd["EBHoverEX"]=(EBHoverE < scWP->cutval("EBHoverEX"));
    cutRcrd["EBSihSihX"]=(EBsigmaIetaIeta < scWP->cutval("EBSihSihX"));
    cutRcrd["EBnBClstX"]=(EBnBasicClusters == scWP->cutval("EBnBClstX"));
    cutRcrd["NoCut"]=1;

    
    std::vector<std::string> cutsApplied = cutParser(cutString);
    bool pass= 1;
    for(int i=0; i!=cutsApplied.size();i++){
      pass=pass&&cutRcrd[cutsApplied[i]];
    }

    xyz* sc=new xyz(EBEt,EBEta,EBPhi,EBHoverE,EBsigmaIetaIeta,EBnBasicClusters);    
    if(trigStuff[trigName] && pass)
      ScVector.push_back(*sc);
  }
  return ScVector;
}


std::vector<std::string> 
trigEffi::cutParser(std::string cutString){
  std::vector<string> substrings;
  std::string::const_iterator itr=cutString.begin();
  stringstream s1;
  for(; itr!=cutString.end(); itr++){
    if((*itr)!=';')
      s1<<*itr;
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
