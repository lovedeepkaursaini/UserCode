//Get the electron/SC variables after each seq. cut of WP80, 95 etc.
//WP's: a cut sequence to be applied to get results for VBTF_ICHEP
//To Do: Put a README to explain the algorithm of this code
//Here is very simple explanation:
// (1)Collect all the electrons into two separate vectors (endcap and barrel)
// (2)Test them for satisfying the cuts.
// (3)Now the leading electrons from both the catagories are used to fill d'butions.      
//since VBTF has defined a new HLT_Trigger, so the implemenation of that is herein as "HLT_Photon10now15_L1R", suggested by Bryan
// The requirement is to get the efficiency curves (turn-on curves), which are obtained after dividing the HLT* histo with NoTrig* histo as a function of ScEt.
// Coversion rejection is also here in
// spike cleaning done
// 
//
// Lovedeep Kaur Saini , Panjab Univesity, Chandigarh.

 
#define EleTrigEffi_cxx
#include "EleTrigEffi.h"
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

void fillSelHisto(eleColl& Collection, histo* h){
  
  eleColl::iterator iter=Collection.begin();
  for(;iter!=Collection.end();iter++){
    std::string name =iter->first;
    std::vector<EleKin> coll=iter->second;
    if(coll.size())(h->hisMap[name+"_Et"])->Fill(coll[0].Et());
    if(coll.size())(h->hisMap[name+"_Eta"])->Fill(coll[0].Eta());
    if(coll.size())(h->hisMap[name+"_Phi"])->Fill(coll[0].Phi());
    if(coll.size())(h->hisMap[name+"_HoE"])->Fill(coll[0].HoverE());
    if(coll.size())(h->hisMap[name+"_Sihih"])->Fill(coll[0].sigmaIetaIeta());
    if(coll.size())(h->hisMap[name+"_dEta"])->Fill(coll[0].dEta());
    if(coll.size())(h->hisMap[name+"_dPhi"])->Fill(coll[0].dPhi());
    if(coll.size())(h->hisMap[name+"_TrkIso"])->Fill(coll[0].TrkIso());
    if(coll.size())(h->hisMap[name+"_EcalIso"])->Fill(coll[0].EcalIso());
    if(coll.size())(h->hisMap[name+"_HcalIso"])->Fill(coll[0].HcalIso());
    if(coll.size())(h->hisMap[name+"_Spikes"])->Fill(coll[0].Spikes());
    if(coll.size())(h->hisMap[name+"_MisHits"])->Fill(coll[0].MisHits());
    if(coll.size())(h->hisMap[name+"_Dist"])->Fill(coll[0].Dist());
    if(coll.size())(h->hisMap[name+"_dCotTheta"])->Fill(coll[0].dCotTheta());
    if(coll.size())(h->hisMap[name+"_HltMatch"])->Fill(coll[0].HltMatch());
    if(coll.size())(h->hisMap[name+"_EcalDriven"])->Fill(coll[0].EcalDriven());
    if(coll.size())(h->hisMap[name+"_TrkDriven"])->Fill(coll[0].TrkDriven());
  }
}


void EleTrigEffi::Loop(Long64_t min, Long64_t max,std::string input)
{

   std::stringstream ss;
   ss<<"Ele_EffiHisto_"<<input<<"_"<<min<<"_"<<max<<"_out.root";
   std::string filenam = ss.str();	
   TFile * file = new TFile(filenam.c_str(),"recreate");
   std::vector<std::string> selVec; 
   selVec.push_back("NoCut");
   selVec.push_back("NoCut;SpikeRemovalX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX;EcalIsoX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX;EcalIsoX;HcalIsoX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX;EcalIsoX;HcalIsoX;SihSihX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX;EcalIsoX;HcalIsoX;SihSihX;dPhiX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX;EcalIsoX;HcalIsoX;SihSihX;dPhiX;dEtaX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX;EcalIsoX;HcalIsoX;SihSihX;dPhiX;dEtaX;HoverEX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX;EcalIsoX;HcalIsoX;SihSihX;dPhiX;dEtaX;HoverEX;HltMatchX");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX;EcalIsoX;HcalIsoX;SihSihX;dPhiX;dEtaX;HoverEX;HltMatchX;EcalDriven");
   selVec.push_back("NoCut;SpikeRemovalX;EtX;MisHitsX;DistX;dCotThetaX;TrkIsoX;EcalIsoX;HcalIsoX;SihSihX;dPhiX;dEtaX;HoverEX;HltMatchX;TrkDriven");

   eleWorkPoint* eleWP80EE=new eleWorkPoint("eleWP80EE",0.9,10.,0,0.02,0.02,0.04,0.05,0.025,0.03,0.03,0.007,0.025,0.5);
   eleWorkPoint* eleWP95EE=new eleWorkPoint("eleWP95EE",0.9,10.,1,0.,0.,0.08,0.06,0.05,0.03,0.7,0.01,0.07,0.5);
 
   eleWorkPoint* eleWP80EB=new eleWorkPoint("eleWP80EB",0.9,10.,0,0.02,0.02,0.09,0.07,0.10,0.01,0.06,0.004,0.04,0.5);
   eleWorkPoint* eleWP95EB=new eleWorkPoint("eleWP95EB",0.9,10.,1,0.,0.,0.15,2.,0.12,0.01,0.8,0.007,0.15,0.5);
   
   UseBarrelEl();
   Loop(selVec,eleWP80EB,file,min,max);
   //Loop(selVec,eleWP95EB,file,min,max);
   UseEndCapEl();
   Loop(selVec,eleWP80EE,file,min,max);
   //Loop(selVec,eleWP95EE,file,min,max);
   
   file->cd();
   file->Write();
   file->Close();
   
}

void 
EleTrigEffi::Loop(std::vector<std::string> selVec, eleWorkPoint* eleWP, TFile* file, Long64_t min, Long64_t max){
  file->cd();
  std::string name="";
  std::string dirNam=eleWP->name();
  TDirectory* dir=gDirectory->mkdir(dirNam.c_str());
  dir->cd();

  histo* test=new histo("NoTrig",selVec,file,name);

  histo* testHltPhot15=new histo("HLT_Photon15_L1R",selVec,file,name);
  histo* testHltPhot10=new histo("HLT_Photon10_L1R",selVec,file,name);
  histo* testHltPhot10now15=new histo("HLT_Photon10now15_L1R",selVec,file,name);

  histo* testL1EG8=new histo("L1_SingleEG8",selVec,file,name);
  histo* testL1EG5=new histo("L1_SingleEG5",selVec,file,name);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   if(min>nentries) return;
   if(max>nentries) max=nentries;


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=min; jentry<max;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //if(!(EvtInfo_BunchXing==1 || EvtInfo_BunchXing==1786))continue;
      
      if(jentry%10000==0)  
	std::cout<<"Jentry: "<<jentry<<std::endl;
      
      trigStuff.clear();
      
      for(unsigned int i=0; i!=trigName->size(); i++)
	trigStuff[trigName->at(i)]=trigResults->at(i);
      
      for(unsigned int i=0; i!=L1trigName->size(); i++)
	trigStuff[L1trigName->at(i)]=L1trigResults->at(i);
      
      trigStuff["NoTrig"]=1;
      trigStuff["HLT_Photon10now15_L1R"]=trigStuff["HLT_Photon10_L1R"];
      
      eleColl noTrigColl=MkCollection("NoTrig",selVec,eleWP);

      eleColl HltPhot10now15TrigColl=MkCollection("HLT_Photon10now15_L1R",selVec,eleWP);
      eleColl HltPhot15TrigColl=MkCollection("HLT_Photon15_L1R",selVec,eleWP);
      eleColl HltPhot10TrigColl=MkCollection("HLT_Photon10_L1R",selVec,eleWP);
      eleColl L1EG8TrigColl=MkCollection("L1_SingleEG8",selVec,eleWP);
      eleColl L1EG5TrigColl=MkCollection("L1_SingleEG5",selVec,eleWP);
      
      fillSelHisto(noTrigColl,test);

      fillSelHisto(HltPhot10now15TrigColl,testHltPhot10now15);
      fillSelHisto(HltPhot15TrigColl,testHltPhot15);
      fillSelHisto(HltPhot10TrigColl,testHltPhot10);
      fillSelHisto(L1EG8TrigColl,testL1EG8);
      fillSelHisto(L1EG5TrigColl,testL1EG5);
   }
}

eleColl 
EleTrigEffi::MkCollection(std::string trigName, std::vector<std::string> selecVec,eleWorkPoint* eleWP){
  eleColl collection;
  for(unsigned int i=0; i!=selecVec.size(); i++){
    std::vector<std::string> cutNams=cutParser(selecVec[i]);
    std::string key=trigName+"_"+ cutNams[cutNams.size()-1];
    std::vector<EleKin> ELEcoll = SelectElectrons(trigName,selecVec[i],eleWP);
    collection[key]=ELEcoll;
  }
  return collection;
}

std::vector<EleKin>
EleTrigEffi::SelectElectrons(std::string trigName, std::string cutString, eleWorkPoint* eleWP){
  std::vector<EleKin> EleVector; 

  for(int i=0; i !=ElectronNum; i++){
    double Et=ElectronscEt->at(i);//changed after meeting : ElectronEt->at(i);
    double Eta=ElectronEta->at(i);
    double Phi=ElectronPhi->at(i);
    double HoverE=ElectronelehcalOverEcal->at(i);
    double sigmaIetaIeta=ElectronelesigmaIetaIeta->at(i);
    double DhSuperClsTrkAtVtx = ElectronDhSuperClsTrkAtVtx->at(i);
    double DphiSuperClsTrkAtVtx=ElectronDphiSuperClsTrkAtVtx->at(i);
    double dr03TkSumPt=Electroneledr03TkSumPt->at(i)/ElectronPt->at(i);
    double dr03EcalRecHitSumEt=Electroneledr03EcalRecHitSumEt->at(i)/ElectronPt->at(i);
    double dr03HcalTowerSumEt=Electroneledr03HcalTowerSumEt->at(i)/ElectronPt->at(i);
    double rmax3x3=Electronelermax3x3->at(i);    
    double MisHits=ElectroneleMissingHits->at(i);
    double Dist=ElectroneleDist->at(i);
    double dCotTheta=ElectroneleDeltaCotTheta->at(i);
    double RelIsoEcal= ElectroneleRelIsoEcal->at(i);
    double RelIsoHcal=ElectroneleRelIsoHcal->at(i);
    double RelIsoTrk=ElectroneleRelIsoTrk->at(i);
    bool ecalDriven= ElectroneleisEcalDriven->at(i);
    bool trkDriven= ElectroneleisTrackerDriven->at(i);

    if(useEndCapEl_ && !(ElectroneleisEndcap->at(i)))continue;
    if(useBrlEl_ && !(ElectroneleisBarrel->at(i)))continue;

    double dRHltMatch=FindHltMatch(trigName, Eta, Phi);

    std::map<std::string, bool> cutRcrd;
    cutRcrd["NoCut"]=1;
    cutRcrd["SpikeRemovalX"]=(rmax3x3 < eleWP->cutval("SpikeRemovalX"));
    cutRcrd["EtX"]=(Et > eleWP->cutval("EtX"));
    cutRcrd["MisHitsX"]=(MisHits <= eleWP->cutval("MisHitsX"));
    cutRcrd["DistX"]=(fabs(Dist) >= eleWP->cutval("DistX"));
    cutRcrd["dCotThetaX"]=(fabs(dCotTheta) >= eleWP->cutval("dCotThetaX"));
    cutRcrd["TrkIsoX"]=(RelIsoTrk)<eleWP->cutval("TrkIsoX");
    cutRcrd["EcalIsoX"]=(RelIsoEcal)<eleWP->cutval("EcalIsoX");
    cutRcrd["HcalIsoX"]=(RelIsoHcal)<eleWP->cutval("HcalIsoX");
    cutRcrd["SihSihX"]=(sigmaIetaIeta < eleWP->cutval("SihSihX"));
    cutRcrd["dPhiX"]=(fabs(DphiSuperClsTrkAtVtx)<eleWP->cutval("dPhiX"));
    cutRcrd["dEtaX"]=(fabs(DhSuperClsTrkAtVtx)<eleWP->cutval("dEtaX"));
    cutRcrd["HoverEX"]=(HoverE < eleWP->cutval("HoverEX"));
    cutRcrd["HltMatchX"]=(dRHltMatch < eleWP->cutval("HltMatchX"));
    cutRcrd["EcalDriven"]=ecalDriven;
    cutRcrd["TrkDriven"]=trkDriven;

  
    std::vector<std::string> cutsApplied = cutParser(cutString);
    bool pass= 1;
    for(unsigned int i=0; i!=cutsApplied.size();i++){
      pass=pass&&cutRcrd[cutsApplied[i]];
    }

    EleKin* ele=new EleKin(Et,Eta,Phi,HoverE,sigmaIetaIeta,DhSuperClsTrkAtVtx, DphiSuperClsTrkAtVtx, RelIsoTrk, RelIsoEcal, RelIsoHcal,rmax3x3,MisHits,Dist,dCotTheta,dRHltMatch,ecalDriven,trkDriven);
    if(trigStuff[trigName] && pass)
      EleVector.push_back(*ele);
  }
  return EleVector;
}


double EleTrigEffi::FindHltMatch(std::string trigname, double Eta, double Phi)
{

  //there was a problem of filling tree, so i am trying ....

  double min=9999.;
  if ( trigname == "HLT_Photon15_L1R" )
  {
    for ( unsigned int i = 0; i < hltPhoton15Pt->size(); ++i )
      {
	double dr=99999.;
	dr = radius(hltPhoton15Eta->at(i), hltPhoton15Phi->at(i), Eta, Phi);
	if ( dr < min )
	  {
	    min= dr;
	  }
      }
  }

  else if ( trigname == "HLT_Photon10_L1R" )
  {
    for ( unsigned int i = 0; i < hltPhoton10Pt->size(); ++i )
      {
	double dr=99999.;
	dr = radius(hltPhoton10Eta->at(i), hltPhoton10Phi->at(i), Eta, Phi);
	if ( dr < min )
	  {
	    min= dr;
	  }
      }
  }
  else if ( trigname == "HLT_Photon10now15_L1R" )
  {
    for ( unsigned int i = 0; i < hltPhoton10now15Pt->size(); ++i )
      {
	double dr=99999.;
	dr = radius(hltPhoton10now15Eta->at(i), hltPhoton10now15Phi->at(i), Eta, Phi);
	if ( dr < min )
	  {
	    min= dr;
	  }
      }
  }
  else if(trigname == "NoTrig") min=0.;

  return min;
}


std::vector<std::string> 
EleTrigEffi::cutParser(std::string cutString){
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
