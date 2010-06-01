#ifndef PANJAB_HISTO_C
#define PANJAB_HISTO_C

#include "TH1D.h"
#include "TFile.h"
#include <string>
#include <map>
#include <sstream>
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
    std::string keyHoE=key+"_HoE";
    std::string keySihih=key+"_Sihih";
    std::string keynBC=key+"_nBC";
    TH1D* hEt  = new TH1D(keyEt.c_str(),keyEt.c_str(),100,10,60);
    TH1D* hEta = new TH1D(keyEta.c_str(),keyEta.c_str(),100,-3,3);
    TH1D* hHoE = new TH1D(keyHoE.c_str(),keyHoE.c_str(),100,0,0.1);
    TH1D* hSihih = new TH1D(keySihih.c_str(),keySihih.c_str(),100,0,0.1);
    TH1D* hnBC = new TH1D(keynBC.c_str(),keynBC.c_str(),10,0,10);
    hisMap[keyEt]=hEt;
    hisMap[keyEta]=hEta;
    hisMap[keyHoE]=hHoE;
    hisMap[keySihih]=hSihih;
    hisMap[keynBC]=hnBC;
  }
}



class scWorkPoint{

public:
  scWorkPoint(std::string name,double ebetx ,double ebhoex,
	      double ebsihsihx ,double ebnbcx)
    :name_(name),ebetx_(ebetx),
     ebhoex_(ebhoex),ebsihsihx_(ebsihsihx),
     ebnbcx_(ebnbcx) {};
  double cutval(std::string);
  std::string name(){return name_;}
private:
  scWorkPoint(){}
  std::string name_;
  double ebetx_;
  double ebhoex_;
  double ebsihsihx_;
  double ebnbcx_;
};

double 
scWorkPoint::cutval(std::string cut)
{
  if(cut=="EBEtX") return ebetx_;
  if(cut=="EBHoverEX") return ebhoex_;
  if(cut=="EBSihSihX") return ebsihsihx_;
  if(cut=="EBnBClstX") return ebnbcx_;
  return -99999;
}



#endif
