#include "PhysicsTools/TagAndProbe/interface/ZGeneratorLineShape.h"

ClassImp(ZGeneratorLineShape)

ZGeneratorLineShape::ZGeneratorLineShape(
					 const char *name, const char *title, 
					 RooAbsReal& _m,
  			                 RooAbsReal& _ptL,
					 RooAbsReal& _ptH,
                                         RooAbsReal& _etaL,
                                         RooAbsReal& _etaH,
                                         RooAbsReal& _pfail,
					 //string ptBin,string etaBin,string pf,
					 char* genfile, char* histoName
					 ): 
  RooAbsPdf(name,title),
  m("m","m", this,_m),  
  dataHist(0)
{
  //std::cout<<"TOKEN-------------------------------------------------------------------------------->:  "<<ptBin<<std::endl;
  TFile *f_gen= TFile::Open(genfile);

  double ptLow = _ptL.getVal();
  double ptHi = _ptH.getVal();
  double etaLw = _etaL.getVal();
  double etaHi = _etaH.getVal();
  double pf = _pfail.getVal();
  
  std::stringstream ss;
  std::string pass = "Pass";
  if(pf==0)pass = "Fail";
  ss<<"hMass_"<<ptLow<<"To"<<ptHi<<"_"<<etaLw<<"To"<<etaHi<<"_"<<pass;
  std::string histName (ss.str());
  histoName = (char*)histName.c_str();
  std::cout<<"hName: "<<histoName<<std::endl;
  TH1F* mass_th1f = (TH1F*)  f_gen->Get(histoName);
  
  dataHist = new RooDataHist("Mass_gen", "Mass_gen", _m, mass_th1f );
  f_gen->Close();
}


ZGeneratorLineShape::ZGeneratorLineShape(const ZGeneratorLineShape& other, const char* name):
  RooAbsPdf(other,name),
  m("m", this,other.m),
  dataHist(other.dataHist)
{
}


Double_t ZGeneratorLineShape::evaluate() const{

  // std::cout<<"gen shape: m, evaluate= "<<m<<", "<<dataHist->weight(m.arg())<<std::endl;
  return dataHist->weight(m.arg()) ;
}
