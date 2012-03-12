#include <TLorentzVector.h>
#include "DQM/Physics/interface/dqmWorkpoint.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"


/*
  Anil Singh
  Panjab University
*/


namespace ewkdqm{
  class Boson;
  class BosonHist;
}


class ewkdqm::Boson{

 public:
  Boson(std::vector<ewkdqm::Electron> dauVec, double metx, double mety);
  TLorentzVector FourVector(){return kine_;}
  double AssociatedMet(){return met_;}
  double IsNull(){return isNull_;}
  double Mass(){return kine_.M();}
  double Pt(){return kine_.Pt();}
  double Y() {return kine_.Rapidity();}
  double isZ(){return isZ_&&(!isW_);}
  double isW(){return isW_&&(!isZ_);}
  double isBB(){return isBB_ && (!isEE_)&&(!isEB_);}
  double isEE(){return isEE_ && (!isBB_)&&(!isEB_);}
  double isEB(){return isEB_ && (!isEE_)&&(!isBB_);}
  

 private:
  Boson();
  bool isZ_;
  bool isEE_;
  bool isEB_;
  bool isBB_;

  bool isW_;
  bool isE_;
  bool isB_;

  bool isNull_;
  double met_;
  TLorentzVector kine_;
};

ewkdqm::Boson::Boson(){
  //not allowed 
}

ewkdqm::Boson::Boson(std::vector<ewkdqm::Electron> dtrs, double metx, double mety){
  
  int nele = dtrs.size();
 
  isZ_=0;
  isW_=0;
  isB_=0;
  isE_=0;

  isBB_=0;
  isEE_=0;
  isEB_=0;
  isNull_=1;
  
  if(nele==1){
    isZ_=0;
    isW_=1;
    isB_=(dtrs[0]).InBarrel();
    isE_=(dtrs[0]).InEndcap();
    isNull_=0;
    }
  
  else if (nele>1){
    isZ_=1;
    isW_=0;
    isBB_ = (dtrs[0]).InBarrel()&&(dtrs[1]).InBarrel();
    isEE_ = (dtrs[0]).InEndcap()&&(dtrs[1]).InEndcap();
    isEB_ = (((dtrs[0]).InBarrel()&&(dtrs[1]).InEndcap())||((dtrs[0]).InEndcap()&&(dtrs[1]).InBarrel()));
    isNull_=0;
  }
 else return; 
  //calculate missing et
  met_ = sqrt(metx*metx+mety*mety);

  TLorentzVector dau1(0,0,0,0); /*Lorentz Vector to hold Z/W decay product*/
  TLorentzVector dau2(0,0,0,0); /*Lorentz Vector to hold Z/W decay product*/

  dau1.SetPxPyPzE((dtrs[0]).px(),(dtrs[0]).py(),(dtrs[0]).pz(),(dtrs[0]).energy());
  if(isZ_)
    dau2.SetPxPyPzE((dtrs[1]).px(),(dtrs[1]).py(),(dtrs[1]).pz(),(dtrs[1]).energy());
  else if(isW_)
    dau2.SetPxPyPzE(metx,metx,0,met_);
  kine_ = dau1+dau2;
}



class ewkdqm::BosonHist{

 public:
  BosonHist(std::string cutLevel,DQMStore* theDbe);
  ~BosonHist();
  void Fill(ewkdqm::Boson& b);

 private:
  BosonHist();
  MonitorElement* BosMass_;
  MonitorElement* BosPt_;
  MonitorElement* BosRap_;
  MonitorElement* AssociatedMet_;
  MonitorElement* InitHist( std::string qty,  std::string cutLevel,
                            std::string title,int nbins, int rlow,
                            int rhigh, DQMStore* theDbe);


};

ewkdqm::BosonHist::BosonHist(std::string cutLevel, DQMStore* theDbe){
  BosMass_= InitHist( "Mass_",  cutLevel,  "Mass [GeV]", 100,0.,200., theDbe);  
  BosPt_  = InitHist( "BosonTransverseMomentum_",  cutLevel,  "p_{T} (Boson) [GeV]", 100,0.,200., theDbe);
  BosRap_ = InitHist( "BosonRapidity_",  cutLevel,  "Boson Rapidity", 100,-5.,5., theDbe);
  AssociatedMet_ = InitHist("AssociatedMet", cutLevel,  "Missing transverse energy [GeV]", 100,0.,200., theDbe);
}

ewkdqm::BosonHist::~BosonHist(){
  delete BosMass_;
  delete BosPt_;
  delete BosRap_;
  delete AssociatedMet_;
}

void
ewkdqm::BosonHist::Fill(Boson& b){
  BosPt_  -> Fill(b.Pt());
  BosRap_ -> Fill(b.Y());
  BosMass_->Fill(b.Mass());
  AssociatedMet_->Fill(b.AssociatedMet());
}

MonitorElement*
ewkdqm::BosonHist::InitHist( std::string qty,  std::string cutLevel,
                                      std::string title,int nbins, int rlow,
                                      int rhigh, DQMStore* theDbe){
  MonitorElement* m;
  std::string name =qty+"_"+cutLevel;
  m = theDbe->book1D(name.c_str(),title.c_str(),nbins,rlow,rhigh);
  m->setAxisTitle("test", 1);
  return m;
}


