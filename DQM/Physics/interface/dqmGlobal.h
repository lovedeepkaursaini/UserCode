#ifndef __EWKDQMGLOBAL_H__
#define __EWKDQMGLOBAL_H__

#include <iostream>
#include <TLorentzVector.h>
#include "DQM/Physics/interface/dqmWorkpoint.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"


namespace ewkdqm{
  class Global;
  class GlobalHist;
}


class ewkdqm::Global{
 public:
  Global();
  Global(double nlep, double nlepgood, double njets, double met, bool trig);
  double Nlep();
  double Njets();
  double Met();
  bool   Trig();

  bool PassTrigX();
  bool PassMetX(ewkdqm::gWorkpoint& g){return 1;}
  bool PassAllX(ewkdqm::gWorkpoint& g);

 private:
  double nlep_;
  double nlepgood_;
  double njets_;
  double met_;
  bool trig_;
  
};

ewkdqm::Global::Global():nlep_(0),nlepgood_(0),njets_(0),met_(0),trig_(0){

}


ewkdqm::Global::Global(double nlep, double nlepgood,
		       double njets, double met, 
		       bool trig):nlep_(nlep),nlepgood_(nlepgood),njets_(njets),met_(met),trig_(trig){}


double ewkdqm::Global::Nlep(){
  return nlep_;
}


double ewkdqm::Global::Njets(){
  return njets_;
}

double ewkdqm::Global::Met(){
  return met_;
}


bool ewkdqm::Global::Trig(){
  return trig_;
}

bool ewkdqm::Global::PassTrigX(){
  return trig_;
}

bool ewkdqm::Global::PassAllX(gWorkpoint& g){

  return (
    PassMetX(g)
    &&PassTrigX()
    );
}


class ewkdqm::GlobalHist{

 public:
  GlobalHist(std::string cutLevel, DQMStore* theDbe);
  ~GlobalHist();
  void Fill(ewkdqm::Global& g); 
  void FillAfterCuts(ewkdqm::Global& g, ewkdqm::gWorkpoint wp); 
  
 private:
  GlobalHist();//dont allow the default constructor.
  MonitorElement * InitHist(std::string qty,  std::string cutLevel, 
			    std::string title,int nbins, int rlow,
			    int rhigh, DQMStore* theDbe);
  
  
  //-------------------------------data members---------
  MonitorElement* trig_;
  MonitorElement* nelectrons_;
  MonitorElement* met_;
  MonitorElement* njets_;
};

//- - - - - - - - - - - - - - - - - - - - - - c o n s t r u c t o r - - - - - - - - - - - - - - - - - -

ewkdqm::GlobalHist::GlobalHist(std::string cutLevel, DQMStore* theDbe){
  trig_=InitHist( "TRIG",  cutLevel,  "Trigger response", 2,0.,2.0, theDbe);
  nelectrons_=InitHist( "NELECTRONS",     cutLevel,  "Number of electrons in event", 10,-0.5,9.5, theDbe);
  met_=InitHist("MET", cutLevel,  "Missing transverse energy [GeV]", 100,0.,200., theDbe);
  njets_=InitHist( "NJETS",  cutLevel, "Number of jets", 10,-0.5,9.5, theDbe);
 
}


//- - - - - - - - - - - - - - - - - - - - - - d e s t r u c t o r - - - - - - - - - - - - - - - - - - - - 
ewkdqm::GlobalHist::~GlobalHist(){
  
  delete trig_;
  delete nelectrons_;
  delete met_;
  delete njets_;
}


//- - - - - - - - - - - - - - - - - - - m e m b e r   f u n c t i o n s - - - - - - - - - - - - - - - - - - -
MonitorElement* 
ewkdqm::GlobalHist::InitHist( std::string qty,  std::string cutLevel, 
				      std::string title,int nbins, int rlow,
				      int rhigh, DQMStore* theDbe){
  MonitorElement* m;
  std::string name =qty+"_"+cutLevel;
  m = theDbe->book1D(name.c_str(),title.c_str(),nbins,rlow,rhigh);
  return m;
}



void
ewkdqm::GlobalHist::Fill(ewkdqm::Global& g){
  trig_->Fill(g.Trig());
  nelectrons_->Fill(g.Nlep());
  met_->Fill(g.Met());
  njets_->Fill(g.Njets());
}

void 
ewkdqm::GlobalHist::FillAfterCuts(ewkdqm::Global& g, ewkdqm::gWorkpoint wp){

  bool metx   = g.PassMetX(wp);
  bool trigx  = g.PassTrigX();

  if(metx && trigx){
    nelectrons_->Fill(g.Nlep());
    njets_->Fill(g.Njets());
    //also fill goodnele
} 
}

//--------------------------
//--------------------------
#endif
