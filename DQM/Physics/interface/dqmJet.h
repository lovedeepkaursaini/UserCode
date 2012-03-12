#ifndef __JET_HIST__
#define __JET_HIST__
#include <TLorentzVector.h>

namespace ewkdqm{
  class JetHist;
}

class ewkdqm::JetHist{
 public:
  JetHist(std::string cutLevel, DQMStore* theDbe);
  ~JetHist();
  void Fill(TLorentzVector& j);
 private:
  JetHist();
  MonitorElement* InitHist( std::string qty,  std::string cutLevel,
                            std::string title,int nbins, int rlow,
                            int rhigh, DQMStore* theDbe);

  MonitorElement* pt_;
  MonitorElement* eta_;

};


ewkdqm::JetHist::JetHist(){
}


ewkdqm::JetHist::~JetHist(){
  delete pt_;
  delete eta_;
}


ewkdqm::JetHist::JetHist(std::string cutLevel, DQMStore* theDbe){

  pt_  =  InitHist( "PT",cutLevel,"Jet transverse momentum [GeV]",100,0.,100.,theDbe);
  eta_ =  InitHist("ETA",cutLevel,"Jet pseudo-rapidity",100,-3,3, theDbe);
}

void ewkdqm::JetHist::Fill(TLorentzVector& j){
  pt_->Fill(j.Pt());
  eta_->Fill(j.Eta());
}


MonitorElement* ewkdqm::JetHist::InitHist(std::string qty,  std::string cutLevel,
                                   std::string title,int nbins, int rlow,
                                   int rhigh, DQMStore* theDbe){
  MonitorElement* m;
  std::string name =qty+"_"+cutLevel;
  m = theDbe->book1D(name.c_str(),title.c_str(),nbins,rlow,rhigh);
  return m;
}


#endif

