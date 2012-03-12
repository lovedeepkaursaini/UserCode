#ifndef _WORKPOINT_H__
#define _WORKPOINT_H__



#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace ewkdqm{
  class  gWorkpoint;
  class  eWorkpoint;
  class  jWorkpoint;
}



//------------------------------------------------------------------
//------------------------------------------------------------------

class ewkdqm::gWorkpoint{
 public:
  gWorkpoint(const edm::ParameterSet& cfg);
  int njetMaxX_;
  int neleMaxX_;
  double metMaxX_;
  double massMaxX_;


  int njetMinX_;
  int neleMinX_;
  double metMinX_;
  double massMinX_;
  bool trigX_;

 private:
  gWorkpoint();
};

ewkdqm::gWorkpoint::gWorkpoint(const edm::ParameterSet& cfg):
  njetMaxX_(cfg.getUntrackedParameter<int>("NJetMax", 999999)),
  neleMaxX_(cfg.getUntrackedParameter<int>("NEleMax", 999999)),
  metMaxX_ (cfg.getUntrackedParameter<double>("MetMax", 999999)),
  massMaxX_(cfg.getUntrackedParameter<double>("MassMax", 999999)),
  
  njetMinX_(cfg.getUntrackedParameter<int>("NJetMin", 0)),
  neleMinX_(cfg.getUntrackedParameter<int>("NEleMin", 0)),
  metMinX_ (cfg.getUntrackedParameter<double>("MetMin", 0)),
  massMinX_(cfg.getUntrackedParameter<double>("MassMin", 0)),
  trigX_(1)
{}




//------------------------------------------------------------------
//------------------------------------------------------------------

class ewkdqm::eWorkpoint{
 public:
  eWorkpoint(const edm::ParameterSet& cfg);
  
  double   ptX_;
  double   etaX_;
  
  double   sIeIeX_barrel_;
  double   dEtaInX_barrel_;
  double   eIsoX_barrel_;
  double   hIsoX_barrel_;
  double   tIsoX_barrel_;
  
  double   sIeIeX_endcap_;
  double   dEtaInX_endcap_;
  double   eIsoX_endcap_;
  double   hIsoX_endcap_;
  double   tIsoX_endcap_;
  
 private:
  eWorkpoint();
};


ewkdqm::eWorkpoint::eWorkpoint(const edm::ParameterSet& cfg):
ptX_(cfg.getUntrackedParameter<double>("ElePtCut", 10.)),
  etaX_(cfg.getUntrackedParameter<double>("EleEtaCut", 2.4)),
  sIeIeX_barrel_(cfg.getUntrackedParameter<double>("SieieBarrel", 0.01)),
  dEtaInX_barrel_(cfg.getUntrackedParameter<double>("DetainBarrel", 0.0071)),
  eIsoX_barrel_(cfg.getUntrackedParameter<double>("EcalIsoCutBarrel", 5.7)),
  hIsoX_barrel_(cfg.getUntrackedParameter<double>("HcalIsoCutBarrel", 8.1)),
  tIsoX_barrel_(cfg.getUntrackedParameter<double>("TrkIsoCutBarrel", 7.2)),
  sIeIeX_endcap_(cfg.getUntrackedParameter<double>("SieieEndcap", 0.028)),
  dEtaInX_endcap_(cfg.getUntrackedParameter<double>("DetainEndcap", 0.0066)),
  eIsoX_endcap_(cfg.getUntrackedParameter<double>("EcalIsoCutEndcap", 5.0)),
  hIsoX_endcap_(cfg.getUntrackedParameter<double>("HcalIsoCutEndcap", 3.4)),
  tIsoX_endcap_(cfg.getUntrackedParameter<double>("TrkIsoCutEndcap", 5.1))
{}




//------------------------------------------------------------------
//------------------------------------------------------------------


//------------------------------------------------------------------

class ewkdqm::jWorkpoint{
 public:
  jWorkpoint(const edm::ParameterSet& cfg);
  double   ptX_;
  double   etaX_;
};

ewkdqm::jWorkpoint::jWorkpoint(const edm::ParameterSet& cfg):
  ptX_(cfg.getUntrackedParameter<double>("JetPtCut", 10.)),
  etaX_(cfg.getUntrackedParameter<double>("JetEtaCut", 2.4))
  {}


#endif
