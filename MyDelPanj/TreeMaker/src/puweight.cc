/*
Lovedeep Kaur Saini
Panjab University, 
Chandigarh
*/
#include <iostream>
#include "DelPanj/TreeMaker/interface/puweight.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

puweight::puweight(std::string name, TTree* tree){
  name_=name;
  tree_=tree;
  SetBranches();
Double_t Fall2011[50] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08
  };
double Data[50]={2.44414e+06,1.25579e+06,3.42957e+06,5.18618e+07,2.54054e+08,4.36586e+08,4.86031e+08,4.63551e+08,4.18993e+08,3.84891e+08,3.65724e+08,3.41505e+08,3.22711e+08,3.0886e+08,2.87693e+08,2.5129e+08,1.99438e+08,1.40551e+08,8.66577e+07,4.6234e+07,2.12058e+07,8.37396e+06,2.88178e+06,882886,246537,63900.8,15571.6,3628.24,840.61,211.248,66.7507,32.1445,26.92,33.3738,47.1181,66.9628,92.522,123.311,158.278,195.606,232.736,266.603,294.026,312.195,319.142,314.095,297.617,271.501,238.455,1.93299e+07};

 for( int i=0; i<50; ++i) {
      Data2011_.push_back(Data[i]);
      Fall2011_.push_back(Fall2011[i]);
   }
 LumiWeights_ = edm::LumiReWeighting(Fall2011_,Data2011_);



}


puweight::~puweight(){
  delete tree_;
}


void
puweight::Fill(const edm::Event& iEvent,bool data){
  Clear();

if(data) PUweight_ = -999999999999.;
else {
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

  std::vector<PileupSummaryInfo>::const_iterator PVI;


  float npT=-1.;
  float npIT=-1.;


  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

    int BX = PVI->getBunchCrossing();

    if(BX == 0) {
      npT = PVI->getTrueNumInteractions();
      npIT = PVI->getPU_NumInteractions();
    }
  }

  double MyWeight = LumiWeights_.weight( npT );
PUweight_ = MyWeight;
}
}

void
puweight::SetBranches(){
  AddBranch(&PUweight_,  "weight");

}



void
puweight::AddBranch(double* x, std::string name){
  std::string brName="PU_"+name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());
}



void 
puweight::Clear(){
PUweight_=-99999.;
}

