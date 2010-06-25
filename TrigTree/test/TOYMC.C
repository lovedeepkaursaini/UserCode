#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

/*
  Lovedeep kaur
  Panjab University
  
  Date: Jun07, 2010
  
 "Toy MC study for estimating the 
  errors on efficiency corr factors
  using RooFit"

  Proceedure Oriented Implementation.

  Usage:
     root -l
     [0] .L TOYMC.C++
     [1] TOYMC("data.root","mc.root","WP","trigname")

*/

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "TH1D.h"
#include "TH1.h"
#include "TFile.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooBifurGauss.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TGraphAsymmErrors.h"
#include "RooFitResult.h"
using namespace RooFit;

//-------------------------------------------------------------------------------------
//----------------------------CLASS: POINT---------------------------------------------
//-------------------------------------------------------------------------------------
//A simple class to hold the data points in a graph.
//A TGraph will simply be a vector of these objects.
class point{
public:
  point(double y, double x, double er, double el):val(y),xloc(x),hErr(er),lErr(el){
  }
  point(const point& p){
    this->xloc=p.xloc;
    this->val=p.val;
    this->hErr=p.hErr;
    this->lErr=p.lErr;
  }
  ~point(){};

  double GetVal(){return val;}
  double GetLoc(){return xloc;}
  double GetHErr(){return hErr;}
  double GetLErr(){return lErr;}

  void SetVal(double val_){val=val_;}
  void SetLoc(double xl){xloc=xl;}
  void SetHErr(double hEr){hEr=hErr;}
  void SetLErr(double lEr){lEr=lErr;}

private:
  point():val(0),xloc(0),hErr(0),lErr(0){};
  double val;
  double xloc;
  double hErr;
  double lErr;
};



//-------------------------------------------------------------------------------------
//------------------------------FXN: READGRAPH ----------------------------------------
//-------------------------------------------------------------------------------------
//Given a TGraph, convert it 
//into a vector of "point" objects
std::vector<point> ReadGraph(TGraph* g1){
  std::vector<point> mypointvec;
  const int num = g1->GetN();
  for(int i=0; i<num; i++){

    double hErr = g1->GetErrorYhigh(i);
    double lErr = g1->GetErrorYlow(i);
    if(hErr<0.000000001)hErr=0.0;
    if(lErr<0.000000001)lErr=0.0;

    double val ; 
    double xloc ;
    g1->GetPoint(i,xloc,val);
    point mypoint(val,xloc,hErr,lErr);
    mypointvec.push_back(mypoint);
  }
  return mypointvec;
}

double NEVTS=20000;







//-------------------------------------------------------------------------------------
//------------------------------FXN: GENASSYGAUSS--------------------------------------
//-------------------------------------------------------------------------------------
//A function that would generate an assymetric gaussian
//distribution, provided u give it the mean and wL, wR
TH1D* GenAssyGauss(double mean, double sigL, double sigR){
  std::stringstream ss;
  ss<<"hist"<<mean<<"_"<<sigL<<"_"<<sigR;
  std::string title=ss.str();

  //RooBiFurGauss seems to have a bug: does not like 0 as sigma(L or R)
  if(!sigL)sigL=0.0000000000000001;
  if(!sigR)sigR=0.0000000000000001;

  RooRealVar effi("effi","effi",1, 0.0,1.3);
  RooRealVar meanA("mean","mean",mean);
  RooRealVar sigLA("sigL","sigL",sigL);
  RooRealVar sigRA("sigR","sigR",sigR);
 
  RooBifurGauss foo("foo","foo",effi,meanA,sigLA,sigRA);
  
  RooDataSet* data = foo.generate(effi, NEVTS);

  TH1D* hdata=(TH1D*) data->createHistogram(title.c_str(),effi);//,RooFit::Binning(100));
   return hdata;
}







//-------------------------------------------------------------------------------------
//-------------------------------FXN: DIVIDEHISTOS-------------------------------------
//-------------------------------------------------------------------------------------
//Here the division means: Generate random numbers from either
//histogram, divide them and fill the divided number into 
//a third histogram. NEVT control the number of such events
//generated

TH1D* DivideHistos(TH1D* h1, TH1D* h2,std::string trigname,std::string WP){
  
  std::string title1(h1->GetTitle());
  std::string title2(h2->GetTitle());

  std::string tit = title1+"_"+title2;

  TH1D* h1rndm = new TH1D("h1rndm","Random events",150,0.4,1.05);
  TH1D* h2rndm = new TH1D("h2rndm","Random events",150,0.4,1.05);

  TH1D* Corr = new TH1D(tit.c_str(),tit.c_str(),100,0,2);
  for(int i=0; i!=NEVTS; i++){
  
    double datai=h1->GetRandom();
    double mci=h2->GetRandom();
    if(datai>1 || mci >1)continue;
    h1rndm->Fill(datai);
    h2rndm->Fill(mci);
    //fout<<"rndm no. : data "<<datai<<" MC: "<<mci<<endl;
    double corrFac=datai/mci;
    Corr->Fill(corrFac);
  }

  std::string title11=title1+"_Random";

  TCanvas* c1=new TCanvas(title11.c_str(),title11.c_str());
  c1->cd(1);
  h2rndm->SetMarkerColor(4); 
  h2rndm->SetLineColor(4);
  h2rndm->SetLineWidth(2);  
  h2rndm->SetMarkerStyle(25);
  h2rndm->SetMarkerSize(0.8);
  h1rndm->SetMarkerColor(1);
  h1rndm->SetLineColor(1);
  h1rndm->SetLineWidth(2);
  h1rndm->SetMarkerStyle(21);  
  h1rndm->SetMarkerSize(0.8);
  h2rndm->Draw();
  h1rndm->Draw("same");
  std::string save2="SAVE/"+WP+"_"+trigname+"_"+title11+".png";
  std::string save29="SAVE/"+WP+"_"+trigname+"_"+title11+".pdf";
  c1->SaveAs(save2.c_str());
  c1->SaveAs(save29.c_str());
  return Corr;
}







//-------------------------------------------------------------------------------------
//---------------------------FXN: TOYMC------------------------------------------------
//-------------------------------------------------------------------------------------
//For a given point p1  in effi graph for data,
//and point p2 in effigraph for MC, generate the
//divided histogram "DivideHistos", containing a
//distribution of corrfactor (for that point).
//Then Fit this distribution.

void toyMc(point p1, point p2,std::string trigname,std::string WP){

  double mean1 =p1.GetVal();
  double siH1=p1.GetHErr();
  double siL1=p1.GetLErr();
  
  double mean2 =p2.GetVal();
  double siH2=p2.GetHErr();
  double siL2=p2.GetLErr();

  double loc1=p1.GetLoc();
  double loc2=p2.GetLoc();
  
  TH1D* test1 = GenAssyGauss(mean1,siL1,siH1);
  TH1D* test2 = GenAssyGauss(mean2,siL2,siH2);
  
  std::stringstream ss1,ss2;
  ss1<<"AGaus_"<<loc1<<"_"<<loc2;
  ss2<<"Divided_Corr_"<<loc1<<"_"<<loc2;
  std::string title1=ss1.str();
  std::string title2=ss2.str();
  test1->SetTitle(title1.c_str());
  TH1D* divided=DivideHistos(test1,test2,trigname,WP);

  RooRealVar corr("corr","Correction Factor",0,2);
  RooDataHist hCorr("hCorr","hCorr",corr,divided);
  double Mean=divided->GetMean();

  RooRealVar meanA("mean","mean",Mean,Mean-0.5*Mean,Mean+0.5*Mean);
  RooRealVar sigLA("sigL","sigL",0.,2.0);
  RooRealVar sigRA("sigR","sigR",0.,2.0);


  double low=0;double high=2.0;

  RooBifurGauss foo("foo","foo",corr,meanA,sigLA,sigRA);
  foo.fitTo(hCorr,RooFit::Range(low,high));
    
  std::stringstream ss3;
  ss3<<"Fitted_Corr_"<<loc1<<"_"<<loc2;
  std::string title3=ss3.str();
  TCanvas* c111=new TCanvas(title3.c_str(),title3.c_str());
  c111->cd();
  RooPlot* plot = corr.frame(RooFit::Name("Correction Factor"),RooFit::Title("Correction Factor"));
  hCorr.plotOn(plot);
  foo.plotOn(plot);
  foo.paramOn(plot,Format("NELU", AutoPrecision(2)), Layout(0.58,1.0,1.0));//0.1, 0.4,0.9));//RooFit::Layout(0.60));
  c111->cd();
  plot->Draw();
  std::string save4="SAVE/"+WP+"_"+trigname+"_"+title3+".png";
  std::string save49="SAVE/"+WP+"_"+trigname+"_"+title3+".pdf";
  c111->SaveAs(save4.c_str());
  c111->SaveAs(save49.c_str());

  cout<<"chi2---> "<<plot->chiSquare()<<endl;
  //meanA.Print("v");
  cout<<meanA.getVal()<<'\t'<<sigLA.getVal()<<'\t'<<sigRA.getVal()<<endl;

}









//-------------------------------------------------------------------------------------
//---------------------------FXN: TOY_MC_OVER_ALLPOINTS--------------------------------
//-------------------------------------------------------------------------------------
//The toyMC function provided above does the study
//for a single point, call it in a loop and do the
//study for all the points in the efficiency graphs
void toyMcOverAllPoints(std::vector<point>& gr1, std::vector<point>& gr2,std::string trigname,std::string WP){
  int steps=0;
  if(gr1.size()>gr2.size())steps=gr2.size();
  else steps=gr1.size();
  for(int i=0; i<steps; i++){
    point p1(gr1[i]);
    point p2(gr2[i]);
    if(gr1[i].GetVal()!=0 && gr2[i].GetVal()!=0 && (gr1[i].GetLoc()==gr2[i].GetLoc()))
      toyMc(p1,p2,trigname,WP);
  }
}





//-------------------------------------------------------------------------------------
//----------------------------FXN: TOYMC-----------------------------------------------
//-------------------------------------------------------------------------------------
//This is the "userDashboard", pick histos from the 
//files, create efficiency graphs...and then do
//MC study for each point in the graph.
void TOYMC(TString data,TString mc,std::string WP, std::string trigname){

  std::string ht=WP+"_"+trigname+".html";
  std::string tx=WP+"_"+trigname+".tex";

  const int BinArraySize = 14;
  Double_t Bins[BinArraySize] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 40,60};
  Int_t nBins = (Int_t) BinArraySize - 1;

  TFile *fileData=new TFile(data);
  TH1D *denomData=(TH1D*)fileData->Get((WP+"/NoTrig_HltMatchX_Et").c_str());
  TH1D *numData=(TH1D*)fileData->Get((WP+"/"+trigname+"_HltMatchX_Et").c_str());
  TH1D * denomDataNew=(TH1D*)denomData->Rebin(nBins,"denomDataNew",Bins);//Clone();//
  TH1D * numDataNew=(TH1D*)numData->Rebin(nBins,"numDataNew",Bins);
  TGraphAsymmErrors* EffData = new TGraphAsymmErrors();
  EffData->BayesDivide(numDataNew,denomDataNew);

  TFile *fileMC=new TFile(mc);
  TH1D *denomMC=(TH1D*)fileMC->Get((WP+"/NoTrig_HltMatchX_Et").c_str());
  TH1D *numMC=(TH1D*)fileMC->Get((WP+"/"+trigname+"_HltMatchX_Et").c_str());
  TH1D * denomMCNew=(TH1D*)denomMC->Rebin(nBins,"denomMCNew",Bins);
  TH1D * numMCNew=(TH1D*)numMC->Rebin(nBins,"numMCNew",Bins);
  TGraphAsymmErrors* EffMC = new TGraphAsymmErrors();
  EffMC->BayesDivide(numMCNew,denomMCNew);

   TCanvas* sh = new TCanvas("sh","effi");
  EffMC->SetMarkerColor(4);
  EffMC->SetLineColor(4);
  EffMC->SetLineWidth(2);
  EffMC->SetMarkerStyle(25);
  EffMC->SetMarkerSize(0.8);
  EffMC->SetMinimum(0.0);
  EffMC->SetMaximum(1.1);
  (EffMC->GetXaxis())->SetTitle("Et (GeV)");
  (EffMC->GetYaxis())->SetTitle("Efficiency");
   EffMC->Draw("AP");
  EffData->SetMarkerColor(1);
  EffData->SetLineColor(1);
  EffData->SetLineWidth(2);
  EffData->SetLineStyle(2);
  EffData->SetMarkerStyle(21);
  EffData->SetMarkerSize(0.8);
  EffData->SetMinimum(0.0);
  EffData->SetMaximum(1.1);
  EffData->Draw("AP");
  EffMC->Draw("Psame");
  std::string effsave="SAVE/Efficiency_"+WP+trigname+"_HltMatchX_Et.png";
  sh->SaveAs(effsave.c_str());
  std::vector<point> gr1 = ReadGraph(EffData);
  std::vector<point> gr2 = ReadGraph(EffMC);
  std::cout<<"========\n\n";
  toyMcOverAllPoints(gr1,gr2,trigname,WP);
}


//-------------------------------------------------------------------------------------
//---------------------****END*****----------------------------------------------------
//-------------------------------------------------------------------------------------
