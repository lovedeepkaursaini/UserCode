#define TrainSample_cxx
#include "TrainSample.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include<sstream>
#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include "../Classes/particle.h"
#include "../Classes/pfjet.h"
typedef pfjet<TrainSample> jet;
namespace zjet{
  TH1D* prettyHistogram(std::string name, std::string title, std::string xtitle, double nbins, double rlow, double rhigh ){
    TH1D* hist = new TH1D(name.c_str(),title.c_str(),nbins,rlow,rhigh);
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    return hist;
  }

TH1D* prettyHistogram(std::string name, std::string title, std::string xtitle, double nbins, float* xbins){
    TH1D* hist = new TH1D(name.c_str(),title.c_str(),nbins,xbins);
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    return hist;
  }


}

void TrainSample::Loop()
{
  TFile* file = new TFile("TrainSample.root","recreate");
  TH1D* hNGenJet = zjet::prettyHistogram("hNGenJet_","N_{GenJet}(Exc) [Z+Jets]","NGenJets_{Exc}",9,1,10);
  TH1D* hNRecjet = zjet::prettyHistogram("hNRecjet_","N_{Recjet}(Exc) [Z+Jets]","NRecjets_{Exc}",9,1,10);
  TH2D* hNjetRes = new TH2D("hNjetRes","hNjetRes",9,1.,10.,9,1.,10.);

/*
  float jpt[10] = {30,50,70,90,120,150,200,250,300,400};
  float jpt2[9] = {30,50,70,90,120,150,200,250,400};
  float jpt3[8] = {30,50,70,90,120,150,200,400};
  float jpt4[6] = {30,50,70,90,120,400};

  float jpt[10] = {30,50,70,90,120,150,200,250,300,400};
  float jpt2[10] = {30,50,70,90,120,150,200,250,300,400};
  float jpt3[6] = {30,50,70,90,120,400};
  float jpt4[6] = {30,50,70,90,120,400};
*/

  float jpt[10] = {30,50,70,90,120,150,180,210,300,400};
  float jpt2[9] = {30,50,70,90,120,150,200,270,400};
  float jpt3[6] = {30,50,70,90,120,400};
  float jpt4[6] = {30,50,70,90,120,400};

 
  TH1D* Jet1Pt_ZJet=zjet::prettyHistogram("JetPt_Z1jet","LeadJet (Z+Jets)","p_{T} [GeV]",9,jpt);
  TH1D* GenJet1Pt_ZJet=zjet::prettyHistogram("GenJetPt_Z1jet","LeadGenJet (Z+Jets)","p_{T} [GeV]",9,jpt);
  TH2D* hJet1PtRes = new TH2D("hJetPtRes","hJetPtRes",9,jpt,9,jpt);

  TH1D* Jet2Pt_ZJet=zjet::prettyHistogram("Jet2Pt_Z1jet","2ndLeadJet (Z+Jets)","p_{T} [GeV]",8,jpt2);
  TH1D* GenJet2Pt_ZJet=zjet::prettyHistogram("GenJet2Pt_Z1jet","2ndLeadGenJet (Z+Jets)","p_{T} [GeV]",8,jpt2);
  TH2D* hJet2PtRes = new TH2D("hJet2PtRes","hJet2PtRes",8,jpt2,8,jpt2);

  TH1D* Jet3Pt_ZJet=zjet::prettyHistogram("Jet3Pt_Z1jet","3rdLeadJet (Z+Jets)","p_{T} [GeV]",5,jpt3);
  TH1D* GenJet3Pt_ZJet=zjet::prettyHistogram("GenJet3Pt_Z1jet","3rdLeadGenJet (Z+Jets)","p_{T} [GeV]",5,jpt3);
  TH2D* hJet3PtRes = new TH2D("hJet3PtRes","hJet3PtRes",5,jpt3,5,jpt3);

  TH1D* Jet4Pt_ZJet=zjet::prettyHistogram("Jet4Pt_Z1jet","4thLeadJet (Z+Jets)","p_{T} [GeV]",5,jpt4);
  TH1D* GenJet4Pt_ZJet=zjet::prettyHistogram("GenJet4Pt_Z1jet","4thLeadGenJet (Z+Jets)","p_{T} [GeV]",5,jpt4);
  TH2D* hJet4PtRes = new TH2D("hJet4PtRes","hJet4PtRes",5,jpt4,5,jpt4);

/*
  TH1D* Jet1Pt_ZJet=zjet::prettyHistogram("JetPt_Z1jet","LeadJet (Z+Jets)","p_{T} [GeV]",25,0.,500.);
  TH1D* GenJet1Pt_ZJet=zjet::prettyHistogram("GenJetPt_Z1jet","LeadGenJet (Z+Jets)","p_{T} [GeV]",25,0.,500.);
  TH2D* hJet1PtRes = new TH2D("hJetPtRes","hJetPtRes",25,0.,500.,25,0.,500.);

  TH1D* Jet2Pt_ZJet=zjet::prettyHistogram("Jet2Pt_Z1jet","2ndLeadJet (Z+Jets)","p_{T} [GeV]",25,0.,500.);
  TH1D* GenJet2Pt_ZJet=zjet::prettyHistogram("GenJet2Pt_Z1jet","2ndLeadGenJet (Z+Jets)","p_{T} [GeV]",25,0.,500.);
  TH2D* hJet2PtRes = new TH2D("hJet2PtRes","hJet2PtRes",25,0.,500.,25,0.,500.);

  TH1D* Jet3Pt_ZJet=zjet::prettyHistogram("Jet3Pt_Z1jet","3rdLeadJet (Z+Jets)","p_{T} [GeV]",25,0.,500.);
  TH1D* GenJet3Pt_ZJet=zjet::prettyHistogram("GenJet3Pt_Z1jet","3rdLeadGenJet (Z+Jets)","p_{T} [GeV]",25,0.,500.);
  TH2D* hJet3PtRes = new TH2D("hJet3PtRes","hJet3PtRes",25,0.,500.,25,0.,500.);

  TH1D* Jet4Pt_ZJet=zjet::prettyHistogram("Jet4Pt_Z1jet","4thLeadJet (Z+Jets)","p_{T} [GeV]",25,0.,500.);
  TH1D* GenJet4Pt_ZJet=zjet::prettyHistogram("GenJet4Pt_Z1jet","4thLeadGenJet (Z+Jets)","p_{T} [GeV]",25,0.,500.);
  TH2D* hJet4PtRes = new TH2D("hJet4PtRes","hJet4PtRes",25,0.,500.,25,0.,500.);
*/
  //-------------------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------

  TH1D* Jet1Eta_ZJet=zjet::prettyHistogram("JetEta_Z1jet","LeadJet (Z+Jets)","#eta (Jet)",30,-3.,3.);
  TH1D* GenJet1Eta_ZJet=zjet::prettyHistogram("GenJetEta_Z1jet","LeadGenJet (Z+Jets)","#eta (GenJet)",30,-3.,3.);
  TH2D* hJet1EtaRes = new TH2D("hJetEtaRes","hJetEtaRes",30,-3.,3.,30,-3.,3.);

  TH1D* Jet2Eta_ZJet=zjet::prettyHistogram("Jet2Eta_Z1jet","2ndLeadJet (Z+Jets)","#eta (Jet)",30,-3.,3.);
  TH1D* GenJet2Eta_ZJet=zjet::prettyHistogram("GenJet2Eta_Z1jet","2ndLeadGenJet (Z+Jets)","#eta (GenJet)",30,-3.,3.);
  TH2D* hJet2EtaRes = new TH2D("hJet2EtaRes","hJet2EtaRes",30,-3.,3.,30,-3.,3.);

  TH1D* Jet3Eta_ZJet=zjet::prettyHistogram("Jet3Eta_Z1jet","3rdLeadJet (Z+Jets)","#eta (Jet)",30,-3.,3.);
  TH1D* GenJet3Eta_ZJet=zjet::prettyHistogram("GenJet3Eta_Z1jet","3rdLeadGenJet (Z+Jets)","#eta (GenJet)",30,-3.,3.);
  TH2D* hJet3EtaRes = new TH2D("hJet3EtaRes","hJet3EtaRes",30,-3.,3.,30,-3.,3.);

  TH1D* Jet4Eta_ZJet=zjet::prettyHistogram("Jet4Eta_Z1jet","4thLeadJet (Z+Jets)","#eta (Jet)",30,-3.,3.);
  TH1D* GenJet4Eta_ZJet=zjet::prettyHistogram("GenJet4Eta_Z1jet","4thLeadGenJet (Z+Jets)","#eta (GenJet)",30,-3.,3.);
  TH2D* hJet4EtaRes = new TH2D("hJet4EtaRes","hJet4EtaRes",30,-3.,3.,30,-3.,3.);





  //TH1D* ZPt_ZJet=zjet::prettyHistogram("ZPt_ZJet","Z p_{T}","p_{T}(Z) [GeV]",25,0.,500.);
  //TH1D* ZRap_ZJet=zjet::prettyHistogram("ZRap_ZJet","Z Y","Y(Z)",25,0.,500.);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%10000 ==0)  std::cout<<"Event Number: "<<jentry<<std::endl; 
      // if (Cut(ientry) < 0) continue;
      if(!(genParPt_->size()))continue;
      std::vector<particle> elec;
      for(int i=0; i!=genParPt_->size();i++){
      double pt = genParPt_->at(i);
      double eta = genParEta_->at(i);
      double phi = genParPhi_->at(i);
      double e = genParE_->at(i);
      double m = 0;
      particle p(pt,eta,phi,m,e);
      if(fabs(genParId_->at(i))==11)
	if(p.PtX(20)&&p.EtaX(2.5))
	  elec.push_back(p);
      }


      int njet = patJetPfAk05Pt_->size();

      //std::cout<<genJetPt_->size()<<"\n";
      if(elec.size()!=2)continue;
    
      particle p = elec[0]+elec[1];
      bool mwinX = p.MX(75,105);

      particle z, d1, d2;
      z=p;d1=elec[0];d2=elec[1];
      if(!mwinX)continue;
   
      std::vector<particle> jets; 
      for(unsigned int i=0; i!=njet; i++){
	jet j(this, i);
	if(!(j.PtX(30)&&j.EtaX(2.4)))continue;
	if(!(j.loosePfJetIdX()))continue;
	if(j.deltaR(d1)<0.3)continue;
	if(j.deltaR(d2)<0.3)continue;
	jets.push_back(j);
       }
      if(jets.size() && jets.size()<5)
        hNRecjet ->Fill(jets.size());
      std::vector<particle> genJets;
      //Collect the gen genJets as well.
      for(int i=0; i!=genJetPt_->size();i++){
	particle jet(
		     genJetPt_->at(i),
		     genJetEta_->at(i),
		     genJetPhi_->at(i),
		     0,
		     genJetE_->at(i)
		     );
	//std::cout<<"jello1\n";
	if(!jet.PtX(30))continue;
	if(!jet.EtaX(2.4))continue;
	if(jet.deltaR(d1)<0.3)continue;
	if(jet.deltaR(d2)<0.3)continue;
//	GenJetPt_ZJet->Fill(jet.Pt());
	genJets.push_back(jet);
      }
      if(genJets.size() && genJets.size()<5) 
       hNGenJet->Fill(genJets.size());
      if(jets.size() && genJets.size()){
       hNjetRes->Fill(jets.size(),genJets.size());
      }
      if(!(jets.size()==genJets.size()))continue;

      //----------------------------------------------------

      if(jets.size()){
	Jet1Pt_ZJet->Fill((jets[0]).Pt());
	Jet1Eta_ZJet->Fill((jets[0]).Eta());
      }
      if(genJets.size()){
	GenJet1Pt_ZJet->Fill((genJets[0]).Pt());
	GenJet1Eta_ZJet->Fill((genJets[0]).Eta());
      }
      if(jets.size()&& genJets.size()){
	hJet1EtaRes->Fill((jets[0]).Eta(), (genJets[0]).Eta());
	hJet1PtRes->Fill((jets[0]).Pt(), (genJets[0]).Pt());
      }


    if(jets.size()>1){
	Jet2Pt_ZJet->Fill((jets[1]).Pt());
	Jet2Eta_ZJet->Fill((jets[1]).Eta());
      }
      if(genJets.size()>1){
	GenJet2Pt_ZJet->Fill((genJets[1]).Pt());
	GenJet2Eta_ZJet->Fill((genJets[1]).Eta());
      }
      if(jets.size()>1&& genJets.size()>1){
	hJet2EtaRes->Fill((jets[1]).Eta(), (genJets[1]).Eta());
	hJet2PtRes->Fill((jets[1]).Pt(), (genJets[1]).Pt());
      }


      if(jets.size()>2){
	Jet3Pt_ZJet->Fill((jets[2]).Pt());
	Jet3Eta_ZJet->Fill((jets[2]).Eta());
      }
      if(genJets.size()>2){
	GenJet3Pt_ZJet->Fill((genJets[2]).Pt());
	GenJet3Eta_ZJet->Fill((genJets[2]).Eta());
      }
      if(jets.size()>2&& genJets.size()>2){
	hJet3EtaRes->Fill((jets[2]).Eta(), (genJets[2]).Eta());
	hJet3PtRes->Fill((jets[2]).Pt(), (genJets[2]).Pt());
      }


   if(jets.size()>3){
	Jet4Pt_ZJet->Fill((jets[3]).Pt());
	Jet4Eta_ZJet->Fill((jets[3]).Eta());
      }
      if(genJets.size()>3){
	GenJet4Pt_ZJet->Fill((genJets[3]).Pt());
	GenJet4Eta_ZJet->Fill((genJets[3]).Eta());
      }
      if(jets.size()>3&& genJets.size()>3){
	hJet4EtaRes->Fill((jets[3]).Eta(), (genJets[3]).Eta());
	hJet4PtRes->Fill((jets[3]).Pt(), (genJets[3]).Pt());
      }

   }
  file->cd();
  file->Write();
  file->Close();
}
