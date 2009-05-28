#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TLegend.h>
#include <TROOT.h>
#include "JetUtilMC.h"

// PYTHIA pT-hat binning
const int nGen=10;
double weight[nGen];
int nEvents[nGen];
float lumi=100;
double crosssection[nGen];
// From DBS
crosssection[0]  =  6430.0;
crosssection[1]  =  230.0 ;
crosssection[2]  =  211.0 ;
crosssection[3]  =  142.0 ;
crosssection[4]  =  56.8  ;
crosssection[5]  =  18.8  ;
crosssection[6]  =  5.4   ;
crosssection[7]  =  1.55  ;
crosssection[8]  =  0.45  ;
crosssection[9]  =  0.20  ;

TString inFileNames[nGen] = {"Summer08-ZeeJets_Pt_0_15.root",
                             "Summer08-ZeeJets_Pt_15_20.root",
                             "Summer08-ZeeJets_Pt_20_30.root",
                             "Summer08-ZeeJets_Pt_30_50.root",
                             "Summer08-ZeeJets_Pt_50_80.root",
                             "Summer08-ZeeJets_Pt_80_120.root",
                             "Summer08-ZeeJets_Pt_120_170.root",
                             "Summer08-ZeeJets_Pt_170_230.root",
                             "Summer08-ZeeJets_Pt_230_300.root",
                             "Summer08-ZeeJets_Pt_300_Inf.root"};

static const unsigned int NUM_JET_MAX= 10 ;

float R1[NUM_JET_MAX] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
float R2[NUM_JET_MAX] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
float p[NUM_JET_MAX] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
 
// branch variables
Float_t eMinusPt;
Float_t eMinusEta;
Float_t eMinusPhi;
Float_t         eMinus_DeltaEtaIn;
Float_t         eMinus_DeltaPhiIn;
Float_t         eMinus_SigmaEtaEta;
Float_t         eMinus_EoverPin;
Float_t eMinusE;

Float_t ePlusPt;
Float_t ePlusEta;
Float_t ePlusPhi;
Float_t         ePlus_DeltaEtaIn;
Float_t         ePlus_DeltaPhiIn;
Float_t         ePlus_SigmaEtaEta;
Float_t         ePlus_EoverPin;

Int_t           NumRecoJetAlgo;
Int_t           NumRecoJets;
Float_t         JetRecoPt[10][10];
Float_t         JetRecoEta[10][10];
Float_t         JetRecoPhi[10][10];

Float_t eMinusPtGen;
Float_t eMinusEtaGen;
Float_t eMinusPhiGen;
Float_t ePlusEtaGen;
Float_t ePlusPtGen;
Float_t ePlusPhiGen;
Int_t           NumGenJetAlgo;
Int_t           NumGenJets;
Float_t         JetGenPt[10][10];
Float_t         JetGenEta[10][10];
Float_t         JetGenPhi[10][10];

Float_t         mZeeGen;
Float_t         Z_PtGen;
Float_t         Z_EtaGen;
Float_t         Z_PhiGen;

Float_t         mZee;
Float_t         Z_Pt;
Float_t         Z_Eta;
Float_t         Z_Phi;

void ZJet_selection()
{
  
  // ****** get the weight value for each pthat bin ******* //
  for(int i=0; i<nGen; i++) 
    {
      TFile f( inFileNames[i], "read");
      TTree* atree = (TTree*) f.Get("ZJet");
      nEvents[i]  = (int) atree->GetEntries();
      cout<<"nEvents_"<<i<<" == "<<nEvents[i]<<endl;
      //if (applyWeights) 
      weight[i] = lumi * crosssection[i]  /  nEvents[i];
      // else weight[i] = 1.0;
      cout << "weight_" << i << " == " << weight[i] << endl;
    }
  
  // Chain all pT bins etc .
  TChain *mychain = new TChain("ZJet") ;
  mychain->Add("Summer08-ZeeJets_Pt_0_15.root") ;
  mychain->Add("Summer08-ZeeJets_Pt_15_20.root") ;
  mychain->Add("Summer08-ZeeJets_Pt_20_30.root") ;
  mychain->Add("Summer08-ZeeJets_Pt_30_50.root") ;
  mychain->Add("Summer08-ZeeJets_Pt_50_80.root") ;
  mychain->Add("Summer08-ZeeJets_Pt_80_120.root") ;
  mychain->Add("Summer08-ZeeJets_Pt_120_170.root") ;
  mychain->Add("Summer08-ZeeJets_Pt_170_230.root") ;
  mychain->Add("Summer08-ZeeJets_Pt_230_300.root") ;
  mychain->Add("Summer08-ZeeJets_Pt_300_Inf.root") ;
  
  mychain->SetBranchAddress("eMinusPt",&eMinusPt);
  mychain->SetBranchAddress("eMinusEta",&eMinusEta);
  mychain->SetBranchAddress("eMinusPhi",&eMinusPhi);
  mychain->SetBranchAddress("eMinus_DeltaEtaIn",&eMinus_DeltaEtaIn);
  mychain->SetBranchAddress("eMinus_DeltaPhiIn",&eMinus_DeltaPhiIn);
  mychain->SetBranchAddress("eMinus_SigmaEtaEta",&eMinus_SigmaEtaEta);
  mychain->SetBranchAddress("eMinus_EoverPin",&eMinus_EoverPin);
  mychain->SetBranchAddress("eMinusE",&eMinusE);
  
  mychain->SetBranchAddress("ePlusPt",&ePlusPt);
  mychain->SetBranchAddress("ePlusEta",&ePlusEta);
  mychain->SetBranchAddress("ePlusPhi",&ePlusPhi);
  mychain->SetBranchAddress("ePlus_DeltaEtaIn",&ePlus_DeltaEtaIn);
  mychain->SetBranchAddress("ePlus_DeltaPhiIn",&ePlus_DeltaPhiIn);
  mychain->SetBranchAddress("ePlus_SigmaEtaEta",&ePlus_SigmaEtaEta);
  mychain->SetBranchAddress("ePlus_EoverPin",&ePlus_EoverPin);
  
  mychain->SetBranchAddress("NumRecoJetAlgo",&NumRecoJetAlgo);
  mychain->SetBranchAddress("NumRecoJets",&NumRecoJets);
  mychain->SetBranchAddress("JetRecoPt",&JetRecoPt);
  mychain->SetBranchAddress("JetRecoEta",&JetRecoEta);
  mychain->SetBranchAddress("JetRecoPhi",&JetRecoPhi);
  
  mychain->SetBranchAddress("eMinusPtGen",&eMinusPtGen);
  mychain->SetBranchAddress("eMinusEtaGen",&eMinusEtaGen);
  mychain->SetBranchAddress("eMinusPhiGen",&eMinusPhiGen);
  mychain->SetBranchAddress("ePlusPtGen",&ePlusPtGen);
  mychain->SetBranchAddress("ePlusEtaGen",&ePlusEtaGen);
  mychain->SetBranchAddress("ePlusPhiGen",&ePlusPhiGen);
  mychain->SetBranchAddress("NumGenJetAlgo",&NumGenJetAlgo);
  mychain->SetBranchAddress("NumGenJets",&NumGenJets);
  mychain->SetBranchAddress("JetGenPt",&JetGenPt);
  mychain->SetBranchAddress("JetGenEta",&JetGenEta);
  mychain->SetBranchAddress("JetGenPhi",&JetGenPhi);
  
  mychain->SetBranchAddress("mZeeGen",&mZeeGen);
  mychain->SetBranchAddress("Z_PtGen",&Z_PtGen);
  mychain->SetBranchAddress("Z_PhiGen",&Z_PhiGen);
  mychain->SetBranchAddress("Z_EtaGen",&Z_EtaGen);
  
  mychain->SetBranchAddress("mZee",&mZee);
  mychain->SetBranchAddress("Z_Pt",&Z_Pt);
  mychain->SetBranchAddress("Z_Phi",&Z_Phi);
  mychain->SetBranchAddress("Z_Eta",&Z_Eta);
  
  mychain->SetBranchStatus("*", 0); 
  mychain->SetBranchStatus("eMinusPt",1);
  mychain->SetBranchStatus("eMinusEta",1);
  mychain->SetBranchStatus("eMinusPhi",1);
  mychain->SetBranchStatus("eMinus_DeltaEtaIn",1);
  mychain->SetBranchStatus("eMinus_DeltaPhiIn",1);
  mychain->SetBranchStatus("eMinus_SigmaEtaEta",1);
  mychain->SetBranchStatus("eMinus_EoverPin",1);
  mychain->SetBranchStatus("eMinusE",1);
    
  mychain->SetBranchStatus("ePlusPt",1);
  mychain->SetBranchStatus("ePlusEta",1);
  mychain->SetBranchStatus("ePlusPhi",1);
  mychain->SetBranchStatus("ePlus_DeltaEtaIn",1);
  mychain->SetBranchStatus("ePlus_DeltaPhiIn",1);
  mychain->SetBranchStatus("ePlus_SigmaEtaEta",1);
  mychain->SetBranchStatus("ePlus_EoverPin",1);
  
  mychain->SetBranchStatus("NumRecoJetAlgo",1);
  mychain->SetBranchStatus("NumRecoJets",1);
  mychain->SetBranchStatus("JetRecoPt",1);
  mychain->SetBranchStatus("JetRecoEta",1);
  mychain->SetBranchStatus("JetRecoPhi",1);
  
  mychain->SetBranchStatus("eMinusPtGen",1);
  mychain->SetBranchStatus("eMinusEtaGen",1);
  mychain->SetBranchStatus("eMinusPhiGen",1);
  mychain->SetBranchStatus("ePlusPtGen",1);
  mychain->SetBranchStatus("ePlusEtaGen",1);
  mychain->SetBranchStatus("ePlusPhiGen",1);
  mychain->SetBranchStatus("NumGenJetAlgo",1);
  mychain->SetBranchStatus("NumGenJets",1);
  mychain->SetBranchStatus("JetGenPt",1);
  mychain->SetBranchStatus("JetGenEta",1);
  mychain->SetBranchStatus("JetGenPhi",1);
  
  mychain->SetBranchStatus("mZeeGen",1);
  mychain->SetBranchStatus("Z_PtGen",1);
  mychain->SetBranchStatus("Z_PhiGen",1);
  mychain->SetBranchStatus("Z_EtaGen",1);
  
  mychain->SetBranchStatus("mZee",1);
  mychain->SetBranchStatus("Z_Pt",1);
  mychain->SetBranchStatus("Z_Phi",1);
  mychain->SetBranchStatus("Z_Eta",1);
  
  const double MIN = 60.0;
  const double MAX = 120.0;
  const int BINz = 60;
  const int BINs = 40;
  
  TFile* rootfile=new TFile("ZeeJet.root","recreate");
  rootfile->mkdir("Reco");
  rootfile->cd("Reco");
  
  TH1F* jet_Pt_reco=new TH1F("jet_Pt_reco","jet_Pt_reco",BINs,0.,400.);
  TH1F* jet_Eta_reco=new TH1F("jet_Eta_reco","jet_Eta_reco",BINs,-5.,5.);
  TH1F* jet_Phi_reco=new TH1F("jet_Phi_reco","jet_Phi_reco",BINs,-3.2,3.2);
  
  TH1F* Z_M_reco=new TH1F("Z_M_reco","Z_M_reco",BINz,60.,120.);
  TH1F* Z_Pt_reco=new TH1F("Z_Pt_reco","Z_Pt_reco",BINs,0.,400.);
  TH1F* Z_Eta_reco=new TH1F("Z_Eta_reco","Z_Eta_reco",BINs,-5.,5.);
  TH1F* Z_Phi_reco=new TH1F("Z_Phi_reco","Z_Phi_reco",BINs,-3.2,3.2);
  
  TH1F* eM_Pt_reco=new TH1F("eM_Pt_reco","eM_Pt_reco",BINs,0.,200.);
  TH1F* eM_Eta_reco=new TH1F("eM_Eta_reco","eM_Eta_reco",BINs,-2.5,2.5);
  TH1F* eM_Phi_reco=new TH1F("eM_Phi_reco","eM_Phi_reco",BINs,-3.2,3.2);
  
  TH1F* eP_Pt_reco=new TH1F("eP_Pt_reco","eP_Pt_reco",BINs,0.,200.);
  TH1F* eP_Eta_reco=new TH1F("eP_Eta_reco","eP_Eta_reco",BINs,-2.5,2.5);
  TH1F* eP_Phi_reco=new TH1F("eP_Phi_reco","eP_Phi_reco",BINs,-3.2,3.2);
  
  TH1F* jet_Pt_reco_id=new TH1F("jet_Pt_reco_id","jet_Pt_reco_id",BINs,0.,400.);
  TH1F* Z_M_reco_id=new TH1F("Z_M_reco_id","Z_M_reco_id",30,75.,105.);
  TH1F* Z_Pt_reco_id=new TH1F("Z_Pt_reco_id","Z_Pt_reco_id",BINs,0.,400.);

  for (Long64_t entry =0; entry < mychain->GetEntries(); entry++) {
    
    mychain->GetEntry(entry);
    double wt = GetWeight( mychain->GetFile()->GetName() );
    if(entry%10000==0)cout<<"entry: "<<entry<<endl;
    int leadJetReco=0; int secJetReco=0; int leadJetGen=0; int secJetGen=0;
    int l=FindLeadIndices( 0, leadJetReco, secJetReco,leadJetGen,secJetGen);
    
    //be careful.....
    if(leadJetReco==-1)continue;
    if(JetRecoPt[0][leadJetReco]==0)continue;
    if(fabs(JetRecoPhi[0][leadJetReco])>10)continue;
    
    //Jets ....
    jet_Pt_reco->Fill(JetRecoPt[0][leadJetReco],wt);
    jet_Eta_reco->Fill(JetRecoEta[0][leadJetReco],wt);
    jet_Phi_reco->Fill(JetRecoPhi[0][leadJetReco],wt);
    
    //Z properties....
    Z_M_reco->Fill(mZee,wt);
    Z_Pt_reco->Fill(Z_Pt,wt);
    Z_Eta_reco->Fill(Z_Eta,wt);
    Z_Phi_reco->Fill(Z_Phi,wt);
    
    //electrons.......
    eM_Pt_reco->Fill(eMinusPt,wt);
    eM_Eta_reco->Fill(eMinusEta,wt);
    eM_Phi_reco->Fill(eMinusPhi,wt);
      
    eP_Pt_reco->Fill(ePlusPt,wt);
    eP_Eta_reco->Fill(ePlusEta,wt);
    eP_Phi_reco->Fill(ePlusPhi,wt);	 
    
    //Defive variables..................
    double dPhiR=dPhi(Z_Phi,JetRecoPhi[0][leadJetReco]);
    double ratioR=(JetRecoPt[0][secJetReco])/(JetRecoPt[0][leadJetReco]);
    double ratio2R=(JetRecoPt[0][secJetReco])/(Z_Pt);
    double dPhiEPMR=dPhi(eMinusPhi,ePlusPhi);
    
    //For my ease, defined all cuts....
    bool ePM_Pt=ePlusPt>25.0 && eMinusPt> 25.0;
    bool eP_Eta= fabs(ePlusEta)<1.4442 ||(fabs(ePlusEta)>1.560 && fabs(ePlusEta)<2.5);
    bool eM_Eta=fabs(eMinusEta)<1.4442 || (fabs(eMinusEta)>1.560 && fabs(eMinusEta)<2.5);
    bool zm= fabs(mZee-91.2) < 10.0 ;
    bool secJ= ratio2R < 0.2 ;
    bool leadJEta= fabs(JetRecoEta[0][leadJetReco]) < 1.3 ;
    bool dPhi= fabs(dPhiR) > 2.94 ;
    bool eM_SEE=eMinus_SigmaEtaEta<0.03 ;
    bool eP_SEE=ePlus_SigmaEtaEta<0.03 ;
    bool eM_EP=eMinus_EoverPin<3 ;
    bool eP_EP=ePlus_EoverPin<3 ;
    bool eM_dPhi=fabs(eMinus_DeltaPhiIn)<0.06 ;
    bool eP_dPhi=fabs(ePlus_DeltaPhiIn)<0.06 ;
    bool eM_dEta=fabs(eMinus_DeltaEtaIn)<0.01 ;
    bool eP_dEta=fabs(ePlus_DeltaEtaIn)<0.01;
    
    //use the function of cuts....
    float leadJetPtR=JetRecoPt[0][leadJetReco];
    float leadJetEtaR=JetRecoEta[0][leadJetReco];
    float seclJetPtR=JetRecoPt[0][secJetReco];
    bool passR=BoolCutResult(Z_Pt, mZee, eMinusPt,
			     eMinusEta, ePlusPt, ePlusEta,
			     leadJetPtR,leadJetEtaR,seclJetPtR,
			     dPhi,
			     eMinus_SigmaEtaEta,ePlus_SigmaEtaEta,
			     eMinus_EoverPin,ePlus_EoverPin,
			     eMinus_DeltaPhiIn,ePlus_DeltaPhiIn,
			     eMinus_DeltaEtaIn,ePlus_DeltaEtaIn);
    
    if(!passR)continue;
  
    Z_M_reco_id->Fill(mZee,wt);
    Z_Pt_reco_id->Fill(Z_Pt,wt);
    jet_Pt_reco_id->Fill(JetRecoPt[0][leadJetReco],wt);
	
  }
  rootfile->Write();
  rootfile->Close();
}

// //  ******************************************** //

double GetWeight(TString filename){
  double wt =0.0;
  for(int i=0; i<nGen; i++) {
    TString st2 = inFileNames[i];
    if(  filename. CompareTo(st2) == 0) wt = weight[i];
  }

  return wt;
}
// // ******************************************************* //

// ////////// Apply event selection cuts ///////////////////

bool BoolCutResult( float zPt, float mZ, float e1Pt,
		    float e1Eta, float e2Pt, float e2Eta,
		    float leadPt, float leadEta, float secondPt,
		    float phidiff, 
		    float eMinus_SigmaEtaEta,float ePlus_SigmaEtaEta,
		    float eMinus_EoverPin,float ePlus_EoverPin,
		    float eMinus_DeltaPhiIn,float ePlus_DeltaPhiIn,
		    float eMinus_DeltaEtaIn,float ePlus_DeltaEtaIn) 
{
  bool result = true;
  //  // Z mass cut
  if( fabs(mZ-91.2) > 10.0 ) result = false;
  //  // electron pT cut
  if( e1Pt < 25.0 ) result = false;
  if( e2Pt < 25.0 ) result = false;
  //  // electron acceptance
  if( !((fabs(e1Eta)<1.4442) ||
	(fabs(e1Eta)>1.560 && fabs(e1Eta)<2.5)) ) result = false;
  if( !((fabs(e2Eta)<1.4442) ||
	(fabs(e2Eta)>1.560 && fabs(e2Eta)<2.5)) ) result = false;
  
  // //   // electron isolation
  // // //   if( e1trackiso > 0.2 )  result = false;
  // // //   if( e2trackiso > 0.2 )  result = false;
  // // //   if( !(e1ecaliso > 0.0 && e1ecaliso < 0.2) ) result = false;
  // // //   if( !(e2ecaliso > 0.0 && e2ecaliso < 0.2) ) result = false;
  
  // electron iD cuts
  if( !(eMinus_SigmaEtaEta<0.03 ))  result = false;
  if( !(ePlus_SigmaEtaEta<0.03 ))  result = false;

  if( !(eMinus_EoverPin<3. ))  result = false;
  if( !(ePlus_EoverPin<3. ))  result = false;
  
  if( !(fabs(eMinus_DeltaPhiIn)<0.06 ))  result = false;
  if( !(fabs(ePlus_DeltaPhiIn)<0.06 ))  result = false;
  
  if( !(fabs(eMinus_DeltaEtaIn)<0.01 ))  result = false;
  if( !(fabs(ePlus_DeltaEtaIn)<0.01 ))  result = false;

  // //   // Minimum jet pT cut on all jets
  // //   if( leadPt < MINPTCUT ) result = false;
  
  // //   // Cut on the second jet pT
  // //   //  if( secondPt / leadPt > 0.2 ) result = false;
  if( (secondPt / zPt) > 0.2 ) result = false;
  
  // //   // Leading jet has to be in the barrel
  if( fabs(leadEta) > 1.3 ) result = false;
  
  // //   // Z and lead jet are back-to-back in phi
  if( fabs(phidiff) < 2.94 ) result = false;
  
  return result;
}

 void FindLeadIndices( int algo, int &leadReco, int &secReco,
		       int &leadGen, int &secGen ) {
   
   FindLeadIndexReco( algo, leadReco, secReco );
   // FindLeadIndexGen( algo, leadGen, secGen );
 }
 
void FindLeadIndexReco( int algo, int &lead, int &sec ) {
  
  // calculate dR of all the jets w.r.t. the two electrons
  for( int j=0; j < NumRecoJets; j++ ) {
    p[j]  = JetRecoPt[algo][j];
    
    // protection:
    if(p[j]==0 || (p[j]/p[j])!=1)continue;
    //this is giving nan as nan is not a number so divison by nan is not defined.
    
    if( p[j] < 0.001 ||p[j]>1000000) { p[j] = 0.001; break; }
    
    R1[j] = radius( JetRecoEta[algo][j], JetRecoPhi[algo][j],
		    eMinusEta,  eMinusPhi );
    R2[j] = radius( JetRecoEta[algo][j], JetRecoPhi[algo][j],
		    ePlusEta,  ePlusPhi );
  }
  
  FindIndex( p, R1, R2, lead, sec, NumRecoJets );
}
 
 void FindLeadIndexGen( int algo, int &lead, int &sec ) {
   
   for( int j=0; j < NumGenJets; j++ ) {
     p[j]  = JetGenPt[algo][j];
     
     // protection:
     if(p[j]==0 ||(p[j]/p[j])!=1)continue;
     //this is giving nan as nan is not a number so divison by nan is not defined.
     
     if( p[j] < 0.001||p[j]>1000000 ) { p[j] = 0.001; break; }
     
     R1[j] = radius( JetGenEta[algo][j], JetGenPhi[algo][j],
		     eMinusEtaGen,  eMinusPhiGen );
     R2[j] = radius( JetGenEta[algo][j], JetGenPhi[algo][j],
		     ePlusEtaGen,  ePlusPhiGen );
   }
   
   FindIndex( p, R1, R2, lead, sec, NumGenJets );
 }
 
 void FindIndex( float pT[], float dr1[], float dr2[], int &lead,
		 int &sec, int maxIndx = NUM_JET_MAX ) {
   float max = 0.0;
   lead = -1;
   sec = -1;
   for (int i=0; i< maxIndx; i++) {
     if( !(dr1[i]>0.1 && dr2[i]>0.1) ) continue;
     if( pT[i]>max ) {  max = pT[i];  lead = i; }
   }

   if( lead == -1 ) return;
   
   max = 0.0;
   for (int i=0; i< maxIndx; i++) {
     if( !(dr1[i]>0.1 && dr2[i]>0.1) || i==lead) continue;
     if( pT[i]>max ) {  max = pT[i];  sec = i; }
   }
 }
 
