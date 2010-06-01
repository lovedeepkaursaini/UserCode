//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 25 16:59:38 2010 by ROOT version 5.22/00d
// from TTree TrigTree/TrigTree
// found on file: TrigTree_data.root
//////////////////////////////////////////////////////////

#ifndef trigEffi_h
#define trigEffi_h
#include <vector>
#include <string>
#include <map>
#include "TLorentzVector.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "histo.C"
#include "TDirectory.h"

using namespace std;





class xyz{
 public:
  xyz(double x_, double y_, double z_,double HoE_, double Sihih_, double nBC_):
    x(x_),y(y_),z(z_),HoE(HoE_),Sihih(Sihih_),nBC(nBC_){}
    ~xyz(){}
    double Pt(){return x;}
    double Eta(){return y;}
    double Phi(){return z;}
    double HoverE() {return HoE;}
    double sigmaIetaIeta() {return Sihih;}
    double nBasicClusters() {return nBC;}
    double x;
    double y;
    double z;
    double HoE;
    double Sihih;
    double nBC;
 private:
    xyz();
};

typedef	std::map<std::string,std::vector<xyz> > scColl;

class trigEffi {
public :
  std::vector<xyz>
    SelectElectrons(std::string trigName, std::string cutString, scWorkPoint* scWP);
  std::vector<std::string> 
    cutParser(std::string cutString);
  
  std::map<std::string, bool> trigStuff;
  //histo* test;
  //histo* test2;
  //histo* test3;
  //histo* test4;
  //  scWorkPoint *scWP1;

  void Loop(std::vector<std::string> selVec, scWorkPoint* scWP1, TFile* file);


  scColl MkCollection(std::string trigName, std::vector<std::string> selecVec,scWorkPoint* scWP);
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leaf types
  Int_t           EvtInfo_EventNum;
  Int_t           EvtInfo_RunNum;
  Int_t           EvtInfo_LumiSection;
  Int_t           EvtInfo_BunchXing;
  vector<bool>    *trigResults;
  vector<string>  *trigName;
  Double_t        scEBNum;
  vector<double>  *scEBEnergy;
  vector<double>  *scEBEt;
  vector<double>  *scEBEta;
   vector<double>  *scEBPhi;
   vector<double>  *scEBRawEnergy;
   vector<double>  *scEBpreshowerEnergy;
   vector<double>  *scEBenergy1;
   vector<double>  *scEBenergy2x2;
   vector<double>  *scEBenergy3x3;
   vector<double>  *scEBenergy1x5;
   vector<double>  *scEBenergy2x5;
   vector<double>  *scEBenergy5x5;
   vector<double>  *scEBHoverE;
   vector<double>  *scEBx;
   vector<double>  *scEBy;
   vector<double>  *scEBz;
   vector<double>  *scEBetaWidth;
   vector<double>  *scEBphiWidth;
   vector<double>  *scEBsigmaetaeta;
   vector<double>  *scEBsigmaIetaIeta;
   vector<double>  *scEBsigmaphiphi;
   vector<double>  *scEBsigmaIphiIphi;
   vector<double>  *scEBsize;
   vector<double>  *scEBnBasicClusters;
   Double_t        scEENum;
   vector<double>  *scEEEnergy;
   vector<double>  *scEEEt;
   vector<double>  *scEEEta;
   vector<double>  *scEEPhi;
   vector<double>  *scEERawEnergy;
   vector<double>  *scEEpreshowerEnergy;
   vector<double>  *scEEenergy1;
   vector<double>  *scEEenergy2x2;
   vector<double>  *scEEenergy3x3;
   vector<double>  *scEEenergy1x5;
   vector<double>  *scEEenergy2x5;
   vector<double>  *scEEenergy5x5;
   vector<double>  *scEEHoverE;
   vector<double>  *scEEx;
   vector<double>  *scEEy;
   vector<double>  *scEEz;
   vector<double>  *scEEetaWidth;
   vector<double>  *scEEphiWidth;
   vector<double>  *scEEsigmaetaeta;
   vector<double>  *scEEsigmaIetaIeta;
   vector<double>  *scEEsigmaphiphi;
   vector<double>  *scEEsigmaIphiIphi;
   vector<double>  *scEEsize;
   vector<double>  *scEEnBasicClusters;
   vector<bool>    *L1trigResults;
   vector<bool>    *L1trigErrCode;
   vector<string>  *L1trigName;

   // List of branches
   TBranch        *b_EvtInfo_EventNum;   //!
   TBranch        *b_EvtInfo_RunNum;   //!
   TBranch        *b_EvtInfo_LumiSection;   //!
   TBranch        *b_EvtInfo_BunchXing;   //!
   TBranch        *b_trigResults;   //!
   TBranch        *b_trigName;   //!
   TBranch        *b_scEBNum;   //!
   TBranch        *b_scEBEnergy;   //!
   TBranch        *b_scEBEt;   //!
   TBranch        *b_scEBEta;   //!
   TBranch        *b_scEBPhi;   //!
   TBranch        *b_scEBRawEnergy;   //!
   TBranch        *b_scEBpreshowerEnergy;   //!
   TBranch        *b_scEBenergy1;   //!
   TBranch        *b_scEBenergy2x2;   //!
   TBranch        *b_scEBenergy3x3;   //!
   TBranch        *b_scEBenergy1x5;   //!
   TBranch        *b_scEBenergy2x5;   //!
   TBranch        *b_scEBenergy5x5;   //!
   TBranch        *b_scEBHoverE;   //!
   TBranch        *b_scEBx;   //!
   TBranch        *b_scEBy;   //!
   TBranch        *b_scEBz;   //!
   TBranch        *b_scEBetaWidth;   //!
   TBranch        *b_scEBphiWidth;   //!
   TBranch        *b_scEBsigmaetaeta;   //!
   TBranch        *b_scEBsigmaIetaIeta;   //!
   TBranch        *b_scEBsigmaphiphi;   //!
   TBranch        *b_scEBsigmaIphiIphi;   //!
   TBranch        *b_scEBsize;   //!
   TBranch        *b_scEBnBasicClusters;   //!
   TBranch        *b_scEENum;   //!
   TBranch        *b_scEEEnergy;   //!
   TBranch        *b_scEEEt;   //!
   TBranch        *b_scEEEta;   //!
   TBranch        *b_scEEPhi;   //!
   TBranch        *b_scEERawEnergy;   //!
   TBranch        *b_scEEpreshowerEnergy;   //!
   TBranch        *b_scEEenergy1;   //!
   TBranch        *b_scEEenergy2x2;   //!
   TBranch        *b_scEEenergy3x3;   //!
   TBranch        *b_scEEenergy1x5;   //!
   TBranch        *b_scEEenergy2x5;   //!
   TBranch        *b_scEEenergy5x5;   //!
   TBranch        *b_scEEHoverE;   //!
   TBranch        *b_scEEx;   //!
   TBranch        *b_scEEy;   //!
   TBranch        *b_scEEz;   //!
   TBranch        *b_scEEetaWidth;   //!
   TBranch        *b_scEEphiWidth;   //!
   TBranch        *b_scEEsigmaetaeta;   //!
   TBranch        *b_scEEsigmaIetaIeta;   //!
   TBranch        *b_scEEsigmaphiphi;   //!
   TBranch        *b_scEEsigmaIphiIphi;   //!
   TBranch        *b_scEEsize;   //!
   TBranch        *b_scEEnBasicClusters;   //!
   TBranch        *b_L1trigResults;   //!
   TBranch        *b_L1trigErrCode;   //!
   TBranch        *b_L1trigName;   //!

   trigEffi(TTree *tree=0);
   virtual ~trigEffi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef trigEffi_cxx
trigEffi::trigEffi(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = TFile::Open("rfio:/castor/cern.ch/user/l/lovedeep/TRIG_Spring10/mrgd_TrigTree_data_MB_GOODCOLLL_v9.root");
     //      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TrigTree_data.root");
     //      if (!f) {
     //  f = new TFile("TrigTree_data.root");
     f->cd("demo");
     //         f->cd("TrigTree_data.root:/demo");
     //      }
      tree = (TTree*)gDirectory->Get("TrigTree");

   }
   Init(tree);
}

trigEffi::~trigEffi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t trigEffi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t trigEffi::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void trigEffi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trigResults = 0;
   trigName = 0;
   scEBEnergy = 0;
   scEBEt = 0;
   scEBEta = 0;
   scEBPhi = 0;
   scEBRawEnergy = 0;
   scEBpreshowerEnergy = 0;
   scEBenergy1 = 0;
   scEBenergy2x2 = 0;
   scEBenergy3x3 = 0;
   scEBenergy1x5 = 0;
   scEBenergy2x5 = 0;
   scEBenergy5x5 = 0;
   scEBHoverE = 0;
   scEBx = 0;
   scEBy = 0;
   scEBz = 0;
   scEBetaWidth = 0;
   scEBphiWidth = 0;
   scEBsigmaetaeta = 0;
   scEBsigmaIetaIeta = 0;
   scEBsigmaphiphi = 0;
   scEBsigmaIphiIphi = 0;
   scEBsize = 0;
   scEBnBasicClusters = 0;
   scEEEnergy = 0;
   scEEEt = 0;
   scEEEta = 0;
   scEEPhi = 0;
   scEERawEnergy = 0;
   scEEpreshowerEnergy = 0;
   scEEenergy1 = 0;
   scEEenergy2x2 = 0;
   scEEenergy3x3 = 0;
   scEEenergy1x5 = 0;
   scEEenergy2x5 = 0;
   scEEenergy5x5 = 0;
   scEEHoverE = 0;
   scEEx = 0;
   scEEy = 0;
   scEEz = 0;
   scEEetaWidth = 0;
   scEEphiWidth = 0;
   scEEsigmaetaeta = 0;
   scEEsigmaIetaIeta = 0;
   scEEsigmaphiphi = 0;
   scEEsigmaIphiIphi = 0;
   scEEsize = 0;
   scEEnBasicClusters = 0;
   L1trigResults = 0;
   L1trigErrCode = 0;
   L1trigName = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EvtInfo_EventNum", &EvtInfo_EventNum, &b_EvtInfo_EventNum);
   fChain->SetBranchAddress("EvtInfo_RunNum", &EvtInfo_RunNum, &b_EvtInfo_RunNum);
   fChain->SetBranchAddress("EvtInfo_LumiSection", &EvtInfo_LumiSection, &b_EvtInfo_LumiSection);
   fChain->SetBranchAddress("EvtInfo_BunchXing", &EvtInfo_BunchXing, &b_EvtInfo_BunchXing);
   fChain->SetBranchAddress("trigResults", &trigResults, &b_trigResults);
   fChain->SetBranchAddress("trigName", &trigName, &b_trigName);
   fChain->SetBranchAddress("scEBNum", &scEBNum, &b_scEBNum);
   fChain->SetBranchAddress("scEBEnergy", &scEBEnergy, &b_scEBEnergy);
   fChain->SetBranchAddress("scEBEt", &scEBEt, &b_scEBEt);
   fChain->SetBranchAddress("scEBEta", &scEBEta, &b_scEBEta);
   fChain->SetBranchAddress("scEBPhi", &scEBPhi, &b_scEBPhi);
   fChain->SetBranchAddress("scEBRawEnergy", &scEBRawEnergy, &b_scEBRawEnergy);
   fChain->SetBranchAddress("scEBpreshowerEnergy", &scEBpreshowerEnergy, &b_scEBpreshowerEnergy);
   fChain->SetBranchAddress("scEBenergy1", &scEBenergy1, &b_scEBenergy1);
   fChain->SetBranchAddress("scEBenergy2x2", &scEBenergy2x2, &b_scEBenergy2x2);
   fChain->SetBranchAddress("scEBenergy3x3", &scEBenergy3x3, &b_scEBenergy3x3);
   fChain->SetBranchAddress("scEBenergy1x5", &scEBenergy1x5, &b_scEBenergy1x5);
   fChain->SetBranchAddress("scEBenergy2x5", &scEBenergy2x5, &b_scEBenergy2x5);
   fChain->SetBranchAddress("scEBenergy5x5", &scEBenergy5x5, &b_scEBenergy5x5);
   fChain->SetBranchAddress("scEBHoverE", &scEBHoverE, &b_scEBHoverE);
   fChain->SetBranchAddress("scEBx", &scEBx, &b_scEBx);
   fChain->SetBranchAddress("scEBy", &scEBy, &b_scEBy);
   fChain->SetBranchAddress("scEBz", &scEBz, &b_scEBz);
   fChain->SetBranchAddress("scEBetaWidth", &scEBetaWidth, &b_scEBetaWidth);
   fChain->SetBranchAddress("scEBphiWidth", &scEBphiWidth, &b_scEBphiWidth);
   fChain->SetBranchAddress("scEBsigmaetaeta", &scEBsigmaetaeta, &b_scEBsigmaetaeta);
   fChain->SetBranchAddress("scEBsigmaIetaIeta", &scEBsigmaIetaIeta, &b_scEBsigmaIetaIeta);
   fChain->SetBranchAddress("scEBsigmaphiphi", &scEBsigmaphiphi, &b_scEBsigmaphiphi);
   fChain->SetBranchAddress("scEBsigmaIphiIphi", &scEBsigmaIphiIphi, &b_scEBsigmaIphiIphi);
   fChain->SetBranchAddress("scEBsize", &scEBsize, &b_scEBsize);
   fChain->SetBranchAddress("scEBnBasicClusters", &scEBnBasicClusters, &b_scEBnBasicClusters);
   fChain->SetBranchAddress("scEENum", &scEENum, &b_scEENum);
   fChain->SetBranchAddress("scEEEnergy", &scEEEnergy, &b_scEEEnergy);
   fChain->SetBranchAddress("scEEEt", &scEEEt, &b_scEEEt);
   fChain->SetBranchAddress("scEEEta", &scEEEta, &b_scEEEta);
   fChain->SetBranchAddress("scEEPhi", &scEEPhi, &b_scEEPhi);
   fChain->SetBranchAddress("scEERawEnergy", &scEERawEnergy, &b_scEERawEnergy);
   fChain->SetBranchAddress("scEEpreshowerEnergy", &scEEpreshowerEnergy, &b_scEEpreshowerEnergy);
   fChain->SetBranchAddress("scEEenergy1", &scEEenergy1, &b_scEEenergy1);
   fChain->SetBranchAddress("scEEenergy2x2", &scEEenergy2x2, &b_scEEenergy2x2);
   fChain->SetBranchAddress("scEEenergy3x3", &scEEenergy3x3, &b_scEEenergy3x3);
   fChain->SetBranchAddress("scEEenergy1x5", &scEEenergy1x5, &b_scEEenergy1x5);
   fChain->SetBranchAddress("scEEenergy2x5", &scEEenergy2x5, &b_scEEenergy2x5);
   fChain->SetBranchAddress("scEEenergy5x5", &scEEenergy5x5, &b_scEEenergy5x5);
   fChain->SetBranchAddress("scEEHoverE", &scEEHoverE, &b_scEEHoverE);
   fChain->SetBranchAddress("scEEx", &scEEx, &b_scEEx);
   fChain->SetBranchAddress("scEEy", &scEEy, &b_scEEy);
   fChain->SetBranchAddress("scEEz", &scEEz, &b_scEEz);
   fChain->SetBranchAddress("scEEetaWidth", &scEEetaWidth, &b_scEEetaWidth);
   fChain->SetBranchAddress("scEEphiWidth", &scEEphiWidth, &b_scEEphiWidth);
   fChain->SetBranchAddress("scEEsigmaetaeta", &scEEsigmaetaeta, &b_scEEsigmaetaeta);
   fChain->SetBranchAddress("scEEsigmaIetaIeta", &scEEsigmaIetaIeta, &b_scEEsigmaIetaIeta);
   fChain->SetBranchAddress("scEEsigmaphiphi", &scEEsigmaphiphi, &b_scEEsigmaphiphi);
   fChain->SetBranchAddress("scEEsigmaIphiIphi", &scEEsigmaIphiIphi, &b_scEEsigmaIphiIphi);
   fChain->SetBranchAddress("scEEsize", &scEEsize, &b_scEEsize);
   fChain->SetBranchAddress("scEEnBasicClusters", &scEEnBasicClusters, &b_scEEnBasicClusters);
   fChain->SetBranchAddress("L1trigResults", &L1trigResults, &b_L1trigResults);
   fChain->SetBranchAddress("L1trigErrCode", &L1trigErrCode, &b_L1trigErrCode);
   fChain->SetBranchAddress("L1trigName", &L1trigName, &b_L1trigName);
   Notify();
}

Bool_t trigEffi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void trigEffi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t trigEffi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef trigEffi_cxx
