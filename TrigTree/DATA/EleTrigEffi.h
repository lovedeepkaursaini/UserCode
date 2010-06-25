//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jun 13 23:12:03 2010 by ROOT version 5.22/00d
// from TTree TrigTree/TrigTree
// found on file: rfio:///castor/cern.ch/user/l/lovedeep/TRIG_Spring10/mrgd_Run2010A-May27thReReco_v1.root
//////////////////////////////////////////////////////////

#ifndef EleTrigEffi_h
#define EleTrigEffi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "histo.C"
#include <iostream>

using namespace std;

class EleKin{
 public:
  EleKin(double et_, double eta_, double phi_,double HoE_, double Sihih_, double deta_,
	 double dphi_,double trkiso_,double ecaliso_,double hcaliso_,
	 double spikes_,unsigned int mishits_, double dist_, 
	 double dcottheta_,double hltmatch_, bool ecaldriven_, bool trkdriven_):
    et(et_),eta(eta_),phi(phi_),HoE(HoE_),Sihih(Sihih_),deta(deta_),dphi(dphi_),trkiso(trkiso_),ecaliso(ecaliso_),hcaliso(hcaliso_),spikes(spikes_),mishits(mishits_),dist(dist_),dcottheta(dcottheta_),hltmatch(hltmatch_),ecaldriven(ecaldriven_),trkdriven(trkdriven_){}
    ~EleKin(){}
    double Et(){return et;}
    double Eta(){return eta;}
    double Phi(){return phi;}
    double HoverE() {return HoE;}
    double sigmaIetaIeta() {return Sihih;}
    double dEta() {return deta;}
    double dPhi() {return dphi;}
    double TrkIso() {return trkiso;}
    double EcalIso() {return ecaliso;}
    double HcalIso() {return hcaliso;}
    double Spikes() {return spikes;}
    double MisHits() {return mishits;}
    double Dist() {return dist;}
    double dCotTheta() {return dcottheta;}
    double HltMatch() {return hltmatch;}
    double EcalDriven() {return ecaldriven;}
    double TrkDriven() {return trkdriven;}
    double et;
    double eta;
    double phi;
    double HoE;
    double Sihih;
    double deta;
    double dphi;
    double trkiso;
    double ecaliso;
    double hcaliso;
    double spikes;
    unsigned int mishits;
    double dist;
    double dcottheta;
    double hltmatch;
    bool ecaldriven;
    bool trkdriven;
  private:
    EleKin();
};

typedef	std::map<std::string,std::vector<EleKin> > eleColl;

class EleTrigEffi {
public :
  bool useEndCapEl_;
  bool useBrlEl_;
  void UseBarrelEl(){useEndCapEl_=0; useBrlEl_=1;}
  void UseEndCapEl(){useEndCapEl_=1; useBrlEl_=0;}
  std::vector<EleKin>
    SelectElectrons(std::string trigName, std::string cutString, eleWorkPoint* eleWP);
  std::vector<std::string> 
    cutParser(std::string cutString);
  
  double FindHltMatch(std::string trigname, double Eta, double Phi);
  std::map<std::string, bool> trigStuff;


  void Loop(std::vector<std::string> selVec, eleWorkPoint* eleWP1, TFile* file, Long64_t min, Long64_t max);
  eleColl MkCollection(std::string trigName, std::vector<std::string> selecVec,eleWorkPoint* eleWP);

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EvtInfo_EventNum;
   Int_t           EvtInfo_RunNum;
   Int_t           EvtInfo_LumiSection;
   Int_t           EvtInfo_BunchXing;
   Double_t        ElectronNum;
   vector<double>  *ElectronPt;
   vector<double>  *ElectronEta;
   vector<double>  *ElectronPhi;
   vector<double>  *ElectronEt;
   vector<double>  *ElectronscEt;
   vector<double>  *ElectronscEta;
   vector<double>  *Electronscphi;
   vector<double>  *Electronscenergy;
   vector<double>  *Electronelermax3x3;
   vector<double>  *ElectroneleisBarrel;
   vector<double>  *ElectroneleisEndcap;
   vector<double>  *ElectroneleisEcalDriven;
   vector<double>  *ElectroneleisTrackerDriven;
   vector<double>  *Electroneleecalenergy;
   vector<double>  *Electroneleecalenergyerror;
   vector<double>  *Electronelescr9;
   vector<double>  *Electronelesceseedoveresupercluster;
   vector<double>  *ElectronDhElecClsTrkAtCalo;
   vector<double>  *ElectronDhSeedClsTrkAtCalo;
   vector<double>  *ElectronDhSuperClsTrkAtVtx;
   vector<double>  *ElectronDphiElecClsTrkAtCalo;
   vector<double>  *ElectronDphiSeedClsTrkAtCalo;
   vector<double>  *ElectronDphiSuperClsTrkAtVtx;
   vector<double>  *ElectronPositionXTrkAtVtx;
   vector<double>  *ElectronPositionYTrkAtVtx;
   vector<double>  *ElectronPositionZTrkAtVtx;
   vector<double>  *ElectronMomentumXTrkAtVtx;
   vector<double>  *ElectronMomentumYTrkAtVtx;
   vector<double>  *ElectronMomentumZTrkAtVtx;
   vector<double>  *ElectronPositionXTrkAtCalo;
   vector<double>  *ElectronPositionYTrkAtCalo;
   vector<double>  *ElectronPositionZTrkAtCalo;
   vector<double>  *ElectronMomentumXTrkAtCalo;
   vector<double>  *ElectronMomentumYTrkAtCalo;
   vector<double>  *ElectronMomentumZTrkAtCalo;
   vector<double>  *ElectroneEleClsOverPout;
   vector<double>  *ElectroneSeedClsOverP;
   vector<double>  *ElectroneSeedClsOverPout;
   vector<double>  *ElectroneSuperClsOverP;
   vector<double>  *ElectroneleMomErr;
   vector<double>  *Electroneledr03EcalRecHitSumEt;
   vector<double>  *Electroneledr03HcalDepth1TowerSumEt;
   vector<double>  *Electroneledr03HcalDepth2TowerSumEt;
   vector<double>  *Electroneledr03HcalTowerSumEt;
   vector<double>  *Electroneledr03TkSumPt;
   vector<double>  *Electroneledr04EcalRecHitSumEt;
   vector<double>  *Electroneledr04HcalDepth1TowerSumEt;
   vector<double>  *Electroneledr04HcalDepth2TowerSumEt;
   vector<double>  *Electroneledr04HcalTowerSumEt;
   vector<double>  *Electroneledr04TkSumPt;
   vector<double>  *ElectroneleRelIsoEcal;
   vector<double>  *ElectroneleRelIsoHcal;
   vector<double>  *ElectroneleRelIsoTrk;
   vector<double>  *ElectroneleRelIsoComb;
   vector<double>  *ElectroneleMissingHits;
   vector<double>  *ElectroneleDist;
   vector<double>  *ElectroneleDeltaCotTheta;
   vector<double>  *ElectroneleConvRadius;
   vector<double>  *Electrone1x5;
   vector<double>  *Electrone2x5Max;
   vector<double>  *Electrone5x5;
   vector<double>  *Electroneler1x5;
   vector<double>  *Electroneler2x5max;
   vector<double>  *Electronscpreshowerenergy;
   vector<double>  *Electronscetawidth;
   vector<double>  *Electronscphiwidth;
   vector<double>  *Electroneleenergy;
   vector<double>  *ElectronelehcalDepth1OverEcal;
   vector<double>  *ElectronelehcalDepth2OverEcal;
   vector<double>  *ElectronelehcalOverEcal;
   vector<double>  *ElectronelesigmaEtaEta;
   vector<double>  *ElectronelesigmaIetaIeta;
   vector<double>  *ElectronelebasicClustersSize;
   vector<double>  *ElectronelenumberOfBrems;
   vector<double>  *Electronelefbrem;
   vector<double>  *ElectronelescPixCharge;
   vector<double>  *Electronelectfcharge;
   vector<double>  *Electronelecharge;
   vector<double>  *ElectroneleisGsfScPixChargeConsistent;
   vector<double>  *ElectroneleisGsfCtfChargeConsistent;
   vector<double>  *ElectroneleisGsfCtfScPixChargeConsistent;
   vector<double>  *ElectronelevertexChi2;
   vector<double>  *ElectronelevertexNdof;
   vector<double>  *ElectronelevertexNormalizedChi2;
   vector<double>  *Electronelevx;
   vector<double>  *Electronelevy;
   vector<double>  *Electronelevz;
   vector<double>  *ElectronelevertexX;
   vector<double>  *ElectronelevertexY;
   vector<double>  *ElectronelevertexZ;
   vector<double>  *ElectronelevertexTIP;
   vector<double>  *Electronelep;
   vector<double>  *Electronelepx;
   vector<double>  *Electronelepy;
   vector<double>  *Electronelepz;
   vector<double>  *Electroneledxy;
   vector<double>  *Electronelegsfcharge;
   vector<double>  *ElectroneleambiguousTracks;
   vector<double>  *ElectronelefoundHits;
   vector<double>  *ElectronelelostHits;
   vector<double>  *Electronelechi2;
   vector<int>    *trigResults;
   vector<string>  *trigName;
   vector<double>  *hltEle10LWPt;
   vector<double>  *hltEle10LWEta;
   vector<double>  *hltEle10LWPhi;
   vector<double>  *hltEle15LWPt;
   vector<double>  *hltEle15LWEta;
   vector<double>  *hltEle15LWPhi;
   vector<double>  *hltPhoton10Pt;
   vector<double>  *hltPhoton10Eta;
   vector<double>  *hltPhoton10Phi;
   vector<double>  *hltPhoton10now15Pt;
   vector<double>  *hltPhoton10now15Eta;
   vector<double>  *hltPhoton10now15Phi;
   vector<double>  *hltPhoton15Pt;
   vector<double>  *hltPhoton15Eta;
   vector<double>  *hltPhoton15Phi;
   vector<int>     *L1trigResults;
   vector<int>     *L1trigErrCode;
   vector<string>  *L1trigName;
   vector<double>  *l1IsoEleEt;
   vector<double>  *l1IsoEleEnergy;
   vector<double>  *l1IsoEleEta;
   vector<double>  *l1IsoElePhi;
   vector<double>  *l1NonIsoEleEt;
   vector<double>  *l1NonIsoEleEnergy;
   vector<double>  *l1NonIsoEleEta;
   vector<double>  *l1NonIsoElePhi;

   // List of branches
   TBranch        *b_EvtInfo_EventNum;   //!
   TBranch        *b_EvtInfo_RunNum;   //!
   TBranch        *b_EvtInfo_LumiSection;   //!
   TBranch        *b_EvtInfo_BunchXing;   //!
   TBranch        *b_ElectronNum;   //!
   TBranch        *b_ElectronPt;   //!
   TBranch        *b_ElectronEta;   //!
   TBranch        *b_ElectronPhi;   //!
   TBranch        *b_ElectronEt;   //!
   TBranch        *b_ElectronscEt;   //!
   TBranch        *b_ElectronscEta;   //!
   TBranch        *b_Electronscphi;   //!
   TBranch        *b_Electronscenergy;   //!
   TBranch        *b_Electronelermax3x3;   //!
   TBranch        *b_ElectroneleisBarrel;   //!
   TBranch        *b_ElectroneleisEndcap;   //!
   TBranch        *b_ElectroneleisEcalDriven;   //!
   TBranch        *b_ElectroneleisTrackerDriven;   //!
   TBranch        *b_Electroneleecalenergy;   //!
   TBranch        *b_Electroneleecalenergyerror;   //!
   TBranch        *b_Electronelescr9;   //!
   TBranch        *b_Electronelesceseedoveresupercluster;   //!
   TBranch        *b_ElectronDhElecClsTrkAtCalo;   //!
   TBranch        *b_ElectronDhSeedClsTrkAtCalo;   //!
   TBranch        *b_ElectronDhSuperClsTrkAtVtx;   //!
   TBranch        *b_ElectronDphiElecClsTrkAtCalo;   //!
   TBranch        *b_ElectronDphiSeedClsTrkAtCalo;   //!
   TBranch        *b_ElectronDphiSuperClsTrkAtVtx;   //!
   TBranch        *b_ElectronPositionXTrkAtVtx;   //!
   TBranch        *b_ElectronPositionYTrkAtVtx;   //!
   TBranch        *b_ElectronPositionZTrkAtVtx;   //!
   TBranch        *b_ElectronMomentumXTrkAtVtx;   //!
   TBranch        *b_ElectronMomentumYTrkAtVtx;   //!
   TBranch        *b_ElectronMomentumZTrkAtVtx;   //!
   TBranch        *b_ElectronPositionXTrkAtCalo;   //!
   TBranch        *b_ElectronPositionYTrkAtCalo;   //!
   TBranch        *b_ElectronPositionZTrkAtCalo;   //!
   TBranch        *b_ElectronMomentumXTrkAtCalo;   //!
   TBranch        *b_ElectronMomentumYTrkAtCalo;   //!
   TBranch        *b_ElectronMomentumZTrkAtCalo;   //!
   TBranch        *b_ElectroneEleClsOverPout;   //!
   TBranch        *b_ElectroneSeedClsOverP;   //!
   TBranch        *b_ElectroneSeedClsOverPout;   //!
   TBranch        *b_ElectroneSuperClsOverP;   //!
   TBranch        *b_ElectroneleMomErr;   //!
   TBranch        *b_Electroneledr03EcalRecHitSumEt;   //!
   TBranch        *b_Electroneledr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electroneledr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_Electroneledr03HcalTowerSumEt;   //!
   TBranch        *b_Electroneledr03TkSumPt;   //!
   TBranch        *b_Electroneledr04EcalRecHitSumEt;   //!
   TBranch        *b_Electroneledr04HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electroneledr04HcalDepth2TowerSumEt;   //!
   TBranch        *b_Electroneledr04HcalTowerSumEt;   //!
   TBranch        *b_Electroneledr04TkSumPt;   //!
   TBranch        *b_ElectroneleRelIsoEcal;   //!
   TBranch        *b_ElectroneleRelIsoHcal;   //!
   TBranch        *b_ElectroneleRelIsoTrk;   //!
   TBranch        *b_ElectroneleRelIsoComb;   //!
   TBranch        *b_ElectroneleMissingHits;   //!
   TBranch        *b_ElectroneleDist;   //!
   TBranch        *b_ElectroneleDeltaCotTheta;   //!
   TBranch        *b_ElectroneleConvRadius;   //!
   TBranch        *b_Electrone1x5;   //!
   TBranch        *b_Electrone2x5Max;   //!
   TBranch        *b_Electrone5x5;   //!
   TBranch        *b_Electroneler1x5;   //!
   TBranch        *b_Electroneler2x5max;   //!
   TBranch        *b_Electronscpreshowerenergy;   //!
   TBranch        *b_Electronscetawidth;   //!
   TBranch        *b_Electronscphiwidth;   //!
   TBranch        *b_Electroneleenergy;   //!
   TBranch        *b_ElectronelehcalDepth1OverEcal;   //!
   TBranch        *b_ElectronelehcalDepth2OverEcal;   //!
   TBranch        *b_ElectronelehcalOverEcal;   //!
   TBranch        *b_ElectronelesigmaEtaEta;   //!
   TBranch        *b_ElectronelesigmaIetaIeta;   //!
   TBranch        *b_ElectronelebasicClustersSize;   //!
   TBranch        *b_ElectronelenumberOfBrems;   //!
   TBranch        *b_Electronelefbrem;   //!
   TBranch        *b_ElectronelescPixCharge;   //!
   TBranch        *b_Electronelectfcharge;   //!
   TBranch        *b_Electronelecharge;   //!
   TBranch        *b_ElectroneleisGsfScPixChargeConsistent;   //!
   TBranch        *b_ElectroneleisGsfCtfChargeConsistent;   //!
   TBranch        *b_ElectroneleisGsfCtfScPixChargeConsistent;   //!
   TBranch        *b_ElectronelevertexChi2;   //!
   TBranch        *b_ElectronelevertexNdof;   //!
   TBranch        *b_ElectronelevertexNormalizedChi2;   //!
   TBranch        *b_Electronelevx;   //!
   TBranch        *b_Electronelevy;   //!
   TBranch        *b_Electronelevz;   //!
   TBranch        *b_ElectronelevertexX;   //!
   TBranch        *b_ElectronelevertexY;   //!
   TBranch        *b_ElectronelevertexZ;   //!
   TBranch        *b_ElectronelevertexTIP;   //!
   TBranch        *b_Electronelep;   //!
   TBranch        *b_Electronelepx;   //!
   TBranch        *b_Electronelepy;   //!
   TBranch        *b_Electronelepz;   //!
   TBranch        *b_Electroneledxy;   //!
   TBranch        *b_Electronelegsfcharge;   //!
   TBranch        *b_ElectroneleambiguousTracks;   //!
   TBranch        *b_ElectronelefoundHits;   //!
   TBranch        *b_ElectronelelostHits;   //!
   TBranch        *b_Electronelechi2;   //!
   TBranch        *b_trigResults;   //!
   TBranch        *b_trigName;   //!
   TBranch        *b_hltEle10LWPt;   //!
   TBranch        *b_hltEle10LWEta;   //!
   TBranch        *b_hltEle10LWPhi;   //!
   TBranch        *b_hltEle15LWPt;   //!
   TBranch        *b_hltEle15LWEta;   //!
   TBranch        *b_hltEle15LWPhi;   //!
   TBranch        *b_hltPhoton10Pt;   //!
   TBranch        *b_hltPhoton10Eta;   //!
   TBranch        *b_hltPhoton10Phi;   //!
   TBranch        *b_hltPhoton10now15Pt;   //!
   TBranch        *b_hltPhoton10now15Eta;   //!
   TBranch        *b_hltPhoton10now15Phi;   //!
   TBranch        *b_hltPhoton15Pt;   //!
   TBranch        *b_hltPhoton15Eta;   //!
   TBranch        *b_hltPhoton15Phi;   //!
   TBranch        *b_L1trigResults;   //!
   TBranch        *b_L1trigErrCode;   //!
   TBranch        *b_L1trigName;   //!
   TBranch        *b_l1IsoEleEt;   //!
   TBranch        *b_l1IsoEleEnergy;   //!
   TBranch        *b_l1IsoEleEta;   //!
   TBranch        *b_l1IsoElePhi;   //!
   TBranch        *b_l1NonIsoEleEt;   //!
   TBranch        *b_l1NonIsoEleEnergy;   //!
   TBranch        *b_l1NonIsoEleEta;   //!
   TBranch        *b_l1NonIsoElePhi;   //!

   EleTrigEffi(TTree *tree=0);
   virtual ~EleTrigEffi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Long64_t min, Long64_t max,std::string input);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EleTrigEffi_cxx
EleTrigEffi::EleTrigEffi(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = TFile::Open("rfio:///castor/cern.ch/user/l/lovedeep/TRIG_Spring10/mrgd_MB_Commissioning10-GOODCOLL-May27thSkim_v5_AND_33_Run2010A-May27thReReco_v1_noPhyDec.root");//mrgd_Wenu_Spring10-START3X_V26_S09-v1.root");
     f->cd("demo");
     tree = (TTree*)gDirectory->Get("TrigTree");
   }
   Init(tree);
}

EleTrigEffi::~EleTrigEffi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EleTrigEffi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EleTrigEffi::LoadTree(Long64_t entry)
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

void EleTrigEffi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ElectronPt = 0;
   ElectronEta = 0;
   ElectronPhi = 0;
   ElectronEt = 0;
   ElectronscEt = 0;
   ElectronscEta = 0;
   Electronscphi = 0;
   Electronscenergy = 0;
   Electronelermax3x3 = 0;
   ElectroneleisBarrel = 0;
   ElectroneleisEndcap = 0;
   ElectroneleisEcalDriven = 0;
   ElectroneleisTrackerDriven = 0;
   Electroneleecalenergy = 0;
   Electroneleecalenergyerror = 0;
   Electronelescr9 = 0;
   Electronelesceseedoveresupercluster = 0;
   ElectronDhElecClsTrkAtCalo = 0;
   ElectronDhSeedClsTrkAtCalo = 0;
   ElectronDhSuperClsTrkAtVtx = 0;
   ElectronDphiElecClsTrkAtCalo = 0;
   ElectronDphiSeedClsTrkAtCalo = 0;
   ElectronDphiSuperClsTrkAtVtx = 0;
   ElectronPositionXTrkAtVtx = 0;
   ElectronPositionYTrkAtVtx = 0;
   ElectronPositionZTrkAtVtx = 0;
   ElectronMomentumXTrkAtVtx = 0;
   ElectronMomentumYTrkAtVtx = 0;
   ElectronMomentumZTrkAtVtx = 0;
   ElectronPositionXTrkAtCalo = 0;
   ElectronPositionYTrkAtCalo = 0;
   ElectronPositionZTrkAtCalo = 0;
   ElectronMomentumXTrkAtCalo = 0;
   ElectronMomentumYTrkAtCalo = 0;
   ElectronMomentumZTrkAtCalo = 0;
   ElectroneEleClsOverPout = 0;
   ElectroneSeedClsOverP = 0;
   ElectroneSeedClsOverPout = 0;
   ElectroneSuperClsOverP = 0;
   ElectroneleMomErr = 0;
   Electroneledr03EcalRecHitSumEt = 0;
   Electroneledr03HcalDepth1TowerSumEt = 0;
   Electroneledr03HcalDepth2TowerSumEt = 0;
   Electroneledr03HcalTowerSumEt = 0;
   Electroneledr03TkSumPt = 0;
   Electroneledr04EcalRecHitSumEt = 0;
   Electroneledr04HcalDepth1TowerSumEt = 0;
   Electroneledr04HcalDepth2TowerSumEt = 0;
   Electroneledr04HcalTowerSumEt = 0;
   Electroneledr04TkSumPt = 0;
   ElectroneleRelIsoEcal = 0;
   ElectroneleRelIsoHcal = 0;
   ElectroneleRelIsoTrk = 0;
   ElectroneleRelIsoComb = 0;
   ElectroneleMissingHits = 0;
   ElectroneleDist = 0;
   ElectroneleDeltaCotTheta = 0;
   ElectroneleConvRadius = 0;
   Electrone1x5 = 0;
   Electrone2x5Max = 0;
   Electrone5x5 = 0;
   Electroneler1x5 = 0;
   Electroneler2x5max = 0;
   Electronscpreshowerenergy = 0;
   Electronscetawidth = 0;
   Electronscphiwidth = 0;
   Electroneleenergy = 0;
   ElectronelehcalDepth1OverEcal = 0;
   ElectronelehcalDepth2OverEcal = 0;
   ElectronelehcalOverEcal = 0;
   ElectronelesigmaEtaEta = 0;
   ElectronelesigmaIetaIeta = 0;
   ElectronelebasicClustersSize = 0;
   ElectronelenumberOfBrems = 0;
   Electronelefbrem = 0;
   ElectronelescPixCharge = 0;
   Electronelectfcharge = 0;
   Electronelecharge = 0;
   ElectroneleisGsfScPixChargeConsistent = 0;
   ElectroneleisGsfCtfChargeConsistent = 0;
   ElectroneleisGsfCtfScPixChargeConsistent = 0;
   ElectronelevertexChi2 = 0;
   ElectronelevertexNdof = 0;
   ElectronelevertexNormalizedChi2 = 0;
   Electronelevx = 0;
   Electronelevy = 0;
   Electronelevz = 0;
   ElectronelevertexX = 0;
   ElectronelevertexY = 0;
   ElectronelevertexZ = 0;
   ElectronelevertexTIP = 0;
   Electronelep = 0;
   Electronelepx = 0;
   Electronelepy = 0;
   Electronelepz = 0;
   Electroneledxy = 0;
   Electronelegsfcharge = 0;
   ElectroneleambiguousTracks = 0;
   ElectronelefoundHits = 0;
   ElectronelelostHits = 0;
   Electronelechi2 = 0;
   trigResults = 0;
   trigName = 0;
   hltEle10LWPt = 0;
   hltEle10LWEta = 0;
   hltEle10LWPhi = 0;
   hltEle15LWPt = 0;
   hltEle15LWEta = 0;
   hltEle15LWPhi = 0;
   hltPhoton10Pt = 0;
   hltPhoton10Eta = 0;
   hltPhoton10Phi = 0;
   hltPhoton10now15Pt = 0;
   hltPhoton10now15Eta = 0;
   hltPhoton10now15Phi = 0;
   hltPhoton15Pt = 0;
   hltPhoton15Eta = 0;
   hltPhoton15Phi = 0;
   L1trigResults = 0;
   L1trigErrCode = 0;
   L1trigName = 0;
   l1IsoEleEt = 0;
   l1IsoEleEnergy = 0;
   l1IsoEleEta = 0;
   l1IsoElePhi = 0;
   l1NonIsoEleEt = 0;
   l1NonIsoEleEnergy = 0;
   l1NonIsoEleEta = 0;
   l1NonIsoElePhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EvtInfo_EventNum", &EvtInfo_EventNum, &b_EvtInfo_EventNum);
   fChain->SetBranchAddress("EvtInfo_RunNum", &EvtInfo_RunNum, &b_EvtInfo_RunNum);
   fChain->SetBranchAddress("EvtInfo_LumiSection", &EvtInfo_LumiSection, &b_EvtInfo_LumiSection);
   fChain->SetBranchAddress("EvtInfo_BunchXing", &EvtInfo_BunchXing, &b_EvtInfo_BunchXing);
   fChain->SetBranchAddress("ElectronNum", &ElectronNum, &b_ElectronNum);
   fChain->SetBranchAddress("ElectronPt", &ElectronPt, &b_ElectronPt);
   fChain->SetBranchAddress("ElectronEta", &ElectronEta, &b_ElectronEta);
   fChain->SetBranchAddress("ElectronPhi", &ElectronPhi, &b_ElectronPhi);
   fChain->SetBranchAddress("ElectronEt", &ElectronEt, &b_ElectronEt);
   fChain->SetBranchAddress("ElectronscEt", &ElectronscEt, &b_ElectronscEt);
   fChain->SetBranchAddress("ElectronscEta", &ElectronscEta, &b_ElectronscEta);
   fChain->SetBranchAddress("Electronscphi", &Electronscphi, &b_Electronscphi);
   fChain->SetBranchAddress("Electronscenergy", &Electronscenergy, &b_Electronscenergy);
   fChain->SetBranchAddress("Electronelermax3x3", &Electronelermax3x3, &b_Electronelermax3x3);
   fChain->SetBranchAddress("ElectroneleisBarrel", &ElectroneleisBarrel, &b_ElectroneleisBarrel);
   fChain->SetBranchAddress("ElectroneleisEndcap", &ElectroneleisEndcap, &b_ElectroneleisEndcap);
   fChain->SetBranchAddress("ElectroneleisEcalDriven", &ElectroneleisEcalDriven, &b_ElectroneleisEcalDriven);
   fChain->SetBranchAddress("ElectroneleisTrackerDriven", &ElectroneleisTrackerDriven, &b_ElectroneleisTrackerDriven);
   fChain->SetBranchAddress("Electroneleecalenergy", &Electroneleecalenergy, &b_Electroneleecalenergy);
   fChain->SetBranchAddress("Electroneleecalenergyerror", &Electroneleecalenergyerror, &b_Electroneleecalenergyerror);
   fChain->SetBranchAddress("Electronelescr9", &Electronelescr9, &b_Electronelescr9);
   fChain->SetBranchAddress("Electronelesceseedoveresupercluster", &Electronelesceseedoveresupercluster, &b_Electronelesceseedoveresupercluster);
   fChain->SetBranchAddress("ElectronDhElecClsTrkAtCalo", &ElectronDhElecClsTrkAtCalo, &b_ElectronDhElecClsTrkAtCalo);
   fChain->SetBranchAddress("ElectronDhSeedClsTrkAtCalo", &ElectronDhSeedClsTrkAtCalo, &b_ElectronDhSeedClsTrkAtCalo);
   fChain->SetBranchAddress("ElectronDhSuperClsTrkAtVtx", &ElectronDhSuperClsTrkAtVtx, &b_ElectronDhSuperClsTrkAtVtx);
   fChain->SetBranchAddress("ElectronDphiElecClsTrkAtCalo", &ElectronDphiElecClsTrkAtCalo, &b_ElectronDphiElecClsTrkAtCalo);
   fChain->SetBranchAddress("ElectronDphiSeedClsTrkAtCalo", &ElectronDphiSeedClsTrkAtCalo, &b_ElectronDphiSeedClsTrkAtCalo);
   fChain->SetBranchAddress("ElectronDphiSuperClsTrkAtVtx", &ElectronDphiSuperClsTrkAtVtx, &b_ElectronDphiSuperClsTrkAtVtx);
   fChain->SetBranchAddress("ElectronPositionXTrkAtVtx", &ElectronPositionXTrkAtVtx, &b_ElectronPositionXTrkAtVtx);
   fChain->SetBranchAddress("ElectronPositionYTrkAtVtx", &ElectronPositionYTrkAtVtx, &b_ElectronPositionYTrkAtVtx);
   fChain->SetBranchAddress("ElectronPositionZTrkAtVtx", &ElectronPositionZTrkAtVtx, &b_ElectronPositionZTrkAtVtx);
   fChain->SetBranchAddress("ElectronMomentumXTrkAtVtx", &ElectronMomentumXTrkAtVtx, &b_ElectronMomentumXTrkAtVtx);
   fChain->SetBranchAddress("ElectronMomentumYTrkAtVtx", &ElectronMomentumYTrkAtVtx, &b_ElectronMomentumYTrkAtVtx);
   fChain->SetBranchAddress("ElectronMomentumZTrkAtVtx", &ElectronMomentumZTrkAtVtx, &b_ElectronMomentumZTrkAtVtx);
   fChain->SetBranchAddress("ElectronPositionXTrkAtCalo", &ElectronPositionXTrkAtCalo, &b_ElectronPositionXTrkAtCalo);
   fChain->SetBranchAddress("ElectronPositionYTrkAtCalo", &ElectronPositionYTrkAtCalo, &b_ElectronPositionYTrkAtCalo);
   fChain->SetBranchAddress("ElectronPositionZTrkAtCalo", &ElectronPositionZTrkAtCalo, &b_ElectronPositionZTrkAtCalo);
   fChain->SetBranchAddress("ElectronMomentumXTrkAtCalo", &ElectronMomentumXTrkAtCalo, &b_ElectronMomentumXTrkAtCalo);
   fChain->SetBranchAddress("ElectronMomentumYTrkAtCalo", &ElectronMomentumYTrkAtCalo, &b_ElectronMomentumYTrkAtCalo);
   fChain->SetBranchAddress("ElectronMomentumZTrkAtCalo", &ElectronMomentumZTrkAtCalo, &b_ElectronMomentumZTrkAtCalo);
   fChain->SetBranchAddress("ElectroneEleClsOverPout", &ElectroneEleClsOverPout, &b_ElectroneEleClsOverPout);
   fChain->SetBranchAddress("ElectroneSeedClsOverP", &ElectroneSeedClsOverP, &b_ElectroneSeedClsOverP);
   fChain->SetBranchAddress("ElectroneSeedClsOverPout", &ElectroneSeedClsOverPout, &b_ElectroneSeedClsOverPout);
   fChain->SetBranchAddress("ElectroneSuperClsOverP", &ElectroneSuperClsOverP, &b_ElectroneSuperClsOverP);
   fChain->SetBranchAddress("ElectroneleMomErr", &ElectroneleMomErr, &b_ElectroneleMomErr);
   fChain->SetBranchAddress("Electroneledr03EcalRecHitSumEt", &Electroneledr03EcalRecHitSumEt, &b_Electroneledr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electroneledr03HcalDepth1TowerSumEt", &Electroneledr03HcalDepth1TowerSumEt, &b_Electroneledr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electroneledr03HcalDepth2TowerSumEt", &Electroneledr03HcalDepth2TowerSumEt, &b_Electroneledr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("Electroneledr03HcalTowerSumEt", &Electroneledr03HcalTowerSumEt, &b_Electroneledr03HcalTowerSumEt);
   fChain->SetBranchAddress("Electroneledr03TkSumPt", &Electroneledr03TkSumPt, &b_Electroneledr03TkSumPt);
   fChain->SetBranchAddress("Electroneledr04EcalRecHitSumEt", &Electroneledr04EcalRecHitSumEt, &b_Electroneledr04EcalRecHitSumEt);
   fChain->SetBranchAddress("Electroneledr04HcalDepth1TowerSumEt", &Electroneledr04HcalDepth1TowerSumEt, &b_Electroneledr04HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electroneledr04HcalDepth2TowerSumEt", &Electroneledr04HcalDepth2TowerSumEt, &b_Electroneledr04HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("Electroneledr04HcalTowerSumEt", &Electroneledr04HcalTowerSumEt, &b_Electroneledr04HcalTowerSumEt);
   fChain->SetBranchAddress("Electroneledr04TkSumPt", &Electroneledr04TkSumPt, &b_Electroneledr04TkSumPt);
   fChain->SetBranchAddress("ElectroneleRelIsoEcal", &ElectroneleRelIsoEcal, &b_ElectroneleRelIsoEcal);
   fChain->SetBranchAddress("ElectroneleRelIsoHcal", &ElectroneleRelIsoHcal, &b_ElectroneleRelIsoHcal);
   fChain->SetBranchAddress("ElectroneleRelIsoTrk", &ElectroneleRelIsoTrk, &b_ElectroneleRelIsoTrk);
   fChain->SetBranchAddress("ElectroneleRelIsoComb", &ElectroneleRelIsoComb, &b_ElectroneleRelIsoComb);
   fChain->SetBranchAddress("ElectroneleMissingHits", &ElectroneleMissingHits, &b_ElectroneleMissingHits);
   fChain->SetBranchAddress("ElectroneleDist", &ElectroneleDist, &b_ElectroneleDist);
   fChain->SetBranchAddress("ElectroneleDeltaCotTheta", &ElectroneleDeltaCotTheta, &b_ElectroneleDeltaCotTheta);
   fChain->SetBranchAddress("ElectroneleConvRadius", &ElectroneleConvRadius, &b_ElectroneleConvRadius);
   fChain->SetBranchAddress("Electrone1x5", &Electrone1x5, &b_Electrone1x5);
   fChain->SetBranchAddress("Electrone2x5Max", &Electrone2x5Max, &b_Electrone2x5Max);
   fChain->SetBranchAddress("Electrone5x5", &Electrone5x5, &b_Electrone5x5);
   fChain->SetBranchAddress("Electroneler1x5", &Electroneler1x5, &b_Electroneler1x5);
   fChain->SetBranchAddress("Electroneler2x5max", &Electroneler2x5max, &b_Electroneler2x5max);
   fChain->SetBranchAddress("Electronscpreshowerenergy", &Electronscpreshowerenergy, &b_Electronscpreshowerenergy);
   fChain->SetBranchAddress("Electronscetawidth", &Electronscetawidth, &b_Electronscetawidth);
   fChain->SetBranchAddress("Electronscphiwidth", &Electronscphiwidth, &b_Electronscphiwidth);
   fChain->SetBranchAddress("Electroneleenergy", &Electroneleenergy, &b_Electroneleenergy);
   fChain->SetBranchAddress("ElectronelehcalDepth1OverEcal", &ElectronelehcalDepth1OverEcal, &b_ElectronelehcalDepth1OverEcal);
   fChain->SetBranchAddress("ElectronelehcalDepth2OverEcal", &ElectronelehcalDepth2OverEcal, &b_ElectronelehcalDepth2OverEcal);
   fChain->SetBranchAddress("ElectronelehcalOverEcal", &ElectronelehcalOverEcal, &b_ElectronelehcalOverEcal);
   fChain->SetBranchAddress("ElectronelesigmaEtaEta", &ElectronelesigmaEtaEta, &b_ElectronelesigmaEtaEta);
   fChain->SetBranchAddress("ElectronelesigmaIetaIeta", &ElectronelesigmaIetaIeta, &b_ElectronelesigmaIetaIeta);
   fChain->SetBranchAddress("ElectronelebasicClustersSize", &ElectronelebasicClustersSize, &b_ElectronelebasicClustersSize);
   fChain->SetBranchAddress("ElectronelenumberOfBrems", &ElectronelenumberOfBrems, &b_ElectronelenumberOfBrems);
   fChain->SetBranchAddress("Electronelefbrem", &Electronelefbrem, &b_Electronelefbrem);
   fChain->SetBranchAddress("ElectronelescPixCharge", &ElectronelescPixCharge, &b_ElectronelescPixCharge);
   fChain->SetBranchAddress("Electronelectfcharge", &Electronelectfcharge, &b_Electronelectfcharge);
   fChain->SetBranchAddress("Electronelecharge", &Electronelecharge, &b_Electronelecharge);
   fChain->SetBranchAddress("ElectroneleisGsfScPixChargeConsistent", &ElectroneleisGsfScPixChargeConsistent, &b_ElectroneleisGsfScPixChargeConsistent);
   fChain->SetBranchAddress("ElectroneleisGsfCtfChargeConsistent", &ElectroneleisGsfCtfChargeConsistent, &b_ElectroneleisGsfCtfChargeConsistent);
   fChain->SetBranchAddress("ElectroneleisGsfCtfScPixChargeConsistent", &ElectroneleisGsfCtfScPixChargeConsistent, &b_ElectroneleisGsfCtfScPixChargeConsistent);
   fChain->SetBranchAddress("ElectronelevertexChi2", &ElectronelevertexChi2, &b_ElectronelevertexChi2);
   fChain->SetBranchAddress("ElectronelevertexNdof", &ElectronelevertexNdof, &b_ElectronelevertexNdof);
   fChain->SetBranchAddress("ElectronelevertexNormalizedChi2", &ElectronelevertexNormalizedChi2, &b_ElectronelevertexNormalizedChi2);
   fChain->SetBranchAddress("Electronelevx", &Electronelevx, &b_Electronelevx);
   fChain->SetBranchAddress("Electronelevy", &Electronelevy, &b_Electronelevy);
   fChain->SetBranchAddress("Electronelevz", &Electronelevz, &b_Electronelevz);
   fChain->SetBranchAddress("ElectronelevertexX", &ElectronelevertexX, &b_ElectronelevertexX);
   fChain->SetBranchAddress("ElectronelevertexY", &ElectronelevertexY, &b_ElectronelevertexY);
   fChain->SetBranchAddress("ElectronelevertexZ", &ElectronelevertexZ, &b_ElectronelevertexZ);
   fChain->SetBranchAddress("ElectronelevertexTIP", &ElectronelevertexTIP, &b_ElectronelevertexTIP);
   fChain->SetBranchAddress("Electronelep", &Electronelep, &b_Electronelep);
   fChain->SetBranchAddress("Electronelepx", &Electronelepx, &b_Electronelepx);
   fChain->SetBranchAddress("Electronelepy", &Electronelepy, &b_Electronelepy);
   fChain->SetBranchAddress("Electronelepz", &Electronelepz, &b_Electronelepz);
   fChain->SetBranchAddress("Electroneledxy", &Electroneledxy, &b_Electroneledxy);
   fChain->SetBranchAddress("Electronelegsfcharge", &Electronelegsfcharge, &b_Electronelegsfcharge);
   fChain->SetBranchAddress("ElectroneleambiguousTracks", &ElectroneleambiguousTracks, &b_ElectroneleambiguousTracks);
   fChain->SetBranchAddress("ElectronelefoundHits", &ElectronelefoundHits, &b_ElectronelefoundHits);
   fChain->SetBranchAddress("ElectronelelostHits", &ElectronelelostHits, &b_ElectronelelostHits);
   fChain->SetBranchAddress("Electronelechi2", &Electronelechi2, &b_Electronelechi2);
   fChain->SetBranchAddress("trigResults", &trigResults, &b_trigResults);
   fChain->SetBranchAddress("trigName", &trigName, &b_trigName);
   fChain->SetBranchAddress("hltEle10LWPt", &hltEle10LWPt, &b_hltEle10LWPt);
   fChain->SetBranchAddress("hltEle10LWEta", &hltEle10LWEta, &b_hltEle10LWEta);
   fChain->SetBranchAddress("hltEle10LWPhi", &hltEle10LWPhi, &b_hltEle10LWPhi);
   fChain->SetBranchAddress("hltEle15LWPt", &hltEle15LWPt, &b_hltEle15LWPt);
   fChain->SetBranchAddress("hltEle15LWEta", &hltEle15LWEta, &b_hltEle15LWEta);
   fChain->SetBranchAddress("hltEle15LWPhi", &hltEle15LWPhi, &b_hltEle15LWPhi);
   fChain->SetBranchAddress("hltPhoton10Pt", &hltPhoton10Pt, &b_hltPhoton10Pt);
   fChain->SetBranchAddress("hltPhoton10Eta", &hltPhoton10Eta, &b_hltPhoton10Eta);
   fChain->SetBranchAddress("hltPhoton10Phi", &hltPhoton10Phi, &b_hltPhoton10Phi);
   fChain->SetBranchAddress("hltPhoton10now15Pt", &hltPhoton10now15Pt, &b_hltPhoton10now15Pt);
   fChain->SetBranchAddress("hltPhoton10now15Eta", &hltPhoton10now15Eta, &b_hltPhoton10now15Eta);
   fChain->SetBranchAddress("hltPhoton10now15Phi", &hltPhoton10now15Phi, &b_hltPhoton10now15Phi);
   fChain->SetBranchAddress("hltPhoton15Pt", &hltPhoton15Pt, &b_hltPhoton15Pt);
   fChain->SetBranchAddress("hltPhoton15Eta", &hltPhoton15Eta, &b_hltPhoton15Eta);
   fChain->SetBranchAddress("hltPhoton15Phi", &hltPhoton15Phi, &b_hltPhoton15Phi);
   fChain->SetBranchAddress("L1trigResults", &L1trigResults, &b_L1trigResults);
   fChain->SetBranchAddress("L1trigErrCode", &L1trigErrCode, &b_L1trigErrCode);
   fChain->SetBranchAddress("L1trigName", &L1trigName, &b_L1trigName);
   fChain->SetBranchAddress("l1IsoEleEt", &l1IsoEleEt, &b_l1IsoEleEt);
   fChain->SetBranchAddress("l1IsoEleEnergy", &l1IsoEleEnergy, &b_l1IsoEleEnergy);
   fChain->SetBranchAddress("l1IsoEleEta", &l1IsoEleEta, &b_l1IsoEleEta);
   fChain->SetBranchAddress("l1IsoElePhi", &l1IsoElePhi, &b_l1IsoElePhi);
   fChain->SetBranchAddress("l1NonIsoEleEt", &l1NonIsoEleEt, &b_l1NonIsoEleEt);
   fChain->SetBranchAddress("l1NonIsoEleEnergy", &l1NonIsoEleEnergy, &b_l1NonIsoEleEnergy);
   fChain->SetBranchAddress("l1NonIsoEleEta", &l1NonIsoEleEta, &b_l1NonIsoEleEta);
   fChain->SetBranchAddress("l1NonIsoElePhi", &l1NonIsoElePhi, &b_l1NonIsoElePhi);
   Notify();
}

Bool_t EleTrigEffi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EleTrigEffi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EleTrigEffi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EleTrigEffi_cxx
