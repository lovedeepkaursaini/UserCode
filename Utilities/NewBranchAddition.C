#include <vector>
#include <cmath>
#include "fstream.h"
/***********************************************************
   ADDING NEW BRANCHES TO THE EXISTING TREE BRANCHES 
   HOW TO RUN:
   .x NewBranchAddition.C("input.txt","analysis",100)
***********************************************************
Arguments are: 
1) input file name with three colomns: 
       root file name, 
       Xsections of datasets, 
       their filterEffeciency (Twiki page helps here)
2) tree name 
3) the required Luminosity for scaling

Lovedeep, Panjab
March, 2009
***********************************************************/
void NewBranchAddition(TString readfile,TString treeName,Float_t ILuminosity) {
  
  ifstream get;
  get.open(readfile);
  if ( !get.good() ) {
    cout << "ur txt file not good " << readfile << endl;
    return;
  }
  
  vector<TString> NtupleFileName;
  vector<Float_t> Xsection;
  vector<Float_t> FilterEff;
  
  TString currentNtupleFileName;
  Float_t currentXsection;
  Float_t currentFilterEff;
  
  cout<<"Input is:"<<endl;
  cout<<"Root File Name, "<<'\t'<<" Xsection, "<<'\t'<<" Filter Eff. "<<endl;
  //READ YOUR *.TXT FILE AND STORE ITS CONTENTS FOR YOUR USE...
  while ( get >> currentNtupleFileName >> currentXsection>>currentFilterEff) { 
    cout<<currentNtupleFileName<<","<<'\t'<<currentXsection<<","<<'\t'<<currentFilterEff<<endl;
    
    if ( !get.good() ) break;
    NtupleFileName.push_back(currentNtupleFileName);
    Xsection.push_back(currentXsection);
    FilterEff.push_back(currentFilterEff);
  }
  cout<<endl;
  //GET THE NUMBER OF INPUT FILES (pT HAT BINNED DATA SETS)
  const int nFiles=NtupleFileName.size();
  int nEvents[nFiles];  
  //LOOP OVER ALL THE BINNED FILES
  for (int i=0; i <nFiles; i++) {    
    //OPEN THE CURRENT FILE AND STORE THE NUMBER OF EVENTS IN nEvents ARRAY
    TFile currentFile(NtupleFileName[i]);
    TTree *mytree = (TTree*)currentFile.Get(treeName);
    nEvents[i]  = (int) mytree->GetEntries();
    Float_t scale=Xsection[i]*ILuminosity*FilterEff[i]/nEvents[i];
    cout<<NtupleFileName[i]<<" has "<<nEvents[i]
	<<" no. of events and weight for this is "<<scale<<endl;
    //CREATE NEW FILES WITH "OLDNAMES+_newTree.root" AND ALL OTHER
    //OLD TREE BRANCHES AND THE NEW ONE ADDED NOW
    TString newFileName = NtupleFileName[i];
    newFileName.ReplaceAll(".root","_newTree.root");
    TFile *newFile = new TFile(newFileName,"recreate");
    TTree *newTree = mytree->CloneTree();
    //ADDING FOLLOWING NEW BRANCHES
    TBranch *X = newTree->Branch("Xsection", &Xsection[i],"Xsection/F");
    TBranch *FE = newTree->Branch("filterEff", &FilterEff[i],"filterEff/F");
    TBranch *W = newTree->Branch("weight", &scale,"weight/F");
    
    Int_t nEntries = (Int_t)newTree->GetEntries();
    for (Int_t j=0; j<nEntries; j++) {
      X->Fill();
      FE->Fill();
      W->Fill();
    }
    
    newFile->cd();
    newTree->Write();
    currentFile.Close();
  }
}


