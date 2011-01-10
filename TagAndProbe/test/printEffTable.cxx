//**********************//
#include <TCanvas.h>
#include <TPad.h>
#include <string>
#include <iostream>
#include <TFile>
//**********************//
void printEffTable(TString meas) {
  std::string dir[4]={"SCToGsf","GsfToIso","IsoToId","IdToHLT"};
    for(int i=0; i!=4; i++){
      std::string filename = dir[i]+".root";
      TFile *file = new TFile(filename.c_str());
      TString vars="";TString method="fit";
      if(meas=="nJets") vars="probe_"+meas+"+probe_"+meas+"_aerr_lo:probe_"+meas+"+probe_"+meas+"_aerr_hi";
      else if(meas=="pt") vars="probe_gsfEle_"+meas+"+probe_gsfEle_"+meas+"_aerr_lo:probe_gsfEle_"+meas+"+probe_gsfEle_"+meas+"_aerr_hi";
      else if(meas=="et_nJets") vars="probe_sc_et+probe_sc_et_aerr_lo:probe_sc_et+probe_sc_et_aerr_hi:probe_nJets+probe_nJets_aerr_lo:probe_nJets+probe_nJets_aerr_hi";
      else if(meas=="eta_nJets") vars="probe_sc_eta+probe_sc_eta_aerr_lo:probe_sc_eta+probe_sc_eta_aerr_hi:probe_nJets+probe_nJets_aerr_lo:probe_nJets+probe_nJets_aerr_hi";
      else if(meas=="abseta_nJets") vars="probe_sc_abseta+probe_sc_abseta_aerr_lo:probe_sc_abseta+probe_sc_abseta_aerr_hi:probe_nJets+probe_nJets_aerr_lo:probe_nJets+probe_nJets_aerr_hi";
      else if(meas=="et_abseta") vars="probe_sc_abseta+probe_sc_abseta_aerr_lo:probe_sc_abseta+probe_sc_abseta_aerr_hi:probe_sc_et+probe_sc_et_aerr_lo:probe_sc_et+probe_sc_et_aerr_hi";
      else if(meas=="et_eta") vars="probe_sc_eta+probe_sc_eta_aerr_lo:probe_sc_eta+probe_sc_eta_aerr_hi:probe_sc_et+probe_sc_et_aerr_lo:probe_sc_et+probe_sc_et_aerr_hi";
      else vars="probe_sc_"+meas+"+probe_sc_"+meas+"_aerr_lo:probe_sc_"+meas+"+probe_sc_"+meas+"_aerr_hi";
      printEffTable(file, dir[i],meas,vars,method);
      //file->Close();
    }
}

void printEffTable(TFile* file, TString dir, TString meas, TString vars="pt:eta", TString method="fit") {
  RooDataSet *ds = (RooDataSet *) file->Get(dir+"/"+meas+"/"+method+"_eff");
  if (ds == 0) {
    std::cerr << "NOT FOUND: " << (dir+"/"+meas+"/"+method+"_eff") << std::endl;
    if (file->GetDirectory(dir) != 0) {
      file->GetDirectory(dir)->ls();
    } else {
      file->ls();
    }
  }
  ds->tree()->Scan(vars+":efficiency:efficiency_aerr_lo:efficiency_aerr_hi");
}

