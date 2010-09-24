//**********************//
/*
  USAGE: 
   root.exe -b -l -q testEff*.root 'PlotEffi.cxx("data")'

Followed Muon TnP tools
*/
//**********************//
#include <TCanvas.h>
#include <TPad.h>
#include <string>
#include <iostream>
//**********************//


//**********************//
//User defined variables//
//**********************//
TString prefix = "plots/";
TString basedir = "IdToHLT";
TString what[4] = { "et", "pt", "eta", "pt_eta" };
TString fit[4] = { "probe_sc_et_PLOT","probe_gsfEle_pt_PLOT","probe_sc_eta_PLOT","probe_sc_et_probe_sc_eta_PLOT" };
TString retitle = "Efficiency";
double yMax = 1.1;
double yMin = 0.0;
bool doSquare = true;
TString preliminary = "CMS Preliminary,   #sqrt{s} = 7 TeV";
bool doPdf = true;
//**********************//


//**********************//
TCanvas *c1 = 0;
void PlotEffi(TString scenario) {
    prefix = prefix+scenario+"/";
    gSystem->mkdir(prefix,true);
    gROOT->ProcessLine(".x tdrstyle.cc");
    gStyle->SetOptStat(0);
    c1 = new TCanvas("c1","c1");
    ((TFile*) gROOT->GetListOfFiles()->At(0))->cd();
    plotTriggerData();
}

void plotTriggerData() {
  for (size_t i = 0; i < 4; ++i) {
        TString idname = what[i];
	std::string fitname = fit[i];
	if(fitname.find("_probe_")!= string::npos) 
	  plotTriggerData2D(idname,fitname);
        else plotTriggerData(idname,fitname);
      }
}

void plotTriggerData(TString idname,TString fitname) {
  TDirectory *fit  = gFile->GetDirectory(basedir+"/"+idname+"/");
  single(fit,  idname, fitname );
}

void plotTriggerData2D(TString idname,TString fitname) {
  TDirectory *fit  = gFile->GetDirectory(basedir+"/"+idname+"/");
  single2D(fit,  idname, fitname );
}

TObject *getFromPrefix(TDirectory *dir, TString prefixName) {
    TIter next(dir->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *) next())) {
        if (strstr(key->GetName(), prefixName.Data()) == key->GetName() ) {
            return dir->Get(key->GetName());
        }
    }
    return 0;
}

void EffPalette()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;

   double mid = 0.5*(yMin+yMax);
   Double_t Red[3]    = { 1.00, 1.00, 0.00 };
   Double_t Green[3]  = { 0.00, 1.00, 1.00 };
   Double_t Blue[3]   = { 0.00, 0.00, 0.00 };
   Double_t Length[3] = { 0,    0.75,  1.0  };

   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
   }
}

void reTitleTAxis(TAxis *ax, TString ytitle, double yoffset=1.0) {
   ax->SetTitle(ytitle); 
   ax->SetTitleOffset(yoffset); 
   ax->SetDecimals(true);
   ax->SetMoreLogLabels(true);
}

void reTitleY(TCanvas *pl, TString ytitle) {
    TH1 *first = (TH1*) pl->GetListOfPrimitives()->At(0);
    TH1 *last = (TH1*) pl->GetListOfPrimitives()->At(pl->GetListOfPrimitives()->GetSize()-1);
    if (first) reTitleTAxis(first->GetYaxis(), ytitle);
    if (last)  reTitleTAxis(last->GetYaxis(), ytitle);
} 

void setRangeY(TCanvas *c, double min=0, double max=1.1) {
    for (size_t i = 0, n = c->GetListOfPrimitives()->GetSize(); i < n; ++i) {
        TObject *o = c->GetListOfPrimitives()->At(i);
        if (o->InheritsFrom("TH1")) ((TH1*) o)->GetYaxis()->SetRangeUser(min,max);
    }
}

void squareCanvas(TCanvas *c) {
    c->SetCanvasSize(600,600);
    gStyle->SetPaperSize(20.,20.);
}

void cmsprelim() {
    TPaveText *cmsprel = new TPaveText(doSquare ? 0.40 : .55,.16,.94,.21,"NDC");
    cmsprel->SetTextSize(doSquare ? 0.040 : 0.05);
    cmsprel->SetFillColor(0);
    cmsprel->SetFillStyle(0);
    cmsprel->SetLineStyle(2);
    cmsprel->SetLineColor(0);
    cmsprel->SetTextAlign(12);
    cmsprel->SetTextFont(42);
    cmsprel->AddText(preliminary);
    cmsprel->Draw("same");
}

void single2D( TDirectory *fit, TString alias, TString fitname) {
    gStyle->SetPaintTextFormat(".4f");
    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    TH2 *hfit = pfit->FindObject(pfit->GetName());
    hfit->SetTitle(retitle == "" ?  "efficiency" : retitle);
    hfit->GetZaxis()->SetTitle("");
    hfit->GetZaxis()->SetRangeUser(yMin, yMax);

    EffPalette();
    c1->cd();
    double orm = c1->GetRightMargin();
    c1->SetRightMargin(0.14);
    hfit->Draw();
    if (doSquare) squareCanvas(c1);
    if (preliminary != "") cmsprelim();
    pfit->Print(prefix+alias+".png"); 
    if (doPdf) pfit->Print(prefix+alias+".pdf"); 
    c1->SetRightMargin(orm);
}

void single( TDirectory *fit, TString alias, TString fitname) {
    TCanvas *pfit = getFromPrefix(fit->GetDirectory("fit_eff_plots"), fitname);
    if (pfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << " in " << fit->GetName() << std::endl;
        return;
    }
    if (retitle != "") reTitleY(pfit, retitle);

    RooHist *hfit = (RooHist *) pfit->FindObject("hxy_fit_eff");
    if (hfit == 0) {
        std::cerr << "NOT FOUND: " << "fit_eff_plots/"+fitname << "/hxy_fit_eff in " << fit->GetName() << std::endl;
        pfit->ls();
        return;
    }
    hfit->SetLineWidth(2);
    hfit->SetLineColor(kBlue);
    hfit->SetMarkerColor(kBlue);
    hfit->SetMarkerStyle(20);
    hfit->SetMarkerSize(1.6);
    setRangeY(pfit, yMin, yMax);
    if (doSquare) squareCanvas(pfit);
    if (preliminary != "") cmsprelim();
    pfit->Print(prefix+alias+".png"); 
    if (doPdf) pfit->Print(prefix+alias+".pdf"); 
}
