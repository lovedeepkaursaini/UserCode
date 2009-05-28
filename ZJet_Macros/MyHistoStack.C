#include "/home/lovedeep/mystyle.C"
void application()
{
  TCanvas *c1;
  THStack *hs1;
  
  MyHistoStackA(c1,hs1,"Reco/jet_Pt_reco");
  MyHistoStackA(c1,hs1,"Reco/Z_M_reco","m_{e^{+}e^{-}}");
  MyHistoStackA(c1,hs1,"Reco/Z_Pt_reco");
  MyHistoStackA(c1,hs1,"Reco/eM_Pt_reco","p_{T}");
  MyHistoStackA(c1,hs1,"Reco/eM_Eta_reco","#eta");
  MyHistoStackA(c1,hs1,"Reco/eM_Phi_reco","#phi");
  MyHistoStackA(c1,hs1,"Reco/eP_Pt_reco","p_{T}");
  MyHistoStackA(c1,hs1,"Reco/eP_Eta_reco","#eta");
  MyHistoStackA(c1,hs1,"Reco/eP_Phi_reco","#phi");
}


void MyHistoStackA(TCanvas *c1, THStack* hs,std::string histname,std::string XTitle="p_{T}",std::string YTitle="Number of Events")
{
  setTDRStyle();
  tdrStyle->SetPadLeftMargin(0.2);
  tdrStyle->SetPadRightMargin(0.1);
  tdrStyle->SetPadTopMargin(0.08);
  tdrStyle->SetLegendBorderSize(0);
  c1= new TCanvas("c1");  
  // c1->SetLogy();
  TFile ww("WW.root");
  TFile wj("Wjets.root");
  TFile tt("TTbar.root");
  TFile bc("BCtoEQCD.root");
  TFile em("emEnrichQCD.root");
  TFile zj("ZeeJet.root");
  
  TH1F *hwj=(TH1F*)wj.Get(histname.c_str());
  TH1F *htt=(TH1F*)tt.Get(histname.c_str());
  TH1F *hbc=(TH1F*)bc.Get(histname.c_str());
  TH1F *hem=(TH1F*)em.Get(histname.c_str());
  TH1F *hzj=(TH1F*)zj.Get(histname.c_str());
  
  hzj->Sumw2();
  hzj->SetMarkerSize(0.7);
  hzj->SetMarkerStyle(20);
  hzj->SetMarkerColor(4);
  hbc->SetFillColor(kGreen);
  hem->SetFillColor(kCyan);
  hwj->SetFillColor(kYellow);
  htt->SetFillColor(kRed);
  hww->SetFillColor(kMagenta);
  
  hs=new THStack(histname.c_str(),histname.c_str());
  hs->Add(hww);
  hs->Add(htt);
  hs->Add(hwj);
  hs->Add(hbc);
  hs->Add(hem);
  hs->Add(hzj);
  hs->Draw();
  hs->GetYaxis()->SetTitle(YTitle.c_str());
  hs->GetYaxis()->SetTitleOffset(1.6);
  hs->SetMinimum(1);
  hs->GetXaxis()->SetTitle(XTitle.c_str());
  hs->GetXaxis()->SetLabelSize(0.03);
  
  leg = new TLegend(0.7031873,0.7029289,0.8685259,0.8953975,NULL,"brNDC"); //coordinates are fractions
  leg->SetBorderSize(0); 
  leg->SetFillColor(0);
  leg->AddEntry(hzj,"ZeeJet","f"); 
  leg->AddEntry(hem,"emEnriched QCD","f"); 
  leg->AddEntry(hbc,"bc to e","f"); 
  leg->AddEntry(hwj,"Wjets","f");  // "l" means line
  leg->AddEntry(htt,"ttbar","f"); 
  leg->AddEntry(hww,"WW","f");  // "l" means line
  // use "f" for a box
  leg->Draw();
  c1->Draw();
  //histnamep=histname+".png";
  histnamee=histname+".eps";
  histnameg=histname+".gif";
  c1->SaveAs(histnamee.c_str());
  //c1->SaveAs(histnamep.c_str());
  c1->SaveAs(histnameg.c_str());
  
}
