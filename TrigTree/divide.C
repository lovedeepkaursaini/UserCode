TCanvas* makeCan(TString h1,TString h2)
{
  TFile *file1=new TFile("EffiHisto_data.root");
  TCanvas* c1=new TCanvas("c1","c1");

  TH1D *hist1=(TH1D*)file1->Get(h1);
  TH1D *hist2=(TH1D*)file1->Get(h2);

  TFile *file2=new TFile("EffiHisto_mc.root");
  TH1D *hist11=(TH1D*)file2->Get(h1);
  TH1D *hist22=(TH1D*)file2->Get(h2);

  //STACK
  THStack * s_1 = new THStack( "s_1", hist2->GetTitle() );
  hist2->SetMarkerColor(1);
  hist2->SetMarkerStyle(8);
  hist2->SetMarkerSize(0.8);
  hist2->Sumw2();
  hist1->Sumw2();
  TGraphAsymmErrors * gr = new TGraphAsymmErrors(); 
  gr->BayesDivide(hist2, hist1,"w");
  gr->SetMarkerColor(1);
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(0.8);

  hist22->SetLineColor(4);
  hist22->SetLineWidth(2);
  hist22->Divide(hist11);
  s_1->Add(hist22);
  //  s_1->GetXaxis()->SetAxisTitle(hist1->GetXaxis()->GetAxisTitle);
  s_1->Draw("nostack");
  TAxis * histXAxis = (TAxis*) hist1->GetXaxis();
  TAxis * stackXAxis = (TAxis*) s_1->GetXaxis();
  stackXAxis->SetTitle(histXAxis->GetTitle());

  TLegend *leg = new TLegend(0.6471774,0.7860169,0.9717742,0.9957627,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(52);
  leg->SetTextSize(0.0545977);
  leg->AddEntry(hist2,"DATA","l");
  leg->AddEntry(hist22,"MC","l");
  leg->Draw();

  s_1->Draw("nostack");
  gr->Draw("zp");
  c1->Modified();
  c1->cd();
  c1->SaveAs("SAVE/eff_"+h2+".gif");
  return c1;
}

void divide()
{
  TCanvas* cc1;
  cc1=makeCan("hEtSCAll","hEtSCTrg");
  cc1=makeCan("hEtSCAllEt10","hEtSCTrgEt10");
  cc1=makeCan("hEtSCAllEt10HoE","hEtSCTrgEt10HoE");
  cc1=makeCan("hEtSCAllEt10HoESihih","hEtSCTrgEt10HoESihih");
  cc1=makeCan("hEtSCAllEt10HoESihihnBC","hEtSCTrgEt10HoESihihnBC");

  cc1=makeCan("hEtSCAll","hEtSCL1Trg");
  cc1=makeCan("hEtSCAllEt10","hEtSCL1TrgEt10");
  cc1=makeCan("hEtSCAllEt10HoE","hEtSCL1TrgEt10HoE");
  cc1=makeCan("hEtSCAllEt10HoESihih","hEtSCL1TrgEt10HoESihih");
  cc1=makeCan("hEtSCAllEt10HoESihihnBC","hEtSCL1TrgEt10HoESihihnBC");

}
