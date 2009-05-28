std::ofstream fout("cut_yield_summary.txt");

void application()
{

  //  for(int i=1;i<=50;i++)
  //    {
      //fout<<"=================>eMinus_DeltaPhiIn_zm<================"<<endl;
      //Intg("Reco/eMinus_DeltaPhiIn_zm",-0.001*i,0.001*i);
  //  }

  fout<<"=================>Z_M_reco<================"<<endl<<endl;
  Intg("Reco/Z_M_reco",60.,120.);
  fout<<endl; 
  
  fout<<"=================>After all cuts<================"<<endl<<endl;
  Intg("Reco/Z_M_reco_id",75.,105.);
  fout<<endl;
  
  fout<<"=================End of above set================"<<endl<<endl;
  fout.close();

}


void Intg(std::string histname,float xmin,float xmax)
{
  double count1=0,count2=0;
  
  double integralfull = 0;
  double integral = 0;

  TFile ww("WW.root");
  TFile wj("Wjets.root");
  TFile tt("TTbar.root");
  TFile bc("BCtoEQCD.root");
  TFile em("emEnrichQCD.root");
  TFile zj("ZeeJet.root");

  fout<<"---------ZJet--------"<<endl; 
  TH1F *hzj=(TH1F*)zj.Get(histname.c_str());

  GetInt(hzj,xmin,xmax,integralfull,integral);
  fout<<endl; 
  
  fout<<"---------QCD emEnrich--------"<<endl; 
  TH1F *hem=(TH1F*)em.Get(histname.c_str());
  GetInt(hem,xmin,xmax,integralfull,integral);
  fout<<endl; 
  GetIntg(hem,xmin,xmax,count1, count2);
  
  fout<<"---------QCD BCtoE--------"<<endl; 
  TH1F *hbc=(TH1F*)bc.Get(histname.c_str());
  GetInt(hbc,xmin,xmax,integralfull,integral);
  fout<<endl; 
  GetIntg(hbc,xmin,xmax,count1, count2);

  fout<<"---------TT bar--------"<<endl; 
  TH1F *htt=(TH1F*)tt.Get(histname.c_str());
  GetInt(htt,xmin,xmax,integralfull,integral);
  fout<<endl; 
  GetIntg(htt,xmin,xmax,count1, count2);
  
  fout<<"---------Wjets--------"<<endl; 
  TH1F *hwj=(TH1F*)wj.Get(histname.c_str());
  GetInt(hwj,xmin,xmax,integralfull,integral);
  fout<<endl; 
  GetIntg(hwj,xmin,xmax,count1, count2);
  
  fout<<"---------WW--------"<<endl; 
  TH1F *hww=(TH1F*)ww.Get(histname.c_str());
  GetInt(hww,xmin,xmax,integralfull,integral);
  fout<<endl; 
  GetIntg(hww,xmin,xmax,count1, count2);

  fout<<endl; 

  fout<<"------------Total BackGround --------"<<endl; 

  fout<<"Full Integral: "<<count1<<endl;
  //  fout<<"Surviving the cut: "<<count2<<endl;
  double eff=count2/count1;
  //  fout<<"Efficiency: "<<eff<<endl;
  //  fout<<endl;
  //  fout<<"Rejection rate: "<<1-eff<<endl;
  
  double M_S=integral/sqrt(integral+count2);
  //  fout<<"Merit of Significance: "<<M_S<<endl;

}


void GetIntg(TH1F *hzj,float xmin,float xmax,double& count1, double& count2)
{
 
  TAxis *axis = hzj->GetXaxis();
  int bmin = axis->FindBin(xmin); 
  int bmax = axis->FindBin(xmax); 
  double integralfull = hzj->Integral();
  double integral = hzj->Integral(bmin,bmax);
  integral -= hzj->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
  integral -= hzj->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/axis->GetBinWidth(bmax);
  count1=count1+integralfull;
  count2=count2+integral;
  // cout<<count1<<'\t'<<count2<<endl;

}

void GetInt(TH1F *hzj,float xmin,float xmax,double& integralfull, double& integral)
{
  TAxis *axis = hzj->GetXaxis();
  int bmin = axis->FindBin(xmin); 
  int bmax = axis->FindBin(xmax); 
  integralfull = hzj->Integral();
  integral = hzj->Integral(bmin,bmax);
  integral -= hzj->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
  integral -= hzj->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/axis->GetBinWidth(bmax);
  fout<<"Full Integral: "<<integralfull<<endl;
}
