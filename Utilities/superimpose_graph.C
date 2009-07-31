{
//created on july 22, 2009 to superimpose graphs, Lovedeep, Anil, Bhawan 
#include "/home/lovedeep/mystyle.C"
setTDRStyle();
tdrStyle->SetPadLeftMargin(0.2);
tdrStyle->SetPadRightMargin(0.1);
tdrStyle->SetPadTopMargin(0.08);
tdrStyle->SetLegendBorderSize(0);
TCanvas* c1= new TCanvas("c1","",7,26,600,600);  

TPad *pad = new TPad("pad","",0,0,1,1);
pad->SetFillColor(0);
//pad->SetGrid();
pad->Draw();
pad->cd();

// draw a frame to define the range
TH1F *hr = c1->DrawFrame(0,0,0.02,1.1);
hr->SetXTitle("#Delta #eta");//#sigma_{i#eta i#eta}");
hr->SetYTitle("efficiency");
hr->GetXaxis()->SetLabelSize(0.03);
hr->GetYaxis()->SetLabelSize(0.03);
hr->GetXaxis()->SetTickLength(0.02);//
hr->GetYaxis()->SetTickLength(0.02);//
pad->GetFrame()->SetFillColor(21);
pad->GetFrame()->SetBorderSize(12);

//create 1 graph
Int_t i;
Float_t x[50];
Float_t y[50];
std::ifstream inputfile("Seff.txt");
float x1,y1;
while(!inputfile.eof())
{
  inputfile>>x1>>y1;
  if(x1==0)continue;
  x[i]=x1;
  y[i]=y1;
  i++;
}

//superimposed graph
Int_t j;
Float_t u[50];
Float_t v[50];
std::ifstream inputfile("Bak_QCD.txt");
float u1,v1;
while(!inputfile.eof())
{
  inputfile>>u1>>v1;
  if(u1==0)continue;
  u[j]=u1;
  v[j]=v1;
  j++;
}

Int_t jt;
Float_t ut[50];
Float_t vt[50];
std::ifstream inputfile("Bak_ttbar.txt");
float ut1,vt1;
while(!inputfile.eof())
{
  inputfile>>ut1>>vt1;
  if(ut1==0)continue;
  ut[jt]=ut1;
  vt[jt]=vt1;
  jt++;
}

Int_t jw;
Float_t uw[50];
Float_t vw[50];
std::ifstream inputfile("Bak_WJets.txt");
float uw1,vw1;
while(!inputfile.eof())
{
  inputfile>>uw1>>vw1;
  if(uw1==0)continue;
  uw[jw]=uw1;
  vw[jw]=vw1;
  jw++;
}



grs=new TGraph(i,x,y);
grq=new TGraph(j,u,v);
grt=new TGraph(j,ut,vt);
grw=new TGraph(j,uw,vw);
grs->SetMarkerColor(kBlue);
grq->SetMarkerColor(kGreen);
grt->SetMarkerColor(kMagenta);
grw->SetMarkerColor(kRed);
grs->SetMarkerStyle(20);
grq->SetMarkerStyle(21);
grt->SetMarkerStyle(23);
grw->SetMarkerStyle(29);
grs->SetLineColor(kBlue);
grq->SetLineColor(kGreen);
grt->SetLineColor(kMagenta);
grw->SetLineColor(kRed);

TMultiGraph *mg = new TMultiGraph();
mg->Add(grs);
mg->Add(grq);
mg->Add(grt);
mg->Add(grw);
mg->Draw("LP");
// gr21->SetMarkerStyle(20);
// gr1->Draw("C*");
// gr21->Draw("CP*");


// create second graph

Int_t jj;
Float_t xx[50];
Float_t yy[50];
std::ifstream inputfile("Ms.txt");
float x11,y11;
while(!inputfile.eof())
    {
      inputfile>>x11>>y11;
      if(x11==0)continue;
      xx[jj]=x11;
      yy[jj]=y11;
      jj++;
    }
 
   //create a transparent pad drawn on top of the main pad


   c1->cd();
   TPad *overlay = new TPad("overlay","",0,0,1,1);
   overlay->SetFillStyle(4000);
   overlay->SetFillColor(0);
   overlay->SetFrameFillStyle(4000);
   overlay->Draw();
   overlay->cd();

    gr2 = new TGraph(jj,xx,yy);
    gr2->SetMarkerColor(kGray);
    gr2->SetLineColor(kGray);
    gr2->SetMarkerStyle(22);
    gr2->SetName("gr2");



   Double_t xmin = pad->GetUxmin();
// Double_t ymin = 115;
   Double_t ymin = 20;
   Double_t xmax = pad->GetUxmax();
//   Double_t ymax = 200;
   Double_t ymax = 45;
   TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
   hframe->GetXaxis()->SetLabelOffset(99);
   hframe->GetYaxis()->SetTickLength(0.00);
   hframe->GetYaxis()->SetLabelOffset(99);
   gr2->Draw("LP");
      
   //Draw an axis on the right side
   TGaxis *axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
   axis->SetLineColor(kGray);
   axis->SetLabelColor(kGray);
   axis->Draw();


  leg = new TLegend(0.7031873,0.7029289,0.8685259,0.8953975,NULL,"brNDC"); //coordinates are fractions
  leg->SetBorderSize(0); 
  leg->SetFillColor(0);
  leg->AddEntry(grs,"Eff. Z_jets","p"); 
  leg->AddEntry(grq,"Eff. QCD","p"); 
  leg->AddEntry(grt,"Eff. ttbar","p"); 
  leg->AddEntry(grw,"Eff. WJets","p"); 
  leg->AddEntry(gr2,"S / #sqrt{S+B}","p"); 
  // use "f" for a box
  leg->Draw();
  c1->Draw();

 }

