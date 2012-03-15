rootlogon(bool compile = true)
{
   gStyle->SetCanvasColor(0);

   gStyle->SetHistLineWidth(2);  
   gStyle->SetOptStat(111111);
   G__loadfile("setTDRStyle.h");
   setTDRStyle();
   gStyle->SetTitleYOffset(1.5);//1.5

   //      gStyle->SetFillColor(0);    // White
   ////   gStyle->SetFillStyle(4000); // Transparent

   gStyle->SetPalette(1,0);
   
   gStyle->SetHistLineWidth(2);   
   gStyle->SetLegendBorderSize(1);
   gStyle->SetFillColor(0);

//gSystem->Load("libMinuit.so");
   //gSystem->Load("libRooFitCore.so");
   //gSystem->Load("libRooFitModels.so");
   //using namespace RooFit;

//    gSystem->Load("libFWCoreFWLite.so");
//    AutoLibraryLoader::enable();
//    gSystem->Load("libDataFormatsFWLite.so");

}
