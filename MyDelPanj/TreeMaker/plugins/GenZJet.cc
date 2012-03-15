
// system include files
#include <memory>
#include "TLorentzVector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"



//Predicate to sort LorentzVectors.

bool DescendingOrder(const TLorentzVector& l1, const TLorentzVector& l2){
  return l1.Pt()>l2.Pt();
}

bool AscendingOrder(const TLorentzVector& l1, const TLorentzVector& l2){
  return l1.Pt()<l2.Pt();
}


namespace genZJet{
  TH1D* prettyHistogram(edm::Service<TFileService> fs, std::string name, std::string title, std::string xtitle, double nbins, double rlow, double rhigh ){
    TH1D* hist = fs->make<TH1D>(name.c_str(),title.c_str(),nbins,rlow,rhigh);
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    return hist;
  }

  
  TH1D* prettyHistogram(edm::Service<TFileService> fs, std::string name, std::string title, std::string xtitle, double nbins, float* xbins){
    TH1D* hist = fs->make<TH1D>(name.c_str(),title.c_str(),nbins,xbins);
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    return hist;
  }

  double Phi_0_2pi(double x) {
    while (x >= 2*M_PI) x -= 2*M_PI;
    while (x <     0.)  x += 2*M_PI;
    return x;
  }

  double radius(double eta1, double phi1,double eta2, double phi2){
    
    const double TWOPI= 2.0*M_PI;
    
    phi1=Phi_0_2pi(phi1);
    phi2=Phi_0_2pi(phi2);
    
    double dphi=Phi_0_2pi(phi1-phi2);
    dphi = TMath::Min(dphi,TWOPI-dphi);
    double deta = eta1-eta2;
    
    return sqrt(deta*deta+dphi*dphi);
  }


}


//
// class declaration
//

class GenZJet : public edm::EDAnalyzer {
   public:
      explicit GenZJet(const edm::ParameterSet&);
      ~GenZJet();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  TH1D* hNjetGoodZX; 
  TH1D* hEta1 ; 
  TH1D* FirstJetPt_Z1jet;
  TH1D* SecondJetPt_Z2jet;
  TH1D* ThirdJetPt_Z3jet;
  TH1D* FourthJetPt_Z4jet;
  TH1D* FifthJetPt_Z5jet;
  TH1D* SixthJetPt_Z6jet;
  
  TH1D* hEffi_ ; 
  TH1D* hNjetGoodIncZX ; 
  TH1D* ZMass_Zinc0jet;
  TH1D* ZMass_Zinc1jet;
  TH1D* ZMass_Zinc2jet;
  TH1D* ZMass_Zinc3jet;
  TH1D* ZMass_Zinc4jet;
  TH1D* ZMass_Zinc5jet;
  TH1D* ZMass_Zinc6jet;
  
  TH1D* ZMass_Zexc0jet;
  TH1D* ZMass_Zexc1jet;
  TH1D* ZMass_Zexc2jet;
  TH1D* ZMass_Zexc3jet;
  TH1D* ZMass_Zexc4jet;
  TH1D* ZMass_Zexc5jet;
  TH1D* ZMass_Zexc6jet;
  
  TH1D* ZPt_Zinc0jet;
  TH1D* ZPt_Zinc1jet;
  TH1D* ZPt_Zinc2jet;
  TH1D* ZPt_Zinc3jet;
  TH1D* ZPt_Zinc4jet;
  TH1D* ZPt_Zinc5jet;
  TH1D* ZPt_Zinc6jet;
  
  
  TH1D* ZPt_Zexc0jet;
  TH1D* ZPt_Zexc1jet;
  TH1D* ZPt_Zexc2jet;
  TH1D* ZPt_Zexc3jet;
  TH1D* ZPt_Zexc4jet;
  TH1D* ZPt_Zexc5jet;
  TH1D* ZPt_Zexc6jet;
  
  
  TH1D* ZRapidity_Zinc0jet;
  TH1D* ZRapidity_Zinc1jet;
  TH1D* ZRapidity_Zinc2jet;
  TH1D* ZRapidity_Zinc3jet;
  TH1D* ZRapidity_Zinc4jet;
  TH1D* ZRapidity_Zinc5jet;
  TH1D* ZRapidity_Zinc6jet;
  
  
  TH1D* ZRapidity_Zexc0jet;
  TH1D* ZRapidity_Zexc1jet;
  TH1D* ZRapidity_Zexc2jet;
  TH1D* ZRapidity_Zexc3jet;
  TH1D* ZRapidity_Zexc4jet;
  TH1D* ZRapidity_Zexc5jet;
  TH1D* ZRapidity_Zexc6jet;
  
  
  TH1D* FirstJetY_Z1jet;
  TH1D* SecondJetY_Z2jet;
  TH1D* ThirdJetY_Z3jet;
  TH1D* FourthJetY_Z4jet;
  TH1D* FifthJetY_Z5jet;
  TH1D* SixthJetY_Z6jet;
  
  TH1D* hAccept ; 
  TH1D* hCounter;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenZJet::GenZJet(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
   hNjetGoodZX = genZJet::prettyHistogram(fs,"h_Njet_GoodZX","N_{Jets}(Exc) [Z+Jets]","NJets_{Exc}",10,0,10);
  
  float jpt[10] = {30,50,70,90,120,150,180,210,300,400};
  float jpt2[9] = {30,50,70,90,120,150,200,270,400};
  float jpt3[6] = {30,50,70,90,120,400};
  float jpt4[6] = {30,50,70,90,120,400};
  
   hEta1 = new TH1D("hEta1","hEta2",30,-3,3);
   FirstJetPt_Z1jet=genZJet::prettyHistogram(fs,"FirstJetPt_Z1jet","First Leading Jet (Z+Jets)","p_{T} [GeV]",9,jpt);
   SecondJetPt_Z2jet=genZJet::prettyHistogram(fs,"SecondJetPt_Z2jet","Second Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",8,jpt2);
   ThirdJetPt_Z3jet=genZJet::prettyHistogram(fs,"ThirdJetPt_Z3jet","Third Leading (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt3);
   FourthJetPt_Z4jet=genZJet::prettyHistogram(fs,"FourthJetPt_Z4jet","Fourth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",5,jpt4);
   FifthJetPt_Z5jet=genZJet::prettyHistogram(fs,"FifthJetPt_Z5jet","Fifth Leading Jet (Z+Jets)","p_{T}(Jet) [GeV]",9,jpt);
   SixthJetPt_Z6jet=genZJet::prettyHistogram(fs,"SixthJetPt_Z6jet","#geq 6th Jets (Z+Jets)","p_{T}(Jet) [GeV]",9,jpt);
  
   hEffi_ = genZJet::prettyHistogram(fs,"hEffi_","Weight","Weight",50,0,1);
   hNjetGoodIncZX = genZJet::prettyHistogram(fs,"h_Njet_GoodIncZX","N_{Jets}(Exc) [Z+Jets]","NJets_{Inc}",10,0,10);
   ZMass_Zinc0jet=genZJet::prettyHistogram(fs,"ZMass_Zinc0jet","Z Invariant Mass (N_{jets} #geq 0)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zinc1jet=genZJet::prettyHistogram(fs,"ZMass_Zinc1jet","Z Invariant Mass (N_{jets} #geq 1)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zinc2jet=genZJet::prettyHistogram(fs,"ZMass_Zinc2jet","Z Invariant Mass (N_{jets} #geq 2)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zinc3jet=genZJet::prettyHistogram(fs,"ZMass_Zinc3jet","Z Invariant Mass (N_{jets} #geq 3)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zinc4jet=genZJet::prettyHistogram(fs,"ZMass_Zinc4jet","Z Invariant Mass (N_{jets} #geq 4)","M_{ee} [GeV]",40,60.,120.);
   ZMass_Zinc5jet=genZJet::prettyHistogram(fs,"ZMass_Zinc5jet","Z Invariant Mass (N_{jets} #geq 5)","M_{ee} [GeV]",40,60.,120.);
   ZMass_Zinc6jet=genZJet::prettyHistogram(fs,"ZMass_Zinc6jet","Z Invariant Mass (N_{jets} #geq 6)","M_{ee} [GeV]",40,60.,120.);

   ZMass_Zexc0jet=genZJet::prettyHistogram(fs,"ZMass_Zexc0jet","Z Invariant Mass (N_{jets} = 0)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zexc1jet=genZJet::prettyHistogram(fs,"ZMass_Zexc1jet","Z Invariant Mass (N_{jets} = 1)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zexc2jet=genZJet::prettyHistogram(fs,"ZMass_Zexc2jet","Z Invariant Mass (N_{jets} = 2)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zexc3jet=genZJet::prettyHistogram(fs,"ZMass_Zexc3jet","Z Invariant Mass (N_{jets} = 3)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zexc4jet=genZJet::prettyHistogram(fs,"ZMass_Zexc4jet","Z Invariant Mass (N_{jets} = 4)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zexc5jet=genZJet::prettyHistogram(fs,"ZMass_Zexc5jet","Z Invariant Mass (N_{jets} = 5)","M_{ee} [GeV]",60,60.,120.);
   ZMass_Zexc6jet=genZJet::prettyHistogram(fs,"ZMass_Zexc6jet","Z Invariant Mass (N_{jets} = 6)","M_{ee} [GeV]",20,60.,120.);

   ZPt_Zinc0jet=genZJet::prettyHistogram(fs,"ZPt_Zinc0jet","Z p_{T} (N_{jets} #geq 0)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zinc1jet=genZJet::prettyHistogram(fs,"ZPt_Zinc1jet","Z p_{T} (N_{jets} #geq 1)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zinc2jet=genZJet::prettyHistogram(fs,"ZPt_Zinc2jet","Z p_{T} (N_{jets} #geq 2)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zinc3jet=genZJet::prettyHistogram(fs,"ZPt_Zinc3jet","Z p_{T} (N_{jets} #geq 3)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zinc4jet=genZJet::prettyHistogram(fs,"ZPt_Zinc4jet","Z p_{T} (N_{jets} #geq 4)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zinc5jet=genZJet::prettyHistogram(fs,"ZPt_Zinc5jet","Z p_{T} (N_{jets} #geq 5)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zinc6jet=genZJet::prettyHistogram(fs,"ZPt_Zinc6jet","Z p_{T} (N_{jets} #geq 6)","p_{T}(Z) [GeV]",25,0.,500.);


   ZPt_Zexc0jet=genZJet::prettyHistogram(fs,"ZPt_Zexc0jet","Z p_{T} (N_{jets} = 0)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zexc1jet=genZJet::prettyHistogram(fs,"ZPt_Zexc1jet","Z p_{T} (N_{jets} = 1)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zexc2jet=genZJet::prettyHistogram(fs,"ZPt_Zexc2jet","Z p_{T} (N_{jets} = 2)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zexc3jet=genZJet::prettyHistogram(fs,"ZPt_Zexc3jet","Z p_{T} (N_{jets} = 3)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zexc4jet=genZJet::prettyHistogram(fs,"ZPt_Zexc4jet","Z p_{T} (N_{jets} = 4)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zexc5jet=genZJet::prettyHistogram(fs,"ZPt_Zexc5jet","Z p_{T} (N_{jets} = 5)","p_{T}(Z) [GeV]",25,0.,500.);
   ZPt_Zexc6jet=genZJet::prettyHistogram(fs,"ZPt_Zexc6jet","Z p_{T} (N_{jets} = 6)","p_{T}(Z) [GeV]",25,0.,500.);


   ZRapidity_Zinc0jet=genZJet::prettyHistogram(fs,"ZRapidity_Zinc0jet","ZRapidity (N_{jets} #geq 0)","Y(Z)",30,-3.,3.);
   ZRapidity_Zinc1jet=genZJet::prettyHistogram(fs,"ZRapidity_Zinc1jet","ZRapidity (N_{jets} #geq 1)","Y(Z)",30,-3.,3.);
   ZRapidity_Zinc2jet=genZJet::prettyHistogram(fs,"ZRapidity_Zinc2jet","ZRapidity (N_{jets} #geq 2)","Y(Z)",30,-3.,3.);
   ZRapidity_Zinc3jet=genZJet::prettyHistogram(fs,"ZRapidity_Zinc3jet","ZRapidity (N_{jets} #geq 3)","Y(Z)",30,-3.,3.);
   ZRapidity_Zinc4jet=genZJet::prettyHistogram(fs,"ZRapidity_Zinc4jet","ZRapidity (N_{jets} #geq 4)","Y(Z)",30,-3.,3.);
   ZRapidity_Zinc5jet=genZJet::prettyHistogram(fs,"ZRapidity_Zinc5jet","ZRapidity (N_{jets} #geq 5)","Y(Z)",30,-3.,3.);
   ZRapidity_Zinc6jet=genZJet::prettyHistogram(fs,"ZRapidity_Zinc6jet","ZRapidity (N_{jets} #geq 6)","Y(Z)",30,-3.,3.);


   ZRapidity_Zexc0jet=genZJet::prettyHistogram(fs,"ZRapidity_Zexc0jet","ZRapidity (N_{jets} = 0)","Y(Z)",10,0.,3.);
   ZRapidity_Zexc1jet=genZJet::prettyHistogram(fs,"ZRapidity_Zexc1jet","ZRapidity (N_{jets} = 1)","Y(Z)",10,0.,3.);
   ZRapidity_Zexc2jet=genZJet::prettyHistogram(fs,"ZRapidity_Zexc2jet","ZRapidity (N_{jets} = 2)","Y(Z)",10,0.,3.);
   ZRapidity_Zexc3jet=genZJet::prettyHistogram(fs,"ZRapidity_Zexc3jet","ZRapidity (N_{jets} = 3)","Y(Z)",10,0.,3.);
   ZRapidity_Zexc4jet=genZJet::prettyHistogram(fs,"ZRapidity_Zexc4jet","ZRapidity (N_{jets} = 4)","Y(Z)",10,0.,3.);
   ZRapidity_Zexc5jet=genZJet::prettyHistogram(fs,"ZRapidity_Zexc5jet","ZRapidity (N_{jets} = 5)","Y(Z)",10,0.,3.);
   ZRapidity_Zexc6jet=genZJet::prettyHistogram(fs,"ZRapidity_Zexc6jet","ZRapidity (N_{jets} = 6)","Y(Z)",10,0.,3.);


   FirstJetY_Z1jet=genZJet::prettyHistogram(fs,"FirstJetY_Z1jet","First Leading Jet (Z+Jets)","Y (Jet)",10,0.,3.);
   SecondJetY_Z2jet=genZJet::prettyHistogram(fs,"SecondJetY_Z2jet","Second Leading Jet (Z+Jets)","Y(Jet) ",10,0.,3.);
   ThirdJetY_Z3jet=genZJet::prettyHistogram(fs,"ThirdJetY_Z3jet","Third Leading (Z+Jets)","Y(Jets) ",10,0.,3.);
   FourthJetY_Z4jet=genZJet::prettyHistogram(fs,"FourthJetY_Z4jet","Fourth Leading Jet (Z+Jets)","Y(Jet) ",10,0.,3.);
   FifthJetY_Z5jet=genZJet::prettyHistogram(fs,"FifthJetY_Z5jet","Fifth Leading Jet (Z+Jets)","Y(Jet) ",10,0.,3.);
   SixthJetY_Z6jet=genZJet::prettyHistogram(fs,"SixthJetY_Z6jet","#geq 6th Jets (Z+Jets)","Y(Jet) ",10,0.,3.);

   hAccept = genZJet::prettyHistogram(fs,"hAccept","Accept","A",50,0.,1.);
   hCounter = genZJet::prettyHistogram(fs,"hCounter","Counter","A",10,0.,10.);


  
}


GenZJet::~GenZJet()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenZJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace edm;
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  if(not iEvent.getByLabel("genParticles", genParticleHandle))
    {
      std::cout<<
        "GenAnalyzer: Generator Level Information not found\n"
	       <<std::endl;
    }
  

  std::vector<TLorentzVector> electrons ;
  const reco::GenParticleCollection* genColl= &(*genParticleHandle);
  reco::GenParticleCollection::const_iterator geni = genColl->begin();
  for(; geni!=genColl->end();geni++){
    reco::GenParticle gen = *geni;
    //Look out for the GenMuons/Electrons
    if(gen.status()==1 && fabs(gen.pdgId())==13){
      TLorentzVector* ele = new TLorentzVector();
      ele->SetPxPyPzE(gen.px(),gen.py(),gen.pz(),gen.energy());
      electrons.push_back(*ele);
    }
  }
  if(electrons.size()<2)return;
  std::sort(electrons.begin(),electrons.end(),DescendingOrder);
  TLorentzVector p((electrons[0]+electrons[1]));
  if((60>p.M())||(120<p.M()))return;

  //Loop over GenJets.

  edm::Handle<reco::GenJetCollection> genJetsHandle;
   if( not iEvent.getByLabel("ak5GenJets",genJetsHandle)){
     edm::LogInfo("GenAnalyzer") << "genJets not found, "
       "skipping event";
     return;
   }

   std::vector<TLorentzVector> jets;
   const reco::GenJetCollection* genJetColl = &(*genJetsHandle);
   reco::GenJetCollection::const_iterator gjeti = genJetColl->begin();
   
   for(; gjeti!=genJetColl->end();gjeti++){
     reco::GenParticle gjet = *gjeti;
     if(gjet.pt()<=30)continue;
     if(fabs(gjet.eta())< 2.4)continue;
     double delR_e1 = genZJet::radius(gjet.eta(),gjet.phi(),(electrons[0]).Eta(), (electrons[0]).Phi());
     double delR_e2= genZJet::radius(gjet.eta(),gjet.phi(),(electrons[1]).Eta(), (electrons[1]).Phi());
     if((delR_e1<0.3)||(delR_e2<0.3))continue;
     TLorentzVector* jet = new TLorentzVector();
     jet->SetPxPyPzE(gjet.px(),gjet.py(),gjet.pz(),gjet.energy());
     jets.push_back(*jet);  
 }
   
   double weight =1 ;
  //std::cout<<"1: "<<std::endl; 
   hNjetGoodIncZX->Fill(0.,weight);
   ZMass_Zinc0jet->Fill(p.M(),weight);
  //std::cout<<"12: "<<std::endl;

   ZPt_Zinc0jet->Fill(p.Pt(),weight);
   ZRapidity_Zinc0jet->Fill(fabs(p.Rapidity()),weight);
  //std::cout<<"13: "<<std::endl;
 
  if(jets.size()==0){
     ZMass_Zexc0jet->Fill(p.M(),weight);
     //std::cout<<"14: "<<std::endl;

     ZPt_Zexc0jet->Fill(p.Pt(),weight);
//std::cout<<"15: "<<std::endl;

     ZRapidity_Zexc0jet->Fill(fabs(p.Rapidity()),weight);
//std::cout<<"16: "<<std::endl;

     hCounter->Fill(0.,weight);
//std::cout<<"17: "<<std::endl;

   }
        //std::cout<<"2: "<<std::endl;

   if(jets.size()>=1){
     
     hNjetGoodIncZX->Fill(1., weight);
     ZMass_Zinc1jet->Fill(p.M(),weight);
     ZPt_Zinc1jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc1jet->Fill(fabs(p.Rapidity()),weight);
     if(jets.size()==1){
       ZMass_Zexc1jet->Fill(p.M(),weight);
       ZPt_Zexc1jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc1jet->Fill(fabs(p.Rapidity()),weight);
       FirstJetPt_Z1jet->Fill(jets[0].Pt(),weight);
       FirstJetY_Z1jet->Fill(fabs(jets[0].Eta()),weight);
       hCounter->Fill(1.,weight);
     }
   }
   
     //std::cout<<"3: "<<std::endl;
   if(jets.size()>=2){
     
     hNjetGoodIncZX->Fill(2., weight);
     ZMass_Zinc2jet->Fill(p.M(),weight);
     ZPt_Zinc2jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc2jet->Fill(fabs(p.Rapidity()),weight);
     if(jets.size()==2){
       ZMass_Zexc2jet->Fill(p.M(),weight);
       ZPt_Zexc2jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc2jet->Fill(fabs(p.Rapidity()),weight);
       SecondJetPt_Z2jet->Fill(jets[1].Pt(),weight);
       SecondJetY_Z2jet->Fill(fabs(jets[1].Eta()),weight);
       hCounter->Fill(2.,weight);
     }
   }
     //std::cout<<"4: "<<std::endl;

   if(jets.size()>=3){
     
     hNjetGoodIncZX->Fill(3., weight);
     ZMass_Zinc3jet->Fill(p.M(),weight);
     ZPt_Zinc3jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc3jet->Fill(fabs(p.Rapidity()),weight);
     if(jets.size()==3){
       ZMass_Zexc3jet->Fill(p.M(),weight);
       ZPt_Zexc3jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc3jet->Fill(fabs(p.Rapidity()),weight);
       ThirdJetPt_Z3jet->Fill(jets[2].Pt(),weight);
       ThirdJetY_Z3jet->Fill(fabs(jets[2].Eta()),weight);
       hCounter->Fill(3.,weight);
     }
   }
   
  //std::cout<<"5: "<<std::endl;

   if(jets.size()>=4){
     
     hNjetGoodIncZX->Fill(4., weight);
     ZMass_Zinc4jet->Fill(p.M(),weight);
     ZPt_Zinc4jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc4jet->Fill(fabs(p.Rapidity()),weight);
     if(jets.size()==4){
       ZMass_Zexc4jet->Fill(p.M(),weight);
       ZPt_Zexc4jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc4jet->Fill(fabs(p.Rapidity()),weight);
       FourthJetPt_Z4jet->Fill(jets[3].Pt(),weight);
       FourthJetY_Z4jet->Fill(fabs(jets[3].Eta()),weight);
       hCounter->Fill(4.,weight);
     }
   }
  //std::cout<<"6: "<<std::endl;

   if(jets.size()>=5){
     
     hNjetGoodIncZX->Fill(5., weight);
     ZMass_Zinc5jet->Fill(p.M(),weight);
     ZPt_Zinc5jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc5jet->Fill(p.Rapidity(),weight);
     if(jets.size()==5){
       ZMass_Zexc5jet->Fill(p.M(),weight);
       ZPt_Zexc5jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc5jet->Fill(p.Rapidity(),weight);
       FifthJetPt_Z5jet->Fill(jets[4].Pt(),weight);
       FifthJetY_Z5jet->Fill(jets[4].Eta(),weight);
       hCounter->Fill(5.,weight);
     }
   }
   
   if(jets.size()>=6){
     
     hNjetGoodIncZX->Fill(6., weight);
     ZMass_Zinc6jet->Fill(p.M(),weight);
     ZPt_Zinc6jet->Fill(p.Pt(),weight);
     ZRapidity_Zinc6jet->Fill(p.Rapidity(),weight);
     SixthJetPt_Z6jet->Fill(jets[5].Pt(),weight);
     SixthJetY_Z6jet->Fill(jets[5].Eta(),weight);
     if(jets.size()==6){
       ZMass_Zexc6jet->Fill(p.M(),weight);
       ZPt_Zexc6jet->Fill(p.Pt(),weight);
       ZRapidity_Zexc6jet->Fill(p.Rapidity(),weight);
       hCounter->Fill(6.,weight);
     }
     
   }
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenZJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenZJet::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
GenZJet::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
GenZJet::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
GenZJet::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
GenZJet::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenZJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenZJet);
