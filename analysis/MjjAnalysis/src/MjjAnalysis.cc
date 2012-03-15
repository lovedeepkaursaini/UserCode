// -*- C++ -*-
//
// Package:    MjjAnalysis
// Class:      MjjAnalysis
// 
/**\class MjjAnalysis MjjAnalysis.cc analysis/MjjAnalysis/src/MjjAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Anil Pratap Singh,,,
//         Created:  Sat Jul  9 10:37:13 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//nclude "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TLorentzVector.h"
#include "analysis/MjjAnalysis/interface/JetUtils.h"

TH1D* prettyHistograms(TH1D* hist,
		       edm::Service<TFileService> fs,
		       std::string token,
		       std::string qty,
		       std::string cut,
		       std::string xtitle,
		       int xbins,
		       double xrangeD,
		       double xrangeU
		       ){
  std::string name = "h_"+token+"_"+qty+"_"+"_"+cut;
  std::string title = token+" "+qty+"("+cut+")";
  hist =fs->make<TH1D>(name.c_str(), title.c_str(),xbins,xrangeD,xrangeU);
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  return hist;
}

class MjjPlots{

public:
  MjjPlots(std::string tag, edm::Service<TFileService> fs);
  void Fill(TLorentzVector j0,TLorentzVector j1, TLorentzVector dau1, 
       TLorentzVector met, TLorentzVector W, TLorentzVector l);
  TH1D* hWPt_;
  TH1D* hWMt_;
  TH1D* hEtj0_;
  TH1D* hEtj1_;
  TH1D* hEtaj0_;
  TH1D* hEtaj1_;

  TH1D* hMj0j1_;
  TH1D* hPtj0j1_;
  TH1D* hdEtaj0j1_; 
  TH1D* hdPhij0j1_ ;
  TH1D* hdRj0j1_ ;
  TH1D* hMj0j1W_ ;
  TH1D* hptRatioj0j1_; 
  TH1D* hdPhij0l0_;
  TH1D* hdPhij1l0_;
  TH1D* hdPhil1Met_;
   };



MjjPlots::MjjPlots(std::string tag,edm::Service<TFileService> fs){
  hWPt_ = prettyHistograms(hWPt_,fs,"W","TransverseMom",tag,"W p_{T} [GeV]",50,0,200.);
  hWMt_ = prettyHistograms(hWMt_,fs,"W","TransverseMass",tag,"W M_{T} [GeV]", 40,0.,200.);
  hMj0j1_  = prettyHistograms(hMj0j1_,fs,"Dijet","InvMass",tag,"M_{j0,j1} [GeV]",60, 0.,300.);
  hMj0j1W_ = prettyHistograms(hMj0j1W_,fs,"DijetPlusW", "InvMass", tag,"M_{j0j1,W} [GeV]",60, 0.,300.);
  hdEtaj0j1_= prettyHistograms(hdEtaj0j1_,fs,"Dijet","DeltaEta_LeadJetSubleadJet",tag,"#Delta#eta_{j0,j1}",40, -5.,5.); 
  hdPhij0j1_ = prettyHistograms(hdPhij0j1_,fs,"Dijet","DeltaPhi_LeadJetSubleadJet", tag,"#Delta#phi_{j0,j1}",40, -5.,5.);
  hdRj0j1_ = prettyHistograms(hdRj0j1_,fs,"Dijet", "DeltaR_LeadJetSubleadJet",tag,"#Delta R_{j0,j1}", 20, 0.,1.);
  hEtj0_ = prettyHistograms(hEtj0_,fs,"Dijet","LeadJet_Et", tag,"E_{T} [GeV]",60, 0.,300.);
  hEtj1_ = prettyHistograms(hEtj1_,fs,"Dijet","SubleadJet_Et",tag,"E_{T} [GeV]", 60, 0.,300.);
  hEtaj0_ = prettyHistograms(hEtaj0_,fs,"Dijet","LeadJet_Eta", tag,"#eta",20, -5.,5.);
  hEtaj1_ = prettyHistograms(hEtaj1_,fs,"Dijet","SubleadJet_Eta", tag,"#eta",20, -5.,5.);
  hdPhij0l0_ = prettyHistograms(hdPhij0l0_,fs,"Dijet","DeltaPhi_LeadJetLepton",tag,"#Delta#phi_{j0,l0}", 20,-5.,5.);
  hdPhij1l0_= prettyHistograms(hdPhij1l0_,fs,"Dijet","DeltaPhi_SubleadJetLepton" , tag,"#Delta#phi_{j1,l0}",20,-5.,5.);
  hdPhil1Met_ = prettyHistograms(hdPhil1Met_,fs,"Dijet","DeltaPhi_leptonMet",tag,"#Delta#phi_{l0,Met}", 20, -5.,5.);
}

void MjjPlots::Fill(TLorentzVector j0,TLorentzVector j1, 
	       TLorentzVector dau1, TLorentzVector met, 
	       TLorentzVector W, TLorentzVector lepPlusMet){
  TLorentzVector dijet = j0+j1;
  TLorentzVector dijetW = dijet+W;
  hWPt_ ->Fill(dijet.Pt());
  hWMt_ ->Fill(lepPlusMet.M());
  hEtj0_->Fill(j0.Pt());
  hEtj1_->Fill(j1.Pt());
  hEtaj0_->Fill(j0.Eta());
  hEtaj1_->Fill(j1.Eta());
  hMj0j1_->Fill(dijet.M());  
  hMj0j1W_ ->Fill( dijetW.M());
  hdEtaj0j1_ ->Fill( dEta(j0.Eta(),j1.Eta()));
  hdPhij0j1_ ->Fill( dPhi(j0.Phi(),j1.Phi()));
  hdRj0j1_ ->Fill( radius(j0.Eta(),j0.Phi(),j1.Eta(),j1.Phi()));
  hdPhij0l0_->Fill(dPhi(j0.Phi(),dau1.Phi()));
  hdPhij1l0_->Fill(dPhi(j1.Phi(),dau1.Phi()));
  hdPhil1Met_->Fill(dPhi(dau1.Phi(),met.Phi()));
}


//
// class declaration
//
double deltaR(){return 0;}

class PtGreater1 {
public:
  template <typename T> bool operator () (const T& i, const T& j) {  return (i.Pt() > j.Pt()); }
};


class MjjAnalysis : public edm::EDAnalyzer {
   public:
      explicit MjjAnalysis(const edm::ParameterSet&);
      ~MjjAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

     edm::InputTag hepMcColl_;
 

      // ----------member data ---------------------------
  TH1D* hWMass_;
  TH1D* hWPt_;
  TH1D* hWRap_;
  TH1D* hWMt_;
  TH1D* hWRecMt_;
  TH1D* hWRecPt_;
  TH1D* hLepPt_;
  TH1D* hLepEta_;
  TH1D* hpt_;
  TH1D* hRecMet_;
  TH1D* hNeutrinoMet_;
  TH1D* hNjetExclusive_;
  TH1D* hNjetInclusive_;
  TH1D* hdelR;
  TH1D* hdelR2;
  TH1D* hNeutrinoPt_;
  TH1D* hNeutrinoSumPt_;

  TH1D* hCounter_;


 //Dijet Plots...
  MjjPlots* dijet_;

  MjjPlots* dijetPtBelow40_;
  MjjPlots* dijetPtAbove40_;
 
  MjjPlots* dijetMassWindow_;
  MjjPlots* dijetSideBand_;
  

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
MjjAnalysis::MjjAnalysis(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  hNeutrinoPt_ = prettyHistograms(hNeutrinoPt_,fs,"Nu","TransverseMom","Inclusive","Neutrino p_{T} [GeV]",40,0,200.);
  hNeutrinoSumPt_ =  prettyHistograms(hNeutrinoSumPt_,fs,"Nu","SumPt","Inclusive","Neutrino SumPt [GeV]", 40, 0.,200.);
  hWMass_   =  prettyHistograms(hWMass_,fs,"W","Mass","Inclusive","M [GeV]",50, 0.,500.);
  hWPt_     =  prettyHistograms(hWPt_,fs,"W"," TransverseMomentum","Inclusive","p_{T} [GeV]",50, 0.,200.);
  hWRecPt_  =  prettyHistograms(hWRecPt_,fs,"W", "TransverseMomentumRec","Inclusive","p_{T} [GeV]",50,0.,200.);
  hWMt_     =  prettyHistograms(hWMt_,fs,"W"," TransverseMass","Inclusive","m_{T} [GeV]",40,0.,200.);
  hWRecMt_  =  prettyHistograms(hWRecMt_,fs,"W","TransverseMassRec","Inclusive", "m_{T} [GeV]",40,0.,200.);
  hRecMet_ = prettyHistograms(hRecMet_,fs,"Associated","MissingEtRec","Inclusive","Missing E{T} [GeV]", 60, 0.,300.);
  hNeutrinoMet_ = prettyHistograms(hNeutrinoMet_,fs,
				   "Associated","MissingEtNeutrino","Inclusive","Missing E{T} [GeV]",60, 0.,300.);
  hNjetExclusive_= prettyHistograms(hNjetExclusive_,fs,"Associated", "JetMulti","Exclusive","N_{jet}",10, 0.,10.);
  hNjetInclusive_= prettyHistograms(hNjetInclusive_,fs,"Associated", "JetMulti","Inclusive","N_{jet}",10, 0.,10.);
  hLepPt_  = prettyHistograms(hLepPt_,fs,"W","LeptonTransverseMomentum","Inclusive","p_{T} [GeV]",50,0.,200.);
  hLepEta_  = prettyHistograms(hLepEta_,fs,"W","LeptonPseudorapidity","Inclusive","#eta [GeV]",20,-5.,5.);

  dijetPtBelow40_=new MjjPlots("PtBelow40",fs);
  dijetPtAbove40_=new MjjPlots("PtAbove40",fs);
  dijetMassWindow_=new MjjPlots("MassWindow",fs);
  dijetSideBand_=new MjjPlots("SideBand",fs);
  dijet_=new MjjPlots("WPlusDijet",fs);


  hdelR = fs->make<TH1D>("delR","delR",400,0,40);
  hdelR2 = fs->make<TH1D>("delR2","delR2",400,0,40);

  hpt_ = fs->make<TH1D>("hpt_", "hpt_", 50, 0.,500.);
  hCounter_= fs->make<TH1D>("hCounter_", "hCounter_", 100, 0.,50.);

}


MjjAnalysis::~MjjAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MjjAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  hCounter_->Fill(0);
  using namespace edm;
  edm::Handle<HepMCProduct> event;
  iEvent.getByLabel("generator", event);
  const HepMC::GenEvent *myGenEvent = event->GetEvent();
  std::vector<const HepMC::GenParticle*> stableGenColl; 
  std::vector<const HepMC::GenParticle*> stableLepColl; 
  std::vector<const HepMC::GenParticle*> allNeutrinoColl;
 
  
  TLorentzVector momsum(0., 0., 0., 0.);  
  TLorentzVector neutrino_momsum(0., 0., 0., 0.);
  for(HepMC::GenEvent::particle_const_iterator giter = myGenEvent->particles_begin(); giter != myGenEvent->particles_end(); ++giter) {
    
    if((*giter)->status()!=1)continue;
    bool isneutrino = false;
    unsigned int pid = fabs((*giter)->pdg_id());   
    if(pid==12 || pid ==14||pid==16){
      isneutrino = true;
      allNeutrinoColl.push_back(*giter);
      TLorentzVector mom((*giter)->momentum().x(), (*giter)->momentum().y(), (*giter)->momentum().z(), (*giter)->momentum().t());
      neutrino_momsum += mom;
    }
    stableGenColl.push_back(*giter);
    if(!(isneutrino)){
      TLorentzVector mom((*giter)->momentum().x(), (*giter)->momentum().y(), (*giter)->momentum().z(), (*giter)->momentum().t());     
      momsum += mom;
      double pt = (*giter)->momentum().perp();
      double eta = fabs((*giter)->momentum().eta());
      if(pid==11 || pid ==13||pid==15){
	if(!(eta>3 ||pt<20))
	  stableLepColl.push_back(*giter);
      }
    }
  }
  //Fill Missing Et stuff.
  
  if(!allNeutrinoColl.size()){std::cout<<"retiring"<<std::endl ;return;}
  hCounter_->Fill(1);
  TLorentzVector met(-1*momsum.Px(),-1*momsum.Py(),0,momsum.Pt());
  
  
  TLorentzVector neutrino(allNeutrinoColl[0]->momentum().x(), allNeutrinoColl[0]->momentum().y(), 
			  allNeutrinoColl[0]->momentum().z(), allNeutrinoColl[0]->momentum().t()); 
  
  TLorentzVector neutrinomet(neutrino_momsum.Px(),neutrino_momsum.Py(),0,neutrino_momsum.Pt());

  TLorentzVector recmet(met);
  
  
  
  hRecMet_->Fill(recmet.Et());
  hNeutrinoPt_->Fill(neutrino.Et());
  hNeutrinoSumPt_->Fill(neutrino_momsum.Pt());
  hNeutrinoMet_->Fill(neutrinomet.Et());

  //first escape from the analysis

  if((stableLepColl.size()+allNeutrinoColl.size())<2)return;
  hCounter_->Fill(2);


  //apply isolation

  std::vector<const HepMC::GenParticle*> stableIsoLepColl;
  for(unsigned int i=0; i!=stableLepColl.size(); i++){
    double etsum = 0;
    for(unsigned int j=0;j!= stableGenColl.size(); j++){
      double delR = radius(stableLepColl[i]->momentum().eta(),stableLepColl[i]->momentum().phi(),
			   stableGenColl[j]->momentum().eta(),stableGenColl[j]->momentum().phi());
      //between lep&gen
      if(delR<0.5 && delR>0.00001){
	etsum += (stableGenColl[j])->momentum().perp();
	if(etsum>0.1*(stableLepColl[i])->momentum().perp());
	break;
      }
    }
    if(etsum<0.1*(stableLepColl[i])->momentum().perp())
      stableIsoLepColl.push_back(stableLepColl[i]);
  }
  
  //2nd escape from the analysis
  if(stableIsoLepColl.size()<1)return ;
  hCounter_->Fill(3);
  
  //3rd escape from the analysis
  if(neutrinomet.Et()<25)return ;
  hCounter_->Fill(4);
  
  
  //reconstruct W boson.
  TLorentzVector dau1(stableIsoLepColl[0]->momentum().x(), stableIsoLepColl[0]->momentum().y(), 
		      stableIsoLepColl[0]->momentum().z(), stableIsoLepColl[0]->momentum().t()); 
  
  
  TLorentzVector w_true = dau1 + neutrino;
  
  //i think that this is same as w_true, just that M() gives Mt now.
  TLorentzVector w_NeutrinoRec = dau1 + neutrinomet;
  
  TLorentzVector w_Rec = dau1+recmet;
  
  if(w_NeutrinoRec.M()<20)return;
    hCounter_->Fill(5);

  hWMass_->Fill(w_true.M());
  hWMt_->Fill(w_NeutrinoRec.M());
  hWRecMt_->Fill(w_Rec.M());
  hWPt_ ->Fill(w_true.Pt());
  hWRecPt_->Fill(w_Rec.Pt());
  hLepPt_->Fill(dau1.Pt());
  hLepEta_->Fill(dau1.Eta());
  
  edm::Handle<reco::GenJetCollection> genJetsHandle;
  if( not iEvent.getByLabel("ak5GenJets",genJetsHandle)){ 
    edm::LogInfo("GenAnalyzer") << "genJets not found, "
      "skipping event"; 
    return;
  }
    hCounter_->Fill(6);

  std::vector<TLorentzVector> selectedGenJets;
  
  const reco::GenJetCollection* genJetColl = &(*genJetsHandle);
  reco::GenJetCollection::const_iterator gjeti = genJetColl->begin();
  for(; gjeti!=genJetColl->end();gjeti++){
    reco::GenJet gjet = *gjeti;
    hpt_->Fill(gjet.pt());
    if(gjeti->pt()<30)continue;
    if(fabs(gjeti->eta())>2.4)continue;
    double deltaR = 9;
    bool keep = 1;
    for(unsigned int ii=0; ii!=stableLepColl.size(); ii++){
      
      double delta = radius(stableLepColl[ii]->momentum().eta(),
			    stableLepColl[ii]->momentum().phi(),
			    gjet.eta(),gjet.phi()
			    );
      if(delta<0.5){keep=0;break;}   
      
    }
    hdelR->Fill(deltaR);
    if (!keep)continue;
    hdelR2->Fill(deltaR);
    selectedGenJets.push_back(TLorentzVector(gjeti->px(), gjeti->py(), gjeti->pz(), gjeti->energy()));
  }
  std::sort(selectedGenJets.begin(),selectedGenJets.end(),PtGreater1());
  

  unsigned int njet =  selectedGenJets.size();
  //std::cout<<njet<<std::endl;   
  hNjetExclusive_->Fill(njet);
    hCounter_->Fill(7);

  if(njet<2)return;
  hCounter_->Fill(8);
  TLorentzVector j0= selectedGenJets[0];
  TLorentzVector j1= selectedGenJets[1];
  TLorentzVector dijet = selectedGenJets[0]+selectedGenJets[1];
  TLorentzVector dijetW = w_true+dijet;
  dijet_->Fill(j0,j1,dau1,met,w_true,w_Rec);  

  if(dijet.Pt()<40)
    dijetPtBelow40_->Fill(j0,j1,dau1,met,w_true,w_Rec);
  else
    dijetPtAbove40_->Fill(j0,j1,dau1,met,w_true,w_Rec);    
  
  
  if ((120<dijet.M())&&(160>dijet.M()))
    dijetMassWindow_->Fill(j0,j1,dau1,met,w_true,w_Rec);//W,lepPlusMet);
  else
    dijetSideBand_->Fill(j0,j1,dau1,met,w_true,w_Rec);//W,lepPlusMet);
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
MjjAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MjjAnalysis::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MjjAnalysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MjjAnalysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MjjAnalysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MjjAnalysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MjjAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MjjAnalysis);
