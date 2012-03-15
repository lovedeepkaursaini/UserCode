#include "DelPanj/TreeMaker/interface/patZJetPlant.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "FWCore/Framework/interface/ESHandle.h"


#include "TMath.h"

typedef math::XYZTLorentzVector LorentzVector;
bool DescendingOrder(const TLorentzVector& l1, const TLorentzVector& l2){                                                                              
    return l1.Pt()>l2.Pt();                                                                                                                                  
  }

bool AscendingOrder(const TLorentzVector& l1, const TLorentzVector& l2){
    return l1.Pt()<l2.Pt();
  }

patZJetPlant::patZJetPlant(
			   std::string desc, TTree* tree, 
			   const edm::ParameterSet& iConfig
			   ):
  baseTree(desc, tree),
  patElecLabel_(iConfig.getParameter<edm::InputTag>("patElectrons")),
  patJetLabel_(iConfig.getParameter<edm::InputTag>("patJetLabel")),
  ewp1_(iConfig.getParameter<edm::ParameterSet>("leadElecPset_")),
  ewp2_(iConfig.getParameter<edm::ParameterSet>("subLeadElecPset_"))
{
  SetBranches();
  Clear();
}



patZJetPlant::~patZJetPlant(){
  delete tree_;
}

void
patZJetPlant::Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup){
  Clear();
  //std::cout<<"Debugger: "<<std::endl;   
  //MC Block.
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  if(not iEvent.getByLabel("genParticles", genParticleHandle))
    std::cout<< "GenAnalyzer: Generator Level Information not found\n" <<std::endl;

  //const reco::GenParticleCollection* genColl= &(*genParticleHandle);
  //reco::GenParticleCollection::const_iterator geni = genColl->begin();
  std::vector<TLorentzVector> electrons ;
  const reco::GenParticleCollection* genColl= &(*genParticleHandle);
  reco::GenParticleCollection::const_iterator geni = genColl->begin();
  for(; geni!=genColl->end();geni++){
    reco::GenParticle gen = *geni;
    //Look out for the GenMuons/Electrons
    if(gen.status()==1 && fabs(gen.pdgId())==11){
      TLorentzVector* ele = new TLorentzVector();
      ele->SetPxPyPzE(gen.px(),gen.py(),gen.pz(),gen.energy());
      electrons.push_back(*ele);
    }
  }
  
  //IF NO: 1 (IF YOU HAVE TWO GEN ELECTRONS)
  if(electrons.size()>1){
    std::sort(electrons.begin(),electrons.end(),DescendingOrder);
    TLorentzVector d1(electrons[0]);
    TLorentzVector d2(electrons[1]);
    TLorentzVector Z(d1+d2);
       pxGenZ = Z.Px();
       pyGenZ = Z.Py();
       pzGenZ = Z.Pz();
       enGenZ = Z.E();
       pxD1GenZ = d1.Px();
       pyD1GenZ = d1.Py();
       pzD1GenZ = d1.Pz();
       enD1GenZ = d1.E();
       pxD2GenZ = d2.Px();
       pyD2GenZ = d2.Py();
       pzD2GenZ = d2.Pz();
       enD2GenZ = d2.E();
  }
  
  
  
  edm::Handle<pat::ElectronCollection> patElecHandle;
  if(not iEvent.getByLabel(patElecLabel_,patElecHandle))
    std::cout<<"Following Not Found: "<<patElecLabel_<<std::endl;
  pat::ElectronCollection eleColl(*(patElecHandle.product()));
  
  //IF NO : 2 (IF YOU HAVE TWO PAT ELECTRONS)
  if(eleColl.size()>1){
    pat::Electron e1 = eleColl[0];
    pat::Electron e2 = eleColl[1];
    std::map<std::string, bool> leadPass = ewp1_.CutRecord(e1);
    std::map<std::string, bool> subLeadPass = ewp2_.CutRecord(e2);
    TLorentzVector l1 = Part2LorVec(e1);
    TLorentzVector l2 = Part2LorVec(e2);
    TLorentzVector patZ(l1+l2);
    pxPatZ = patZ.Px();
    pyPatZ = patZ.Py();
    pzPatZ = patZ.Pz();
    enPatZ = patZ.E();
    pxD1PatZ = l1.Px();
    pyD1PatZ = l1.Py();
    pzD1PatZ = l1.Pz();
    enD1PatZ = l1.E();
    pxD2PatZ = l2.Px();
    pyD2PatZ = l2.Py();
    pzD2PatZ = l2.Pz();
    enD2PatZ = l2.E();
    passD1PatZ = (int)PassAll(leadPass);
    passD2PatZ = (int)PassAll(subLeadPass);
    //std::cout<<passD1PatZ<<"\t"<<passD2PatZ<<std::endl;
  }
   //if(pxD1GenZ!=-9999)
   //std::cout<<"NEle: pat, gen"<<eleColl.size()<<"\t"<<electrons.size()<<std::endl;
   bool verbose = 0;
   if(verbose)
   if(eleColl.size()>1 && electrons.size()>1)
   if((pxD1GenZ!=-9999)&&(pxD2GenZ!=-9999)){
   std::cout<<"============================================\n";
   std::cout<<"\tGenZDau1: "<<pxD1GenZ<<"\t"<<pyD1GenZ<<"\t"<<pzD1GenZ<<"\t"<<enD1GenZ<<"\n";
   std::cout<<"\tGenZDau2: "<<pxD2GenZ<<"\t"<<pyD2GenZ<<"\t"<<pzD2GenZ<<"\t"<<enD2GenZ<<"\n";
   std::cout<<"\tPatZDau1: "<<pxD1PatZ<<"\t"<<pyD1PatZ<<"\t"<<pzD1PatZ<<"\t"<<enD1PatZ<<"\n";
   std::cout<<"\tPatZDau2: "<<pxD2PatZ<<"\t"<<pyD2PatZ<<"\t"<<pzD2PatZ<<"\t"<<enD2PatZ<<"\n";
   std::cout<<"============================================\n";
   }
 
  //Fill The GenJet pt>20 GeV
  edm::Handle<reco::GenJetCollection> genJetsHandle;
  if( not iEvent.getByLabel("ak5GenJets",genJetsHandle)){
    edm::LogInfo("GenAnalyzer") << "genJets not found, "
      "skipping event";
    return;
  }

  const reco::GenJetCollection* genJetColl = &(*genJetsHandle);
  reco::GenJetCollection::const_iterator gjeti = genJetColl->begin();
  for(; gjeti!=genJetColl->end();gjeti++){
    reco::GenParticle gjet = *gjeti;
    if(gjet.pt()<=20)continue;
    if(fabs(gjet.eta())> 3)continue;
    genJetPx.push_back(gjet.px());
    genJetPy.push_back(gjet.py());
    genJetPz.push_back(gjet.pz());
    genJetEn.push_back(gjet.energy());
      }
  
  //Fill The RecoJets pt>20 GeV
  edm::Handle<std::vector<pat::Jet> > patJetHandle;
  if(not iEvent.getByLabel("selectedPatJetsPFlow",patJetHandle))
    std::cout<<"PAT Jets Not Found" <<std::endl;
  const std::vector<pat::Jet>*  patJets = patJetHandle.product();
  std::vector<pat::Jet>::const_iterator pjet =patJets->begin();
  for(;pjet!=patJets->end();pjet++){
    if(pjet->pt()<20)continue;
    if(abs(pjet->eta())>3)continue;
    patJetPx.push_back(pjet->px());
    patJetPy.push_back(pjet->py());
    patJetPz.push_back(pjet->pz());
    patJetEn.push_back(pjet->energy());
  }
}


void
patZJetPlant::SetBranches(){
  AddBranch(&pxGenZ,  "pxGenZ");  
  AddBranch(&pyGenZ  ,"pyGenZ");
  AddBranch(&pzGenZ  ,"pzGenZ");
  AddBranch(&enGenZ  ,"enGenZ");
  AddBranch(&pxD1GenZ  ,"pxD1GenZ");
  AddBranch(&pyD1GenZ  ,"pyD1GenZ");
  AddBranch(&pzD1GenZ  ,"pzD1GenZ");
  AddBranch(&enD1GenZ  ,"enD1GenZ");
  AddBranch(&pxD2GenZ  ,"pxD2GenZ");
  AddBranch(&pyD2GenZ  ,"pyD2GenZ");
  AddBranch(&pzD2GenZ  ,"pzD2GenZ");
  AddBranch(&enD2GenZ  ,"enD2GenZ");
  AddBranch(&pxPatZ  ,"pxPatZ");
  AddBranch(&pyPatZ  ,"pyPatZ");
  AddBranch(&pzPatZ  ,"pzPatZ");
  AddBranch(&enPatZ  ,"enPatZ");
  AddBranch(&pxD1PatZ  ,"pxD1PatZ");
  AddBranch(&pyD1PatZ  ,"pyD1PatZ");
  AddBranch(&pzD1PatZ  ,"pzD1PatZ");
  AddBranch(&enD1PatZ  ,"enD1PatZ");
  AddBranch(&passD1PatZ  ,"passD1PatZ");
  AddBranch(&pxD2PatZ  ,"pxD2PatZ");
  AddBranch(&pyD2PatZ  ,"pyD2PatZ");
  AddBranch(&pzD2PatZ  ,"pzD2PatZ");
  AddBranch(&enD2PatZ  ,"enD2PatZ");
  AddBranch(&passD2PatZ,  "passD2PatZ");  
  AddBranch(&genJetPx,"genJetPx");
  AddBranch(&genJetPy,"genJetPy");
  AddBranch(&genJetPz,"genJetPz");
  AddBranch(&genJetEn,"genJetEn");
  AddBranch(&patJetPx,"patJetPx");
  AddBranch(&patJetPy,"patJetPy");
  AddBranch(&patJetPz,"patJetPz");
  AddBranch(&patJetEn,"patJetEn");
}


void
patZJetPlant::Clear(){
  pxGenZ = -9999;
  pyGenZ = -9999;
  pzGenZ = -9999;
  enGenZ = -9999;
  pxD1GenZ = -9999;
  pyD1GenZ = -9999;
  pzD1GenZ = -9999;
  enD1GenZ = -9999;
  pxD2GenZ = -9999;
  pyD2GenZ = -9999;
  pzD2GenZ = -9999;
  enD2GenZ = -9999;
  pxPatZ = -9999;
  pyPatZ = -9999;
  pzPatZ = -9999;
  enPatZ = -9999;
  pxD1PatZ = -9999;
  pyD1PatZ = -9999;
  pzD1PatZ = -9999;
  enD1PatZ = -9999;
  passD1PatZ = 0;
  pxD2PatZ = -9999;
  pyD2PatZ = -9999;
  pzD2PatZ = -9999;
  enD2PatZ = -9999;
  passD2PatZ = 0;
  genJetPx.clear();
  genJetPy.clear();
  genJetPz.clear();
  genJetEn.clear();
  patJetPx.clear();
  patJetPy.clear();
  patJetPz.clear();
  patJetEn.clear();
}
