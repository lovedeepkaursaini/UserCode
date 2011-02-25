#include "UserCode/TrigTree/interface/genInfoTree.h"

//---------------------------------------------------------------
// Add Branches to the genTree
//---------------------------------------------------------------
void
genInfoTree::AddBranch(double* x, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),x,(brName+"/D").c_str());

}


//---------------------------------------------------------------
//---------------------------------------------------------------
void
genInfoTree::AddBranch(std::vector<double>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);


}

//---------------------------------------------------------------
//---------------------------------------------------------------
void
genInfoTree::AddBranch(std::vector<int>* vec, std::string name){
  std::string brName=name;
  tree_->Branch(brName.c_str(),vec);
}

genInfoTree::genInfoTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig)
{
  tree_=tree; 
  genPartLabel_ = iConfig.getParameter<edm::InputTag>("genPartLabel");
  SetBranches();
}


genInfoTree::~genInfoTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete tree_;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
genInfoTree::Fill(const edm::Event& iEvent)
{
  Clear();
   using namespace edm;
   edm::Handle<reco::GenParticleCollection> genParticleHandle;
   if(not iEvent.getByLabel(genPartLabel_, genParticleHandle))
     {
       std::cout<<
	 "GenAnalyzer: Generator Level Information not found\n"
		<<std::endl;
     }


   const reco::GenParticleCollection* genColl= &(*genParticleHandle);
   reco::GenParticleCollection::const_iterator geni = genColl->begin();
   for(; geni!=genColl->end();geni++){
     reco::GenParticle gen = *geni;
     
     //Look out for the GenMuons
     if((abs(gen.pdgId())==13 && gen.status()==1)){
       genMuPx_.push_back(gen.px());
       genMuPy_.push_back(gen.py());
       genMuPz_.push_back(gen.pz());
       genMuE_.push_back(gen.energy());
       genMuP_.push_back(gen.p());
       genMuPt_.push_back(gen.pt());
       genMuEta_.push_back(gen.eta());
       genMuPhi_.push_back(gen.phi());
       genMuTheta_.push_back(gen.theta());
       genMuEt_.push_back(gen.et());
       genMuQ_.push_back(gen.charge());
       
       //       std::cout<<gen.pt()<<std::endl;
       //       std::cout<<gen.eta()<<std::endl;

     }
 

     //Look out for the GenElectrons
     if((abs(gen.pdgId())==11 && gen.status()==1)){
       genElePx_.push_back(gen.px());
       genElePy_.push_back(gen.py());
       genElePz_.push_back(gen.pz());
       genEleE_.push_back(gen.energy());
       genEleP_.push_back(gen.p());
       genElePt_.push_back(gen.pt());
       genEleEta_.push_back(gen.eta());
       genElePhi_.push_back(gen.phi());
       genEleTheta_.push_back(gen.theta());
       genEleEt_.push_back(gen.et());
       genEleQ_.push_back(gen.charge());
       //genMuChar_.push_back(gen.)
     }
   


   }

  

   edm::Handle<reco::GenJetCollection> genJetsHandle;
   if( not iEvent.getByLabel("iterativeCone5GenJets",genJetsHandle)){ 
     edm::LogInfo("GenAnalyzer") << "genJets not found, "
       "skipping event"; 
     return;
   }
   
   
   const reco::GenJetCollection* genJetColl = &(*genJetsHandle);
   reco::GenJetCollection::const_iterator gjeti = genJetColl->begin();
   

   for(; gjeti!=genJetColl->end();gjeti++){
     reco::GenParticle gjet = *gjeti;
     genJetPx_.push_back(gjet.px());
     genJetPy_.push_back(gjet.py()); 
     genJetPz_.push_back(gjet.pz()); 
     genJetE_.push_back(gjet.energy()); 
     genJetP_.push_back(gjet.p()); 
     genJetPt_.push_back(gjet.pt());
     genJetEta_.push_back(gjet.eta());
     genJetPhi_.push_back(gjet.phi());
     genJetTheta_.push_back(gjet.theta());
     genJetEt_.push_back(gjet.et());
     //genJetQ_.push_back(gjet.charge());
     //std::cout<<gjet.pt()<<std::endl;

   }
}



void  
genInfoTree::SetBranches(){
  AddBranch(&genMuPx_, "genMuPx");
  AddBranch(&genMuPy_,"genMuPy");
  AddBranch(&genMuPz_,"genMuPz");
  AddBranch(&genMuE_, "genMuE");
  AddBranch(&genMuP_,"genMuP");
  AddBranch(&genMuTheta_,"genMuTheta");
  AddBranch(&genMuPt_, "genMuPt");
  AddBranch(&genMuEta_,"genMuEta");
  AddBranch(&genMuPhi_,"genMuPhi");
  AddBranch(&genMuEt_,"genMuEt");
  AddBranch(&genMuQ_,"genMuQ");

  
  AddBranch(&genElePx_, "genElePx");
  AddBranch(&genElePy_,"genElePy");
  AddBranch(&genElePz_,"genElePz");
  AddBranch(&genEleE_, "genEleE");
  AddBranch(&genEleP_,"genEleP");
  AddBranch(&genEleTheta_,"genEleTheta");
  AddBranch(&genElePt_, "genElePt");
  AddBranch(&genEleEta_,"genEleEta");
  AddBranch(&genElePhi_,"genElePhi");
  AddBranch(&genEleEt_,"genEleEt"); 
  AddBranch(&genEleQ_,"genEleQ");

  AddBranch(&genJetPx_, "genJetPx");
  AddBranch(&genJetPy_,"genJetPy");
  AddBranch(&genJetPz_,"genJetPz");
  AddBranch(&genJetE_, "genJetE");
  AddBranch(&genJetP_,"genJetP");
  AddBranch(&genJetTheta_,"genJetTheta");
  AddBranch(&genJetPt_,"genJetPt");
  AddBranch(&genJetEta_,"genJetEta");
  AddBranch(&genJetPhi_,"genJetPhi");
  AddBranch(&genJetEt_,"genJetEt"); 
  AddBranch(&genJetQ_,"genJetQ");

}


void  
genInfoTree::Clear(){
  genMuPx_.clear();
  genMuPy_.clear();
  genMuPz_.clear();  
  genMuE_.clear();
  genMuP_.clear();
  genMuPt_.clear();
  genMuEta_.clear();
  genMuPhi_.clear();  
  genMuEt_.clear();
  genMuQ_.clear();
  genMuTheta_.clear();

  genElePx_.clear();
  genElePy_.clear();
  genElePz_.clear();  
  genEleE_.clear();
  genEleP_.clear();
  genElePt_.clear();
  genEleEta_.clear();
  genElePhi_.clear();
  genEleEt_.clear();
  genEleQ_.clear();
  genEleTheta_.clear();

  genJetPx_.clear();
  genJetPy_.clear();
  genJetPz_.clear();  
  genJetE_.clear();
  genJetP_.clear();
  genJetPt_.clear();
  genJetEta_.clear();
  genJetPhi_.clear();
  genJetEt_.clear();
  genJetQ_.clear();
  genJetTheta_.clear();


  //  genTQuark_.clear();
}
