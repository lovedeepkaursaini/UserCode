// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include <boost/lexical_cast.hpp>

// ATLAS=0, CMS=1, CDF=2
#define EXPERIMENT 1


namespace Rivet {

#if EXPERIMENT==0
  class MC_LES_HOUCHES_SYSTEMATICS_ATLAS : public Analysis {
#elif EXPERIMENT==1
  class MC_LES_HOUCHES_SYSTEMATICS_CMS : public Analysis {
#elif EXPERIMENT==2
  class MC_LES_HOUCHES_SYSTEMATICS_CDF : public Analysis {
#endif
  public:

    /// Default constructor
#if EXPERIMENT==0
    MC_LES_HOUCHES_SYSTEMATICS_ATLAS() : Analysis("MC_LES_HOUCHES_SYSTEMATICS_ATLAS")
#elif EXPERIMENT==1
    MC_LES_HOUCHES_SYSTEMATICS_CMS() : Analysis("MC_LES_HOUCHES_SYSTEMATICS_CMS")
#elif EXPERIMENT==2
    MC_LES_HOUCHES_SYSTEMATICS_CDF() : Analysis("MC_LES_HOUCHES_SYSTEMATICS_CDF")
#endif
    {
      setNeedsCrossSection(true);
    }


    void init() {

#if EXPERIMENT==0
      VisibleFinalState fs(-5.0, 5.0, 0.*GeV);
      IdentifiedFinalState electrons(-5.0, 5.0, 20.*GeV);
      electrons.acceptIdPair(ELECTRON);
      IdentifiedFinalState muons(-5.0, 5.0, 20.*GeV);
      muons.acceptIdPair(MUON);
#elif EXPERIMENT==1
      VisibleFinalState fs(-3.0, 3.0, 0.*GeV);
      IdentifiedFinalState electrons(-3.0, 3.0, 20.*GeV);
      electrons.acceptIdPair(ELECTRON);
      IdentifiedFinalState muons(-3.0, 3.0, 20.*GeV);
      muons.acceptIdPair(MUON);
#elif EXPERIMENT==2
      VisibleFinalState fs(-3.0, 3.0, 0.*GeV);
      IdentifiedFinalState electrons(-1.0, 1.0, 20.*GeV);
      electrons.acceptIdPair(ELECTRON);
      IdentifiedFinalState muons(-1.0, 1.0, 20.*GeV);
      muons.acceptIdPair(MUON);
#endif
      MergedFinalState leptons(electrons, muons);

      addProjection(fs, "FS");

      VetoedFinalState vfs(fs);
      vfs.addVetoPairDetail(ELECTRON, 20.*GeV, 7000.*GeV);
      vfs.addVetoPairDetail(MUON,     20.*GeV, 7000.*GeV);

      VisibleFinalState missing(-10.0, 10.0, 0.*GeV);

      addProjection(electrons, "ELECTRONS");
      addProjection(muons, "MUONS");
      addProjection(leptons, "LEPTONS");
      addProjection(vfs, "VFS");
      addProjection(missing, "MISSING");

      _hnjets = bookHistogram1D("njet", 7, -0.5, 6.5);            // inclusive jet multiplicity (0..6)
      _hnjetsratio = bookDataPointSet("njetratio", 6, 0.5, 6.5);  // n/(n-1)

      for (int i=0 ; i<6 ; i++) {
        _hjetpt.push_back(bookHistogram1D("jetpt"+boost::lexical_cast<string>(i), 50, 0, 250));  // jet pT  (1..6)
        _hjeteta.push_back(bookHistogram1D("jeteta"+boost::lexical_cast<string>(i), 16, -4, 4));  // jet eta (1..6)
      }

      for (int i=0 ; i<7 ; i++) {
        _hHTjet.push_back(bookHistogram1D("HTjet"+boost::lexical_cast<string>(i), 50, 0, 1000));  // HT from jets
        _hHTall.push_back(bookHistogram1D("HTall"+boost::lexical_cast<string>(i), 50, 0, 1000));  // HT from jets + lepton + missing
        _hsumET.push_back(bookHistogram1D("sumET"+boost::lexical_cast<string>(i), 50, 0, 1000));  // sum ET of visible particles
      }

      _hdEtaj0j1    = bookHistogram1D("dEtaj0j1", 20, 0, 5);     // deltaEta(leading jets)
      _hdPhij0j1    = bookHistogram1D("dPhij0j1", 20, 0, M_PI);  // deltaPhi(leading jets)
      _hdRj0j1      = bookHistogram1D("dRj0j1", 20, 0, 5);       // deltaR(leading jets)
      _hmj0j1       = bookHistogram1D("mj0j1", 60, 0, 300);      // mass(jet0 + jet1)
      _hptratioj1j0 = bookHistogram1D("ptratioj1j0", 20, 0, 1);  // pT(jet1)/pT(jet0)

      _hdEtaj0l = bookHistogram1D("dEtaj0l", 20, 0, 5);          // deltaEta(leading jet, lepton)
      _hdPhij0l = bookHistogram1D("dPhij0l", 20, 0, M_PI);       // deltaPhi(leading jet, lepton)
      _hdRj0l   = bookHistogram1D("dRj0l", 20, 0, 5);            // deltaR(leading jet, lepton)
      _hmj0l    = bookHistogram1D("mj0l", 60, 0, 300);           // mass(jet0 + lepton)

      _hmj0j1W = bookHistogram1D("mj0j1W", 50, 0, 500);          // mass(jet0 + jet1 + W)

      _hbeamthrustjets      = bookHistogram1D("beamthrustjets", 50, 0, 100);      // Whatever-1.
      _hbeamthrustparticles = bookHistogram1D("beamthrustparticles", 50, 0, 100); // Whatever-2.

      _hleptonpt  = bookHistogram1D("leptonpt", 50, 0, 200);     // lepton pT
      _hleptoneta = bookHistogram1D("leptoneta", 12, -3, 3);     // lepton eta

      _hWpt   = bookHistogram1D("Wpt", 50, 0, 200);              // W pT
      _hWeta  = bookHistogram1D("Weta", 20, -5, 5);              // W eta
      _hWmass = bookHistogram1D("Wmass", 50, 0, 500);           // W mass
      _hWInvmass = bookHistogram1D("WInvmass", 50, 0, 300);           // W mass
      _hWmt   = bookHistogram1D("Wmt", 40, 0, 200);              // W transverse mass

      _hsigmatot = bookDataPointSet("sigmatot", 1, 0, 1);        // sigma_tot as reported by the generator
      _hsigmacut = bookHistogram1D("sigmacut", 6, -0.5, 5.5);    // sigma after each cut 
    }


    void analyze(const Event & event) {
      const double weight = event.weight();

      _hsigmacut->fill(0, weight);

      const FinalState& allleptonsfs = applyProjection<FinalState>(event, "LEPTONS");
      ParticleVector allleptons = allleptonsfs.particlesByPt();
      if (allleptons.size() < 1) vetoEvent;

      _hsigmacut->fill(1, weight);

      // Isolation cut
      Particle lepton;
      const FinalState& fullfs = applyProjection<FinalState>(event, "MISSING");
      bool found_lepton = false;
      for (size_t i=0 ; i<allleptons.size() ; i++) {
        FourMomentum testmom = allleptons[i].momentum();
#if EXPERIMENT==0
        if (fabs(testmom.eta())>2.5) continue;
#elif EXPERIMENT==1
        if ((abs(allleptons[i].pdgId())==MUON && fabs(testmom.eta())>2.1) ||
            (abs(allleptons[i].pdgId())==ELECTRON && fabs(testmom.eta())>2.5)) continue;
#elif EXPERIMENT==2
        if (fabs(testmom.eta())>1.0) continue;
#endif
        double etsum(-testmom.Et());
        foreach (Particle hit, fullfs.particles()) {
          FourMomentum trackmom = hit.momentum();
          if (deltaR(testmom,trackmom)<0.5) {
            etsum += trackmom.Et();
            if (etsum>0.1*testmom.Et())
              break;
          }
        }
        if (etsum<0.1*testmom.Et()) {
          lepton = allleptons[i];
          allleptons.erase(allleptons.begin()+i);
          found_lepton = true;
          break;
        }
      }
      if (!found_lepton) vetoEvent;
      _hsigmacut->fill(2, weight);


      // Missing ET cut
      FourMomentum missingmom;
      foreach (Particle hit, fullfs.particles()) {
        missingmom += hit.momentum();
      }
      missingmom *= -1; // missing is "minus visible"
      missingmom.setE(missingmom.vector3().mod()); // assume neutrinos are massless
#if EXPERIMENT==0
      if (missingmom.Et()<25.*GeV) vetoEvent;
#elif EXPERIMENT==2
      if (missingmom.Et()<25.*GeV) vetoEvent;
#endif
      _hsigmacut->fill(3, weight);



      // Create a W
      FourMomentum Wmom = missingmom + lepton.momentum();



      // Transverse mass cut
      double mT2 = pow(lepton.momentum().pT()+missingmom.pT(),2)-Wmom.pT2();
#if EXPERIMENT==0
      if (sqrt(mT2) < 40.*GeV) vetoEvent;
#elif EXPERIMENT==1
      if (sqrt(mT2) < 20.*GeV) vetoEvent;
#elif EXPERIMENT==2
      if (sqrt(mT2) < 30.*GeV) vetoEvent;
#endif
      _hsigmacut->fill(4, weight);


      // Reconstruct jets
      const FinalState& vfs = applyProjection<FinalState>(event, "VFS");
#if EXPERIMENT==0
      FastJets jetsproj(vfs, FastJets::ANTIKT, 0.4);
      jetsproj.calc(vfs.particles()+allleptons);
      Jets alljets = jetsproj.jetsByPt(25.*GeV);
#elif EXPERIMENT==1
      FastJets jetsproj(vfs, FastJets::ANTIKT, 0.5);
      jetsproj.calc(vfs.particles()+allleptons);
      Jets alljets = jetsproj.jetsByPt(30.*GeV);
#elif EXPERIMENT==2
      FastJets jetsproj(vfs, FastJets::CDFJETCLU, 0.4);
      jetsproj.calc(vfs.particles()+allleptons);
      Jets alljets = jetsproj.jetsByPt(30.*GeV);
#endif
      Jets jets;
      foreach (Jet jet, alljets) {
#if EXPERIMENT==0
        if (fabs(jet.momentum().eta())<4.4)
#elif EXPERIMENT==1
        if (fabs(jet.momentum().eta())<2.4)
#elif EXPERIMENT==2
        if (fabs(jet.momentum().eta())<2.4)
#endif
          jets.push_back(jet);
      }


#if EXPERIMENT==2
      if (jets.size()<2 || (jets[0].momentum()+jets[1].momentum()).pT()<40.*GeV) vetoEvent;
#endif
      _hsigmacut->fill(5, weight);


      // Fill the histograms
      _hWpt->fill(Wmom.pT(),weight);
      _hWeta->fill(Wmom.eta(),weight);
      _hWInvmass->fill(sqrt(Wmom.invariant()),weight);
      _hWmass->fill(Wmom.mass(),weight);
      _hWmt->fill(sqrt(mT2),weight);
//cout<<"W mass"<<Wmom.mass()<<"  inv. "<<sqrt(Wmom.invariant())<<endl;
      _hleptonpt->fill(lepton.momentum().pT(),weight);
      _hleptoneta->fill(lepton.momentum().eta(),weight);

      double HTjet = 0.;
      double HTall = 0.;
      double sumET = 0.;
      double beamthrustjets = 0.;
      double beamthrustparticles = 0.;
      foreach (Jet jet, jets) {
        HTjet += jet.momentum().Et();
        HTall += jet.momentum().Et();
        beamthrustjets += jet.momentum().E() - fabs(jet.momentum().z());
      }
      HTall += Wmom.Et();

      foreach (Particle p, vfs.particles()+allleptons) {
        sumET += p.momentum().Et();
        beamthrustparticles += p.momentum().E() - fabs(p.momentum().z());
      }

      _hbeamthrustjets->fill(beamthrustjets, weight);
      _hbeamthrustparticles->fill(beamthrustparticles, weight);
      _hnjets->fill(0, weight);  // I guess we always have at least 0 jets
      _hHTjet[0]->fill(HTjet, weight);
      _hHTall[0]->fill(HTall, weight);
      _hsumET[0]->fill(sumET, weight);
      for (size_t i=0 ; i<jets.size() ; i++) {
        if (i==6) break;
        _hnjets->fill(i+1, weight);  // njets is inclusive
        _hjetpt[i]->fill(jets[i].momentum().pT(), weight);
        _hjeteta[i]->fill(jets[i].momentum().eta(), weight);
        _hHTjet[i+1]->fill(HTjet, weight);
        _hHTall[i+1]->fill(HTall, weight);
        _hsumET[i+1]->fill(sumET, weight);
      }
      if (jets.size() >= 1) {
        _hdEtaj0l->fill(deltaEta(jets[0],lepton), weight);
        _hdPhij0l->fill(deltaPhi(jets[0],lepton),weight);
        _hdRj0l->fill(deltaR(jets[0],lepton),weight);
        _hmj0l->fill((jets[0].momentum() + lepton.momentum()).mass(),weight);
      }

      if (jets.size() >= 2) {
        _hdEtaj0j1->fill(deltaEta(jets[0],jets[1]),weight);
        _hdPhij0j1->fill(deltaPhi(jets[0],jets[1]),weight);
        _hdRj0j1->fill(deltaR(jets[0],jets[1]),weight);
        _hmj0j1->fill((jets[0].momentum() + jets[1].momentum()).mass(),weight);
        _hmj0j1W->fill((jets[0].momentum() + jets[1].momentum() + Wmom).mass(),weight);
        _hptratioj1j0->fill(jets[1].momentum().pT()/jets[0].momentum().pT(),weight);
      }

    }

    /// Finalize
    void finalize() {
      AIDA::IHistogramFactory& hf = histogramFactory();
      const string dir = histoDir();

      for (int i=0 ; i<6 ; i++) {
        hf.divide(dir + "/HTjet" + boost::lexical_cast<string>(i+1) + "over" + boost::lexical_cast<string>(i), *_hHTjet[i+1], *_hHTjet[i]);
        hf.divide(dir + "/HTall" + boost::lexical_cast<string>(i+1) + "over" + boost::lexical_cast<string>(i), *_hHTall[i+1], *_hHTall[i]);
        hf.divide(dir + "/sumET" + boost::lexical_cast<string>(i+1) + "over" + boost::lexical_cast<string>(i), *_hsumET[i+1], *_hsumET[i]);
      }

      std::vector<double> y, yerr;
      for (int i=0; i<_hnjets->axis().bins()-1; i++) {
        double val = 0.;
        double err = 0.;
        if (!fuzzyEquals(_hnjets->binHeight(i), 0)) {
          val = _hnjets->binHeight(i+1) / _hnjets->binHeight(i);
          err = val * sqrt(  pow(_hnjets->binError(i+1)/_hnjets->binHeight(i+1), 2)
                           + pow(_hnjets->binError(i)  /_hnjets->binHeight(i)  , 2) );
        }
        y.push_back(val);
        yerr.push_back(err);
      }
      _hnjetsratio->setCoordinate(1, y, yerr);

      std::vector<double> sigma, sigmaerr;
      sigma.push_back(crossSection());
      sigmaerr.push_back(0);
      _hsigmatot->setCoordinate(1, sigma, sigmaerr);
      //      scale(_hnjets, crossSection()/sumOfWeights());  original style to balance
      double csection = 1; 
      //double csection = crossSection();
      scale(_hnjets, csection/sumOfWeights());
      scale(_hdEtaj0j1, csection/sumOfWeights());
      scale(_hdPhij0j1, csection/sumOfWeights());
      scale(_hdRj0j1, csection/sumOfWeights());
      scale(_hmj0j1, csection/sumOfWeights());
      scale(_hptratioj1j0, csection/sumOfWeights());
      scale(_hdEtaj0l, csection/sumOfWeights());
      scale(_hdPhij0l, csection/sumOfWeights());
      scale(_hdRj0l, csection/sumOfWeights());
      scale(_hmj0l, csection/sumOfWeights());
      scale(_hmj0j1W, csection/sumOfWeights());
      scale(_hbeamthrustjets, csection/sumOfWeights());
      scale(_hbeamthrustparticles, csection/sumOfWeights());
      scale(_hleptonpt, csection/sumOfWeights());
      scale(_hleptoneta, csection/sumOfWeights());
      scale(_hWpt, csection/sumOfWeights());
      scale(_hWeta, csection/sumOfWeights());
      scale(_hWmass, csection/sumOfWeights());
      scale(_hWInvmass, csection/sumOfWeights());
      scale(_hWmt, csection/sumOfWeights());
      scale(_hsigmacut, csection/sumOfWeights());

      for (int i=0 ; i<6 ; i++) {
        scale(_hjetpt[i], csection/sumOfWeights());
        scale(_hjeteta[i], csection/sumOfWeights());
      }

      for (int i=0 ; i<7 ; i++) {
        scale(_hHTjet[i], csection/sumOfWeights());
        scale(_hHTall[i], csection/sumOfWeights());
        scale(_hsumET[i], csection/sumOfWeights());
      }
    }

  private:

    AIDA::IHistogram1D * _hnjets;               // inclusive jet multiplicity (0..6)
    AIDA::IDataPointSet* _hnjetsratio;          // n/(n-1)

    std::vector<AIDA::IHistogram1D*> _hjetpt;   // jet pT  (1..6)
    std::vector<AIDA::IHistogram1D*> _hjeteta;  // jet eta (1..6)

    std::vector<AIDA::IHistogram1D*> _hHTjet;   // HT from jets
    std::vector<AIDA::IHistogram1D*> _hHTall;   // HT from jets + lepton + missing
    std::vector<AIDA::IHistogram1D*> _hsumET;   // sum ET of visible particles

    AIDA::IHistogram1D * _hdEtaj0j1;            // deltaEta(leading jets)
    AIDA::IHistogram1D * _hdPhij0j1;            // deltaPhi(leading jets)
    AIDA::IHistogram1D * _hdRj0j1;              // deltaR(leading jets)
    AIDA::IHistogram1D * _hmj0j1;               // mass(jet0 + jet1)
    AIDA::IHistogram1D * _hptratioj1j0;         // pT(jet1)/pT(jet0)

    AIDA::IHistogram1D * _hdEtaj0l;             // deltaEta(leading jet, lepton)
    AIDA::IHistogram1D * _hdPhij0l;             // deltaPhi(leading jet, lepton)
    AIDA::IHistogram1D * _hdRj0l;               // deltaR(leading jet, lepton)
    AIDA::IHistogram1D * _hmj0l;                // mass(jet0 + lepton)

    AIDA::IHistogram1D * _hmj0j1W;              // mass(jet0 + jet1 + W)

    AIDA::IHistogram1D * _hbeamthrustjets;      // Whatever-1.
    AIDA::IHistogram1D * _hbeamthrustparticles; // Whatever-2.

    AIDA::IHistogram1D * _hleptonpt;            // lepton pT
    AIDA::IHistogram1D * _hleptoneta;           // lepton eta

    AIDA::IHistogram1D * _hWpt;                 // W pT
    AIDA::IHistogram1D * _hWeta;                // W eta
    AIDA::IHistogram1D * _hWmass;               // W mass
    AIDA::IHistogram1D * _hWInvmass;               // W mass
    AIDA::IHistogram1D * _hWmt;                 // W tranverse mass

    AIDA::IDataPointSet* _hsigmatot;            // sigma_tot as reported by the generator
    AIDA::IHistogram1D * _hsigmacut;            // sigma after each cut
  };

  // This global object acts as a hook for the plugin system
#if EXPERIMENT==0
  AnalysisBuilder<MC_LES_HOUCHES_SYSTEMATICS_ATLAS> plugin_MC_LES_HOUCHES_SYSTEMATICS_ATLAS;
#elif EXPERIMENT==1
  AnalysisBuilder<MC_LES_HOUCHES_SYSTEMATICS_CMS> plugin_MC_LES_HOUCHES_SYSTEMATICS_CMS;
#elif EXPERIMENT==2
  AnalysisBuilder<MC_LES_HOUCHES_SYSTEMATICS_CDF> plugin_MC_LES_HOUCHES_SYSTEMATICS_CDF;
#endif
}

