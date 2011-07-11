// -*- C++ -*-
// AUTHOR:  Anil Singh (anil@cern.ch), Lovedeep Saini (lovedeep@cern.ch)

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"


namespace Rivet {
  
  
  class CMS_EWK_10_012 : public Analysis {
  public:
    
    CMS_EWK_10_012()
      : Analysis("CMS_EWK_10_012")
    {
      setBeams(PROTON, PROTON);
      setNeedsCrossSection(true);
    }
    
    
    /// Book histograms and initialise projections before the run
    void init() {
      
      const FinalState fs(-MAXRAPIDITY,MAXRAPIDITY);
      addProjection(fs, "FS");
      
      vector<pair<PdgId,PdgId> > vidsZ;
      vidsZ.push_back(make_pair(ELECTRON, POSITRON));
      vidsZ.push_back(make_pair(MUON, ANTIMUON));

      FinalState fsZ(-MAXRAPIDITY,MAXRAPIDITY);
      InvMassFinalState invfsZ(fsZ, vidsZ, 60*GeV, 120*GeV);
      addProjection(invfsZ, "INVFSZ");
      
      vector<pair<PdgId,PdgId> > vidsW;
      vidsW.push_back(make_pair(ELECTRON, NU_EBAR));
      vidsW.push_back(make_pair(POSITRON, NU_E));
      vidsW.push_back(make_pair(MUON, NU_MUBAR));
      vidsW.push_back(make_pair(ANTIMUON, NU_MU));
      
      FinalState fsW(-MAXRAPIDITY,MAXRAPIDITY);
      InvMassFinalState invfsW(fsW, vidsW, 20*GeV, 99999*GeV);
      addProjection(invfsW, "INVFSW");
      
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(invfsZ);
      vfs.addVetoOnThisFinalState(invfsW);
      addProjection(vfs, "VFS");
      addProjection(FastJets(vfs, FastJets::ANTIKT, 0.5), "Jets");
      
      for (int i = 0; i < 4; ++i) {
	_histJetMultNoverNm1Welec[i]   = bookDataPointSet(2,i+1,i+1);
	_histJetMultNoverN0Welec[i] = bookDataPointSet(1, i+1, i+1);
        _histJetMultNoverNm1Wmu[i]   = bookDataPointSet(4,i+1,i+1);
        _histJetMultNoverN0Wmu[i] = bookDataPointSet(3, i+1, i+1);
        _histJetMultRatioWmuPlusMinus[i] = bookDataPointSet(10, i+1, i+1);
	_histJetMultRatioWelPlusMinus[i] = bookDataPointSet(9, i+1, i+1);
	_histJetMultNoverNm1Zelec[i]   = bookDataPointSet(6,i+1,i+1);
	_histJetMultNoverN0Zelec[i] = bookDataPointSet(5, i+1, i+1);
	_histJetMultNoverNm1Zmu[i]     = bookDataPointSet(8,i+1,i+1);
	_histJetMultNoverN0Zmu[i]   = bookDataPointSet(7, i+1, i+1);
      }

      _histJetMultWelec  = bookHistogram1D("njetWenu", 7, -0.5, 6.5);
      _histJetMultWmu    = bookHistogram1D("njetWmunu", 7, -0.5, 6.5);
      _histJetMultZelec  = bookHistogram1D("njetZee", 7, -0.5, 6.5);
      _histJetMultZmu    = bookHistogram1D("njetZmumu", 7, -0.5, 6.5);

      _histJetMultWmuPlus = bookHistogram1D("njetWmuPlus", 7, -0.5, 6.5);
      _histJetMultWmuMinus = bookHistogram1D("njetWmuMinus", 7, -0.5, 6.5);
      _histJetMultWelPlus = bookHistogram1D("njetWePlus", 7, -0.5, 6.5);
      _histJetMultWelMinus = bookHistogram1D("njetWeMinus", 7, -0.5, 6.5);

    } 

    
    bool ApplyElectronCutsForZee(double pt1, double pt2, double eta1, double eta2){
      bool isFid1 = ((fabs(eta1)<1.4442)||((fabs(eta1)>1.566)&&(fabs(eta1)<2.5)));
      bool isFid2 = ((fabs(eta2)<1.4442)||((fabs(eta2)>1.566)&&(fabs(eta2)<2.5)));
      if( isFid1 && isFid2 && pt1>20 && pt2 >10) return true;
      else return false;
    }

    
    bool ApplyMuonCutsForZmm(double pt1, double pt2, double eta1, double eta2){
      bool isFid1 = ((fabs(eta1)<2.1));
      bool isFid2 = ((fabs(eta2)<2.4));
      if( isFid1 && isFid2 && pt1>20 && pt2 >10) return true;
      else return false;
    }


    bool ApplyElectronCutsForWen(double pt1, double eta1){
      bool isFid1 = ((fabs(eta1)<1.4442)||((fabs(eta1)>1.566)&&(fabs(eta1)<2.5)));
      if( isFid1 && pt1>20 ) return true;
      return 0;
    }
 
   
    bool ApplyMuonCutsForWmn(double pt1, double eta1){
      bool isFid1 = ((fabs(eta1)<2.1));
      if( isFid1 && pt1>20) return true;
      return 0;
    }
    

    void Fill(AIDA::IHistogram1D*& _histJetMult, const double& weight, std::vector<FourMomentum>& finaljet_list){
      _histJetMult->fill(0, weight);
      for (size_t i=0 ; i<finaljet_list.size() ; ++i) {
        if (i==6) break;
        _histJetMult->fill(i+1, weight);  // inclusive
      }
    }  
    
    
    void FillHistogramSet(  AIDA::IHistogram1D*& _histJetMult, AIDA::IDataPointSet* _histJetMultNoverNm1[4], AIDA::IDataPointSet* _histJetMultNoverN0[4] ){
      for (int i = 0; i < 4; ++i) {
	std::vector<double> xval; xval.push_back(i+1);
	std::vector<double> xerr; xerr.push_back(.5);      
	double ratioNtoNminus1=0;       
	double ratioNtoNzero=0;      
	double frac_errNminus1 = 0.;   
	double frac_errNtoNzero = 0.;
	CalculateRatioAndErrors(_histJetMult, i, ratioNtoNminus1, ratioNtoNzero, frac_errNminus1, frac_errNtoNzero);       
/*
	std::cout<<"(njet >= "<<i+1<<")/(njet >="<<i<<") ratio: "<<ratioNtoNminus1<<"\t";	
	std::cout<<"(njet >= "<<i+1<<")/(njet >="<<0<<") ratio: "<<ratioNtoNzero<<"\n";	
*/	
	std::vector<double> yval;yval.push_back(ratioNtoNminus1);        
	std::vector<double> yerr;yerr.push_back(ratioNtoNminus1*frac_errNminus1);    
	std::vector<double> yval2;yval2.push_back(ratioNtoNzero);            
	std::vector<double> yerr2;yerr2.push_back(ratioNtoNzero*frac_errNtoNzero);    
        
	_histJetMultNoverNm1[i]->setCoordinate(0,xval,xerr);              
	_histJetMultNoverNm1[i]->setCoordinate(1,yval,yerr);           
	_histJetMultNoverN0[i]->setCoordinate(0,xval,xerr);         
	_histJetMultNoverN0[i]->setCoordinate(1,yval2,yerr2);        
        }                      
    }

    
   void FillChargeAssymHistogramSet(  AIDA::IHistogram1D*& _histJetMult1,AIDA::IHistogram1D*& _histJetMult2, AIDA::IDataPointSet* _histJetMultRatio12[4] ){
      for (int i = 0; i < 4; ++i) {
        std::vector<double> xval; xval.push_back(i);
        std::vector<double> xerr; xerr.push_back(.5);
        double ratio = 0;
        double err = 0.;
        double num = _histJetMult1->binHeight(i)-_histJetMult2->binHeight(i);
	double den = _histJetMult1->binHeight(i)+_histJetMult2->binHeight(i);
	double errNum = 0;
	errNum = std::pow(_histJetMult1->binError(i),2)+std::pow(_histJetMult2->binError(i),2);
	double errDen = 0;
	errDen = std::pow(_histJetMult1->binError(i),2)+std::pow(_histJetMult2->binError(i),2); 

        if (den)ratio = num/den;

        if(num)
	  errNum = errNum/(num*num); 
        if(den) 
	  errDen = errDen/(den*den);

        err = std::sqrt(errDen+errNum);
	if(!(err==err))err=0;
//	std::cout<<ratio<<"\t"<<err<<std::endl;
        std::vector<double> yval;yval.push_back(ratio);
        std::vector<double> yerr;yerr.push_back(ratio*err);
        _histJetMultRatio12[i]->setCoordinate(0,xval,xerr);
        _histJetMultRatio12[i]->setCoordinate(1,yval,yerr);
      }
    }


    void CalculateRatioAndErrors( AIDA::IHistogram1D*& _histJetMult, int i, double& ratioNtoNminus1, double& ratioNtoNzero, 
				  double& frac_errNminus1, double& frac_errNtoNzero ){       
      double errorTermForZeroJet = 0;
      if (_histJetMult->binHeight(i) > 0.)      
	ratioNtoNminus1 = _histJetMult->binHeight(i+1)/_histJetMult->binHeight(i);             
      if (_histJetMult->binHeight(0) > 0.) {
	ratioNtoNzero = _histJetMult->binHeight(i+1)/_histJetMult->binHeight(0);       
      }
      double errCont0 = 0, errCont1=0;      
      
      if (_histJetMult->binHeight(i) > 0.) {                                                                                                  
	double delMult0 = _histJetMult->binError(i);
	double mult0 = _histJetMult->binHeight(i);    
	errCont0 = delMult0/mult0;   
      }
      if(i==0)errorTermForZeroJet =  errCont0;
      
      if (_histJetMult->binHeight(i+1)>0 ){     
	double delMult = _histJetMult->binError(i+1);     
	double mult = _histJetMult->binHeight(i+1);     
	errCont1 = delMult/mult;      
      }     
      
      frac_errNminus1 = std::sqrt(errCont0*errCont0 + errCont1* errCont1);         
      frac_errNtoNzero = std::sqrt(errorTermForZeroJet*errorTermForZeroJet + errCont1* errCont1);      
    }    


    void analyze(const Event& event) {
      //some flag definitions.
      bool isZmm =false;
      bool isZee =false;
      bool isWmn =false;
      bool isWen =false;
      bool isWmnMinus =false;
      bool isWmnPlus  =false;
      bool isWenMinus =false;
      bool isWenPlus  =false;

      const double weight = event.weight();
      
      const InvMassFinalState& invMassFinalStateZ = applyProjection<InvMassFinalState>(event, "INVFSZ");
      const InvMassFinalState& invMassFinalStateW = applyProjection<InvMassFinalState>(event, "INVFSW");
      
      bool isW(false); bool isZ(false);
      
      isW  = (invMassFinalStateZ.empty() && !(invMassFinalStateW.empty()));
      isZ  = (!(invMassFinalStateZ.empty()) && invMassFinalStateW.empty());

      const ParticleVector&  ZDecayProducts =  invMassFinalStateZ.particles();
      const ParticleVector&  WDecayProducts =  invMassFinalStateW.particles();

      if (ZDecayProducts.size() < 2 && WDecayProducts.size() <2) vetoEvent;
      
      double pt1=-9999.,  pt2=-9999.;
      double phi1=-9999., phi2=-9999.;
      double eta1=-9999., eta2=-9999.;
      
      double mt = 999999;
      if(isZ){
	pt1  = ZDecayProducts[0].momentum().pT();
	pt2  = ZDecayProducts[1].momentum().pT();
	eta1 = ZDecayProducts[0].momentum().eta();
	eta2 = ZDecayProducts[1].momentum().eta();
	phi1 = ZDecayProducts[0].momentum().phi();
	phi2 = ZDecayProducts[1].momentum().phi();
      }
      
      if(isW){
	if(
	   (fabs(WDecayProducts[1].pdgId()) == NU_MU) || (fabs(WDecayProducts[1].pdgId()) == NU_E)){
	  pt1  = WDecayProducts[0].momentum().pT();
	  pt2  = WDecayProducts[1].momentum().Et();
          eta1 = WDecayProducts[0].momentum().eta();
          eta2 = WDecayProducts[1].momentum().eta();
	  phi1 = WDecayProducts[0].momentum().phi();
	  phi2 = WDecayProducts[1].momentum().phi();
          mt=sqrt(2.0*pt1*pt2*(1.0-cos(phi1-phi2)));
	}
	else {
	  pt1  = WDecayProducts[1].momentum().pT();
	  pt2  = WDecayProducts[0].momentum().Et();
          eta1 = WDecayProducts[1].momentum().eta();
          eta2 = WDecayProducts[0].momentum().eta();
	  phi1 = WDecayProducts[1].momentum().phi();
	  phi2 = WDecayProducts[0].momentum().phi();
          mt=sqrt(2.0*pt1*pt2*(1.0-cos(phi1-phi2)));
	}
      }

      if(isW && mt<20)vetoEvent;
            
      isZmm = isZ && ((fabs(ZDecayProducts[0].pdgId()) == 13) && (fabs(ZDecayProducts[1].pdgId()) == 13));
      isZee = isZ && ((fabs(ZDecayProducts[0].pdgId()) == 11) && (fabs(ZDecayProducts[1].pdgId()) == 11));
      isWmn  = isW && ((fabs(WDecayProducts[0].pdgId()) == 14) || (fabs(WDecayProducts[1].pdgId()) == 14));
      isWen  = isW && ((fabs(WDecayProducts[0].pdgId()) == 12) || (fabs(WDecayProducts[1].pdgId()) == 12));
      
      if(isWmn){
        if((WDecayProducts[0].pdgId()==-13)|| (WDecayProducts[1].pdgId()==-13)){
	isWmnMinus = false;
	isWmnPlus = true;
	  }	
       else{
	isWmnMinus = true;
	isWmnPlus  = false;
      } 
     }

      if(isWen){
       if((WDecayProducts[0].pdgId()==11)|| (WDecayProducts[1].pdgId()==11)){
	isWenMinus = true;
	isWenPlus = false;
      }
      else{
	isWenMinus = false;
	isWenPlus  = true;
       }
      }

      if(!((isZmm||isZee)||(isWmn||isWen)))vetoEvent;
              
      bool passBosonConditions = false;
      if(isZmm)passBosonConditions = ApplyMuonCutsForZmm(pt1,pt2,eta1,eta2);
      if(isZee)passBosonConditions = ApplyElectronCutsForZee(pt1,pt2,eta1,eta2);
      if(isWen)passBosonConditions = ApplyElectronCutsForWen(pt1,eta1);  
      if(isWmn)passBosonConditions = ApplyMuonCutsForWmn(pt1,eta1);  
      
      if(!passBosonConditions)vetoEvent;
          
      //Obtain the jets.
      vector<FourMomentum> finaljet_list;
      foreach (const Jet& j, applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV)) {
	const double jeta = j.momentum().eta();
	const double jphi = j.momentum().phi();
	const double jpt = j.momentum().pT();
	if (fabs(jeta) < 2.4) 
	  if(jpt>30){
	      if(isZee){
		  if (deltaR(pt1, phi1, jeta, jphi) > 0.3 && deltaR(pt2, phi2, jeta, jphi) > 0.3)
		    finaljet_list.push_back(j.momentum());
		  continue;
		}
	      else if(isWen){
		  if (deltaR(pt1, phi1, jeta, jphi) > 0.3)
		    finaljet_list.push_back(j.momentum());
		  continue;
	      }
	      
	      else  finaljet_list.push_back(j.momentum());
	  }
      }
      //Multiplicity plots.	
      if(isWen)Fill(_histJetMultWelec, weight, finaljet_list);
      if(isWmn)Fill(_histJetMultWmu, weight, finaljet_list);
      if(isWmnPlus)Fill(_histJetMultWmuPlus, weight, finaljet_list);
      if(isWmnMinus)Fill(_histJetMultWmuMinus, weight, finaljet_list);
      if(isWenPlus)Fill(_histJetMultWelPlus, weight, finaljet_list);
      if(isWenMinus)Fill(_histJetMultWelMinus, weight, finaljet_list);
      if(isZee)Fill(_histJetMultZelec, weight, finaljet_list);
      if(isZmm)Fill(_histJetMultZmu, weight, finaljet_list);
    }
    
    
    /// Normalise histograms etc., after the run
    void finalize() {
//      std::cout<<"For the Z->mm case: \n";
      FillHistogramSet(_histJetMultZmu,  _histJetMultNoverNm1Zmu,  _histJetMultNoverN0Zmu);
//      std::cout<<"For the Z->ee case: \n";
      FillHistogramSet(_histJetMultZelec,_histJetMultNoverNm1Zelec,_histJetMultNoverN0Zelec);
//      std::cout<<"For the W->en case: \n";
      FillHistogramSet(_histJetMultWelec,_histJetMultNoverNm1Welec,_histJetMultNoverN0Welec);
//      std::cout<<"For the W->mn case: \n";
      FillHistogramSet(_histJetMultWmu,_histJetMultNoverNm1Wmu,_histJetMultNoverN0Wmu);
      FillChargeAssymHistogramSet(_histJetMultWmuPlus,_histJetMultWmuMinus, _histJetMultRatioWmuPlusMinus);
      FillChargeAssymHistogramSet(_histJetMultWelPlus,_histJetMultWelMinus, _histJetMultRatioWelPlusMinus);
    }

  private:

    AIDA::IDataPointSet* _histJetMultNoverNm1Welec[4];
    AIDA::IDataPointSet* _histJetMultNoverN0Welec[4];  
    AIDA::IHistogram1D*  _histJetMultWelec;
    
    AIDA::IDataPointSet* _histJetMultNoverNm1Wmu[4];
    AIDA::IDataPointSet* _histJetMultNoverN0Wmu[4];  
    AIDA::IHistogram1D*  _histJetMultWmu;

    AIDA::IHistogram1D*  _histJetMultWelMinus;
    AIDA::IHistogram1D*  _histJetMultWelPlus;
    AIDA::IDataPointSet* _histJetMultRatioWelPlusMinus[4];
    
    AIDA::IHistogram1D*  _histJetMultWmuMinus;
    AIDA::IHistogram1D*  _histJetMultWmuPlus;
    AIDA::IDataPointSet* _histJetMultRatioWmuPlusMinus[4];
   
    AIDA::IDataPointSet* _histJetMultNoverNm1Zelec[4];
    AIDA::IDataPointSet* _histJetMultNoverN0Zelec[4];  
    AIDA::IHistogram1D*  _histJetMultZelec;

    AIDA::IDataPointSet* _histJetMultNoverNm1Zmu[4];
    AIDA::IDataPointSet* _histJetMultNoverN0Zmu[4];  
    AIDA::IHistogram1D*  _histJetMultZmu;
  };
  
  AnalysisBuilder<CMS_EWK_10_012> plugin_CMS_EWK_10_012;
  
}

